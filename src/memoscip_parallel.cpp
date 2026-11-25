#include "memoscip_parallel.h"
#include "scip/cons_linear.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <unordered_map>

namespace {

const double PARALLEL_TAU = 1e-4;   // rejected if |cos| > 1 - PARALLEL_TAU

bool initialised = false;
int ncoord = 0;
std::unordered_map<SCIP_VAR*, int> coordOf;                        // original (non-slack) variable -> coordinate
std::unordered_map<SCIP_VAR*, std::vector<double>> slackDef;       // slack -> definition over the coordinates
std::vector<std::vector<int>> blocks;                              // domain blocks
std::vector<std::vector<double>> kept;                             // kept learned cuts, canonical and normalised

std::string pendingName;
std::vector<double> pending;
bool hasPending = false;

long long nChecked = 0, nRejected = 0, nSkipped = 0;
double checkTime = 0.0;

bool isModelSlack(SCIP_VAR* var) { return strstr(SCIPvarGetName(var), "slack_") != nullptr; }

/** out += c * var, the slacks being replaced by their definition. Returns false for an unknown variable. */
bool addTerm(std::vector<double>& out, SCIP_VAR* var, double c)
{
    auto it = coordOf.find(var);
    if( it != coordOf.end() )
    {
        out[it->second] += c;
        return true;
    }
    auto s = slackDef.find(var);
    if( s == slackDef.end() )
        return false;
    const std::vector<double>& d = s->second;
    for( int j = 0; j < ncoord; ++j )
        if( d[j] != 0.0 )
            out[j] += c * d[j];
    return true;
}

/** Centres each domain block and normalises. Returns false for a zero vector. */
bool canonicalise(std::vector<double>& v)
{
    for( const std::vector<int>& b : blocks )
    {
        double mean = 0.0;
        for( int j : b )
            mean += v[j];
        mean /= b.size();
        for( int j : b )
            v[j] -= mean;
    }
    double n2 = 0.0;
    for( double x : v )
        n2 += x * x;
    if( n2 <= 1e-24 )
        return false;
    const double inv = 1.0 / std::sqrt(n2);
    for( double& x : v )
        x *= inv;
    return true;
}

/** Coordinates, domain blocks and model slack definitions, from the transformed problem (no presolving). */
void initialise(SCIP* scip)
{
    if( initialised )
        return;
    initialised = true;

    SCIP_VAR** vars = SCIPgetVars(scip);
    for( int v = 0; v < SCIPgetNVars(scip); ++v )
        if( !isModelSlack(vars[v]) )
            coordOf[vars[v]] = ncoord++;

    SCIP_CONS** conss = SCIPgetConss(scip);
    const int nconss = SCIPgetNConss(scip);
    std::vector<char> inBlock(ncoord, 0);
    for( int c = 0; c < nconss; ++c )
    {
        if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
            continue;
        const SCIP_Real lhs = SCIPgetLhsLinear(scip, conss[c]);
        const SCIP_Real rhs = SCIPgetRhsLinear(scip, conss[c]);
        if( !SCIPisEQ(scip, lhs, rhs) )
            continue;
        const int nvars = SCIPgetNVarsLinear(scip, conss[c]);
        SCIP_VAR** cvars = SCIPgetVarsLinear(scip, conss[c]);
        SCIP_Real* cvals = SCIPgetValsLinear(scip, conss[c]);

        // Domain block: sum_v x_v = 1, binaries, coefficient 1, each variable in one block only
        if( SCIPisEQ(scip, rhs, 1.0) && nvars >= 2 )
        {
            std::vector<int> b;
            bool ok = true;
            for( int k = 0; k < nvars && ok; ++k )
            {
                auto it = coordOf.find(cvars[k]);
                ok = it != coordOf.end() && cvals[k] == 1.0 && SCIPvarIsBinary(cvars[k]) && !inBlock[it->second];
                if( ok )
                    b.push_back(it->second);
            }
            if( ok )
            {
                for( int j : b )
                    inBlock[j] = 1;
                blocks.push_back(std::move(b));
                continue;
            }
        }

        // Model slack: a.x + sigma*s = rhs  ->  s = -a.x / sigma + constant (constants are irrelevant for the cosine)
        int slackPos = -1;
        for( int k = 0; k < nvars; ++k )
            if( isModelSlack(cvars[k]) )
                slackPos = k;
        if( slackPos < 0 )
            continue;
        std::vector<double> def(ncoord, 0.0);
        bool ok = true;
        for( int k = 0; k < nvars && ok; ++k )
            if( k != slackPos )
                ok = addTerm(def, cvars[k], -cvals[k] / cvals[slackPos]);
        if( ok )
            slackDef[cvars[slackPos]] = std::move(def);
    }
}

}  // namespace

bool parallelCheck(SCIP* scip, const std::string& cutname, const std::vector<SCIP_VAR*>& consvars,
                   const std::vector<SCIP_Real>& consvals)
{
    initialise(scip);
    const auto t0 = std::chrono::steady_clock::now();
    hasPending = false;

    std::vector<double> v(ncoord, 0.0);
    bool ok = true;
    for( size_t k = 0; k < consvars.size() && ok; ++k )
        ok = addTerm(v, consvars[k], consvals[k]);

    bool reject = false;
    if( !ok || !canonicalise(v) )
        nSkipped++;  // unknown variable, or constant in the canonical space: kept, not compared
    else
    {
        nChecked++;
        for( const std::vector<double>& p : kept )
        {
            double d = 0.0;
            for( int j = 0; j < ncoord; ++j )
                d += v[j] * p[j];
            if( std::fabs(d) > 1.0 - PARALLEL_TAU )
            {
                reject = true;
                break;
            }
        }
        if( reject )
            nRejected++;
        else
        {
            pending = std::move(v);
            pendingName = cutname;
            hasPending = true;
        }
    }
    checkTime += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return reject;
}

void parallelRegister(SCIP* scip, const std::string& cutname, SCIP_VAR* slack,
                      const std::vector<SCIP_VAR*>& consvars, const std::vector<SCIP_Real>& consvals)
{
    initialise(scip);
    const auto t0 = std::chrono::steady_clock::now();

    // consvals . consvars = rhs with the slack coefficient c_s: s = -(a.x)/c_s + constant
    SCIP_Real slackCoef = 0.0;
    for( size_t k = 0; k < consvars.size(); ++k )
        if( consvars[k] == slack )
            slackCoef = consvals[k];
    if( slackCoef != 0.0 )
    {
        std::vector<double> def(ncoord, 0.0);
        bool ok = true;
        for( size_t k = 0; k < consvars.size() && ok; ++k )
            if( consvars[k] != slack )
                ok = addTerm(def, consvars[k], -consvals[k] / slackCoef);
        if( ok )
            slackDef[slack] = std::move(def);
    }

    if( hasPending && pendingName == cutname )
    {
        kept.push_back(std::move(pending));
        hasPending = false;
    }
    checkTime += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void parallelPrintStatistics()
{
    std::cout << "Parallel check : " << nChecked << " cuts checked, " << nRejected << " rejected, " << nSkipped
              << " not checked, " << kept.size() << " kept, " << blocks.size() << " domain blocks, time " << checkTime
              << " s" << std::endl;
}
