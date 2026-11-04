#include "memoscip_reader.h"

#include "scip/scip.h"
#include "scip/scip_cons.h"
#include "scip/cons_linear.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <lzma.h>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace std;

static inline string trim(const string& s)
{
    size_t a = 0;
    while (a < s.size() && isspace((unsigned char)s[a])) ++a;
    size_t b = s.size();
    while (b > a && isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

// Custom streambuf for LZMA decompression
class LZMAStreambuf : public std::streambuf {
private:
    FILE* file_;
    lzma_stream stream_;
    std::vector<uint8_t> in_buf_;
    std::vector<uint8_t> out_buf_;
    static const size_t BUF_SIZE = 65536;

    int fill_buffer() {
        if (in_buf_.empty()) in_buf_.resize(BUF_SIZE);
        if (out_buf_.empty()) out_buf_.resize(BUF_SIZE);

        size_t nread = fread(in_buf_.data(), 1, BUF_SIZE, file_);
        if (nread == 0 && ferror(file_)) return -1;

        stream_.next_in = in_buf_.data();
        stream_.avail_in = nread;
        stream_.next_out = out_buf_.data();
        stream_.avail_out = out_buf_.size();

        lzma_ret ret = lzma_code(&stream_, nread == 0 ? LZMA_FINISH : LZMA_RUN);
        if (ret != LZMA_OK && ret != LZMA_STREAM_END) return -1;

        size_t produced = out_buf_.size() - stream_.avail_out;
        setg((char*)out_buf_.data(), (char*)out_buf_.data(), (char*)out_buf_.data() + produced);

        return produced > 0 ? *eback() : EOF;
    }

public:
    LZMAStreambuf(FILE* f) : file_(f), stream_(LZMA_STREAM_INIT) {
        lzma_ret ret = lzma_stream_decoder(&stream_, UINT64_MAX, LZMA_CONCATENATED);
        if (ret != LZMA_OK) throw std::runtime_error("lzma_stream_decoder failed");
        fill_buffer();
    }

    ~LZMAStreambuf() {
        lzma_end(&stream_);
    }

    int underflow() override {
        if (gptr() < egptr()) return *gptr();
        return fill_buffer();
    }
};

class LZMAIStream : public std::istream {
private:
    std::unique_ptr<LZMAStreambuf> buf_;
public:
    LZMAIStream(FILE* f) : std::istream(nullptr) {
        buf_ = std::make_unique<LZMAStreambuf>(f);
        rdbuf(buf_.get());
    }
};

/** Opens a plain or .xz compressed instance file. */
static std::unique_ptr<std::istream> open_instance_file(const std::string& path, std::string& err) {
    if (path.length() >= 3 && path.substr(path.length() - 3) == ".xz") {
        FILE* f = fopen(path.c_str(), "rb");
        if (!f) {
            err = "Failed to open .xz file: " + path;
            return nullptr;
        }
        try {
            return std::make_unique<LZMAIStream>(f);
        } catch (const std::exception& e) {
            fclose(f);
            err = std::string("LZMA error: ") + e.what();
            return nullptr;
        }
    } else {
        auto ifs = std::make_unique<std::ifstream>(path);
        if (!ifs->is_open()) {
            err = "Failed to open file: " + path;
            return nullptr;
        }
        return ifs;
    }
}

/** Creates and adds a binary variable; on error sets err and returns false. */
static bool addBinaryVar(SCIP* scip, const char* name, SCIP_Real obj, SCIP_VAR*& var, std::string& err)
{
    if (SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, obj, SCIP_VARTYPE_BINARY) != SCIP_OKAY) { err = "SCIP error creating var"; return false; }
    if (SCIPaddVar(scip, var) != SCIP_OKAY) { err = "SCIP error adding var"; return false; }
    return true;
}

/** Creates, adds and releases the linear constraint lhs <= coefs.vars <= rhs; on error sets err
 *  ("SCIP error creating/adding/releasing <what>") and returns false. */
static bool addLinearCons(SCIP* scip, const char* name, int nvars, SCIP_VAR** vars, SCIP_Real* coefs,
    SCIP_Real lhs, SCIP_Real rhs, const std::string& what, std::string& err)
{
    SCIP_CONS* cons = nullptr;
    if (SCIPcreateConsBasicLinear(scip, &cons, name, nvars, vars, coefs, lhs, rhs) != SCIP_OKAY) { err = "SCIP error creating " + what; return false; }
    if (SCIPaddCons(scip, cons) != SCIP_OKAY) { err = "SCIP error adding " + what; return false; }
    if (SCIPreleaseCons(scip, &cons) != SCIP_OKAY) { err = "SCIP error releasing " + what; return false; }
    return true;
}


bool readKPCGInstance(SCIP* scip, const std::string& path, int encoding, std::string& err)
{
    if (scip == nullptr) {
        err = "SCIP pointer is null";
        return false;
    }

    auto stream = open_instance_file(path, err);
    if (!stream) {
        return false;
    }
    std::istream& ifs = *stream;

    int n_items = -1;
    int capacity = -1;
    vector<tuple<int,int,int>> items;  // (index, profit, weight)
    vector<pair<int,int>> edges;       // conflicts (i, j) with i < j
    string mode;

    regex param_re(R"(^\s*param\s+(\w+)\s*:=\s*(\d+)\s*;\s*$)");

    string raw;
    while (std::getline(ifs, raw)) {
        string line = trim(raw);
        if (line.empty() || line[0] == '#')
            continue;

        smatch m;
        if (regex_match(line, m, param_re)) {
            string name = m[1];
            int val = stoi(m[2]);
            if (name == "n") n_items = val;
            else if (name == "c") capacity = val;
            continue;
        }

        if (line.rfind("param :", 0) == 0) {
            mode = "items";
            continue;
        }
        if (line.rfind("set E", 0) == 0) {
            mode = "edges";
            continue;
        }

        if (line == ";") {
            if (mode == "items") mode.clear();
            else if (mode == "edges") break;
            continue;
        }

        if (mode == "items") {
            // parse: idx profit weight
            istringstream iss(line);
            int idx, profit, weight;
            if (!(iss >> idx >> profit >> weight))
                continue;
            items.emplace_back(idx, profit, weight);
        }
        else if (mode == "edges") {
            istringstream iss(line);
            int i, j;
            if (!(iss >> i >> j))
                continue;
            if (i > j) swap(i, j);
            edges.emplace_back(i, j);
        }
    }

    if (n_items < 0 || capacity < 0) {
        err = "Missing param n and/or param c in instance file";
        return false;
    }

    if ((int)items.size() != n_items) {
        err = "Number of items does not match param n";
        return false;
    }

    // validate indices
    vector<int> indices;
    for (auto &t : items) indices.push_back(get<0>(t));
    sort(indices.begin(), indices.end());
    for (int i = 0; i < n_items; ++i) {
        if (indices[i] != i) {
            err = "Item indices must be 0..n-1";
            return false;
        }
    }

    // Create SCIP variables: x<i>_0 = item i not taken (pays its profit), x<i>_1 = item i taken
    vector<SCIP_VAR*> var0(n_items, nullptr);
    vector<SCIP_VAR*> var1(n_items, nullptr);

    int total_weight = 0;
    for (auto &t : items) total_weight += get<2>(t);
    int knapsack_rhs = total_weight - capacity;

    for (int idx = 0; idx < n_items; ++idx) {
        int profit = get<1>(items[idx]);

        char name0[64];
        char name1[64];
        snprintf(name0, sizeof(name0), "x%d_0", idx);
        snprintf(name1, sizeof(name1), "x%d_1", idx);

        if (!addBinaryVar(scip, name0, (SCIP_Real)profit, var0[idx], err)) return false;
        if (!addBinaryVar(scip, name1, 0.0, var1[idx], err)) return false;

        // the problem keeps a reference to the variables, the pointers stay valid
        SCIP_VAR* v0 = var0[idx];
        SCIP_VAR* v1 = var1[idx];
        SCIPreleaseVar(scip, &v0);
        SCIPreleaseVar(scip, &v1);
    }

    // knapsack constraint: sum weight * x0 >= knapsack_rhs
    {
        vector<SCIP_Real> coefs;
        coefs.reserve(n_items);
        for (int idx = 0; idx < n_items; ++idx)
            coefs.push_back((SCIP_Real)get<2>(items[idx]));
        if (!addLinearCons(scip, "knapsack", n_items, var0.data(), coefs.data(),
                           (SCIP_Real)knapsack_rhs, SCIPinfinity(scip), "knapsack constraint", err))
            return false;
    }

    // Conflict edges use either the original direct formulation or a tuple formulation.
    for (size_t k = 0; k < edges.size(); ++k) {
        int i = edges[k].first;
        int j = edges[k].second;
        if (i < 0 || i >= n_items || j < 0 || j >= n_items || i == j) {
            err = "Conflict edge contains an invalid item index";
            return false;
        }

        char cname[128];
        if (encoding == 0) {
            // x<i>_1 + x<j>_1 <= 1
            SCIP_VAR* consvars[2] = { var1[i], var1[j] };
            SCIP_Real coefs[2] = { 1.0, 1.0 };
            snprintf(cname, sizeof(cname), "conflict_%zu", k);
            if (!addLinearCons(scip, cname, 2, consvars, coefs, -SCIPinfinity(scip), 1.0, "conflict constraint", err))
                return false;
            continue;
        }

        // Tuple formulation: one variable per assignment of (x<i>, x<j>), the forbidden tuple 11 is fixed to 0
        SCIP_VAR* tuples[4] = { nullptr, nullptr, nullptr, nullptr };
        const char* tuple_labels[4] = { "00", "11", "10", "01" };
        for (int tuple = 0; tuple < 4; ++tuple) {
            char name[128];
            snprintf(name, sizeof(name), "x%d_%s_%d", i, tuple_labels[tuple], j);

            SCIP_VAR* tuple_var = nullptr;
            SCIP_RETCODE ret = SCIPcreateVarBasic(scip, &tuple_var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY);
            if (ret != SCIP_OKAY) { err = "SCIP error creating tuple variable"; return false; }
            ret = SCIPaddVar(scip, tuple_var);
            if(tuple == 1) {
                SCIP_CALL( SCIPchgVarUb(scip, tuple_var, 0) );
            }
            if (ret != SCIP_OKAY) { err = "SCIP error adding tuple variable"; return false; }
            tuples[tuple] = tuple_var;
            ret = SCIPreleaseVar(scip, &tuple_var);
            if (ret != SCIP_OKAY) { err = "SCIP error releasing tuple variable"; return false; }
        }

        // Marginal constraints: each literal is the sum of the tuples containing it
        SCIP_Real marginal_coefs[3] = { 1.0, 1.0, -1.0 };
        SCIP_VAR* marginal_vars[4][3] = {
            { tuples[0], tuples[3], var0[i] },
            { tuples[2], tuples[1], var1[i] },
            { tuples[0], tuples[2], var0[j] },
            { tuples[3], tuples[1], var1[j] }
        };
        const char* marginal_labels[4] = { "i0", "i1", "j0", "j1" };
        for (int marginal = 0; marginal < 4; ++marginal) {
            // i1 and j1 are implied by i0/j0, the tuple sum and the domain constraints (i1 = sum - i0 - domain_i):
            // keeping them makes the equality rows linearly dependent and the LP duals non-unique (huge dual values)
            if (marginal == 1 || marginal == 3)
                continue;
            snprintf(cname, sizeof(cname), "tuple_%zu_%s", k, marginal_labels[marginal]);
            if (!addLinearCons(scip, cname, 3, marginal_vars[marginal], marginal_coefs, 0.0, 0.0, "tuple marginal constraint", err))
                return false;
        }

        // Exactly one tuple is selected
        SCIP_Real tuple_sum_coefs[4] = { 1.0, 1.0, 1.0, 1.0 };
        snprintf(cname, sizeof(cname), "tuple_%zu_sum", k);
        if (!addLinearCons(scip, cname, 4, tuples, tuple_sum_coefs, 1.0, 1.0, "tuple sum constraint", err))
            return false;
    }

    // domain constraints: x0 + x1 = 1
    for (int idx = 0; idx < n_items; ++idx) {
        SCIP_VAR* consvars[2] = { var0[idx], var1[idx] };
        SCIP_Real coefs[2] = { 1.0, 1.0 };
        char cname[128];
        snprintf(cname, sizeof(cname), "domain_constraint_%d", idx);
        if (!addLinearCons(scip, cname, 2, consvars, coefs, 1.0, 1.0, "domain constraint", err))
            return false;
    }

    return true;
}

/** Returns the value of attribute `name` in the XML start tag `tag`, or an empty string. */
static string xmlAttr(const string& tag, const string& name)
{
    string key = " " + name + "=\"";
    size_t a = tag.find(key);
    if (a == string::npos) return "";
    a += key.size();
    size_t b = tag.find('"', a);
    return b == string::npos ? "" : tag.substr(a, b - a);
}

/** Calls f(starttag, body) for each element <elem ...>body</elem> or <elem .../> of the document. */
template <class F>
static void forEachXmlElement(const string& doc, const string& elem, F&& f)
{
    string open = "<" + elem + " ";
    string close = "</" + elem + ">";
    size_t pos = 0;
    while ((pos = doc.find(open, pos)) != string::npos) {
        size_t end = doc.find('>', pos);
        if (end == string::npos) return;
        string tag = doc.substr(pos, end - pos + 1);
        string body;
        if (doc[end - 1] != '/') {
            size_t cend = doc.find(close, end);
            if (cend == string::npos) return;
            body = doc.substr(end + 1, cend - end - 1);
            end = cend;
        }
        f(tag, body);
        pos = end + 1;
    }
}

bool readXCSPInstance(SCIP* scip, const std::string& path, int encoding, std::string& err,
    std::vector<std::string>* varnames_out)
{
    if (scip == nullptr) {
        err = "SCIP pointer is null";
        return false;
    }

    auto stream = open_instance_file(path, err);
    if (!stream) {
        return false;
    }
    string doc((std::istreambuf_iterator<char>(*stream)), std::istreambuf_iterator<char>());

    // domains: list of integers and/or intervals a..b
    unordered_map<string, vector<int>> domains;
    forEachXmlElement(doc, "domain", [&](const string& tag, const string& body) {
        vector<int>& values = domains[xmlAttr(tag, "name")];
        istringstream iss(body);
        string tok;
        while (iss >> tok) {
            size_t dots = tok.find("..");
            if (dots == string::npos) values.push_back(stoi(tok));
            else for (int v = stoi(tok.substr(0, dots)); v <= stoi(tok.substr(dots + 2)); ++v) values.push_back(v);
        }
    });

    // variables
    vector<string> varnames;
    vector<const vector<int>*> vardom;
    unordered_map<string, int> varidx;
    bool ok = true;
    forEachXmlElement(doc, "variable", [&](const string& tag, const string&) {
        auto it = domains.find(xmlAttr(tag, "domain"));
        if (it == domains.end()) { err = "Unknown domain for variable " + xmlAttr(tag, "name"); ok = false; return; }
        varidx[xmlAttr(tag, "name")] = (int)varnames.size();
        varnames.push_back(xmlAttr(tag, "name"));
        vardom.push_back(&it->second);
    });
    if (!ok) return false;
    if (varnames.empty()) {
        err = "No variables found in instance file";
        return false;
    }
    if (varnames_out) *varnames_out = varnames;

    // relations: list of value pairs, either supports or conflicts
    struct Relation { bool supports; vector<pair<int,int>> tuples; };
    unordered_map<string, Relation> relations;
    forEachXmlElement(doc, "relation", [&](const string& tag, const string& body) {
        if (xmlAttr(tag, "arity") != "2") { err = "Only binary relations are supported"; ok = false; return; }
        Relation& rel = relations[xmlAttr(tag, "name")];
        rel.supports = xmlAttr(tag, "semantics") == "supports";
        string tuples = body;
        replace(tuples.begin(), tuples.end(), '|', ' ');
        istringstream iss(tuples);
        int a, b;
        while (iss >> a >> b) rel.tuples.emplace_back(a, b);
    });
    if (!ok) return false;

    // x<i>_<v> = 1 iff variable i takes value v
    vector<vector<SCIP_VAR*>> vars(varnames.size());
    vector<unordered_map<int,int>> valpos(varnames.size());
    for (size_t i = 0; i < varnames.size(); ++i) {
        for (int v : *vardom[i]) {
            char name[64];
            snprintf(name, sizeof(name), "x%zu_%d", i, v);
            SCIP_VAR* var = nullptr;
            if (!addBinaryVar(scip, name, 0.0, var, err)) return false;
            valpos[i][v] = (int)vars[i].size();
            vars[i].push_back(var);
            // the problem keeps a reference to the variable, the pointer stays valid
            SCIPreleaseVar(scip, &var);
        }
    }

    // constraints: encoding 0 adds x<i>_<a> + x<j>_<b> <= 1 for each forbidden pair (a, b) (hard CSP), otherwise
    // one tuple variable per pair of values is created and the forbidden ones cost 1 (Max-CSP: minimise the number
    // of violated tuples)
    int k = 0;
    forEachXmlElement(doc, "constraint", [&](const string& tag, const string&) {
        if (!ok) return;
        istringstream scope(xmlAttr(tag, "scope"));
        string sx, sy, extra;
        if (!(scope >> sx >> sy) || (scope >> extra)) { err = "Only binary constraints are supported"; ok = false; return; }
        auto itx = varidx.find(sx), ity = varidx.find(sy);
        auto itr = relations.find(xmlAttr(tag, "reference"));
        if (itx == varidx.end() || ity == varidx.end() || itx->second == ity->second || itr == relations.end()) {
            err = "Invalid scope or reference in constraint " + xmlAttr(tag, "name"); ok = false; return;
        }
        int i = itx->second, j = ity->second;
        const Relation& rel = itr->second;
        size_t di = vars[i].size(), dj = vars[j].size();

        // listed[a][b]: pair (value a of i, value b of j) appears in the relation; values outside the domains are ignored
        vector<vector<bool>> listed(di, vector<bool>(dj, false));
        for (auto& t : rel.tuples) {
            auto pa = valpos[i].find(t.first), pb = valpos[j].find(t.second);
            if (pa != valpos[i].end() && pb != valpos[j].end()) listed[pa->second][pb->second] = true;
        }
        auto forbidden = [&](size_t a, size_t b) { return listed[a][b] != rel.supports; };

        char cname[128];
        if (encoding == 0) {
            for (size_t a = 0; a < di; ++a)
                for (size_t b = 0; b < dj; ++b) {
                    if (!forbidden(a, b)) continue;
                    SCIP_VAR* consvars[2] = { vars[i][a], vars[j][b] };
                    SCIP_Real coefs[2] = { 1.0, 1.0 };
                    snprintf(cname, sizeof(cname), "conflict_%d_%zu_%zu", k, a, b);
                    if (!addLinearCons(scip, cname, 2, consvars, coefs, -SCIPinfinity(scip), 1.0, "conflict constraint", err)) { ok = false; return; }
                }
            ++k;
            return;
        }

        // Tuple formulation: t[a][b] = 1 iff x<i> takes its a-th value and x<j> its b-th value, forbidden tuples cost 1
        vector<vector<SCIP_VAR*>> t(di, vector<SCIP_VAR*>(dj, nullptr));
        vector<SCIP_VAR*> alltuples;
        for (size_t a = 0; a < di; ++a)
            for (size_t b = 0; b < dj; ++b) {
                char name[128];
                snprintf(name, sizeof(name), "x%d_%d_x%d_%d", i, (*vardom[i])[a], j, (*vardom[j])[b]);
                if (!addBinaryVar(scip, name, forbidden(a, b) ? 1.0 : 0.0, t[a][b], err)) { ok = false; return; }
                alltuples.push_back(t[a][b]);
                SCIP_VAR* tv = t[a][b];
                if (SCIPreleaseVar(scip, &tv) != SCIP_OKAY) { err = "SCIP error releasing tuple variable"; ok = false; return; }
            }

        // Marginal constraints: each literal is the sum of the tuples containing it. The last value of each variable
        // is implied by the other marginals, the tuple sum and the domain constraint: keeping it makes the equality
        // rows linearly dependent and the LP duals non-unique
        for (size_t a = 0; a + 1 < di; ++a) {
            vector<SCIP_VAR*> mvars(t[a].begin(), t[a].end());
            mvars.push_back(vars[i][a]);
            vector<SCIP_Real> mcoefs(dj, 1.0);
            mcoefs.push_back(-1.0);
            snprintf(cname, sizeof(cname), "tuple_%d_i%zu", k, a);
            if (!addLinearCons(scip, cname, (int)mvars.size(), mvars.data(), mcoefs.data(), 0.0, 0.0, "tuple marginal constraint", err)) { ok = false; return; }
        }
        for (size_t b = 0; b + 1 < dj; ++b) {
            vector<SCIP_VAR*> mvars;
            for (size_t a = 0; a < di; ++a) mvars.push_back(t[a][b]);
            mvars.push_back(vars[j][b]);
            vector<SCIP_Real> mcoefs(di, 1.0);
            mcoefs.push_back(-1.0);
            snprintf(cname, sizeof(cname), "tuple_%d_j%zu", k, b);
            if (!addLinearCons(scip, cname, (int)mvars.size(), mvars.data(), mcoefs.data(), 0.0, 0.0, "tuple marginal constraint", err)) { ok = false; return; }
        }

        // Exactly one tuple is selected
        vector<SCIP_Real> sumcoefs(alltuples.size(), 1.0);
        snprintf(cname, sizeof(cname), "tuple_%d_sum", k);
        if (!addLinearCons(scip, cname, (int)alltuples.size(), alltuples.data(), sumcoefs.data(), 1.0, 1.0, "tuple sum constraint", err)) { ok = false; return; }
        ++k;
    });
    if (!ok) return false;

    // domain constraints: each variable takes exactly one value
    for (size_t i = 0; i < vars.size(); ++i) {
        vector<SCIP_Real> coefs(vars[i].size(), 1.0);
        char cname[128];
        snprintf(cname, sizeof(cname), "domain_constraint_%zu", i);
        if (!addLinearCons(scip, cname, (int)vars[i].size(), vars[i].data(), coefs.data(), 1.0, 1.0, "domain constraint", err))
            return false;
    }

    return true;
}
