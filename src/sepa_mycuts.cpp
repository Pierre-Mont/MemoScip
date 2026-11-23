#include "sepa_mycuts.h"

#include "memoscip_option.h"
#include "memoscip_types.h"

#include "scip/cons_linear.h"
#include "scip/cuts.h"
#include "memoscip_dual.h"
#include "memoscip_events.h"
#include "memoscip_parallel.h"
#include "scip/var.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#define SEPA_NAME              "mycuts"
#define SEPA_DESC              "separator for queued cuts from events"
#define SEPA_PRIORITY          1000000
#define SEPA_FREQ              1
#define SEPA_MAXBOUNDDIST      1.0
#define SEPA_USESSUBSCIP       FALSE
#define SEPA_DELAY             FALSE

static
SCIP_DECL_SEPAEXECLP(sepaExeclpMycuts)
{
   if( g_memoscip.CutQueue.queue.empty() )
   {
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }
   if(g_options.verbosity > 1)
      std::cout<<"enter sepaExeclpMycuts with "<<g_memoscip.CutQueue.queue.size()<<" cuts in queue"<<std::endl;
   SCIP_CALL(AddCutstoProb(scip, &g_memoscip.CutQueue, &g_memoscip, double(g_memoscip.TOP), true, &g_memoscip.nbctr, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)),sepa));
   SCIP_CALL( SCIPflushLP(scip) );
   *result = SCIP_SEPARATED;
   return SCIP_OKAY;
}

void includeMycutsSeparator(SCIP* scip)
{
   SCIP_SEPA* sepa = nullptr;

   SCIP_CALL_ABORT(SCIPincludeSepaBasic(
      scip,
      &sepa,
      SEPA_NAME,
      SEPA_DESC,
      SEPA_PRIORITY,
      SEPA_FREQ,
      SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP,
      SEPA_DELAY,
      sepaExeclpMycuts,
      nullptr,
      nullptr));
   SCIP_CALL_ABORT(SCIPsetIntParam(scip,"separating/" SEPA_NAME "/expbackoff",1));
}

static double roundToDecimals(double val, int decimals) {
   double factor = std::pow(10.0, decimals);
   return std::round(val * factor) / factor;
}

/** Ceiling that tolerates values slightly above an integer. */
SCIP_Longint Myceil(double val){
   if(val - floor(val) < 1e-6)
      return (SCIP_Longint)floor(val);
   else
      return (SCIP_Longint)ceil(val);
}

SCIP_Real MyPrimalbound(SCIP* scip){
   SCIP_Real ub = SCIPgetPrimalbound(scip);
   if(SCIPisInfinity(scip, fabs(ub)))
      return ub;
   return floor(ub + 1e-6);
}

/** Sets the negative coefficients of slack variables to 0 and returns the sum of the remaining negative coefficients. */
SCIP_Real WeakenSlack(std::vector<SCIP_Real>& vals, const std::vector<SCIP_VAR*>& vars){
   SCIP_Real Negterm=0;
   for(size_t i = 0; i < vals.size(); ++i )
   {
      if(isSlackVar(vars[i]))
         vals[i] = std::max(0.0,vals[i]);
      if(vals[i]<0)
         Negterm+=vals[i];
   }
   assert(vals[vals.size()-1]!=-1);
   return Negterm;
}

/** Sets all negative coefficients to 0. */
void Weaken(std::vector<SCIP_Real>& vals){
   for(size_t i = 0; i < vals.size(); ++i )
   {
      vals[i] = std::max(0.0,vals[i]);
   }
   assert(vals[vals.size()-1]!=-1);
}

/** Caps the coefficients of non-slack variables to rhs. */
void Saturate(std::vector<SCIP_Real>& vals, const std::vector<SCIP_VAR*>& vars, SCIP_Real rhs){
   for(size_t i = 0; i < vals.size(); ++i )
   {
      if(!isSlackVar(vars[i]))
         vals[i] = std::min(vals[i],rhs);
   }
   assert(vals[vals.size()-1]!=-1);
}

/** queues a cut to be added in the next separation round */
SCIP_RETCODE SCIPqueueCut(
   SCIP*                 scip,
   const char*           name,
   int                   nvars,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_Real             lhs,
   SCIP_Real             rhs,
   MyQueue*              q,
   const std::vector<SCIP_Real>* refsol
   )
{
   QueuedCut cut;
   cut.name = (name != NULL) ? name : "queued_cut";
   cut.vars.assign(vars, vars + nvars);
   cut.vals.assign(vals, vals + nvars);
   cut.lhs = lhs;
   cut.rhs = rhs;
   if( refsol != nullptr && !refsol->empty() )
   {
      cut.refvals.reserve(cut.vars.size());
      for( SCIP_VAR* v : cut.vars )
      {
         const int idx = allVarsIndex(v);
         cut.refvals.push_back(idx < (int)refsol->size() ? (*refsol)[idx] : 0.0);
      }
   }
   q->queue.push_back(cut);
   return SCIP_OKAY;
}

/** prints information about queued cuts */
void SCIPprintQueuedCuts(
   SCIP*                 scip,
   MyQueue*              q
   )
{
   int ncuts = (int)q->queue.size();

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "==================================================\n");
   SCIPinfoMessage(scip, NULL, "Queued cuts: %d\n", ncuts);
   SCIPinfoMessage(scip, NULL, "==================================================\n");

   if( ncuts == 0 )
   {
      SCIPinfoMessage(scip, NULL, "  (queue is empty)\n");
      SCIPinfoMessage(scip, NULL, "==================================================\n\n");
      return;
   }

   for( size_t i = 0; i < q->queue.size(); ++i )
   {
      const QueuedCut& cut = q->queue[i];

      SCIPinfoMessage(scip, NULL, "\nCut #%d: %s\n", (int)i + 1, cut.name.c_str());
      SCIPinfoMessage(scip, NULL, "  LHS: %g\n", cut.lhs);
      SCIPinfoMessage(scip, NULL, "  RHS: %g\n", cut.rhs);
      SCIPinfoMessage(scip, NULL, "  Variables (%d):\n", (int)cut.vars.size());

      for( size_t j = 0; j < cut.vars.size(); ++j )
      {
         if( j < 10 || j >= cut.vars.size() - 3 )  // Show the first 10 and the last 3
         {
            SCIPinfoMessage(scip, NULL, "    %+.6g * %s\n", cut.vals[j], SCIPvarGetName(cut.vars[j]));
         }
         else if( j == 10 )
         {
            SCIPinfoMessage(scip, NULL, "    ... (%d more variables)\n",
                          (int)cut.vars.size() - 13);
         }
      }
   }

   SCIPinfoMessage(scip, NULL, "\n==================================================\n\n");
}


/** Option -subslack: replaces in the cut every learned slack s_k by its definition a_k x - b_k, then reduces the
 *  complementary literals and drops the negligible coefficients. Returns false if the cut becomes trivial. */
static bool substituteSlacks(SCIP* scip, QueuedCut& cut, MemoSCIPContext* ctx)
{
   const size_t n = ctx->AllVars.size();
   std::vector<SCIP_Real> coefs(n, 0.0);
   for( size_t i = 0; i < cut.vars.size(); ++i )
      coefs[allVarsIndex(cut.vars[i])] += cut.vals[i];

   SCIP_Real rhs = cut.rhs;
   bool changed = false;
   for( size_t j = 0; j < n; ++j )
   {
      if( coefs[j] == 0.0 || j >= ctx->SlackHasDef.size() || !ctx->SlackHasDef[j] )
         continue;
      // A negative slack term is dropped (s_j >= 0), as the weakening does without the option: substituting it would
      // lower the rhs (possibly below 0) and only add negative coefficients that the weakening removes anyway
      const bool modelSlack = j < ctx->SlackIsModel.size() && ctx->SlackIsModel[j];
      if( coefs[j] < 0.0 && !modelSlack )
      {
         coefs[j] = 0.0;
         changed = true;
         continue;
      }
      // c s_j = c a_j x - c b_j: the definitions only contain non learned-slack variables (a model slack is
      // substituted for both signs of c, its definition being an equality of the model)
      const SCIP_Real c = coefs[j];
      coefs[j] = 0.0;
      for( const auto& t : ctx->SlackDefTerms[j] )
         coefs[t.first] += c * t.second;
      rhs += c * ctx->SlackDefRhs[j];
      changed = true;
   }
   if( !changed )
      return true;
   ctx->nbSubstCuts++;

   for( size_t i = 0; i < n; ++i )
      if( coefs[i] != 0.0 )
         reduceComplementaryCoefs(scip, coefs, i, rhs);
   dropNegligibleCoefs(scip, coefs, rhs);

   cut.vars.clear();
   cut.vals.clear();
   bool allNonneg = true;
   for( size_t i = 0; i < n; ++i )
   {
      if( coefs[i] == 0.0 )
         continue;
      cut.vars.push_back(ctx->AllVars[i]);
      cut.vals.push_back(coefs[i]);
      allNonneg = allNonneg && coefs[i] > 0.0;
   }
   cut.rhs = rhs;
   if( g_options.verbosity > 2 )
   {
      std::cout << "Cut " << cut.name << " after slack substitution:";
      for( size_t i = 0; i < cut.vars.size(); ++i )
         std::cout << " " << (cut.vals[i] > 0 ? "+" : "") << cut.vals[i] << " " << SCIPvarGetName(cut.vars[i]);
      std::cout << " >= " << cut.rhs << std::endl;
   }

   // Trivial: 0 >= rhs, or nonnegative coefficients (on nonnegative variables) with rhs <= 0
   const SCIP_Real tol = 1e-6 * std::max(1.0, fabs(rhs));
   if( cut.vars.empty() || (allNonneg && rhs <= tol) )
   {
      if( cut.vars.empty() && rhs > tol )
         std::cout << "WARNING: cut " << cut.name << " reduces to 0 >= " << rhs << " after slack substitution" << std::endl;
      if( g_options.verbosity > 2 )
         std::cout << "Cut " << cut.name << " trivial after slack substitution (" << cut.vars.size() << " vars, rhs " << rhs << ")" << std::endl;
      ctx->nbSubstTrivial++;
      return false;
   }
   return true;
}

/** Adds the cut consvals . consvars >= rhs with a new slack variable sc_<k>: consvals . consvars - sc_<k> = rhs.
 *  With sepa, it is added as a separator row (and to the cut pool), otherwise as a linear constraint. */
static SCIP_RETCODE addCutWithSlack(SCIP* scip, MemoSCIPContext* g_memoscip, const std::string& cutname,
   std::vector<SCIP_VAR*> consvars, std::vector<SCIP_Real> consvals, SCIP_Real rhs, bool sepa, int* nbctr,
   SCIP_Longint nodenumber, SCIP_SEPA* sepaname)
{
   // Option -subslack: definition of the new slack, s = consvals . consvars - rhs (no learned slack left in it)
   std::vector<std::pair<int, SCIP_Real>> slackDef;
   if( g_options.subSlack )
   {
      slackDef.reserve(consvars.size());
      for( size_t i = 0; i < consvars.size(); ++i )
      {
         const int idx = allVarsIndex(consvars[i]);
         assert(idx >= (int)g_memoscip->SlackHasDef.size() || !g_memoscip->SlackHasDef[idx]);
         slackDef.emplace_back(idx, consvals[i]);
      }
   }

   // New slack variable sc_<k> of the cut
   SCIP_VAR* slack=nullptr;
   char slackname[50];
   snprintf(slackname, 50, "sc_%d", g_memoscip->nbAddedctr);
   SCIP_CALL( SCIPcreateVarBasic(scip, &slack, slackname, 0.0, SCIPinfinity(scip),0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL(SCIPvarSetRemovable(slack,FALSE));
   SCIP_CALL(SCIPvarSetInitial(slack,TRUE));

   SCIP_CALL( SCIPaddVar(scip, slack) );
   consvars.push_back(slack);
   consvals.push_back(-1.0);
   g_memoscip->AllVars.push_back(slack);
   if( g_options.parallel )
      parallelRegister(scip, cutname, slack, consvars, consvals);
   g_memoscip->VarUb.push_back(SCIPinfinity(scip));
   g_memoscip->VarLb.push_back(0);
   const int slackAllIdx = (int)g_memoscip->AllVars.size()-1;
   g_memoscip->SlackImpliedUb.resize(g_memoscip->AllVars.size(), std::numeric_limits<SCIP_Real>::quiet_NaN());
   g_memoscip->SlackCons.resize(g_memoscip->AllVars.size(), nullptr);
   if( g_options.subSlack )
   {
      g_memoscip->SlackHasDef.resize(g_memoscip->AllVars.size(), 0);
      g_memoscip->SlackDefTerms.resize(g_memoscip->AllVars.size());
      g_memoscip->SlackDefRhs.resize(g_memoscip->AllVars.size(), 0.0);
      g_memoscip->SlackHasDef[slackAllIdx] = 1;
      g_memoscip->SlackDefTerms[slackAllIdx] = std::move(slackDef);
      g_memoscip->SlackDefRhs[slackAllIdx] = rhs;
   }

   const int slackIdx = SCIPvarGetIndex(slack)-g_memoscip->nOriginalVar;
   g_memoscip->CurrentVarIdx.resize(slackIdx+1);
   g_memoscip->AuxVar.resize(slackIdx+1);
   g_memoscip->CurrentVarIdx[slackIdx] = g_memoscip->AllVars.size()-1;
   g_memoscip->AuxVar[slackIdx] = { (int)g_memoscip->AllVars.size()-1 };
   g_memoscip->nbAddedctr++;

   // slack's column has no valid LP position until flushed; needed before evaluating cut efficacy against the LP solution
   SCIP_CALL( SCIPflushLP(scip) );

   if(sepa){
      SCIP_ROW* row;
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepaname, cutname.c_str(),
                                     rhs,rhs,
                                     FALSE,   // local
                                     FALSE,   // modifiable
                                     TRUE) ); // removable
      for (size_t i = 0; i < consvars.size(); ++i) {
         SCIP_CALL( SCIPaddVarToRow(scip, row, consvars[i], consvals[i]) );
      }
      // Implied upper bound of the slack from the exact row (SCIPaddVarToRow may transform non-active variables)
      g_memoscip->SlackImpliedUb[slackAllIdx] = slackUbFromRow(scip, row, SCIProwGetRhs(row), slack);
      SCIP_Bool infeasible = FALSE;
      SCIP_CALL( SCIPaddPoolCut(scip, row) );   // persistent, re-injected when violated

      SCIP_CALL( SCIPaddRow(scip, row, FALSE /*forcecut*/, &infeasible) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }else{
      SCIP_CONS* cons = NULL;
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname.c_str(),
                                 (int)consvars.size(),
                                 consvars.data(),
                                 consvals.data(),
                                 rhs,
                                 rhs,
                                 TRUE,    // initial
                                 TRUE,    // separate
                                 TRUE,    // enforce
                                 TRUE,    // check
                                 TRUE,    // propagate
                                 FALSE,   // local
                                 FALSE,   // modifiable
                                 FALSE,   // dynamic
                                 FALSE,   // removable
                                 FALSE) );// stickingatnode
      SCIP_CALL( SCIPaddCons(scip, cons) );
      // The implied upper bound of the slack is computed from the row of the constraint once it is in the LP
      g_memoscip->SlackCons[slackAllIdx] = cons;
      g_memoscip->nodeInfoMap[nodenumber].AddedConstraints.push_back(cons);
      g_memoscip->Unconstructed_cons.push_back(cons);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   (*nbctr)++;

   SCIP_CALL(SCIPreleaseVar(scip, &slack));
   return SCIP_OKAY;
}

/** Efficacy (violation / euclidean norm) of the cut vals . vars >= rhs at the point refvals (aligned with vars). */
static SCIP_Real cutEfficacyAt(const std::vector<SCIP_Real>& vals, SCIP_Real rhs, const std::vector<SCIP_Real>& refvals)
{
   SCIP_Real act = 0.0, norm = 0.0;
   for( size_t i = 0; i < vals.size(); ++i )
   {
      act += vals[i] * (i < refvals.size() ? refvals[i] : 0.0);
      norm += vals[i] * vals[i];
   }
   return norm > 0.0 ? (rhs - act) / std::sqrt(norm) : 0.0;
}

/** Option -mir: MIR cut (SCIPcalcMIR, global bounds only) of the learned cut vals . vars >= rhs.
 *  The node LP solution has at most a few fractional variables, so the c-MIR heuristic of SCIP (whose deltas are the
 *  coefficients of the fractional integer variables) finds nothing. The deltas tried are the distinct absolute
 *  coefficients of the integer variables of the cut (and the best one divided by 2, 4, 8); the cut kept is the most
 *  efficacious one on the reference point refvals (aligned with vars; the current LP solution if empty).
 *  Returns false if no MIR cut was found or if it is numerically unsafe. */
static bool computeMIRCut(SCIP* scip, const std::vector<SCIP_VAR*>& vars, const std::vector<SCIP_Real>& vals,
   SCIP_Real rhs, const std::vector<SCIP_Real>& refvals,
   std::vector<SCIP_VAR*>& mirVars, std::vector<SCIP_Real>& mirVals, SCIP_Real& mirRhs, SCIP_Real& mirEfficacy)
{
   const int nprobvars = SCIPgetNVars(scip);
   SCIP_VAR** probvars = SCIPgetVars(scip);
   std::vector<int> inds;
   std::vector<SCIP_Real> negvals;
   std::vector<SCIP_Real> deltas;
   for( size_t i = 0; i < vars.size(); ++i )
   {
      const int pi = SCIPvarGetProbindex(vars[i]);
      if( pi < 0 )
         return false;
      inds.push_back(pi);
      negvals.push_back(-vals[i]);  // the aggregation row is in <= form: -vals . vars <= -rhs
      if( SCIPvarIsIntegral(vars[i]) && fabs(vals[i]) > 1e-6 )
         deltas.push_back(fabs(vals[i]));
   }
   std::sort(deltas.begin(), deltas.end());
   deltas.erase(std::unique(deltas.begin(), deltas.end(), [](SCIP_Real x, SCIP_Real y){ return fabs(x - y) <= 1e-9 * std::max(1.0, y); }), deltas.end());
   if( deltas.empty() )
      return false;

   SCIP_SOL* sol = nullptr;
   if( !refvals.empty() )
   {
      SCIP_CALL_ABORT( SCIPcreateSol(scip, &sol, nullptr) );
      for( size_t i = 0; i < vars.size(); ++i )
         SCIP_CALL_ABORT( SCIPsetSolVal(scip, sol, vars[i], refvals[i]) );
   }
   SCIP_AGGRROW* aggr = nullptr;
   SCIP_CALL_ABORT( SCIPaggrRowCreate(scip, &aggr) );
   SCIP_CALL_ABORT( SCIPaggrRowAddCustomCons(scip, aggr, inds.data(), negvals.data(), (int)inds.size(), -rhs, 1.0, 0, FALSE) );

   std::vector<SCIP_Real> cutcoefs(nprobvars), bestcoefs;
   std::vector<int> cutinds(nprobvars), bestinds;
   SCIP_Real bestrhs = 0.0;
   SCIP_Real bestdelta = -1.0;
   mirEfficacy = -SCIPinfinity(scip);
   auto tryDelta = [&](SCIP_Real delta) {
      SCIP_Real cutrhs = 0.0, efficacy = 0.0;
      int cutnnz = 0, cutrank = 0;
      SCIP_Bool cutislocal = FALSE, success = FALSE;
      SCIP_CALL_ABORT( SCIPcalcMIR(scip, sol, TRUE, 0.9999, 0, FALSE, FALSE, nullptr, nullptr, 0.05, 0.999, 1.0 / delta,
         aggr, cutcoefs.data(), &cutrhs, cutinds.data(), &cutnnz, &efficacy, &cutrank, &cutislocal, &success) );
      if( !success || cutislocal || cutnnz == 0 || efficacy <= mirEfficacy )
         return;
      SCIP_Real maxabs = fabs(cutrhs);
      for( int k = 0; k < cutnnz; ++k )
         maxabs = std::max(maxabs, fabs(cutcoefs[k]));
      if( maxabs >= 1e6 )
         return;
      mirEfficacy = efficacy;
      bestdelta = delta;
      bestrhs = cutrhs;
      bestcoefs.assign(cutcoefs.begin(), cutcoefs.begin() + cutnnz);
      bestinds.assign(cutinds.begin(), cutinds.begin() + cutnnz);
   };
   for( SCIP_Real delta : deltas )
      tryDelta(delta);
   if( bestdelta > 0.0 )
   {
      const SCIP_Real d = bestdelta;
      for( SCIP_Real div : {2.0, 4.0, 8.0} )
         tryDelta(d / div);
   }
   SCIPaggrRowFree(scip, &aggr);
   if( sol != nullptr )
      SCIP_CALL_ABORT( SCIPfreeSol(scip, &sol) );
   if( bestdelta <= 0.0 )
      return false;

   mirVars.clear();
   mirVals.clear();
   for( size_t k = 0; k < bestinds.size(); ++k )
   {
      mirVars.push_back(probvars[bestinds[k]]);
      mirVals.push_back(-bestcoefs[k]);  // back to >= form
   }
   mirRhs = -bestrhs;
   return true;
}

/** True if the two sparse cuts are equal up to a positive scaling (relative tolerance 1e-9). */
static bool sameCutUpToScaling(const std::vector<SCIP_VAR*>& v1, const std::vector<SCIP_Real>& c1, SCIP_Real r1,
   const std::vector<SCIP_VAR*>& v2, const std::vector<SCIP_Real>& c2, SCIP_Real r2)
{
   if( v1.size() != v2.size() || v1.empty() )
      return false;
   std::vector<std::pair<SCIP_VAR*, SCIP_Real>> a, b;
   for( size_t i = 0; i < v1.size(); ++i ) a.emplace_back(v1[i], c1[i]);
   for( size_t i = 0; i < v2.size(); ++i ) b.emplace_back(v2[i], c2[i]);
   std::sort(a.begin(), a.end());
   std::sort(b.begin(), b.end());
   if( a[0].second == 0.0 )
      return false;
   const SCIP_Real ratio = b[0].second / a[0].second;
   if( ratio <= 0.0 )
      return false;
   for( size_t i = 0; i < a.size(); ++i )
      if( a[i].first != b[i].first || fabs(b[i].second - ratio * a[i].second) > 1e-9 * std::max(1.0, fabs(b[i].second)) )
         return false;
   return fabs(r2 - ratio * r1) <= 1e-9 * std::max(1.0, fabs(r2));
}

/** Adds the queued cuts to the problem, each one with a new slack variable: sum coef*x - slack = rhs.
 *  With sepa, cuts are added as separator rows (and to the cut pool), otherwise as linear constraints. */
SCIP_RETCODE AddCutstoProb(SCIP* scip, MyQueue* q, MemoSCIPContext* g_memoscip, double Top,bool sepa, int* nbctr, SCIP_Longint nodenumber, SCIP_SEPA* sepaname)
{
   (void)Top;
   assert(!SCIPinDive(scip));
   while( !q->queue.empty() )
   {
      QueuedCut cut = std::move(q->queue.back());
      q->queue.pop_back();

      assert(!cut.vars.empty());

      // Substitution on the added cut only: the node proofs passed to the parents keep the slacks
      if( (g_options.subSlack || g_options.subModelSlack) && !substituteSlacks(scip, cut, g_memoscip) )
         continue;

      SCIP_Real rhs = roundToDecimals(cut.rhs, 9);

      std::vector<SCIP_VAR*> consvars;
      std::vector<SCIP_Real> consvals;
      consvars.reserve(cut.vars.size()+1);
      consvals.reserve(cut.vals.size()+1);
      SCIP_Real Negterm =0;
      if(g_options.Weaken==0)
         Negterm=WeakenSlack(cut.vals,cut.vars);
      else
         Weaken(cut.vals);

      // rhs <= Negterm (the minimum activity): the cut is trivially satisfied, and saturating at a cap <= 0 would be invalid
      if( rhs - Negterm <= 1e-6 * std::max(1.0, fabs(rhs)) )
      {
         if( g_options.subSlack )
            g_memoscip->nbSubstTrivial++;
         continue;
      }
      Saturate(cut.vals, cut.vars, rhs-Negterm);
      for( size_t i = 0; i < cut.vars.size(); ++i )
      {
         SCIP_Real coef = roundToDecimals(cut.vals[i] , 9);
         if( SCIPisZero(scip, coef) )
            continue;
         consvars.push_back(cut.vars[i]);
         consvals.push_back(coef);
         assert(fabs(coef)<1000000 || coef >= 1e+08);
      }

      if(consvars.empty())
         continue;

      // Option -parallel: a cut nearly parallel to an earlier learned cut is not added
      if( g_options.parallel && parallelCheck(scip, cut.name, consvars, consvals) )
         continue;

      // Option -mir: the c-MIR cut of the learned cut is added next to it (computed before its slack is created)
      std::vector<SCIP_VAR*> mirVars;
      std::vector<SCIP_Real> mirVals;
      SCIP_Real mirRhs = 0.0;
      bool addMir = false;
      if( g_options.mir )
      {
         g_memoscip->nbMirTried++;
         // Reference values of the variables kept in the cut
         std::vector<SCIP_Real> refvals;
         if( cut.refvals.size() == cut.vars.size() )
         {
            for( SCIP_VAR* v : consvars )
            {
               const size_t k = std::find(cut.vars.begin(), cut.vars.end(), v) - cut.vars.begin();
               refvals.push_back(k < cut.vars.size() ? cut.refvals[k] : 0.0);
            }
         }
         SCIP_Real mirEff = 0.0;
         if( computeMIRCut(scip, consvars, consvals, rhs, refvals, mirVars, mirVals, mirRhs, mirEff) )
         {
            const SCIP_Real cutEff = refvals.empty() ? -SCIPinfinity(scip) : cutEfficacyAt(consvals, rhs, refvals);
            if( !refvals.empty() && cutEff > 1e-6 )
               g_memoscip->nbMirCutViolated++;
            if( !refvals.empty() && mirEff > 1e-6 )
               g_memoscip->nbMirViolated++;
            if( g_options.verbosity > 2 )
               std::cout << "MIR of " << cut.name << ": efficacy on the node LP solution " << mirEff << " versus cut " << cutEff << std::endl;
            if( sameCutUpToScaling(consvars, consvals, rhs, mirVars, mirVals, mirRhs) )
               g_memoscip->nbMirSame++;
            else if( mirEff <= cutEff + 1e-9 )
               g_memoscip->nbMirWeaker++;  // not stronger than the learned cut on the reference point
            else
               addMir = true;
         }
         else
            g_memoscip->nbMirFailed++;
      }

      // With -mir=2 a stronger MIR cut replaces the learned cut, otherwise the learned cut is always added (first)
      if( g_options.mir != 2 )
         SCIP_CALL( addCutWithSlack(scip, g_memoscip, cut.name, consvars, consvals, rhs, sepa, nbctr, nodenumber, sepaname) );
      bool mirAdded = false;
      if( addMir )
      {
         // Same filtering as the learned cut: rounded coefficients, zeros removed
         std::vector<SCIP_VAR*> mv;
         std::vector<SCIP_Real> mc;
         for( size_t i = 0; i < mirVars.size(); ++i )
         {
            const SCIP_Real coef = roundToDecimals(mirVals[i], 9);
            if( SCIPisZero(scip, coef) )
               continue;
            mv.push_back(mirVars[i]);
            mc.push_back(coef);
         }
         // The rhs is rounded down so that rounding never makes the cut stronger than the computed one
         const SCIP_Real mrhs = std::floor(mirRhs * 1e9) / 1e9;
         if( !mv.empty() )
         {
            if( g_options.verbosity > 2 )
            {
               std::cout << "MIR cut of " << cut.name << ":";
               for( size_t i = 0; i < mv.size(); ++i )
                  std::cout << " " << mc[i] << " " << SCIPvarGetName(mv[i]);
               std::cout << " >= " << mrhs << std::endl;
            }
            SCIP_CALL( addCutWithSlack(scip, g_memoscip, "cut_mir_" + cut.name, mv, mc, mrhs, sepa, nbctr, nodenumber, sepaname) );
            g_memoscip->nbMirAdded++;
            mirAdded = true;
         }
      }
      if( g_options.mir == 2 && !mirAdded )
         SCIP_CALL( addCutWithSlack(scip, g_memoscip, cut.name, consvars, consvals, rhs, sepa, nbctr, nodenumber, sepaname) );
   }

   return SCIP_OKAY;
}
