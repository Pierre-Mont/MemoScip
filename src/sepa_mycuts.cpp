#include "sepa_mycuts.h"

#include "memoscip_option.h"
#include "memoscip_types.h"

#include "scip/cons_linear.h"
#include "memoscip_dual.h"
#include "memoscip_events.h"
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
   MyQueue*              q
   )
{
   QueuedCut cut;
   cut.name = (name != NULL) ? name : "queued_cut";
   cut.vars.assign(vars, vars + nvars);
   cut.vals.assign(vals, vals + nvars);
   cut.lhs = lhs;
   cut.rhs = rhs;
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
      if( coefs[j] < 0.0 )
      {
         coefs[j] = 0.0;
         changed = true;
         continue;
      }
      // c s_j = c a_j x - c b_j: the definitions only contain non learned-slack variables
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
      if( g_options.subSlack && !substituteSlacks(scip, cut, g_memoscip) )
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
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &row, sepaname, cut.name.c_str(),
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
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cut.name.c_str(),
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
   }

   return SCIP_OKAY;
}
