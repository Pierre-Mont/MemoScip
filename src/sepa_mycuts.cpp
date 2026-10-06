#include "memoscip_types.h"

#include "sepa_mycuts.h"

#include "scip/cons_linear.h"
#include "scip/var.h"
#include "memoscip_option.h"

#include <algorithm>
#include <cmath>
#include <iostream>
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

double roundToDecimals(double val, int decimals) {
    double factor = std::pow(10.0, decimals);
    return std::round(val * factor) / factor;
}


SCIP_Longint Myceil(double val){
      if(val - floor(val) < 1e-5)
         return (SCIP_Longint)floor(val);
      else
         return (SCIP_Longint)ceil(val);
}

SCIP_Real WeakenSlack(std::vector<SCIP_Real>& vals, const std::vector<SCIP_VAR*>& vars){
      SCIP_Real Negterm=0;
      for(size_t i = 0; i < vals.size(); ++i )
      {
         if(SCIPvarGetName(vars[i])[2]=='s' || SCIPvarGetName(vars[i])[0]=='s')
            vals[i] = std::max(0.0,vals[i]);
         if(vals[i]<0)
            Negterm+=vals[i];
      }
      assert(vals[vals.size()-1]!=-1);
      return Negterm;
}


void Weaken(std::vector<SCIP_Real>& vals){
      for(size_t i = 0; i < vals.size(); ++i )
      {
         vals[i] = std::max(0.0,vals[i]);
      }
      assert(vals[vals.size()-1]!=-1);
}

void Saturate(std::vector<SCIP_Real>& vals, const std::vector<SCIP_VAR*>& vars, SCIP_Real rhs){
      for(size_t i = 0; i < vals.size(); ++i )
      {
         if(SCIPvarGetName(vars[i])[2]!='s' && SCIPvarGetName(vars[i])[0]!='s')
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
   MyQueue*                q
   )
{   
   // Créer le cut en attente
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
   MyQueue*               q
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
         const char* varname = SCIPvarGetName(cut.vars[j]);
         SCIP_Real coef = cut.vals[j];
         
         if( j < 10 || j >= cut.vars.size() - 3 )  // Affiche les 10 premiers et 3 derniers
         {
            SCIPinfoMessage(scip, NULL, "    %+.6g * %s\n", coef, varname);
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



SCIP_RETCODE AddCutstoProb(SCIP* scip, MyQueue* q, MemoSCIPContext* g_memoscip, double Top,bool sepa, int* nbctr, SCIP_Longint nodenumber, SCIP_SEPA* sepaname)
{
   (void)Top;
   assert(!SCIPinDive(scip));
   while( !q->queue.empty() )
   {
      //std::cout<<"Adding cut "<<q->queue.back().name<<" with rhs "<<q->queue.back().rhs<<" and "<<q->queue.back().vars.size()<<" variables"<<std::endl;
      QueuedCut cut = std::move(q->queue.back());
      q->queue.pop_back();
      
      // Créer la contrainte linéaire
      
      //std::cout << "Adding cut " << cut.name << " with rhs " << cut.rhs << " and " << cut.vars.size() << " variables" <<cut.vals.size() <<std::endl;
      //std::cout<< " d "<< cut.vals.back() <<std::endl;
      assert(!cut.vars.empty());
      //std::cout<<SCIPvarGetName(cut.vars.back())<<std::endl;
      //assert(SCIPvarGetName(cut.vars.back())[0]=='s');
      //assert(cut.vals.back()!=-1);
      /*SCIP_Real lhs = ceil(cut.lhs);
      SCIP_Real rhs = ceil(cut.rhs);*/

      SCIP_Real rhs = roundToDecimals(cut.rhs, 9);
      // Ajouter les variables

      std::vector<SCIP_VAR*> consvars;
      std::vector<SCIP_Real> consvals;
      consvars.reserve(cut.vars.size()+1);
      consvals.reserve(cut.vals.size()+1);
      SCIP_Real Negterm =0;
      if(g_options.Weaken==0)
         Negterm=WeakenSlack(cut.vals,cut.vars);
      else
         Weaken(cut.vals);
      
      Saturate(cut.vals, cut.vars, rhs-Negterm);
      for( size_t i = 0; i < cut.vars.size(); ++i )
      {
         SCIP_Real coef = roundToDecimals(cut.vals[i] , 9);
         if( SCIPisZero(scip, coef) )
            continue;
         consvars.push_back(cut.vars[i]);
         consvals.push_back(coef);
         assert(fabs(coef)<1000000);
      }

      if(consvars.empty())
         continue;
      //SCIP_VAR* slack = g_memoscip->slackVar[0];
      //g_memoscip->slackVar.erase(g_memoscip->slackVar.begin());
      SCIP_VAR* slack=nullptr;
      char slackname[50];
      snprintf(slackname, 50, "sc_%d", g_memoscip->nbAddedctr);
      SCIP_CALL( SCIPcreateVarBasic(scip, &slack, slackname, 0.0, SCIPinfinity(scip),0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL(SCIPvarSetRemovable(slack,FALSE));
      SCIP_CALL(SCIPvarSetInitial(slack,TRUE));

      //SCIP_CALL( SCIPaddPricedVar(scip, slack, 100) );
      SCIP_CALL( SCIPaddVar(scip, slack) );
      consvars.push_back(slack);
      consvals.push_back(-1.0);
      g_memoscip->AllVars.push_back(slack);
      g_memoscip->VarUb.push_back(SCIPinfinity(scip));
      g_memoscip->VarLb.push_back(0);

      g_memoscip->CurrentVarIdx.resize(SCIPvarGetIndex(slack)+1-g_memoscip->nOriginalVar);
      g_memoscip->AuxVar.resize(SCIPvarGetIndex(slack)+1-g_memoscip->nOriginalVar);
      g_memoscip->TupleAuxVar.resize(SCIPvarGetIndex(slack)+1-g_memoscip->nOriginalVar);
      g_memoscip->CurrentVarIdx[SCIPvarGetIndex(slack)-g_memoscip->nOriginalVar] = g_memoscip->AllVars.size()-1;
      g_memoscip->AuxVar[SCIPvarGetIndex(slack)-g_memoscip->nOriginalVar] = g_memoscip->AllVars.size()-1;
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
            SCIP_Bool infeasible = FALSE;
            SCIP_CALL( SCIPaddPoolCut(scip, row) );   // persistant, réinjecté quand violé

            SCIP_CALL( SCIPaddRow(scip, row, FALSE /*forcecut*/, &infeasible) );
         //SCIP_Real efficacy = SCIPgetCutEfficacy(scip, nullptr, row);
         //printf("Cut efficacy: %f\n", efficacy);
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }else{
         SCIP_CONS* cons = NULL;
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cut.name.c_str(),
                                    (int)consvars.size(),
                                    consvars.empty() ? NULL : consvars.data(),
                                    consvals.empty() ? NULL : consvals.data(),
                                    rhs,
                                    rhs,
                                    TRUE,    // initial
                                    TRUE,   // separate
                                    TRUE,    // enforce
                                    TRUE,    // check
                                    TRUE,    // propagate
                                    FALSE,   // local
                                    FALSE,   // modifiable
                                    FALSE,   // dynamic
                                    FALSE,   // removable
                                    FALSE) );// stickingatnode
         SCIP_CALL( SCIPaddCons(scip, cons) );
         g_memoscip->nodeInfoMap[nodenumber].AddedConstraints.push_back(cons);
         g_memoscip->Unconstructed_cons.push_back(cons);
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      }
      //assert(SCIPvarIsInLP(consvars.back()));
      // Ajouter la contrainte au problème
      //std::cout<<"ADD CONS "<<SCIPconsGetName(cons)<<" with rhs "<<rhs<<" and "<<consvars.size()<<" variables"<<std::endl;
      


      /*std::cout<<"SEPA MYCUTS Adding cut "<<cut.name<<" with rhs "<<rhs<<" and "<<consvars.size()<<" variables"<<std::endl;
      //Print the constraint 
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      std::cout << "active = " << SCIPconsIsActive(cons) << std::endl;
      std::cout << "deleted = " << SCIPconsIsDeleted(cons) << std::endl; 
      std::cout << "enabled = " << SCIPconsIsEnabled(cons) << std::endl;
      std::cout<<std::endl<<"Cut added to the problem"<<std::endl; */

      (*nbctr)++;
                 
      SCIP_CALL(SCIPreleaseVar(scip, &slack));
   }
   
   return SCIP_OKAY;
}