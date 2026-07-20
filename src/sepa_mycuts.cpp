#include "sepa_mycuts.h"

#include "scip/cons_linear.h"
#include "scip/var.h"

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

double roundToDecimals(double val, int decimals) {
    double factor = std::pow(10.0, decimals);
    return std::ceil(val * factor) / factor;
}


SCIP_Longint Myceil(double val){
      if(val - floor(val) < 1e-7)
         return (SCIP_Longint)floor(val);
      else
         return (SCIP_Longint)ceil(val);
}

void Weaken(std::vector<SCIP_Real>& vals){
      for(size_t i = 0; i < vals.size()-1; ++i )
      {
         vals[i] = std::max(0.0,vals[i]);
      }
      assert(vals[vals.size()-1]==-1);
}

void Saturate(std::vector<SCIP_Real>& vals, const std::vector<SCIP_VAR*>& vars, double rhs){
      for(size_t i = 0; i < vals.size()-1; ++i )
      {
         if(SCIPvarGetName(vars[i])[2]!='s')
            vals[i] = std::min(vals[i],rhs);
      }
      assert(vals[vals.size()-1]==-1);
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


/*

SCIP_RETCODE AddCutstoProb(SCIP* scip, MyQueue* q, double Top)
{
   (void)Top;
   
   while( !q->queue.empty() )
   {
      //std::cout<<"Adding cut "<<q->queue.back().name<<" with rhs "<<q->queue.back().rhs<<" and "<<q->queue.back().vars.size()<<" variables"<<std::endl;
      QueuedCut cut = std::move(q->queue.back());
      q->queue.pop_back();
   // Build and add a row directly to the LP (lighter than creating SCIP_CONS).
   SCIP_ROW* row = NULL;

      SCIP_Real lhs = ceil(cut.lhs);
      SCIP_Real rhs = ceil(cut.rhs);
les variables

      Weaken(cut.vals);
      //Saturate(cut.vals, roundToDecimals(cut.lhs, 5));
      Saturate(cut.vals, lhs);

      std::vector<SCIP_VAR*> consvars;
      std::vector<SCIP_Real> consvals;
      consvars.reserve(cut.vars.size());
      consvals.reserve(cut.vals.size());

      for( size_t i = 0; i < cut.vars.size(); ++i )
      {
         SCIP_Real coef = ceil(cut.vals[i]);
         if( SCIPisZero(scip, coef) )
            continue;
         consvars.push_back(cut.vars[i]);
         consvals.push_back(coef);
      }

      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, cut.name.c_str(),
                                          lhs, rhs,
                                          FALSE,  // local
                                          FALSE,  // modifiable
                                          FALSE) ); // removable

      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      for( size_t i = 0; i < consvars.size(); ++i )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, consvars[i], consvals[i]) );
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      SCIP_Bool infeasible = FALSE;
      SCIP_CALL( SCIPaddRow(scip, row, TRUE, &infeasible) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );   
   }
   
   return SCIP_OKAY;
}

*/



SCIP_RETCODE AddCutstoProb(SCIP* scip, MyQueue* q, double Top,bool sepa, int* nbctr)
{
   (void)Top;
   
   while( !q->queue.empty() )
   {
      //std::cout<<"Adding cut "<<q->queue.back().name<<" with rhs "<<q->queue.back().rhs<<" and "<<q->queue.back().vars.size()<<" variables"<<std::endl;
      QueuedCut cut = std::move(q->queue.back());
      q->queue.pop_back();
      // Créer la contrainte linéaire
      SCIP_CONS* cons = NULL;

      /*SCIP_Real lhs = ceil(cut.lhs);
      SCIP_Real rhs = ceil(cut.rhs);*/

      SCIP_Real rhs = roundToDecimals(cut.rhs, 5);
      // Ajouter les variables

      Weaken(cut.vals);
      Saturate(cut.vals, cut.vars, rhs);
      //Saturate(cut.vals, lhs);

      std::vector<SCIP_VAR*> consvars;
      std::vector<SCIP_Real> consvals;
      consvars.reserve(cut.vars.size());
      consvals.reserve(cut.vals.size());

      for( size_t i = 0; i < cut.vars.size(); ++i )
      {
         SCIP_Real coef = roundToDecimals(cut.vals[i] , 5);
         if( SCIPisZero(scip, coef) )
            continue;
         consvars.push_back(cut.vars[i]);
         consvals.push_back(coef);
      }
      if(sepa){      
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cut.name.c_str(),
                                    (int)consvars.size(),
                                    consvars.empty() ? NULL : consvars.data(),
                                    consvals.empty() ? NULL : consvals.data(),
                                    rhs,
                                    rhs,
                                    FALSE,    // initial
                                    FALSE,   // separate
                                    FALSE,    // enforce
                                    FALSE,    // check
                                    TRUE,    // propagate
                                    FALSE,   // local
                                    FALSE,   // modifiable
                                    TRUE,   // dynamic
                                    TRUE,   // removable
                                    FALSE) );// stickingatnode

         SCIP_ROW* row;
         SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, cut.name.c_str(),
                                             rhs, rhs, FALSE, FALSE, TRUE) ); // removable=TRUE
         for (size_t i = 0; i < consvars.size(); ++i) {
            SCIP_CALL( SCIPaddVarToRow(scip, row, consvars[i], consvals[i]) );
         }
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );

      }else{
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cut.name.c_str(),
                                    (int)consvars.size(),
                                    consvars.empty() ? NULL : consvars.data(),
                                    consvals.empty() ? NULL : consvals.data(),
                                    rhs,
                                    rhs,
                                    TRUE,    // initial
                                    FALSE,   // separate
                                    TRUE,    // enforce
                                    FALSE,    // check
                                    TRUE,    // propagate
                                    FALSE,   // local
                                    FALSE,   // modifiable
                                    TRUE,   // dynamic
                                    TRUE,   // removable
                                    FALSE) );// stickingatnode
      }
      
      // Ajouter la contrainte au problème
      SCIP_CALL( SCIPaddCons(scip, cons) );

      // Libérer la contrainte (SCIP en garde une référence)
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      (*nbctr)++;
   
      //std::cout<<"Adding cut "<<cut.name<<" with rhs "<<rhs<<" and "<<consvars.size()<<" variables"<<std::endl;
                 
      
   }
   
   return SCIP_OKAY;
}