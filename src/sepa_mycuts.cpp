#include "sepa_mycuts.h"

#include "scip/cons_linear.h"
#include "scip/var.h"

#include <algorithm>
#include <cmath>

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
void Weaken(std::vector<SCIP_Real>& vals){
      for(size_t i = 0; i < vals.size()-1; ++i )
      {
         vals[i] = std::max(0.0,vals[i]);
      }
      assert(vals[vals.size()-1]==-1);
}

void Saturate(std::vector<SCIP_Real>& vals, double rhs){
      for(size_t i = 0; i < vals.size()-1; ++i )
      {
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




SCIP_RETCODE AddCutstoProb(SCIP* scip, MyQueue* q, double Top)
{
   //SCIPprintQueuedCuts(scip,q);
     
   //SCIPinfoMessage(scip, NULL, "Forcing %d queued cuts into LP...\n", nqueued);
   
   // Forcer tous les cuts dans le LP
   while( !q->queue.empty() )
   {
      QueuedCut cut = q->queue.back();
      q->queue.pop_back();
      /*SCIP_ROW* row = NULL;
      
      // Créer la row
      SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, cut.name.c_str(),
                                        ceil(cut.lhs), ceil(cut.rhs),
                                        FALSE,   // local
                                        TRUE,   // modifiable
                                        FALSE) ); // removable
      
      // Ajouter les variables
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      for( size_t i = 0; i < cut.vars.size(); ++i )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, cut.vars[i], ceil(cut.vals[i]) ));
      }
      
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      
      // Ajouter au LP
      SCIP_Bool infeasible = FALSE;
      //SCIP_CALL( SCIPaddRowProbing(scip, row) );
      SCIP_CALL( SCIPForceaddRow(scip, row) );

      // Ajouter au cut pool global (optionnel)
      //SCIP_CALL( SCIPaddPoolCut(scip, row) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      if( infeasible )
      {
         return SCIP_OKAY;
      }*/

      
      
      // Créer la contrainte linéaire
      SCIP_CONS* cons = NULL;

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cut.name.c_str(),
                                    0,       // nvars 
                                    NULL,    // vars
                                    NULL,    // vals
                                    roundToDecimals(cut.lhs, 5), // lhs
                                    roundToDecimals(cut.rhs, 5), // rhs
                                    TRUE,    // initial
                                    FALSE,   // separate
                                    TRUE,    // enforce
                                    TRUE,    // check
                                    TRUE,    // propagate
                                    FALSE,   // local
                                    FALSE,    // modifiable
                                    FALSE,   // dynamic
                                    FALSE,   // removable
                                    FALSE) );// stickingatnode

      /*SCIP_Real NegCoef=0;
      for( size_t i = 0; i < cut.vars.size()-1; ++i )
      {
         if(cut.vals[i]<0 && SCIPvarGetName(cut.vars[i])[2]!='s')
            NegCoef+=cut.vals[i];
         if(SCIPvarGetName(cut.vars[i])[2]=='s')
            NegCoef+=-9000;
      }*/
      // Ajouter les variables

      Weaken(cut.vals);
      Saturate(cut.vals, roundToDecimals(cut.lhs, 5));
      for( size_t i = 0; i < cut.vars.size(); ++i )
      {
         if(SCIPvarGetName(cut.vars[i])[2]!='s'){
            //SCIP_CALL( SCIPaddCoefLinear(scip, cons, cut.vars[i], std::min(cut.lhs - NegCoef, cut.vals[i])) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, cut.vars[i],  roundToDecimals(cut.vals[i], 5)));
         }
         else
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, cut.vars[i], roundToDecimals(cut.vals[i], 5)) );
      }
      

      // Ajouter la contrainte au problème
      SCIP_CALL( SCIPaddCons(scip, cons) );

      // Libérer la contrainte (SCIP en garde une référence)
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      
                 
      
   }
   
   return SCIP_OKAY;
}