#ifndef SEPA_MYCUTS_H
#define SEPA_MYCUTS_H

#include "scip/scip.h"

#include <string>
#include <vector>

struct QueuedCut
{
   std::string name;
   std::vector<SCIP_VAR*> vars;
   std::vector<SCIP_Real> vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
};

struct MyQueue
{
   std::vector<QueuedCut> queue;
};

void SCIPprintQueuedCuts(SCIP* scip, MyQueue* q);
void Weaken(std::vector<SCIP_Real>& vals);
void Saturate(std::vector<SCIP_Real>& vals, double rhs);
/** queues a cut to be added in the next separation round */
SCIP_EXPORT
SCIP_RETCODE SCIPqueueCut(
   SCIP* scip,
   const char* name,
   int nvars,
   SCIP_VAR** vars,
   SCIP_Real* vals,
   SCIP_Real lhs,
   SCIP_Real rhs,
   MyQueue* q);

SCIP_RETCODE AddCutstoProb(SCIP* scip, MyQueue* q, double Top);

#endif
