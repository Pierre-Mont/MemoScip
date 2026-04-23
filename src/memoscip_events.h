#ifndef MEMOSCIP_EVENTS_H
#define MEMOSCIP_EVENTS_H

#include "memoscip_option.h"

#include "scip/scip.h"

void includeNodeSolvedEventHandler(SCIP* scip);
void includeBoundChangedEventHandler(SCIP* scip);
void includeNodeFocusedEventHandler(SCIP* scip);
void includeNodeInfeasibleEventHandler(SCIP* scip);
void includeLPSolveEventHandler(SCIP* scip);
void includeNodeDeleteEventHandler(SCIP* scip);
void register_memoscip_branchrule(SCIP* scip);

SCIP_BDCHGINFO* getBoundChangeReasonDetailed(SCIP* scip, SCIP_VAR* var,
    SCIP_BOUNDTYPE boundtype,
    SCIP_Real oldbound,
    SCIP_Real newbound);

#endif
