#ifndef MEMOSCIP_EVENTS_H
#define MEMOSCIP_EVENTS_H

#include "memoscip_option.h"

#include "scip/scip.h"

#include <vector>

void includeNodeSolvedEventHandler(SCIP* scip);
/** Removes from coefs the common part of the coefficient of literal i and of its complementary literal(s), decreasing rhs. */
void reduceComplementaryCoefs(SCIP* scip, std::vector<SCIP_Real>& coefs, size_t i, SCIP_Real& rhs);
void includeBoundChangedEventHandler(SCIP* scip);
void includeNodeFocusedEventHandler(SCIP* scip);
void includeNodeInfeasibleEventHandler(SCIP* scip);
void includeLPSolveEventHandler(SCIP* scip);
void includeBestSolFoundEventHandler(SCIP* scip);
void includeNodeDeleteEventHandler(SCIP* scip);
void register_memoscip_branchrule(SCIP* scip);

SCIP_BDCHGINFO* getBoundChangeReasonDetailed(SCIP* scip, SCIP_VAR* var,
    SCIP_BOUNDTYPE boundtype,
    SCIP_Real oldbound,
    SCIP_Real newbound);

#endif
