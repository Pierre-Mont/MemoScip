#ifndef MEMOSCIP_DUAL_H
#define MEMOSCIP_DUAL_H

#include "scip/scip.h"
#include "memoscip_types.h"

#include <vector>

void Weaken(std::vector<SCIP_Real>& coef);
void Saturate(std::vector<SCIP_Real>& coef, SCIP_Real rhs);
std::vector<SCIP_Real> FusionResolution(SCIP* scip, SCIP_Longint child0, SCIP_Longint child1);
SCIP_RETCODE FinishLearn(SCIP* scip, SCIP_Longint node);
SCIP_RETCODE AnalyseDualSol(SCIP* scip, SCIP_Longint nodeNumber);


SCIP_RETCODE gatherLPInfeasibilityCertificate(
    SCIP* scip,
    SCIP_Longint nodeNumber);

SCIP_RETCODE gatherFarkasCertificate(
    SCIP* scip,
    SCIP_Longint nodeNumber);
SCIP_RETCODE getDualProofFromPseudoObjCutoff(
    SCIP* scip,
    SCIP_Longint nodeNumber, SCIP_Longint parentNumber, int decVar);
#endif
