#ifndef MEMOSCIP_DUAL_H
#define MEMOSCIP_DUAL_H

#include "scip/scip.h"
#include "memoscip_types.h"

#include <vector>

std::vector<SCIP_Real> FusionResolution(SCIP* scip, SCIP_Longint child0, SCIP_Longint child1);
SCIP_RETCODE FinishLearn(SCIP* scip, SCIP_Longint node);
SCIP_RETCODE AnalyseDualSol(SCIP* scip, SCIP_Longint nodeNumber);

SCIP_RETCODE gatherLPInfeasibilityCertificate(
    SCIP* scip,
    SCIP_Longint nodeNumber);

SCIP_RETCODE getDualProofFromPseudoObjCutoff(
    SCIP* scip,
    SCIP_Longint nodeNumber, SCIP_Longint parentNumber, int decVar);

/** Writes the LP and the transformed problem of the node when option -o=N is given with N > 1 and node >= N. */
SCIP_RETCODE writeNodeLPFilesIfRequested(SCIP* scip, SCIP_Longint nodenumber);

/* Steps shared by FinishLearn and the node-deletion event to build the proof of a node from the proofs of its children. */

/** Debug checks on the proof of the only child of a node. */
void checkSingleChildProof(SCIP_Longint nodenumber, SCIP_Longint child);

/** If both children increased the lower bound, fuses their proofs into ResFus with right-hand side rhsFus. */
SCIP_RETCODE fuseChildrenProofs(SCIP* scip, SCIP_Longint nodenumber, const std::vector<long>& children,
    std::vector<SCIP_Real>& ResFus, SCIP_Real& rhsFus);

/** Adds the node's own dual proof (if it increased the lower bound) to ResFus and rhs_summed. */
void addNodeProof(SCIP_Longint nodenumber, std::vector<SCIP_Real>& ResFus, SCIP_Real& rhs_summed);

/** Replaces variables removed by bound propagation with their explanation. */
void applyBoundPropagationRemovals(SCIP_Longint nodenumber, std::vector<SCIP_Real>& ResFus, SCIP_Real rhs_summed);

/** Replaces variables removed by reduced cost filtering (node consistency) with the parent reduced costs. */
void applyNodeConsistency(SCIP* scip, SCIP_Longint nodenumber, SCIP_Longint parentNumber,
    std::vector<SCIP_Real>& ResFus, SCIP_Real& rhs_summed);

/** Debug check (option -d) that rhs_summed matches the lower bound increases of the node and its children. */
void checkSummedRhs(SCIP_Longint nodenumber, const std::vector<long>& children, SCIP_Real rhs_summed);

/** Removes the children of nodenumber from nodeInfoMap. */
void eraseChildren(SCIP_Longint nodenumber, const std::vector<long>& children);

/** Prints the nonzero coefficients of coefs followed by "<relation> rhs".
 *  If mapThroughCurrentVarIdx, variable i is printed as AllVars[CurrentVarIdx[i]] instead of AllVars[i]. */
void printProof(SCIP* scip, const std::vector<SCIP_Real>& coefs, const char* relation, SCIP_Real rhs,
    bool mapThroughCurrentVarIdx);

#endif
