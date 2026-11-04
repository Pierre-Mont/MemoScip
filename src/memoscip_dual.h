#ifndef MEMOSCIP_DUAL_H
#define MEMOSCIP_DUAL_H

#include "scip/scip.h"
#include "memoscip_types.h"

#include <vector>

/** Sets to 0 the coefficients of the proof sum coefs[i] x_i >= rhs that are below 1e-9 times the largest one.
 *  rhs is relaxed with the global bound of each dropped variable, so the proof stays valid. */
void dropNegligibleCoefs(SCIP* scip, std::vector<SCIP_Real>& coefs, SCIP_Real& rhs);

/** Upper bound of the slack s of the cut row a x - s (+ row constant) = rhs, from the exact LP row: max activity of
 *  a x + constant - rhs over the global bounds (the bounds of other cut slacks are their implied ones).
 *  rhs is passed explicitly since the sides of a row are changed in dives. +infinity if unbounded. */
SCIP_Real slackUbFromRow(SCIP* scip, SCIP_ROW* row, SCIP_Real rhs, SCIP_VAR* slack);

/** Implied upper bound of AllVars[idx] if it is the slack of a learned cut, +infinity otherwise. */
SCIP_Real getSlackImpliedUb(SCIP* scip, int idx);

std::vector<SCIP_Real> FusionResolution(SCIP* scip, SCIP_Longint child0, SCIP_Longint child1, SCIP_Real& rhs);
SCIP_RETCODE FinishLearn(SCIP* scip, SCIP_Longint node);
SCIP_RETCODE AnalyseDualSol(SCIP* scip, SCIP_Longint nodeNumber);

SCIP_RETCODE gatherLPInfeasibilityCertificate(
    SCIP* scip,
    SCIP_Longint nodeNumber);

/** Proof of a node cut off by its pseudo objective: the objective of the variables fixedTo1 (AllVars indices),
 *  relative to the dual solution of the parent. */
SCIP_RETCODE getDualProofFromPseudoObjCutoff(
    SCIP* scip,
    SCIP_Longint nodeNumber, SCIP_Longint parentNumber, const std::vector<int>& fixedTo1);

/** Writes the LP and the transformed problem of the node when option -o=N is given with N > 1 and node >= N. */
SCIP_RETCODE writeNodeLPFilesIfRequested(SCIP* scip, SCIP_Longint nodenumber);

/* Steps shared by FinishLearn and the node-deletion event to build the proof of a node from the proofs of its children. */

/** Debug checks on the proof of the only child of a node. */
void checkSingleChildProof(SCIP_Longint nodenumber, SCIP_Longint child);

/** True if rhs is the rhs of an infeasibility proof (TOP, possibly lowered by dropNegligibleCoefs). */
bool isInfeasibilityRhs(SCIP_Real rhs);

/** Queues the clause sum x_i >= 1 over the TOP literals of the infeasibility proof sum coefs[i] x_i >= rhs. */
SCIP_RETCODE queueInfeasibilityClause(SCIP* scip, SCIP_Longint nodenumber, std::vector<SCIP_Real> coefs, SCIP_Real rhs);

/** If both children increased the lower bound, fuses their proofs into ResFus with right-hand side rhsFus.
 *  infeasibleFusion is set if the fused proof proves the node infeasible; its clause version is then learned. */
SCIP_RETCODE fuseChildrenProofs(SCIP* scip, SCIP_Longint nodenumber, const std::vector<long>& children,
    std::vector<SCIP_Real>& ResFus, SCIP_Real& rhsFus, bool& infeasibleFusion);

/** Adds the node's own dual proof (if it increased the lower bound) to ResFus and rhs_summed. */
void addNodeProof(SCIP_Longint nodenumber, std::vector<SCIP_Real>& ResFus, SCIP_Real& rhs_summed);

/** Replaces variables removed by bound propagation with their explanation. */
void applyBoundPropagationRemovals(SCIP_Longint nodenumber, std::vector<SCIP_Real>& ResFus, SCIP_Real rhs_summed);

/** Replaces variables removed by reduced cost filtering (node consistency) with the parent reduced costs. */
void applyNodeConsistency(SCIP* scip, SCIP_Longint nodenumber, SCIP_Longint parentNumber,
    std::vector<SCIP_Real>& ResFus, SCIP_Real& rhs_summed);

/** Debug check (option -d) that rhs_summed matches the lower bound increases of the node and its children. */
void checkSummedRhs(SCIP_Longint nodenumber, const std::vector<long>& children, SCIP_Real rhs_summed,
    bool infeasibleFusion);

/** Removes the children of nodenumber from nodeInfoMap. */
void eraseChildren(SCIP_Longint nodenumber, const std::vector<long>& children);

/** Prints the nonzero coefficients of coefs followed by "<relation> rhs".
 *  If mapThroughCurrentVarIdx, variable i is printed as AllVars[CurrentVarIdx[i]] instead of AllVars[i]. */
void printProof(SCIP* scip, const std::vector<SCIP_Real>& coefs, const char* relation, SCIP_Real rhs,
    bool mapThroughCurrentVarIdx);

#endif
