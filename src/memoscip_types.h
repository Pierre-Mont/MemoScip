#ifndef MEMOSCIP_TYPES_H
#define MEMOSCIP_TYPES_H

#include "memoscip_option.h"
#include "scip/scip.h"
#include "sepa_mycuts.h"

#include <chrono>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct SolCollector {
    std::vector<SCIP_Real> dualsol;
    std::vector<SCIP_Real> redcost;
    SCIP_Real objval = 0.0;

    void init(SCIP* scip);
    SCIP_RETCODE collect(SCIP* scip);
};

/** Reduced costs d = c - y^T A of the LP columns for the dual solution y (indexed by LP row position). */
void computeRedcosts(SCIP* scip, const std::vector<SCIP_Real>& dualsol, std::vector<SCIP_Real>& redcost);

/** A piece of a dual proof (option -split): sparse coefficients (AllVars index) and its own rhs. A proof is a
 *  combination of equality rows (bound duals go through the domain equality x + aux = 1), so splitting it row by row
 *  gives pieces that are valid on their own and sum to the whole proof. */
struct ProofPart {
    std::vector<int> inds;
    std::vector<SCIP_Real> coefs;
    SCIP_Real rhs = 0.0;
    bool orphan = false;  // node increment piece with no unfixed variable at the node (newly fixed variables, fully fixed rows)
};

struct NodeInfo {
    int number = 0;
    SCIP_Real lb_increase = 0.0;
    SCIP_Real lb = 0.0;
    SCIP_Real ub = 0.0;
    bool visited = false;
    bool LPsolved = false;
    std::vector<long> children;
    int BranchVar = 0;  // shifted index of the branching variable of this node
    bool BranchVal = true;  // true if BranchVar is fixed to 1 at this node, false if it is fixed to 0
    std::string Status = "Unknown";
    bool LbgUb = false;

    std::vector<int> Var1;  // AllVars indices of the variables fixed to 1 at this node
    std::pair<std::vector<SCIP_Real>, SCIP_Real> DualProofCoef;  // (coefficients, rhs) of the node dual proof
    std::vector<int> DualProofInds;  // AllVars index of each coefficient of DualProofCoef

    SCIP_Longint parent = -1;
    std::vector<SCIP_Real> DualSol;  // dual solution indexed by SCIP row index

    std::vector<SCIP_Real> y_i;   // dual solution of constraints x_i<=1. Based on AllVars index. y_i[i] = dual solution of constraint AllVars[i]<=1
    std::vector<SCIP_Real> redcosts; // redcosts based on AllVars index. redcosts[i] = redcost of AllVars[i]
    std::vector<SCIP_Real> LPLb, LPUb;  // with -comp: local bounds (AllVars index) of the LP whose duals are DualSol
    // Option -split
    std::vector<std::pair<int,int>> CompLabel;  // sorted (AllVars index, component id) of the unfixed columns of the node LP
    std::vector<ProofPart> IncParts;  // node increment (DualProofCoef) split by component of the node LP
    std::vector<ProofPart> Parts;     // proof of the node subtree split by component, set when the node is deleted
    bool HasParts = false;

    std::vector<std::pair<SCIP_Longint, SCIP_Real>> NodeConsistency;  // (variable removed by reduced cost filtering, UB - parent lb)

    std::vector<std::pair<int, std::vector<int>>> BPRemoval;  // (variable removed by bound propagation, explanation)

    std::vector<SCIP_CONS*> AddedConstraints;
};


/** Shared solver / learning state (single instance, see g_memoscip). */
struct MemoSCIPContext {
    int nbctr = 0;
    int nbAddedctr = 0;
    int nOriginalVar = 0;
    std::unordered_map<SCIP_Longint, NodeInfo> nodeInfoMap;
    std::vector<std::vector<int>> AuxVar;  //Auxiliary variable indices. AuxVar[i] = AllVars indices of the other variables of the exactly-one constraint of i (domain constraint, or tuple sum in encoding 1: sum_j x_j = 1).
    std::vector<SCIP_VAR*> AllVars; //Contains all variables including auxiliary and slack variables. Always same ordering.
    std::vector<int> CurrentVarIdx; //Transforms SCIP index to AllVars index. CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar] = index in AllVars
    MyQueue CutQueue;
    SolCollector collector;
    std::vector<SCIP_Real> VarLb; //Index are based on AllVars. VarLb[i] = lower bound of AllVars[i]
    std::vector<SCIP_Real> VarUb; //Index are based on AllVars. VarUb[i] = upper bound of AllVars[i]
    std::vector<SCIP_VAR*> OrderVar;
    std::string FileExt;
    std::vector<SCIP_CONS*> Unconstructed_cons; //Constraint added to the problem but not yet added in the LP
    std::vector<SCIP_ROW*> DiveRelaxedRows; //Learned cut rows relaxed to one side in the current dive because their slack column is not in the LP
    // Learned cut slacks s_k >= 0 (row a_k x - s_k = b_k) have no upper bound, but s_k <= max activity of a_k x - b_k
    // over the global bounds. Indexed by AllVars: the implied upper bound (NaN while not computed from the LP row).
    std::vector<SCIP_Real> SlackImpliedUb;
    std::vector<SCIP_CONS*> SlackCons; //AllVars index -> linear constraint of the cut of the slack (nullptr for separator rows, whose bound is computed at creation)


    // Option -comp: nodes whose duals were mixed per connected component, kept unchanged, or not mixed (bound check failed)
    long long DualMixApplied = 0;
    long long DualMixUnchanged = 0;
    long long DualMixRejected = 0;
    long long SplitNodes = 0;       // node deletions handled per component (-split)
    long long SplitFallback = 0;    // node deletions that fell back to the single proof fusion
    long long SplitPassengers = 0;  // parts of other components passed from a child to the node
    long long SplitMerged = 0;      // child parts that joined another component at the node

    long long TOP = 1e+9;
    std::chrono::duration<double> AnalyseDualSolTime{0.0}; // cumulated time spent in AnalyseDualSol
    std::chrono::duration<double> NodeDeleteTime{0.0}; // cumulated time spent in eventExecNodeDelete
    std::chrono::duration<double> NodeSolvedTime{0.0}; // cumulated time spent in eventExecNodeSolved
    std::chrono::duration<double> LPSolvedTime{0.0}; // cumulated time spent in eventExecLPSolved
    std::chrono::duration<double> NodeFocusedTime{0.0}; // cumulated time spent in eventExecNodeFocused
    std::chrono::duration<double> BoundChangedTime{0.0}; // cumulated time spent in eventExecBoundChanged
    std::chrono::duration<double> GatherLPInfeasibilityCertificateTime{0.0}; // cumulated time spent in gatherLPInfeasibilityCertificate
    std::chrono::duration<double> GetDualProofFromPseudoObjCutoffTime{0.0}; // cumulated time spent in getDualProofFromPseudoObjCutoff
};

/** Single shared instance; inline (C++17) so no separate .cpp definition is required. */
inline MemoSCIPContext g_memoscip{};


// Accumulates elapsed time into acc for the whole lifetime of the enclosing scope, covering every early return.
struct ScopedTimer {
    std::chrono::duration<double>& acc;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    explicit ScopedTimer(std::chrono::duration<double>& acc_) : acc(acc_) {}
    ~ScopedTimer() { acc += std::chrono::steady_clock::now() - start; }
};

/** SCIP index of var shifted by nOriginalVar: position of var in CurrentVarIdx and AuxVar. */
inline int shiftedIndex(SCIP_VAR* var)
{
    return SCIPvarGetIndex(var) - g_memoscip.nOriginalVar;
}

/** Index of var in AllVars (and VarLb, VarUb, ...). */
inline int allVarsIndex(SCIP_VAR* var)
{
    return g_memoscip.CurrentVarIdx[shiftedIndex(var)];
}

/** True if row is a learned cut row relaxed to one side in the current dive (its slack column is not in the LP). */
inline bool isDiveRelaxedRow(SCIP_ROW* row)
{
    for (SCIP_ROW* r : g_memoscip.DiveRelaxedRows)
        if (r == row)
            return true;
    return false;
}

/** Value of the equality "row = value" used in the proofs: the rhs, or for a row relaxed in the dive its remaining
 *  finite side, which is the value of the original equality (the slack term is kept in the aggregation). */
inline SCIP_Real equalityRowValue(SCIP* scip, SCIP_ROW* row)
{
    if (isDiveRelaxedRow(row) && SCIPisInfinity(scip, SCIProwGetRhs(row)))
        return SCIProwGetLhs(row);
    return SCIProwGetRhs(row);
}

/** Slack variables are recognised by their name: "slack_*", "sc_*", or the transformed "t_slack_*", "t_sc_*". */
inline bool isSlackVar(SCIP_VAR* var)
{
    const char* name = SCIPvarGetName(var);
    return name[0] == 's' || name[2] == 's';
}

/** Calls f on the AllVars index of each literal complementary to the variable at shifted index pos: its auxiliary
 *  variables, i.e. the other variables of its domain constraint or of its tuple sum constraint. */
template <class F>
inline void forEachAuxVar(int pos, F&& f)
{
    for (int aux : g_memoscip.AuxVar[pos])
        f(aux);
}

/** True if pred holds for at least one auxiliary variable (AllVars index) of the variable at shifted index pos. */
template <class P>
inline bool anyAuxVar(int pos, P&& pred)
{
    for (int aux : g_memoscip.AuxVar[pos])
        if (pred(aux))
            return true;
    return false;
}

/** True if the variable at shifted index pos has auxiliary variables and pred holds for all of them. */
template <class P>
inline bool allAuxVars(int pos, P&& pred)
{
    const std::vector<int>& auxs = g_memoscip.AuxVar[pos];
    for (int aux : auxs)
        if (!pred(aux))
            return false;
    return !auxs.empty();
}

#endif
