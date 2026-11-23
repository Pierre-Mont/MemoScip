#ifndef MEMOSCIP_TYPES_H
#define MEMOSCIP_TYPES_H

#include "memoscip_option.h"
#include "scip/scip.h"
#include "sepa_mycuts.h"

#include <chrono>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

/** Statistics on the variables of a node proof (non negligible coefficients only), by sign of the coefficient:
 *  [0] positive, [1] negative. "atNode" = fixed at the depth of the node (branching or propagation at the node),
 *  "before" = fixed at an ancestor, globally or in presolving, "unfixed" = not fixed at the node, "cutSlack" = slack
 *  of a learned cut, "modelSlack" = slack of a model inequality. A variable is fixed if its MemoSCIP bounds
 *  (VarLb/VarUb, which include the fixings implied by the domain constraints) are equal; its depth is the earliest
 *  depth at which it or one of its auxiliary variables was fixed in SCIP. */
struct ProofVarStats {
    long long proofs = 0;
    long long emptyProofs = 0;
    long long negligible = 0;          // coefficients dropped by dropNegligibleCoefs
    long long noBefore = 0;            // non empty proofs without any variable fixed before the node
    long long unfixed[2] = {0, 0};
    long long atNode[2] = {0, 0};
    long long atNode1[2] = {0, 0};     // part of atNode assigned to 1
    long long before[2] = {0, 0};
    long long before1[2] = {0, 0};     // part of before assigned to 1
    long long cutSlack[2] = {0, 0};
    long long modelSlack[2] = {0, 0};
    long long unknownDepth[2] = {0, 0};  // fixed in the MemoSCIP bounds without a SCIP fixing (should stay 0)

    void add(const ProofVarStats& o, int sign)
    {
        proofs += sign * o.proofs;
        emptyProofs += sign * o.emptyProofs;
        negligible += sign * o.negligible;
        noBefore += sign * o.noBefore;
        for (int k = 0; k < 2; ++k) {
            unfixed[k] += sign * o.unfixed[k];
            atNode[k] += sign * o.atNode[k];
            atNode1[k] += sign * o.atNode1[k];
            before[k] += sign * o.before[k];
            before1[k] += sign * o.before1[k];
            cutSlack[k] += sign * o.cutSlack[k];
            modelSlack[k] += sign * o.modelSlack[k];
            unknownDepth[k] += sign * o.unknownDepth[k];
        }
    }
};

/** Use of the learned cuts by the node LPs, from the full dual solution of each analysed LP (not in dives/probing). */
struct CutUsageStats {
    long long lps = 0;               // analysed node LPs
    long long lpsWithCuts = 0;       // ... with at least one learned cut row in the LP
    long long cutRowsInLP = 0;       // sum over the analysed LPs of the number of learned cut rows
    long long lpsUsingCuts = 0;      // ... with at least one learned cut row with a nonzero dual
    long long cutDuals = 0;          // sum over the analysed LPs of the number of learned cut rows with a nonzero dual
    double boundShare = 0.0;         // sum over the LPs using cuts of (sum of dual * rhs of the cut rows) / LP value
};

struct NodeInfo {
    int number = 0;
    // lb_increase is the proven bound increase (LP gain, or subtree gain once the node is deleted) and
    // DualProofCoef.second the rhs of the stored proof. dropNegligibleCoefs lowers the rhs to stay valid, so the two
    // differ: drop_loss is an upper bound on lb_increase - DualProofCoef.second (>= 0), used by the debug checks
    SCIP_Real lb_increase = 0.0;
    SCIP_Real drop_loss = 0.0;
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
    std::vector<SCIP_Real> LPSol;  // option -mir: primal LP solution of the node (AllVars index), reference point of its cuts
    std::vector<SCIP_Real> LPLb, LPUb;  // with -comp: local bounds (AllVars index) of the LP whose duals are DualSol
    // Option -split
    std::vector<std::pair<int,int>> CompLabel;  // sorted (AllVars index, component id) of the unfixed columns of the node LP
    std::vector<ProofPart> IncParts;  // node increment (DualProofCoef) split by component of the node LP
    std::vector<ProofPart> Parts;     // proof of the node subtree split by component, set when the node is deleted
    bool HasParts = false;

    std::vector<std::pair<SCIP_Longint, SCIP_Real>> NodeConsistency;  // (variable removed by reduced cost filtering, UB - parent lb)

    std::vector<std::pair<int, std::vector<int>>> BPRemoval;  // (variable removed by bound propagation, explanation)

    std::vector<SCIP_CONS*> AddedConstraints;

    // Statistics of the last proof computed at the node (its LP can be analysed twice), see MemoSCIPContext::ProofStats
    ProofVarStats ProofStats;
    int ProofStatsKind = -1;  // -1: none, 0: LP dual proof, 1: infeasibility proof
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
    // Option -subslack: definition of the learned slack s_k = a_k x - b_k, indexed by AllVars (a_k over AllVars indices
    // of non learned-slack variables). Only filled for the cuts added with the option.
    std::vector<char> SlackHasDef;
    std::vector<std::vector<std::pair<int, SCIP_Real>>> SlackDefTerms;
    std::vector<SCIP_Real> SlackDefRhs;
    // Option -submodel: slacks of model inequalities (slack_<i>) whose definition is in SlackDef*. They are substituted
    // whatever the sign of their coefficient (an identity), unlike learned slacks whose negative terms are dropped.
    std::vector<char> SlackIsModel;
    long long nbSubstCuts = 0;    // cuts in which at least one learned slack was substituted
    long long nbSubstTrivial = 0; // substituted cuts skipped because trivial
    // Option -mir: learned cuts for which a c-MIR cut was computed, added, equal to the cut, or not found
    long long nbMirTried = 0;
    long long nbMirAdded = 0;
    long long nbMirSame = 0;
    long long nbMirFailed = 0;
    long long nbMirWeaker = 0;       // MIR cut not more efficacious than the learned cut on the node LP solution
    long long nbMirCutViolated = 0;  // learned cuts violated by the LP solution of their node
    long long nbMirViolated = 0;     // MIR cuts violated by the LP solution of their node


    // Option -comp: nodes whose duals were mixed per connected component, kept unchanged, or not mixed (bound check failed)
    long long DualMixApplied = 0;
    long long DualMixUnchanged = 0;
    long long DualMixRejected = 0;
    long long SplitNodes = 0;       // node deletions handled per component (-split)
    long long SplitFallback = 0;    // node deletions that fell back to the single proof fusion
    long long SplitPassengers = 0;  // parts of other components passed from a child to the node
    long long SplitMerged = 0;      // child parts that joined another component at the node

    // Statistics on the node proofs: [0] LP dual proofs (AnalyseDualSol), [1] infeasibility proofs
    ProofVarStats ProofStats[2];
    CutUsageStats CutUsage;
    std::unordered_set<int> UsedCutRows;  // SCIP row index of the learned cuts that had a nonzero dual in a node LP

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
 *  finite side, which is the value of the original equality (the slack term is kept in the aggregation). A row made
 *  free in the dive (slacks of both signs missing from the LP) gets +infinity, so that the proofs skip it as a row
 *  that is not in the LP. */
inline SCIP_Real equalityRowValue(SCIP* scip, SCIP_ROW* row)
{
    if (isDiveRelaxedRow(row) && SCIPisInfinity(scip, SCIProwGetRhs(row)))
        return SCIPisInfinity(scip, -SCIProwGetLhs(row)) ? SCIPinfinity(scip) : SCIProwGetLhs(row);
    return SCIProwGetRhs(row);
}

/** Learned cut rows are recognised by the names given to the cuts queued with SCIPqueueCut ("cut_*", "l2cut_*",
 *  "clause_*"): sepa_mycuts appends their slack sc_<k> as the last variable. Model rows can have as many nonzeros
 *  (e.g. a tuple marginal has d+1, a tuple sum d*d with encoding 1), so the size of a row does not identify a cut. */
inline bool isLearnedCutRow(SCIP_ROW* row)
{
    const char* name = SCIProwGetName(row);
    return std::strncmp(name, "cut_", 4) == 0 || std::strncmp(name, "l2cut_", 6) == 0 || std::strncmp(name, "clause_", 7) == 0;
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
