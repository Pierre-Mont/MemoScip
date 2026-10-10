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

struct NodeInfo {
    int number = 0;
    SCIP_Real lb_increase = 0.0;
    SCIP_Real lb = 0.0;
    SCIP_Real ub = 0.0;
    bool visited = false;
    bool LPsolved = false;
    std::vector<long> children;
    int BranchVar = 0;
    std::string Status = "Unknown";
    bool LbgUb = false;

    std::vector<int> Var1;  // AllVars indices of the variables fixed to 1 at this node
    std::pair<std::vector<SCIP_Real>, SCIP_Real> DualProofCoef;  // (coefficients, rhs) of the node dual proof
    std::vector<int> DualProofInds;  // AllVars index of each coefficient of DualProofCoef

    SCIP_Longint parent = -1;
    std::vector<SCIP_Real> DualSol;  // dual solution indexed by SCIP row index

    std::vector<SCIP_Real> y_i;   // dual solution of constraints x_i<=1. Based on AllVars index. y_i[i] = dual solution of constraint AllVars[i]<=1
    std::vector<SCIP_Real> redcosts; // redcosts based on AllVars index. redcosts[i] = redcost of AllVars[i]

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
    std::vector<int> AuxVar;  //Auxiliary variable index. AuxVar[i] = index of the auxilliary variable in AllVars.
    std::vector<std::vector<int>> TupleAuxVar; //Auxilliary variables indices for tuple variables (Encoding=1)
    std::vector<SCIP_VAR*> AllVars; //Contains all variables including auxiliary and slack variables. Always same ordering.
    std::vector<int> CurrentVarIdx; //Transforms SCIP index to AllVars index. CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar] = index in AllVars
    MyQueue CutQueue;
    SolCollector collector;
    std::vector<SCIP_Real> VarLb; //Index are based on AllVars. VarLb[i] = lower bound of AllVars[i]
    std::vector<SCIP_Real> VarUb; //Index are based on AllVars. VarUb[i] = upper bound of AllVars[i]
    std::vector<SCIP_VAR*> OrderVar;
    std::string FileExt;
    std::vector<SCIP_CONS*> Unconstructed_cons; //Constraint added to the problem but not yet added in the LP


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

/** SCIP index of var shifted by nOriginalVar: position of var in CurrentVarIdx, AuxVar and TupleAuxVar. */
inline int shiftedIndex(SCIP_VAR* var)
{
    return SCIPvarGetIndex(var) - g_memoscip.nOriginalVar;
}

/** Index of var in AllVars (and VarLb, VarUb, ...). */
inline int allVarsIndex(SCIP_VAR* var)
{
    return g_memoscip.CurrentVarIdx[shiftedIndex(var)];
}

/** Slack variables are recognised by their name: "slack_*", "sc_*", or the transformed "t_slack_*", "t_sc_*". */
inline bool isSlackVar(SCIP_VAR* var)
{
    const char* name = SCIPvarGetName(var);
    return name[0] == 's' || name[2] == 's';
}

/** Calls f on the AllVars index of each literal complementary to the variable at shifted index pos:
 *  its auxiliary variable, or the other variables of its tuple constraint when encoding 1 is used.
 *  If checkVarTuple, a variable that belongs to no tuple constraint falls back to its auxiliary variable;
 *  otherwise the fallback only happens when the whole TupleAuxVar table is empty. */
template <class F>
inline void forEachAuxVar(int pos, bool checkVarTuple, F&& f)
{
    const bool useTuple = g_options.encoding != 0
        && (checkVarTuple ? !g_memoscip.TupleAuxVar[pos].empty() : !g_memoscip.TupleAuxVar.empty());
    if (!useTuple)
        f(g_memoscip.AuxVar[pos]);
    else
        for (int aux : g_memoscip.TupleAuxVar[pos])
            f(aux);
}

#endif
