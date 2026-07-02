#ifndef MEMOSCIP_TYPES_H
#define MEMOSCIP_TYPES_H

#include "scip/scip.h"
#include "sepa_mycuts.h"

#include <cmath>
#include <map>
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
    std::vector<long> children;
    int BranchVar = 0;
    std::vector<SCIP_Real> redcosts;
    int nvarsparent = 0;
    SCIP_CONS* objcons = nullptr;
    std::string Status = "Unknown";
    bool LbgUb = false;
    
    std::pair<std::vector<SCIP_Real>, SCIP_Real> DualProofCoef;
    std::vector<int> DualProofInds;

    std::pair<std::vector<SCIP_Real>, SCIP_Real> FirstDualProofCoef;
    std::vector<int> FirstDualProofInds;

    SCIP_Longint parent = 0;
    std::vector<SCIP_Real> DualSol;
    std::vector<SCIP_Real> y_i;

    std::vector<std::pair<SCIP_Longint, SCIP_Real>> NodeConsistency;

    std::vector<int> farkasRowLPIndex;
    std::vector<SCIP_Real> farkasRowCoef;
    std::vector<std::string> farkasRowNames;
    std::vector<SCIP_Real> farkasColCoef;
    std::vector<SCIP_VAR*> farkasColVars;

    std::vector<std::pair<int, std::vector<int>>> BPRemoval;

    std::vector<int> VarLbDebug;
    std::vector<int> VarUbDebug;
};


/** Shared solver / learning state (single instance, see g_memoscip). */
struct MemoSCIPContext {
    int counter = 0;
    int nbAddedctr = 0;
    int nOriginalVar = 0;
    std::unordered_map<SCIP_Longint, NodeInfo> nodeInfoMap;
    std::vector<int> AuxVar;
    std::vector<int> slackVar;
    std::vector<SCIP_VAR*> AllVars;
    std::vector<int> CurrentVarIdx;
    MyQueue CutQueue;
    std::map<int, int> eqvid;
    SolCollector collector;
    std::vector<SCIP_Real> VarLb;
    std::vector<SCIP_Real> VarUb;
    std::vector<SCIP_VAR*> OrderVar;
    std::string FileExt;
    long long TOP = 1e+9;
};

/** Single shared instance; inline (C++17) so no separate .cpp definition is required. */
inline MemoSCIPContext g_memoscip{};

#endif
