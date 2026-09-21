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
    bool LPsolved=false;
    std::vector<long> children;
    int BranchVar = 0;
    int nvarsparent = 0;
    SCIP_CONS* objcons = nullptr;
    std::string Status = "Unknown";
    bool LbgUb = false;
    
    std::vector<int> Var1;
    std::pair<std::vector<SCIP_Real>, SCIP_Real> DualProofCoef;
    std::vector<int> DualProofInds;

    std::pair<std::vector<SCIP_Real>, SCIP_Real> FirstDualProofCoef;
    std::vector<int> FirstDualProofInds;

    SCIP_Longint parent = -1;
    std::vector<SCIP_Real> DualSol;

    std::vector<SCIP_Real> y_i;   // dual solution of constraints x_i<=1. Based on AllVars index. y_i[i] = dual solution of constraint AllVars[i]<=1
    std::vector<SCIP_Real> redcosts; // redcosts based on AllVars index. redcosts[i] = redcost of AllVars[i]

    std::vector<std::pair<SCIP_Longint, SCIP_Real>> NodeConsistency;

    std::vector<int> farkasRowLPIndex;
    std::vector<SCIP_Real> farkasRowCoef;
    std::vector<std::string> farkasRowNames;
    std::vector<SCIP_Real> farkasColCoef;
    std::vector<SCIP_VAR*> farkasColVars;

    std::vector<std::pair<int, std::vector<int>>> BPRemoval;

    std::vector<SCIP_CONS*> AddedConstraints;

};


/** Shared solver / learning state (single instance, see g_memoscip). */
struct MemoSCIPContext {
    int counter = 0;
    int nbctr = 0;
    int nbAddedctr = 0;
    int nOriginalVar = 0;
    int addedslack=0;
    std::unordered_map<SCIP_Longint, NodeInfo> nodeInfoMap;
    std::vector<int> AuxVar;  //Auxiliary variable index. AuxVar[i] = index of the auxilliary variable in AllVars. 
    std::vector<SCIP_VAR*> slackVar;
    std::vector<SCIP_VAR*> AllVars; //Contains all variables including auxiliary and slack variables. Always same ordering.
    std::vector<int> CurrentVarIdx; //Transforms SCIP index to AllVars index. CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar] = index in AllVars
    MyQueue CutQueue;
    std::map<int, int> eqvid;
    SolCollector collector;
    std::vector<SCIP_Real> VarLb; //Index are based on AllVars. VarLb[i] = lower bound of AllVars[i]
    std::vector<SCIP_Real> VarUb; //Index are based on AllVars. VarUb[i] = upper bound of AllVars[i]
    std::vector<SCIP_VAR*> OrderVar;
    std::vector<SCIP_CONS*> SlackNonNegConss;
    //std::vector<SCIP_VAR*> Slackvars;
    std::string FileExt;
    long long TOP = 1e+9;
};

/** Single shared instance; inline (C++17) so no separate .cpp definition is required. */
inline MemoSCIPContext g_memoscip{};

#endif
