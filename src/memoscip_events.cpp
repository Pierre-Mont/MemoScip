#include "memoscip_events.h"

#include "memoscip_dual.h"
#include "memoscip_option.h"
#include "memoscip_types.h"
#include "sepa_mycuts.h"

#include "scip/cons_linear.h"
#include "scip/scip.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace {

struct ReplayDecision {
    std::string varname;
    int value = 0;
};

std::vector<ReplayDecision> g_replayDecisions;
size_t g_replayDecisionPos = 0;
bool g_replayLoaded = false;
std::vector<SCIP_VAR*> g_boundChangedEventVars;

const char* kReplayBranchingFile = "/home/kisouke/MemoSCIP/Instance/branch_replay_sequence.txt";

std::string trim_copy(const std::string& s)
{
    const size_t first = s.find_first_not_of(" \t\n\r");
    if (first == std::string::npos)
        return "";
    const size_t last = s.find_last_not_of(" \t\n\r");
    return s.substr(first, last - first + 1);
}

/** Parses "var1=0; var2=1; ..." into a list of branching decisions. */
bool parse_replay_sequence(const std::string& text, std::vector<ReplayDecision>& out)
{
    out.clear();
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ';'))
    {
        token = trim_copy(token);
        if (token.empty())
            continue;

        const size_t eqpos = token.find('=');
        if (eqpos == std::string::npos)
            return false;

        std::string lhs = trim_copy(token.substr(0, eqpos));
        std::string rhs = trim_copy(token.substr(eqpos + 1));
        if (lhs.empty() || rhs.empty())
            return false;

        if (rhs != "0" && rhs != "1")
            return false;

        ReplayDecision d;
        d.varname = lhs;
        d.value = (rhs == "1") ? 1 : 0;
        out.push_back(d);
    }
    return !out.empty();
}

/** Finds a variable by name, also trying the transformed name "t_<name>". */
SCIP_VAR* find_replay_var(SCIP* scip, const std::string& name)
{
    SCIP_VAR** vars = SCIPgetVars(scip);
    const int nvars = SCIPgetNVars(scip);

    for (int i = 0; i < nvars; ++i)
    {
        if (name == SCIPvarGetName(vars[i]))
            return vars[i];
    }

    if (name.rfind("t_", 0) != 0)
    {
        const std::string alt = "t_" + name;
        for (int i = 0; i < nvars; ++i)
        {
            if (alt == SCIPvarGetName(vars[i]))
                return vars[i];
        }
    }

    return nullptr;
}

bool ensure_replay_loaded()
{
    if (g_replayLoaded)
        return !g_replayDecisions.empty();

    g_replayLoaded = true;
    std::ifstream in(kReplayBranchingFile);
    if (!in)
        return false;

    std::ostringstream buf;
    buf << in.rdbuf();
    return parse_replay_sequence(buf.str(), g_replayDecisions);
}

}

/* ---------------------------------------------------------------------------------------------------------------- */
/* Node solved                                                                                                       */
/* ---------------------------------------------------------------------------------------------------------------- */

/** Initialises AllVars, CurrentVarIdx, VarLb, VarUb, AuxVar and TupleAuxVar from the columns of the first LP. */
static void initVarIndexMaps(SCIP* scip)
{
    int nconss = SCIPgetNConss(scip);
    SCIP_CONS** conss = SCIPgetConss(scip);
    g_memoscip.nbctr=nconss;
    int ncols = SCIPgetNLPCols(scip);
    SCIP_COL** cols = SCIPgetLPCols(scip);
    auto pos = [ncols](SCIP_VAR* v) { return SCIPvarGetIndex(v)-ncols; };

    int maxindex = SCIPvarGetIndex(SCIPcolGetVar(cols[ncols-1]));
    assert(maxindex>=2*ncols-1);
    assert(g_memoscip.nOriginalVar==ncols);
    g_memoscip.CurrentVarIdx.resize(maxindex+1-ncols,0);
    g_memoscip.VarLb.resize(ncols,0);
    g_memoscip.VarUb.resize(ncols,1);
    g_memoscip.AuxVar.resize(maxindex+1-ncols,-1);
    g_memoscip.TupleAuxVar.resize(maxindex+1-ncols);

    for( int i = 0; i < ncols; i++ )
    {
        SCIP_VAR* var = SCIPcolGetVar(cols[i]);
        g_memoscip.AllVars.push_back(var);
        g_memoscip.CurrentVarIdx[pos(var)]=i;
        g_memoscip.VarLb[i]=SCIPvarGetLbLocal(var);
        if(!isSlackVar(var) && SCIPvarGetUbLocal(var)==SCIPinfinity(scip))
            g_memoscip.VarUb[i]=1;
        else
            g_memoscip.VarUb[i]=SCIPvarGetUbLocal(var);

        // By default a variable is its own auxiliary variable
        if(g_memoscip.AuxVar[pos(var)]==-1)
            g_memoscip.AuxVar[pos(var)]=pos(var);
    }

    // Domain constraints x_i0 + x_i1 = 1 link each variable to its auxiliary variable
    bool foundAuxVars=false;
    for(int i = 0; i < nconss; ++i)
    {
        SCIP_CONS* cons = conss[i];
        SCIP_VAR** vars = SCIPgetVarsLinear(scip,cons);
        int nvars = SCIPgetNVarsLinear(scip,cons);
        const char* consname = SCIPconsGetName(cons);
        if((g_memoscip.FileExt==".opb" && consname[0]=='u' && consname[1]=='b')
            || ((g_memoscip.FileExt==".lp" || g_options.kpcg) && SCIPgetRhsLinear(scip, cons)==1 && nvars==2)){
            g_memoscip.AuxVar[pos(vars[0])]=g_memoscip.CurrentVarIdx[pos(vars[1])];
            g_memoscip.AuxVar[pos(vars[1])]=g_memoscip.CurrentVarIdx[pos(vars[0])];
            foundAuxVars=true;
        }
        // Encoding 1: each tuple variable of a constraint sum of 4 tuples = 1 is linked to the 3 other ones
        if(g_options.encoding==1 && SCIPgetRhsLinear(scip, cons)==1 && nvars==4){
            for(int k = 0; k < 4; ++k){
                vector<int>& others = g_memoscip.TupleAuxVar[pos(vars[k])];
                others.clear();
                for(int j = 0; j < 4; ++j)
                    if(j != k)
                        others.push_back(g_memoscip.CurrentVarIdx[pos(vars[j])]);
            }
        }
    }
    if(!foundAuxVars)
        cout<<"-----------------No auxiliary variables found !!!!--------------"<<endl;
}

/** If the lower bound didn't change then we copy the parent information. */
static void copyParentInfo(SCIP_Longint nodenumber)
{
    NodeInfo& node = g_memoscip.nodeInfoMap[nodenumber];
    const NodeInfo& parent = g_memoscip.nodeInfoMap[node.parent];
    node.redcosts = parent.redcosts;
    node.lb=parent.lb;
    node.y_i=parent.y_i;
    node.DualSol=parent.DualSol;
    node.DualProofCoef.second=0;
    node.Status="FEASIBLE";
    node.lb_increase=0;
    for (size_t i = 0; i < parent.y_i.size(); i++){
        if(g_memoscip.VarLb[i]==1){
            node.Var1.push_back(allVarsIndex(g_memoscip.AllVars[i]));
        }
    }
}

/** Builds the proof of the node from its solved LP. Sets handled to false if the LP status is not
 *  infeasible, optimal or objective limit. */
static SCIP_RETCODE processNodeLP(SCIP* scip, SCIP_LPSOLSTAT stat, SCIP_Longint nodenumber, SCIP_Longint parentnumber,
    bool& handled)
{
    handled = true;
    if(stat==SCIP_LPSOLSTAT_INFEASIBLE){
        if(g_options.verbosity >= 1)
            cout<<"LP INFEASIBLE "<<endl;
        SCIP_CALL(gatherLPInfeasibilityCertificate(scip, nodenumber));
    }else if(stat==SCIP_LPSOLSTAT_OPTIMAL || stat==SCIP_LPSOLSTAT_OBJLIMIT){
        SCIP_Real lpobj=SCIPgetLPObjval(scip);
        if(g_options.verbosity >= 1)
            cout<<"LP OPTIMAL: "<<lpobj<<endl;
        // Only analyse the LP if the lower bound increased enough
        SCIP_Real parentlb = g_memoscip.nodeInfoMap[parentnumber].lb;
        if(!SCIPisDualfeasZero(scip,lpobj) && (fabs(parentlb - lpobj)>= g_options.minInc || Myceil(lpobj)>=SCIPgetPrimalbound(scip) || lpobj - Myceil(parentlb) > 1e-5))
            AnalyseDualSol(scip,nodenumber );
        else
            copyParentInfo(nodenumber);
    }else{
        handled = false;
    }
    return SCIP_OKAY;
}

/** AllVars index of the variable fixed to 1 by the branching decision that created node,
 *  or -1 if node was not created by a single branching x >= 1. */
static int getBranchingVarFixedTo1(SCIP_NODE* node)
{
    int nvars;
    // First query the required size with branchvarssize = 0
    SCIPnodeGetParentBranchings(node, NULL, NULL, NULL, &nvars, 0);
    if(nvars != 1)
        return -1;
    SCIP_VAR* var;
    SCIP_Real bound;
    SCIP_BOUNDTYPE boundtype;
    SCIPnodeGetParentBranchings(node, &var, &bound, &boundtype, &nvars, 1);
    if(boundtype != SCIP_BOUNDTYPE_LOWER || bound != 1)
        return -1;
    return allVarsIndex(var);
}

/** Pseudo objective of the child of parentnumber where AllVars[decVar] is fixed to 1: objective of the variables
 *  fixed to 1 at the parent plus the one of decVar (same value as the proof of getDualProofFromPseudoObjCutoff). */
static SCIP_Real pseudoObjAfterBranching(SCIP_Longint parentnumber, int decVar)
{
    SCIP_Real obj = 0;
    for (int idx : g_memoscip.nodeInfoMap[parentnumber].Var1)
        obj += SCIPvarGetObj(g_memoscip.AllVars[idx]);
    obj += SCIPvarGetObj(g_memoscip.AllVars[decVar]);
    return obj;
}

static
SCIP_DECL_EVENTEXEC(eventExecNodeSolved)
{
    ScopedTimer NodeSolvedTimer(g_memoscip.NodeSolvedTime);
    /* Event provides the node that was solved */
    SCIP_NODE* node = SCIPeventGetNode(event);
    if (node == nullptr)
        return SCIP_OKAY;
    SCIP_Longint nodenumber = SCIPnodeGetNumber(node);
    SCIP_Longint parentnumber = g_memoscip.nodeInfoMap[nodenumber].parent;

    if(g_options.verbosity > 1){
        cout<<"-----NODE SOLVED EVENT -------"<<endl;
        SCIP_STATUS status = SCIPgetStatus(scip);
        cout << "Node " << nodenumber<<" status "<<status<<" "<<g_memoscip.nodeInfoMap[nodenumber].Status<<endl;
    }

    if(g_memoscip.CurrentVarIdx.empty())
        initVarIndexMaps(scip);

    SCIP_Bool cutoff, lperror;
    SCIP_LPSOLSTAT stat = SCIPgetLPSolstat(scip);
    bool handled;

    // A child whose pseudo objective (objective of the variables fixed to 1) reaches the incumbent is cut off by SCIP
    // as soon as it is focused, without propagation nor LP: its proof is that pseudo objective, as for the
    // never-focused nodes handled in eventExecNodeDelete.
    const int decVar = (stat==SCIP_LPSOLSTAT_NOTSOLVED && nodenumber != 1) ? getBranchingVarFixedTo1(node) : -1;
    if(g_options.pseudoObjCut != 0 && decVar >= 0 && pseudoObjAfterBranching(parentnumber, decVar) >= SCIPgetPrimalbound(scip)){
        if(g_options.verbosity > 3)
            cout<<"Node "<<nodenumber<<" cut off by its pseudo objective, decision "<<SCIPvarGetName(g_memoscip.AllVars[decVar])<<endl;
        // The LP of a node cut off at focus is not constructed, but the proof needs its rows (weighted by the parent duals)
        if(!SCIPisLPConstructed(scip)){
            SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
        }
        SCIP_CALL( SCIPflushLP(scip) );
        SCIP_CALL( getDualProofFromPseudoObjCutoff(scip, nodenumber, parentnumber, decVar) );
    }else if(stat==SCIP_LPSOLSTAT_NOTSOLVED){
        // The node was decided before its LP was solved (e.g. infeasible by propagation): solve it in a dive to get a proof
        if(g_options.verbosity > 3){
            cout<<"Construct and solve dive LP"<<endl;
        }
        if(!SCIPisLPConstructed(scip)){
            SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
        }
        SCIP_CALL( SCIPstartDive(scip) );
        // Must keep in case problem decided infeasible (by bound propagation for example) before LP is constructed and solved.
        // Learned constraints not yet in the LP are relaxed in the dive.
        for(size_t i =0; i <g_memoscip.Unconstructed_cons.size(); i++){
            SCIP_CONS* cons = g_memoscip.Unconstructed_cons[i];
            SCIP_ROW* row = SCIPgetRowLinear(scip, cons);
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            if(row!=nullptr){
                SCIP_CALL( SCIPchgRowLhsDive(scip, row, -SCIPinfinity(scip)) );
                SCIP_CALL( SCIPchgRowRhsDive(scip, row, SCIPinfinity(scip)) );
            }
        }
        SCIP_CALL( SCIPflushLP(scip) );
        SCIP_CALL( SCIPsolveDiveLP(scip, -1, &lperror, &cutoff) );
        SCIP_LPSOLSTAT stat2 = SCIPgetLPSolstat(scip);

        SCIP_CALL( processNodeLP(scip, stat2, nodenumber, parentnumber, handled) );
        if(!handled){
            cout<<"lperror "<<lperror<<" cutoff "<<cutoff<<" stat2 "<<stat2<<endl;
            cout<<"LP SOLVED EVENT with status "<<stat2<<" node "<<nodenumber<<endl;
            char nameprobsol[50];
            snprintf(nameprobsol, 50, "infeasible%s_%lld.lp", g_options.instancePath.c_str(), nodenumber);
            SCIP_CALL( SCIPwriteTransProblem(scip, nameprobsol,"lp",FALSE) );
            SCIP_CALL( SCIPwriteLP(scip, nameprobsol) );
            assert(1==0);
        }

        SCIP_CALL( SCIPendDive(scip) );
    }else{
        SCIP_CALL( processNodeLP(scip, stat, nodenumber, parentnumber, handled) );
        if(!handled){
            cout<<"LP SOLVED EVENT with status "<<stat<<" node "<<nodenumber<<endl;
            assert(1==0);
        }
    }

    NodeInfo& ninfo = g_memoscip.nodeInfoMap[nodenumber];
    if(stat==SCIP_LPSOLSTAT_OPTIMAL){
        assert(ninfo.Status=="FEASIBLE");
    }

    if(ninfo.lb>=SCIPgetPrimalbound(scip))
        ninfo.LbgUb=true;

    if(g_options.verbosity > 2)
        cout<<"END OF NODE SOLVED EVENT  node "<<nodenumber<<" my lb "<<ninfo.lb<<" SCIP node lb "<<SCIPgetNodeLowerbound(scip, node)<<" ub "<<SCIPgetPrimalbound(scip)<< " diff between lb and scip node lb "<<fabs(SCIPgetNodeLowerbound(scip, node)-ninfo.lb)<<endl;
    SCIP_CALL( writeNodeLPFilesIfRequested(scip, nodenumber) );

    return SCIP_OKAY;
}

/* ---------------------------------------------------------------------------------------------------------------- */
/* LP solved                                                                                                         */
/* ---------------------------------------------------------------------------------------------------------------- */

static
SCIP_DECL_EVENTEXEC(eventExecLPSolved)
{
    ScopedTimer LPSolvedTimer(g_memoscip.LPSolvedTime);
    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    if( eventtype != SCIP_EVENTTYPE_LPSOLVED )
        return SCIP_OKAY;

    SCIP_NODE* node = SCIPeventGetNode(event);
    if(node == NULL){
        return SCIP_OKAY;
    }

    SCIP_Longint nodenumber= SCIPnodeGetNumber(node);
    g_memoscip.nodeInfoMap[nodenumber].LPsolved=true;
    g_memoscip.Unconstructed_cons.clear();
    if(g_options.verbosity > 1)
        cout<<"-----LP SOLVED EVENT ------- "<<nodenumber<< " depth "<<SCIPgetDepth(scip)<<" parent "<<g_memoscip.nodeInfoMap[nodenumber].parent<<endl;
    SCIP_LPSOLSTAT lpsolstat = SCIPgetLPSolstat(scip);
    if(g_options.verbosity >= 1){
        cout << "LP Solved - Status: ";
        switch(lpsolstat) {
            case SCIP_LPSOLSTAT_OPTIMAL:
                cout << "OPTIMAL "<<SCIPgetLPObjval(scip)<<" "<<SCIPgetPrimalbound(scip)<<" node status "<<g_memoscip.nodeInfoMap[nodenumber].Status<<endl;
                cout << " LP size rows=" << SCIPgetNLPRows(scip)
                        << " cols=" << SCIPgetNLPCols(scip)
                        << endl;
                break;
            case SCIP_LPSOLSTAT_INFEASIBLE:
                cout << "INFEASIBLE";
                cout << " LP size rows=" << SCIPgetNLPRows(scip)
                        << " cols=" << SCIPgetNLPCols(scip)
                        << endl;
                break;
            case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
                cout << "UNBOUNDED";
                break;
            case SCIP_LPSOLSTAT_OBJLIMIT:
                cout << "OBJECTIVE LIMIT REACHED";
                break;
            case SCIP_LPSOLSTAT_ITERLIMIT:
                cout << "ITERATION LIMIT";
                break;
            case SCIP_LPSOLSTAT_TIMELIMIT:
                cout << "TIME LIMIT";
                break;
            case SCIP_LPSOLSTAT_NOTSOLVED:
                cout << "NOT SOLVED";
                break;
            case SCIP_LPSOLSTAT_ERROR:
                cout << "ERROR";
                break;
            default:
                cout << "UNKNOWN (" << lpsolstat << ")";
                break;
        }
        cout << endl;
    }
    return SCIP_OKAY;
}

/* ---------------------------------------------------------------------------------------------------------------- */
/* Node focused                                                                                                      */
/* ---------------------------------------------------------------------------------------------------------------- */

static
SCIP_DECL_EVENTEXEC(eventExecNodeFocused)
{
    ScopedTimer NodeFocusedTimer(g_memoscip.NodeFocusedTime);
    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    if( eventtype != SCIP_EVENTTYPE_NODEFOCUSED )
        return SCIP_OKAY;

    SCIP_NODE* node = SCIPeventGetNode(event);
    if(node == NULL){
        return SCIP_OKAY;
    }

    SCIP_Longint nodenumber=SCIPnodeGetNumber(node);

    SCIP_NODE* parent = nullptr;
    SCIP_Longint parentNumber=0;
    NodeInfo &ninfo = g_memoscip.nodeInfoMap[nodenumber]; // creates if not existing
    ninfo.visited = true;
    ninfo.number = nodenumber;
    ninfo.ub = SCIPgetPrimalbound(scip);
    if(nodenumber!=1)  // skip root
    {
        parent = SCIPnodeGetParent(node);
        parentNumber = SCIPnodeGetNumber(parent);
        if(g_options.verbosity > 1)
            cout<<"-----------FOCUSING NODE "<<nodenumber<<" depth "<<SCIPgetDepth(scip)<<" with Parent "<<parentNumber<<"---------------------and "<<g_memoscip.nbAddedctr<<" CONSTRAINT added "<<endl<<"UB :"<< SCIPgetPrimalbound(scip)<<endl<<" LB : "<<SCIPnodeGetLowerbound(node)<<endl;
        ninfo.parent=parentNumber;
    }

    // The branching variable is the variable fixed to 1 at this node (the auxiliary variable of a variable fixed to 0)
    SCIP_DOMCHG* domchg = SCIPnodeGetDomchg(node);
    int nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
    for(int i = 0; i < nboundchgs; ++i )
    {
        SCIP_BOUNDCHG* boundchg = SCIPdomchgGetBoundchg(domchg, i);
        SCIP_VAR* var = SCIPboundchgGetVar(boundchg);
        SCIP_Real newbound = SCIPboundchgGetNewbound(boundchg);
        if(newbound==1.0)
            ninfo.BranchVar = shiftedIndex(var);
        else
            ninfo.BranchVar = g_memoscip.AuxVar[shiftedIndex(var)];
    }
    assert(g_memoscip.VarLb.size()==g_memoscip.VarUb.size());

    if(g_options.debug){
        for (size_t i = 0; i < g_memoscip.VarLb.size(); i++)
        {
            SCIP_VAR* var= g_memoscip.AllVars[i];
            assert(SCIPvarGetLbLocal(var)==g_memoscip.VarLb[i] || SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[shiftedIndex(var)]]])==0);
            assert(SCIPvarGetUbLocal(var)==g_memoscip.VarUb[i] ||  SCIPvarGetUbLocal(var)==SCIPinfinity(scip) || SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[shiftedIndex(var)]])==1);
            (void)var;
        }
    }
    int qsize=g_memoscip.CutQueue.queue.size();
    if(qsize>g_options.maxCutQueueSize && !g_options.sepa){
        SCIP_CALL(AddCutstoProb(scip,&g_memoscip.CutQueue,&g_memoscip,double(g_memoscip.TOP),g_options.sepa,&g_memoscip.nbctr,nodenumber,nullptr));
    }
    // link this node as child of its parent (if parent exists)
    if (parent != nullptr)
    {
        g_memoscip.nodeInfoMap[parentNumber].children.push_back(SCIPnodeGetNumber(node));
    }
    if(nodenumber!=1 && g_options.verbosity >= 1)
        cout<<"Depth "<<SCIPnodeGetDepth(node)<<" Node "<<nodenumber<<" lb: "<<g_memoscip.nodeInfoMap[parentNumber].lb<<" ub: "<<SCIPgetPrimalbound(scip)<<" Dec Var: "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[ninfo.BranchVar]])<<" nb added ctrs : "<<g_memoscip.nbAddedctr<<" nb nodes " <<SCIPgetNTotalNodes(scip)<<endl;
    if(g_options.verbosity > 1)
        cout<<"-------------------------------------ENDOFCUS EVENT ----------------"<<endl;
    if(g_options.verbosity > 3){
        printf("Total constraints: %d\n", SCIPgetNConss(scip));
        printf("Original constraints: %d\n", SCIPgetNOrigConss(scip));
        printf("Active constraints: %d\n", SCIPgetNActiveConss(scip));
    }
    return SCIP_OKAY;
}

/* ---------------------------------------------------------------------------------------------------------------- */
/* Node deleted: bottom-up aggregation of the children proofs                                                        */
/* ---------------------------------------------------------------------------------------------------------------- */

/** Removes from coefs the common part of the coefficient of literal i and of its complementary literal(s)
 *  (auxiliary variable, then tuple variables in encoding 1), decreasing rhs accordingly.
 *  If checkVarTuple, the tuple reduction is only done when literal i belongs to a tuple constraint. */
static void reduceComplementaryCoefs(SCIP* scip, vector<SCIP_Real>& coefs, size_t i, SCIP_Real& rhs, bool checkVarTuple)
{
    const int aux = g_memoscip.AuxVar[i];
    if((int)i!= aux && !SCIPisDualfeasZero(scip,coefs[aux])){
        if(coefs[i]!= 0 && coefs[aux]!=0){
            SCIP_Real minval = min(coefs[i],coefs[aux]);
            rhs -= minval;
            coefs[i]-=minval;
            coefs[aux]-=minval;
            assert(coefs[i]==0 || coefs[aux]==0);
        }
    }

    if(g_options.encoding==1 && (checkVarTuple ? !g_memoscip.TupleAuxVar[i].empty() : !g_memoscip.TupleAuxVar.empty())){
        const vector<int>& tuple = g_memoscip.TupleAuxVar[i];
        SCIP_Real minval = coefs[i];
        for (int j = 0; j < (int)tuple.size(); j++)
            minval = min(minval, coefs[tuple[j]]);
        rhs -= minval;
        coefs[i]-=minval;
        coefs[tuple[0]]-=minval;
        coefs[tuple[1]]-=minval;
        coefs[tuple[2]]-=minval;
    }
}

static
SCIP_DECL_EVENTEXEC(eventExecNodeDelete)
{
    ScopedTimer NodeDeleteTimer(g_memoscip.NodeDeleteTime);
    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);

    if( eventtype != SCIP_EVENTTYPE_NODEDELETE )
        return SCIP_OKAY;

    SCIP_NODE* node = SCIPeventGetNode(event);
    if(node == NULL){
        return SCIP_OKAY;
    }
    SCIP_Longint nodenumber=SCIPnodeGetNumber(node);
    if(nodenumber==1) return SCIP_OKAY; // skip root
    SCIP_NODE* parent = SCIPnodeGetParent(node);
    NodeInfo& nodeInfo = g_memoscip.nodeInfoMap[nodenumber];
    SCIP_Longint parentNumber = nodeInfo.parent;
    if(g_options.verbosity > 1){
        cout<<"Deleting node "<<nodenumber<<" from nodeInfoMap, scip Node Lb "<< SCIPgetNodeLowerbound(scip,node)  <<" my node lb : "<<nodeInfo.lb<< " "<<" depth : "<<SCIPgetDepth(scip)<< " parent: "<<parentNumber<<endl;
        cout<<"UB "<<SCIPgetPrimalbound(scip)<<endl;
        cout<<"number of children "<<nodeInfo.children.size()<<endl;
        cout<<"Bracnh var "<<nodeInfo.BranchVar<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[nodeInfo.BranchVar]])<<endl;
        cout<<"Curnnet node lb increase "<<nodeInfo.lb_increase<<" versus children "<<endl;
    }

    //Case node was not processed. Mostly because it was pruned immediately after creation, since a brother node found a better solution than the pseudo obj of this node.
    if(parentNumber==-1){
        SCIP_Longint SCIP_parentNumber = SCIPnodeGetNumber(parent);
        if(SCIPisLT(scip, ceil(g_memoscip.nodeInfoMap[SCIP_parentNumber].lb), SCIPgetPrimalbound(scip))){
            int nvars;
            // First query the required size with branchvarssize = 0
            SCIPnodeGetParentBranchings(node, NULL, NULL, NULL, &nvars, 0);
            assert(nvars==1);
            int decvar;
            if(g_options.verbosity > 3)
                cout<<"case node cutted of parentNumber "<<SCIP_parentNumber<<endl;

            if (nvars > 0)
            {
                std::vector<SCIP_VAR*> vars(nvars);
                std::vector<SCIP_Real> bounds(nvars);
                std::vector<SCIP_BOUNDTYPE> boundtypes(nvars);
                SCIPnodeGetParentBranchings(node, vars.data(), bounds.data(),
                                                boundtypes.data(), &nvars, nvars);

                for (int i = 0; i < nvars; ++i)
                {
                    decvar = allVarsIndex(vars[i]);
                    assert(SCIPvarGetName(vars[i]) == SCIPvarGetName(g_memoscip.AllVars[decvar]));
                    assert(bounds[i]==1);
                    nodeInfo.BranchVar= shiftedIndex(vars[i]);
                }
            }

            assert(nodeInfo.lb==0);
            assert(SCIP_parentNumber>0);
            SCIP_CALL(getDualProofFromPseudoObjCutoff(scip,nodenumber, SCIP_parentNumber, decvar));
            g_memoscip.nodeInfoMap[SCIP_parentNumber].children.push_back(nodenumber);
            return SCIP_OKAY;
        }
    }

    vector<long> children = nodeInfo.children;
    if (g_memoscip.nbAddedctr<(int)g_options.nbslack && (nodeInfo.lb_increase != 0 || children.size()==2) ) {

        vector<SCIP_Real> ResFus(g_memoscip.AllVars.size(),0.0);
        SCIP_Real rhsFus=0;
        if(children.size()==1){
            checkSingleChildProof(nodenumber, children[0]);
        }
        else if(children.size()==2){
            if(g_options.verbosity > 2){
                for(long child : children)
                    cout<<"Two children "<<child<<" lb "<<g_memoscip.nodeInfoMap[child].lb<<" versus parent lb "<<nodeInfo.lb<<" "<<g_memoscip.nodeInfoMap[child].DualProofCoef.first.size()<<" "<<g_memoscip.nodeInfoMap[child].DualProofCoef.second<<" incr lb "<<g_memoscip.nodeInfoMap[child].lb_increase<<endl;
            }
            SCIP_CALL( fuseChildrenProofs(scip, nodenumber, children, ResFus, rhsFus) );
        }else{
            if(g_options.verbosity > 2)
                cout<<"No child "<<endl;
        }
        assert(ResFus.size()==g_memoscip.AllVars.size());

        if(g_options.verbosity > 2){
            cout<<"ResFus "<<endl;
            printProof(scip, ResFus, ">=", rhsFus, false);
        }
        SCIP_Real rhs_summed=rhsFus;
        vector<int> NewInds;
        vector<SCIP_Real> NewCoefs;
        vector<SCIP_VAR*> NewVar;

        // Learning level 2: learn the fused proof of the children alone
        if(rhs_summed>0 && g_options.learning==2 && rhs_summed < g_memoscip.TOP){
            char name[100];
            snprintf(name, 100, "l2cut_%lld_%d", nodenumber,SCIPgetDepth(scip));
            SCIP_Real rhs_ctr=rhs_summed;
            vector<SCIP_Real> ResFus2=ResFus;

            // Reduce first, then collect: a tuple reduction can raise the coefficient of an already visited variable
            for (size_t i = 0; i < ResFus2.size(); i++)
                if(!SCIPisDualfeasZero(scip,ResFus2[i]))
                    reduceComplementaryCoefs(scip, ResFus2, i, rhs_ctr, true);
            for (size_t i = 0; i < ResFus2.size(); i++)
            {
                if(!SCIPisDualfeasZero(scip,ResFus2[i])){
                    NewInds.push_back(i);
                    NewCoefs.push_back(ResFus2[i]);
                    NewVar.push_back(g_memoscip.AllVars[i]);
                }
            }

            if(g_options.verbosity > 2)
                cout<<"NODE DELETE Adding cut "<<name<<" with rhs "<<rhs_ctr<<" and "<<NewVar.size()<< " variables vs Original vars "<< g_memoscip.nOriginalVar<< endl;

            assert(NewVar.size()==NewCoefs.size());
            assert(NewVar.size()!=0);
            SCIP_CALL( SCIPqueueCut(scip, name,
                                (int)NewVar.size(),
                                NewVar.data(),
                                NewCoefs.data(),
                                rhs_ctr,  // lhs (equality constraint)
                                rhs_ctr,&g_memoscip.CutQueue) );

            NewInds.clear();
            NewCoefs.clear();
            NewVar.clear();
        }

        assert(fabs(nodeInfo.DualProofCoef.second-nodeInfo.lb_increase)<1e-5);

        addNodeProof(nodenumber, ResFus, rhs_summed);
        applyBoundPropagationRemovals(nodenumber, ResFus, rhs_summed);
        applyNodeConsistency(scip, nodenumber, parentNumber, ResFus, rhs_summed);

        if(g_options.verbosity > 2){
            cout<<"final res "<<endl;
            printProof(scip, ResFus, ">=", rhs_summed, true);
        }

        checkSummedRhs(nodenumber, children, rhs_summed);

        // The proof of the node (before reduction) is kept to be fused at its parent
        vector<int> Inds;
        vector<SCIP_Real> Coefs;
        for (size_t i = 0; i < ResFus.size(); i++)
        {
            if(!SCIPisDualfeasZero(scip,ResFus[i])){
                Inds.push_back(i);
                Coefs.push_back(ResFus[i]);
                assert( i >= g_memoscip.nodeInfoMap[parentNumber].redcosts.size() || ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<1e-3 || g_memoscip.VarUb[i] == 0);
            }
        }
        if(rhs_summed>0  && rhs_summed < g_memoscip.TOP){
            char name[100];
            snprintf(name, 100, "cut_%lld_%d", nodenumber,SCIPgetDepth(scip));
            SCIP_Real rhs_ctr=rhs_summed;
            // Reduce first, then collect: a tuple reduction can raise the coefficient of an already visited variable
            for (size_t i = 0; i < ResFus.size(); i++)
            {
                if(!SCIPisDualfeasZero(scip,ResFus[i])){
                    assert(g_options.encoding==0 || g_memoscip.TupleAuxVar[i].empty() || g_memoscip.TupleAuxVar[i].size() == 3);
                    reduceComplementaryCoefs(scip, ResFus, i, rhs_ctr, true);
                }
            }
            for (size_t i = 0; i < ResFus.size(); i++)
            {
                if(!SCIPisDualfeasZero(scip,ResFus[i])){
                    if(i<g_memoscip.nodeInfoMap[parentNumber].redcosts.size() && g_options.verbosity > 4){
                      cout<<SCIPvarGetName(g_memoscip.AllVars[i])<< " i "<< i<<" resfus "<<ResFus[i]<<" parentrc "<<g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<" diff "<<ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<" Bounds "<<g_memoscip.VarUb[i]<<" "<< g_memoscip.VarLb[i]<<endl;
                    }
                    NewInds.push_back(i);
                    NewCoefs.push_back(ResFus[i]);
                    NewVar.push_back(g_memoscip.AllVars[i]);
                }
            }

            if(g_options.verbosity > 2)
                cout<<"NODE DELETE Adding cut "<<name<<" with rhs "<<rhs_ctr<<" and "<<NewVar.size()<< " variables vs Original vars "<< g_memoscip.nOriginalVar<< endl;

            assert(NewVar.size()==NewCoefs.size());
            assert(NewVar.size()!=0);
            // Learning level 4 learns at every node; otherwise only at nodes with two children whose proof is useful
            if(g_options.learning==4 || (children.size()==2 && (g_options.learning==3 || nodeInfo.lb_increase>0 || (g_memoscip.nodeInfoMap[parentNumber].lb_increase==0 && (g_memoscip.nodeInfoMap[children[0]].lb_increase>0 || g_memoscip.nodeInfoMap[children[1]].lb_increase>0 ))) ))
                SCIP_CALL( SCIPqueueCut(scip, name,
                                (int)NewVar.size(),
                                NewVar.data(),
                                NewCoefs.data(),
                                rhs_ctr,  // lhs (equality constraint)
                                rhs_ctr,&g_memoscip.CutQueue) );
            if(g_options.verbosity > 2)
                SCIPprintQueuedCuts(scip, &g_memoscip.CutQueue);
        }

        assert(abs(nodeInfo.DualProofCoef.second-nodeInfo.lb_increase) < 1e-9);
        assert(NewCoefs.size()<=g_memoscip.CurrentVarIdx.size());
        nodeInfo.DualProofCoef.second=rhs_summed;
        nodeInfo.DualProofInds=Inds;
        nodeInfo.DualProofCoef.first=Coefs;
        assert(nodeInfo.lb_increase<=rhs_summed || !nodeInfo.NodeConsistency.empty());
        nodeInfo.lb_increase=rhs_summed;
        assert(Myceil(g_memoscip.nodeInfoMap[parentNumber].lb + nodeInfo.lb_increase) >= SCIPgetPrimalbound(scip)-1 ||SCIPgetPrimalbound(scip)==SCIPinfinity(scip) );
        eraseChildren(nodenumber, children);
    }
    return SCIP_OKAY;
}

/* ---------------------------------------------------------------------------------------------------------------- */
/* Bound changed: explanations of bound propagation and reduced cost filtering                                       */
/* ---------------------------------------------------------------------------------------------------------------- */

/** Walks the lb/ub change stack of an active variable (a single SCIP_VAR representation). */
static SCIP_BDCHGINFO* searchBdchgStackForEvent(
    SCIP* scip,
    SCIP_VAR* var,
    SCIP_BOUNDTYPE boundtype,
    SCIP_Real eventOldbound,
    SCIP_Real eventNewbound)
{
    if (var == NULL || !SCIPvarIsActive(var))
        return NULL;

    const int n = (boundtype == SCIP_BOUNDTYPE_LOWER)
        ? SCIPvarGetNBdchgInfosLb(var)
        : SCIPvarGetNBdchgInfosUb(var);
    if (n <= 0)
        return NULL;

    auto getAtPos = [&](int pos) -> SCIP_BDCHGINFO* {
        return (boundtype == SCIP_BOUNDTYPE_LOWER)
            ? SCIPvarGetBdchgInfoLb(var, pos)
            : SCIPvarGetBdchgInfoUb(var, pos);
    };

    auto pairMatches = [&](SCIP_BDCHGINFO* info) -> bool {
        return SCIPisFeasEQ(scip, SCIPbdchginfoGetOldbound(info), eventOldbound)
            && SCIPisFeasEQ(scip, SCIPbdchginfoGetNewbound(info), eventNewbound);
    };
    auto newMatches = [&](SCIP_BDCHGINFO* info) -> bool {
        return SCIPisFeasEQ(scip, SCIPbdchginfoGetNewbound(info), eventNewbound);
    };

    /* 1) Same (old, new) bound pair as the event, non-redundant. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (pairMatches(info) && !SCIPbdchginfoIsRedundant(info))
            return info;
    }
    /* 2) Same pair, even if redundant. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (pairMatches(info))
            return info;
    }
    /* 3) Same new bound only (previous behaviour), non-redundant. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (newMatches(info) && !SCIPbdchginfoIsRedundant(info))
            return info;
    }
    /* 4) Same new bound only, redundant allowed. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (newMatches(info))
            return info;
    }
    /* 5) Last non-redundant entry. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (!SCIPbdchginfoIsRedundant(info))
            return info;
    }
    return getAtPos(n - 1);
}

/** Finds the SCIP_BDCHGINFO matching the event; tries the event variable, then SCIPvarGetProbvar (negation/aggregation chain). */
SCIP_BDCHGINFO* getBoundChangeReasonDetailed(SCIP* scip, SCIP_VAR* var,
    SCIP_BOUNDTYPE boundtype,
    SCIP_Real oldbound,
    SCIP_Real newbound)
{
    SCIP_VAR* prob = SCIPvarGetProbvar(var);

    SCIP_BDCHGINFO* info = searchBdchgStackForEvent(scip, var, boundtype, oldbound, newbound);
    if (info != NULL)
        return info;
    if (prob != NULL && prob != var)
        info = searchBdchgStackForEvent(scip, prob, boundtype, oldbound, newbound);
    return info;
}

/** Asserts that the stored bounds of v match its local bounds, unless its auxiliary variable is fixed. */
static void assertStoredBounds(SCIP_VAR* v)
{
    SCIP_VAR* auxv = g_memoscip.AllVars[g_memoscip.AuxVar[shiftedIndex(v)]];
    assert(g_memoscip.VarLb[allVarsIndex(v)]==SCIPvarGetLbLocal(v) || (SCIPvarGetLbLocal(auxv)==SCIPvarGetUbLocal(auxv)));
    assert(g_memoscip.VarUb[allVarsIndex(v)]==SCIPvarGetUbLocal(v) || (SCIPvarGetLbLocal(auxv)==SCIPvarGetUbLocal(auxv)));
    (void)auxv;
}

static void printConstraintReason(SCIP_CONS* cons)
{
    if(g_options.verbosity > 2){
        std::string constype = SCIPconshdlrGetName(SCIPconsGetHdlr(cons));
        std::string consname = SCIPconsGetName(cons);
        cout<< "CONSTRAINT (" + constype + ": " + consname + ")"<<endl;
    }
}

/** Explanation of a lower bound change of var inferred by a linear constraint (vars, coefs):
 *  AllVars indices of the literals assigned to 0 that forced it. */
static vector<int> explainLbChange(SCIP_VAR* var, SCIP_VAR** vars, SCIP_Real* coefs, int nvars)
{
    vector<int> explanation;
    bool varsign;
    for(int i=0; i < nvars; i++) {
        if(vars[i] == var){
            varsign = coefs[i] > 0;
            break;
        }
    }
    for(int i = 0; i < nvars; i++) {
        SCIP_VAR* v = vars[i];
        if(v == var)
            continue;
        if(g_options.verbosity > 3 )
            cout<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
        assertStoredBounds(v);
        const bool sameSign = varsign ? coefs[i]>1e-7 : coefs[i]<-1e-7;
        const bool oppositeSign = varsign ? coefs[i]<-1e-7 : coefs[i]>1e-7;
        // A variable fixed to 0 is its own explanation, a variable fixed to 1 is explained by its auxiliary variable(s)
        if(sameSign && g_memoscip.VarUb[allVarsIndex(v)]==0)
            explanation.push_back(allVarsIndex(v));
        if(oppositeSign && g_memoscip.VarLb[allVarsIndex(v)]==1)
            forEachAuxVar(shiftedIndex(v), true, [&](int aux){ explanation.push_back(aux); });
    }
    return explanation;
}

/** Explanation of an upper bound change of var inferred by a linear constraint (vars, coefs):
 *  AllVars indices of the literals assigned to 0 that forced it (mirror of explainLbChange). */
static vector<int> explainUbChange(SCIP_VAR* var, SCIP_VAR** vars, SCIP_Real* coefs, int nvars)
{
    vector<int> explanation;
    bool varsign;
    for(int i=0; i < nvars; i++) {
        if(vars[i] == var){
            varsign = coefs[i] > 0;
            break;
        }
    }
    for(int i = 0; i < nvars; i++) {
        SCIP_VAR* v = vars[i];
        if(v == var)
            continue;
        assertStoredBounds(v);
        const bool sameSign = varsign ? coefs[i]>1e-7 : coefs[i]<-1e-7;
        const bool oppositeSign = varsign ? coefs[i]<-1e-7 : coefs[i]>1e-7;
        // A variable of opposite sign fixed to 0 is its own explanation,
        // a variable of same sign fixed to 1 is explained by its auxiliary variable(s)
        if(oppositeSign && g_memoscip.VarUb[allVarsIndex(v)]==0){
            explanation.push_back(allVarsIndex(v));
            if(g_options.verbosity > 3 )
                cout<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
        }
        if(sameSign && g_memoscip.VarLb[allVarsIndex(v)]==1){
            forEachAuxVar(shiftedIndex(v), false, [&](int aux){ explanation.push_back(aux); });
            if(g_options.verbosity > 3 )
                cout<<"Auxvar "<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
        }
    }
    return explanation;
}

/** Returns true if prop is reduced cost filtering (used for node consistency). */
static bool isReducedCostFiltering(SCIP_PROP* prop)
{
    std::string propname = SCIPpropGetName(prop);
    return propname == "redcost" || propname == "rootredcost";
}

static void printReducedCostFiltering(SCIP* scip, SCIP_NODE* node, SCIP_Longint ParentNumber, SCIP_VAR* var, SCIP_BDCHGINFO* reason)
{
    if(g_options.verbosity > 2){
        std::string propname = SCIPpropGetName(SCIPbdchginfoGetInferProp(reason));
        cout<<"PROPAGATION (" + propname + ", info=" + std::to_string(SCIPbdchginfoGetInferInfo(reason)) + "): reduced-cost filtering"<<endl;
        cout<<"Parent node "<<ParentNumber<<" lb "<<g_memoscip.nodeInfoMap[ParentNumber].lb<<" redcost "<<g_memoscip.nodeInfoMap[ParentNumber].redcosts[allVarsIndex(var)]<<" new lb "<<SCIPnodeGetLowerbound(node)<<" "<<SCIPgetPrimalbound(scip)<<endl;
    }
}

static void handleLbChanged(SCIP* scip, SCIP_NODE* node, SCIP_Longint nodenumber, SCIP_VAR* var, SCIP_Real oldbound, SCIP_Real newbound)
{
    if(g_options.verbosity > 1){
        cout<<"LB CHANGED "<<SCIPvarGetName(var)<<" "<<oldbound<<" ->"<<newbound<<endl;
    }
    assert(!isSlackVar(var));
    const int pos = shiftedIndex(var);
    g_memoscip.VarLb[allVarsIndex(var)]=newbound;
    bool auxvarassigned=false;
    if(g_options.encoding==0 || g_memoscip.TupleAuxVar[pos].empty()){
        // The auxiliary variable is the complement of var
        const int aux = g_memoscip.AuxVar[pos];
        auxvarassigned =  g_memoscip.VarUb[aux] == g_memoscip.VarLb[aux];
        if(newbound==1)
            g_memoscip.VarUb[aux]=0;
        else if(newbound==0)
            g_memoscip.VarUb[aux]=1;
        else
            g_memoscip.VarUb[aux]=newbound;
    }

    SCIP_BDCHGINFO* reason = getBoundChangeReasonDetailed(scip, var, SCIP_BOUNDTYPE_LOWER, oldbound, newbound);

    assert(!SCIPinProbing(scip));
    if(reason!=NULL && !auxvarassigned){
        if(g_options.verbosity > 2)
            cout<<"reason found LB Change "<<SCIPvarGetName(var)<<" "<<oldbound<<" "<<newbound<<endl;
        switch (SCIPbdchginfoGetChgtype(reason)) {
            case SCIP_BOUNDCHGTYPE_CONSINFER: {
                SCIP_CONS* cons = SCIPbdchginfoGetInferCons(reason);
                if (cons != NULL) {
                    printConstraintReason(cons);
                    int nvars = SCIPgetNVarsLinear(scip,cons);
                    if(nvars>2){
                        vector<int> explanation = explainLbChange(var, SCIPgetVarsLinear(scip,cons), SCIPgetValsLinear(scip,cons), nvars);
                        if(g_options.verbosity > 3)
                            cout<<"Explanation LB for "<<SCIPvarGetName(var)<<" "<<explanation.size()<<endl;
                        // The auxiliary variable(s) of var are removed by bound propagation
                        NodeInfo& ninfo = g_memoscip.nodeInfoMap[nodenumber];
                        forEachAuxVar(pos, true, [&](int aux){ ninfo.BPRemoval.push_back({aux, explanation}); });
                    }
                }
                break;
            }
            case SCIP_BOUNDCHGTYPE_PROPINFER: {
                SCIP_PROP* prop = SCIPbdchginfoGetInferProp(reason);
                if (prop != NULL && isReducedCostFiltering(prop)) {
                    NodeInfo& ninfo = g_memoscip.nodeInfoMap[nodenumber];
                    SCIP_Longint ParentNumber=ninfo.parent;
                    printReducedCostFiltering(scip, node, ParentNumber, var, reason);
                    // The auxiliary variable(s) of var are removed by reduced cost filtering
                    forEachAuxVar(pos, true, [&](int aux){
                        ninfo.NodeConsistency.push_back({aux,SCIPgetPrimalbound(scip)-g_memoscip.nodeInfoMap[ParentNumber].lb});
                    });
                    assert(!isSlackVar(var));
                }
                break;
            }
            default:
                break;
        }
    }
    if(g_options.verbosity > 1)
        cout << "[" << reason << "] "<< SCIPvarGetName(var) << " LB: "<< oldbound << " → " << newbound << endl;
}

static void handleUbChanged(SCIP* scip, SCIP_NODE* node, SCIP_Longint nodenumber, SCIP_VAR* var, SCIP_Real oldbound, SCIP_Real newbound)
{
    if(g_options.verbosity > 1)
        cout<<"UB CHANGED "<<SCIPvarGetName(var)<<" "<<oldbound<<" ->"<<newbound<<endl;

    const int pos = shiftedIndex(var);
    bool auxvarassigned=false;
    if(g_memoscip.AuxVar[pos] != -1)
        auxvarassigned =  g_memoscip.VarUb[g_memoscip.AuxVar[pos]] == g_memoscip.VarLb[g_memoscip.AuxVar[pos]];
    SCIP_BDCHGINFO* reason = getBoundChangeReasonDetailed(scip, var, SCIP_BOUNDTYPE_UPPER, oldbound, newbound);
    if(reason!=NULL && !auxvarassigned){
        switch (SCIPbdchginfoGetChgtype(reason)) {
            case SCIP_BOUNDCHGTYPE_CONSINFER: {
                SCIP_CONS* cons = SCIPbdchginfoGetInferCons(reason);
                if (cons != NULL) {
                    printConstraintReason(cons);
                    int nvars = SCIPgetNVarsLinear(scip,cons);
                    if(nvars>2){
                        vector<int> explanation = explainUbChange(var, SCIPgetVarsLinear(scip,cons), SCIPgetValsLinear(scip,cons), nvars);
                        NodeInfo& ninfo = g_memoscip.nodeInfoMap[nodenumber];
                        if(g_options.verbosity > 3)
                            cout<<"ExplanationUB for "<<SCIPvarGetName(var)<<" "<<explanation.size()<<" "<<ninfo.BPRemoval.size()<<endl;
                        // var itself is removed by bound propagation
                        ninfo.BPRemoval.push_back({allVarsIndex(var), explanation});
                    }
                }
                break;
            }
            case SCIP_BOUNDCHGTYPE_PROPINFER: {
                SCIP_PROP* prop = SCIPbdchginfoGetInferProp(reason);
                if (prop == NULL) {
                    cout<<"PROPAGATION (unknown, info=" + std::to_string(SCIPbdchginfoGetInferInfo(reason)) + ")"<<endl;
                }
                else if (isReducedCostFiltering(prop)) {
                    NodeInfo& ninfo = g_memoscip.nodeInfoMap[nodenumber];
                    SCIP_Longint ParentNumber=ninfo.parent;
                    printReducedCostFiltering(scip, node, ParentNumber, var, reason);
                    assert(!isSlackVar(var));
                    // var itself is removed by reduced cost filtering
                    ninfo.NodeConsistency.push_back({allVarsIndex(var),SCIPgetPrimalbound(scip)-g_memoscip.nodeInfoMap[ParentNumber].lb});
                }
                break;
            }
            default:
                break;
        }
    }
    if(g_options.verbosity > 1 )
        cout << "[" << reason << "] "<< SCIPvarGetName(var) << " UB: "<< oldbound << " → " << newbound << endl;

    const int idx = allVarsIndex(var);
    if(isSlackVar(var)){
        g_memoscip.VarUb[idx]=newbound;
    }else if(newbound==SCIPinfinity(scip) || newbound==1){
        g_memoscip.VarUb[idx]=1;
        if(g_options.encoding==0 || g_memoscip.TupleAuxVar[pos].empty())
            g_memoscip.VarLb[g_memoscip.AuxVar[pos]]=0;
    }else{
        g_memoscip.VarUb[idx]=newbound;
        assert(newbound==0);
        if(g_options.encoding==0 || g_memoscip.TupleAuxVar[pos].empty())
            g_memoscip.VarLb[g_memoscip.AuxVar[pos]]=1;
    }
}

static SCIP_DECL_EVENTEXEC(eventExecBoundChanged)
{
    ScopedTimer BoundChangedTimer(g_memoscip.BoundChangedTime);
    SCIP_VAR* var = SCIPeventGetVar(event);
    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    SCIP_NODE* node = SCIPgetCurrentNode(scip);
    if(node == NULL || g_memoscip.CurrentVarIdx.empty()){
        return SCIP_OKAY;
    }
    SCIP_Longint nodenumber=SCIPnodeGetNumber(node);

    if(g_options.verbosity > 1)
        cout<<"-----BOUND CHANGED EVENT ------- "<< nodenumber<<" "<<" "<<SCIPvarGetName(var)<<" "<<SCIPgetDepth(scip)<<endl;

    SCIP_Real oldbound = SCIPeventGetOldbound(event);
    SCIP_Real newbound = SCIPeventGetNewbound(event);

    if (eventtype & SCIP_EVENTTYPE_LBCHANGED)
        handleLbChanged(scip, node, nodenumber, var, oldbound, newbound);

    if (eventtype & SCIP_EVENTTYPE_UBCHANGED)
        handleUbChanged(scip, node, nodenumber, var, oldbound, newbound);

    return SCIP_OKAY;
}

// Init callback: attach the bound change events to every variable
static SCIP_DECL_EVENTINIT(eventInitBoundChanged)
{
    SCIP_VAR** vars = SCIPgetVars(scip);
    int nvars = SCIPgetNVars(scip);
    g_boundChangedEventVars.clear();

    for (int i = 0; i < nvars; i++) {
        // Variable events must be caught with SCIPcatchVarEvent (not SCIPcatchEvent)
        SCIP_CALL( SCIPcatchVarEvent(scip, vars[i],
                                     SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
                                     eventhdlr,
                                     NULL,    // eventdata
                                     NULL) ); // filterpos
        g_boundChangedEventVars.push_back(vars[i]);
    }

    return SCIP_OKAY;
}

static SCIP_DECL_EVENTEXIT(eventExitBoundChanged)
{
    for (SCIP_VAR* var : g_boundChangedEventVars) {
        // Variable events must be dropped with SCIPdropVarEvent (not SCIPdropEvent)
        SCIP_CALL( SCIPdropVarEvent(scip, var,
                                    SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
                                    eventhdlr,
                                    NULL,    // eventdata
                                    -1) );   // filterpos (-1 = auto)
    }
    g_boundChangedEventVars.clear();

    return SCIP_OKAY;
}

/* ---------------------------------------------------------------------------------------------------------------- */
/* Other events                                                                                                      */
/* ---------------------------------------------------------------------------------------------------------------- */

static
SCIP_DECL_EVENTEXEC(eventExecNodeInfeasible)
{
    // Nothing to do: infeasible nodes are handled when the node is solved (eventExecNodeSolved)
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXEC(eventExecBestSolFound)
{
    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    if( eventtype != SCIP_EVENTTYPE_SOLFOUND )
        return SCIP_OKAY;

    SCIP_SOL* sol = SCIPeventGetSol(event);
    if( sol == nullptr )
        return SCIP_OKAY;

    if( g_options.verbosity > 0 )
    {
        cout << "-----INTEGER SOLUTION FOUND EVENT -------" << endl;
        cout << "node=" << SCIPgetNNodes(scip)
             << " obj=" << SCIPgetSolOrigObj(scip, sol)
             << " nsols=" << SCIPgetNSols(scip)
             << " primalbound=" << SCIPgetPrimalbound(scip)
             << endl;
    }

    return SCIP_OKAY;
}

/* ---------------------------------------------------------------------------------------------------------------- */
/* Branching rule                                                                                                    */
/* ---------------------------------------------------------------------------------------------------------------- */

/** Replays the branching decisions read from kReplayBranchingFile (option -b=4). */
static SCIP_RETCODE branchReplay(SCIP* scip, SCIP_RESULT* result)
{
    if (!ensure_replay_loaded())
    {
        if (g_options.verbosity > 0)
            cout << "Replay branching sequence not loaded from " << kReplayBranchingFile << endl;
        *result = SCIP_DIDNOTFIND;
        return SCIP_OKAY;
    }

    while (g_replayDecisionPos < g_replayDecisions.size())
    {
        const ReplayDecision& step = g_replayDecisions[g_replayDecisionPos];
        SCIP_VAR* v = find_replay_var(scip, step.varname);
        if (v == nullptr)
        {
            if (g_options.verbosity > 2)
                cout << "Replay skip: unknown variable " << step.varname << endl;
            ++g_replayDecisionPos;
            continue;
        }

        if (SCIPvarGetType(v) == SCIP_VARTYPE_CONTINUOUS)
        {
            ++g_replayDecisionPos;
            continue;
        }

        const SCIP_Real lb = SCIPvarGetLbLocal(v);
        const SCIP_Real ub = SCIPvarGetUbLocal(v);
        const SCIP_Real target = static_cast<SCIP_Real>(step.value);

        if (target < lb - 1e-9 || target > ub + 1e-9)
        {
            if (g_options.verbosity > 2)
                cout << "Replay skip: incompatible bounds for " << SCIPvarGetName(v) << " target=" << step.value << " lb=" << lb << " ub=" << ub << endl;
            ++g_replayDecisionPos;
            continue;
        }

        SCIP_Real x = SCIPgetSolVal(scip, nullptr, v);
        if (SCIPisFeasIntegral(scip, x) && SCIPisFeasEQ(scip, x, target))
        {
            ++g_replayDecisionPos;
            continue;
        }

        SCIP_NODE* downchild = nullptr;
        SCIP_NODE* eqchild   = nullptr;
        SCIP_NODE* upchild   = nullptr;

        SCIP_CALL( SCIPbranchVarVal(scip, v, target, &downchild, &eqchild, &upchild) );
        *result = SCIP_BRANCHED;

        if (g_options.verbosity > 0)
        {
            cout << "Replay branching step " << (g_replayDecisionPos + 1) << "/" << g_replayDecisions.size()
                 << " var=" << SCIPvarGetName(v) << " target=" << step.value << " sol=" << x << endl;
        }

        ++g_replayDecisionPos;
        return SCIP_OKAY;
    }

    if (g_options.verbosity > 0)
        cout << "Replay branching sequence exhausted" << endl;

    *result = SCIP_DIDNOTFIND;
    return SCIP_OKAY;
}

static
SCIP_DECL_BRANCHEXECLP(branchExecLP)
{
    if(g_options.verbosity > 1)
        cout<<"branching on fractional variables"<<endl;

    // If domains were tightened, let SCIP re-solve the LP first
    if( *result == SCIP_REDUCEDDOM )
        return SCIP_OKAY;

    if (g_options.branchingHeuristic == 4)
        return branchReplay(scip, result);

    if(g_options.branchingHeuristic < 1){
        // Branch on the first fractional variable (in g_memoscip.OrderVar order)
        for( SCIP_VAR* v : g_memoscip.OrderVar )
        {
            // Only branch on discrete vars
            if( SCIPvarGetType(v) == SCIP_VARTYPE_CONTINUOUS )
                continue;

            SCIP_Real x = SCIPgetSolVal(scip, nullptr, v);  // current LP solution value
            if( !SCIPisFeasIntegral(scip, x) )
            {
                SCIP_NODE* downchild = nullptr;
                SCIP_NODE* eqchild   = nullptr;
                SCIP_NODE* upchild   = nullptr;

                SCIP_CALL( SCIPbranchVarVal(scip, v, 1, &downchild, &eqchild, &upchild) );
                *result = SCIP_BRANCHED;

                if(g_options.verbosity > 1)
                    cout << "Branching on first fractional var " << SCIPvarGetName(v) << " x=" << x << " probindex=" << SCIPvarGetProbindex(v) << endl;
                return SCIP_OKAY;
            }
        }
    }

    // Nothing fractional found -> nothing to do
    *result = SCIP_DIDNOTFIND;
    return SCIP_OKAY;
}

/* ---------------------------------------------------------------------------------------------------------------- */
/* Registration                                                                                                      */
/* ---------------------------------------------------------------------------------------------------------------- */

/** Init / exit callbacks catching / dropping the global event type EventType. */
template <SCIP_EVENTTYPE EventType>
static SCIP_DECL_EVENTINIT(eventInitCatch)
{
    SCIP_CALL(SCIPcatchEvent(scip, EventType, eventhdlr, nullptr, nullptr));
    return SCIP_OKAY;
}

template <SCIP_EVENTTYPE EventType>
static SCIP_DECL_EVENTEXIT(eventExitDrop)
{
    SCIP_CALL(SCIPdropEvent(scip, EventType, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}

static void includeEventHandler(SCIP* scip, const char* name, const char* desc,
    SCIP_DECL_EVENTEXEC((*eventexec)), SCIP_DECL_EVENTINIT((*eventinit)), SCIP_DECL_EVENTEXIT((*eventexit)))
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;
    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(scip, &eventhdlr, name, desc, eventexec, nullptr));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventinit));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventexit));
}

/** Includes an event handler that catches the global event type EventType during the solve. */
template <SCIP_EVENTTYPE EventType>
static void includeEventHandler(SCIP* scip, const char* name, const char* desc, SCIP_DECL_EVENTEXEC((*eventexec)))
{
    includeEventHandler(scip, name, desc, eventexec, eventInitCatch<EventType>, eventExitDrop<EventType>);
}

void includeNodeSolvedEventHandler(SCIP* scip)
{
    includeEventHandler<SCIP_EVENTTYPE_NODESOLVED>(scip,
        "nodesolved_event", "Notify when a node is solved", eventExecNodeSolved);
}

void includeBoundChangedEventHandler(SCIP* scip)
{
    includeEventHandler(scip,
        "boundchanged_event", "Capture bound changes",
        eventExecBoundChanged, eventInitBoundChanged, eventExitBoundChanged);
}

void includeNodeFocusedEventHandler(SCIP* scip)
{
    includeEventHandler<SCIP_EVENTTYPE_NODEFOCUSED>(scip,
        "nodefocused_event", "Capture LB increases on node focus", eventExecNodeFocused);
}

void includeNodeInfeasibleEventHandler(SCIP* scip)
{
    includeEventHandler<SCIP_EVENTTYPE_NODEINFEASIBLE>(scip,
        "nodeinf_event", "detect infeasible nodes", eventExecNodeInfeasible);
}

void includeLPSolveEventHandler(SCIP* scip)
{
    includeEventHandler<SCIP_EVENTTYPE_LPSOLVED>(scip,
        "lpsolve_event", "Print duals after each LP solve", eventExecLPSolved);
}

void includeBestSolFoundEventHandler(SCIP* scip)
{
    includeEventHandler<SCIP_EVENTTYPE_SOLFOUND>(scip,
        "bestsolfound_event", "Notify when a new incumbent integer solution is found", eventExecBestSolFound);
}

void includeNodeDeleteEventHandler(SCIP* scip)
{
    includeEventHandler<SCIP_EVENTTYPE_NODEDELETE>(scip,
        "event_nodedelete", "handles node deletion events for bottom-up aggregation", eventExecNodeDelete);
}

void register_memoscip_branchrule(SCIP* scip)
{
    SCIP_CALL_ABORT(SCIPincludeBranchrule(scip,
        "mybranchrule",
        "My custom branch rule",
        6000000,  // priority
        -1,       // maxdepth
        1.0,      // maxbounddist
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        branchExecLP,
        NULL,
        NULL,
        NULL));
}
