#include "memoscip_events.h"

#include "memoscip_dual.h"
#include "memoscip_option.h"
#include "memoscip_types.h"
#include "sepa_mycuts.h"
#include "scip/pub_var.h"
#include "scip/scip.h"
#include "scip/scip_cons.h"
#include "scip/scip_lp.h"
#include "scip/lp.h"
#include "lpi/lpi.h"
#include "scip/scip_prob.h"
#include "scip/set.h"
#include "scip/struct_lp.h"
#include "scip/struct_scip.h"
#include "scip/var.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cuts.h"
#include "scip/struct_cuts.h"
#include "scip/conflict_general.h"
#include "lpi/lpi.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/scip_lp.h"
#include "scip/scip_prob.h"
#include "scip/var.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
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

// Accumulates elapsed time into acc for the whole lifetime of the enclosing scope, covering every early return.
struct ScopedTimer {
    std::chrono::duration<double>& acc;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    explicit ScopedTimer(std::chrono::duration<double>& acc_) : acc(acc_) {}
    ~ScopedTimer() { acc += std::chrono::steady_clock::now() - start; }
};

std::string trim_copy(const std::string& s)
{
    const size_t first = s.find_first_not_of(" \t\n\r");
    if (first == std::string::npos)
        return "";
    const size_t last = s.find_last_not_of(" \t\n\r");
    return s.substr(first, last - first + 1);
}

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
    SCIP_Real lpobj;
    
    if(g_options.verbosity > 1){
        cout<<"-----NODE SOLVED EVENT -------"<<endl;
        SCIP_STATUS status = SCIPgetStatus(scip);
        cout << "Node " << nodenumber<<" status "<<status<<" "<<g_memoscip.nodeInfoMap[nodenumber].Status<<endl;
    }

    //initializaion of CurrentVarIdx, VarLb, VarUb and AuxVar
    if(g_memoscip.CurrentVarIdx.empty()){
        int nconss = SCIPgetNConss(scip);
        SCIP_CONS** conss = SCIPgetConss(scip);
        g_memoscip.nbctr=nconss;
        int ncols = SCIPgetNLPCols(scip);
        SCIP_COL** cols = SCIPgetLPCols(scip);
        
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
            SCIP_COL* col = cols[i];
            SCIP_VAR* var = SCIPcolGetVar(col);  // ← variable correspondante            
            g_memoscip.AllVars.push_back(var);
            g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-ncols]=i;
            g_memoscip.VarLb[i]=SCIPvarGetLbLocal(var);
            if(SCIPvarGetName(var)[0]!='s' && SCIPvarGetName(var)[2]!='s'){
                if(SCIPvarGetUbLocal(var)==SCIPinfinity(scip))
                    g_memoscip.VarUb[i]=1;
                else
                    g_memoscip.VarUb[i]=SCIPvarGetUbLocal(var);
            }else {
                //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetUbLocal(var)<<" "<<SCIPinfinity(scip)<<endl;
                g_memoscip.VarUb[i]=SCIPvarGetUbLocal(var);
            }
            
            if(g_memoscip.AuxVar[SCIPvarGetIndex(var)-ncols]==-1){
                g_memoscip.AuxVar[SCIPvarGetIndex(var)-ncols]=SCIPvarGetIndex(var)-ncols;
                //assert(SCIPvarGetName(var)[2]=='s');
            }
        }
        bool test=false;
        for(int i = 0; i < nconss; ++i)
        {  
            SCIP_CONS* cons = conss[i];

            // Get variables in this constraint
            SCIP_VAR** vars = SCIPgetVarsLinear(scip,cons);
            int nvars = SCIPgetNVarsLinear(scip,cons);
            if(g_memoscip.FileExt==".opb" && SCIPconsGetName(cons)[0]=='u' and SCIPconsGetName(cons)[1]=='b' ){
                g_memoscip.AuxVar[SCIPvarGetIndex(vars[0])-ncols]=g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[1])-ncols];
                g_memoscip.AuxVar[SCIPvarGetIndex(vars[1])-ncols]=g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[0])-ncols];
                test=true;
            }
            else if((g_memoscip.FileExt==".lp" || g_options.kpcg) && SCIPgetRhsLinear(scip, cons)==1 && nvars==2){
                g_memoscip.AuxVar[SCIPvarGetIndex(vars[0])-ncols]=g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[1])-ncols];
                g_memoscip.AuxVar[SCIPvarGetIndex(vars[1])-ncols]=g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[0])-ncols];
                test=true;
            }
            if(g_options.encoding==1 && SCIPgetRhsLinear(scip, cons)==1 && nvars==4){
                g_memoscip.TupleAuxVar[SCIPvarGetIndex(vars[0])-ncols] = {g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[1])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[2])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[3])-ncols]};
                g_memoscip.TupleAuxVar[SCIPvarGetIndex(vars[1])-ncols] = {g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[0])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[2])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[3])-ncols]};
                g_memoscip.TupleAuxVar[SCIPvarGetIndex(vars[2])-ncols] = {g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[0])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[1])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[3])-ncols]};
                g_memoscip.TupleAuxVar[SCIPvarGetIndex(vars[3])-ncols] = {g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[0])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[1])-ncols], g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[2])-ncols]};
            }
        }
        if(!test)
            cout<<"-----------------No auxiliary variables found !!!!--------------"<<endl;
    }




    SCIP_Bool cutoff, lperror;
    SCIP_LPSOLSTAT stat = SCIPgetLPSolstat(scip);

    if(stat==SCIP_LPSOLSTAT_NOTSOLVED){
       if(g_options.verbosity > 3){
        cout<<"Construct and solve dive LP"<<endl;
       }
        if(!SCIPisLPConstructed(scip)){
            SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
        }
        SCIP_CALL( SCIPstartDive(scip) );
        // Must keep in case problem decided infeasible (by bound propagationf for exampel) before LP is constructed and solved.
        if(!g_memoscip.Unconstructed_cons.empty()){
            for(size_t i =0; i <g_memoscip.Unconstructed_cons.size(); i++){
                SCIP_CONS* cons = g_memoscip.Unconstructed_cons[i];
                SCIP_ROW* row = SCIPgetRowLinear(scip, cons);
                SCIP_CALL( SCIPdelConsLocal(scip, cons) );
                if(row!=nullptr){
                    SCIP_CALL( SCIPchgRowLhsDive(scip, row, -SCIPinfinity(scip)) );
                    SCIP_CALL( SCIPchgRowRhsDive(scip, row, SCIPinfinity(scip)) );
                }
            }
        }
        SCIP_CALL( SCIPflushLP(scip) );
        SCIP_CALL( SCIPsolveDiveLP(scip, -1, &lperror, &cutoff) );
        //SCIP_CALL( SCIPwriteTransProblem(scip, "infeasibleproblem.lp","lp",FALSE) );
        SCIP_LPSOLSTAT stat2 = SCIPgetLPSolstat(scip);
        
        if(stat2==SCIP_LPSOLSTAT_INFEASIBLE){
            if(g_options.verbosity >= 1)
                cout<<"LP INFEASIBLE "<<endl;
            SCIP_CALL(gatherLPInfeasibilityCertificate(scip, nodenumber));
            //SCIP_CALL( SCIPwriteLP(scip, "infeasible.lp") );
            lpobj=100000;
        }else if(stat2==SCIP_LPSOLSTAT_OPTIMAL || stat2==SCIP_LPSOLSTAT_OBJLIMIT){
            lpobj=SCIPgetLPObjval(scip);
            if(g_options.verbosity >= 1)
                cout<<"LP OPTIMAL: "<<lpobj<<endl;        
            //cout<<"lpobj "<<lpobj<<" parent lb "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb<<endl;
            if(!SCIPisDualfeasZero(scip,lpobj) && (fabs(g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb - lpobj)>= g_options.minInc || Myceil(lpobj)>=SCIPgetPrimalbound(scip) || lpobj - Myceil(g_memoscip.nodeInfoMap[parentnumber].lb) > 1e-5))
                AnalyseDualSol(scip,nodenumber );
            else{
                //If lb didn't change then we copy the parent information
                g_memoscip.nodeInfoMap[nodenumber].redcosts = g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].redcosts;
                g_memoscip.nodeInfoMap[nodenumber].lb=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb;
                g_memoscip.nodeInfoMap[nodenumber].y_i=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].y_i;
                g_memoscip.nodeInfoMap[nodenumber].DualSol=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].DualSol;
                g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second=0;
                g_memoscip.nodeInfoMap[nodenumber].Status="FEASIBLE"; 
                g_memoscip.nodeInfoMap[nodenumber].lb_increase=0;
                for (size_t i = 0; i < g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].y_i.size(); i++){
                    if(g_memoscip.VarLb[i]==1){
                        g_memoscip.nodeInfoMap[nodenumber].Var1.push_back(g_memoscip.CurrentVarIdx[SCIPvarGetIndex(g_memoscip.AllVars[i])-g_memoscip.nOriginalVar]);
                    }
                }
            }
        }else{
            cout<<"lperror "<<lperror<<" cutoff "<<cutoff<<" stat2 "<<stat2<<endl;
            cout<<"LP SOLVED EVENT with status "<<stat2<<" node "<<nodenumber<<endl;
            char nameprobsol[50];
            snprintf(nameprobsol, 50, "infeasible%s_%lld.lp", g_options.instancePath.c_str(), nodenumber);
            SCIP_CALL( SCIPwriteTransProblem(scip, nameprobsol,"lp",FALSE) );
            snprintf(nameprobsol, 50, "infeasible%s_%lld.lp", g_options.instancePath.c_str(), nodenumber);
            SCIP_CALL( SCIPwriteLP(scip, nameprobsol) );
             assert(1==0);
        }
           
        SCIP_CALL( SCIPendDive(scip) );
        
        //assert(ceil(lpobj)-SCIPgetPrimalbound(scip)>=-1e-9);
    }else if(stat==SCIP_LPSOLSTAT_INFEASIBLE){
        if(g_options.verbosity >= 1)
            cout<<"LP INFEASIBLE "<<endl;
        SCIP_CALL(gatherLPInfeasibilityCertificate(scip, nodenumber));
        //SCIP_CALL( SCIPwriteLP(scip, "infeasible.lp") );
        lpobj=100000;
    }else if(stat==SCIP_LPSOLSTAT_OPTIMAL || stat==SCIP_LPSOLSTAT_OBJLIMIT  ) {
       
        lpobj=SCIPgetLPObjval(scip);
        if(g_options.verbosity >= 1)
            cout<<"LP OPTIMAL: "<<lpobj<<endl;
        //cout<<"lpobj "<<lpobj<<" parent lb "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb<<endl;
        if(!SCIPisDualfeasZero(scip,lpobj) && (fabs(g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb - lpobj)>= g_options.minInc || Myceil(lpobj)>=SCIPgetPrimalbound(scip) || lpobj - Myceil(g_memoscip.nodeInfoMap[parentnumber].lb) > 1e-5))
            AnalyseDualSol(scip,nodenumber );
        else{
                //If lb didn't change then we copy the parent information
                g_memoscip.nodeInfoMap[nodenumber].redcosts = g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].redcosts;
                g_memoscip.nodeInfoMap[nodenumber].lb=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb;
                g_memoscip.nodeInfoMap[nodenumber].y_i=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].y_i;
                g_memoscip.nodeInfoMap[nodenumber].DualSol=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].DualSol;
                g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second=0;
                g_memoscip.nodeInfoMap[nodenumber].Status="FEASIBLE"; 
                g_memoscip.nodeInfoMap[nodenumber].lb_increase=0;
                for (size_t i = 0; i < g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].y_i.size(); i++){
                    if(g_memoscip.VarLb[i]==1){
                        g_memoscip.nodeInfoMap[nodenumber].Var1.push_back(g_memoscip.CurrentVarIdx[SCIPvarGetIndex(g_memoscip.AllVars[i])-g_memoscip.nOriginalVar]);
                    }
                }
            }
    }else{
        cout<<"LP SOLVED EVENT with status "<<stat<<" node "<<nodenumber<<endl;
        assert(1==0);
    }
    
    //cout<<g_memoscip.nodeInfoMap[SCIPnodeGetNumber(node)].Status<<endl;
    if(stat==SCIP_LPSOLSTAT_OPTIMAL){
        assert(g_memoscip.nodeInfoMap[nodenumber].Status=="FEASIBLE");
    }
   
    //cout<<SCIPgetNodeLowerbound(scip, node)<<" - "<<g_memoscip.nodeInfoMap[SCIPnodeGetNumber(node)].lb<<" = "<<fabs(SCIPgetNodeLowerbound(scip, node)-g_memoscip.nodeInfoMap[SCIPnodeGetNumber(node)].lb)<<endl;
    if(g_memoscip.nodeInfoMap[nodenumber].lb>=SCIPgetPrimalbound(scip))
        g_memoscip.nodeInfoMap[nodenumber].LbgUb=true;
    
    if(g_options.verbosity > 2)
        cout<<"END OF NODE SOLVED EVENT  node "<<nodenumber<<" my lb "<<g_memoscip.nodeInfoMap[nodenumber].lb<<" SCIP node lb "<<SCIPgetNodeLowerbound(scip, node)<<" ub "<<SCIPgetPrimalbound(scip)<< " diff between lb and scip node lb "<<fabs(SCIPgetNodeLowerbound(scip, node)-g_memoscip.nodeInfoMap[nodenumber].lb)<<endl;
    if((g_options.output>1 && nodenumber >= g_options.output)){
        char nameprobsol[50];
        snprintf(nameprobsol, 50, "solved_%lld.lp", nodenumber);
        SCIP_CALL( SCIPwriteLP(scip, nameprobsol) ); 
        snprintf(nameprobsol, 50, "problem_%lld.lp", nodenumber);
        SCIP_CALL( SCIPwriteTransProblem(scip, nameprobsol,"lp",FALSE) );
    }
    
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTFREE(eventFreeNodeSolved)
{
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTINIT(eventInitNodeSolved)
{
    /* 'eventhdlr' parameter is provided by SCIP (do not redeclare) */
    SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, nullptr, nullptr));
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXIT(eventExitNodeSolved)
{
    SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}

void includeNodeSolvedEventHandler(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;
    //std::cout << "Include NODESOLVED Event Handler" << std::endl;

    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(scip, &eventhdlr,
        "nodesolved_event", "Notify when a node is solved",
        eventExecNodeSolved, nullptr));

    SCIP_CALL_ABORT(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeNodeSolved));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitNodeSolved));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitNodeSolved));
}
// Event handler execution callback for LP solved (existing)
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
    //SCIP_Longint parentNumber = g_memoscip.nodeInfoMap[nodenumber].parent;
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
    //cout<<"END LP SOLVED EVENT "<<nodenumber<<" "<<g_memoscip.nodeInfoMap[nodenumber].Status<<" "<<g_memoscip.nodeInfoMap[nodenumber].lb<<" "<<SCIPgetPrimalbound(scip)<<endl;
    return SCIP_OKAY;
}


// ---- New: NODEFOCUSED event handler ----
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

    if(!g_memoscip.SlackNonNegConss.empty()){
        for(SCIP_CONS* cons : g_memoscip.SlackNonNegConss){
            cout<<"Deleting slack nonnegativity constraint "<<SCIPconsGetName(cons)<<endl;
            //SCIP_CALL(SCIPdelCons(scip, cons));
            //SCIP_CALL(SCIPreleaseCons(scip, &cons));
        }
        g_memoscip.SlackNonNegConss.clear();
    }

    SCIP_Longint nodenumber=SCIPnodeGetNumber(node);

    /*if(nodenumber==87){
        SCIP_VAR* slack = nullptr;
            char slackname[50];
            snprintf(slackname, 50, "sc_%d", g_memoscip.nbAddedctr);
            SCIP_CALL( SCIPcreateVarBasic(scip, &slack, slackname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL(SCIPvarSetRemovable(slack,FALSE));
            SCIP_CALL(SCIPvarSetInitial(slack,TRUE));
            SCIP_CALL( SCIPaddVar(scip, slack) );

            SCIP_CONS* slackcons = nullptr;
            SCIP_Real slackcoef = 1.0;
            char slackconsname[64];
            snprintf(slackconsname, 64, "slack_nonneg_%d", g_memoscip.nbAddedctr);
            SCIP_CALL( SCIPcreateConsLinear(scip, &slackcons, slackconsname,
                                    1,
                                    &slack,
                                    &slackcoef,
                                    5,
                                    SCIPinfinity(scip),
                                    TRUE,    // initial
                                    FALSE,   // separate
                                    TRUE,    // enforce
                                    TRUE,    // check
                                    TRUE,    // propagate
                                    FALSE,   // local
                                    FALSE,   // modifiable
                                    TRUE,   // dynamic
                                    FALSE,   // removable
                                    FALSE) );// stickingatnode
            SCIP_CALL(SCIPaddCons(scip, slackcons));
            g_memoscip.SlackNonNegConss.push_back(slackcons);
            SCIP_CALL(SCIPreleaseCons(scip, &slackcons));
            g_memoscip.nbAddedctr+=1;
    }*/
    SCIP_NODE* parent = nullptr;
    SCIP_Longint parentNumber=0;
    NodeInfo &ninfo = g_memoscip.nodeInfoMap[nodenumber]; // creates if not existing
    ninfo.visited = true;
    //ninfo.lb = node_lb;
    ninfo.number = nodenumber;
    ninfo.ub = SCIPgetPrimalbound(scip);
    if(nodenumber!=1)  // skip root
    {
        parent = SCIPnodeGetParent(node);
        parentNumber = SCIPnodeGetNumber(parent);
        if(g_options.verbosity > 1)
            cout<<"-----------FOCUSING NODE "<<nodenumber<<" depth "<<SCIPgetDepth(scip)<<" with Parent "<<parentNumber<<"---------------------and "<<g_memoscip.nbAddedctr<<" CONSTRAINT added "<<endl<<"UB :"<< SCIPgetPrimalbound(scip)<<endl<<" LB : "<<SCIPnodeGetLowerbound(node)<<endl;
        g_memoscip.nodeInfoMap[nodenumber].parent=parentNumber;

    }else{
        //CHECK ROOT NODE ASSIGNMENT 
    }
   
    SCIP_DOMCHG* domchg;
    domchg = SCIPnodeGetDomchg(node);
    int nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
    for(int i = 0; i < nboundchgs; ++i )
    {
        SCIP_BOUNDCHG* boundchg;
        SCIP_VAR* var;
        SCIP_Real newbound;
        
        boundchg = SCIPdomchgGetBoundchg(domchg, i);
        var = SCIPboundchgGetVar(boundchg);
        newbound = SCIPboundchgGetNewbound(boundchg);
        //cout<<"At node "<<nodenumber<<" Bound change for var "<<SCIPvarGetName(var)<<" to "<<newbound<<" with obj coef "<<SCIPvarGetObj(var)<<endl;
        //cout<<"Global lb of var "<<SCIPvarGetName(var)<<" is "<<SCIPvarGetLbGlobal(var)<< " "<< SCIPvarGetUbGlobal(var)<<endl;
        if(newbound==1.0){
            ninfo.BranchVar = SCIPvarGetIndex(var)-g_memoscip.nOriginalVar;
            
        }else{
            ninfo.BranchVar = g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar];
        }        
    }
    assert(g_memoscip.VarLb.size()==g_memoscip.VarUb.size());

    if(g_options.debug){
        for (size_t i = 0; i < g_memoscip.VarLb.size(); i++)
        {   
            SCIP_VAR* var= g_memoscip.AllVars[i];
            //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetLbLocal(var)<<" "<<g_memoscip.VarLb[i]<<endl;
            //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetUbLocal(var)<<" "<<g_memoscip.VarUb[i]<<endl;
            assert(SCIPvarGetLbLocal(var)==g_memoscip.VarLb[i] || SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]])==0);
            assert(SCIPvarGetUbLocal(var)==g_memoscip.VarUb[i] ||  SCIPvarGetUbLocal(var)==SCIPinfinity(scip) || SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]])==1);
        }
    }
    int qsize=g_memoscip.CutQueue.queue.size();
    if(qsize>g_options.maxCutQueueSize && !g_options.sepa){
        //cout<<"add "<<qsize<<" cuts to problem at node "<<nodenumber<<" with parent "<<parentNumber<<" and depth "<<SCIPgetDepth(scip)<<" and "<<g_memoscip.nbAddedctr<<" constraints added"<<endl;
        SCIP_CALL(AddCutstoProb(scip,&g_memoscip.CutQueue,&g_memoscip,double(g_memoscip.TOP),g_options.sepa,&g_memoscip.nbctr,nodenumber,nullptr)); 
        //SCIP_CALL( SCIPwriteLP(scip, "queue.lp") );    
        //SCIP_CALL( SCIPflushLP(scip) );
    }
    // link this node as child of its parent (if parent exists)
    if (parent != nullptr)
    {
        // add child pointer to parent's children list if not already present
        g_memoscip.nodeInfoMap[parentNumber].children.push_back(SCIPnodeGetNumber(node));
        //cout<<SCIPnodeGetLowerbound(node)<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb -SCIPnodeGetLowerbound(node)<<" "<<endl;
        //assert(fabs(g_memoscip.nodeInfoMap[parentNumber].lb -SCIPnodeGetLowerbound(node)) <= 1 || Myceil(g_memoscip.nodeInfoMap[parentNumber].lb) - SCIPgetPrimalbound(scip)>= -1e-9 || g_memoscip.nodeInfoMap[parentNumber].Status=="INFEASIBLE");
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
    if( node == nullptr )
        return SCIP_OKAY;
    if(nodenumber==1) return SCIP_OKAY; // skip root
    SCIP_NODE* parent = SCIPnodeGetParent(node);
    SCIP_Longint parentNumber = g_memoscip.nodeInfoMap[nodenumber].parent;
    if(g_options.verbosity > 1){
        cout<<"Deleting node "<<nodenumber<<" from nodeInfoMap, scip Node Lb "<< SCIPgetNodeLowerbound(scip,node)  <<" my node lb : "<<g_memoscip.nodeInfoMap[nodenumber].lb<< " "<<" depth : "<<SCIPgetDepth(scip)<< " parent: "<<parentNumber<<endl;
        cout<<"UB "<<SCIPgetPrimalbound(scip)<<endl;
        cout<<"number of children "<<g_memoscip.nodeInfoMap[nodenumber].children.size()<<endl;
        cout<<"Bracnh var "<<g_memoscip.nodeInfoMap[nodenumber].BranchVar<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[nodenumber].BranchVar]])<<endl;
        cout<<"Curnnet node lb increase "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" versus children "<<endl;
    }

    //Case node was not processed. Mostly because it was pruned immediately after creation, since a brother node found a better solution than the pseudo obj of this node. 
    if(parentNumber==-1){
        SCIP_Longint SCIP_parentNumber = SCIPnodeGetNumber(parent);
        if(SCIPisLT(scip, ceil(g_memoscip.nodeInfoMap[SCIP_parentNumber].lb), SCIPgetPrimalbound(scip))){
            int nvars;
            // D'abord interroger la taille nécessaire avec branchvarssize = 0
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
                    decvar =g_memoscip.CurrentVarIdx[SCIPvarGetIndex(vars[i])-g_memoscip.nOriginalVar];
                    assert(SCIPvarGetName(vars[i]) == SCIPvarGetName(g_memoscip.AllVars[decvar]));
                    assert(bounds[i]==1);
                    g_memoscip.nodeInfoMap[nodenumber].BranchVar= SCIPvarGetIndex(vars[i])-g_memoscip.nOriginalVar;
                }
            }
            
            assert(g_memoscip.nodeInfoMap[nodenumber].lb==0);
            assert(SCIP_parentNumber>0);
            SCIP_CALL(getDualProofFromPseudoObjCutoff(scip,nodenumber, SCIP_parentNumber, decvar));
            g_memoscip.nodeInfoMap[SCIP_parentNumber].children.push_back(nodenumber);
            return SCIP_OKAY;
        }
    }

    vector<long> children = g_memoscip.nodeInfoMap[nodenumber].children;
    //cout<< "NodeInfoMap size before deletion: " << g_memoscip.nodeInfoMap.size() << endl;
    //SCIPinfoMessage(scip, NULL, "[DEBUG] Total LP rows now: %d\n", SCIPgetNLPRows(scip));
    //Print children of the node being deleted
    if (g_memoscip.nodeInfoMap.contains(nodenumber) && g_memoscip.nbAddedctr<(int)g_options.nbslack && (g_memoscip.nodeInfoMap[nodenumber].lb_increase || children.size()==2) ) {
        
        vector<SCIP_Real> ResFus(g_memoscip.AllVars.size(),0.0);
        SCIP_Real rhsFus=0;
        if(children.size()==1){
            assert(g_memoscip.nodeInfoMap.contains(children[0]));
            if(g_options.verbosity > 2)
                cout<<"Single child "<<children[0]<<" lb "<<g_memoscip.nodeInfoMap[children[0]].lb<<" versus parent lb "<<g_memoscip.nodeInfoMap[nodenumber].lb<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.first.size()<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second<<" incr lb "<<g_memoscip.nodeInfoMap[children[0]].lb_increase<<endl;
            assert(g_memoscip.nodeInfoMap[children[0]].lb<g_memoscip.TOP);
            assert(fabs(g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second-g_memoscip.nodeInfoMap[children[0]].lb_increase)<1e-5);
        }
        else if(children.size()==2){
            if(g_options.verbosity > 2){
                cout<<"Two children "<<children[0]<<" lb "<<g_memoscip.nodeInfoMap[children[0]].lb<<" versus parent lb "<<g_memoscip.nodeInfoMap[nodenumber].lb<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.first.size()<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second<<" incr lb "<<g_memoscip.nodeInfoMap[children[0]].lb_increase<<endl;
                cout<<"Two children "<<children[1]<<" lb "<<g_memoscip.nodeInfoMap[children[1]].lb<<" versus parent lb "<<g_memoscip.nodeInfoMap[nodenumber].lb<<" "<<g_memoscip.nodeInfoMap[children[1]].DualProofCoef.first.size()<<" "<<g_memoscip.nodeInfoMap[children[1]].DualProofCoef.second<<" incr lb "<<g_memoscip.nodeInfoMap[children[1]].lb_increase<<endl;
            }
            if(g_memoscip.nodeInfoMap[children[0]].lb_increase>0 && g_memoscip.nodeInfoMap[children[1]].lb_increase>0){
                assert(g_memoscip.nodeInfoMap.contains(children[0]) && g_memoscip.nodeInfoMap.contains(children[1]));
                
                //cout<<children[0]<<" lb "<<g_memoscip.nodeInfoMap[children[0]].lb<<"  "<<children[1]<<" lb "<<g_memoscip.nodeInfoMap[children[1]].lb<<" versus parent lb "<<g_memoscip.nodeInfoMap[nodenumber].lb<<endl;
                if(g_options.verbosity > 1){
                    cout<<"Fusing proofs of children "<<children[0]<<" and "<<children[1]<<" to recover parent node "<<nodenumber<<" proof "<<endl;
                    cout<<"proofs of children "<<children[0]<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second<<" lb increase: "<<g_memoscip.nodeInfoMap[children[0]].lb_increase<<endl;
                    cout<<"proofs of children "<<children[1]<<" "<<g_memoscip.nodeInfoMap[children[1]].DualProofCoef.second<<" lb increase: "<<g_memoscip.nodeInfoMap[children[1]].lb_increase<<endl;
                }
                
                assert(fabs(g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second-g_memoscip.nodeInfoMap[children[0]].lb_increase)<1e-5 || g_memoscip.nodeInfoMap[children[0]].LbgUb);
                assert(fabs(g_memoscip.nodeInfoMap[children[1]].DualProofCoef.second-g_memoscip.nodeInfoMap[children[1]].lb_increase)<1e-5 || g_memoscip.nodeInfoMap[children[1]].LbgUb);

                assert(g_memoscip.nodeInfoMap[children[1]].BranchVar==g_memoscip.AuxVar[SCIPvarGetIndex(g_memoscip.AllVars[g_memoscip.nodeInfoMap[children[0]].BranchVar])-g_memoscip.nOriginalVar]);
                ResFus=FusionResolution(scip,children[0],children[1]);
                
                if(g_memoscip.nodeInfoMap[children[0]].Status!="INFEASIBLE" && g_memoscip.nodeInfoMap[children[1]].Status!="INFEASIBLE"){
                    rhsFus=min(g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second,g_memoscip.nodeInfoMap[children[1]].DualProofCoef.second);
                }else if(g_memoscip.nodeInfoMap[children[0]].Status!="INFEASIBLE"){
                    rhsFus=g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second;
                }else if(g_memoscip.nodeInfoMap[children[1]].Status!="INFEASIBLE"){
                    rhsFus=g_memoscip.nodeInfoMap[children[1]].DualProofCoef.second;
                }else{
                    rhsFus=g_memoscip.TOP;
                    cout<<"This should not happen "<<endl;
                    SCIP_CALL( SCIPwriteTransProblem(scip, "model.lp","lp",FALSE) );
                    assert(1==0);
                }
            }
        }else{
            if(g_options.verbosity > 2)
                cout<<"No child "<<endl;
        }
        assert(ResFus.size()==g_memoscip.AllVars.size());
        //cout<<"Result of Fusion Resolution: rhs "<< rhsFus<<endl;
        //cout<<g_memoscip.nodeInfoMap[nodenumber].redcosts.size()<<"' "<<ResFus.size()<<endl;

        //PRINT RESFUS
        if(g_options.verbosity > 2){
            cout<<"ResFus "<<endl;
            for (size_t i = 0; i < ResFus.size(); i++)
            {
                if(!SCIPisDualfeasZero(scip,ResFus[i])){
                    if(ResFus[i]>0)
                        cout<<"+"<<ResFus[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[i])<<" ";
                    else
                        cout<<ResFus[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[i])<<" ";
                }

            }
            cout<<" >= "<<rhsFus<<endl;
        }
        SCIP_Real rhs_summed=rhsFus;
        vector<int> NewInds;
        vector<SCIP_Real> NewCoefs;
        vector<SCIP_VAR*> NewVar;

        if(rhs_summed>0 && g_options.learning==2 && rhs_summed < g_memoscip.TOP){
                char name[100];

            snprintf(name, 100, "l2cut_%lld_%d", nodenumber,SCIPgetDepth(scip));
            SCIP_Real rhs_ctr=rhs_summed;
            vector<SCIP_Real> ResFus2=ResFus;

            for (size_t i = 0; i < ResFus2.size(); i++)
            {
                if(!SCIPisDualfeasZero(scip,ResFus2[i])){
                    if((int)i!= g_memoscip.AuxVar[i] && !SCIPisDualfeasZero(scip,ResFus2[g_memoscip.AuxVar[i]])){
                        if(ResFus2[i]!= 0 && ResFus2[g_memoscip.AuxVar[i]]!=0){
                            SCIP_Real minval = min(ResFus2[i],ResFus2[g_memoscip.AuxVar[i]]);
                            rhs_ctr -= minval;
                            ResFus2[i]-=minval;
                            ResFus2[g_memoscip.AuxVar[i]]-=minval;
                            assert(ResFus2[i]==0 || ResFus2[g_memoscip.AuxVar[i]]==0);
                        }
                    }
                    
                    if(g_options.encoding==1 && !g_memoscip.TupleAuxVar.empty()){
                        SCIP_Real minval = ResFus2[i];
                        for (int j = 0; j < (int)g_memoscip.TupleAuxVar[i].size(); j++)
                            minval = min(minval, ResFus2[g_memoscip.TupleAuxVar[i][j]]);
                        rhs_ctr -= minval;
                        ResFus2[i]-=minval;
                        ResFus2[g_memoscip.TupleAuxVar[i][0]]-=minval;
                        ResFus2[g_memoscip.TupleAuxVar[i][1]]-=minval;
                        ResFus2[g_memoscip.TupleAuxVar[i][2]]-=minval;
                    }
                    NewInds.push_back(i);
                    NewCoefs.push_back(ResFus2[i]);
                    NewVar.push_back(g_memoscip.AllVars[i]);
                    //assert(ResFus[i]<g_memoscip.TOP);
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
                                rhs_ctr,  // lhs (pour inégalité)
                                rhs_ctr,&g_memoscip.CutQueue) ); 
                            
            //g_memoscip.nbAddedctr++;
            NewInds.clear();
            NewCoefs.clear();
            NewVar.clear();
        }

        assert(fabs(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second-g_memoscip.nodeInfoMap[nodenumber].lb_increase)<1e-5);
        
        if(g_memoscip.nodeInfoMap[nodenumber].lb_increase>0){
            //cout<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.size()<<endl;
            for (size_t i = 0; i < g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.size(); i++)
            {
                int idx =g_memoscip.nodeInfoMap[nodenumber].DualProofInds[i];
                ResFus[idx]+=g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i];
                //cout<<idx<<" "<<SCIPvarGetName(g_memoscip.AllVars[idx])<<" "<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i]<<" "<<ResFus[idx]<<endl;
                //cout<<"red cost "<<g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<" FIncal ResFus "<<ResFus[idx]<<endl;
                //assert(ResFus[idx]<g_memoscip.TOP);
                
            }
            rhs_summed+=g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second;
        }
        
        if(g_options.verbosity > 4)
            cout<<"BOund propagation "<<g_memoscip.nodeInfoMap[nodenumber].BPRemoval.size()<<endl;

        for(int i = (int)g_memoscip.nodeInfoMap[nodenumber].BPRemoval.size()-1; i >= 0; i--) {
            if(g_options.verbosity > 4)
                cout<<"removed "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].first])<<" "<<g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second.size()<<endl;
            int idxi = g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].first;
            if(ResFus[idxi]>1e-5){
                assert(ResFus[idxi]>1e-5);
                for(size_t j = 0; j < g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second.size(); j++) {
                    if(g_options.verbosity > 4)
                        cout<<"explanation "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second[j]])<<endl;
                    int idxj =g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second[j];
                    //Temporary fix because negative coef are hard to handle.
                    ResFus[idxj]=rhs_summed*100;

                }
            }
            ResFus[idxi]=0;
        }
        if(g_options.verbosity > 4)
            cout<<"Node consistency "<<g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.size()<<endl;
        for (size_t i = 0; i < g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.size(); i++)
        {
            int idx =g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].first;
            ResFus[idx]=g_memoscip.nodeInfoMap[parentNumber].redcosts[idx];
            rhs_summed=min(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second,rhs_summed);
            if(g_options.verbosity > 4)
                cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[idx]])<< " "<< idx<<" "<<g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second<<" "<<Myceil(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second+g_memoscip.nodeInfoMap[nodenumber].lb)<<" "<<SCIPgetPrimalbound(scip)<<endl;
            //assert(Myceil(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second+g_memoscip.nodeInfoMap[nodenumber].lb)>=SCIPgetPrimalbound(scip));
        }
        
        
        if(g_options.verbosity > 2){
            cout<<"final res "<<endl;
            for (size_t i = 0; i < ResFus.size(); i++)
            {
                if(!SCIPisDualfeasZero(scip,ResFus[i])){
                    if(ResFus[i]>0)
                        cout<<"+"<<ResFus[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
                    else
                        cout<<ResFus[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
                }

            }
            cout<<" >= "<<rhs_summed<<endl;
        }

        if(g_options.debug){
            if(children.size()==2){
                //cout<<"Adding cut with rhs "<<rhs_summed<<" vs "<<(min(g_memoscip.nodeInfoMap[children[1]].lb_increase,g_memoscip.nodeInfoMap[children[0]].lb_increase) + g_memoscip.nodeInfoMap[nodenumber].lb_increase )<<endl;
                //cout<<"Min children "<<min(g_memoscip.nodeInfoMap[children[1]].lb_increase,g_memoscip.nodeInfoMap[children[0]].lb_increase) <<endl;
                assert(abs(rhs_summed - (min(g_memoscip.nodeInfoMap[children[1]].lb_increase,g_memoscip.nodeInfoMap[children[0]].lb_increase) + g_memoscip.nodeInfoMap[nodenumber].lb_increase )<1e-9)|| !g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.empty()); 
            }if(children.size()==1){
                //cout<<"Adding cut with rhs "<<rhs_summed<<" vs "<<g_memoscip.nodeInfoMap[children[0]].lb_increase + g_memoscip.nodeInfoMap[nodenumber].lb_increase<<endl;
                //cout<<"Min children "<<g_memoscip.nodeInfoMap[children[0]].lb_increase <<endl;
                assert(abs(rhs_summed - (g_memoscip.nodeInfoMap[children[0]].lb_increase + g_memoscip.nodeInfoMap[nodenumber].lb_increase )<1e-9) || !g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.empty()); 
            }
        }
    
        

        //cout<<"rhs summed "<<rhs_summed<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb<<endl;
        //cout<<Myceil(rhs_summed+g_memoscip.nodeInfoMap[parentNumber].lb)<<" "<<SCIPgetPrimalbound(scip)<<endl;
        //assert(SCIPgetPrimalbound(scip) == SCIPinfinity(scip)  || Myceil(rhs_summed+g_memoscip.nodeInfoMap[parentNumber].lb) - SCIPgetPrimalbound(scip) > -1e-9);

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
            for (size_t i = 0; i < ResFus.size(); i++)
            {

                if(!SCIPisDualfeasZero(scip,ResFus[i])){
                    assert(g_options.encoding==0 || g_memoscip.TupleAuxVar[i].empty() || g_memoscip.TupleAuxVar[i].size() == 3);
                    if((int)i!= g_memoscip.AuxVar[i] && !SCIPisDualfeasZero(scip,ResFus[g_memoscip.AuxVar[i]])){
                        if(ResFus[i]!= 0 && ResFus[g_memoscip.AuxVar[i]]!=0){
                            SCIP_Real minval = min(ResFus[i],ResFus[g_memoscip.AuxVar[i]]);
                            rhs_ctr -= minval;
                            ResFus[i]-=minval;
                            ResFus[g_memoscip.AuxVar[i]]-=minval;
                            assert(ResFus[i]==0 || ResFus[g_memoscip.AuxVar[i]]==0);
                        }
                    }

                    if(g_options.encoding==1 && !g_memoscip.TupleAuxVar[i].empty()){
                        SCIP_Real minval = ResFus[i];
                        for (int j = 0; j < (int)g_memoscip.TupleAuxVar[i].size(); j++)
                            minval = min(minval, ResFus[g_memoscip.TupleAuxVar[i][j]]);
                        rhs_ctr -= minval;
                        ResFus[i]-=minval;
                        ResFus[g_memoscip.TupleAuxVar[i][0]]-=minval;
                        ResFus[g_memoscip.TupleAuxVar[i][1]]-=minval;
                        ResFus[g_memoscip.TupleAuxVar[i][2]]-=minval;
                    }

                    if(i<g_memoscip.nodeInfoMap[parentNumber].redcosts.size() && g_options.verbosity > 4){
                      cout<<SCIPvarGetName(g_memoscip.AllVars[i])<< " i "<< i<<" resfus "<<ResFus[i]<<" parentrc "<<g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<" diff "<<ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<" Bounds "<<g_memoscip.VarUb[i]<<" "<< g_memoscip.VarLb[i]<<endl;
                    }
                    NewInds.push_back(i);
                    NewCoefs.push_back(ResFus[i]);
                    NewVar.push_back(g_memoscip.AllVars[i]);
                    
                    //assert(ResFus[i]<g_memoscip.TOP);
                }
            }
            
            if(g_options.verbosity > 2)
                cout<<"NODE DELETE Adding cut "<<name<<" with rhs "<<rhs_ctr<<" and "<<NewVar.size()<< " variables vs Original vars "<< g_memoscip.nOriginalVar<< endl;

            assert(NewVar.size()==NewCoefs.size());
            assert(NewVar.size()!=0);
            if(g_options.learning==4 || (children.size()==2 && (g_options.learning==3 || g_memoscip.nodeInfoMap[nodenumber].lb_increase>0 || (g_memoscip.nodeInfoMap[parentNumber].lb_increase==0 && (g_memoscip.nodeInfoMap[children[0]].lb_increase>0 || g_memoscip.nodeInfoMap[children[1]].lb_increase>0 ))) ))
                SCIP_CALL( SCIPqueueCut(scip, name,
                                (int)NewVar.size(),
                                NewVar.data(),
                                NewCoefs.data(),
                                rhs_ctr,  // lhs (pour inégalité)
                                rhs_ctr,&g_memoscip.CutQueue) ); 
            if(g_options.verbosity > 2)
                SCIPprintQueuedCuts(scip, &g_memoscip.CutQueue);
            //g_memoscip.nbAddedctr++;
        }
        
        //cout<<"Dual proof coef second "<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second<<" lb increase "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" rhs_summed "<<rhs_summed<<endl;
        assert(abs(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second-g_memoscip.nodeInfoMap[nodenumber].lb_increase) < 1e-9);
        assert(NewCoefs.size()<=g_memoscip.CurrentVarIdx.size());
        g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second=rhs_summed;
        g_memoscip.nodeInfoMap[nodenumber].DualProofInds=Inds;
        g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first=Coefs;
        assert(g_memoscip.nodeInfoMap[nodenumber].lb_increase<=rhs_summed || !g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.empty());
        g_memoscip.nodeInfoMap[nodenumber].lb_increase=rhs_summed;
        //cout<<" plb "<<g_memoscip.nodeInfoMap[parentNumber].lb<<" lbinc "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" ub "<<SCIPgetPrimalbound(scip)<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb + g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" "<<Myceil(g_memoscip.nodeInfoMap[parentNumber].lb + g_memoscip.nodeInfoMap[nodenumber].lb_increase)<<endl;
        assert(Myceil(g_memoscip.nodeInfoMap[parentNumber].lb + g_memoscip.nodeInfoMap[nodenumber].lb_increase) >= SCIPgetPrimalbound(scip)-1 ||SCIPgetPrimalbound(scip)==SCIPinfinity(scip) );
        if(!children.empty()){
            g_memoscip.nodeInfoMap.erase(children[0]);
            if(children.size()==2)  
                g_memoscip.nodeInfoMap.erase(children[1]);
            g_memoscip.nodeInfoMap[nodenumber].children.clear();
        }  
        
    }
    return SCIP_OKAY;
}


/** Parcourt la pile lb/ub d'une variable active (une seule représentation SCIP_VAR). */
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

    /* 1) Même (ancienne, nouvelle) borne — aligné sur l'événement, non redondant. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (!pairMatches(info))
            continue;
        if (!SCIPbdchginfoIsRedundant(info))
            return info;
    }
    /* 2) Même paire même si redondant. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (pairMatches(info))
            return info;
    }
    /* 3) Même nouvelle borne seule (comportement précédent), non redondant. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (!newMatches(info))
            continue;
        if (!SCIPbdchginfoIsRedundant(info))
            return info;
    }
    /* 4) Même nouvelle borne seule, redondant autorisé. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (newMatches(info))
            return info;
    }
    /* 5) Dernière entrée non redondante. */
    for (int pos = n - 1; pos >= 0; --pos) {
        SCIP_BDCHGINFO* info = getAtPos(pos);
        if (!SCIPbdchginfoIsRedundant(info))
            return info;
    }
    return getAtPos(n - 1);
}

/** Résout le SCIP_BDCHGINFO aligné sur l'événement ; essaie la var de l'événement puis SCIPvarGetProbvar (chaîne neg/agg). */
static SCIP_BDCHGINFO* getBdchgInfoForEventBound(
    SCIP* scip,
    SCIP_VAR* var,
    SCIP_BOUNDTYPE boundtype,
    SCIP_Real eventOldbound,
    SCIP_Real eventNewbound)
{
    SCIP_VAR* prob = SCIPvarGetProbvar(var);

    SCIP_BDCHGINFO* info = searchBdchgStackForEvent(scip, var, boundtype, eventOldbound, eventNewbound);
    if (info != NULL)
        return info;
    if (prob != NULL && prob != var)
        info = searchBdchgStackForEvent(scip, prob, boundtype, eventOldbound, eventNewbound);
    return info;
}

SCIP_BDCHGINFO* getBoundChangeReasonDetailed(SCIP* scip, SCIP_VAR* var,
    SCIP_BOUNDTYPE boundtype,
    SCIP_Real oldbound,
    SCIP_Real newbound)
{
    SCIP_BDCHGINFO* bdchginfo = getBdchgInfoForEventBound(scip, var, boundtype, oldbound, newbound);
    return bdchginfo;
    if (bdchginfo == NULL) {
        return NULL;
    }
}

static SCIP_DECL_EVENTEXEC(eventExecBoundChanged)
{
    ScopedTimer BoundChangedTimer(g_memoscip.BoundChangedTime);
    SCIP_VAR* var = SCIPeventGetVar(event);
    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    SCIP_Real newbound;
    SCIP_BDCHGINFO* reason;
    SCIP_NODE* node = SCIPgetCurrentNode(scip);
    if(node == NULL || g_memoscip.CurrentVarIdx.empty()){
        return SCIP_OKAY;
    }
    SCIP_Longint nodenumber=SCIPnodeGetNumber(node);
        
    if(g_options.verbosity > 1)
        cout<<"-----BOUND CHANGED EVENT ------- "<< nodenumber<<" "<<" "<<SCIPvarGetName(var)<<" "<<SCIPgetDepth(scip)<<endl;

    SCIP_Real oldbound = SCIPeventGetOldbound(event);
    newbound = SCIPeventGetNewbound(event);
   
    if (eventtype & SCIP_EVENTTYPE_LBCHANGED) {
        if(g_options.verbosity > 1){
            cout<<"LB CHANGED "<<SCIPvarGetName(var)<<" "<<oldbound<<" ->"<<newbound<<endl;
        }
        assert(SCIPvarGetName(var)[0]!='s' && SCIPvarGetName(var)[2]!='s');
        g_memoscip.VarLb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=newbound;
        bool auxvarassigned=false;
        if(g_options.encoding==0 || g_memoscip.TupleAuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar].empty()){
            auxvarassigned =  g_memoscip.VarUb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]] == g_memoscip.VarLb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]];
            if(newbound==1)
                g_memoscip.VarUb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=0;
            else if(newbound==0)
                g_memoscip.VarUb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=1;
            else
                g_memoscip.VarUb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=newbound ;
        }
        
        reason = getBoundChangeReasonDetailed(scip, var, SCIP_BOUNDTYPE_LOWER, oldbound, newbound);

         
        
        assert(!SCIPinProbing(scip));
        //if reason is COnstraint :
        if(reason!=NULL && !auxvarassigned){
            SCIP_BOUNDCHGTYPE chgtype = SCIPbdchginfoGetChgtype(reason);
            if(g_options.verbosity > 2)
                cout<<"reason found LB Change "<<SCIPvarGetName(var)<<" "<<oldbound<<" "<<newbound<<endl;
            switch (chgtype) {
                case SCIP_BOUNDCHGTYPE_CONSINFER: {
                    SCIP_CONS* cons = SCIPbdchginfoGetInferCons(reason);
                    if (cons != NULL) {
                        std::string constype = SCIPconshdlrGetName(SCIPconsGetHdlr(cons));
                        std::string consname = SCIPconsGetName(cons);
                        if(g_options.verbosity > 2)
                            cout<< "CONSTRAINT (" + constype + ": " + consname + ")"<<endl;
                        //Loop over variables in constraint and coef
                        SCIP_VAR** vars = SCIPgetVarsLinear(scip,cons);
                        int nvars = SCIPgetNVarsLinear(scip,cons);
                        SCIP_Real* coefs = SCIPgetValsLinear(scip,cons);                        
                        if(nvars>2){
                            vector<int> explanation;
                            bool varsign;
                            for(int i=0; i < nvars; i++) {
                                SCIP_VAR* v = vars[i];
                                if(v == var){
                                    varsign = coefs[i] > 0;
                                    break;
                                }
                            }
                            for(int i = 0; i < nvars; i++) {
                                SCIP_VAR* v = vars[i];
                                if(v != var){
                                    if(g_options.verbosity > 3 )
                                        cout<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;                                    
                                    assert(g_memoscip.VarLb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==SCIPvarGetLbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]])));
                                    assert(g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==SCIPvarGetUbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)- g_memoscip.nOriginalVar]])));
                                    if(varsign && coefs[i]>1e-7 && g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==0){
                                        explanation.push_back(g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                        assert(var!=v);
                                    }
                                    if(varsign && coefs[i]<-1e-7 && g_memoscip.VarLb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==1){
                                        if(g_options.encoding==0 || g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar].empty())
                                            explanation.push_back(g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                        else{
                                            for(size_t j = 0; j < g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar].size(); j++){
                                                explanation.push_back(g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar][j]);
                                            }
                                        }
                                        assert(var!=v);
                                    }
                                    if(!varsign && coefs[i]<-1e-7 && g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==0){
                                        explanation.push_back(g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                        assert(var!=v);
                                    }
                                    if(!varsign && coefs[i]>1e-7 && g_memoscip.VarLb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==1){
                                        if(g_options.encoding==0 || g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar].empty())
                                            explanation.push_back(g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                        else{
                                            for(size_t j = 0; j < g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar].size(); j++){
                                                explanation.push_back(g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar][j]);
                                            }
                                        }
                                        assert(var!=v);
                                    }
                                }
                            }
                            if(g_options.verbosity > 3)
                                cout<<"Explanation LB for "<<SCIPvarGetName(var)<<" "<<explanation.size()<<endl;
                            if(g_options.encoding==0 || g_memoscip.TupleAuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar].empty())
                                g_memoscip.nodeInfoMap[nodenumber].BPRemoval.push_back({g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar], explanation});
                            else{
                                for(size_t j = 0; j < g_memoscip.TupleAuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar].size(); j++){
                                    g_memoscip.nodeInfoMap[nodenumber].BPRemoval.push_back({g_memoscip.TupleAuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar][j], explanation});
                                }
                            }
                        }
                    }
                    break;
                }
                case SCIP_BOUNDCHGTYPE_PROPINFER: {
                    SCIP_PROP* prop = SCIPbdchginfoGetInferProp(reason);
                    int inferinfo = SCIPbdchginfoGetInferInfo(reason);
                    if (prop != NULL) {
                        std::string propname = SCIPpropGetName(prop);
                        bool fromRedcost = (propname == "redcost") || (propname == "rootredcost");
                        if(fromRedcost){
                            SCIP_Longint ParentNumber=g_memoscip.nodeInfoMap[nodenumber].parent;
                            if(g_options.verbosity > 2){
                                cout<<"PROPAGATION (" + propname + ", info=" + std::to_string(inferinfo) + "): reduced-cost filtering"<<endl;
                                cout<<"Parent node "<<ParentNumber<<" lb "<<g_memoscip.nodeInfoMap[ParentNumber].lb<<" redcost "<<g_memoscip.nodeInfoMap[ParentNumber].redcosts[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]<<" new lb "<<SCIPnodeGetLowerbound(node)<<" "<<SCIPgetPrimalbound(scip)<<endl;
                            }
                            if(g_options.encoding==0 || g_memoscip.TupleAuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar].empty())
                                g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.push_back({g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar],SCIPgetPrimalbound(scip)-g_memoscip.nodeInfoMap[ParentNumber].lb});
                            else{
                                for(size_t j = 0; j < g_memoscip.TupleAuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar].size(); j++){
                                    g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.push_back({g_memoscip.TupleAuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar][j],SCIPgetPrimalbound(scip)-g_memoscip.nodeInfoMap[ParentNumber].lb});
                                }
                            }
                            assert(SCIPvarGetName(var)[0]!='s' && SCIPvarGetName(var)[2]!='s');
                        }
                    }
                }
                default:
                    //cout<<"OTHER"<<endl;
                    break;
            }
        }else{
            //cout<<"no reason found"<<endl;
        }
        if(g_options.verbosity > 1)
            cout << "[" << reason << "] "<< SCIPvarGetName(var) << " LB: "<< oldbound << " → " << newbound << endl;
    }
    
    if (eventtype & SCIP_EVENTTYPE_UBCHANGED) {
        if(g_options.verbosity > 1)
            cout<<"UB CHANGED "<<SCIPvarGetName(var)<<" "<<oldbound<<" ->"<<newbound<<endl;

        bool auxvarassigned=false;
        if(g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar] != -1)
            auxvarassigned =  g_memoscip.VarUb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]] == g_memoscip.VarLb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]];
        reason = getBoundChangeReasonDetailed(scip, var, SCIP_BOUNDTYPE_UPPER, oldbound, newbound);
        if(reason!=NULL && !auxvarassigned){
            SCIP_BOUNDCHGTYPE chgtype = SCIPbdchginfoGetChgtype(reason);
            //cout<<chgtype<<endl;
            switch (chgtype) {
                case SCIP_BOUNDCHGTYPE_BRANCHING:
                    //cout<< "BRANCHING"<<endl;
                    break;
                    
                    case SCIP_BOUNDCHGTYPE_CONSINFER: {
                        SCIP_CONS* cons = SCIPbdchginfoGetInferCons(reason);
                        if (cons != NULL) {
                            std::string constype = SCIPconshdlrGetName(SCIPconsGetHdlr(cons));
                            std::string consname = SCIPconsGetName(cons);
                            if(g_options.verbosity > 2)
                                cout<< "CONSTRAINT (" + constype + ": " + consname + ")"<<endl;
                            //Loop over variables in constraint and coef
                            SCIP_VAR** vars = SCIPgetVarsLinear(scip,cons);
                            int nvars = SCIPgetNVarsLinear(scip,cons);
                            SCIP_Real* coefs = SCIPgetValsLinear(scip,cons);
                            if(nvars>2){                                             
                                vector<int> explanation;
                                bool varsign;
                                for(int i=0; i < nvars; i++) {
                                    SCIP_VAR* v = vars[i];
                                    if(v == var){
                                        varsign = coefs[i] > 0;
                                        break;
                                    }
                                }
                                for(int i = 0; i < nvars; i++) {
                                    SCIP_VAR* v = vars[i];
                                    if(v != var){
                                        assert(g_memoscip.VarLb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==SCIPvarGetLbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]])));
                                        assert(g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==SCIPvarGetUbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]])));
                                        if(varsign && coefs[i]<-1e-7 && g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==0){
                                            explanation.push_back(g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                            assert(var!=v);
                                            if(g_options.verbosity > 3 )
                                                cout<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
                                        }
                                        if(varsign && coefs[i]>1e-7 && g_memoscip.VarLb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==1){
                                            if(g_options.encoding==0 || g_memoscip.TupleAuxVar.empty())
                                                explanation.push_back(g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                            else{
                                                for(size_t j = 0; j < g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar].size(); j++){
                                                    explanation.push_back(g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar][j]);
                                                }
                                            }
                                            assert(var!=v);
                                            if(g_options.verbosity > 3 )
                                                cout<<"Auxvar "<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
                                        }


                                        if(!varsign && coefs[i]<-1e-7 && g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==0){
                                            explanation.push_back(g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                            assert(var!=v);
                                            if(g_options.verbosity > 3 )
                                                cout<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
                                        }
                                        if(!varsign && coefs[i]>1e-7 && g_memoscip.VarLb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]]==1){
                                            if(g_options.encoding==0 || g_memoscip.TupleAuxVar.empty())
                                                explanation.push_back(g_memoscip.AuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar]);
                                            else{
                                                for(size_t j = 0; j < g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar].size(); j++){
                                                    explanation.push_back(g_memoscip.TupleAuxVar[SCIPvarGetIndex(v)-g_memoscip.nOriginalVar][j]);
                                                }
                                            }
                                            assert(var!=v);
                                            if(g_options.verbosity > 3 )
                                                cout<<"Auxvar "<<SCIPvarGetName(v)<<" "<<coefs[i]<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
                                        }
                                    }
                                }
                                if(g_options.verbosity > 3)
                                    cout<<"ExplanationUB for "<<SCIPvarGetName(var)<<" "<<explanation.size()<<" "<<g_memoscip.nodeInfoMap[nodenumber].BPRemoval.size()<<endl;
                                
                                g_memoscip.nodeInfoMap[nodenumber].BPRemoval.push_back({g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar], explanation});
                            }
                        } else {
                            //cout << "CONSTRAINT (unknown)" << endl;
                        }
                        break;
                    }
                
                case SCIP_BOUNDCHGTYPE_PROPINFER: {
                    SCIP_PROP* prop = SCIPbdchginfoGetInferProp(reason);
                    int inferinfo = SCIPbdchginfoGetInferInfo(reason);
                    if (prop != NULL) {
                        std::string propname = SCIPpropGetName(prop);
                        bool fromRedcost = (propname == "redcost") || (propname == "rootredcost");
                        if(fromRedcost){
                            SCIP_Longint ParentNumber=g_memoscip.nodeInfoMap[nodenumber].parent;
                            if(g_options.verbosity > 2){
                                cout<<"PROPAGATION (" + propname + ", info=" + std::to_string(inferinfo) + "): reduced-cost filtering"<<endl;
                                cout<<"Parent node "<<ParentNumber<<" lb "<<g_memoscip.nodeInfoMap[ParentNumber].lb<<" redcost "<<g_memoscip.nodeInfoMap[ParentNumber].redcosts[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]<<" new lb "<<SCIPnodeGetLowerbound(node)<<" "<<SCIPgetPrimalbound(scip)<<endl;
                            }
                            assert(SCIPvarGetName(var)[0]!='s' && SCIPvarGetName(var)[2]!='s');
                            g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.push_back({g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar],SCIPgetPrimalbound(scip)-g_memoscip.nodeInfoMap[ParentNumber].lb});
                        }
                    }
                    else {
                        cout<<"PROPAGATION (unknown, info=" + std::to_string(inferinfo) + ")"<<endl;
                    }
                    break;
                }
                
                default:
                    //cout<<"OTHER"<<endl;
                    break;
            }
        }
        if(g_options.verbosity > 1 )
            cout << "[" << reason << "] "<< SCIPvarGetName(var) << " UB: "<< oldbound << " → " << newbound << endl;
        if(SCIPvarGetName(var)[0]=='s' || SCIPvarGetName(var)[2]=='s'){
            //cout<<"SLack var"<<endl;
            g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=newbound;
        }else{
            //cout<<SCIPinfinity(scip) <<" "<<newbound<<endl;
            if(newbound==SCIPinfinity(scip) || newbound==1){
                g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=1;
                if(g_options.encoding==0 || g_memoscip.TupleAuxVar.empty())
                    g_memoscip.VarLb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=0;
            }
            else{
                g_memoscip.VarUb[g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=newbound;
                assert(newbound==0);
                if(g_options.encoding==0 || g_memoscip.TupleAuxVar.empty())
                    g_memoscip.VarLb[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]=1;
            }
        }
        
    }
    //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetLbGlobal(var)<<" "<<SCIPvarGetUbGlobal(var)<<endl;
    //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetLbLocal(var)<<" "<<SCIPvarGetUbLocal(var)<<endl;
    return SCIP_OKAY;
}




// Callback d'initialisation - ATTACHER les événements aux variables
static SCIP_DECL_EVENTINIT(eventInitBoundChanged)
{
    SCIP_VAR** vars = SCIPgetVars(scip);
    int nvars = SCIPgetNVars(scip);
    g_boundChangedEventVars.clear();
    
    //std::cout << "Attaching bound change events to " << nvars << " variables" << std::endl;
    
    for (int i = 0; i < nvars; i++) {
        // ✅ Utiliser SCIPcatchVarEvent (PAS SCIPcatchEvent!)
        SCIP_CALL( SCIPcatchVarEvent(scip, vars[i],
                                     SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
                                     eventhdlr, 
                                     NULL,    // eventdata
                                     NULL) ); // filterpos
                        g_boundChangedEventVars.push_back(vars[i]);
    }
    
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXEC(eventExecNodeInfeasible)
{
    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    if( eventtype != SCIP_EVENTTYPE_NODEINFEASIBLE )
        return SCIP_OKAY;

    SCIP_NODE* node = SCIPeventGetNode(event);
    if(node == NULL){
        return SCIP_OKAY;
    }
    //cout<<"-----NODE INFEASIBLE EVENT -------"<<endl;
    
    //SCIP_Longint nodenumber= SCIPnodeGetNumber(node);

    /*if( g_memoscip.nodeInfoMap[nodenumber].Status!="INFEASIBLE"){
        cout<<"Node "<<nodenumber<<" is infeasible"<<endl;
        g_memoscip.nodeInfoMap[nodenumber].Status="INFEASIBLE";
        SCIP_CALL(gatherLPInfeasibilityCertificate(scip, nodenumber));
        //SCIP_CALL(SCIPwriteLP(scip, "infeasible.lp"));
        //SCIP_CALL(SCIPwriteTransProblem(scip, "infeasibleproblem.lp","lp",FALSE));
    }*/
    return SCIP_OKAY;
}

static
SCIP_DECL_BRANCHEXECLP(branchExecLP)
{
   // If you tightened domains, let SCIP re-solve LP first
    if(g_options.verbosity > 1)
        cout<<"branching on fractional variables"<<endl;
    SCIP_NODE* node = SCIPgetCurrentNode(scip);
    SCIP_Longint nodenumber=SCIPnodeGetNumber(node);

    (void)nodenumber;

    if( *result == SCIP_REDUCEDDOM )
        return SCIP_OKAY;

    if (g_options.branchingHeuristic == 4)
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

    if(g_options.branchingHeuristic < 1){
        //cout<<"Using default branching heuristic" << endl;
        // Branch on the first fractional variable (in g_memoscip.CurrentVarIdx order)
        for( size_t i = 0; i < g_memoscip.OrderVar.size(); ++i )
        {
            SCIP_VAR* v = g_memoscip.OrderVar[i];
            // Only branch on discrete vars
            if( SCIPvarGetType(v) == SCIP_VARTYPE_CONTINUOUS )
                continue;

            SCIP_Real x = SCIPgetSolVal(scip, nullptr, v);  // current LP solution value
            //cout << "Current LP solution value for var " << SCIPvarGetName(v) << " is " << x << endl;
            if( !SCIPisFeasIntegral(scip, x) )
            {
                SCIP_NODE* downchild = nullptr;
                SCIP_NODE* eqchild   = nullptr;
                SCIP_NODE* upchild   = nullptr;

                SCIP_CALL( SCIPbranchVarVal(scip, v, 1, &downchild, &eqchild, &upchild) );
                *result = SCIP_BRANCHED;

                // optional debug
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



static
SCIP_DECL_EVENTFREE(eventFreeInf)
{
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTINIT(eventInitInf)
{
    SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEINFEASIBLE, eventhdlr, nullptr, nullptr));
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXIT(eventExitInf)
{
    SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEINFEASIBLE, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}

// Event handler free callback (existing)
static
SCIP_DECL_EVENTFREE(eventFreeLP)
{
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

static
SCIP_DECL_EVENTFREE(eventFreeBestSolFound)
{
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTINIT(eventInitBestSolFound)
{
    SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, nullptr, nullptr));
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXIT(eventExitBestSolFound)
{
    SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}

// Event handler init callback (existing)
static
SCIP_DECL_EVENTINIT(eventInitLP)
{
    SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, nullptr, nullptr));
    return SCIP_OKAY;
}

// Event handler exit callback (existing)
static
SCIP_DECL_EVENTEXIT(eventExitLP)
{
    SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}



static
SCIP_DECL_EVENTFREE(eventFreeNode)
{
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTINIT(eventInitNode)
{
    SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, nullptr, nullptr));
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXIT(eventExitNode)
{
    SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}


static
SCIP_DECL_EVENTFREE(eventFreeDel)
{
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTINIT(eventInitDel)
{
    SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEDELETE, eventhdlr, nullptr, nullptr));
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXIT(eventExitDel)
{
    SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEDELETE, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTFREE(eventFreeBoundChanged)
{
    //SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, nullptr, -1));
    return SCIP_OKAY;
}

static SCIP_DECL_EVENTEXIT(eventExitBoundChanged)
{
    for (SCIP_VAR* var : g_boundChangedEventVars) {
        // ✅ Utiliser SCIPdropVarEvent (PAS SCIPdropEvent!)
        SCIP_CALL( SCIPdropVarEvent(scip, var,
                                    SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
                                    eventhdlr, 
                                    NULL,    // eventdata
                                    -1) );   // filterpos (-1 = auto)
    }
    g_boundChangedEventVars.clear();
    
    return SCIP_OKAY;
}

void includeBoundChangedEventHandler(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;

    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(scip, &eventhdlr,
        "boundchanged_event", "Capture bound changes",
        eventExecBoundChanged, nullptr));
    
    SCIP_CALL_ABORT(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeBoundChanged));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitBoundChanged));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitBoundChanged));
}

// Include the node-focused event handler
void includeNodeFocusedEventHandler(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;

    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(scip, &eventhdlr,
        "nodefocused_event", "Capture LB increases on node focus",
        eventExecNodeFocused, nullptr));
    
    SCIP_CALL_ABORT(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeNode));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitNode));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitNode));
}

void includeNodeInfeasibleEventHandler(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;

    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(scip, &eventhdlr,
        "nodeinf_event", "detect infeasible nodes",
        eventExecNodeInfeasible, nullptr));

    SCIP_CALL_ABORT(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeInf));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitInf));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitInf));
}

// Include the LP event handler (existing)
void includeLPSolveEventHandler(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;

    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(scip, &eventhdlr,
        "lpsolve_event", "Print duals after each LP solve",
        eventExecLPSolved, nullptr));

    // Set the callbacks - these will be called automatically by SCIP
    SCIP_CALL_ABORT(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeLP));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitLP));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitLP));
}

void includeBestSolFoundEventHandler(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;

    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(scip, &eventhdlr,
        "bestsolfound_event", "Notify when a new incumbent integer solution is found",
        eventExecBestSolFound, nullptr));

    SCIP_CALL_ABORT(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeBestSolFound));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitBestSolFound));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitBestSolFound));
}

void includeNodeDeleteEventHandler(SCIP* scip)
{
    SCIP_EVENTHDLR* eventhdlr = nullptr;

    SCIP_CALL_ABORT(SCIPincludeEventhdlrBasic(
        scip, &eventhdlr,
        "event_nodedelete",
        "handles node deletion events for bottom-up aggregation",
        eventExecNodeDelete, nullptr));

    SCIP_CALL_ABORT(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeDel));
    SCIP_CALL_ABORT(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitDel));
    SCIP_CALL_ABORT(SCIPsetEventhdlrExit(scip, eventhdlr, eventExitDel));
}

void register_memoscip_branchrule(SCIP* scip)
{
    SCIP_CALL_ABORT(SCIPincludeBranchrule(scip,
        "mybranchrule",
        "My custom branch rule",
        6000000,
        -1,
        1.0,
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
