#include "memoscip_events.h"

#include "memoscip_dual.h"
#include "memoscip_option.h"
#include "memoscip_types.h"
#include "sepa_mycuts.h"

#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/scip_lp.h"
#include "scip/scip_prob.h"
#include "scip/var.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

static
SCIP_DECL_EVENTEXEC(eventExecNodeSolved)
{
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
    SCIP_Bool cutoff, lperror;
    SCIP_LPSOLSTAT stat = SCIPgetLPSolstat(scip);


    if(g_memoscip.nodeInfoMap[nodenumber].Status=="Unknown" || stat==SCIP_LPSOLSTAT_NOTSOLVED){   
        if(stat==SCIP_LPSOLSTAT_NOTSOLVED){
            SCIP_CALL( SCIPconstructLP(scip,&cutoff));
            SCIP_CALL( SCIPstartProbing(scip) );
            SCIP_CALL( SCIPflushLP(scip) );
            SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) ); 
            SCIP_CALL( SCIPwriteTransProblem(scip, "infeasibleproblem.lp","lp",FALSE) );
 
            if(stat==SCIP_LPSOLSTAT_INFEASIBLE){
                SCIP_CALL(gatherLPInfeasibilityCertificate(scip, nodenumber));
                //SCIP_CALL( SCIPwriteLP(scip, "infeasible.lp") );
                lpobj=100000;
            }else{
                lpobj=SCIPgetLPObjval(scip);
                //cout<<"lpobj "<<lpobj<<" parent lb "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb<<endl;
                if(!SCIPisDualfeasZero(scip,lpobj) && (fabs(g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb - lpobj)>= g_options.minInc || ceil(lpobj)>=SCIPgetPrimalbound(scip) || lpobj - ceil(g_memoscip.nodeInfoMap[parentnumber].lb) > 1e-5))
                    AnalyseDualSol(scip,nodenumber );
            }
            SCIP_CALL( SCIPendProbing(scip));

        }else{
            lpobj=SCIPgetLPObjval(scip);
            //cout<<"lpobj "<<lpobj<<" parent lb "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb<<endl;
            if(!SCIPisDualfeasZero(scip,lpobj) && (fabs(g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodenumber].parent].lb - lpobj)>= g_options.minInc || ceil(lpobj)>=SCIPgetPrimalbound(scip) || lpobj - ceil(g_memoscip.nodeInfoMap[parentnumber].lb) > 1e-5))
                AnalyseDualSol(scip,nodenumber );
        }
        //cout<<"lb "<<SCIPgetLPObjval(scip)<<" ceil "<<ceil(SCIPgetLPObjval(scip))<<" ub"<<SCIPgetPrimalbound(scip)<<" "<<ceil(SCIPgetLPObjval(scip))-SCIPgetPrimalbound(scip)<<endl;

        //assert(ceil(lpobj)-SCIPgetPrimalbound(scip)>=-1e-9);
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

    //print redcosts.
    
    /*for (size_t i = 0; i < g_memoscip.nodeInfoMap[nodenumber].redcosts.size(); i++)
    {
        if(!SCIPisDualfeasZero(scip, g_memoscip.nodeInfoMap[nodenumber].redcosts[i]))
            cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" "<<g_memoscip.nodeInfoMap[nodenumber].redcosts[i]<<" "<<endl;
    }*/
    
    
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

    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    if( eventtype != SCIP_EVENTTYPE_LPSOLVED )
        return SCIP_OKAY;

    SCIP_NODE* node = SCIPeventGetNode(event);
    if(node == NULL){
        return SCIP_OKAY;
    }
    SCIP_Longint nodeNumber= SCIPnodeGetNumber(node);
    SCIP_Longint parentNumber = g_memoscip.nodeInfoMap[nodeNumber].parent;
    if(g_options.verbosity > 1)
        cout<<"-----LP SOLVED EVENT ------- "<<nodeNumber<< " depth "<<SCIPgetDepth(scip)<<" parent "<<g_memoscip.nodeInfoMap[nodeNumber].parent<<endl;
    SCIP_LPSOLSTAT lpsolstat = SCIPgetLPSolstat(scip);
    if(g_memoscip.CurrentVarIdx.empty()){
        int nconss = SCIPgetNConss(scip);
        SCIP_CONS** conss = SCIPgetConss(scip);
        static unordered_map<SCIP_VAR*, SCIP_VAR*> TempAuxVar;
        for(int i = 0; i < nconss; ++i)
        {  
            SCIP_CONS* cons = conss[i];

            // Get variables in this constraint
            SCIP_VAR** vars = SCIPgetVarsLinear(scip,cons);
            int nvars = SCIPgetNVarsLinear(scip,cons);
            if(g_memoscip.FileExt==".opb" && SCIPconsGetName(cons)[0]=='u' and SCIPconsGetName(cons)[1]=='b' ){
                TempAuxVar[vars[0]]=vars[1];
                TempAuxVar[vars[1]]=vars[0];
            }
            else if((g_memoscip.FileExt==".lp" || g_options.kpcg) && SCIPgetRhsLinear(scip, cons)==1 && nvars==2){
                TempAuxVar[vars[0]]=vars[1];
                TempAuxVar[vars[1]]=vars[0];
            }
        }

        int ncols = SCIPgetNLPCols(scip);
        SCIP_COL** cols = SCIPgetLPCols(scip);
        for( int i = 0; i < ncols; i++ )
        {
            SCIP_COL* col = cols[i];
            SCIP_VAR* var = SCIPcolGetVar(col);  // ← variable correspondante            
            g_memoscip.AllVars.push_back(var);
            if(SCIPvarGetName(var)[2]=='s' && SCIPvarGetName(var)[3]=='c')
                g_memoscip.slackVar.push_back(i);
            else{
                g_memoscip.eqvid[SCIPvarGetIndex(var)]=g_memoscip.CurrentVarIdx.size();
                g_memoscip.CurrentVarIdx.push_back(i);
                g_memoscip.VarLb.push_back(SCIPvarGetLbLocal(var));

                if(SCIPvarGetName(var)[2]!='s'){
                    if(SCIPvarGetUbLocal(var)==SCIPinfinity(scip))
                        g_memoscip.VarUb.push_back(1);
                    else
                        g_memoscip.VarUb.push_back(SCIPvarGetUbLocal(var));
                }else {
                    //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetUbLocal(var)<<" "<<SCIPinfinity(scip)<<endl;
                    g_memoscip.VarUb.push_back(SCIPvarGetUbLocal(var));
                }
            }
        }

        for(size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); i++){
            if(TempAuxVar.contains(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])){
                g_memoscip.AuxVar.push_back(g_memoscip.eqvid[SCIPvarGetIndex(TempAuxVar[g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]]])]);
            }
            else
                g_memoscip.AuxVar.push_back(i);
        }       
    }


    if(g_options.verbosity > 1){
        cout << "LP Solved - Status: ";
        switch(lpsolstat) {
            case SCIP_LPSOLSTAT_OPTIMAL:
                cout << "OPTIMAL "<<SCIPgetLPObjval(scip)<<" "<<SCIPgetPrimalbound(scip)<<" node status "<<g_memoscip.nodeInfoMap[nodeNumber].Status<<endl;
                break;
            case SCIP_LPSOLSTAT_INFEASIBLE:
                cout << "INFEASIBLE";
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
    if(g_options.learning > 0){
        if(lpsolstat==SCIP_LPSOLSTAT_INFEASIBLE){
            g_memoscip.nodeInfoMap[nodeNumber].Status="INFEASIBLE";
            SCIP_CALL(gatherLPInfeasibilityCertificate(scip, nodeNumber));
        }
        else {
            SCIP_Real lpobj = SCIPgetLPObjval(scip);
            //We analyse the dual solution only if the lower bound increased by at least minInc or if the lower bound is greater than the upper bound.
            if(fabs(g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb - lpobj) >= g_options.minInc || ceil(lpobj) >= SCIPgetPrimalbound(scip) || lpobj - ceil(g_memoscip.nodeInfoMap[parentNumber].lb)> 1e-5)
                AnalyseDualSol(scip, nodeNumber);
            else{
                //If lb didn't change then we copy the parent information
                g_memoscip.nodeInfoMap[nodeNumber].redcosts = g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].redcosts;
                g_memoscip.nodeInfoMap[nodeNumber].lb=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb;
                g_memoscip.nodeInfoMap[nodeNumber].y_i=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i;
                g_memoscip.nodeInfoMap[nodeNumber].DualSol=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol;
                g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=0;
                g_memoscip.nodeInfoMap[nodeNumber].Status="FEASIBLE"; 
                g_memoscip.nodeInfoMap[nodeNumber].lb_increase=0;
            }
        }
    }
    
    //cout<<"END LP SOLVED EVENT "<<nodeNumber<<" "<<g_memoscip.nodeInfoMap[nodeNumber].Status<<" "<<g_memoscip.nodeInfoMap[nodeNumber].lb<<" "<<SCIPgetPrimalbound(scip)<<endl;
    return SCIP_OKAY;
}


// ---- New: NODEFOCUSED event handler ----
static
SCIP_DECL_EVENTEXEC(eventExecNodeFocused)
{    

    SCIP_EVENTTYPE eventtype = SCIPeventGetType(event);
    if( eventtype != SCIP_EVENTTYPE_NODEFOCUSED )
        return SCIP_OKAY;

    SCIP_NODE* node = SCIPeventGetNode(event);
    if(node == NULL){
        return SCIP_OKAY;
    }

    SCIP_Longint nodeNumber=SCIPnodeGetNumber(node);
    SCIP_NODE* parent = nullptr;
    SCIP_Longint parentNumber=0;
    NodeInfo &ninfo = g_memoscip.nodeInfoMap[nodeNumber]; // creates if not existing
    ninfo.visited = true;
    //ninfo.lb = node_lb;
    ninfo.number = nodeNumber;
    ninfo.ub = SCIPgetPrimalbound(scip);
    if(nodeNumber!=1)  // skip root
    {
        parent = SCIPnodeGetParent(node);
        parentNumber = SCIPnodeGetNumber(parent);
        if(g_options.verbosity > 1)
            cout<<"-----------FOCUSING NODE "<<nodeNumber<<" depth "<<SCIPgetDepth(scip)<<" with Parent "<<parentNumber<<"---------------------and "<<g_memoscip.nbAddedctr<<" CONSTRAINT added "<<endl<<"UB :"<< SCIPgetPrimalbound(scip)<<endl<<" LB : "<<SCIPnodeGetLowerbound(node)<<endl;
        g_memoscip.nodeInfoMap[nodeNumber].parent=parentNumber;

    }else{
        //CHECK ROOT NODE ASSIGNMENT 
    }
   
    SCIP_DOMCHG* domchg;
    SCIP_Bool cutoff;
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
        //cout<<"At node "<<nodeNumber<<" Bound change for var "<<SCIPvarGetName(var)<<" to "<<newbound<<" with obj coef "<<SCIPvarGetObj(var)<<endl;
        //cout<<"Global lb of var "<<SCIPvarGetName(var)<<" is "<<SCIPvarGetLbGlobal(var)<< " "<< SCIPvarGetUbGlobal(var)<<endl;
        if(newbound==1.0){
            ninfo.BranchVar = g_memoscip.eqvid[SCIPvarGetIndex(var)];
            
        }else{
            ninfo.BranchVar = g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]];
        }        
    }
    assert(g_memoscip.VarLb.size()==g_memoscip.VarUb.size());

    if(g_options.debug){
        for (size_t i = 0; i < g_memoscip.VarLb.size(); i++)
        {   
            g_memoscip.nodeInfoMap[nodeNumber].VarLbDebug.push_back(g_memoscip.VarLb[i]);
            g_memoscip.nodeInfoMap[nodeNumber].VarUbDebug.push_back(g_memoscip.VarUb[i]);
            SCIP_VAR* var= g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]];
            //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetLbLocal(var)<<" "<<g_memoscip.VarLb[i]<<endl;
            //cout<<SCIPvarGetName(var)<<" "<<SCIPvarGetUbLocal(var)<<" "<<g_memoscip.VarUb[i]<<endl;
            assert(SCIPvarGetLbLocal(var)==g_memoscip.VarLb[i] || SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[i]]])==0);
            assert(SCIPvarGetUbLocal(var)==g_memoscip.VarUb[i] ||  SCIPvarGetUbLocal(var)==SCIPinfinity(scip) || SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[i]]])==1);
        }
        
        //SCIPprintQueuedCuts(scip,&g_memoscip.CutQueue);
        for (int i = (int)g_memoscip.CutQueue.queue.size()-1; i >= 0; i--)
        {
            //cout<<SCIPvarGetName(g_memoscip.CutQueue.queue[i].vars.back())<<" "<<SCIPvarGetProbindex(g_memoscip.CutQueue.queue[i].vars.back())<<" "<<SCIPvarIsInLP(g_memoscip.CutQueue.queue[i].vars.back())<<endl;
            //assert(SCIPvarGetProbindex(g_memoscip.CutQueue.queue[i].vars.back())!=-1);
            assert(SCIPvarIsInLP(g_memoscip.CutQueue.queue[i].vars.back()));
            //assert(SCIPvarGetName(g_memoscip.CutQueue.queue[i].vars.back())[2]=='s' && SCIPvarGetName(g_memoscip.CutQueue.queue[i].vars.back())[3]=='c');
            //g_memoscip.CurrentVarIdx.push_back(g_memoscip.slackVar[g_memoscip.nbAddedctr-i]);
        }
    }
    int qsize=g_memoscip.CutQueue.queue.size();
    if(qsize>0){
        SCIP_CALL(AddCutstoProb(scip,&g_memoscip.CutQueue,double(g_memoscip.TOP)));    
    }
    SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
    SCIP_CALL(SCIPflushLP(scip));
    if( g_options.output>1 && nodeNumber >= g_options.output){
        char nameprobsol[50];

        snprintf(nameprobsol, 50, "solved_%lld.lp", nodeNumber);
        SCIP_CALL( SCIPwriteLP(scip, nameprobsol) ); 
        snprintf(nameprobsol, 50, "problem_%lld.lp", nodeNumber);
        SCIP_CALL( SCIPwriteTransProblem(scip, nameprobsol,"lp",FALSE) );
    }
    
    // link this node as child of its parent (if parent exists)
    if (parent != nullptr)
    {
        // add child pointer to parent's children list if not already present
        g_memoscip.nodeInfoMap[parentNumber].children.push_back(SCIPnodeGetNumber(node));
        //cout<<SCIPnodeGetLowerbound(node)<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb -SCIPnodeGetLowerbound(node)<<" "<<endl;
        //assert(fabs(g_memoscip.nodeInfoMap[parentNumber].lb -SCIPnodeGetLowerbound(node)) <= 1 || ceil(g_memoscip.nodeInfoMap[parentNumber].lb) - SCIPgetPrimalbound(scip)>= -1e-9 || g_memoscip.nodeInfoMap[parentNumber].Status=="INFEASIBLE");
    }
    
    if(g_options.verbosity > 1)
        cout<<"-------------------------------------ENDOFCUS EVENT ----------------"<<endl;
    
    return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXEC(eventExecNodeDelete)
{
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
    if(g_options.verbosity > 1){
        cout<<"Deleting node "<<nodenumber<<" from nodeInfoMap, scip Node Lb "<< SCIPgetNodeLowerbound(scip,node)  <<" my node lb : "<<g_memoscip.nodeInfoMap[nodenumber].lb<< " "<<" depth : "<<SCIPgetDepth(scip)<<endl;
        cout<<"UB "<<SCIPgetPrimalbound(scip)<<endl;
        cout<<"number of children "<<g_memoscip.nodeInfoMap[nodenumber].children.size()<<endl;
        cout<<"Bracnh var "<<g_memoscip.nodeInfoMap[nodenumber].BranchVar<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[nodenumber].BranchVar]])<<endl;
        cout<<"Curnnet node lb increase "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" versus children "<<endl;
        
    }
    SCIP_Longint parentNumber = g_memoscip.nodeInfoMap[nodenumber].parent;

    //cout<< "NodeInfoMap size before deletion: " << g_memoscip.nodeInfoMap.size() << endl;
    //SCIPinfoMessage(scip, NULL, "[DEBUG] Total LP rows now: %d\n", SCIPgetNLPRows(scip));
    //Print children of the node being deleted
    if (g_memoscip.nodeInfoMap.contains(nodenumber) && g_memoscip.nbAddedctr<(int)g_memoscip.slackVar.size()) {
        vector<long> children = g_memoscip.nodeInfoMap[nodenumber].children;
        vector<SCIP_Real> ResFus(g_memoscip.CurrentVarIdx.size(),0.0);
        SCIP_Real rhsFus=0;
        if(children.size()>0){
            if(children.size()==1){
                
                assert(g_memoscip.nodeInfoMap.contains(children[0]));
                if(g_options.verbosity > 2)
                    cout<<"Single child "<<children[0]<<" lb "<<g_memoscip.nodeInfoMap[children[0]].lb<<" versus parent lb "<<g_memoscip.nodeInfoMap[nodenumber].lb<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.first.size()<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second<<" incr lb "<<g_memoscip.nodeInfoMap[children[0]].lb_increase<<endl;
                assert(g_memoscip.nodeInfoMap[children[0]].lb<g_memoscip.TOP);
                assert(fabs(g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second-g_memoscip.nodeInfoMap[children[0]].lb_increase)<1e-5);
            }
            if(children.size()==2){
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

                    assert(g_memoscip.nodeInfoMap[children[1]].BranchVar==g_memoscip.AuxVar[g_memoscip.nodeInfoMap[children[0]].BranchVar]);
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
            }
            assert(ResFus.size()==g_memoscip.CurrentVarIdx.size());
            //cout<<"Result of Fusion Resolution: rhs "<< rhsFus<<endl;
            //cout<<g_memoscip.nodeInfoMap[nodenumber].redcosts.size()<<"' "<<ResFus.size()<<endl;

            //PRINT RESFUS
            if(g_options.verbosity > 2){
                cout<<"ResFus "<<endl;
                for (size_t i = 0; i < ResFus.size(); i++)
                {
                    if(!SCIPisDualfeasZero(scip,ResFus[i])){
                        if(ResFus[i]>0)
                            cout<<"+"<<ResFus[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
                        else
                            cout<<ResFus[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
                    }

                }
                cout<<" >= "<<rhsFus<<endl;
            }
            SCIP_Real rhs_summed=rhsFus;
            vector<int> NewInds;
            vector<SCIP_Real> NewCoefs;
            vector<SCIP_VAR*> NewVar;

            if(rhs_summed>0 && g_options.learning==2){
                 char name[100];

                snprintf(name, 100, "cut_%lld_%d", nodenumber,SCIPgetDepth(scip));
                

                for (size_t i = 0; i < ResFus.size(); i++)
                {
                    if(!SCIPisDualfeasZero(scip,ResFus[i])){
                        NewInds.push_back(i);
                        NewCoefs.push_back(ResFus[i]);
                        NewVar.push_back(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]]);
                        //assert(ResFus[i]<g_memoscip.TOP);
                    }
                }

                SCIP_VAR* slack=g_memoscip.AllVars[g_memoscip.slackVar[g_memoscip.nbAddedctr]];     
                NewVar.push_back(slack);
                NewInds.push_back(g_memoscip.CurrentVarIdx.size());
                NewCoefs.push_back(-1.0);

                assert(NewVar.size()==NewCoefs.size());
                assert(NewVar.size()!=0);
                SCIP_CALL( SCIPqueueCut(scip, name,
                                (int)NewVar.size(),
                                NewVar.data(),
                                NewCoefs.data(),
                                rhs_summed,  // lhs (pour inégalité)
                                rhs_summed,&g_memoscip.CutQueue) ); 


                g_memoscip.AuxVar.push_back(g_memoscip.CurrentVarIdx.size());
                g_memoscip.eqvid[SCIPvarGetIndex(slack)]=g_memoscip.CurrentVarIdx.size();
                g_memoscip.CurrentVarIdx.push_back(g_memoscip.slackVar[g_memoscip.nbAddedctr]);
                g_memoscip.VarUb.push_back(SCIPinfinity(scip));
                g_memoscip.VarLb.push_back(0);

                g_memoscip.nbAddedctr++;


                if(g_options.verbosity > 2)
                    cout<<"Adding cut "<<name<<" with rhs "<<rhs_summed<<" and "<<NewVar.size()<<" variables"<<endl;
            }

            assert(fabs(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second-g_memoscip.nodeInfoMap[nodenumber].lb_increase)<1e-5);
            
            if(g_memoscip.nodeInfoMap[nodenumber].lb_increase>0){
                //cout<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.size()<<endl;
                for (size_t i = 0; i < g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.size(); i++)
                {
                    int idx =g_memoscip.nodeInfoMap[nodenumber].DualProofInds[i];
                    ResFus[idx]+=g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i];
                    //cout<<idx<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[idx]])<<" "<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i]<<" "<<ResFus[idx]<<endl;
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
                    cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[idx]])<< " "<< idx<<" "<<g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second<<" "<<ceil(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second+g_memoscip.nodeInfoMap[nodenumber].lb)<<" "<<SCIPgetPrimalbound(scip)<<endl;
                //assert(ceil(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second+g_memoscip.nodeInfoMap[nodenumber].lb)>=SCIPgetPrimalbound(scip));
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
                }else{
                    //cout<<"Adding cut with rhs "<<rhs_summed<<" vs "<<g_memoscip.nodeInfoMap[children[0]].lb_increase + g_memoscip.nodeInfoMap[nodenumber].lb_increase<<endl;
                    //cout<<"Min children "<<g_memoscip.nodeInfoMap[children[0]].lb_increase <<endl;
                    assert(abs(rhs_summed - (g_memoscip.nodeInfoMap[children[0]].lb_increase + g_memoscip.nodeInfoMap[nodenumber].lb_increase )<1e-9) || !g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.empty()); 

                }
            }
        
            

            //cout<<"rhs summed "<<rhs_summed<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb<<endl;
            //cout<<ceil(rhs_summed+g_memoscip.nodeInfoMap[parentNumber].lb)<<" "<<SCIPgetPrimalbound(scip)<<endl;
            //assert(SCIPgetPrimalbound(scip) == SCIPinfinity(scip)  || ceil(rhs_summed+g_memoscip.nodeInfoMap[parentNumber].lb) - SCIPgetPrimalbound(scip) > -1e-9);

            if(rhs_summed>0 && g_options.learning==1){
                 char name[100];

                snprintf(name, 100, "cut_%lld_%d", nodenumber,SCIPgetDepth(scip));
                

                for (size_t i = 0; i < ResFus.size(); i++)
                {
                    if(!SCIPisDualfeasZero(scip,ResFus[i])){
                        NewInds.push_back(i);
                        NewCoefs.push_back(ResFus[i]);
                        NewVar.push_back(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]]);
                        //cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<< " "<< i<<" "<<ResFus[i]<<" "<<g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<endl;
                        assert(ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<1e-5 || g_memoscip.VarUb[i] == g_memoscip.VarLb[i]);
                        //assert(ResFus[i]<g_memoscip.TOP);
                    }
                }

                SCIP_VAR* slack=g_memoscip.AllVars[g_memoscip.slackVar[g_memoscip.nbAddedctr]];     
                NewVar.push_back(slack);
                NewInds.push_back(g_memoscip.CurrentVarIdx.size());
                NewCoefs.push_back(-1.0);
                g_memoscip.AuxVar.push_back(g_memoscip.CurrentVarIdx.size());
                g_memoscip.eqvid[SCIPvarGetIndex(slack)]=g_memoscip.CurrentVarIdx.size();
                g_memoscip.CurrentVarIdx.push_back(g_memoscip.slackVar[g_memoscip.nbAddedctr]);
                g_memoscip.VarUb.push_back(SCIPinfinity(scip));
                g_memoscip.VarLb.push_back(0);
                if(g_options.verbosity > 2)
                    cout<<"Adding cut "<<name<<" with rhs "<<rhs_summed<<" and "<<NewVar.size()<<" variables"<<endl;

                assert(NewVar.size()==NewCoefs.size());
                assert(NewVar.size()!=0);
                SCIP_CALL( SCIPqueueCut(scip, name,
                                (int)NewVar.size(),
                                NewVar.data(),
                                NewCoefs.data(),
                                rhs_summed,  // lhs (pour inégalité)
                                rhs_summed,&g_memoscip.CutQueue) ); 
               
                g_memoscip.nbAddedctr++;
            }
            
            //cout<<"Dual proof coef second "<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second<<" lb increase "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<endl;
            assert(abs(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second-g_memoscip.nodeInfoMap[nodenumber].lb_increase) < 1e-9);
            assert(NewCoefs.size()<=g_memoscip.CurrentVarIdx.size());
            g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second=rhs_summed;
            g_memoscip.nodeInfoMap[nodenumber].DualProofInds=NewInds;
            g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first=NewCoefs;
            assert(g_memoscip.nodeInfoMap[nodenumber].lb_increase<=rhs_summed);
            g_memoscip.nodeInfoMap[nodenumber].lb_increase=rhs_summed;
            //cout<<" plb "<<g_memoscip.nodeInfoMap[parentNumber].lb<<" lbinc "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" ub "<<SCIPgetPrimalbound(scip)<<endl;
            assert(ceil(g_memoscip.nodeInfoMap[parentNumber].lb + g_memoscip.nodeInfoMap[nodenumber].lb_increase) >= SCIPgetPrimalbound(scip));
            g_memoscip.nodeInfoMap.erase(children[0]);
            if(children.size()==2)  
                g_memoscip.nodeInfoMap.erase(children[1]);
            g_memoscip.nodeInfoMap[nodenumber].children.clear();
           
        }else{

            //Case no child.
            if(g_options.verbosity > 3)
                cout<<"BOund propagation "<<g_memoscip.nodeInfoMap[nodenumber].BPRemoval.size()<<endl;
            for(int i = (int)g_memoscip.nodeInfoMap[nodenumber].BPRemoval.size()-1; i >= 0; i--) {
                //cout<<"removed "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].first])<<" "<<g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second.size()<<endl;
                int idxi = g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].first;
                auto it = find(g_memoscip.nodeInfoMap[nodenumber].DualProofInds.begin(), g_memoscip.nodeInfoMap[nodenumber].DualProofInds.end(), idxi);
                //cout<<"it "<<it - g_memoscip.nodeInfoMap[nodenumber].DualProofInds.begin()<<endl;
                assert(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.size()==g_memoscip.nodeInfoMap[nodenumber].DualProofInds.size());
                int pos = it - g_memoscip.nodeInfoMap[nodenumber].DualProofInds.begin();
                if(it != g_memoscip.nodeInfoMap[nodenumber].DualProofInds.end()){
                    //cout<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[pos]<<endl;
                    if( g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[pos] > 1e-5){
                        assert(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[pos]>1e-5);
                        for(size_t j = 0; j < g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second.size(); j++) {
                            //cout<<"explanation "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second[j]])<<endl;
                            int idxj =g_memoscip.nodeInfoMap[nodenumber].BPRemoval[i].second[j];
                            auto itj = find(g_memoscip.nodeInfoMap[nodenumber].DualProofInds.begin(), g_memoscip.nodeInfoMap[nodenumber].DualProofInds.end(), idxj);
                            if(itj != g_memoscip.nodeInfoMap[nodenumber].DualProofInds.end()){
                                g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[itj - g_memoscip.nodeInfoMap[nodenumber].DualProofInds.begin()]=g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second;
                            }
                            else {
                                g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.push_back(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second);
                                g_memoscip.nodeInfoMap[nodenumber].DualProofInds.push_back(idxj);
                            }
                        }
                    }
                    g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[pos]=0;
                }
            }
            if(g_options.verbosity > 3)
                cout<<"Node consistency "<<g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.size()<<" "<<nodenumber<<endl;
            for (size_t i = 0; i < g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.size(); i++)
            {
                int idx =g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].first;
                // search idx in g_memoscip.nodeInfoMap[nodenumber].DualProofInds
                auto it = find(g_memoscip.nodeInfoMap[nodenumber].DualProofInds.begin(), g_memoscip.nodeInfoMap[nodenumber].DualProofInds.end(), idx);
                if(it != g_memoscip.nodeInfoMap[nodenumber].DualProofInds.end()){
                    g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[it - g_memoscip.nodeInfoMap[nodenumber].DualProofInds.begin()]=g_memoscip.nodeInfoMap[parentNumber].redcosts[idx];
                }else {
                    g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.push_back(g_memoscip.nodeInfoMap[parentNumber].redcosts[idx]);
                    g_memoscip.nodeInfoMap[nodenumber].DualProofInds.push_back(idx);
                }
                g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second=min(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second,g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second);
                
                //cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[idx]])<< " "<< idx<<" "<<g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second<<" "<<g_memoscip.nodeInfoMap[parentNumber].redcosts[idx]<<endl;
                assert(ceil(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second+g_memoscip.nodeInfoMap[nodenumber].lb)>=SCIPgetPrimalbound(scip));
            }
            

            //Print dual proof
            if(g_options.verbosity > 2){
                cout<<"Dual proof after no child "<<g_memoscip.nodeInfoMap[nodenumber].DualProofInds.size()<<endl;
                for(size_t i = 0; i < g_memoscip.nodeInfoMap[nodenumber].DualProofInds.size(); i++)
                {
                    if(!SCIPisDualfeasZero(scip, g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i])){
                        if( g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i]<0){
                            if(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i]>-1e3)
                                cout<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[nodenumber].DualProofInds[i]]])<<" ";
                        }else{
                            cout<<"+"<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[nodenumber].DualProofInds[i]]])<<" ";
                        }
                    }
                }
                cout<<"<= "<<g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second<<endl;
            }
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
        assert(SCIPvarGetName(var)[2]!='s');
        g_memoscip.VarLb[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=newbound;
        if(SCIPvarGetName(var)[2]!='s'){
            if(newbound==1)
                g_memoscip.VarUb[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]]]=0;
            else if(newbound==0)
                g_memoscip.VarUb[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]]]=1;
            else
                g_memoscip.VarUb[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]]]=newbound ;
        }
        reason = getBoundChangeReasonDetailed(scip, var, SCIP_BOUNDTYPE_LOWER, oldbound, newbound);

         
        
        assert(!SCIPinProbing(scip));
        //if reason is COnstraint :
        if(reason!=NULL){
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
                            for(int i = 0; i < nvars; i++) {
                                SCIP_VAR* v = vars[i];
                                if(v != var){
                                    assert(g_memoscip.VarLb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==SCIPvarGetLbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])));
                                    assert(g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==SCIPvarGetUbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])));
                                    if(coefs[i]>1e-7 && g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==0){
                                        explanation.push_back(g_memoscip.eqvid[SCIPvarGetIndex(v)]);
                                        assert(var!=v);
                                    }
                                    if(coefs[i]<-1e-7 && g_memoscip.VarLb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==1){
                                        explanation.push_back(g_memoscip.eqvid[SCIPvarGetIndex(v)]);
                                        assert(var!=v);
                                    }
                                }
                            }
                            //cout<<"Explanation for "<<SCIPvarGetName(var)<<" "<<explanation.size()<<endl;
                            g_memoscip.nodeInfoMap[nodenumber].BPRemoval.push_back({g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]], explanation});
                            
                            if(g_options.debug){
                                g_memoscip.nodeInfoMap[nodenumber].VarLbDebug[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=0;
                                g_memoscip.nodeInfoMap[nodenumber].VarUbDebug[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=0;
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
                                cout<<"Parent node "<<ParentNumber<<" lb "<<g_memoscip.nodeInfoMap[ParentNumber].lb<<" redcost "<<g_memoscip.nodeInfoMap[ParentNumber].redcosts[g_memoscip.eqvid[SCIPvarGetIndex(var)]]<<" new lb "<<SCIPnodeGetLowerbound(node)<<" "<<SCIPgetPrimalbound(scip)<<endl;
                            }
                            g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.push_back({g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]],SCIPgetPrimalbound(scip)-g_memoscip.nodeInfoMap[ParentNumber].lb});
                            assert(SCIPvarGetName(var)[2]!='s');
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

        reason = getBoundChangeReasonDetailed(scip, var, SCIP_BOUNDTYPE_UPPER, oldbound, newbound);
        if(reason!=NULL){
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
                                for(int i = 0; i < nvars; i++) {
                                    SCIP_VAR* v = vars[i];
                                    if(v != var){
                                        //cout<<SCIPvarGetName(v)<<" "<<SCIPvarGetLbLocal(v)<<" "<<SCIPvarGetUbLocal(v)<<endl;
                                        //cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])<<endl;
                                        //cout<<" "<<g_memoscip.VarLb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]<<" "<<g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]<<endl;
                                        assert(g_memoscip.VarLb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==SCIPvarGetLbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])));
                                        assert(g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==SCIPvarGetUbLocal(v) || (SCIPvarGetLbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])==SCIPvarGetUbLocal(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(v)]]]])));
                                        if(coefs[i]>1e-7 && g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==0){
                                            explanation.push_back(g_memoscip.eqvid[SCIPvarGetIndex(v)]);
                                            assert(var!=v);
                                        }
                                        if(coefs[i]<-1e-7 && g_memoscip.VarLb[g_memoscip.eqvid[SCIPvarGetIndex(v)]]==1){
                                            explanation.push_back(g_memoscip.eqvid[SCIPvarGetIndex(v)]);
                                            assert(var!=v);
                                        }
                                    }
                                }
                                //cout<<"Explanation for "<<SCIPvarGetName(var)<<" "<<explanation.size()<<endl;
                                
                                    g_memoscip.nodeInfoMap[nodenumber].BPRemoval.push_back({g_memoscip.eqvid[SCIPvarGetIndex(var)], explanation});
                                
                                if(g_options.debug){
                                    g_memoscip.nodeInfoMap[nodenumber].VarLbDebug[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=0;
                                    g_memoscip.nodeInfoMap[nodenumber].VarUbDebug[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=0;
                                }
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
                                cout<<"Parent node "<<ParentNumber<<" lb "<<g_memoscip.nodeInfoMap[ParentNumber].lb<<" redcost "<<g_memoscip.nodeInfoMap[ParentNumber].redcosts[g_memoscip.eqvid[SCIPvarGetIndex(var)]]<<" new lb "<<SCIPnodeGetLowerbound(node)<<" "<<SCIPgetPrimalbound(scip)<<endl;
                            }
                            assert(SCIPvarGetName(var)[2]!='s');
                            g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.push_back({g_memoscip.eqvid[SCIPvarGetIndex(var)],SCIPgetPrimalbound(scip)-g_memoscip.nodeInfoMap[ParentNumber].lb});
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
        if(SCIPvarGetName(var)[2]=='s'){
            //cout<<"SLack var"<<endl;
            g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=newbound;
            //cout<<g_memoscip.eqvid[SCIPvarGetIndex(var)]<<endl;
        }else{
            //cout<<SCIPinfinity(scip) <<" "<<newbound<<endl;
            if(newbound==SCIPinfinity(scip) || newbound==1){
                g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=1;
                g_memoscip.VarLb[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]]]=0;
            }
            else{
                g_memoscip.VarUb[g_memoscip.eqvid[SCIPvarGetIndex(var)]]=newbound;
                assert(newbound==0);
                g_memoscip.VarLb[g_memoscip.AuxVar[g_memoscip.eqvid[SCIPvarGetIndex(var)]]]=1;
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
    
    //std::cout << "Attaching bound change events to " << nvars << " variables" << std::endl;
    
    for (int i = 0; i < nvars; i++) {
        // ✅ Utiliser SCIPcatchVarEvent (PAS SCIPcatchEvent!)
        SCIP_CALL( SCIPcatchVarEvent(scip, vars[i],
                                     SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
                                     eventhdlr, 
                                     NULL,    // eventdata
                                     NULL) ); // filterpos
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
    
    //SCIP_Longint nodeNumber= SCIPnodeGetNumber(node);

    /*if( g_memoscip.nodeInfoMap[nodeNumber].Status!="INFEASIBLE"){
        cout<<"Node "<<nodeNumber<<" is infeasible"<<endl;
        g_memoscip.nodeInfoMap[nodeNumber].Status="INFEASIBLE";
        SCIP_CALL(gatherLPInfeasibilityCertificate(scip, nodeNumber));
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
    SCIP_Longint nodeNumber=SCIPnodeGetNumber(node);
    if(g_options.verbosity >=1)
        cout<<"Node "<<nodeNumber<<" ub "<<SCIPgetPrimalbound(scip)<<" lb "<<SCIPgetNodeLowerbound(scip, node)<<" GAP "<<SCIPgetGap(scip)<<endl;

    if( *result == SCIP_REDUCEDDOM )
        return SCIP_OKAY;

    if(g_options.branchingHeuristic < 1){
    // Branch on the first fractional variable (in g_memoscip.CurrentVarIdx order)
        for( size_t i = 0; i < g_memoscip.OrderVar.size(); ++i )
        {
            SCIP_VAR* v = g_memoscip.OrderVar[i];
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
    SCIP_VAR** vars = SCIPgetVars(scip);
    int nvars = SCIPgetNVars(scip);
        
    for (int i = 0; i < nvars; i++) {
        // ✅ Utiliser SCIPdropVarEvent (PAS SCIPdropEvent!)
        SCIP_CALL( SCIPdropVarEvent(scip, vars[i],
                                    SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED,
                                    eventhdlr, 
                                    NULL,    // eventdata
                                    -1) );   // filterpos (-1 = auto)
    }
    
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
