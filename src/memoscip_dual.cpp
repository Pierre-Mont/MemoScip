#include "memoscip_dual.h"

#include "memoscip_option.h"
#include "memoscip_types.h"
#include "sepa_mycuts.h"

#include "lpi/lpi.h"
#include "scip/cuts.h"
#include "scip/scip.h"
#include "scip/struct_lp.h"
#include "scip/struct_scip.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace std;


/** Adds weight * row to the proof row. */
static SCIP_RETCODE addRowToProof(SCIP* scip, SCIP_AGGRROW* proofrow, SCIP_ROW* row, SCIP_Real weight)
{
    if( weight < 0.0 )
    {
        assert(-row->lhs<=g_memoscip.TOP);
        SCIP_CALL( SCIPaggrRowAddRow(scip, proofrow, row, weight, -1) );
    }
    else
    {
        assert(row->lhs<=g_memoscip.TOP);
        SCIP_CALL( SCIPaggrRowAddRow(scip, proofrow, row, weight, +1) );
    }
    return SCIP_OKAY;
}

/** Writes the LP and the transformed problem of the node when option -o=N is given with N > 1 and node >= N. */
SCIP_RETCODE writeNodeLPFilesIfRequested(SCIP* scip, SCIP_Longint nodenumber)
{
    if(g_options.output>1 && nodenumber >= g_options.output){
        char nameprobsol[50];
        snprintf(nameprobsol, 50, "solved_%lld.lp", nodenumber);
        SCIP_CALL( SCIPwriteLP(scip, nameprobsol) );
        snprintf(nameprobsol, 50, "problem_%lld.lp", nodenumber);
        SCIP_CALL( SCIPwriteTransProblem(scip, nameprobsol,"lp",FALSE) );
    }
    return SCIP_OKAY;
}

void printProof(SCIP* scip, const vector<SCIP_Real>& coefs, const char* relation, SCIP_Real rhs,
    bool mapThroughCurrentVarIdx)
{
    for (size_t i = 0; i < coefs.size(); i++)
    {
        if(SCIPisDualfeasZero(scip,coefs[i]))
            continue;
        SCIP_VAR* var = mapThroughCurrentVarIdx ? g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]] : g_memoscip.AllVars[i];
        if(coefs[i]>0)
            cout<<"+";
        cout<<coefs[i]<<" "<<SCIPvarGetName(var)<<" ";
    }
    cout<<" "<<relation<<" "<<rhs<<endl;
}


vector<SCIP_Real> FusionResolution(SCIP* scip, SCIP_Longint child0, SCIP_Longint child1)
{
    const NodeInfo& info0 = g_memoscip.nodeInfoMap[child0];
    const NodeInfo& info1 = g_memoscip.nodeInfoMap[child1];
    vector<SCIP_Real> FusCoef(g_memoscip.AllVars.size(),0.0);

    // Coefficients of the branching literal of each child are dropped: the fused proof is valid at the parent.
    assert(info0.DualProofCoef.first.size()==info0.DualProofInds.size());
    for( size_t i = 0; i < info0.DualProofCoef.first.size(); ++i )
    {
        int idx = info0.DualProofInds[i];
        SCIP_Real coef = info0.DualProofCoef.first[i];
        if(!SCIPisDualfeasZero(scip,coef) && idx != g_memoscip.AuxVar[shiftedIndex(g_memoscip.AllVars[info0.BranchVar])])
            FusCoef[idx]=max(coef,0.0);
    }

    assert(info1.DualProofCoef.first.size()==info1.DualProofInds.size());
    for( size_t i = 0; i < info1.DualProofCoef.first.size(); ++i )
    {
        int idx = info1.DualProofInds[i];
        SCIP_Real coef = info1.DualProofCoef.first[i];
        if(!SCIPisDualfeasZero(scip,coef) && idx != g_memoscip.AuxVar[shiftedIndex(g_memoscip.AllVars[info1.BranchVar])])
            FusCoef[idx]=max(FusCoef[idx],coef);
    }

    return FusCoef;
}


void checkSingleChildProof(SCIP_Longint nodenumber, SCIP_Longint child)
{
    auto& nodeInfoMap = g_memoscip.nodeInfoMap;
    assert(nodeInfoMap.contains(child));
    if(g_options.verbosity > 2)
        cout<<"Single child "<<child<<" lb "<<nodeInfoMap[child].lb<<" versus parent lb "<<nodeInfoMap[nodenumber].lb<<" "<<nodeInfoMap[child].DualProofCoef.first.size()<<" "<<nodeInfoMap[child].DualProofCoef.second<<" incr lb "<<nodeInfoMap[child].lb_increase<<endl;
    assert(nodeInfoMap[child].lb<g_memoscip.TOP);
    assert(fabs(nodeInfoMap[child].DualProofCoef.second-nodeInfoMap[child].lb_increase)<1e-5);
}

SCIP_RETCODE fuseChildrenProofs(SCIP* scip, SCIP_Longint nodenumber, const vector<long>& children,
    vector<SCIP_Real>& ResFus, SCIP_Real& rhsFus)
{
    auto& nodeInfoMap = g_memoscip.nodeInfoMap;
    NodeInfo& child0 = nodeInfoMap[children[0]];
    if(!(child0.lb_increase>0 && nodeInfoMap[children[1]].lb_increase>0))
        return SCIP_OKAY;
    NodeInfo& child1 = nodeInfoMap[children[1]];

    assert(nodeInfoMap.contains(children[0]) && nodeInfoMap.contains(children[1]));
    if(g_options.verbosity > 1){
        cout<<"Fusing proofs of children "<<children[0]<<" and "<<children[1]<<" to recover parent node "<<nodenumber<<" proof "<<endl;
        cout<<"proofs of children "<<children[0]<<" "<<child0.DualProofCoef.second<<" lb increase: "<<child0.lb_increase<<endl;
        cout<<"proofs of children "<<children[1]<<" "<<child1.DualProofCoef.second<<" lb increase: "<<child1.lb_increase<<endl;
    }

    assert(fabs(child0.DualProofCoef.second-child0.lb_increase)<1e-5 || child0.LbgUb);
    assert(fabs(child1.DualProofCoef.second-child1.lb_increase)<1e-5 || child1.LbgUb);

    assert(child1.BranchVar==g_memoscip.AuxVar[shiftedIndex(g_memoscip.AllVars[child0.BranchVar])]);
    ResFus=FusionResolution(scip,children[0],children[1]);

    // An infeasible child has no finite bound: the fused bound is the one of the other child.
    const bool feasible0 = child0.Status!="INFEASIBLE";
    const bool feasible1 = child1.Status!="INFEASIBLE";
    if(feasible0 && feasible1){
        rhsFus=min(child0.DualProofCoef.second,child1.DualProofCoef.second);
    }else if(feasible0){
        rhsFus=child0.DualProofCoef.second;
    }else if(feasible1){
        rhsFus=child1.DualProofCoef.second;
    }else{
        rhsFus=g_memoscip.TOP;
        cout<<"This should not happen "<<endl;
        SCIP_CALL( SCIPwriteTransProblem(scip, "model.lp","lp",FALSE) );
        assert(1==0);
    }
    return SCIP_OKAY;
}

void addNodeProof(SCIP_Longint nodenumber, vector<SCIP_Real>& ResFus, SCIP_Real& rhs_summed)
{
    const NodeInfo& node = g_memoscip.nodeInfoMap[nodenumber];
    if(node.lb_increase>0){
        for (size_t i = 0; i < node.DualProofCoef.first.size(); i++)
            ResFus[node.DualProofInds[i]]+=node.DualProofCoef.first[i];
        rhs_summed+=node.DualProofCoef.second;
    }
}

void applyBoundPropagationRemovals(SCIP_Longint nodenumber, vector<SCIP_Real>& ResFus, SCIP_Real rhs_summed)
{
    const auto& BPRemoval = g_memoscip.nodeInfoMap[nodenumber].BPRemoval;
    if(g_options.verbosity > 4)
        cout<<"BOund propagation "<<BPRemoval.size()<<endl;

    for(int i = (int)BPRemoval.size()-1; i >= 0; i--) {
        int idxi = BPRemoval[i].first;
        const vector<int>& explanation = BPRemoval[i].second;
        if(g_options.verbosity > 4)
            cout<<"removed "<<SCIPvarGetName(g_memoscip.AllVars[idxi])<<" "<<explanation.size()<<endl;
        if(ResFus[idxi]>1e-5){
            for(int idxj : explanation) {
                if(g_options.verbosity > 4)
                    cout<<"explanation "<<SCIPvarGetName(g_memoscip.AllVars[idxj])<<endl;
                //Temporary fix because negative coef are hard to handle.
                ResFus[idxj]=rhs_summed*100;
            }
        }
        ResFus[idxi]=0;
    }
}

void applyNodeConsistency(SCIP* scip, SCIP_Longint nodenumber, SCIP_Longint parentNumber,
    vector<SCIP_Real>& ResFus, SCIP_Real& rhs_summed)
{
    const NodeInfo& node = g_memoscip.nodeInfoMap[nodenumber];
    if(g_options.verbosity > 4)
        cout<<"Node consistency "<<node.NodeConsistency.size()<<endl;
    for (const auto& [idx, ubGap] : node.NodeConsistency)
    {
        ResFus[idx]=g_memoscip.nodeInfoMap[parentNumber].redcosts[idx];
        rhs_summed=min(ubGap,rhs_summed);
        if(g_options.verbosity > 4)
            cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[idx]])<< " "<< idx<<" "<<ubGap<<" "<<Myceil(ubGap+node.lb)<<" "<<SCIPgetPrimalbound(scip)<<endl;
    }
}

void checkSummedRhs(SCIP_Longint nodenumber, const vector<long>& children, SCIP_Real rhs_summed)
{
    if(!g_options.debug)
        return;
    auto& nodeInfoMap = g_memoscip.nodeInfoMap;
    if(children.size()==2){
        assert(abs(rhs_summed - (min(nodeInfoMap[children[1]].lb_increase,nodeInfoMap[children[0]].lb_increase) + nodeInfoMap[nodenumber].lb_increase )<1e-9)|| !nodeInfoMap[nodenumber].NodeConsistency.empty());
    }
    if(children.size()==1){
        assert(abs(rhs_summed - (nodeInfoMap[children[0]].lb_increase + nodeInfoMap[nodenumber].lb_increase )<1e-9) || !nodeInfoMap[nodenumber].NodeConsistency.empty());
    }
    (void)nodenumber;
    (void)rhs_summed;
}

void eraseChildren(SCIP_Longint nodenumber, const vector<long>& children)
{
    if(!children.empty()){
        g_memoscip.nodeInfoMap.erase(children[0]);
        if(children.size()==2)
            g_memoscip.nodeInfoMap.erase(children[1]);
        g_memoscip.nodeInfoMap[nodenumber].children.clear();
    }
}


SCIP_RETCODE FinishLearn(SCIP* scip, SCIP_Longint nodenumber){
    NodeInfo& node = g_memoscip.nodeInfoMap[nodenumber];
    vector<long> children = node.children;
    vector<SCIP_Real> ResFus(g_memoscip.CurrentVarIdx.size(),0.0);
    SCIP_Longint parentNumber = node.parent;

    SCIP_Real rhsFus;
    cout<<"Finishing node "<<nodenumber<<endl;

    cout<<"Deleting node "<<nodenumber<<" "<<node.lb<<endl;
    cout<<"number of children "<<node.children.size()<<endl;
    cout<<"Bracnh var "<<node.BranchVar<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[node.BranchVar]])<<endl;
    cout<<"Curnnet node lb increase "<<node.lb_increase<<" versus children "<<endl;

    if(children.empty())
        return SCIP_OKAY;

    if(children.size()==1)
        checkSingleChildProof(nodenumber, children[0]);
    if(children.size()==2){
        FinishLearn(scip,children[0]);
        FinishLearn(scip,children[1]);
        SCIP_CALL( fuseChildrenProofs(scip, nodenumber, children, ResFus, rhsFus) );
    }

    if(g_options.verbosity > 2){
        cout<<"ResFus "<<endl;
        printProof(scip, ResFus, ">=", rhsFus, false);
    }
    SCIP_Real rhs_summed=rhsFus;
    vector<int> NewInds;
    vector<SCIP_Real> NewCoefs;
    vector<SCIP_VAR*> NewVar;

    addNodeProof(nodenumber, ResFus, rhs_summed);
    applyBoundPropagationRemovals(nodenumber, ResFus, rhs_summed);
    applyNodeConsistency(scip, nodenumber, parentNumber, ResFus, rhs_summed);

    cout<<"final res "<<endl;
    printProof(scip, ResFus, ">=", rhs_summed, true);

    checkSummedRhs(nodenumber, children, rhs_summed);

    if(rhs_summed>0 && g_options.learning==1 && rhs_summed < g_memoscip.TOP){
        char name[100];
        snprintf(name, 100, "cut_%lld_%d", nodenumber,SCIPgetDepth(scip));

        for (size_t i = 0; i < ResFus.size(); i++)
        {
            if(!SCIPisDualfeasZero(scip,ResFus[i])){
                NewInds.push_back(i);
                NewCoefs.push_back(ResFus[i]);
                NewVar.push_back(g_memoscip.AllVars[i]);
                assert( i >= g_memoscip.nodeInfoMap[parentNumber].redcosts.size() || ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<1e-3 || g_memoscip.VarUb[i] == g_memoscip.VarLb[i]);
            }
        }

        if(g_options.verbosity > 2)
            cout<<"Adding cut "<<name<<" with rhs "<<rhs_summed<<" and "<<NewVar.size()<< " variables"<<endl;

        assert(NewVar.size()==NewCoefs.size());
        assert(NewVar.size()!=0);
        SCIP_CALL( SCIPqueueCut(scip, name,
                        (int)NewVar.size(),
                        NewVar.data(),
                        NewCoefs.data(),
                        rhs_summed,  // lhs (equality constraint)
                        rhs_summed,&g_memoscip.CutQueue) );
        if(g_options.verbosity > 2)
            SCIPprintQueuedCuts(scip, &g_memoscip.CutQueue);
    }

    assert(abs(node.DualProofCoef.second-node.lb_increase) < 1e-9);
    assert(NewCoefs.size()<=g_memoscip.CurrentVarIdx.size());
    node.DualProofCoef.second=rhs_summed;
    node.DualProofInds=NewInds;
    node.DualProofCoef.first=NewCoefs;
    assert(node.lb_increase<=rhs_summed || !node.NodeConsistency.empty());
    node.lb_increase=rhs_summed;
    assert(ceil(g_memoscip.nodeInfoMap[parentNumber].lb + node.lb_increase+1e-5) >= SCIPgetPrimalbound(scip));
    eraseChildren(nodenumber, children);

    return SCIP_OKAY;
}


//Analyse the dual solution of a general LP to produce a dual proof constraint
SCIP_RETCODE AnalyseDualSol(SCIP* scip, SCIP_Longint nodeNumber){
    ScopedTimer AnalyseDualSolTimer(g_memoscip.AnalyseDualSolTime);
    SolCollector& collector = g_memoscip.collector;
    const size_t nAllVars = g_memoscip.AllVars.size();

    int ncols = SCIPgetNLPCols(scip);
    int nrows = SCIPgetNLPRows(scip);
    SCIP_CALL( writeNodeLPFilesIfRequested(scip, nodeNumber) );

    NodeInfo& node = g_memoscip.nodeInfoMap[nodeNumber];
    NodeInfo& parentInfo = g_memoscip.nodeInfoMap[node.parent];

    vector<SCIP_Real> IndexedDualSol(g_memoscip.nbctr,0.0);
    vector<SCIP_Real> newy_i(nAllVars,0.0);
    vector<SCIP_Real> redcost(nAllVars,0.0);

    node.DualProofInds.clear();
    node.DualProofCoef.first.clear();
    node.DualProofCoef.second=0.0;

    if( collector.dualsol.empty() || nrows != (int)collector.dualsol.size() || ncols != (int)collector.redcost.size() )
        collector.init(scip);  // rebuild mapping + resize buffers

    SCIP_CALL(collector.collect(scip));

    SCIP_COL** cols = SCIPgetLPCols(scip);
    SCIP_ROW** rows = SCIPgetLPRows(scip);
    SCIP_LPSOLSTAT stat = SCIPgetLPSolstat(scip);
    assert(stat==SCIP_LPSOLSTAT_OPTIMAL || stat==SCIP_LPSOLSTAT_OBJLIMIT);
    assert((int)collector.dualsol.size()==nrows);
    assert((int)collector.redcost.size()==ncols);
    (void)stat;
    if(g_options.debug){
        for(int i =0; i < ncols; i++){
            assert(SCIPinDive(scip) || SCIPisFeasEQ(scip, SCIPgetVarRedcost(scip, SCIPcolGetVar(cols[i])), collector.redcost[i]));
            if(!isSlackVar(g_memoscip.AllVars[i])){
                assert(g_options.encoding!=0 || g_memoscip.VarLb[i] == SCIPvarGetLbLocal(g_memoscip.AllVars[i]) || g_memoscip.VarUb[g_memoscip.AuxVar[shiftedIndex(g_memoscip.AllVars[i])]] ==0 );
                assert(g_options.encoding!=0 || g_memoscip.VarUb[i] == SCIPvarGetUbLocal(g_memoscip.AllVars[i]) || g_memoscip.VarLb[g_memoscip.AuxVar[shiftedIndex(g_memoscip.AllVars[i])]] ==1 );
            }
        }
    }

    const vector<SCIP_Real>& parentredcost = (nodeNumber==1) ? g_memoscip.nodeInfoMap[0].redcosts : parentInfo.redcosts;

    //We obtain the dual values corresponding to bound constraints by looking at the redcosts
    //If the variable is fixed to 0 then the dual value is 0.
    //If the variable is fixed to 1 then the dual value is the red cost (bound constraint is equality constaraint and dual value is unsigned)
    //If the variable is not fixed then the dual value is != 0 only if redcost is negative (bound constraint is <= constraint and dual value is negative)
    for(int i=0; i < ncols; i++){
        SCIP_VAR* var = SCIPcolGetVar(cols[i]);
        SCIP_Real redc=collector.redcost[i];
        int idx = allVarsIndex(var);

        redcost[idx]=redc;
        if(!SCIPisDualfeasZero(scip,redc)){
            if(g_memoscip.VarUb[idx]==0)
                newy_i[idx]=0;
            else if(g_memoscip.VarLb[idx]==1)
                newy_i[idx]=redc;
            else
                newy_i[idx]=min(redc,0.0);
            assert(!isSlackVar(g_memoscip.AllVars[idx]) || newy_i[idx]==0);
        }
        assert(g_options.encoding==0 || g_memoscip.TupleAuxVar[shiftedIndex(var)].empty() || g_memoscip.TupleAuxVar[shiftedIndex(var)].size() == 3);
    }

    SCIP_Real sumdual=0.0;
    SCIP_Real Truesumdual=0.0;
    SCIP_AGGRROW* proofrow = nullptr;
    SCIP_AGGRROW* Trueproofrow = nullptr;  // with option -d: proof built from the full dual solution (not the difference with the parent)

    SCIP_CALL(SCIPaggrRowCreate(scip, &proofrow));

    if(g_options.debug){
        SCIP_CALL(SCIPaggrRowCreate(scip, &Trueproofrow));
    }

    //We iterate over all rows and add the rows multiplied by the dual solution (minus the parent one) to the proof row
    for (int i = 0; i < nrows; ++i)
    {
        SCIP_Real dual = collector.dualsol[i];
        int idx = SCIProwGetIndex(rows[i]);
        IndexedDualSol[idx] = dual;
        const bool inParent = idx < (int)parentInfo.DualSol.size();
        if(inParent){
            dual-=parentInfo.DualSol[idx];
        }

        // Rows with a single nonzero are bound constraints, handled through the reduced costs
        if(SCIProwGetNNonz(rows[i])>1){
            if(g_options.verbosity > 2  && dual!=0){
                if(inParent)
                    cout << "Row " << SCIProwGetName(rows[i]) << " Idew "<< SCIProwGetIndex(rows[i]) <<" i "<<i<<" Age: " << SCIProwGetAge(rows[i])<<" "<<" dualsol "<<collector.dualsol[i]<<" parent "<<parentInfo.DualSol[idx]<<" "<<dual<<" sumdual "<<sumdual<<endl;
                else
                    cout << "Row " << SCIProwGetName(rows[i]) <<" Idew "<< SCIProwGetIndex(rows[i])<<" i "<<i<<" dualsol "<<collector.dualsol[i]<<" "<<dual<<" Age: " << SCIProwGetAge(rows[i])<<endl;
            }
            assert(fabs(dual)<1e6);
            if(fabs(dual)  > 1e-7   && SCIProwGetRhs(rows[i]) != SCIPinfinity(scip)){
                sumdual+=dual*SCIProwGetRhs(rows[i]);
                SCIP_CALL( addRowToProof(scip, proofrow, rows[i], dual) );
            }
            if(g_options.debug){
                if( fabs(SCIProwGetRhs(rows[i])*collector.dualsol[i])  > 1e-5  && SCIProwGetRhs(rows[i]) != SCIPinfinity(scip)){
                    Truesumdual+=collector.dualsol[i]*SCIProwGetRhs(rows[i]);
                    SCIP_CALL( SCIPaggrRowAddRow(scip, Trueproofrow, rows[i], collector.dualsol[i], collector.dualsol[i]<0 ? -1 : +1) );
                }
            }
        }else {
            assert(collector.dualsol[i]>=-1e-5);
        }
    }

    if(g_options.debug){
        // Learned cut rows end with their slack variable
        for(int i = 5; i < nrows; ++i)
        {
            SCIP_COL** rowcols = SCIProwGetCols(rows[i]);
            int nnonz = SCIProwGetNNonz(rows[i]);
            if( nnonz > 4 )
            {
                SCIP_VAR* var = SCIPcolGetVar(rowcols[nnonz-1]);
                assert(SCIPvarGetName(var)[0] == 's');
                assert(SCIPvarIsInLP(var) || collector.dualsol[i]==0);
                (void)var;
            }
        }
    }

    vector<SCIP_Real> diffy_i=newy_i;
    assert(g_memoscip.AuxVar.size()==g_memoscip.CurrentVarIdx.size());
    if(nodeNumber>1){
        //If g_memoscip.nodeInfoMap[nodeNumber].y_i.empty() then it is the first time we analyse the node
        //Otherwise it is the second time aka we solved the LP applied NC and BP and resolved the LP
        //If variable is fixed to 0 then dual value is 0.
        for (size_t i = 0; i < parentInfo.y_i.size(); i++){
            if(parentInfo.y_i[i]!=0)
            {
                if(g_memoscip.VarUb[i]==0){
                    assert(newy_i[i]==0);
                    diffy_i[i]=-parentInfo.y_i[i];
                }
                else{
                    diffy_i[i]-=parentInfo.y_i[i];
                    assert(g_memoscip.VarLb[g_memoscip.AuxVar[shiftedIndex(g_memoscip.AllVars[i])]]==0 || SCIPisDualfeasZero(scip,diffy_i[i]) || (g_options.encoding==1 && !g_memoscip.TupleAuxVar.empty()));
                }
            }
            if(g_memoscip.VarLb[i]==1){
                node.Var1.push_back(allVarsIndex(g_memoscip.AllVars[i]));
            }
            
        }
    }
    vector<SCIP_Real> TempCoef(nAllVars,0.0);
    vector<SCIP_Real> TrueTempCoef;
    if(g_options.debug){
        TrueTempCoef=vector<SCIP_Real>(nAllVars,0.0);
    }

    //Create the dual proof coef: aggregated rows plus the change of the bound constraint duals x_i<=1
    for(size_t i = 0; i < nAllVars; ++i )
    {
        SCIP_VAR* var = g_memoscip.AllVars[i];
        TempCoef[i]+=SCIPaggrRowGetProbvarValue(proofrow, SCIPvarGetProbindex(var));
        if(g_options.debug){
            TrueTempCoef[i]+=SCIPaggrRowGetProbvarValue(Trueproofrow, SCIPvarGetProbindex(var));
        }

        if(!SCIPisDualfeasZero(scip,diffy_i[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
            TempCoef[i]+=diffy_i[i];
            forEachAuxVar(shiftedIndex(var), true, [&](int aux){ TempCoef[aux]+=diffy_i[i]; });
            sumdual+=diffy_i[i];
        }
        if(g_options.debug){
            if(!SCIPisDualfeasZero(scip,newy_i[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
                Truesumdual+=newy_i[i];
                TrueTempCoef[i]+=newy_i[i];
                forEachAuxVar(shiftedIndex(var), true, [&](int aux){ TrueTempCoef[aux]+=newy_i[i]; });
            }
        }
    }

    assert(nAllVars == TempCoef.size());
    if(g_options.verbosity > 2){
        cout<<"DUAL proof"<<endl;
        printProof(scip, TempCoef, ">=", sumdual, true);
    }

    if(g_options.verbosity > 2 && g_options.debug){
        cout<<"TRUEDUAL proof"<<endl;
        printProof(scip, TrueTempCoef, ">=", Truesumdual, true);
    }
    assert(!parentredcost.empty());
    for(size_t i = 0; i < nAllVars; ++i )
    {
        SCIP_Real coef = TempCoef[i];
        //compute the true redcost
        if(i<parentredcost.size())
            redcost[i]=parentredcost[i]-coef;
        else
            redcost[i]=-coef;
        if(g_memoscip.VarLb[i]==g_memoscip.VarUb[i]){
            redcost[i]=0;
        }
        if(SCIPisDualfeasZero(scip,coef))
            continue;
        node.DualProofCoef.first.push_back(coef);
        node.DualProofInds.push_back(i);

        if(g_options.verbosity > 3){
            if(i<parentredcost.size()){
                cout<<"COEF IN PROOF "<<SCIPvarGetName(g_memoscip.AllVars[i])<<" "<<coef<<" vs parent red "<<parentredcost[i]<<" new red "<<parentredcost[i]-coef<< " LB "<<g_memoscip.VarLb[i]<<" UB "<<g_memoscip.VarUb[i]<<endl;
                cout<<parentredcost[i]-coef<<endl;
                assert(parentredcost[i]-coef>-9.9e-1|| 0==g_memoscip.VarUb[i]);
            }else {
                cout<<SCIPvarGetName(g_memoscip.AllVars[i])<<" "<<coef<<endl;
                assert(isSlackVar(g_memoscip.AllVars[i]));
                assert(coef<1e-9);
            }
        }
    }
    node.DualSol=IndexedDualSol;
    node.y_i=newy_i;
    node.redcosts = redcost;
    assert(!redcost.empty());
    node.lb=collector.objval;
    node.DualProofCoef.second=sumdual;
    node.Status="FEASIBLE";
    node.lb_increase=sumdual;
    if(g_options.verbosity > 1){
        cout<<collector.objval<<" "<<parentInfo.lb<<" Truesumdual "<<Truesumdual<<" "<<Truesumdual-collector.objval<<endl;
        cout<<"DIFF Objval "<<collector.objval-parentInfo.lb<<" versus sumdiffdual "<<sumdual<<endl;
    }
    assert(fabs(collector.objval-sumdual-parentInfo.lb) < 1e-1|| collector.objval>= SCIPgetPrimalbound(scip));

    SCIPaggrRowFree(scip, &proofrow);
    if(Trueproofrow != nullptr)
        SCIPaggrRowFree(scip, &Trueproofrow);

    // With option -d: check that the reduced costs are correct, i.e. that the LP with objective redcost has value 0
    if(g_options.debug && !SCIPinDive(scip))
    {
        SCIP_CALL( SCIPstartProbing(scip) );
        SCIP_Bool lperror, cutoff;
        for (size_t i = 0; i < nAllVars; ++i)
        {
            SCIP_CALL( SCIPchgVarObjProbing(scip, g_memoscip.AllVars[i], redcost[i]) );
        }
        SCIP_CALL( SCIPflushLP(scip) );
        SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );

        SCIP_CALL( SCIPwriteTransProblem(scip, "probling.lp","lp",FALSE) );
        assert(SCIPgetLPSolstat(scip)==SCIP_LPSOLSTAT_OPTIMAL);

        assert(SCIPgetLPObjval(scip)<1e-2);
        SCIP_CALL( SCIPendProbing(scip));
    }

    return SCIP_OKAY;
}

/** Marks the node as infeasible: its lower bound and proof right-hand side become TOP. */
static void markNodeInfeasible(NodeInfo& node)
{
    node.lb=g_memoscip.TOP;
    node.DualProofCoef.second=g_memoscip.TOP;
    node.Status="INFEASIBLE";
    node.lb_increase=g_memoscip.TOP;
}

//Get an infeasibility certificate for the current LP relaxation of the node
SCIP_RETCODE gatherLPInfeasibilityCertificate(
    SCIP* scip,
    SCIP_Longint nodeNumber)
{
    ScopedTimer GatherLPInfeasibilityCertificateTimer(g_memoscip.GatherLPInfeasibilityCertificateTime);
    NodeInfo& node = g_memoscip.nodeInfoMap[nodeNumber];
    const size_t nAllVars = g_memoscip.AllVars.size();
    int nrows = SCIPgetNLPRows(scip);
    SCIP_ROW** rows = SCIPgetLPRows(scip);

    //If INF is true then we collect the infeasibility certificate from farkas rows
    //Otherwise return a basic infeasibility certificate from the variables that are fixed to 0
    if(g_options.INF){
        SCIP_LPI* lpi = scip->lp->lpi;
        SCIP_AGGRROW* proofrow = nullptr;
        SCIP_CALL(SCIPaggrRowCreate(scip, &proofrow));
        int ncols = SCIPgetNLPCols(scip);
        vector<SCIP_Real> rowFarkas(nrows+ncols);
        SCIP_CALL( SCIPlpiGetDualfarkas(lpi, rowFarkas.data()));
        for(int i = 0; i < nrows+ncols; ++i)
        {
            if(!SCIPisDualfeasZero(scip,rowFarkas[i])){
                if(g_options.verbosity > 2){
                    if(i<nrows)
                        cout<<"Row "<<i<<" "<<SCIProwGetName(rows[i])<<" "<<rowFarkas[i]<<endl;
                    else
                        cout<<"Col "<<i-nrows<<" "<<SCIPvarGetName(g_memoscip.AllVars[i-nrows])<<" "<<rowFarkas[i]<<endl;
                }
                SCIP_Bool negated = rowFarkas[i] < 0.0;

                if(SCIProwGetRhs(rows[i])==SCIProwGetLhs(rows[i]))
                    SCIP_CALL( SCIPaggrRowAddRow(scip, proofrow, rows[i], rowFarkas[i], negated ? -1 : 1) );
                else
                    assert(SCIProwGetName(rows[i])[0]=='s' && SCIProwGetName(rows[i])[1]=='c');
            }
        }
        vector<SCIP_Real> Coefs(nAllVars,0.0);
        bool negrhs = SCIPaggrRowGetRhs(proofrow)<0;
        [[maybe_unused]] SCIP_Real testup=0;  // maximum activity of the proof under the current bounds
        for (size_t i = 0; i < nAllVars; ++i)
        {
            SCIP_VAR* var = g_memoscip.AllVars[i];
            SCIP_Real coef = SCIPaggrRowGetProbvarValue(proofrow, SCIPvarGetProbindex(var));

            if(g_memoscip.VarLb[i]==1 || (coef >0 && g_memoscip.VarUb[i]!=0))
                testup+=coef;

            if(negrhs)
                coef = -coef;

            if(!SCIPisDualfeasZero(scip,coef)){
                //Only keep variables assigned to 0 in the proof.
                //A variable assigned to 1 is replaced by its auxiliary variable(s), which are assigned to 0.
                if(g_memoscip.VarUb[i]==0 && coef > 0){
                    Coefs[i]+=coef;
                }else if(g_memoscip.VarLb[i]==1 && coef < 0){
                    forEachAuxVar(shiftedIndex(var), true, [&](int aux){ Coefs[aux]-=coef; });
                }
            }
        }
        assert(testup<SCIPaggrRowGetRhs(proofrow));
        if(g_options.debug){
            // Learned cut rows end with their slack variable
            for(int i = 5; i < nrows; ++i)
            {
                SCIP_COL** rowcols = SCIProwGetCols(rows[i]);
                int nnonz = SCIProwGetNNonz(rows[i]);
                if( nnonz > 5)
                    assert(SCIPvarGetName(SCIPcolGetVar(rowcols[nnonz-1]))[0] == 's');
                (void)rowcols;
            }
        }

        for(size_t i = 0; i < nAllVars; ++i )
        {
            if(SCIPisDualfeasZero(scip,Coefs[i]))
                continue;
            assert(g_memoscip.VarUb[i]==0);
            node.DualProofCoef.first.push_back(g_memoscip.TOP);
            node.DualProofInds.push_back(i);
        }
        markNodeInfeasible(node);
        SCIPaggrRowFree(scip, &proofrow);

    }else{
        for (size_t i = 0; i < nAllVars; ++i)
        {
            if(g_memoscip.VarUb[i]==0){
                node.DualProofCoef.first.push_back(g_memoscip.TOP);
                node.DualProofInds.push_back(i);
            }
        }
        markNodeInfeasible(node);
        if(g_options.debug){
            // Learned cut rows end with their slack variable
            for(int i = 5; i < nrows; ++i)
            {
                SCIP_ROW* row = rows[i];
                SCIP_COL** rowcols = SCIProwGetCols(row);
                int nnonz = SCIProwGetNNonz(row);

                if( nnonz > 5 )
                {
                    SCIP_VAR* var = SCIPcolGetVar(rowcols[nnonz-1]);
                    if(g_options.verbosity > 7){
                        cout << "Row " << SCIProwGetName(row) << " Idew "<< SCIProwGetIndex(row) <<" last var "<<SCIPvarGetName(var) << " status "<<SCIPvarGetStatus(var)<<endl;
                        cout<<"pose "<<SCIPcolGetLPPos(SCIPvarGetCol(var))<<endl;

                        printf("status       = %d\n", SCIPvarGetStatus(var));
                        printf("inLP         = %d\n", SCIPvarIsInLP(var));
                        printf("column       = %p\n", (void*)SCIPvarGetCol(var));
                        printf("lp position  = %d\n", SCIPcolGetLPPos(SCIPvarGetCol(var)));
                        printf("nLPCols      = %d\n", SCIPgetNLPCols(scip));
                    }
                    assert(SCIPvarGetName(var)[0] == 's');
                }
            }
        }
    }

    return SCIP_OKAY;
}



//Generate a dual proof when a node was cut off due to pseudo objective higher than UB
SCIP_RETCODE getDualProofFromPseudoObjCutoff(
    SCIP* scip,
    SCIP_Longint nodeNumber, SCIP_Longint parentNumber, int decVar)
{
    ScopedTimer GetDualProofFromPseudoObjCutoffTimer(g_memoscip.GetDualProofFromPseudoObjCutoffTime);
    if(g_options.verbosity > 3){
        cout << "Generating dual proof for cut off node " << nodeNumber << " with parent " << parentNumber << endl;
    }
    const size_t nAllVars = g_memoscip.AllVars.size();
    NodeInfo& parentInfo = g_memoscip.nodeInfoMap[parentNumber];
    vector<SCIP_Real> y_i(nAllVars,0.0);
    vector<SCIP_Real> diffy(nAllVars,0.0);

    // The pseudo objective of the node is the objective of the variables fixed to 1 (parent ones plus decVar)
    SCIP_Real obj=0;
    for (int idx : parentInfo.Var1){
        y_i[idx] = SCIPvarGetObj(g_memoscip.AllVars[idx]);
        obj += y_i[idx];
    }
    y_i[decVar] = SCIPvarGetObj(g_memoscip.AllVars[decVar]);
    obj += y_i[decVar];
    for(size_t i = 0; i < nAllVars; ++i){
        if(i<parentInfo.y_i.size()){
            diffy[i]=y_i[i]-parentInfo.y_i[i];
        }
    }
    assert(obj >= SCIPgetPrimalbound(scip));

    SCIP_Real sumdual=0;
    SCIP_AGGRROW* proofrow = nullptr;

    int nrows = SCIPgetNLPRows(scip);
    SCIP_ROW** rows = SCIPgetLPRows(scip);

    SCIP_CALL(SCIPaggrRowCreate(scip, &proofrow));

    //We iterate over all rows and remove the rows multiplied by the parent dual solution from the proof row
    for (int i = 0; i < nrows; ++i)
    {
        int idx = SCIProwGetIndex(rows[i]);
        if(idx >= (int)parentInfo.DualSol.size())
            continue;
        SCIP_Real dual=-parentInfo.DualSol[idx];

        // Rows with a single nonzero are bound constraints, handled through y_i
        if(SCIProwGetNNonz(rows[i])>1){
            if(g_options.verbosity > 2  && dual!=0){
                cout << "Row " << SCIProwGetName(rows[i]) << " Idew "<< SCIProwGetIndex(rows[i]) <<" i "<<i<<" Age: " << SCIProwGetAge(rows[i])<<" "<<" dualsol "<<g_memoscip.collector.dualsol[i]<<" parent "<<parentInfo.DualSol[idx]<<" "<<dual<<" sumdual "<<sumdual<<endl;
            }
            assert(fabs(dual)<1e6);
            if(fabs(SCIProwGetRhs(rows[i])*dual)  > 1e-7   && SCIProwGetRhs(rows[i]) != SCIPinfinity(scip)){
                sumdual+=dual*SCIProwGetRhs(rows[i]);
                SCIP_CALL( addRowToProof(scip, proofrow, rows[i], dual) );
            }
        }else {
            assert(g_memoscip.collector.dualsol[i]>=-1e-5);
        }
    }

    vector<SCIP_Real> TempCoef(nAllVars,0.0);
    //Create the dual proof coef
    for(size_t i = 0; i < nAllVars; ++i )
    {
        SCIP_VAR* var = g_memoscip.AllVars[i];
        TempCoef[i]+=SCIPaggrRowGetProbvarValue(proofrow, SCIPvarGetProbindex(var));

        if(!SCIPisDualfeasZero(scip,diffy[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
            TempCoef[i]+=diffy[i];
            forEachAuxVar(shiftedIndex(var), true, [&](int aux){ TempCoef[aux]+=diffy[i]; });
            sumdual+=diffy[i];
        }
    }
    SCIPaggrRowFree(scip, &proofrow);
    assert(SCIPisDualfeasZero(scip,obj-sumdual-parentInfo.lb));

    NodeInfo& node = g_memoscip.nodeInfoMap[nodeNumber];
    for(size_t i = 0; i < nAllVars; ++i )
    {
        SCIP_Real coef = TempCoef[i];
        if(SCIPisDualfeasZero(scip,coef))
            continue;
        node.DualProofCoef.first.push_back(coef);
        node.DualProofInds.push_back(i);
        if(g_options.verbosity > 3){
            cout<<"COEF IN PROOF "<<SCIPvarGetName(g_memoscip.AllVars[i])<<" "<<coef<<" LB "<<g_memoscip.VarLb[i]<<" UB "<<g_memoscip.VarUb[i]<<endl;
        }
    }
    node.lb=obj;
    node.DualProofCoef.second=sumdual;
    node.Status="FEASIBLE";
    node.lb_increase=sumdual;

    return SCIP_OKAY;
}
