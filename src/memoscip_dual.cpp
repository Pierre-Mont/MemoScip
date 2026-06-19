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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace std;


vector<SCIP_Real> FusionResolution(SCIP* scip, SCIP_Longint child0, SCIP_Longint child1)
{
    vector<SCIP_Real> CDP1 = g_memoscip.nodeInfoMap[child0].DualProofCoef.first;
    vector<SCIP_Real> CDP2 = g_memoscip.nodeInfoMap[child1].DualProofCoef.first;
    size_t nnz1 = CDP1.size();
    vector<int> inds1 = g_memoscip.nodeInfoMap[child0].DualProofInds;
    size_t i,idx;
    int nvars = g_memoscip.CurrentVarIdx.size();
    vector<SCIP_Real> FusCoef(nvars,0.0);
    
    //printf("AGGRROW: rhs=%g, nnz=%d\n", SCIPaggrRowGetRhs(DP1), nnz1);
    /*cout<<SCIPaggrRowGetRhs(DP1)<<" "<<g_memoscip.nodeInfoMap[child0].DualProofCoef.second<<endl;
    cout<<SCIPaggrRowGetRhs(DP2)<<" "<<g_memoscip.nodeInfoMap[child1].DualProofCoef.second<<endl;
    assert(g_memoscip.nodeInfoMap[child0].DualProofCoef.second==SCIPaggrRowGetRhs(DP1));
    assert(g_memoscip.nodeInfoMap[child1].DualProofCoef.second==SCIPaggrRowGetRhs(DP2));*/
    //cout<<"START CDP1 "<<endl;
    /*cout<<"FUSION RESOLUTION between node "<<child0<<" and "<<child1<<" with parent "<<g_memoscip.nodeInfoMap[child0].parent<<endl;
    cout<< "Child0 proof rhs: " << g_memoscip.nodeInfoMap[child0].DualProofCoef.second << endl;
    cout<< "Child1 proof rhs: " << g_memoscip.nodeInfoMap[child1].DualProofCoef.second << endl;
    cout<<CDP1.size()<<" "<<g_memoscip.nodeInfoMap[child0].DualProofVars.size()<<endl;*/
    assert(nnz1==g_memoscip.nodeInfoMap[child0].DualProofVars.size());
    assert(nnz1==inds1.size());
    //cout<<"child 0 with "<<nnz1<<" variables and rhs "<< g_memoscip.nodeInfoMap[child0].DualProofCoef.second<<" branch var "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[child0].BranchVar]])<<endl;
    for( i = 0; i < nnz1; ++i )
    {
        idx = inds1[i];
        SCIP_Real coef = CDP1[i];
        if(!SCIPisDualfeasZero(scip,coef) && g_memoscip.nodeInfoMap[child0].DualProofInds[i] != g_memoscip.AuxVar[g_memoscip.nodeInfoMap[child0].BranchVar])
            FusCoef[idx]=max(coef,0.0);
        //cout<<i<<" "<<coef<<" * "<<SCIPvarGetName(g_memoscip.nodeInfoMap[child0].DualProofVars[i])<<" "<<FusCoef[idx]<<endl;
    }

    
    size_t nnz2 = CDP2.size();
    vector<int> inds2 = g_memoscip.nodeInfoMap[child1].DualProofInds;
    //printf("AGGRROW: rhs=%g, nnz=%d\n", SCIPaggrRowGetRhs(DP2), nnz2);
    //cout<<"START CDP2 "<<endl;
    //cout<<"child 1 with "<<nnz2<<" variables and rhs "<<g_memoscip.nodeInfoMap[child1].DualProofCoef.second<<" branch var "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[child1].BranchVar]])<<endl;
    assert(nnz2==g_memoscip.nodeInfoMap[child1].DualProofVars.size());
    assert(nnz2==inds2.size());
    for( i = 0; i < nnz2; ++i )
    {
        idx = inds2[i];
        SCIP_Real coef = CDP2[i];
        if(!SCIPisDualfeasZero(scip,coef) && g_memoscip.nodeInfoMap[child1].DualProofInds[i] != g_memoscip.AuxVar[g_memoscip.nodeInfoMap[child1].BranchVar]){
            FusCoef[idx]=max(FusCoef[idx],coef);
        }
        //cout<<i<<" "<<coef<<" * "<<SCIPvarGetName(g_memoscip.nodeInfoMap[child1].DualProofVars[i])<<" "<<FusCoef[idx]<<endl;
    }
    
    return FusCoef;
}
SCIP_RETCODE FinishLearn(SCIP* scip, SCIP_Longint node){
    vector<long> children = g_memoscip.nodeInfoMap[node].children;
    vector<SCIP_Real> ResFus(g_memoscip.CurrentVarIdx.size(),0.0);
    SCIP_Real rhsFus;
    cout<<"Finishing node "<<node<<endl;
    if(children.empty())
        return SCIP_OKAY;
    if(children.size()==1){
        FinishLearn(scip,children[0]);
        assert(g_memoscip.nodeInfoMap[children[0]].lb<g_memoscip.TOP);
        vector<SCIP_Real> CDP1 = g_memoscip.nodeInfoMap[children[0]].DualProofCoef.first;
        size_t nnz1 = CDP1.size();
        vector<int> inds1 = g_memoscip.nodeInfoMap[children[0]].DualProofInds;
        assert(nnz1==g_memoscip.nodeInfoMap[children[0]].DualProofVars.size());
        assert(nnz1==inds1.size());
        for(size_t i = 0; i < nnz1; ++i )
        {
            SCIP_Real coef = CDP1[i];
            if(!SCIPisDualfeasZero(scip,coef) && g_memoscip.nodeInfoMap[children[0]].DualProofInds[i] != g_memoscip.AuxVar[g_memoscip.nodeInfoMap[children[0]].BranchVar])
                ResFus[inds1[i]]=max(coef,0.0);
            //cout<<i<<" "<<coef<<" * "<<SCIPvarGetName(g_memoscip.nodeInfoMap[child0].DualProofInds[i])<<endl;
        }
        rhsFus=g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second;
    }
    if(children.size()==2){
        FinishLearn(scip,children[0]);
        FinishLearn(scip,children[1]);

        assert(g_memoscip.nodeInfoMap.contains(children[0]) && g_memoscip.nodeInfoMap.contains(children[1]));
    
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
    for (size_t i = 0; i < g_memoscip.nodeInfoMap[node].DualProofCoef.first.size(); i++)
    {
        int idx =g_memoscip.nodeInfoMap[node].DualProofInds[i];
        //cout<<i<<" "<<SCIPvarGetName(g_memoscip.nodeInfoMap[node].DualProofVars[i])<<" ResFus "<<ResFus[idx]<<" dual "<<g_memoscip.nodeInfoMap[node].DualProofCoef.first[i]<<endl;
        ResFus[idx]+=g_memoscip.nodeInfoMap[node].DualProofCoef.first[i];
        //cout<<"red cost "<<g_memoscip.nodeInfoMap[node].redcosts[i]<<" FIncal ResFus "<<ResFus[idx]<<endl;
        assert(ResFus[idx]<g_memoscip.TOP);
    }
    SCIP_Real rhs_summed=rhsFus+g_memoscip.nodeInfoMap[node].DualProofCoef.second;
    SCIP_Longint parentNumber=g_memoscip.nodeInfoMap[node].parent;
    vector<int> NewInds;
    vector<SCIP_Real> NewCoefs;
    vector<SCIP_VAR*> NewVar;
    //cout<<"New cut with "<<ResFus.size()<<" variables "<< g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first.size()<<endl;
    for (size_t i = 0; i < ResFus.size(); i++)
    {
        
        if(!SCIPisDualfeasZero(scip,ResFus[i])){
            NewInds.push_back(i);
            NewCoefs.push_back(ResFus[i]);
            NewVar.push_back(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]]);
            
            assert(ResFus[i]<g_memoscip.TOP);
        }
        if(i < g_memoscip.nodeInfoMap[parentNumber].redcosts.size()){
            //cout<<i<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" "<<ResFus[i]<<" "<<g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<endl;
                //cout<<g_memoscip.VarLb[i]<<" "<<g_memoscip.VarUb[i]<<" "<<g_memoscip.VarLb[g_memoscip.AuxVar[i]]<<" "<<g_memoscip.VarUb[g_memoscip.AuxVar[i]]<<endl;
            //assert(g_memoscip.nodeInfoMap[parentNumber].redcosts[i] - ResFus[i] > -1e-9 ||( i> g_memoscip.AuxVar.size() || (g_memoscip.VarUb[i]==0 || g_memoscip.VarLb[g_memoscip.AuxVar[i]]==1)));
        }
        else
            assert(SCIPisDualfeasZero(scip,ResFus[i]));
    }



    SCIP_VAR* slack=g_memoscip.AllVars[g_memoscip.slackVar[g_memoscip.nbAddedctr]]; 
   
    NewVar.push_back(slack);
    NewInds.push_back(g_memoscip.CurrentVarIdx.size());
    NewCoefs.push_back(-1.0);
    g_memoscip.VarUb.push_back(SCIPinfinity(scip));
    g_memoscip.VarLb.push_back(0);
    g_memoscip.AuxVar.push_back(g_memoscip.CurrentVarIdx.size());
    g_memoscip.CurrentVarIdx.push_back(g_memoscip.slackVar[g_memoscip.nbAddedctr]);
    //cout<<"MY lb incr "<<g_memoscip.nodeInfoMap[node].lb_increase<<endl;

    if(children.size()==2){
        //cout<<"Adding cut with rhs "<<rhs_summed<<" vs "<<(min(g_memoscip.nodeInfoMap[children[1]].lb_increase,g_memoscip.nodeInfoMap[children[0]].lb_increase) + g_memoscip.nodeInfoMap[node].lb_increase )<<endl;
        //cout<<"Min children "<<min(g_memoscip.nodeInfoMap[children[1]].lb_increase,g_memoscip.nodeInfoMap[children[0]].lb_increase) <<endl;
        assert(abs(rhs_summed - (min(g_memoscip.nodeInfoMap[children[1]].lb_increase,g_memoscip.nodeInfoMap[children[0]].lb_increase) + g_memoscip.nodeInfoMap[node].lb_increase )<1e-9)); 

    }else{
        assert(abs(rhs_summed - (g_memoscip.nodeInfoMap[children[0]].lb_increase + g_memoscip.nodeInfoMap[node].lb_increase )<1e-9)); 
    }

    assert(abs(g_memoscip.nodeInfoMap[node].DualProofCoef.second-g_memoscip.nodeInfoMap[node].lb_increase) < 1e-9);

    g_memoscip.nodeInfoMap[node].DualProofCoef.second=rhs_summed;
    g_memoscip.nodeInfoMap[node].DualProofInds=NewInds;
    g_memoscip.nodeInfoMap[node].DualProofCoef.first=NewCoefs;
    g_memoscip.nodeInfoMap[node].DualProofVars=NewVar;
    g_memoscip.nodeInfoMap[node].lb_increase=rhs_summed;
    //cout<<"rhs summed "<<rhs_summed<<" "<<g_memoscip.nodeInfoMap[parentNumber].lb<<endl;
    //cout<<ceil(rhs_summed+g_memoscip.nodeInfoMap[parentNumber].lb)<<" "<<SCIPgetPrimalbound(scip)<<endl;
    assert(SCIPgetPrimalbound(scip) == SCIPinfinity(scip)  || ceil(rhs_summed+g_memoscip.nodeInfoMap[parentNumber].lb) - SCIPgetPrimalbound(scip) > -1e-9);

    char name[100];
    snprintf(name, 100, "cut_%lld", node);
    cout<<"Adding cut "<<name<<" with rhs "<<rhs_summed<<" and "<<NewVar.size()<<" variables"<<endl;
    
    SCIP_CALL( SCIPqueueCut(scip, name,
                    (int)NewVar.size(),
                    NewVar.data(),
                    NewCoefs.data(),
                    rhs_summed,  // lhs (pour inégalité)
                    rhs_summed,&g_memoscip.CutQueue) ); 
    
    for (size_t i = 0; i <  g_memoscip.nodeInfoMap[node].DualProofCoef.first.size(); i++)
    {
        if(!SCIPisDualfeasZero(scip, g_memoscip.nodeInfoMap[node].DualProofCoef.first[i])){
            if( g_memoscip.nodeInfoMap[node].DualProofCoef.first[i]<0)
                cout<<g_memoscip.nodeInfoMap[node].DualProofCoef.first[i]<<" "<<SCIPvarGetName(NewVar[i])<<" ";
            else
                cout<<"+"<<g_memoscip.nodeInfoMap[node].DualProofCoef.first[i]<<" "<<SCIPvarGetName(NewVar[i])<<" ";
        }
        
    }
    cout<<" >= "<<rhs_summed<<endl;

    g_memoscip.nbAddedctr++;

    return SCIP_OKAY;

}



SCIP_RETCODE AnalyseDualSol(SCIP* scip, SCIP_Longint nodeNumber){
    
    SCIP_LPI* lpi;
    SCIP_CALL(SCIPgetLPI(scip, &lpi));

    int ncols = SCIPgetNLPCols(scip);
    int nrows = SCIPgetNLPRows(scip);

    
    vector<SCIP_Real> parentredcost;
    
    vector<SCIP_Real> diffdualsol(nrows);
    
    vector<SCIP_Real> diffredcost(ncols);
    vector<SCIP_Real> diffy_i(g_memoscip.CurrentVarIdx.size());
    vector<SCIP_Real> newy_i(g_memoscip.CurrentVarIdx.size(),0.0);

    vector<SCIP_Real> redcost(g_memoscip.CurrentVarIdx.size(),0.0);
    if( g_memoscip.collector.primsol.empty() )
        g_memoscip.collector.init(scip);

    if( SCIPgetNLPCols(scip) != (int)g_memoscip.collector.primsol.size() || SCIPgetNLPRows(scip) != (int)g_memoscip.collector.dualsol.size() )
    {
        g_memoscip.collector.init(scip);  // rebuild mapping + resize buffers
    }

    SCIP_CALL(g_memoscip.collector.collect(scip));
    
    SCIP_COL** cols = SCIPgetLPCols(scip);
    SCIP_ROW** rows = SCIPgetLPRows(scip);

    if(g_options.debug){
        for( size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); i++ )
        {
            int currv=g_memoscip.CurrentVarIdx[i];
            SCIP_VAR* var = SCIPcolGetVar(cols[currv]);
            SCIP_Real rc  = g_memoscip.collector.redcost[currv];  
            redcost[i]=rc;

            //cout<<SCIPvarGetName(g_memoscip.AllVars[currv])<<" "<<SCIPvarGetName(var)<<endl;
            //cout<<SCIPgetVarRedcost(scip,g_memoscip.AllVars[currv])<<" "<<g_memoscip.collector.redcost[currv]<<endl;
            assert(var == g_memoscip.AllVars[currv]);

            assert(SCIPisFeasEQ(scip, SCIPgetVarRedcost(scip, g_memoscip.AllVars[currv]), g_memoscip.collector.redcost[currv]));
        }
    }




    if (nodeNumber==1)
    {
        parentredcost=g_memoscip.nodeInfoMap[0].redcosts;
    }else{
        if(g_memoscip.nodeInfoMap[nodeNumber].DualSol.empty())
            parentredcost=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].redcosts;
        else
            parentredcost=g_memoscip.nodeInfoMap[nodeNumber].redcosts;
    }


    //We obtain the dual values corresponding to bound constraints by looking at the redcosts
    //If the variable is fixed to 0 then the dual value is 0. 
    //If the variable is fixed to 1 then the dual value is the red cost (bound constraint is equality constaraint and dual value is unsigned)
    //If the variable is not fixed then the dual value is != 0 only if redcost is negative (bound constraint is <= constraint and dual value is negative)
    SCIP_Real redc;
    for (size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); i++)
    {
        redc=g_memoscip.collector.redcost[g_memoscip.CurrentVarIdx[i]];
        redcost[g_memoscip.CurrentVarIdx[i]]=redc;
        if(!SCIPisDualfeasZero(scip,redc)){
            //cout<<i<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" "<<g_memoscip.VarLb[i]<<" "<<g_memoscip.VarUb[i]<<" "<<redc<<endl;
            if(g_memoscip.VarUb[i]==0)
                newy_i[i]=0;
            else if(g_memoscip.VarLb[i]==1){
                newy_i[i]=redc;
            }
            else{
                newy_i[i]=min(redc,0.0);
            }
            assert(SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])[2]!='s' || newy_i[i]==0);
        }
    }
    
    SCIP_Real sumdual=0.0;
    SCIP_Real Truesumdual=0.0;
    SCIP_SET* set= scip->set;  
    SCIP_AGGRROW* proofrow = nullptr;
    SCIP_AGGRROW* Trueproofrow = nullptr;

    SCIP_CALL(SCIPaggrRowCreate(scip, &proofrow));

    if(g_options.debug){
        SCIP_CALL(SCIPaggrRowCreate(scip, &Trueproofrow));
    }

    diffdualsol=g_memoscip.collector.dualsol;
    if(nodeNumber>1){
        //If g_memoscip.nodeInfoMap[nodeNumber].DualSol.empty() then it is the first time we analyse the node 
        //Otherwise it is the second time aka we solved the LP applied NC and BP and resolved the LP 
        if(g_memoscip.nodeInfoMap[nodeNumber].DualSol.empty()){
            //cout<<"First time analysing"<<endl;

            //We simply do the difference between the parent dual solution and the current dual solution
            assert(diffdualsol.size()>= g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol.size());
            for (size_t i = 0; i < g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol.size(); i++){
                diffdualsol[i]-=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol[i];
            }
        }
        else{
            //cout<<"Second time analysing"<<endl;
            //We save the first dual proof coef and inds and vars
            //We do the difference between the current dual solution and the first dual solution
            g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofCoef=g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef;
            g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofInds=g_memoscip.nodeInfoMap[nodeNumber].DualProofInds;
            g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofVars=g_memoscip.nodeInfoMap[nodeNumber].DualProofVars;
            g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.clear();
            g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.clear();
            g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=0;
            g_memoscip.nodeInfoMap[nodeNumber].DualProofVars.clear();
            for (size_t i = 0; i < g_memoscip.nodeInfoMap[nodeNumber].DualSol.size(); i++){
                diffdualsol[i]-=g_memoscip.nodeInfoMap[nodeNumber].DualSol[i];
            }
        }
    }

    //We iterate over all rows and add the rows multiplied by the dual solution to the proof row
    for (int i = 0; i < nrows; ++i)
    {
        SCIP_Real dual = diffdualsol[i];
        //if size of a row is greater than 1
        if(SCIProwGetNNonz(rows[i])>1){
            if(g_options.verbosity > 1 && dual!=0){
                if(i<(int)g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol.size())
                    cout << "Row " << SCIProwGetName(rows[i]) <<" "<<g_memoscip.collector.dualsol[i]<<" "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol[i]<<" "<<dual<<endl;
                else 
                    cout << "Row " << SCIProwGetName(rows[i]) <<" "<<g_memoscip.collector.dualsol[i]<<" "<<dual<<endl;
            }
            if(fabs(SCIProwGetRhs(rows[i])*dual)  > 1e-5   && SCIProwGetRhs(rows[i]) != SCIPinfinity(scip)){
                sumdual+=dual*SCIProwGetRhs(rows[i]);
                if( dual < 0.0 )
                {
                    assert(-rows[i]->lhs<=g_memoscip.TOP);
                    SCIP_CALL( SCIPaggrRowAddRow(set->scip, proofrow, rows[i], dual, -1) );
                }
                else
                {
                    assert(rows[i]->lhs<=g_memoscip.TOP);
                    SCIP_CALL( SCIPaggrRowAddRow(set->scip, proofrow, rows[i], dual, +1) );
                }
            }
            if(g_options.debug){
                if( fabs(SCIProwGetRhs(rows[i])*g_memoscip.collector.dualsol[i])  > 1e-5  && SCIProwGetRhs(rows[i]) != SCIPinfinity(scip)){
                    Truesumdual+=g_memoscip.collector.dualsol[i]*SCIProwGetRhs(rows[i]);
                    if(g_memoscip.collector.dualsol[i]<0){
                        SCIP_CALL( SCIPaggrRowAddRow(set->scip, Trueproofrow, rows[i], g_memoscip.collector.dualsol[i], -1) );
                    }
                    else{
                        SCIP_CALL( SCIPaggrRowAddRow(set->scip, Trueproofrow, rows[i], g_memoscip.collector.dualsol[i], +1) );
                    }
                }
            }
        }else {
            //cout<<"Row "<<SCIProwGetName(rows[i])<<" has dual value "<<g_memoscip.collector.dualsol[i]<<" and rhs "<<SCIProwGetRhs(rows[i])<<endl;
            assert(g_memoscip.collector.dualsol[i]>=-1e-5);
        }
    }
   

    diffy_i=newy_i;
    assert(g_memoscip.AuxVar.size()==g_memoscip.VarUb.size());
    assert(g_memoscip.AuxVar.size()==g_memoscip.CurrentVarIdx.size());
    if(nodeNumber>1){
        //If g_memoscip.nodeInfoMap[nodeNumber].y_i.empty() then it is the first time we analyse the node 
        //Otherwise it is the second time aka we solved the LP applied NC and BP and resolved the LP 
        if(g_memoscip.nodeInfoMap[nodeNumber].y_i.empty()){
            //If variable is fixed to 0 then dual value is 0. 
            for (size_t i = 0; i < g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i.size(); i++){
                if(g_memoscip.VarUb[i]==0){
                    assert(newy_i[i]==0);
                    diffy_i[i]=-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i[i];
                } 
                else{
                    diffy_i[i]-=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i[i];
                    assert(g_memoscip.VarLb[g_memoscip.AuxVar[i]]==0 || SCIPisDualfeasZero(scip,diffy_i[i]));
                }
            }
        }else{
            for (size_t i = 0; i < g_memoscip.nodeInfoMap[nodeNumber].y_i.size(); i++){
                if(g_memoscip.VarUb[i]==0){
                    assert(newy_i[i]==0);
                    diffy_i[i]=-g_memoscip.nodeInfoMap[nodeNumber].y_i[i];
                } 
                else{
                    diffy_i[i]-=g_memoscip.nodeInfoMap[nodeNumber].y_i[i];
                }
            }
        }
    }

    vector<SCIP_Real> TempCoef(g_memoscip.CurrentVarIdx.size(),0.0);
    vector<SCIP_Real> TrueTempCoef;
    if(g_options.debug){
        TrueTempCoef=vector<SCIP_Real>(g_memoscip.CurrentVarIdx.size(),0.0);
    }

    //Create the dual proof coef
    for(size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i )
    {
        int currv = g_memoscip.CurrentVarIdx[i];
        SCIP_Real coef = SCIPaggrRowGetProbvarValue(proofrow, currv);
        TempCoef[i]+=coef;
        if(g_options.debug){
            TrueTempCoef[i]+=SCIPaggrRowGetProbvarValue(Trueproofrow, currv);
        }
        if(!SCIPisDualfeasZero(scip,diffy_i[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
            TempCoef[i]+=diffy_i[i];
            TempCoef[g_memoscip.AuxVar[i]]+=diffy_i[i];
            sumdual+=diffy_i[i];
        }
        if(g_options.debug){
            if(!SCIPisDualfeasZero(scip,newy_i[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
                    Truesumdual+=newy_i[i];
                    TrueTempCoef[i]+=newy_i[i];
                    TrueTempCoef[g_memoscip.AuxVar[i]]+=newy_i[i];
                }
        }
    } 
    if(g_options.verbosity > 1){
        cout<<"DUAL proof"<<endl;
        for(size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i )
        {
            SCIP_Real coef = TempCoef[i];
            if(SCIPisDualfeasZero(scip,coef))
                continue;
            if(coef<0){
                cout<<coef<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
            }else
                cout<<"+"<<coef<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
        }
        cout<<" <= "<<sumdual<<endl;
    }

    for(size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i )
    {
        int currv = g_memoscip.CurrentVarIdx[i];
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
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(coef);
        g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.push_back(i);
        g_memoscip.nodeInfoMap[nodeNumber].DualProofVars.push_back(g_memoscip.AllVars[currv]);
        if(g_options.debug){
            if(i<parentredcost.size()){
                if( g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofInds.empty()){
                    if(g_options.verbosity > 2){
                        cout<<"COEF IN PROOF "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" "<<coef<<" vs parent red "<<parentredcost[i]<<" new red "<<parentredcost[i]-coef<< " LB "<<g_memoscip.VarLb[i]<<" UB "<<g_memoscip.VarUb[i]<<endl;
                        cout<<parentredcost[i]-coef<<endl;
                    }
                    assert(parentredcost[i]-coef>-9.9e-1|| 0==g_memoscip.VarUb[i]);
                }
            }else {
                assert(SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])[2]=='s');
                assert(coef<1e-9);
            }
        }
    }
    
    g_memoscip.nodeInfoMap[nodeNumber].DualSol=g_memoscip.collector.dualsol;
    g_memoscip.nodeInfoMap[nodeNumber].y_i=newy_i;
    g_memoscip.nodeInfoMap[nodeNumber].redcosts = redcost;
    g_memoscip.nodeInfoMap[nodeNumber].lb=g_memoscip.collector.objval;
    g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=sumdual;
    g_memoscip.nodeInfoMap[nodeNumber].Status="FEASIBLE"; 
    g_memoscip.nodeInfoMap[nodeNumber].lb_increase=sumdual;
    if(g_options.verbosity > 0){
        cout<<g_memoscip.collector.objval<<" "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb<<" Truesumdual "<<Truesumdual<<" "<<Truesumdual-g_memoscip.collector.objval<<endl;
    }
    if(g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofInds.empty()){
        if(g_options.verbosity > 0){
            cout<<"DIFF Objval "<<g_memoscip.collector.objval-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb<<" versus sumdiffdual "<<sumdual<<endl;
        }
        assert(fabs(g_memoscip.collector.objval-sumdual-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb) < 1e-2|| g_memoscip.collector.objval>= SCIPgetPrimalbound(scip));
    }else{
        if(g_options.verbosity > 0){
            cout<<"DIFF Objval "<<g_memoscip.collector.objval-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb<<" versus sumdiffdual "<<sumdual+g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofCoef.second<<endl;
        }
        assert(fabs(g_memoscip.collector.objval-sumdual-g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofCoef.second-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb) < 1e-2 || g_memoscip.collector.objval>= SCIPgetPrimalbound(scip));
    }

    SCIPaggrRowFree(scip, &proofrow);
        
    //TO TEST ::
    if(g_options.debug && !SCIPinProbing(scip))
    { 
        SCIP_CALL( SCIPstartProbing(scip) );
        SCIP_Bool lperror, cutoff;
        for (size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i)
        {
            SCIP_CALL( SCIPchgVarObjProbing(scip, g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]], redcost[i]) );
        }
        SCIP_CALL( SCIPflushLP(scip) );
        SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
        char nameprob[50];
        //SCIP_CALL( SCIPwriteLP(scip, "mort.lp") ); 

        snprintf(nameprob, 50, "probing_%lld.lp", nodeNumber);
        //SCIP_CALL( SCIPwriteTransProblem(scip, nameprob,"lp",FALSE) );
        //cout<<"Diving LP solved "<<cutoff<<" "<<lperror<<" "<<SCIPgetLPObjval(scip)<< endl;
        SCIP_LPSOLSTAT st = SCIPgetLPSolstat(scip);
        assert(st==SCIP_LPSOLSTAT_OPTIMAL);
        assert(SCIPgetLPObjval(scip)<1e-3);
        SCIP_CALL( SCIPendProbing(scip));
    }
    
    return SCIP_OKAY;
}


SCIP_RETCODE gatherLPInfeasibilityCertificate(
    SCIP* scip,
    SCIP_Longint nodeNumber)
{
    if(!g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.empty()){
        g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofCoef=g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef;
        g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofInds=g_memoscip.nodeInfoMap[nodeNumber].DualProofInds;
        g_memoscip.nodeInfoMap[nodeNumber].FirstDualProofVars=g_memoscip.nodeInfoMap[nodeNumber].DualProofVars;
        g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.clear();
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.clear();
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=0;
        g_memoscip.nodeInfoMap[nodeNumber].DualProofVars.clear();
        g_memoscip.nodeInfoMap[nodeNumber].DualSol.clear();
        g_memoscip.nodeInfoMap[nodeNumber].y_i.clear();
        g_memoscip.nodeInfoMap[nodeNumber].redcosts.clear();
    }
    //If INF is true then we collect the infeasibility certificate from farkas rows
    //Otherwise return a basic infeasibility certificate from the variables that are fixed to 0
    if(g_options.INF){
        SCIP_SET* set= scip->set;  
        SCIP_LPI* lpi = scip->lp->lpi;
        SCIP_ROW** rows = SCIPgetLPRows(scip);
        SCIP_AGGRROW* proofrow = nullptr;
        SCIP_CALL(SCIPaggrRowCreate(scip, &proofrow));
        int nrows = SCIPgetNLPRows(scip);
        int ncols = SCIPgetNLPCols(scip);
        vector<SCIP_Real> rowFarkas(nrows+ncols);
        SCIP_CALL( SCIPlpiGetDualfarkas(lpi, rowFarkas.data()));
        for(int i = 0; i < nrows+ncols; ++i)
        {
            if(!SCIPisDualfeasZero(scip,rowFarkas[i])){
                if(g_options.verbosity > 1){
                    if(i<nrows)
                        cout<<"Row "<<i<<" "<<SCIProwGetName(rows[i])<<" "<<rowFarkas[i]<<endl;
                    else
                        cout<<"Col "<<i-nrows<<" "<<SCIPvarGetName(g_memoscip.AllVars[i-nrows])<<" "<<rowFarkas[i]<<endl;
                }
                if(SCIProwGetRhs(rows[i])==SCIProwGetLhs(rows[i]))
                    SCIP_CALL( SCIPaggrRowAddRow(set->scip, proofrow, rows[i], rowFarkas[i], +1) );
                else
                    assert(SCIProwGetName(rows[i])[0]=='s' && SCIProwGetName(rows[i])[1]=='c');
            }
        }
        
        vector<SCIP_Real> Coefs(g_memoscip.CurrentVarIdx.size(),0.0);
        SCIP_Real rhs = SCIPaggrRowGetRhs(proofrow);
        SCIP_Real coef;
        for (size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i)
        {
            int currv = g_memoscip.CurrentVarIdx[i];
            int curraux= g_memoscip.CurrentVarIdx[g_memoscip.AuxVar[i]];
            coef= SCIPaggrRowGetProbvarValue(proofrow, currv);
            if(!SCIPisDualfeasZero(scip,coef)){
                //cout<<i<<" "<<SCIPvarGetName(g_memoscip.AllVars[currv])<<" "<<coef<<" bounds "<<SCIPvarGetLbLocal(g_memoscip.AllVars[currv])<<" "<<SCIPvarGetUbLocal(g_memoscip.AllVars[currv])<<endl;
                //cout<<i<<" "<<SCIPvarGetName(g_memoscip.AllVars[curraux])<<" bounds "<<SCIPvarGetLbLocal(g_memoscip.AllVars[curraux])<<" "<<SCIPvarGetUbLocal(g_memoscip.AllVars[curraux])<<endl;
                //assert(g_memoscip.VarUb[i]==0 || g_memoscip.VarLb[i]==1 || coef < 0);
                //If variable is assigned to 0 then add the the var to the explanation 
                //If variable is assigned to 1 then add the auxiliate var to the explanation
                //If variable is unassigned add the var if the coef is <0 otherwise add the auxiliate var
                if(g_memoscip.VarUb[i]==0 || (coef < 0 && g_memoscip.VarLb[i]==0)){
                    Coefs[i]+=coef;
                }else{
                    Coefs[g_memoscip.AuxVar[i]]-=coef;
                    rhs-=coef;
                }
            }
        }
        if(g_options.verbosity > 1){
            cout<<"INFEASIBLE PROOF "<< rhs<<endl;
            for(size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i)
            {
                int currv = g_memoscip.CurrentVarIdx[i];
                if(!SCIPisDualfeasZero(scip,Coefs[i])){
                    cout<<Coefs[i]<<" "<<SCIPvarGetName(g_memoscip.AllVars[currv])<<endl;
                    assert(Coefs[i]<g_memoscip.TOP);
                }
            }
        }

        
        for(size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i )
        {
            int currv = g_memoscip.CurrentVarIdx[i];
            SCIP_Real coef = Coefs[i];
            if(SCIPisDualfeasZero(scip,coef))
                continue;
            if(coef>0){
                assert(g_memoscip.VarUb[i]==0);
                g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(g_memoscip.TOP);
            }
            else{
                g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(-g_memoscip.TOP);
            }
            g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.push_back(i);
            g_memoscip.nodeInfoMap[nodeNumber].DualProofVars.push_back(g_memoscip.AllVars[currv]);
        }
        g_memoscip.nodeInfoMap[nodeNumber].lb=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].Status="INFEASIBLE"; 
        g_memoscip.nodeInfoMap[nodeNumber].lb_increase=g_memoscip.TOP;
        assert(rhs>=0);
        SCIPaggrRowFree(scip, &proofrow);

    }else{
        for (size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); ++i)
        {
            int currv = g_memoscip.CurrentVarIdx[i];
            if(g_memoscip.VarUb[i]==0){
                g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(g_memoscip.TOP);
                g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.push_back(i);
                g_memoscip.nodeInfoMap[nodeNumber].DualProofVars.push_back(g_memoscip.AllVars[currv]);
            }
        }
        g_memoscip.nodeInfoMap[nodeNumber].lb=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].Status="INFEASIBLE"; 
        g_memoscip.nodeInfoMap[nodeNumber].lb_increase=g_memoscip.TOP;
    }
   
    return SCIP_OKAY;
}
