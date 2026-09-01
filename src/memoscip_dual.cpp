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
    int nvars = g_memoscip.AllVars.size();
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
    cout<<CDP1.size()<<" "<<g_memoscip.nodeInfoMap[child0].DualProofInds.size()<<endl;*/
    assert(nnz1==g_memoscip.nodeInfoMap[child0].DualProofInds.size());
    assert(nnz1==inds1.size());
    //cout<<"child 0 with "<<nnz1<<" variables and rhs "<< g_memoscip.nodeInfoMap[child0].DualProofCoef.second<<" branch var "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[child0].BranchVar]])<<endl;
    for( i = 0; i < nnz1; ++i )
    {
        idx = inds1[i];
        SCIP_Real coef = CDP1[i];
        if(!SCIPisDualfeasZero(scip,coef) && g_memoscip.nodeInfoMap[child0].DualProofInds[i] != g_memoscip.AuxVar[SCIPvarGetIndex(g_memoscip.AllVars[g_memoscip.nodeInfoMap[child0].BranchVar])-g_memoscip.nOriginalVar])
            FusCoef[idx]=max(coef,0.0);
    }

    
    size_t nnz2 = CDP2.size();
    vector<int> inds2 = g_memoscip.nodeInfoMap[child1].DualProofInds;
    //printf("AGGRROW: rhs=%g, nnz=%d\n", SCIPaggrRowGetRhs(DP2), nnz2);
    //cout<<"START CDP2 "<<endl;
    //cout<<"child 1 with "<<nnz2<<" variables and rhs "<<g_memoscip.nodeInfoMap[child1].DualProofCoef.second<<" branch var "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[child1].BranchVar]])<<endl;
    assert(nnz2==g_memoscip.nodeInfoMap[child1].DualProofInds.size());
    assert(nnz2==inds2.size());
    for( i = 0; i < nnz2; ++i )
    {
        idx = inds2[i];
        SCIP_Real coef = CDP2[i];
        if(!SCIPisDualfeasZero(scip,coef) && g_memoscip.nodeInfoMap[child1].DualProofInds[i] != g_memoscip.AuxVar[SCIPvarGetIndex(g_memoscip.AllVars[g_memoscip.nodeInfoMap[child1].BranchVar])-g_memoscip.nOriginalVar]){
            FusCoef[idx]=max(FusCoef[idx],coef);
        }
    }
    
    return FusCoef;
}
SCIP_RETCODE FinishLearn(SCIP* scip, SCIP_Longint nodenumber){
    vector<long> children = g_memoscip.nodeInfoMap[nodenumber].children;
    vector<SCIP_Real> ResFus(g_memoscip.CurrentVarIdx.size(),0.0);
    SCIP_Longint parentNumber = g_memoscip.nodeInfoMap[nodenumber].parent;
    
    SCIP_Real rhsFus;
    cout<<"Finishing node "<<nodenumber<<endl;
    
    cout<<"Deleting node "<<nodenumber<<" "<<g_memoscip.nodeInfoMap[nodenumber].lb<<endl;
    cout<<"number of children "<<g_memoscip.nodeInfoMap[nodenumber].children.size()<<endl;
    cout<<"Bracnh var "<<g_memoscip.nodeInfoMap[nodenumber].BranchVar<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[g_memoscip.nodeInfoMap[nodenumber].BranchVar]])<<endl;
    cout<<"Curnnet node lb increase "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" versus children "<<endl;
        
    
    if(children.empty())
        return SCIP_OKAY;

    if(children.size()==1){
        assert(g_memoscip.nodeInfoMap.contains(children[0]));
        if(g_options.verbosity > 2)
            cout<<"Single child "<<children[0]<<" lb "<<g_memoscip.nodeInfoMap[children[0]].lb<<" versus parent lb "<<g_memoscip.nodeInfoMap[nodenumber].lb<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.first.size()<<" "<<g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second<<" incr lb "<<g_memoscip.nodeInfoMap[children[0]].lb_increase<<endl;
        assert(g_memoscip.nodeInfoMap[children[0]].lb<g_memoscip.TOP);
        assert(fabs(g_memoscip.nodeInfoMap[children[0]].DualProofCoef.second-g_memoscip.nodeInfoMap[children[0]].lb_increase)<1e-5);
    }
    if(children.size()==2){
        FinishLearn(scip,children[0]);
        FinishLearn(scip,children[1]);
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
    }

    
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
            cout<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[idx]])<< " "<< idx<<" "<<g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second<<" "<<ceil(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second+g_memoscip.nodeInfoMap[nodenumber].lb)<<" "<<SCIPgetPrimalbound(scip)<<endl;
        //assert(ceil(g_memoscip.nodeInfoMap[nodenumber].NodeConsistency[i].second+g_memoscip.nodeInfoMap[nodenumber].lb)>=SCIPgetPrimalbound(scip));
    }

   
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

    if(rhs_summed>0 && g_options.learning==1 && rhs_summed < g_memoscip.TOP){
        char name[100];

        snprintf(name, 100, "cut_%lld_%d", nodenumber,SCIPgetDepth(scip));
        

        for (size_t i = 0; i < ResFus.size(); i++)
        {
            if(!SCIPisDualfeasZero(scip,ResFus[i])){
                NewInds.push_back(i);
                NewCoefs.push_back(ResFus[i]);
                NewVar.push_back(g_memoscip.AllVars[i]);
                //if(i<g_memoscip.nodeInfoMap[parentNumber].redcosts.size())
                //  cout<<SCIPvarGetName(g_memoscip.AllVars[i])<< " i "<< i<<" resfus "<<ResFus[i]<<" parentrc "<<g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<" diff "<<ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<<endl;
                assert( i >= g_memoscip.nodeInfoMap[parentNumber].redcosts.size() || ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<1e-3 || g_memoscip.VarUb[i] == g_memoscip.VarLb[i]);
                //assert(ResFus[i]<g_memoscip.TOP);
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
                        rhs_summed,  // lhs (pour inégalité)
                        rhs_summed,&g_memoscip.CutQueue) ); 
        if(g_options.verbosity > 2)
            SCIPprintQueuedCuts(scip, &g_memoscip.CutQueue);
    }

    
    assert(abs(g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second-g_memoscip.nodeInfoMap[nodenumber].lb_increase) < 1e-9);
    assert(NewCoefs.size()<=g_memoscip.CurrentVarIdx.size());
    g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.second=rhs_summed;
    g_memoscip.nodeInfoMap[nodenumber].DualProofInds=NewInds;
    g_memoscip.nodeInfoMap[nodenumber].DualProofCoef.first=NewCoefs;
    assert(g_memoscip.nodeInfoMap[nodenumber].lb_increase<=rhs_summed || !g_memoscip.nodeInfoMap[nodenumber].NodeConsistency.empty());
    g_memoscip.nodeInfoMap[nodenumber].lb_increase=rhs_summed;
    //cout<<" plb "<<g_memoscip.nodeInfoMap[parentNumber].lb<<" lbinc "<<g_memoscip.nodeInfoMap[nodenumber].lb_increase<<" ub "<<SCIPgetPrimalbound(scip)<<endl;
    assert(ceil(g_memoscip.nodeInfoMap[parentNumber].lb + g_memoscip.nodeInfoMap[nodenumber].lb_increase+1e-5) >= SCIPgetPrimalbound(scip));
    if(!children.empty()){
        g_memoscip.nodeInfoMap.erase(children[0]);
        if(children.size()==2)  
            g_memoscip.nodeInfoMap.erase(children[1]);
        g_memoscip.nodeInfoMap[nodenumber].children.clear();
    } 
    
    return SCIP_OKAY;

}



SCIP_RETCODE AnalyseDualSol(SCIP* scip, SCIP_Longint nodeNumber){
    SCIP_LPI* lpi;
    SCIP_CALL(SCIPgetLPI(scip, &lpi));

    int ncols = SCIPgetNLPCols(scip);
    int nrows = SCIPgetNLPRows(scip);
     if((g_options.output>1 && nodeNumber >= g_options.output)){
        cout<<"WRITE PORVLME"<<endl;
        char nameprobsol[50];
        snprintf(nameprobsol, 50, "solved_%lld.lp", nodeNumber);
        SCIP_CALL( SCIPwriteLP(scip, nameprobsol) ); 
        snprintf(nameprobsol, 50, "problem_%lld.lp", nodeNumber);
        SCIP_CALL( SCIPwriteTransProblem(scip, nameprobsol,"lp",FALSE) );
    }
    //cout<<"SSolving LP "<<nodeNumber<<" with "<<ncols<<" columns and "<<nrows<<" rows"<<endl;
    
    vector<SCIP_Real> parentredcost;
    
    vector<SCIP_Real> IndexedDualSol(g_memoscip.nbctr,0.0);
    
    vector<SCIP_Real> diffy_i(g_memoscip.AllVars.size(),0.0);
    vector<SCIP_Real> newy_i(g_memoscip.AllVars.size(),0.0);

    vector<SCIP_Real> redcost(g_memoscip.AllVars.size(),0.0);

    g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.clear();
    g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.clear();
    g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=0.0;
    if( g_memoscip.collector.dualsol.empty() )
        g_memoscip.collector.init(scip);

    if( SCIPgetNLPRows(scip) != (int)g_memoscip.collector.dualsol.size() || SCIPgetNLPCols(scip) != (int)g_memoscip.collector.redcost.size() )
    {
        g_memoscip.collector.init(scip);  // rebuild mapping + resize buffers
    }

    SCIP_CALL(g_memoscip.collector.collect(scip));

    SCIP_COL** cols = SCIPgetLPCols(scip);
    SCIP_ROW** rows = SCIPgetLPRows(scip);
    SCIP_LPSOLSTAT stat = SCIPgetLPSolstat(scip);
    assert(stat==SCIP_LPSOLSTAT_OPTIMAL || stat==SCIP_LPSOLSTAT_OBJLIMIT);
    assert(g_memoscip.collector.dualsol.size()==nrows);
    assert(g_memoscip.collector.redcost.size()==ncols);
    if(g_options.debug){
        cout<<"Size of  g_memoscip.collector.redcost "<<g_memoscip.collector.redcost.size()<<" vs "<<ncols<<endl;
        for(int i =0; i < ncols; i++){
            SCIP_VAR* var = SCIPcolGetVar(cols[i]);
            SCIP_Real rc  = g_memoscip.collector.redcost[i];  
            //cout<<"i "<<i<<" varindex "<<SCIPvarGetIndex(var)<<" varname "<<SCIPvarGetName(var)<<" redcost "<<SCIPgetVarRedcost(scip,var)<<" rc "<<rc<<endl;
            assert(SCIPinDive(scip) || SCIPisFeasEQ(scip, SCIPgetVarRedcost(scip, var), rc));
        }


        /*for( size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); i++ )
        {
            int currv=g_memoscip.CurrentVarIdx[i];
            SCIP_VAR* var = SCIPcolGetVar(cols[currv]);
            SCIP_Real rc  = g_memoscip.collector.redcost[currv];  
            //if(rc!=0)
            //    cout<<" "<<SCIPvarGetName(var)<<" "<<g_memoscip.collector.redcost[currv]<<endl;
            //cout<<SCIPgetVarRedcost(scip,g_memoscip.AllVars[currv])<<" "<<g_memoscip.collector.redcost[currv]<<endl;
            cout<<i<<" "<<SCIPvarGetIndex(var)<<" "<<SCIPvarGetName(var)<<endl;
            cout<<" "<<SCIPvarGetName(var)<<" "<<SCIPvarGetName(g_memoscip.AllVars[currv])<<endl;
            assert(var == g_memoscip.AllVars[currv]);

            assert(SCIPisFeasEQ(scip, SCIPgetVarRedcost(scip, g_memoscip.AllVars[currv]), g_memoscip.collector.redcost[currv]));
        }*/
    }




    if (nodeNumber==1)
    {
        parentredcost=g_memoscip.nodeInfoMap[0].redcosts;
    }else{
        parentredcost=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].redcosts;
    }


    //We obtain the dual values corresponding to bound constraints by looking at the redcosts
    //If the variable is fixed to 0 then the dual value is 0. 
    //If the variable is fixed to 1 then the dual value is the red cost (bound constraint is equality constaraint and dual value is unsigned)
    //If the variable is not fixed then the dual value is != 0 only if redcost is negative (bound constraint is <= constraint and dual value is negative)
    SCIP_Real redc;

    for(int i=0; i < ncols; i++){
        SCIP_VAR* var = SCIPcolGetVar(cols[i]);
        redc=g_memoscip.collector.redcost[i];
        int idx =g_memoscip.CurrentVarIdx[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar];
        redcost[idx]=redc;
        if(!SCIPisDualfeasZero(scip,redc)){
            //cout<<"idx " << idx <<" SCIPvarGetIndex(var) "<<SCIPvarGetIndex(var)<<" "<<SCIPvarGetName(g_memoscip.AllVars[idx])<<" "<<g_memoscip.VarLb[idx]<<" "<<g_memoscip.VarUb[idx]<<" redc "<<redc<<endl;
            if(g_memoscip.VarUb[idx]==0)
                newy_i[idx]=0;
            else if(g_memoscip.VarLb[idx]==1){
                newy_i[idx]=redc;
            }
            else{
                newy_i[idx]=min(redc,0.0);
            }
            assert((SCIPvarGetName(g_memoscip.AllVars[idx])[0]!='s' && SCIPvarGetName(g_memoscip.AllVars[idx])[2]!='s') || newy_i[idx]==0);
        }
    }

    /*for (size_t i = 0; i < g_memoscip.CurrentVarIdx.size(); i++)
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
    }*/
    
    SCIP_Real sumdual=0.0;
    SCIP_Real Truesumdual=0.0;
    SCIP_SET* set= scip->set;  
    SCIP_AGGRROW* proofrow = nullptr;
    SCIP_AGGRROW* Trueproofrow = nullptr;

    SCIP_CALL(SCIPaggrRowCreate(scip, &proofrow));

    if(g_options.debug){
        SCIP_CALL(SCIPaggrRowCreate(scip, &Trueproofrow));
    }


    //We iterate over all rows and add the rows multiplied by the dual solution to the proof row
    for (int i = 0; i < nrows; ++i)
    {
        SCIP_Real dual = g_memoscip.collector.dualsol[i];
        int idx = SCIProwGetIndex(rows[i]);
        IndexedDualSol[idx] = dual;
        if(idx < (int)g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol.size()){
            dual-=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol[idx];
        }
        
        //if size of a row is greater than 1
        if(SCIProwGetNNonz(rows[i])>1){            
            if(g_options.verbosity > 2  && dual!=0){
                if(idx<(int)g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol.size())
                    cout << "Row " << SCIProwGetName(rows[i]) << " Idew "<< SCIProwGetIndex(rows[i]) <<" i "<<i<<" Age: " << SCIProwGetAge(rows[i])<<" "<<" dualsol "<<g_memoscip.collector.dualsol[i]<<" parent "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].DualSol[idx]<<" "<<dual<<" sumdual "<<sumdual<<endl;
                else 
                    cout << "Row " << SCIProwGetName(rows[i]) <<" Idew "<< SCIProwGetIndex(rows[i])<<" i "<<i<<" dualsol "<<g_memoscip.collector.dualsol[i]<<" "<<dual<<" Age: " << SCIProwGetAge(rows[i])<<endl;
            }
            assert(fabs(dual)<1e6);
            if(fabs(SCIProwGetRhs(rows[i])*dual)  > 1e-7   && SCIProwGetRhs(rows[i]) != SCIPinfinity(scip)){
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



    if(g_options.debug){
        for(int i = 5; i < nrows; ++i)
        {
            SCIP_ROW* row;
            SCIP_COL** cols;
            int nnonz;

            row = rows[i];
            cols = SCIProwGetCols(row);
            nnonz = SCIProwGetNNonz(row);

            if( nnonz > 2 )
            {
                SCIP_VAR* var;

                var = SCIPcolGetVar(cols[nnonz-1]);
                //if(g_options.verbosity > 2)
                    //cout << "Row " << SCIProwGetName(row) << " Idew "<< SCIProwGetIndex(row) <<" last var "<<SCIPvarGetName(var) << " status "<<SCIPvarGetStatus(var)<<endl;
                assert(SCIPvarGetName(var)[0] == 's');
                assert(SCIPvarIsInLP(var) || g_memoscip.collector.dualsol[i]==0);
            }
        }
    }

    diffy_i=newy_i;
    assert(g_memoscip.AuxVar.size()==g_memoscip.CurrentVarIdx.size());
    if(nodeNumber>1){
        //If g_memoscip.nodeInfoMap[nodeNumber].y_i.empty() then it is the first time we analyse the node 
        //Otherwise it is the second time aka we solved the LP applied NC and BP and resolved the LP 
        //If variable is fixed to 0 then dual value is 0. 
        for (size_t i = 0; i < g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i.size(); i++){
            if(g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i[i]!=0)
            {
                //cout<<"Diff y_i "<<diffy_i[i]<<" Parent dual y_i "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i[i]<<" for variable "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<endl;
                if(g_memoscip.VarUb[i]==0){
                    assert(newy_i[i]==0);
                    diffy_i[i]=-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i[i];
                } 
                else{
                    diffy_i[i]-=g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].y_i[i];
                    assert(g_memoscip.VarLb[g_memoscip.AuxVar[SCIPvarGetIndex(g_memoscip.AllVars[i])-g_memoscip.nOriginalVar]]==0 || SCIPisDualfeasZero(scip,diffy_i[i]));
                }
            }
        }
       
    }

    vector<SCIP_Real> TempCoef(g_memoscip.AllVars.size(),0.0);
    vector<SCIP_Real> TrueTempCoef;
    if(g_options.debug){
        TrueTempCoef=vector<SCIP_Real>(g_memoscip.AllVars.size(),0.0);
    }


    //Print proofrow
    /*int nnz = SCIPaggrRowGetNNz(proofrow);
    int* inds = SCIPaggrRowGetInds(proofrow);
    if(g_options.verbosity > 2){
        cout<<"Proof row "<<endl;
        for(int i = 0; i < nnz; ++i )
        {
            cout<<inds[i]<<" "<<SCIPaggrRowGetProbvarValue(Trueproofrow, inds[i])<<endl;
        }
    }*/
    //Create the dual proof coef
    for(size_t i = 0; i < g_memoscip.AllVars.size(); ++i )
    {
        SCIP_VAR* var = g_memoscip.AllVars[i];
        SCIP_Real coef = SCIPaggrRowGetProbvarValue(proofrow, SCIPvarGetProbindex(var));
        //cout<<"proof coef "<<i <<" "<<SCIPvarGetName(var)<<" "<<SCIPvarGetIndex(var)<<" "<<coef<<" "<<endl;

        TempCoef[i]+=coef;
        if(g_options.debug){
            TrueTempCoef[i]+=SCIPaggrRowGetProbvarValue(Trueproofrow, SCIPvarGetProbindex(var));
        }
        if(!SCIPisDualfeasZero(scip,diffy_i[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
            TempCoef[i]+=diffy_i[i];
            TempCoef[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]+=diffy_i[i];
            sumdual+=diffy_i[i];
        }
        if(g_options.debug){
            if(!SCIPisDualfeasZero(scip,newy_i[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
                //cout<<SCIPvarGetName(var)<<" "<<newy_i[i]<<" "<<Truesumdual<<endl;
                Truesumdual+=newy_i[i];
                TrueTempCoef[i]+=newy_i[i];
                TrueTempCoef[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]+=newy_i[i];
            }
        }
    } 

    assert(g_memoscip.AllVars.size() == TempCoef.size());
    if(g_options.verbosity > 2){
        cout<<"DUAL proof"<<endl;
        for(size_t i = 0; i < g_memoscip.AllVars.size(); ++i )
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

    if(g_options.verbosity > 2 && g_options.debug){
        cout<<"TRUEDUAL proof"<<endl;
        for(size_t i = 0; i < g_memoscip.AllVars.size(); ++i )
        {
            SCIP_Real coef = TrueTempCoef[i];
            if(SCIPisDualfeasZero(scip,coef))
                continue;
            if(coef<0){
                cout<<coef<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
            }else
                cout<<"+"<<coef<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[i]])<<" ";
        }
        cout<<" <= "<<Truesumdual<<endl;
    }
    assert(!parentredcost.empty());
    for(size_t i = 0; i < g_memoscip.AllVars.size(); ++i )
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
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(coef);
        g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.push_back(i);
        
        if(g_options.verbosity > 3){
            if(i<parentredcost.size()){
                cout<<"COEF IN PROOF "<<SCIPvarGetName(g_memoscip.AllVars[i])<<" "<<coef<<" vs parent red "<<parentredcost[i]<<" new red "<<parentredcost[i]-coef<< " LB "<<g_memoscip.VarLb[i]<<" UB "<<g_memoscip.VarUb[i]<<endl;
                cout<<parentredcost[i]-coef<<endl;
                assert(parentredcost[i]-coef>-9.9e-1|| 0==g_memoscip.VarUb[i] || !g_memoscip.nodeInfoMap[nodeNumber].DualSol.empty());
            }else {
                cout<<SCIPvarGetName(g_memoscip.AllVars[i])<<" "<<coef<<endl;
                assert(SCIPvarGetName(g_memoscip.AllVars[i])[0]=='s' || SCIPvarGetName(g_memoscip.AllVars[i])[2]=='s');
                assert(coef<1e-9);
            }
        }
    }
    g_memoscip.nodeInfoMap[nodeNumber].DualSol=IndexedDualSol;
    g_memoscip.nodeInfoMap[nodeNumber].y_i=newy_i;
    g_memoscip.nodeInfoMap[nodeNumber].redcosts = redcost;
    assert(!redcost.empty());
    g_memoscip.nodeInfoMap[nodeNumber].lb=g_memoscip.collector.objval;
    g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=sumdual;
    g_memoscip.nodeInfoMap[nodeNumber].Status="FEASIBLE"; 
    g_memoscip.nodeInfoMap[nodeNumber].lb_increase=sumdual;
    if(g_options.verbosity > 1){
        cout<<g_memoscip.collector.objval<<" "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb<<" Truesumdual "<<Truesumdual<<" "<<Truesumdual-g_memoscip.collector.objval<<endl;
    }

    if(g_options.verbosity > 1){
        cout<<"DIFF Objval "<<g_memoscip.collector.objval-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb<<" versus sumdiffdual "<<sumdual<<endl;
    }
    assert(fabs(g_memoscip.collector.objval-sumdual-g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].lb) < 1e-1|| g_memoscip.collector.objval>= SCIPgetPrimalbound(scip));
    

    SCIPaggrRowFree(scip, &proofrow);
        
    //TO TEST ::
    if(g_options.debug && !SCIPinDive(scip))
    { 
        SCIP_CALL( SCIPstartProbing(scip) );
        SCIP_Bool lperror, cutoff;
        for (size_t i = 0; i < g_memoscip.AllVars.size(); ++i)
        {
            SCIP_CALL( SCIPchgVarObjProbing(scip, g_memoscip.AllVars[i], redcost[i]) );
        }
        SCIP_CALL( SCIPflushLP(scip) );
        SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
        //SCIP_CALL( SCIPwriteLP(scip, "mort.lp") ); 

        SCIP_CALL( SCIPwriteTransProblem(scip, "probling.lp","lp",FALSE) );
        //cout<<"Diving LP solved "<<cutoff<<" "<<lperror<<" "<<SCIPgetLPObjval(scip)<< endl;
        SCIP_LPSOLSTAT st = SCIPgetLPSolstat(scip);
        assert(st==SCIP_LPSOLSTAT_OPTIMAL);
        assert(SCIPgetLPObjval(scip)<1e-2);
        SCIP_CALL( SCIPendProbing(scip));
    }
    
    return SCIP_OKAY;
}


SCIP_RETCODE gatherLPInfeasibilityCertificate(
    SCIP* scip,
    SCIP_Longint nodeNumber)
{
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
        SCIP_Bool negated;
        for(int i = 0; i < nrows+ncols; ++i)
        {
            if(!SCIPisDualfeasZero(scip,rowFarkas[i])){
                if(g_options.verbosity > 2){
                    if(i<nrows)
                        cout<<"Row "<<i<<" "<<SCIProwGetName(rows[i])<<" "<<rowFarkas[i]<<endl;
                    else
                        cout<<"Col "<<i-nrows<<" "<<SCIPvarGetName(g_memoscip.AllVars[i-nrows])<<" "<<rowFarkas[i]<<endl;
                }
                negated = rowFarkas[i] < 0.0;

                if(SCIProwGetRhs(rows[i])==SCIProwGetLhs(rows[i]))
                    SCIP_CALL( SCIPaggrRowAddRow(set->scip, proofrow, rows[i], rowFarkas[i], negated ? -1 : 1) );
                else
                    assert(SCIProwGetName(rows[i])[0]=='s' && SCIProwGetName(rows[i])[1]=='c');
            }
        }
        
        vector<SCIP_Real> Coefs(g_memoscip.AllVars.size(),0.0);
        //SCIP_Real rhs = SCIPaggrRowGetRhs(proofrow);
        SCIP_Real coef;
        for (size_t i = 0; i < g_memoscip.AllVars.size(); ++i)
        {   
            SCIP_VAR* var = g_memoscip.AllVars[i];
            coef = SCIPaggrRowGetProbvarValue(proofrow, SCIPvarGetProbindex(var));
            //cout<<"SCIPFARKAS "<<SCIPgetVarFarkasCoef(scip,var)<<" vs coef "<<coef<<endl;
            /*if(!SCIPisDualfeasZero(scip,coef)){
                //cout<<i<<" "<<SCIPvarGetName(g_memoscip.AllVars[currv])<<" "<<coef<<" bounds "<<SCIPvarGetLbLocal(g_memoscip.AllVars[currv])<<" "<<SCIPvarGetUbLocal(g_memoscip.AllVars[currv])<<endl;
                //cout<<i<<" "<<SCIPvarGetName(g_memoscip.AllVars[curraux])<<" bounds "<<SCIPvarGetLbLocal(g_memoscip.AllVars[curraux])<<" "<<SCIPvarGetUbLocal(g_memoscip.AllVars[curraux])<<endl;
                //assert(g_memoscip.VarUb[i]==0 || g_memoscip.VarLb[i]==1 || coef < 0);
                //If variable is assigned to 0 then add the the var to the explanation 
                //If variable is assigned to 1 then add the auxiliate var to the explanation
                //If variable is unassigned add the var if the coef is <0 otherwise add the auxiliate var
                if(g_memoscip.VarUb[i]==0 || (coef < 0 && g_memoscip.VarLb[i]==0)){
                    Coefs[i]+=coef;
                }else{
                    if(coef >0 && g_memoscip.VarLb[i]==1)
                        cout<<SCIPvarGetName(g_memoscip.AllVars[i])<<" "<<coef<<" "<<g_memoscip.nodeInfoMap[g_memoscip.nodeInfoMap[nodeNumber].parent].redcosts[i]<<endl;
                    Coefs[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]-=coef;
                    rhs-=coef;
                }
            }*/


            if(!SCIPisDualfeasZero(scip,coef)){
                //Only keep variables assigned to 0 in the proof. 
                if(g_memoscip.VarUb[i]==0 && coef < 0){
                    Coefs[i]+=coef;
                }else if(g_memoscip.VarLb[i]==1 && coef > 0)
                    Coefs[g_memoscip.AuxVar[SCIPvarGetIndex(var)-g_memoscip.nOriginalVar]]-=coef;
                
            }
        }
        if(g_options.debug){
            for(int i = 5; i < nrows; ++i)
            {
                SCIP_ROW* row;
                SCIP_COL** cols;
                int nnonz;

                row = rows[i];
                cols = SCIProwGetCols(row);
                nnonz = SCIProwGetNNonz(row);

                if( nnonz > 2 )
                {
                    SCIP_VAR* var;

                    var = SCIPcolGetVar(cols[nnonz-1]);
                    //if(g_options.verbosity > 2)
                        //cout << "Row " << SCIProwGetName(row) << " Idew "<< SCIProwGetIndex(row) <<" last var "<<SCIPvarGetName(var) << endl;
                    assert(SCIPvarGetName(var)[0] == 's');
                }
        }
    }
        
        for(size_t i = 0; i < g_memoscip.AllVars.size(); ++i )
        {
            SCIP_Real coef = Coefs[i];
            if(SCIPisDualfeasZero(scip,coef))
                continue;
            //if(coef>0){
            assert(g_memoscip.VarUb[i]==0);
            g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(g_memoscip.TOP);
            g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.push_back(i);
            //}
            /*else{
                g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(-g_memoscip.TOP);
            }
             g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.push_back(i);*/
           
        }
        g_memoscip.nodeInfoMap[nodeNumber].lb=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].Status="INFEASIBLE"; 
        g_memoscip.nodeInfoMap[nodeNumber].lb_increase=g_memoscip.TOP;
        //assert(rhs>=0);
        SCIPaggrRowFree(scip, &proofrow);

    }else{
        for (size_t i = 0; i < g_memoscip.AllVars.size(); ++i)
        {
            if(g_memoscip.VarUb[i]==0){
                g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.first.push_back(g_memoscip.TOP);
                g_memoscip.nodeInfoMap[nodeNumber].DualProofInds.push_back(i);
            }
        }
        g_memoscip.nodeInfoMap[nodeNumber].lb=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].DualProofCoef.second=g_memoscip.TOP;
        g_memoscip.nodeInfoMap[nodeNumber].Status="INFEASIBLE"; 
        g_memoscip.nodeInfoMap[nodeNumber].lb_increase=g_memoscip.TOP;
        int nrows = SCIPgetNLPRows(scip);
        SCIP_ROW** rows = SCIPgetLPRows(scip);
        if(g_options.debug){
            for(int i = 5; i < nrows; ++i)
            {
                SCIP_ROW* row;
                SCIP_COL** cols;
                int nnonz;

                row = rows[i];
                cols = SCIProwGetCols(row);
                nnonz = SCIProwGetNNonz(row);

                if( nnonz > 2 )
                {
                    SCIP_VAR* var;

                    var = SCIPcolGetVar(cols[nnonz-1]);
                    if(g_options.verbosity > 7){
                        cout << "Row " << SCIProwGetName(row) << " Idew "<< SCIProwGetIndex(row) <<" last var "<<SCIPvarGetName(var) << " status "<<SCIPvarGetStatus(var)<<endl;
                        cout<<"pose "<<SCIPcolGetLPPos(SCIPvarGetCol(var))<<endl;

                        printf("status       = %d\n", SCIPvarGetStatus(var));
                        printf("inLP         = %d\n", SCIPvarIsInLP(var));
                        printf("column       = %p\n", SCIPvarGetCol(var));
                        printf("lp position  = %d\n", SCIPcolGetLPPos(SCIPvarGetCol(var)));
                        printf("nLPCols      = %d\n", SCIPgetNLPCols(scip));
                    }
                    assert(SCIPvarGetName(var)[0] == 's');
                    //assert(SCIPvarIsInLP(var));
                }
            }
        }
    }
   
    return SCIP_OKAY;
}

