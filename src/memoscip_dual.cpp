#include "memoscip_dual.h"

#include "memoscip_option.h"
#include "memoscip_types.h"
#include "sepa_mycuts.h"

#include "lpi/lpi.h"
#include "scip/cons_linear.h"
#include "scip/cuts.h"
#include "scip/scip.h"
#include "scip/struct_lp.h"
#include "scip/struct_scip.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace std;


/** Adds weight * row to the proof row. */
static SCIP_RETCODE addRowToProof(SCIP* scip, SCIP_AGGRROW* proofrow, SCIP_ROW* row, SCIP_Real weight)
{
    // A row relaxed to one side in a dive is added through its finite side
    if( isDiveRelaxedRow(row) )
    {
        SCIP_CALL( SCIPaggrRowAddRow(scip, proofrow, row, weight, SCIPisInfinity(scip, row->rhs) ? -1 : +1) );
        return SCIP_OKAY;
    }
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


/** Max activity of sum vals[k] vars[k] + constant - rhs over the global bounds, without the term of slack (coef -1). */
static SCIP_Real slackUbFromTerms(SCIP* scip, int nvars, SCIP_VAR** vars, const SCIP_Real* vals, SCIP_Real constant,
    SCIP_Real rhs, SCIP_VAR* slack)
{
    if(SCIPisInfinity(scip, fabs(rhs)) || SCIPisInfinity(scip, fabs(constant)))
        return SCIPinfinity(scip);
    SCIP_Real maxact = constant - rhs;
    SCIP_Real sumabs = fabs(constant) + fabs(rhs);
    bool slackFound = false;
    for(int k = 0; k < nvars; ++k)
    {
        SCIP_VAR* var = vars[k];
        const SCIP_Real a = vals[k];
        if(var == slack){
            assert(a == -1.0);
            slackFound = true;
            continue;
        }
        if(a == 0.0)
            continue;
        SCIP_Real bound;
        if(a > 0){
            bound = SCIPvarGetUbGlobal(var);
            if(SCIPisInfinity(scip, bound) && isSlackVar(var))
                bound = getSlackImpliedUb(scip, allVarsIndex(var));
        }else
            bound = SCIPvarGetLbGlobal(var);
        if(SCIPisInfinity(scip, fabs(bound)))
            return SCIPinfinity(scip);
        maxact += a * bound;
        sumabs += fabs(a * bound);
    }
    if(!slackFound)
        return SCIPinfinity(scip);
    // Safety margin for the rounding errors of the sum
    return maxact + 1e-9 * max(1.0, sumabs);
}

SCIP_Real slackUbFromRow(SCIP* scip, SCIP_ROW* row, SCIP_Real rhs, SCIP_VAR* slack)
{
    const int ncols = SCIProwGetNNonz(row);
    SCIP_COL** cols = SCIProwGetCols(row);
    vector<SCIP_VAR*> vars(ncols);
    for(int k = 0; k < ncols; ++k)
        vars[k] = SCIPcolGetVar(cols[k]);
    return slackUbFromTerms(scip, ncols, vars.data(), SCIProwGetVals(row), SCIProwGetConstant(row), rhs, slack);
}

SCIP_Real getSlackImpliedUb(SCIP* scip, int idx)
{
    vector<SCIP_Real>& ubs = g_memoscip.SlackImpliedUb;
    if(idx < 0 || idx >= (int)ubs.size())
        return SCIPinfinity(scip);
    if(!std::isnan(ubs[idx]))
        return ubs[idx];
    SCIP_CONS* cons = idx < (int)g_memoscip.SlackCons.size() ? g_memoscip.SlackCons[idx] : nullptr;
    if(cons == nullptr)
        return SCIPinfinity(scip);
    SCIP_VAR* slack = g_memoscip.AllVars[idx];
    // Guard against a cycle between cuts (a cut only contains the slacks of earlier cuts)
    ubs[idx] = SCIPinfinity(scip);
    // The rhs of the constraint, not of the row, whose sides are relaxed in dives
    const SCIP_Real rhs = SCIPgetRhsLinear(scip, cons);
    SCIP_ROW* row = SCIPgetRowLinear(scip, cons);
    SCIP_Real ub;
    if(row != nullptr){
        // Exact row of the LP: cached
        ub = slackUbFromRow(scip, row, rhs, slack);
        ubs[idx] = ub;
    }else{
        // Row not built yet: the constraint data gives a valid bound (global bounds of possibly non-active
        // variables), but it is not cached so that the exact row is used once it exists
        ub = slackUbFromTerms(scip, SCIPgetNVarsLinear(scip, cons), SCIPgetVarsLinear(scip, cons),
            SCIPgetValsLinear(scip, cons), 0.0, rhs, slack);
        ubs[idx] = std::numeric_limits<SCIP_Real>::quiet_NaN();
    }
    return ub;
}

void dropNegligibleCoefs(SCIP* scip, vector<SCIP_Real>& coefs, SCIP_Real& rhs)
{
    const SCIP_Real tol = 1e-6 * max(1.0, fabs(rhs));
    for (size_t i = 0; i < coefs.size(); ++i)
    {
        const SCIP_Real c = coefs[i];
        if(c == 0.0 || fabs(c) >= tol)
        {
            assert(fabs(c) > 1e-9 || c==0);
            continue;
        }
    
        // c*x_i >= min(c*lb, c*ub) globally, so dropping c*x_i requires rhs -= max(c*lb, c*ub)
        SCIP_VAR* var = g_memoscip.AllVars[i];
        SCIP_Real bound = c > 0 ? SCIPvarGetUbGlobal(var) : SCIPvarGetLbGlobal(var);
        // A learned cut slack has no upper bound, but its row a_k x - s_k = b_k bounds it
        if(c > 0 && SCIPisInfinity(scip, bound) && isSlackVar(var))
        {
            bound = getSlackImpliedUb(scip, (int)i);
            if(g_options.verbosity > 3)
                cout<<"Drop slack coef "<<SCIPvarGetName(var)<<" "<<c<<" with implied ub "<<bound<<endl;
        }
        if(SCIPisInfinity(scip, fabs(bound)))
            continue;
        rhs -= c * bound;
        coefs[i] = 0.0;
    }
}

/** Dense coefficients (indexed by AllVars) of the dual proof of child. */
static vector<SCIP_Real> denseProof(const NodeInfo& child)
{
    assert(child.DualProofCoef.first.size()==child.DualProofInds.size());
    vector<SCIP_Real> coef(g_memoscip.AllVars.size(),0.0);
    for( size_t i = 0; i < child.DualProofInds.size(); ++i )
        coef[child.DualProofInds[i]]=child.DualProofCoef.first[i];
    return coef;
}

/** Fuses the proofs of the two children of a branching on ResolveVar. In the proof of the child where ResolveVar is
 *  fixed to 1 its auxiliary variables (fixed to 0 by the domain constraint) are dropped, in the proof of the child
 *  where it is fixed to 0 ResolveVar itself is dropped. If one proof had no coefficient on the dropped variables it is
 *  valid at the parent and returned alone, otherwise the fused proof takes the max coefficients and the min rhs. */
vector<SCIP_Real> FusionResolution(SCIP* scip, SCIP_Longint child0, SCIP_Longint child1, SCIP_Real& rhs)
{
    const NodeInfo& info0 = g_memoscip.nodeInfoMap[child0];
    const NodeInfo& info1 = g_memoscip.nodeInfoMap[child1];
    const NodeInfo& infoOne = info0.BranchVal ? info0 : info1;
    const NodeInfo& infoZero = info0.BranchVal ? info1 : info0;
    const int ResolveVar = g_memoscip.CurrentVarIdx[infoOne.BranchVar];
    const vector<int>& ResolveVarOpposite = g_memoscip.AuxVar[infoOne.BranchVar];

    vector<SCIP_Real> ctr1 = denseProof(infoOne);
    vector<SCIP_Real> ctr2 = denseProof(infoZero);

    bool iszero = false;
    bool iszero2 = true;
    if(ctr2[ResolveVar]==0)
        iszero=true;
    else
        ctr2[ResolveVar]=0;
    for(int i : ResolveVarOpposite){
        if(ctr1[i]!=0)
            iszero2=false;
        ctr1[i]=0;
    }
    if(iszero && !iszero2){
        rhs=infoZero.DualProofCoef.second;
        return ctr2;
    }else if(iszero2){
        rhs=infoOne.DualProofCoef.second;
        return ctr1;
    }
    for(size_t i = 0; i < ctr1.size(); ++i)
        ctr2[i]=max(ctr1[i],ctr2[i]);
    // An infeasible child has rhs TOP: the min is the bound of the other child
    rhs=min(infoOne.DualProofCoef.second,infoZero.DualProofCoef.second);
    (void)scip;
    return ctr2;
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

bool isInfeasibilityRhs(SCIP_Real rhs)
{
    return rhs >= 0.5*g_memoscip.TOP;
}

SCIP_RETCODE queueInfeasibilityClause(SCIP* scip, SCIP_Longint nodenumber, vector<SCIP_Real> coefs, SCIP_Real rhs)
{
    assert(isInfeasibilityRhs(rhs));
    // sum TOP*x + sum c*y >= rhs with sum c*ub(y) < rhs implies that some x is 1. Slacks (unbounded) are assumed
    // to stay far below TOP
    [[maybe_unused]] SCIP_Real finitemax=0;
    for (size_t i = 0; i < coefs.size(); ++i)
    {
        if(coefs[i] >= 0.5*g_memoscip.TOP){
            coefs[i] = 1.0;
            continue;
        }
        SCIP_VAR* var = g_memoscip.AllVars[i];
        if(coefs[i] > 0 && !isSlackVar(var))
            finitemax += coefs[i]*SCIPvarGetUbGlobal(var);
        coefs[i] = 0.0;
    }
    assert(finitemax < rhs);
    rhs = 1.0;

    vector<SCIP_VAR*> ClauseVars;
    vector<SCIP_Real> ClauseCoefs;
    for (size_t i = 0; i < coefs.size(); ++i)
    {
        if(coefs[i] != 0.0){
            ClauseVars.push_back(g_memoscip.AllVars[i]);
            ClauseCoefs.push_back(coefs[i]);
        }
    }
    char name[100];
    snprintf(name, 100, "clause_%lld_%d", nodenumber,SCIPgetDepth(scip));
    if(g_options.verbosity > 2)
        cout<<"Adding clause "<<name<<" with "<<ClauseVars.size()<<" variables"<<endl;
    if(!ClauseVars.empty())
        SCIP_CALL( SCIPqueueCut(scip, name, (int)ClauseVars.size(), ClauseVars.data(), ClauseCoefs.data(),
                    rhs, rhs, &g_memoscip.CutQueue) );
    return SCIP_OKAY;
}

SCIP_RETCODE fuseChildrenProofs(SCIP* scip, SCIP_Longint nodenumber, const vector<long>& children,
    vector<SCIP_Real>& ResFus, SCIP_Real& rhsFus, bool& infeasibleFusion)
{
    auto& nodeInfoMap = g_memoscip.nodeInfoMap;
    infeasibleFusion = false;
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

    assert(child1.BranchVar==child0.BranchVar && child1.BranchVal!=child0.BranchVal);
    ResFus=FusionResolution(scip,children[0],children[1],rhsFus);
    // The fused proof proves the node infeasible: its clause version is learned. The proof summed with the node dual
    // proof keeps an infeasibility rhs and is not learned
    infeasibleFusion = isInfeasibilityRhs(rhsFus);
    if(infeasibleFusion)
        SCIP_CALL( queueInfeasibilityClause(scip, nodenumber, ResFus, rhsFus) );
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

void checkSummedRhs(SCIP_Longint nodenumber, const vector<long>& children, SCIP_Real rhs_summed,
    bool infeasibleFusion)
{
    if(!g_options.debug)
        return;
    auto& nodeInfoMap = g_memoscip.nodeInfoMap;
    if(children.size()==2 && infeasibleFusion){
        // The rhs is TOP plus the node lb increase
        assert(rhs_summed - nodeInfoMap[nodenumber].lb_increase > -1e-9 || !nodeInfoMap[nodenumber].NodeConsistency.empty());
    }
    else if(children.size()==2){
        // The fusion may keep a single child proof, whose rhs is at least the min of the children
        assert(rhs_summed - (min(nodeInfoMap[children[1]].lb_increase,nodeInfoMap[children[0]].lb_increase) + nodeInfoMap[nodenumber].lb_increase ) > -1e-9 || !nodeInfoMap[nodenumber].NodeConsistency.empty());
    }
    if(children.size()==1){
        // As in addNodeProof, a node increase that is not positive (e.g. -1e-9 of dual noise) is not added
        assert(abs(rhs_summed - (nodeInfoMap[children[0]].lb_increase + max(0.0, nodeInfoMap[nodenumber].lb_increase) )<1e-9) || !nodeInfoMap[nodenumber].NodeConsistency.empty());
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
    bool infeasibleFusion=false;
    cout<<"Finishing node "<<nodenumber<<endl;

    cout<<"Deleting node "<<nodenumber<<" "<<node.lb<<endl;
    cout<<"number of children "<<node.children.size()<<endl;
    cout<<"Bracnh var "<<node.BranchVar<<" "<<SCIPvarGetName(g_memoscip.AllVars[g_memoscip.CurrentVarIdx[node.BranchVar]])<<" = "<<node.BranchVal<<endl;
    cout<<"Curnnet node lb increase "<<node.lb_increase<<" versus children "<<endl;

    if(children.empty())
        return SCIP_OKAY;

    if(children.size()==1)
        checkSingleChildProof(nodenumber, children[0]);
    if(children.size()==2){
        FinishLearn(scip,children[0]);
        FinishLearn(scip,children[1]);
        SCIP_CALL( fuseChildrenProofs(scip, nodenumber, children, ResFus, rhsFus, infeasibleFusion) );
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

    checkSummedRhs(nodenumber, children, rhs_summed, infeasibleFusion);
    dropNegligibleCoefs(scip, ResFus, rhs_summed);

    if(rhs_summed>0 && g_options.learning==1 && !isInfeasibilityRhs(rhs_summed)){
        char name[100];
        snprintf(name, 100, "cut_%lld_%d", nodenumber,SCIPgetDepth(scip));

        for (size_t i = 0; i < ResFus.size(); i++)
        {
            if(ResFus[i] != 0.0){
                NewInds.push_back(i);
                NewCoefs.push_back(ResFus[i]);
                NewVar.push_back(g_memoscip.AllVars[i]);
                assert( i >= g_memoscip.nodeInfoMap[parentNumber].redcosts.size() || ResFus[i]- g_memoscip.nodeInfoMap[parentNumber].redcosts[i]<1e-3 || g_memoscip.VarUb[i] == g_memoscip.VarLb[i]);
            }
        }

        if(g_options.verbosity > 2)
            cout<<"Adding cut "<<name<<" with rhs "<<rhs_summed<<" and "<<NewVar.size()<< " variables"<<endl;

        assert(NewVar.size()==NewCoefs.size());
        // A proof whose coefficients all cancel (e.g. c*(x + aux) >= c reduced to 0 >= 0) is trivial: its rhs is only rounding noise
        assert(!NewVar.empty() || rhs_summed <= 1e-6);
        if(!NewVar.empty())
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
    assert(node.lb_increase<=rhs_summed+1e-6*max(1.0,fabs(rhs_summed)) || !node.NodeConsistency.empty());
    node.lb_increase=rhs_summed;
    assert(Myceil(g_memoscip.nodeInfoMap[parentNumber].lb + node.lb_increase) >= MyPrimalbound(scip));
    eraseChildren(nodenumber, children);

    return SCIP_OKAY;
}


/** Dual bound sum_r y_r (side_r - const_r) + sum_j min(d_j lb_j, d_j ub_j) of the LP for the duals y and reduced
 *  costs d (local bounds); -infinity if a multiplier beyond the dual feasibility tolerance needs an infinite side or
 *  bound (within the tolerance it is LP noise, e.g. on the unbounded slack of a learned cut, and is ignored). */
static SCIP_Real lagrangianBound(SCIP* scip, const vector<SCIP_Real>& dualsol, const vector<SCIP_Real>& redcost)
{
    SCIP_ROW** rows = SCIPgetLPRows(scip);
    SCIP_COL** cols = SCIPgetLPCols(scip);
    SCIP_Real bound = 0.0;
    for (size_t r = 0; r < dualsol.size(); ++r)
    {
        const SCIP_Real y = dualsol[r];
        if (y == 0.0)
            continue;
        const SCIP_Real lhs = SCIProwGetLhs(rows[r]);
        const SCIP_Real rhs = SCIProwGetRhs(rows[r]);
        const SCIP_Real side = (lhs == rhs || y < 0.0) ? rhs : lhs;
        if (SCIPisInfinity(scip, fabs(side)))
        {
            if (SCIPisDualfeasZero(scip, y))
                continue;
            return -SCIPinfinity(scip);
        }
        bound += y * (side - SCIProwGetConstant(rows[r]));
    }
    for (size_t j = 0; j < redcost.size(); ++j)
    {
        const SCIP_Real d = redcost[j];
        if (d == 0.0)
            continue;
        const SCIP_Real b = d > 0.0 ? SCIPcolGetLb(cols[j]) : SCIPcolGetUb(cols[j]);
        if (SCIPisInfinity(scip, fabs(b)))
        {
            if (SCIPisDualfeasZero(scip, d))
                continue;
            return -SCIPinfinity(scip);
        }
        bound += d * b;
    }
    return bound;
}

static int ufFind(vector<int>& uf, int a)
{
    while (uf[a] != a)
        a = uf[a] = uf[uf[a]];
    return a;
}

static bool isUnfixedCol(SCIP_COL* col)
{
    return SCIPcolGetLb(col) != SCIPcolGetUb(col);
}

/** Connected components of the unfixed columns of the LP, joined by the rows that contain them: uf is a union-find
 *  over the LP column positions, rowComp[r] a column of the component of row r (-1 if all its columns are fixed). */
static void computeLPComponents(SCIP* scip, vector<int>& uf, vector<int>& rowComp)
{
    SCIP_ROW** rows = SCIPgetLPRows(scip);
    const int nrows = SCIPgetNLPRows(scip);
    const int ncols = SCIPgetNLPCols(scip);
    uf.resize(ncols);
    for (int j = 0; j < ncols; ++j)
        uf[j] = j;
    rowComp.assign(nrows, -1);
    for (int r = 0; r < nrows; ++r)
    {
        SCIP_COL** rcols = SCIProwGetCols(rows[r]);
        for (int k = 0; k < SCIProwGetNNonz(rows[r]); ++k)
        {
            const int j = SCIPcolGetLPPos(rcols[k]);
            if (j < 0 || !isUnfixedCol(rcols[k]))
                continue;
            if (rowComp[r] < 0)
                rowComp[r] = j;
            else
                uf[ufFind(uf, j)] = ufFind(uf, rowComp[r]);
        }
    }
}

/** Option -comp. The columns of the node LP that are unfixed are joined by the rows that contain them into connected
 *  components. A component is affected by the node if one of its rows contains a variable whose bounds changed at this
 *  node (branching or propagation) or is not in the parent LP; the other components are the same LP blocks as in the
 *  parent, whose duals are optimal for them. Their rows, and the rows with only fixed columns that do not touch a
 *  changed variable, get the parent duals, so that the difference with the parent is 0 on them instead of the noise of
 *  a degenerate LP. The mixed duals replace the collected ones only if they give the same LP bound. */
static bool mixDualsByComponents(SCIP* scip, const NodeInfo& parentInfo, SolCollector& collector)
{
    // Without the bounds of the parent LP (e.g. its duals come from a dive) the unchanged components are unknown
    if (parentInfo.DualSol.empty() || parentInfo.LPLb.empty() || SCIPinDive(scip) || SCIPinProbing(scip))
        return false;
    SCIP_ROW** rows = SCIPgetLPRows(scip);
    SCIP_COL** cols = SCIPgetLPCols(scip);
    const int nrows = SCIPgetNLPRows(scip);
    const int ncols = SCIPgetNLPCols(scip);

    // Components of the unfixed columns
    vector<int> uf, rowComp;
    computeLPComponents(scip, uf, rowComp);
    auto unfixed = [&](SCIP_COL* col){ return isUnfixedCol(col); };

    // Rows touching a variable whose bounds differ from the ones of the parent LP, and rows not in the parent LP.
    // Comparing with the bounds of the LP of the parent duals (rather than using the domchg of the focus node) also
    // catches the bound changes made after that LP: at the parent after its LP (reduced cost fixing, propagation),
    // or at the parent and earlier ancestors whose duals were copied because their bound did not change.
    vector<char> rowTouched(nrows, 0);
    auto touchRowsOf = [&](SCIP_COL* col){
        SCIP_ROW** crows = SCIPcolGetRows(col);
        for (int k = 0; k < SCIPcolGetNNonz(col); ++k)
        {
            const int r = SCIProwGetLPPos(crows[k]);
            if (r >= 0)
                rowTouched[r] = 1;
        }
    };
    for (int j = 0; j < ncols; ++j)
    {
        const size_t idx = allVarsIndex(SCIPcolGetVar(cols[j]));
        if (idx >= parentInfo.LPLb.size() || SCIPcolGetLb(cols[j]) != parentInfo.LPLb[idx] || SCIPcolGetUb(cols[j]) != parentInfo.LPUb[idx])
            touchRowsOf(cols[j]);
    }
    for (int r = 0; r < nrows; ++r)
        if (SCIProwGetIndex(rows[r]) >= (int)parentInfo.DualSol.size())
            rowTouched[r] = 1;
    vector<char> affected(ncols, 0);  // indexed by component root
    for (int r = 0; r < nrows; ++r)
        if (rowTouched[r] && rowComp[r] >= 0)
            affected[ufFind(uf, rowComp[r])] = 1;

    // Mixed duals: the parent dual on the rows of unaffected components and on untouched rows with only fixed columns
    vector<SCIP_Real> mixed = collector.dualsol;
    bool changed = false;
    for (int r = 0; r < nrows; ++r)
    {
        const bool keepChild = rowComp[r] >= 0 ? affected[ufFind(uf, rowComp[r])] : rowTouched[r];
        if (keepChild)
            continue;
        mixed[r] = parentInfo.DualSol[SCIProwGetIndex(rows[r])];
        changed = changed || mixed[r] != collector.dualsol[r];
    }
    if (!changed)
    {
        g_memoscip.DualMixUnchanged++;
        return false;
    }

    // Accept the mixed duals only if their LP bound is the one of the collected duals
    vector<SCIP_Real> mixedRedcost(collector.redcost.size());
    computeRedcosts(scip, mixed, mixedRedcost);
    const SCIP_Real childBound = lagrangianBound(scip, collector.dualsol, collector.redcost);
    const SCIP_Real mixedBound = lagrangianBound(scip, mixed, mixedRedcost);
    const SCIP_Real tol = 1e-7 * max(1.0, fabs(childBound));
    if (g_options.verbosity > 1)
    {
        int ncomp = 0, naffected = 0;
        for (int j = 0; j < ncols; ++j)
            if (unfixed(cols[j]) && !isSlackVar(SCIPcolGetVar(cols[j])) && ufFind(uf, j) == j)
            {
                ncomp++;
                naffected += affected[j];
            }
        cout<<"Dual components: "<<ncomp<<" components, "<<naffected<<" affected, child bound "<<childBound
            <<" mixed bound "<<mixedBound<<" LP obj "<<collector.objval<<endl;
    }
    if (SCIPisInfinity(scip, -mixedBound) || mixedBound < childBound - tol || fabs(childBound - collector.objval) > tol)
    {
        g_memoscip.DualMixRejected++;
        return false;
    }
    collector.dualsol = mixed;
    collector.redcost = mixedRedcost;
    g_memoscip.DualMixApplied++;
    return true;
}

/** dropNegligibleCoefs on a sparse proof part. */
static void dropNegligiblePart(SCIP* scip, ProofPart& part)
{
    const SCIP_Real tol = 1e-6 * max(1.0, fabs(part.rhs));
    size_t w = 0;
    for (size_t k = 0; k < part.inds.size(); ++k)
    {
        const int i = part.inds[k];
        const SCIP_Real c = part.coefs[k];
        bool drop = (c == 0.0);
        if (!drop && fabs(c) < tol)
        {
            SCIP_VAR* var = g_memoscip.AllVars[i];
            SCIP_Real bound = c > 0 ? SCIPvarGetUbGlobal(var) : SCIPvarGetLbGlobal(var);
            if (c > 0 && SCIPisInfinity(scip, bound) && isSlackVar(var))
                bound = getSlackImpliedUb(scip, i);
            if (!SCIPisInfinity(scip, fabs(bound)))
            {
                part.rhs -= c * bound;
                drop = true;
            }
        }
        if (!drop)
        {
            part.inds[w] = i;
            part.coefs[w] = c;
            ++w;
        }
    }
    part.inds.resize(w);
    part.coefs.resize(w);
}

/** Option -split: splits the node increment proof (rows weighted by the dual difference with the parent, plus the
 *  bound dual differences diffy) by connected component of the node LP, and stores the component of each unfixed
 *  column. A row goes to the component of its unfixed columns. A bound dual goes to the component of its column if
 *  it is unfixed; for a variable fixed to 1 at the node and at the parent, diffy = -sum_r a_r (y_r - y_r^P) is split
 *  over its rows so that it cancels the row terms of each component; the rest (variables newly fixed at the node,
 *  fully fixed rows) goes to an orphan part. The parts sum to the increment (checked with -d). */
static void buildIncrementParts(SCIP* scip, NodeInfo& node, const NodeInfo& parentInfo, const SolCollector& collector,
    const vector<SCIP_Real>& diffy, const vector<SCIP_Real>& total, SCIP_Real totalRhs)
{
    SCIP_ROW** rows = SCIPgetLPRows(scip);
    SCIP_COL** cols = SCIPgetLPCols(scip);
    const int nrows = SCIPgetNLPRows(scip);
    const int ncols = SCIPgetNLPCols(scip);
    vector<int> uf, rowComp;
    computeLPComponents(scip, uf, rowComp);
    auto compOfCol = [&](int j){ return allVarsIndex(SCIPcolGetVar(cols[ufFind(uf, j)])); };

    node.CompLabel.clear();
    for (int j = 0; j < ncols; ++j)
        if (isUnfixedCol(cols[j]))
            node.CompLabel.push_back({allVarsIndex(SCIPcolGetVar(cols[j])), compOfCol(j)});
    sort(node.CompLabel.begin(), node.CompLabel.end());

    // group -1: orphan
    std::map<int, std::map<int, SCIP_Real>> coef;
    std::map<int, SCIP_Real> rhs;
    auto rowGroup = [&](int r){ return rowComp[r] >= 0 ? compOfCol(rowComp[r]) : -1; };

    vector<SCIP_Real> rowDelta(nrows, 0.0);
    for (int r = 0; r < nrows; ++r)
    {
        if (SCIProwGetNNonz(rows[r]) <= 1)
            continue;
        const int idx = SCIProwGetIndex(rows[r]);
        SCIP_Real dual = collector.dualsol[r];
        if (idx < (int)parentInfo.DualSol.size())
            dual -= parentInfo.DualSol[idx];
        rowDelta[r] = dual;
        const SCIP_Real rowValue = equalityRowValue(scip, rows[r]);
        if (dual == 0.0 || rowValue == SCIPinfinity(scip))
            continue;
        const int g = rowGroup(r);
        SCIP_COL** rcols = SCIProwGetCols(rows[r]);
        SCIP_Real* rvals = SCIProwGetVals(rows[r]);
        for (int k = 0; k < SCIProwGetNNonz(rows[r]); ++k)
            coef[g][allVarsIndex(SCIPcolGetVar(rcols[k]))] += dual * rvals[k];
        rhs[g] += dual * rowValue;
    }

    auto addBound = [&](int g, int i, SCIP_Real t){
        coef[g][i] += t;
        forEachAuxVar(shiftedIndex(g_memoscip.AllVars[i]), [&](int aux){ coef[g][aux] += t; });
        rhs[g] += t;
    };
    for (size_t i = 0; i < diffy.size(); ++i)
    {
        const SCIP_Real d = diffy[i];
        if (g_memoscip.VarUb[i] == SCIPinfinity(scip))
            continue;
        SCIP_VAR* var = g_memoscip.AllVars[i];
        SCIP_COL* col = SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN ? SCIPvarGetCol(var) : nullptr;
        const int j = col != nullptr ? SCIPcolGetLPPos(col) : -1;
        // Even with d = 0 the row pieces of a variable fixed to 1 at both nodes must be cancelled per component
        const bool fixedOneAtBoth = j >= 0 && g_memoscip.VarLb[i] == 1 && i < parentInfo.LPLb.size() && parentInfo.LPLb[i] == 1;
        if (d == 0.0 && !fixedOneAtBoth)
            continue;
        if (j >= 0 && isUnfixedCol(cols[j]))
        {
            addBound(compOfCol(j), (int)i, d);
            continue;
        }
        SCIP_Real rest = d;
        if (fixedOneAtBoth)
        {
            SCIP_ROW** crows = SCIPcolGetRows(col);
            SCIP_Real* cvals = SCIPcolGetVals(col);
            for (int k = 0; k < SCIPcolGetNNonz(col); ++k)
            {
                const int r = SCIProwGetLPPos(crows[k]);
                if (r < 0 || SCIProwGetNNonz(crows[k]) <= 1 || rowDelta[r] == 0.0)
                    continue;
                const SCIP_Real t = -cvals[k] * rowDelta[r];
                addBound(rowGroup(r), (int)i, t);
                rest -= t;
            }
        }
        if (rest != 0.0)
            addBound(-1, (int)i, rest);
    }

    node.IncParts.clear();
    vector<SCIP_Real> check;
    SCIP_Real checkRhs = 0.0;
    if (g_options.debug)
        check.assign(total.size(), 0.0);
    for (auto& [g, m] : coef)
    {
        ProofPart part;
        part.orphan = (g == -1);
        part.rhs = rhs[g];
        for (auto& [i, v] : m)
        {
            if (g_options.debug)
                check[i] += v;
            if (v != 0.0)
            {
                part.inds.push_back(i);
                part.coefs.push_back(v);
            }
        }
        if (g_options.debug)
            checkRhs += part.rhs;
        dropNegligiblePart(scip, part);
        if (!part.inds.empty() || part.rhs != 0.0)
            node.IncParts.push_back(std::move(part));
    }
    for (auto& [g, v] : rhs)
        if (!coef.count(g) && g_options.debug)
            checkRhs += v;
    if (g_options.debug)
    {
        for (size_t i = 0; i < total.size(); ++i)
            assert(fabs(check[i] - total[i]) <= 1e-6 * max(1.0, fabs(total[i])));
        assert(fabs(checkRhs - totalRhs) <= 1e-6 * max(1.0, fabs(totalRhs)));
    }
    (void)totalRhs;
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
    // With mixed duals the reduced costs are no longer the ones of the LP solver
    const bool mixedDuals = g_options.dualComponents && nodeNumber>1 && mixDualsByComponents(scip, parentInfo, collector);
    (void)mixedDuals;

    SCIP_COL** cols = SCIPgetLPCols(scip);
    SCIP_ROW** rows = SCIPgetLPRows(scip);
    SCIP_LPSOLSTAT stat = SCIPgetLPSolstat(scip);
    assert(stat==SCIP_LPSOLSTAT_OPTIMAL || stat==SCIP_LPSOLSTAT_OBJLIMIT);
    assert((int)collector.dualsol.size()==nrows);
    assert((int)collector.redcost.size()==ncols);
    (void)stat;
    if(g_options.debug){
        for(int i =0; i < ncols; i++){
            assert(SCIPinDive(scip) || mixedDuals || SCIPisFeasEQ(scip, SCIPgetVarRedcost(scip, SCIPcolGetVar(cols[i])), collector.redcost[i]));
            if(!isSlackVar(g_memoscip.AllVars[i])){
                const int pos = shiftedIndex(g_memoscip.AllVars[i]);
                assert(g_options.encoding!=0 || g_memoscip.VarLb[i] == SCIPvarGetLbLocal(g_memoscip.AllVars[i]) || anyAuxVar(pos, [](int aux){ return g_memoscip.VarUb[aux]==0; }) );
                assert(g_options.encoding!=0 || g_memoscip.VarUb[i] == SCIPvarGetUbLocal(g_memoscip.AllVars[i]) || anyAuxVar(pos, [](int aux){ return g_memoscip.VarLb[aux]==1; }) );
                (void)pos;
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
        if(redc != 0.0){
            if(g_memoscip.VarUb[idx]==0)
                newy_i[idx]=0;
            else if(g_memoscip.VarLb[idx]==1)
                newy_i[idx]=redc;
            else if(g_memoscip.VarUb[idx]==SCIPinfinity(scip)){
                // No upper bound constraint, hence no bound dual: a negative reduced cost can only be dual infeasibility noise
                assert(redc > -1e-6);
                newy_i[idx]=0;
            }
            else
                newy_i[idx]=min(redc,0.0);
            assert(!isSlackVar(g_memoscip.AllVars[idx]) || newy_i[idx]==0);
        }
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
            const SCIP_Real rowValue = equalityRowValue(scip, rows[i]);
            if(dual != 0.0 && rowValue != SCIPinfinity(scip)){
                sumdual+=dual*rowValue;
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
                    assert(!anyAuxVar(shiftedIndex(g_memoscip.AllVars[i]), [](int aux){ return g_memoscip.VarLb[aux]!=0; }) || SCIPisDualfeasZero(scip,diffy_i[i]) || g_options.encoding==1);
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

        if(diffy_i[i] != 0.0 && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
            TempCoef[i]+=diffy_i[i];
            forEachAuxVar(shiftedIndex(var), [&](int aux){ TempCoef[aux]+=diffy_i[i]; });
            sumdual+=diffy_i[i];
        }
        if(g_options.debug){
            if(!SCIPisDualfeasZero(scip,newy_i[i]) && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
                Truesumdual+=newy_i[i];
                TrueTempCoef[i]+=newy_i[i];
                forEachAuxVar(shiftedIndex(var), [&](int aux){ TrueTempCoef[aux]+=newy_i[i]; });
            }
        }
    }

    assert(nAllVars == TempCoef.size());
    node.IncParts.clear();
    if(g_options.splitProofs){
        node.CompLabel.clear();
        if(!SCIPinDive(scip) && !SCIPinProbing(scip))
            buildIncrementParts(scip, node, parentInfo, collector, diffy_i, TempCoef, sumdual);
    }
    if(g_options.verbosity > 2){
        cout<<"DUAL proof"<<endl;
        printProof(scip, TempCoef, ">=", sumdual, true);
    }

    if(g_options.verbosity > 2 && g_options.debug){
        cout<<"TRUEDUAL proof"<<endl;
        printProof(scip, TrueTempCoef, ">=", Truesumdual, true);
    }
    // The reduced costs use the exact coefficients, the stored proof the safely cleaned ones
    vector<SCIP_Real> ProofCoef=TempCoef;
    dropNegligibleCoefs(scip, ProofCoef, sumdual);
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
        if(ProofCoef[i] == 0.0)
            continue;
        coef = ProofCoef[i];
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
                // Slack s >= 0: its new redcost -coef may be dual infeasible up to the LP dual feasibility tolerance
                assert(SCIPisDualfeasLE(scip, coef, 0.0));
            }
        }
    }
    node.DualSol=IndexedDualSol;
    node.LPLb.clear();
    node.LPUb.clear();
    if(g_options.dualComponents && !SCIPinDive(scip) && !SCIPinProbing(scip)){
        node.LPLb.assign(nAllVars, 0.0);
        node.LPUb.assign(nAllVars, -1.0);  // lb > ub: never equal for a variable not in the LP
        for(int i = 0; i < ncols; i++){
            const int idx = allVarsIndex(SCIPcolGetVar(cols[i]));
            node.LPLb[idx] = SCIPcolGetLb(cols[i]);
            node.LPUb[idx] = SCIPcolGetUb(cols[i]);
        }
    }
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
    assert(fabs(collector.objval-sumdual-parentInfo.lb) <= 1e-6*max(1.0,fabs(collector.objval)) || Myceil(collector.objval)>= MyPrimalbound(scip));

    SCIPaggrRowFree(scip, &proofrow);
    if(Trueproofrow != nullptr)
        SCIPaggrRowFree(scip, &Trueproofrow);

    // With option -d: check that the reduced costs are correct, i.e. that the LP with objective redcost has value 0
    /*if(g_options.debug && !SCIPinDive(scip))
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
    }*/

    return SCIP_OKAY;
}

/** Marks the node as infeasible: its lower bound and proof right-hand side become TOP. */
static void markNodeInfeasible(NodeInfo& node)
{
    node.IncParts.clear();
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
                else if(isDiveRelaxedRow(rows[i]))  // equality whose slack column is not in the LP: use its finite side
                    SCIP_CALL( SCIPaggrRowAddRow(scip, proofrow, rows[i], rowFarkas[i], SCIPisInfinity(scip, SCIProwGetRhs(rows[i])) ? -1 : 1) );
                else
                    assert(SCIProwGetName(rows[i])[0]=='s' && SCIProwGetName(rows[i])[1]=='c');
            }
        }
        vector<SCIP_Real> Coefs(nAllVars,0.0);
        // The aggregated row a.x = b is infeasible under the current bounds: either maxact(a.x) < b, and the proof is
        // a.x >= b, or minact(a.x) > b, and the proof is -a.x >= -b.
        const SCIP_Real proofrhs = SCIPaggrRowGetRhs(proofrow);
        SCIP_Real maxact = 0.0, minact = 0.0;
        for (size_t i = 0; i < nAllVars; ++i)
        {
            SCIP_Real coef = SCIPaggrRowGetProbvarValue(proofrow, SCIPvarGetProbindex(g_memoscip.AllVars[i]));
            if(SCIPisDualfeasZero(scip,coef))
                continue;
            maxact += coef * (coef > 0 ? g_memoscip.VarUb[i] : g_memoscip.VarLb[i]);
            minact += coef * (coef > 0 ? g_memoscip.VarLb[i] : g_memoscip.VarUb[i]);
        }
        assert(maxact < proofrhs || minact > proofrhs);
        const bool negrhs = !(maxact < proofrhs);
        for (size_t i = 0; i < nAllVars; ++i)
        {
            SCIP_VAR* var = g_memoscip.AllVars[i];
            SCIP_Real coef = SCIPaggrRowGetProbvarValue(proofrow, SCIPvarGetProbindex(var));

            if(negrhs)
                coef = -coef;

            if(!SCIPisDualfeasZero(scip,coef)){
                //Only keep variables assigned to 0 in the proof.
                //A variable assigned to 1 is replaced by its auxiliary variable(s), which are assigned to 0.
                if(g_memoscip.VarUb[i]==0 && coef > 0){
                    Coefs[i]+=coef;
                }else if(g_memoscip.VarLb[i]==1 && coef < 0){
                    forEachAuxVar(shiftedIndex(var), [&](int aux){ 
                        if(g_memoscip.VarUb[aux]==0)
                            Coefs[aux]-=coef;
                        else{
                            assert(g_memoscip.VarLb[aux]==1);
                            //find index of aux in BPRemoval and get the explanation
                            int auxIdx = std::distance(node.BPRemoval.begin(), std::find_if(node.BPRemoval.begin(), node.BPRemoval.end(),[aux](const std::pair<int, vector<int>>& removal){ return removal.first == aux; }));
                            assert(auxIdx < (int)node.BPRemoval.size());
                            const vector<int>& explanation = node.BPRemoval[auxIdx].second;
                            for(int idxj : explanation) 
                                Coefs[idxj]-=coef;
                        }
                    });
                }
            }
        }
        /*char nameprobsol[50];
        snprintf(nameprobsol, 50, "inf.lp");
        SCIP_CALL( SCIPwriteLP(scip, nameprobsol) );
        snprintf(nameprobsol, 50, "inflp.lp");
        SCIP_CALL( SCIPwriteTransProblem(scip, nameprobsol,"lp",FALSE) );*/

        for(size_t i = 0; i < nAllVars; ++i )
        {
            if(SCIPisDualfeasZero(scip,Coefs[i]))
                continue;
            int idx=i;
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
    SCIP_Longint nodeNumber, SCIP_Longint parentNumber, const vector<int>& fixedTo1)
{
    ScopedTimer GetDualProofFromPseudoObjCutoffTimer(g_memoscip.GetDualProofFromPseudoObjCutoffTime);
    if(g_options.verbosity > 3){
        cout << "Generating dual proof for cut off node " << nodeNumber << " with parent " << parentNumber << endl;
    }
    const size_t nAllVars = g_memoscip.AllVars.size();
    NodeInfo& parentInfo = g_memoscip.nodeInfoMap[parentNumber];
    vector<SCIP_Real> y_i(nAllVars,0.0);
    vector<SCIP_Real> diffy(nAllVars,0.0);

    // The pseudo objective of the node is the objective of the variables fixed to 1
    SCIP_Real obj=0;
    for (int idx : fixedTo1){
        y_i[idx] = SCIPvarGetObj(g_memoscip.AllVars[idx]);
        obj += y_i[idx];
    }
    for(size_t i = 0; i < nAllVars; ++i){
        if(i<parentInfo.y_i.size()){
            diffy[i]=y_i[i]-parentInfo.y_i[i];
        }
    }
    assert(Myceil(obj) >= MyPrimalbound(scip));

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
            const SCIP_Real rowValue = equalityRowValue(scip, rows[i]);
            if(dual != 0.0 && rowValue != SCIPinfinity(scip)){
                sumdual+=dual*rowValue;
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

        if(diffy[i] != 0.0 && g_memoscip.VarUb[i]!=SCIPinfinity(scip) ){
            TempCoef[i]+=diffy[i];
            forEachAuxVar(shiftedIndex(var), [&](int aux){ TempCoef[aux]+=diffy[i]; });
            sumdual+=diffy[i];
        }
    }
    SCIPaggrRowFree(scip, &proofrow);
    assert(fabs(obj-sumdual-parentInfo.lb) <= 1e-6*max(1.0,fabs(obj)));
    dropNegligibleCoefs(scip, TempCoef, sumdual);

    NodeInfo& node = g_memoscip.nodeInfoMap[nodeNumber];
    for(size_t i = 0; i < nAllVars; ++i )
    {
        SCIP_Real coef = TempCoef[i];
        if(coef == 0.0)
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
