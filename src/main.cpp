#include "scip/scipdefplugins.h"

#include "memoscip_dual.h"
#include "memoscip_events.h"
#include "memoscip_option.h"
#include "memoscip_reader.h"
#include "memoscip_types.h"
#include "sepa_mycuts.h"

#include "scip/scip_cons.h"
#include "scip/scip_prob.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>


using namespace std;

/** Sets the SCIP parameters: pure branch-and-bound with DFS, no heuristics, presolving, conflict analysis,
 *  native separators nor propagators (except the ones enabled by the options). */
static SCIP_RETCODE configureScip(SCIP* scip)
{
    SCIP_CALL(SCIPsetIntParam(scip, "lp/refactorinterval", 200));
    //Do DFS node selection
    SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/dfs/stdpriority", 1000000) );
    //Pscost branching (disable strong branching), or relpscost with -b=3
    if(g_options.branchingHeuristic==3){
        SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", -2000000) );
        SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/priority", 2000000) );
    }
    else
        SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", 2000000) );

    //FULL VERBOSITY
    if(g_options.verbosity > 5){
        SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 2) );
        SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) );
    }else{
        SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
        SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", FALSE) );
    }

    // disable heuristics, presolving
    SCIP_CALL(SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
    SCIP_CALL(SCIPsetBoolParam(scip,"symmetries/enabled",FALSE));


    //Disable all separators for the moment
    SCIP_CALL( SCIPsetIntParam(scip, "constraints/knapsack/sepacardfreq", -1) );
    SCIP_CALL( SCIPsetIntParam(scip, "constraints/knapsack/sepafreq", -1) );
    SCIP_CALL( SCIPsetIntParam(scip, "constraints/knapsack/maxsepacuts", 0) );
    SCIP_CALL( SCIPsetIntParam(scip, "constraints/knapsack/maxsepacutsroot", 0) );
    SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/maxsepacutsroot", 0) );
    SCIP_CALL( SCIPsetIntParam(scip, "separating/mcf/freq", -1) );
    SCIP_CALL( SCIPsetBoolParam(scip, "separating/mcf/separateknapsack", FALSE) );
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/upgrade/knapsack", FALSE));

    // Disable EVERY native SCIP separator individually (freq = -1),
    // without using SCIPsetSeparating, so that maxrounds is left untouched
    int nsepas = SCIPgetNSepas(scip);
    SCIP_SEPA** sepas = SCIPgetSepas(scip);
    for (int i = 0; i < nsepas; ++i) {
        SCIPsepaSetFreq(sepas[i], -1);   // -1 = never (no SCIP_CALL, the function returns void)
    }


    //TO VERIFY
    SCIP_CALL( SCIPsetIntParam(scip, "lp/disablecutoff", 1) );  // no objective limit: LPs are solved to optimality, so proofs certify the full bound
    SCIP_CALL( SCIPsetBoolParam(scip, "lp/alwaysgetduals", TRUE) );
    // No LP presolving: re-solves only need a few iterations and presolving the (dense) learned rows dominated the LP time
    SCIP_CALL( SCIPsetBoolParam(scip, "lp/presolving", FALSE) );

    SCIP_CALL( SCIPsetIntParam(scip, "pricing/maxvars", 10000) );
    SCIP_CALL( SCIPsetIntParam(scip, "pricing/maxvarsroot", 20000) );

    if(g_options.sepa){
        cout<<"Enabling custom separation"<<endl;
        //Enable separators
        SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", -1) );        // unlimited per node (default behaviour)
        SCIP_CALL( SCIPsetIntParam(scip, "separating/maxroundsroot", -1) );    // unlimited at the root
        SCIP_CALL( SCIPsetIntParam(scip, "separating/poolfreq", 1) );    // at every node
        includeMycutsSeparator(scip);
        SCIP_CALL(SCIPsetBoolParam(scip, "lp/cleanuprows", FALSE)); //Disable Row removal (TO BE TESTED)
        SCIP_CALL(SCIPsetBoolParam(scip, "lp/cleanuprowsroot", FALSE)); //Disable Row removal (TO BE TESTED)
        SCIP_CALL(SCIPsetIntParam(scip, "lp/rowagelimit", -1)); //Disable row aging : TO TEST DIFFERENT VALUES
        SCIP_CALL(SCIPsetIntParam(scip, "separating/maxcuts", g_options.nbslack));      // outside the root (default 100)
        SCIP_CALL(SCIPsetIntParam(scip, "separating/maxcutsroot", g_options.nbslack));
        SCIP_CALL(SCIPsetRealParam(scip, "cutselection/hybrid/minortho", 0.0));
        SCIP_CALL(SCIPsetRealParam(scip, "cutselection/hybrid/minorthoroot", 0.0));
    }else{
        SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    }

    //Multi threading disabled
    SCIP_CALL( SCIPsetIntParam(scip, "lp/threads", 1) );
    SCIP_CALL( SCIPsetIntParam(scip, "parallel/maxnthreads", 1) );

    if(g_options.BC){
        cout<<"Bound propagation enabled"<<endl;
        //Enable basic propagation
        SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", -1) );
        SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", -1) );
    }else{
        cout<<"Bound propagation disabled"<<endl;
        SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", 0) );
        SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );
    }

    if(g_options.NC)
        cout<<"Node consistency enabled"<<endl;
    else
        cout<<"Node consistency disabled"<<endl;

    if(g_options.INF)
        cout<<"Infeasibility certificate collection enabled"<<endl;
    else
        cout<<"Infeasibility certificate collection disabled"<<endl;

    //disable other propagation (reduced cost filtering is kept for node consistency)
    vector<const char*> disabledProps = {"dualfix", "genvbounds", "nlobbt", "obbt", "probing", "pseudoobj", "rootredcost"};
    if(!g_options.NC)
        disabledProps.push_back("redcost");
    disabledProps.push_back("vbounds");
    for(const char* prop : disabledProps){
        SCIP_CALL(SCIPsetIntParam(scip, ("propagating/" + string(prop) + "/freq").c_str(), -1));
        SCIP_CALL(SCIPsetIntParam(scip, ("propagating/" + string(prop) + "/maxprerounds").c_str(), 0));
    }
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/rangedrowartcons", FALSE));


    // Collect Farkas proofs
    SCIP_CALL( SCIPsetBoolParam(scip, "lp/checkfarkas", TRUE) );
    //Disable conflict
    SCIP_CALL(SCIPsetBoolParam(scip, "conflict/enable", FALSE));


    //Fix numerics
    SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol",     1e-6) );
    SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-7) );
    SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon",     1e-9) );

    SCIP_CALL( SCIPsetBoolParam(scip, "misc/scaleobj", FALSE) );

    return SCIP_OKAY;
}

/** Turns every inequality constraint into an equality by adding a nonnegative slack variable slack_<i>.
 *  Returns the number of converted constraints in converted. */
static SCIP_RETCODE addSlackVariables(SCIP* scip, int& converted)
{
    int nconss = SCIPgetNConss(scip);
    SCIP_CONS** conss = SCIPgetConss(scip);
    converted = 0;

    for (int i = 0; i < nconss; ++i)
    {
        SCIP_CONS* cons = conss[i];
        if (cons == nullptr) continue;

        SCIP_Real lhs = SCIPgetLhsLinear(scip,cons);
        SCIP_Real rhs = SCIPgetRhsLinear(scip,cons);
        if (SCIPisEQ(scip, lhs, rhs))
            continue;

        char sname[256];
        snprintf(sname, sizeof(sname), "slack_%d", i);

        SCIP_VAR* slack = nullptr;
        SCIP_CALL( SCIPcreateVarBasic(scip, &slack, sname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
        SCIP_CALL( SCIPaddVar(scip, slack) );
        if (!SCIPisInfinity(scip, rhs))
        {
            // a.x <= rhs  ->  a.x + slack = rhs
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, slack, 1.0) );
            SCIP_CALL( SCIPchgLhsLinear(scip, cons, rhs) );
            SCIP_CALL( SCIPchgRhsLinear(scip, cons, rhs) );
            converted++;
        }
        else if (!SCIPisInfinity(scip, lhs)){
            // a.x >= lhs  ->  a.x - slack = lhs
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, slack, -1.0) );
            SCIP_CALL( SCIPchgLhsLinear(scip, cons, lhs) );
            SCIP_CALL( SCIPchgRhsLinear(scip, cons, lhs) );
            converted++;
        }
        else
        {
            SCIPdebugMessage("Skipping ranged/general linear constraint (index %d)\n", i);
        }
        SCIP_CALL(SCIPsetConsRemovable(scip,cons,FALSE));
        SCIP_CALL( SCIPreleaseVar(scip, &slack) );
    }

    std::cout << "Converted " << converted << " linear constraints to equalities by adding slack vars\n";
    return SCIP_OKAY;
}

/** For .opb instances: adds a binary variable <x>_0 and the constraint x + <x>_0 = 1 for the first nvars variables. */
static SCIP_RETCODE addOpbDomainConstraints(SCIP* scip, int nvars)
{
    for (int i = 0; i < nvars; ++i)
    {
        SCIP_VAR* v = SCIPgetVars(scip)[i];  // the variable array changes as variables are added
        char consname[256];
        snprintf(consname, sizeof(consname), "ub1_%d_%s", i, SCIPvarGetName(v));
        char sname[256];
        snprintf(sname, sizeof(sname), "%s_0", SCIPvarGetName(v));
        SCIP_VAR* slack = nullptr;
        cout<<"Creating slack variable "<<sname<<" for upper bound constraint of variable "<<SCIPvarGetName(v)<<endl;
        SCIP_CALL( SCIPcreateVarBasic(scip, &slack, sname, 0.0, 1, 0.0, SCIP_VARTYPE_BINARY) );
        SCIP_CALL( SCIPaddVar(scip, slack) );

        SCIP_VAR* consvars[2] = { v, slack };
        SCIP_Real coefs[2] = { 1.0, 1.0 };

        SCIP_CONS* ubcons = nullptr;

        SCIP_CALL(SCIPcreateConsLinear(scip, &ubcons, consname,
            2, consvars, coefs,
            1.0, 1.0,   /* lhs, rhs: equality 1 <= ... <= 1 */
            TRUE,   /* initial */
            FALSE,  /* separate */
            TRUE,   /* enforce */
            TRUE,   /* check */
            TRUE,   /* propagate */
            FALSE,  /* local */
            FALSE,  /* modifiable */
            FALSE,  /* dynamic */
            FALSE,  /* removable: not removed by LP aging */
            FALSE   /* stickingatnode */
        ));

        SCIP_CALL( SCIPaddCons(scip, ubcons) );
        SCIP_CALL( SCIPreleaseCons(scip, &ubcons) );
        SCIP_CALL( SCIPreleaseVar(scip, &slack) );
    }
    std::cout << "Added " << nvars << " upper-bound-1 constraints (var <= 1)." << std::endl;
    return SCIP_OKAY;
}

/** Variables x_{i0} / x_{i1} are named "*_0" / "*_1". The preferred side (-bo) is put in OrderVar for the custom
 *  branching rule; with SCIP branching (-b>=1) the other side gets a low branching priority. */
static SCIP_RETCODE setBranchingOrder(SCIP* scip, vector<SCIP_Real>& objcoefs)
{
    SCIP_VAR** vars = SCIPgetVars(scip);
    int nvars = SCIPgetNVars(scip);
    for (int i = 0; i < nvars; ++i)
    {
        SCIP_VAR* v=vars[i];
        objcoefs.push_back(SCIPvarGetObj(v));
        string varname = SCIPvarGetName(v);
        const bool ends0 = varname.ends_with("_0");
        const bool ends1 = varname.ends_with("_1");
        const bool preferred = (ends0 && g_options.branchingOrder==0) || (ends1 && g_options.branchingOrder==1);
        const bool other = (ends0 && g_options.branchingOrder==1) || (ends1 && g_options.branchingOrder==0);
        if (g_options.branchingHeuristic >= 1 && other)
        {
            SCIP_CALL(SCIPchgVarBranchPriority(scip, v, -1000000));
        }
        if (preferred) {
           g_memoscip.OrderVar.push_back(v);
        }
    }
    return SCIP_OKAY;
}

static void printStatistics(SCIP* scip)
{
    cout<<"Solving Status : "<<SCIPgetStatus(scip)<<endl;
    cout<<"Solving time : "<<SCIPgetSolvingTime(scip)<<endl;
    cout<<"Nodes : "<<SCIPgetNNodes(scip)<<endl;
    cout<<"Primal bound : "<<SCIPgetPrimalbound(scip)<<endl;
    cout<<"Dual bound : "<<SCIPgetDualbound(scip)<<endl;
    cout<<"Number of learned cuts : "<<g_memoscip.nbAddedctr<<endl;
    cout<<"Gap : "<<SCIPgetGap(scip)<<endl;
    cout<<"Time in AnalyseDualSol : "<<g_memoscip.AnalyseDualSolTime.count()<<" seconds"<<endl;
    cout<<"Time in eventExecNodeDelete : "<<g_memoscip.NodeDeleteTime.count()<<" seconds"<<endl;
    cout<<"Time in eventExecNodeSolved : "<<g_memoscip.NodeSolvedTime.count()<<" seconds"<<endl;
    cout<<"Time in eventExecLPSolved : "<<g_memoscip.LPSolvedTime.count()<<" seconds"<<endl;
    cout<<"Time in eventExecNodeFocused : "<<g_memoscip.NodeFocusedTime.count()<<" seconds"<<endl;
    cout<<"Time in eventExecBoundChanged : "<<g_memoscip.BoundChangedTime.count()<<" seconds"<<endl;
    cout<<"Time in gatherLPInfeasibilityCertificate : "<<g_memoscip.GatherLPInfeasibilityCertificateTime.count()<<" seconds"<<endl;
    cout<<"Time in getDualProofFromPseudoObjCutoff : "<<g_memoscip.GetDualProofFromPseudoObjCutoffTime.count()<<" seconds"<<endl;
}

/** With -v>2: prints the best solution and writes it to solution.sol. */
static SCIP_RETCODE printSolution(SCIP* scip)
{
    SCIP_SOL* sol = SCIPgetBestSol(scip);
    if (sol != nullptr)
    {
        std::cout << "Solution found with objective: " << SCIPgetSolOrigObj(scip, sol) << std::endl;
        cout << "Solving time : " << SCIPgetSolvingTime(scip) << endl;
        cout << "Nodes : " << SCIPgetNNodes(scip) << endl;

        SCIP_CALL(SCIPprintSol(scip, sol, NULL, FALSE));

        FILE* file = fopen("solution.sol", "w");
        if (file != nullptr)
        {
            SCIP_CALL(SCIPprintSol(scip, sol, file, FALSE));
            fclose(file);
        }
    }
    else
    {
        std::cout << "No solution found." << std::endl;
    }
    cout<<"Should have added " << g_memoscip.nbAddedctr << " constraints." <<endl;
    if(!g_memoscip.CutQueue.queue.empty())
        cout<<"Still in queue "<<g_memoscip.CutQueue.queue.size()<<endl;
    return SCIP_OKAY;
}

int main(int argc, char** argv)
{
    std::string parseErr;
    if (!parse_memoscip_cli(argc, argv, g_options, parseErr))
    {
        std::cerr << "Error: " << parseErr << "\n";
        print_memoscip_usage(argv[0] ? argv[0] : "memobound");
        return 1;
    }

    string filepath = g_options.instancePath;
    string fileext = filesystem::path(filepath).extension().string();


    //If extension is .wcsp then we need to convert it to .lp using wcsp2lp-support_learn.py
    if(fileext==".wcsp"){
        cout<<"Converting "<<filepath<<" to .lp using wcsp2lp-support_learn.py"<<endl;
        string lpfile = filepath.substr(0, filepath.size()-5)+".lp";
        string command = "python3 /home/kisouke/MemoSCIP/wcsp2lp-support_learn.py "+filepath+" "+lpfile;
        system(command.c_str());
        filepath = lpfile;
        fileext = ".lp";
    }
    g_memoscip.FileExt=fileext;
    SCIP* scip = nullptr;
    SCIP_CALL(SCIPcreate(&scip));
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL(SCIPcreateProbBasic(scip, "pure_bab"));
    // read problem filepath
    if (g_options.kpcg) {
        std::string rerr;
        if (!readKPCGInstance(scip, filepath, g_options.encoding, rerr)) {
            std::cerr << "Error reading KPCG instance: " << rerr << std::endl;
            return 1;
        }
    } else {
        SCIP_CALL(SCIPreadProb(scip, filepath.c_str(), NULL));
    }

    if (g_options.verbosity < 0)
        g_options.verbosity = 0;

    SCIP_CALL( configureScip(scip) );

    includeLPSolveEventHandler(scip);
    register_memoscip_branchrule(scip);

    int nvarsBeforeSlacks = SCIPgetNVars(scip);
    int converted = 0;
    SCIP_CALL( addSlackVariables(scip, converted) );

    if(fileext==".opb"){
        //If problem is in .opb format, add x_{i0}+ x_{i1} = 1 for all variables
        SCIP_CALL( addOpbDomainConstraints(scip, nvarsBeforeSlacks-converted) );
    }
    g_memoscip.nOriginalVar=SCIPgetNVars(scip);

    // Learning event handlers when -l=N (N>=1)
    if (g_options.learning > 0)
    {
        includeNodeFocusedEventHandler(scip);
        includeNodeSolvedEventHandler(scip);
        includeNodeDeleteEventHandler(scip);
        includeNodeInfeasibleEventHandler(scip);
        includeBoundChangedEventHandler(scip);
    }

    SCIP_CALL(SCIPtransformProb(scip));
    if (g_options.output > 0)
        SCIP_CALL(SCIPwriteTransProblem(scip, "begin.lp", "lp", FALSE));

    // The objective coefficients are the reduced costs of the virtual node 0 (parent of the root)
    vector<SCIP_Real> objcoefs;
    SCIP_CALL( setBranchingOrder(scip, objcoefs) );
    g_memoscip.nodeInfoMap[0].redcosts=objcoefs;

    if (g_options.profiling)
        std::cerr << "[memobound] profiling flag (-p) set (reserved for future timing / SCIP hooks)\n";
    if (g_options.verbosity > 1 || g_options.output)
        cout << "Start Solving " << endl;
    if (g_options.timelimit > 0)
    {
        SCIP_CALL( SCIPsetRealParam(scip, "limits/time", (SCIP_Real)g_options.timelimit) );
        if (g_options.verbosity > 1)
            cout << "Time limit set to " << g_options.timelimit << " seconds" << endl;
    }

    if (g_options.hasInitialUB)
        SCIP_CALL( SCIPsetObjlimit(scip, g_options.initialUB) );

    SCIP_CALL(SCIPsolve(scip));
    if(g_options.finish_Learn)
        FinishLearn(scip, 1);

    if (g_options.output > 0){
        SCIP_CALL(SCIPwriteTransProblem(scip, "end.lp", "lp", FALSE));
    }

    printStatistics(scip);
    if(g_options.verbosity > 2)
        SCIP_CALL( printSolution(scip) );
    //Print SCIP statistics
    //SCIP_CALL( SCIPprintStatistics(scip, NULL) );
    // Cleanup: clear per-node bookkeeping map
    g_memoscip.nodeInfoMap.clear();
    SCIP_CALL(SCIPfree(&scip));
    return 0;
}
