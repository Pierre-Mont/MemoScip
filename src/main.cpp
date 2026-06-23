#include "scip/scipdefplugins.h"

#include "memoscip_dual.h"
#include "memoscip_events.h"
#include "memoscip_option.h"
#include "memoscip_reader.h"
#include "memoscip_types.h"
#include "sepa_mycuts.h"

#include "scip/scip_cons.h"
#include "scip/scip_prob.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

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
        if (!readKPCGInstance(scip, filepath, rerr)) {
            std::cerr << "Error reading KPCG instance: " << rerr << std::endl;
            return 1;
        }
    } else {
        SCIP_CALL(SCIPreadProb(scip, filepath.c_str(), NULL));
    }
    
    if (g_options.verbosity < 0)
        g_options.verbosity = 0;
    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", g_options.verbosity));


    //Do DFS node selection
    SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/dfs/stdpriority", 1000000) );
    //Pscost branching (disable strong branching)
    if(g_options.branchingHeuristic==3){
        SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", -2000000) );
        SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/priority", 2000000) );
        }    // register handlers
    else
        SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", 2000000) );    // register handlers

    //FULL VERBOSITY
     if(g_options.verbosity > 3){
        SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 2) );
        SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) );
    }

    // disable heuristics, presolving, separation
    SCIP_CALL(SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
    //SCIP_CALL(SCIPsetIntParam(scip,"misc/usesymmetry",0));
    SCIP_CALL(SCIPsetBoolParam(scip,"symmetries/enabled",FALSE));

    SCIP_CALL( SCIPsetIntParam(scip, "lp/disablecutoff", 1) );

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

     //disable other propagation
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/dualfix/freq", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/dualfix/maxprerounds", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/genvbounds/freq", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/genvbounds/maxprerounds", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/nlobbt/freq", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/nlobbt/maxprerounds", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/obbt/freq", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/obbt/maxprerounds", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/probing/freq", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/probing/maxprerounds", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/pseudoobj/freq", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/pseudoobj/maxprerounds", 0));
    //if(g_options.learning > 0){
        SCIP_CALL(SCIPsetIntParam(scip, "propagating/redcost/freq", -1));
        SCIP_CALL(SCIPsetIntParam(scip, "propagating/redcost/maxprerounds", 0));
        SCIP_CALL(SCIPsetIntParam(scip, "propagating/rootredcost/freq", -1));
        SCIP_CALL(SCIPsetIntParam(scip, "propagating/rootredcost/maxprerounds", 0));
    //}
    //SCIP_CALL(SCIPsetIntParam(scip, "propagating/symmetry/freq", -1));
    //SCIP_CALL(SCIPsetIntParam(scip, "propagating/symmetry/maxprerounds", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/vbounds/freq", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "propagating/vbounds/maxprerounds", 0));
    SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/rangedrowartcons", FALSE));


    // COllect FArkas
    SCIP_CALL( SCIPsetBoolParam(scip, "lp/checkfarkas", TRUE) );
    SCIP_CALL( SCIPsetBoolParam(scip, "lp/alwaysgetduals", TRUE) );
    //Disable conflict
    SCIP_CALL(SCIPsetBoolParam(scip, "conflict/enable", FALSE));

    
    //Fix numerics
    SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol",     9e-5) );
    SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 9e-5) );
    SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon",     1e-6) );
    SCIP_CALL( SCIPsetIntParam(scip, "lp/scaling", 0) );

    

    includeLPSolveEventHandler(scip);
    register_memoscip_branchrule(scip);

    //print branch rules
    /*SCIP_BRANCHRULE** rules = SCIPgetBranchrules(scip);
    int nrules = SCIPgetNBranchrules(scip);

    for( int i = 0; i < nrules; i++ )
    {
      cout<<
            SCIPbranchruleGetName(rules[i])<<" "<<
            SCIPbranchruleGetPriority(rules[i])<<endl;
    }*/
    SCIP_VAR* slack = nullptr;

    int nconss = SCIPgetNConss(scip);
    SCIP_CONS** conss = SCIPgetConss(scip);
    int converted = 0;
    string pronbame=filepath;
    SCIP_VAR* v;
    int nvars = SCIPgetNVars(scip);

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

        SCIP_CALL( SCIPcreateVarBasic(scip, &slack, sname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
        SCIP_CALL( SCIPaddVar(scip, slack) );
        if ((!SCIPisInfinity(scip, rhs)))
        {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, slack, 1.0) );
            SCIP_CALL( SCIPchgLhsLinear(scip, cons, rhs) );
            SCIP_CALL( SCIPchgRhsLinear(scip, cons, rhs) );
            converted++;
        }
        else if ((!SCIPisInfinity(scip, lhs))){
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
    SCIP_VAR** vars = SCIPgetVars(scip);
    vars = SCIPgetVars(scip);
    bool ok = false;
    while(!ok){
        ok=true;
        ok = true;
        vars = SCIPgetVars(scip);  
        nvars = SCIPgetNVars(scip);
        for (int i = 0; i < nvars; ++i)
        {
            SCIP_VAR* v = vars[i];
            if(SCIPvarGetType(v) != SCIP_VARTYPE_INTEGER && SCIPvarGetType(v) != SCIP_VARTYPE_CONTINUOUS){
                ok=false;
                SCIP_Bool infeasible;
                SCIP_CALL( SCIPchgVarType(scip, v, SCIP_VARTYPE_INTEGER, &infeasible) );
                assert(!infeasible);
                SCIP_CALL( SCIPchgVarUb(scip, v, 1) );
            }
        }
    }
    //Creat a pool of slack variables with a dummy constraint to keep them in the problem
    for (int i = 0; i < g_options.nbslack; i++)
    {
        slack = nullptr;
        char slackname[50];
        snprintf(slackname, 50, "sc_%d", i);
        SCIP_CALL( SCIPcreateVarBasic(scip, &slack, slackname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
        SCIP_CALL(SCIPvarSetRemovable(slack,FALSE));
        SCIP_CALL( SCIPaddVar(scip, slack) );

        char consname[256];
        snprintf(consname, sizeof(consname), "scc_%d_%s", i, SCIPvarGetName(slack));
        SCIP_VAR* consvars[1];
        SCIP_Real coefs[1];
        consvars[0] = slack;
        coefs[0] = 1.0;
        SCIP_CONS* ubcons = nullptr;

        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &ubcons,
            consname,
            1,
            consvars,
            coefs,
            0,
            SCIPinfinity(scip)) );

        SCIP_CALL( SCIPaddCons(scip, ubcons) );
        SCIP_CALL( SCIPsetConsModifiable(scip,ubcons,TRUE));
        SCIP_CALL( SCIPreleaseCons(scip, &ubcons) );
        SCIP_CALL( SCIPreleaseVar(scip, &slack) );
    }

    if(fileext==".opb"){
        //If problem is in .opb format, add x_{i0}+ x_{i1} = 1 for all variables
        for (int i = 0; i < nvars-converted; ++i)
        {
            vars = SCIPgetVars(scip);
            v = vars[i];
            char consname[256];
            snprintf(consname, sizeof(consname), "ub1_%d_%s", i, SCIPvarGetName(v));
            SCIP_VAR* consvars[2];
            SCIP_Real coefs[2];
            consvars[0] = v;
            coefs[0] = 1.0;
            char sname[256];
            snprintf(sname, sizeof(sname), "%s_0", SCIPvarGetName(v));
            slack = nullptr;
            cout<<"Creating slack variable "<<sname<<" for upper bound constraint of variable "<<SCIPvarGetName(v)<<endl;
            SCIP_CALL( SCIPcreateVarBasic(scip, &slack, sname, 0.0, 1, 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, slack) );

            consvars[1] = slack;
            coefs[1] = 1.0;

            SCIP_CONS* ubcons = nullptr;

            SCIP_CALL(SCIPcreateConsLinear(scip, &ubcons, consname,
                2, consvars, coefs,
                1.0, 1.0,   /* lhs, rhs : égalité 1 <= ... <= 1 */
                TRUE,   /* initial */
                FALSE,  /* separate */
                TRUE,   /* enforce */
                TRUE,   /* check */
                TRUE,   /* propagate */
                FALSE,  /* local */
                FALSE,  /* modifiable — ou TRUE si tu veux pouvoir la modifier plus tard */
                FALSE,  /* dynamic */
                FALSE,  /* removable : non supprimable par vieillissement LP */
                FALSE   /* stickingatnode */
            ));

            SCIP_CALL( SCIPaddCons(scip, ubcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &ubcons) );
            SCIP_CALL( SCIPreleaseVar(scip, &slack) );
        }
        std::cout << "Added " << nvars-converted << " upper-bound-1 constraints (var <= 1)." << std::endl;

    }
    nvars = SCIPgetNVars(scip);
    g_memoscip.nOriginalVar=nvars-g_options.nbslack;

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


    vector<SCIP_Real> objcoefs;
    vars = SCIPgetVars(scip);
    nvars = SCIPgetNVars(scip);
    //Always prefer branching on varible x_{i1} and not x_{i0}.
    for (int i = 0; i < nvars-g_options.nbslack; ++i)
    {
        v=vars[i];
        objcoefs.push_back(SCIPvarGetObj(v));
        string varname = SCIPvarGetName(v);
        if (g_options.branchingHeuristic >= 1 && varname.size() >= 2 && ((varname.substr(varname.size() - 2) == "_0" && g_options.branchingOrder==0) || (varname.substr(varname.size() - 2) != "_0" && g_options.branchingOrder==1)))
        {
            SCIP_CALL(SCIPchgVarBranchPriority(scip, v, -1000000));
        }
        if (varname.size() >= 2 && ((varname.substr(varname.size() - 2) == "_0" && g_options.branchingOrder==0) || (varname.substr(varname.size() - 2) != "_0" && g_options.branchingOrder==1)) ) {
          g_memoscip.OrderVar.push_back((v));
        }
    }

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
    SCIP_CALL(SCIPsolve(scip));

    /*if (g_options.output > 1)
        FinishLearn(scip, 1);*/

    if (g_options.output > 0)
        SCIP_CALL(SCIPwriteTransProblem(scip, "end.lp", "lp", FALSE));

    cout<<"Solving Status : "<<SCIPgetStatus(scip)<<endl;
    cout<<"Solving time : "<<SCIPgetSolvingTime(scip)<<endl;
    cout<<"Nodes : "<<SCIPgetNNodes(scip)<<endl;
    cout<<"Primal bound : "<<SCIPgetPrimalbound(scip)<<endl;
    cout<<"Dual bound : "<<SCIPgetDualbound(scip)<<endl;
    cout<<"Number of learned cuts : "<<g_memoscip.nbAddedctr<<endl;
    cout<<"Gap : "<<SCIPgetGap(scip)<<endl;
    if(g_options.verbosity > 2){
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
    }
    // Cleanup: clear per-node bookkeeping map
    g_memoscip.nodeInfoMap.clear();
    
    SCIP_CALL(SCIPfree(&scip));
    return 0;
}
