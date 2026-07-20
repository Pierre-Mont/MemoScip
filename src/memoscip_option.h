#ifndef MEMOSCIP_OPTION_H
#define MEMOSCIP_OPTION_H

#include <string>

/** Parsed command-line options for memobound. */
struct MemoSCIPRunOptions {
    std::string instancePath;
    /** SCIP `display/verblevel` (0 = quiet). */
    int verbosity = 0;
    /** Learning pipeline: 0 = off, >= 1 = register learning event handlers (level for future use). */
    int learning = 0;
    /** Extra outputs: intermediate LP dumps, etc. */
    int output = 0;
    /** Reserved / future profiling hooks. */
    bool profiling = false;

    bool debug = false;
    /*Node consistency*/
    bool NC = false;
    /** Bound propagation: 0 = off, >= 1 = enable branching consistency. */
    bool BC = false;

    bool INF= true; // True = collect infeasibility certificate, False = do not collect infeasibility certificate
    float minInc = 1.0f;
    /** Treat the input instance as KPCG format and use internal reader. */
    bool kpcg = false;

    int nbslack = 1000; // Number of slack variables 
    /** Maximum time limit in seconds (0 = no limit). */
    int timelimit = 0;

    int branchingHeuristic = 0;  //0: first fraction, 1+: fallback SCIP branching, 3: relpscost, 4: replay from Instance/branch_replay_sequence.txt

    int branchingOrder = 1;    //0: branch on x_{i0} first, 1: branch on x_{i1} first

    int maxCutQueueSize = 0;

    bool sepa = false; // True = enable separation, False = disable separation
};

/** Parse argv. Expects one instance file path (any non-flag argument) and optional flags. */
bool parse_memoscip_cli(int argc, char** argv, MemoSCIPRunOptions& out, std::string& errMsg);

void print_memoscip_usage(const char* programName);


/** Single shared instance; inline (C++17) so no separate .cpp definition is required. */
inline MemoSCIPRunOptions g_options{};


#endif
