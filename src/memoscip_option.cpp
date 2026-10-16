#include "memoscip_option.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

static bool parse_int_suffix(const char* flagName, const char* rest, int& outVal, std::string& errMsg)
{
    if (rest == nullptr || *rest == '\0')
    {
        errMsg = std::string(flagName) + " requires a value (e.g. " + flagName + "=1)";
        return false;
    }
    if (*rest == '=')
        ++rest;
    char* end = nullptr;
    long v = std::strtol(rest, &end, 10);
    if (end == rest || (end != nullptr && *end != '\0'))
    {
        errMsg = std::string("invalid integer for ") + flagName + ": " + rest;
        return false;
    }
    outVal = static_cast<int>(v);
    return true;
}

static bool parse_float_suffix(const char* flagName, const char* rest, float& outVal, std::string& errMsg)
{
    if (rest == nullptr || *rest == '\0')
    {
        errMsg = std::string(flagName) + " requires a value (e.g. " + flagName + "=0.5)";
        return false;
    }
    if (*rest == '=')
        ++rest;
    char* end = nullptr;
    double v = std::strtod(rest, &end);
    if (end == rest || (end != nullptr && *end != '\0'))
    {
        errMsg = std::string("invalid float for ") + flagName + ": " + rest;
        return false;
    }
    outVal = static_cast<float>(v);
    return true;
}

/** True if s is the flag `name` followed by nothing, '=' or a digit (also '.', '+' or '-' if allowSign). */
static bool matches_value_flag(const char* s, const char* name, bool allowSign = false)
{
    const size_t len = std::strlen(name);
    if (std::strncmp(s, name, len) != 0)
        return false;
    const char c = s[len];
    return c == '\0' || c == '=' || (c >= '0' && c <= '9') || (allowSign && (c == '.' || c == '+' || c == '-'));
}

/** Parses the value of an integer flag "-x=N" / "-xN"; missingMsg is the error when no value is given. */
static bool parse_int_flag(const char* s, const char* name, const char* missingMsg, int& outVal, std::string& errMsg)
{
    const char* rest = s + std::strlen(name);
    if (*rest == '\0')
    {
        errMsg = missingMsg;
        return false;
    }
    return parse_int_suffix(name, rest, outVal, errMsg);
}

static bool handle_short_flag(const char* s, MemoSCIPRunOptions& out, std::string& errMsg)
{
    if (s == nullptr || s[0] != '-' || s[1] == '\0')
        return false;

    /* boolean flags */
    if (std::strcmp(s, "-p") == 0)
    {
        out.profiling = true;
        return true;
    }
    if (std::strcmp(s, "-d") == 0)
    {
        out.debug = true;
        return true;
    }
    if (std::strcmp(s, "-NC") == 0)
    {
        out.NC = true;
        return true;
    }
    if (std::strcmp(s, "-BC") == 0)
    {
        out.BC = true;
        return true;
    }
    if (std::strcmp(s, "-INF") == 0)
    {
        out.INF = false;
        return true;
    }
    if (std::strcmp(s, "-kpcg") == 0)
    {
        out.kpcg = true;
        return true;
    }
    if (std::strncmp(s, "-encoding", 9) == 0)
    {
        return parse_int_suffix("-encoding", s + 9, out.encoding, errMsg);
    }
    if (std::strcmp(s, "-sep") == 0)
    {
        out.sepa = true;
        return true;
    }
    if (std::strcmp(s, "-fl") == 0)
    {
        out.finish_Learn = true;
        return true;
    }

    /* flags with a value: -x=N or -xN (checked after the boolean flags, e.g. -sep before -s) */
    if (matches_value_flag(s, "-v"))
        return parse_int_flag(s, "-v", "-v requires a value (e.g. -v=1 or -v1)", out.verbosity, errMsg);

    if (matches_value_flag(s, "-l"))
        return parse_int_flag(s, "-l", "-l requires a value (e.g. -l=2 or -l2)", out.learning, errMsg);

    if (matches_value_flag(s, "-o"))
    {
        const char* msg = "-o requires a positive integer value (e.g. -o=2 or -o2)";
        if (!parse_int_flag(s, "-o", msg, out.output, errMsg))
            return false;
        if (out.output <= 0)
        {
            errMsg = msg;
            return false;
        }
        return true;
    }

    if (matches_value_flag(s, "-s"))
        return parse_int_flag(s, "-s", "-s requires an integer value (e.g. -s=10 or -s10)", out.nbslack, errMsg);

    if (matches_value_flag(s, "-t"))
    {
        if (!parse_int_flag(s, "-t", "-t requires an integer value (e.g. -t=60 or -t60)", out.timelimit, errMsg))
            return false;
        if (out.timelimit < 0)
        {
            errMsg = "-t requires a non-negative integer value (e.g. -t=60)";
            return false;
        }
        return true;
    }

    if (matches_value_flag(s, "-m", true))
    {
        if (s[2] == '\0')
        {
            errMsg = "-m requires a value (e.g. -m=0.5 or -m0.5)";
            return false;
        }
        return parse_float_suffix("-m", s + 2, out.minInc, errMsg);
    }

    if (matches_value_flag(s, "-ub", true))
    {
        const char* rest = s + 3;
        if (*rest == '=')
            ++rest;
        char* end = nullptr;
        out.initialUB = std::strtod(rest, &end);
        if (end == rest || *end != '\0')
        {
            errMsg = "-ub requires a numeric value (e.g. -ub=9996)";
            return false;
        }
        out.hasInitialUB = true;
        return true;
    }

    if (matches_value_flag(s, "-bo"))
        return parse_int_flag(s, "-bo", "-bo requires an integer value (e.g. -bo=1 or -bo1)", out.branchingOrder, errMsg);

    if (matches_value_flag(s, "-c"))
        return parse_int_flag(s, "-c", "-c requires an integer value (e.g. -c=1 or -c1)", out.maxCutQueueSize, errMsg);

    if (matches_value_flag(s, "-pc"))
        return parse_int_flag(s, "-pc", "-pc requires an integer value (e.g. -pc=1 or -pc1)", out.pseudoObjCut, errMsg);

    if (matches_value_flag(s, "-w"))
        return parse_int_flag(s, "-w", "-w requires an integer value (e.g. -w=1 or -w1)", out.Weaken, errMsg);

    if (matches_value_flag(s, "-b"))
        return parse_int_flag(s, "-b", "-b requires an integer value (e.g. -b=1 or -b1)", out.branchingHeuristic, errMsg);

    return false;
}

bool parse_memoscip_cli(int argc, char** argv, MemoSCIPRunOptions& out, std::string& errMsg)
{
    out = MemoSCIPRunOptions{};
    errMsg.clear();

    for (int i = 1; i < argc; ++i)
    {
        const char* s = argv[i];
        if (s == nullptr || s[0] == '\0')
            continue;

        if (s[0] != '-')
        {
            if (!out.instancePath.empty())
            {
                errMsg = "multiple instance paths (expected exactly one): " + std::string(s);
                return false;
            }
            out.instancePath = s;
            continue;
        }

        if (!handle_short_flag(s, out, errMsg))
            return false;
    }

    if (out.instancePath.empty())
    {
        errMsg = "missing instance file path";
        return false;
    }

    return true;
}

void print_memoscip_usage(const char* programName)
{
    std::cerr
        << "Usage: " << programName << " <instance.lp|opb|wcsp> [options...]\n"
        << "\n"
        << "Options:\n"
        << "  -v=N       SCIP display verbosity (integer, default 0)\n"
        << "  -l=N       Learning mode level (0 = off, >= 1 enables learning handlers)\n"
        << "  -o         Extra outputs: intermediate LP dumps, etc.\n"
        << "  -m=F       Minimum increment (float, default 1.0)\n"
        << "  -t=N       Maximum time in seconds before timeout (0 = no limit)\n"
        << "  -s=N       Number of slack variables\n"
        << "  -p         Profiling / reserved for future use\n"
        << "  -kpcg      Treat input as KPCG instance and use internal reader\n"
        << "  -encoding=<n>  KPCG encoding: 0 = direct, nonzero = tuple\n"
        << "  -b=N       Branching heuristic (0=custom, 3=relpscost, 4=replay from Instance/branch_replay_sequence.txt)\n"
        << "  -bo=N      Branching order (integer, default 1)\n"
        << "  -w=N       Weaken constraints (0 = no weakening, 1 = weaken constraints)\n"
        << "  -ub=F      Initial primal bound: objective limit F (only solutions strictly better are accepted)\n"
        << "  -pc=N      Pseudo objective cutoff proof (0 = always dive, 1 = use pseudo objective proof, default 1)\n"
        << "\n"
        << "Examples:\n"
        << "  " << programName << " problem.lp -v=1\n"
        << "  " << programName << " problem.lp -v1 -l=2 -o -p\n";
}
