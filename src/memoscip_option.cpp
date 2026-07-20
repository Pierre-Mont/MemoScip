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

static bool handle_short_flag(const char* s, MemoSCIPRunOptions& out, std::string& errMsg)
{
    if (s == nullptr || s[0] != '-' || s[1] == '\0')
        return false;

    /* -o / -p  */
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

    if(std::strcmp(s, "-NC") == 0)
    {
        out.NC = true;
        return true;
    }
    if(std::strcmp(s, "-BC") == 0)
    {
        out.BC = true;
        return true;
    }
    if(std::strcmp(s, "-INF") == 0)
    {
        out.INF = false;
        return true;
    }
    if (std::strcmp(s, "-kpcg") == 0)
    {
        out.kpcg = true;
        return true;
    }
    if (std::strcmp(s, "-sep") == 0)
    {
        out.sepa = true;
        return true;
    }
    /* -v... / -l... */
    if (s[1] == 'v' && (s[2] == '\0' || s[2] == '=' || (s[2] >= '0' && s[2] <= '9')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-v requires a value (e.g. -v=1 or -v1)";
            return false;
        }
        return parse_int_suffix("-v", rest, out.verbosity, errMsg);
    }
    if (s[1] == 'l' && (s[2] == '\0' || s[2] == '=' || (s[2] >= '0' && s[2] <= '9')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-l requires a value (e.g. -l=2 or -l2)";
            return false;
        }
        if (!parse_int_suffix("-l", rest, out.learning, errMsg))
            return false;
        return true;
    }

    if (s[1] == 'o' && (s[2] == '\0' || s[2] == '=' || (s[2] >= '0' && s[2] <= '9')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-o requires a positive integer value (e.g. -o=2 or -o2)";
            return false;
        }
        if (!parse_int_suffix("-o", rest, out.output, errMsg))
            return false;
        if (out.output <= 0)
        {
            errMsg = "-o requires a positive integer value (e.g. -o=2 or -o2)";
            return false;
        }
        return true;
    }

    if (s[1] == 's' && (s[2] == '\0' || s[2] == '=' || (s[2] >= '0' && s[2] <= '9')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-s requires an integer value (e.g. -s=10 or -s10)";
            return false;
        }
        if (!parse_int_suffix("-s", rest, out.nbslack, errMsg))
            return false;
        return true;
    }
    if (s[1] == 't' && (s[2] == '\0' || s[2] == '=' || (s[2] >= '0' && s[2] <= '9')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-t requires an integer value (e.g. -t=60 or -t60)";
            return false;
        }
        if (!parse_int_suffix("-t", rest, out.timelimit, errMsg))
            return false;
        if (out.timelimit < 0)
        {
            errMsg = "-t requires a non-negative integer value (e.g. -t=60)";
            return false;
        }
        return true;
    }
    
    if (s[1] == 'm' && (s[2] == '\0' || s[2] == '=' || ( (s[2] >= '0' && s[2] <= '9') || s[2] == '.' || s[2] == '+' || s[2] == '-')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-m requires a value (e.g. -m=0.5 or -m0.5)";
            return false;
        }
        if (!parse_float_suffix("-m", rest, out.minInc, errMsg))
            return false;
        return true;
    }

    if (s[1] == 'b' && s[2] == 'o' && (s[3] == '\0' || s[3] == '=' || (s[3] >= '0' && s[3] <= '9')))
    {
        const char* rest = s + 3;
        if (*rest == '\0')
        {
            errMsg = "-bo requires an integer value (e.g. -bo=1 or -bo1)";
            return false;
        }
        if (!parse_int_suffix("-bo", rest, out.branchingOrder, errMsg))
            return false;
        return true;
    }

    if (s[1] == 'c' && (s[2] == '\0' || s[2] == '=' || (s[2] >= '0' && s[2] <= '9')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-c requires an integer value (e.g. -c=1 or -c1)";
            return false;
        }
        if (!parse_int_suffix("-c", rest, out.maxCutQueueSize, errMsg))
            return false;
        return true;
    }
    
    if (s[1] == 'b' && (s[2] == '\0' || s[2] == '=' || (s[2] >= '0' && s[2] <= '9')))
    {
        const char* rest = s + 2;
        if (*rest == '\0')
        {
            errMsg = "-b requires an integer value (e.g. -b=1 or -b1)";
            return false;
        }
        if (!parse_int_suffix("-b", rest, out.branchingHeuristic, errMsg))
            return false;
        return true;
    }
    

    
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
        << "  -b=N       Branching heuristic (0=custom, 3=relpscost, 4=replay from Instance/branch_replay_sequence.txt)\n"
        << "  -bo=N      Branching order (integer, default 1)\n"
        << "  -c=N       Custom integer option (default 0)\n"
        << "\n"
        << "Examples:\n"
        << "  " << programName << " problem.lp -v=1\n"
        << "  " << programName << " problem.lp -v1 -l=2 -o -p\n";
}
