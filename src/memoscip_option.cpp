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
            errMsg = "-o requires a value (e.g. -o=2 or -o2)";
            return false;
        }
        if (!parse_int_suffix("-o", rest, out.output, errMsg))
            return false;
        return true;
    }

    errMsg = std::string("unknown option: ") + s;
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
        << "  -o         Extra output (write begin.lp and end.lp)\n"
        << "  -p         Profiling / reserved for future use\n"
        << "\n"
        << "Examples:\n"
        << "  " << programName << " problem.lp -v=1\n"
        << "  " << programName << " problem.lp -v1 -l=2 -o -p\n";
}
