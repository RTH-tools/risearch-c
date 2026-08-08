// Runs RIsearch end to end inside the test process.
//
// The library target compiles risearch.c with main() renamed to risearch1_main,
// so a test can invoke the whole program the way the shell would and inspect what
// it printed -- no subprocess, no temp binary, no shell script.

#pragma once

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <unistd.h> // optind

// risearch.cpp is compiled as C++ with main renamed, so this has C++ linkage.
int risearch1_main(int argc, char* argv[]);

// getopt keeps its parse position in a global. Setting it to 0 makes glibc
// re-initialise completely on the next call, which is what lets a second test
// parse its own arguments correctly.

namespace risearch_test {

// Runs RIsearch with `args` (argv[0] is supplied) and returns everything it
// wrote to stdout.
inline std::string Run(const std::vector<std::string>& args)
{
    std::vector<std::string> owned;
    owned.reserve(args.size() + 1);
    owned.push_back("RIsearch");
    for (const auto& a : args)
        owned.push_back(a);

    std::vector<char*> argv;
    argv.reserve(owned.size() + 1);
    for (auto& s : owned)
        argv.push_back(s.data());
    argv.push_back(nullptr);

    optind = 0; // rewind getopt before each run
    testing::internal::CaptureStdout();
    risearch1_main(static_cast<int>(argv.size()) - 1, argv.data());
    return testing::internal::GetCapturedStdout();
}

// Number of non-empty lines, i.e. the hit count in the -p2 output.
inline int CountLines(const std::string& out)
{
    int n = 0;
    for (size_t i = 0, start = 0; i <= out.size(); ++i)
        if (i == out.size() || out[i] == '\n') {
            if (i > start)
                ++n;
            start = i + 1;
        }
    return n;
}

} // namespace risearch_test
