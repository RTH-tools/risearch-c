#pragma once

#include <vector>
#include <string>


namespace {
const char* kQuery = RISEARCH_TEST_DATA "/query.fa";
const char* kTarget = RISEARCH_TEST_DATA "/target.fa";


std::vector<std::string> BaseArgs(const char* min_score)
{
    return {"-q", kQuery, "-t",        kTarget, "-s", min_score, "-d",
            "30", "-m",   "su95_noGU", "-n",    "0",  "-R",      "-p2"};
}
} // namespace
