#include <gtest/gtest.h>

#include <vector>
#include <string>

#include "base_args.h"
#include "risearch_runner.h"


TEST(ForceStartDeathTest, RejectsWeightsWithoutAForceStart)
{
    // -w on its own: the weights only mean anything relative to a forced start,
    // so the search refuses rather than silently ignoring them.
    std::vector<std::string> args = {"-q", kQuery,      "-t", kTarget,
                                     "-m", "su95_noGU", "-w", "noweights"};

    EXPECT_EXIT(risearch_test::Run(args), ::testing::ExitedWithCode(1),
                "Parameter -f must be set when using weights");
}

TEST(ForceStartDeathTest, RejectsAForceStartWithoutWeights)
{
    // The mirror of the -w-without-f case: -f alone is refused rather than
    // defaulting the weights to something the caller did not ask for.
    std::vector<std::string> args = {"-q", kQuery,      "-t", kTarget,
                                     "-m", "su95_noGU", "-f", "20000"};

    EXPECT_EXIT(risearch_test::Run(args), ::testing::ExitedWithCode(1),
                "Parameter -w must be set when using force start");
}

TEST(ForceStartDeathTest, RejectsAnExtensionPenaltyAlongsideForceStart)
{
    // -d is one of six options the -f/-w path cannot honour; it is rejected
    // rather than silently ignored.
    std::vector<std::string> args = {"-q", kQuery, "-t",        kTarget,     "-m", "su95_noGU",
                                     "-f", "20000", "-w",       "noweights", "-d", "30"};

    EXPECT_EXIT(risearch_test::Run(args), ::testing::ExitedWithCode(1),
                "are not available in combination");
}

TEST(ForceStartDeathTest, RejectsAnOutputModeAlongsideForceStart)
{
    // printShort is checked for being non-zero, not for a particular value, so
    // this reaches the same message through a different term of the condition.
    std::vector<std::string> args = {"-q", kQuery, "-t",        kTarget,     "-m", "su95_noGU",
                                     "-f", "20000", "-w",       "noweights", "-p2"};

    EXPECT_EXIT(risearch_test::Run(args), ::testing::ExitedWithCode(1),
                "are not available in combination");
}