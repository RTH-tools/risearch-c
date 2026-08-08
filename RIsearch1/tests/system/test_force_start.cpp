// End-to-end tests for the -f / -w path (forced start and position weights).
//
// The seven parameter sets used to check the search against a reference binary
// never touch this path, so it went through the header split, two signature
// changes, an energy-formula rewrite and an allocation rewrite without once
// being executed. These tests keep that from happening again.

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "risearch_runner.h"

namespace {

const char *kQuery = RISEARCH_TEST_DATA "/query.fa";
const char *kTarget = RISEARCH_TEST_DATA "/target.fa";

// -f must exceed 200 * max(query, target) or the search reports that it could
// not force the start; the committed target is 95 nt, so 20000 is ample.
std::vector<std::string> ForceStartArgs(const char *weights, const char *matrix)
{
    return {"-q", kQuery, "-t", kTarget, "-m", matrix, "-f", "20000", "-w", weights};
}

TEST(ForceStart, RunsWithoutWeights)
{
    const std::string out = risearch_test::Run(ForceStartArgs("noweights", "su95_noGU"));
    EXPECT_NE(out.find("***Structures and Energies***"), std::string::npos);
    EXPECT_NE(out.find("Free energy [kcal/mol]"), std::string::npos);
}

TEST(ForceStart, RunsWithCrisprWeights)
{
    const std::string out = risearch_test::Run(ForceStartArgs("CRISPR_20nt_5p_3p", "su95_noGU"));
    EXPECT_NE(out.find("***Structures and Energies***"), std::string::npos);
}

TEST(ForceStart, WeightsChangeTheResult)
{
    // The weights scale the per-position contributions, so the two must differ.
    EXPECT_NE(risearch_test::Run(ForceStartArgs("noweights", "su95_noGU")),
              risearch_test::Run(ForceStartArgs("CRISPR_20nt_5p_3p", "su95_noGU")));
}

TEST(ForceStart, WorksForEveryMatrix)
{
    for (const char *matrix : {"su95_noGU", "t04", "slh04_noGU"}) {
        const std::string out = risearch_test::Run(ForceStartArgs("noweights", matrix));
        EXPECT_NE(out.find("***Structures and Energies***"), std::string::npos)
            << "matrix " << matrix;
    }
}

TEST(ForceStart, ReportsEveryQueryAgainstTheTarget)
{
    const std::string out = risearch_test::Run(ForceStartArgs("noweights", "su95_noGU"));
    for (const char *q : {"aso1", "aso2", "aso3"})
        EXPECT_NE(out.find(q), std::string::npos) << "missing query " << q;
}

TEST(ForceStart, IsDeterministic)
{
    EXPECT_EQ(risearch_test::Run(ForceStartArgs("noweights", "su95_noGU")),
              risearch_test::Run(ForceStartArgs("noweights", "su95_noGU")));
}

TEST(ForceStartDeathTest, RejectsAForceValueThatIsTooSmall)
{
    // Below 200 * max(m, n) the search cannot pin the start, and says so rather
    // than reporting an interaction that does not begin where it was asked to.
    std::vector<std::string> args = {"-q",        kQuery, "-t", kTarget, "-m",
                                     "su95_noGU", "-f",   "10", "-w",    "noweights"};
    EXPECT_EXIT(risearch_test::Run(args), ::testing::ExitedWithCode(1),
                "Force start option did not work");
}

} // namespace
