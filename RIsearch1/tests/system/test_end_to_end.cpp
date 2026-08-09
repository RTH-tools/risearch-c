#include <algorithm>
#include <cstdlib>
#include <vector>
// End-to-end tests: run the whole program in-process and check what it prints.

#include <gtest/gtest.h>

#include <string>

#include "risearch_runner.h"
#include "base_args.h"


namespace {

TEST(EndToEnd, ProducesTabSeparatedHits)
{
    const std::string out = risearch_test::Run(BaseArgs("400"));
    ASSERT_FALSE(out.empty());

    // -p2 is: query, qbeg, qend, target, tbeg, tend, score, energy
    const std::string first = out.substr(0, out.find('\n'));
    EXPECT_EQ(std::count(first.begin(), first.end(), '\t'), 7);
    EXPECT_EQ(first.substr(0, 4), "aso1");
}

TEST(EndToEnd, RaisingTheScoreCutoffReportsFewerHits)
{
    const int loose = risearch_test::CountLines(risearch_test::Run(BaseArgs("400")));
    const int strict = risearch_test::CountLines(risearch_test::Run(BaseArgs("800")));

    EXPECT_GT(loose, 0);
    EXPECT_GT(strict, 0);
    EXPECT_LT(strict, loose) << "a higher -s must not report more hits";
}

TEST(EndToEnd, IsDeterministic)
{
    EXPECT_EQ(risearch_test::Run(BaseArgs("400")), risearch_test::Run(BaseArgs("400")));
}

TEST(EndToEnd, ReportsEveryQueryAgainstEveryTarget)
{
    const std::string out = risearch_test::Run(BaseArgs("400"));
    for (const char* q : {"aso1", "aso2"})
        EXPECT_NE(out.find(q), std::string::npos) << "missing query " << q;
    for (const char* t : {"gene1", "gene2"})
        EXPECT_NE(out.find(t), std::string::npos) << "missing target " << t;
}

TEST(EndToEnd, TransposeChangesTheResult)
{
    // -R swaps the roles of the two strands, so for a non-symmetric matrix the
    // reported hits must differ.
    auto without = BaseArgs("400");
    without.erase(std::find(without.begin(), without.end(), "-R"));

    EXPECT_NE(risearch_test::Run(BaseArgs("400")), risearch_test::Run(without));
}

TEST(EndToEnd, HitCountModeAgreesWithTheHitLines)
{
    auto counted = BaseArgs("400");
    counted.back() = "-p3"; // one line per pair: query, target, hit-count

    int total = 0;
    const std::string out = risearch_test::Run(counted);
    for (size_t pos = 0; (pos = out.find('\t', pos)) != std::string::npos;) {
        const size_t last = out.find('\t', pos + 1);
        if (last == std::string::npos)
            break;
        total += std::atoi(out.c_str() + last + 1);
        pos = out.find('\n', last);
        if (pos == std::string::npos)
            break;
    }
    EXPECT_EQ(total, risearch_test::CountLines(risearch_test::Run(BaseArgs("400"))));
}

TEST(EndToEnd, ReportsAnUnreadableTargetFile)
{
    auto args = BaseArgs("400");
    *std::find(args.begin(), args.end(), kTarget) = RISEARCH_TEST_DATA "/does_not_exist.fa";

    testing::internal::CaptureStderr();
    const std::string out = risearch_test::Run(args);
    const std::string err = testing::internal::GetCapturedStderr();

    EXPECT_TRUE(out.empty()) << "a run that could not open its target printed hits";
    EXPECT_NE(err.find("is not readable"), std::string::npos) << err;
}

TEST(EndToEnd, ReportsAnUnreadableQueryFile)
{
    auto args = BaseArgs("400");
    *std::find(args.begin(), args.end(), kQuery) = RISEARCH_TEST_DATA "/does_not_exist.fa";

    testing::internal::CaptureStderr();
    const std::string out = risearch_test::Run(args);
    const std::string err = testing::internal::GetCapturedStderr();

    EXPECT_TRUE(out.empty()) << "a run that could not open its target printed hits";
    EXPECT_NE(err.find("is not readable"), std::string::npos) << err;
}


} // namespace
