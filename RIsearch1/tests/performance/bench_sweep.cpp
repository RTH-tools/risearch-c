// Quick throughput checks, so a change can be judged in a second rather than by
// running the full regression against a reference binary.
//
// Sequences are generated in process from a fixed seed, so there is no file I/O
// and the work is identical on every run and every machine. What is reported is
// cell updates per second: one cell is one (query position, target position)
// pair, which the recursion evaluates in three states.
//
// The assertions are deliberately loose -- they exist to catch an order of
// magnitude, not a percent. Read the printed numbers for anything finer, and
// compare against a reference binary for anything you intend to act on.

#include <gtest/gtest.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "cli.h"
#include "dsm.h"
#include "align/linspace.h"

namespace {

constexpr int kQueryLength = 20;    // an ASO
constexpr int kTargetLength = 4000; // roughly one human mRNA

// A fixed-seed generator, so every run scores the same sequences. std::mt19937
// would do, but this keeps the benchmark free of <random>'s setup cost and is
// reproducible across standard library versions.
class Lcg {
public:
    explicit Lcg(std::uint64_t seed) : m_state(seed)
    {
    }

    unsigned char next_base()
    {
        m_state = m_state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<unsigned char>((m_state >> 33) & 3); // A, C, G or U
    }

private:
    std::uint64_t m_state;
};

std::vector<unsigned char> make_sequence(int length, std::uint64_t seed)
{
    Lcg rng(seed);
    std::vector<unsigned char> seq(length);
    for (int i = 0; i < length; i++)
        seq[i] = rng.next_base();
    return seq;
}

// The production invocation: -d 30 -m su95_noGU -n 0 -p2.
config_st bench_config(int min_score)
{
    config_st config{};
    config.mat_name = "su95_noGU";
    config.extension_penalty = 30;
    config.min_score = min_score;
    config.doSubopt = 1;
    config.max_energy = INT_MAX;
    config.printShort = 2;
    config.tblen = 40;
    config.vicinity = 0;
    config.all_vs_all = 1;
    return config;
}

// Runs the search `repeats` times and returns cell updates per second. Output is
// swallowed: what is being timed is the search, not the terminal.
double measure(int min_score, int repeats)
{
    const auto query = make_sequence(kQueryLength, 0x5eed);
    const auto target = make_sequence(kTargetLength, 0xf00d);
    const config_st config = bench_config(min_score);

    short dsm[6][6][6][6];
    char matname[] = "su95_noGU";
    getMat(matname, &dsm[0][0][0][0], config.extension_penalty, 0);

    std::fflush(stdout);
    testing::internal::CaptureStdout();
    const auto start = std::chrono::steady_clock::now();
    for (int r = 0; r < repeats; r++) {
        RIs_linSpace(const_cast<unsigned char*>(query.data()),
                     const_cast<unsigned char*>(target.data()), kQueryLength, kTargetLength, dsm,
                     config.extension_penalty, config.min_score, "q", "t", &config);
    }
    const auto elapsed = std::chrono::steady_clock::now() - start;
    testing::internal::GetCapturedStdout();

    const double seconds = std::chrono::duration<double>(elapsed).count();
    const double cells = static_cast<double>(repeats) * kQueryLength * kTargetLength;
    return cells / seconds;
}

void report(const char* what, double cells_per_second)
{
    std::printf("[ PERF     ] %-22s %8.1f M cell-updates/s\n", what, cells_per_second / 1e6);
}

// The sweep on its own: the threshold is unreachable, so no hit is ever
// re-aligned or printed. This is the number the DP work shows up in, and about
// 85% of a production run.
TEST(Performance, SweepThroughput)
{
    const double rate = measure(INT_MAX, 40);
    report("sweep only", rate);
    EXPECT_GT(rate, 20e6) << "the sweep has slowed by more than an order of magnitude";
}

// Sweep plus traceback, at the cutoff the pipeline actually uses. The difference
// against the sweep-only number is what reporting hits costs.
TEST(Performance, SweepAndReportThroughput)
{
    const double rate = measure(800, 40);
    report("sweep + traceback", rate);
    EXPECT_GT(rate, 10e6) << "reporting hits has slowed by more than an order of magnitude";
}

// Traceback cost as a share of the whole, which is what caps any speedup aimed
// only at the sweep. Reported, not asserted on: the two arms are timed
// separately and the random sequences produce few hits, so the difference
// between them is within run-to-run noise and can come out either sign.
TEST(Performance, TracebackShareOfRuntime)
{
    const double sweep = measure(INT_MAX, 40);
    const double both = measure(800, 40);
    const double share = 1.0 - both / sweep;
    std::printf("[ PERF     ] %-22s %8.1f %% of runtime\n", "traceback", share * 100.0);
}

} // namespace
