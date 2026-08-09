// Tests for the FASTA reader and the RAII wrapper around its handle.

#include <gtest/gtest.h>

#include <cstdlib>
#include <string>

#include "FastaRAII.h"

namespace {

// tests/system/data/query.fa: three 20-mers named aso1..aso3.
const char* kQuery = RISEARCH_TEST_DATA "/query.fa";


TEST(Fasta, ReadsNameAndSequenceOfTheFirstRecord)
{
    FASTAFILE* ffp = OpenFASTA(kQuery);
    ASSERT_NE(ffp, nullptr);

    char* seq = nullptr;
    char* name = nullptr;
    std::uint32_t len = 0;
    ASSERT_TRUE(ReadFASTA(ffp, &seq, &name, &len));

    EXPECT_STREQ(name, "aso1");
    EXPECT_STREQ(seq, "GCGCUGUACGAUCGAUCGAU");
    EXPECT_EQ(len, 20u);

    free(seq);
    free(name);
    CloseFASTA(ffp);
}

TEST(Fasta, ReadsEveryRecordThenStops)
{
    FASTAFILE* ffp = OpenFASTA(kQuery);
    ASSERT_NE(ffp, nullptr);

    char* seq = nullptr;
    char* name = nullptr;
    std::uint32_t len = 0;
    int count = 0;
    while (ReadFASTA(ffp, &seq, &name, &len)) {
        EXPECT_EQ(len, 20u) << "record " << count;
        free(seq);
        free(name);
        ++count;
    }

    EXPECT_EQ(count, 3);
    CloseFASTA(ffp);
}

TEST(FastaRAIITest, HandleIsUsableForAReadableFile)
{
    FastaRAII fasta(kQuery);
    ASSERT_NE(fasta.handle(), nullptr);

    char* seq = nullptr;
    char* name = nullptr;
    std::uint32_t len = 0;
    ASSERT_TRUE(ReadFASTA(fasta.handle(), &seq, &name, &len));
    EXPECT_STREQ(name, "aso1");
    free(seq);
    free(name);
}

TEST(FastaRAIITest, HandleIsNullForAMissingFile)
{
    // The destructor must cope with this: it is the case that would double-free
    // or crash if it did not check.
    FastaRAII fasta(RISEARCH_TEST_DATA "/does_not_exist.fa");
    EXPECT_EQ(fasta.handle(), nullptr);
}

TEST(FastaRAIITest, ClosesTheHandleWhenItGoesOutOfScope)
{
    // Nothing observable is returned by CloseFASTA, so this is a leak check: it
    // is meaningful under ASan/valgrind, and a crash here means double close.
    for (int i = 0; i < 100; ++i) {
        FastaRAII fasta(kQuery);
        ASSERT_NE(fasta.handle(), nullptr);
    }
}

} // namespace
