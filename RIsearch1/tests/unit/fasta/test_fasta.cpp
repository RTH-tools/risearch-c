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

    ByteBuffer seq;
    ByteBuffer name;
    ASSERT_TRUE(ReadFASTA(ffp, seq, name));

    EXPECT_STREQ(name.c_str(), "aso1");
    EXPECT_STREQ(seq.c_str(), "GCGCUGUACGAUCGAUCGAU");
    EXPECT_EQ(seq.size(), 20u);

    CloseFASTA(ffp);
}

TEST(Fasta, ReadsEveryRecordThenStops)
{
    FASTAFILE* ffp = OpenFASTA(kQuery);
    ASSERT_NE(ffp, nullptr);

    ByteBuffer seq;
    ByteBuffer name;
    int count = 0;
    while (ReadFASTA(ffp, seq, name)) {
        EXPECT_EQ(seq.size(), 20u) << "record " << count;
        ++count;
    }

    EXPECT_EQ(count, 3);
    CloseFASTA(ffp);
}

TEST(FastaRAIITest, HandleIsUsableForAReadableFile)
{
    FastaRAII fasta(kQuery);
    ASSERT_NE(fasta.handle(), nullptr);

    ByteBuffer seq;
    ByteBuffer name;
    ASSERT_TRUE(ReadFASTA(fasta.handle(), seq, name));
    EXPECT_STREQ(name.c_str(), "aso1");
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
