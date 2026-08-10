

#include <gtest/gtest.h>

#include "fasta.h"
#include <cstdlib>
#include <string>

std::string WriteFasta(const char* stem, const std::string& content)
{
    const std::string path = testing::TempDir() + "/risearch_" + stem + ".fa";
    FILE* f = std::fopen(path.c_str(), "wb");
    EXPECT_NE(f, nullptr) << path;
    std::fwrite(content.data(), 1, content.size(), f);
    std::fclose(f);
    return path;
}

TEST(Fasta, OpenReturnsNullForAMissingFile)
{
    EXPECT_EQ(OpenFASTA(RISEARCH_TEST_DATA "/does_not_exist.fa"), nullptr);
}

TEST(Fasta, OpenEmptyFileReturnsNull)
{
    EXPECT_EQ(OpenFASTA(WriteFasta("empty", "").c_str()), nullptr);
}

// Reads every record, for the death tests below.
static void ReadAll(const std::string& path)
{
    FASTAFILE* ffp = OpenFASTA(path.c_str());
    if (ffp == nullptr)
        return;
    char* seq = nullptr;
    char* name = nullptr;
    std::uint32_t len = 0;
    while (ReadFASTA(ffp, &seq, &name, &len)) {
        free(seq);
        free(name);
    }
    CloseFASTA(ffp);
}

// Each of these was dropped silently, leaving a shorter sequence that still
// looked valid and scored as if it were the real one. Note the split string
// literals: a hex escape swallows every hex digit that follows it, so
// "\xffACGU" is one escape rather than a byte and four letters.

TEST(FastaDeathTest, RejectsAByteAbove127)
{
    const std::string path = WriteFasta("high_bit", ">gene1\nACGU\xff" "ACGU\n");
    EXPECT_EXIT(ReadAll(path), ::testing::ExitedWithCode(1), "Corrupt byte 0xff");
}

TEST(FastaDeathTest, RejectsUtf8InASequence)
{
    // A file that went through an editor and picked up an accented character.
    const std::string path = WriteFasta("utf8_seq", ">gene1\nACGU\xc3\xa9" "ACGU\n");
    EXPECT_EXIT(ReadAll(path), ::testing::ExitedWithCode(1), "Corrupt byte 0xc3");
}

TEST(FastaDeathTest, RejectsAControlCharacter)
{
    // What the tail of a truncated file tends to look like.
    const std::string path = WriteFasta("control", ">gene1\nACGU\x01" "ACGU\n");
    EXPECT_EXIT(ReadAll(path), ::testing::ExitedWithCode(1), "Corrupt byte 0x01");
}

TEST(FastaDeathTest, NamesTheRecordItRejected)
{
    // A 40k-record file gives nothing to go on without the name.
    const std::string path = WriteFasta("named", ">first\nACGU\n>second\nAC\xff" "GU\n");
    EXPECT_EXIT(ReadAll(path), ::testing::ExitedWithCode(1), "of 'second'");
}

// The other half of the rule: the formatting a FASTA file is allowed to carry
// is still accepted, so the check cannot be tightened by accident.
TEST(Fasta, AcceptsGapsCoordinateColumnsAndCarriageReturns)
{
    const std::string path =
        WriteFasta("tolerated", ">gaps\r\nAC--GU..AC\r\n>coords\r\n   1 ACGU GGCC\r\n");
    FASTAFILE* ffp = OpenFASTA(path.c_str());
    ASSERT_NE(ffp, nullptr);

    char* seq = nullptr;
    char* name = nullptr;
    std::uint32_t len = 0;

    ASSERT_TRUE(ReadFASTA(ffp, &seq, &name, &len));
    EXPECT_STREQ(seq, "ACGUAC");
    free(seq);
    free(name);

    ASSERT_TRUE(ReadFASTA(ffp, &seq, &name, &len));
    EXPECT_STREQ(seq, "ACGUGGCC");
    free(seq);
    free(name);

    CloseFASTA(ffp);
}
