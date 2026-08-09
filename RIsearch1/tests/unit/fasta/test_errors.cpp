

#include <gtest/gtest.h>

#include "fasta.h"
#include <fstream>
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
