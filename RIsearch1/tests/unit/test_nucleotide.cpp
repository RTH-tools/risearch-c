#include <gtest/gtest.h>

#include <cstring>
#include <string>
#include <string_view>
#include <vector>

#include "nucleotide.h"

TEST(NucleotideCoding, MapsEachBaseToItsIndex)
{
    EXPECT_EQ(nt2index('A'), 0);
    EXPECT_EQ(nt2index('C'), 1);
    EXPECT_EQ(nt2index('G'), 2);
    EXPECT_EQ(nt2index('T'), 3);
    EXPECT_EQ(nt2index('N'), 4);
    EXPECT_EQ(nt2index('-'), 5);
}

TEST(NucleotideCoding, MapsEachIndexBackToItsBase)
{
    // index2nt is the RNA-facing inverse of nt2index, so index 3 is U, not T.
    EXPECT_EQ(index2nt(0), 'A');
    EXPECT_EQ(index2nt(1), 'C');
    EXPECT_EQ(index2nt(2), 'G');
    EXPECT_EQ(index2nt(3), 'U');
    EXPECT_EQ(index2nt(4), 'N');
    EXPECT_EQ(index2nt(5), '-');
}

TEST(NucleotideCoding, RoundTripsThroughIndexAndBack)
{
    // T is absent: it codes to the same index as U and comes back as U.
    const std::string_view bases = "ACGUN-";
    for (char base : bases) {
        EXPECT_EQ(index2nt(nt2index(base)), base) << "base " << base;
    }
}

TEST(NucleotideCoding, Seq2ixEncodesASequenceAndReportsGapsRemoved)
{
    // seq2ix writes one index per non-gap character, compacting gaps out, and
    // returns how many it dropped. Callers subtract that from the length.
    char seq[] = "AC-GU";
    std::vector<unsigned char> out(sizeof(seq), 0xFF);

    const int gaps = seq2ix(5, seq, out.data(), (char*)"test", (char*)"query");

    EXPECT_EQ(gaps, 1);
    EXPECT_EQ(out[0], 0); // A
    EXPECT_EQ(out[1], 1); // C
    EXPECT_EQ(out[2], 2); // G, shifted down over the gap
    EXPECT_EQ(out[3], 3); // U
}
