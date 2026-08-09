// Tests for set_alignment_symbols: the character shown under each base pair in
// the alignment output.

#include <gtest/gtest.h>

#include "align/symbols.h"

namespace {

char Symbol(char query_nt, char target_nt)
{
    char q = '\0';
    char t = '\0';
    set_alignment_symbols(query_nt, target_nt, &q, &t);
    EXPECT_EQ(q, t) << "outputs disagree for " << query_nt << "/" << target_nt;
    return q;
}

TEST(AlignmentSymbols, WatsonCrickPairsArePipes)
{
    EXPECT_EQ(Symbol('A', 'U'), '|');
    EXPECT_EQ(Symbol('G', 'C'), '|');
}

TEST(AlignmentSymbols, WobblePairIsW)
{
    EXPECT_EQ(Symbol('G', 'U'), 'W');
}

TEST(AlignmentSymbols, NonPairingBasesAreMismatches)
{
    EXPECT_EQ(Symbol('A', 'C'), 'M');
    EXPECT_EQ(Symbol('A', 'A'), 'M');
    EXPECT_EQ(Symbol('C', 'C'), 'M');
    EXPECT_EQ(Symbol('N', 'A'), 'M');
}

TEST(AlignmentSymbols, IsSymmetricInItsTwoBases)
{
    // Each branch tests both arguments, so order must not matter.
    for (char a : {'A', 'C', 'G', 'U', 'N'})
        for (char b : {'A', 'C', 'G', 'U', 'N'})
            EXPECT_EQ(Symbol(a, b), Symbol(b, a)) << a << "/" << b;
}

} // namespace
