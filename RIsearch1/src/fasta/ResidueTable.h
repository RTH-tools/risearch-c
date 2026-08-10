#pragma once

/* Which bytes of a sequence line are residues, i.e. what the reader keeps.
 * isalpha() answers the same question, but it is a locale-dependent call per
 * character; this is a load from a table that stays in L1. It also indexes on
 * unsigned char, where isalpha() on a plain char passes a negative value for
 * any byte above 127, which is undefined.
 */
struct ResidueTable {
    unsigned char is[256] = {};

    constexpr ResidueTable()
    {
        for (int c = 'A'; c <= 'Z'; c++)
            is[c] = 1;
        for (int c = 'a'; c <= 'z'; c++)
            is[c] = 1;
    }
};

inline constexpr ResidueTable kResidue{};
