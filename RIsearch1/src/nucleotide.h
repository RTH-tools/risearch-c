#pragma once

#include <climits>
#include <cctype>
#include <cstdio>

#include <cstdint>

#include "memory/ByteBuffer.hpp"

#define GAP 5 /* position of '-' in alphabet, not as define if read from matrix... */

#define NEGINF INT_MIN / 2

[[maybe_unused]] static unsigned char nt2index(char nt)
{
    switch (nt) {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
    case 'U':
    case 'u':
        return 3;
    case 'N':
    case 'n':
        return 4;
    case '-':
        return 5; /*any case!? */
    default:
        fprintf(stderr, "Nonstandard nucleotide code: %c\n", nt);
        return 4;
    }
}

static char index2nt(unsigned char ix)
{
    switch (ix) {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'U';
    case 4:
        return 'N';
    case 5:
        return '-'; /*any case!? */
    default:
        fprintf(stderr, "\nUnknown symbol >>%d<< found.\n", ix);
        exit(1); /*just skip it!? not fail? */
    }
}


/* Encodes a sequence into matrix indices, dropping gap characters. The buffer
 * ends up exactly as long as the encoding, so there is no gap count to hand back
 * and no length for the caller to correct. False means the sequence held a
 * character that is not a nucleotide code at all, and should be skipped.
 */
static bool seq2ix(std::uint32_t len, const char* seq, ByteBuffer& retIx, const char* name,
                   const char* type)
{
    retIx.clear();

    for (auto i = 0u; i < len; i++) {
        switch (seq[i]) {
        case 'A':
        case 'a':
            retIx.push_back(0);
            break;
        case 'C':
        case 'c':
            retIx.push_back(1);
            break;
        case 'G':
        case 'g':
            retIx.push_back(2);
            break;
        case 'T':
        case 't':
        case 'U':
        case 'u':
            retIx.push_back(3);
            break;
        case 'N':
        case 'n':
            retIx.push_back(4);
            break;
        case '-':
        case '.': /*discard gaps from input*/
            break;
        default:
            if (isalpha(seq[i])) {
                fprintf(stderr,
                        "Nonstandard nucleotide code '%c' in %s sequence '%s'. Replaced with 'N'\n",
                        seq[i], type, name);
                retIx.push_back(4);
                break;
            } else { /*skip sequence!? */
                fprintf(stderr,
                        "Unexpected character '%c' in %s sequence '%s'. Skipping sequence.\n",
                        seq[i], type, name);
                return false;
            }
        }
    }
    return true;
}
