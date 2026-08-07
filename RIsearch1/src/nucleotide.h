#pragma once

#include <stdio.h>

#define NEGINF INT_MIN/2

static unsigned char nt2index(char nt)
{
    switch (nt)
    {
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
    switch (ix)
    {
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


static int seq2ix (int len, const char *seq, unsigned char *retIx, const char *name, const char *type)
{
    int i;
    int gapcnt = 0;
    for (i = 0; i < len; i++)
    {
        switch (seq[i])
        {
        case 'A':
        case 'a':
            *(retIx + i - gapcnt) = 0;
            break;
        case 'C':
        case 'c':
            *(retIx + i - gapcnt) = 1;
            break;
        case 'G':
        case 'g':
            *(retIx + i - gapcnt) = 2;
            break;
        case 'T':
        case 't':
        case 'U':
        case 'u':
            *(retIx + i - gapcnt) = 3;
            break;
        case 'N':
        case 'n':
            *(retIx + i - gapcnt) = 4;
            break;
        case '-':
        case '.':		/*discard gaps from input  --  also add "case ' ' :"??? */
            gapcnt++;
            break;
        default:
            if (isalpha (seq[i]))
            {
                fprintf (stderr,
                     "Nonstandard nucleotide code '%c' in %s sequence '%s'. Replaced with 'N'\n",
                     seq[i], type, name);
                *(retIx + i - gapcnt) = 4;
                break;
            }
            else
            {			/*skip sequence!? */
                fprintf (stderr,
                     "Unexpected character '%c' in %s sequence '%s'. Skipping sequence.\n",
                     seq[i], type, name);
                return -1;
            }
        }
    }
    return gapcnt;
}
