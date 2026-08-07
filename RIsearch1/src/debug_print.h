#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "nucleotide.h"


static void printMat(int** mat, int rows, int cols, unsigned char* seq1,
                     unsigned char* seq2)
{
    int i, j;
    printf("\t-");
    for (j = 0; j < cols - 1; j++)
    {
        printf("\t%c", index2nt(*(seq2 + j)));
    }
    for (i = 0; i < rows; i++)
    {
        printf("\n%c", (i == 0 ? '-' : index2nt(*(seq1 + i - 1))));
        for (j = 0; j < cols; j++)
        {
            printf("\t%d", (mat[i][j] == NEGINF ? -8 : mat[i][j])); /* -8 as dummy for -inf */
        }
    }
    printf("\n\n");
}

static void printfloatMat(float** mat, int rows, int cols, unsigned char* seq1,
                          unsigned char* seq2)
{
    int i, j;
    printf("\t-");
    for (j = 0; j < cols - 1; j++)
    {
        printf("\t%c", index2nt(*(seq2 + j)));
    }
    for (i = 0; i < rows; i++)
    {
        printf("\n%c", (i == 0 ? '-' : index2nt(*(seq1 + i - 1))));
        for (j = 0; j < cols; j++)
        {
            printf("\t%f", (mat[i][j] == NEGINF ? -8 : mat[i][j])); /* -8 as dummy for -inf */
        }
    }
    printf("\n\n");
}
