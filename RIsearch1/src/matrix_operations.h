#pragma once

#include <stdio.h>

#include <stdlib.h>

#include "memory/MallocRAII.hpp"

static int** allocIntMatrix(int rows, int cols)
{
    /* function to allocate an Integer Matrix of size rows x cols */
    int** m;
    int i;

    m = malloc(rows * sizeof(int*));

    if (!m)
    {
        printf("Cannot allocate integer matrix with %d rows\n", rows);
        exit(1);
    }

    for (i = 0; i < rows; i++)
    {
        m[i] = malloc(cols * sizeof(int));
        if (!m[i])
        {
            printf
            ("Cannot allocate column %d of matrix with %d rows and %d cols\n",
             i, rows, cols);
            exit(1);
        }
    }
    return m;
}

static float** allocFloatMatrix(int rows, int cols)
{
    /* function to allocate an Integer Matrix of size rows x cols */
    float** m;
    int i;

    m = malloc(rows * sizeof(float*));

    if (!m)
    {
        printf("Cannot allocate integer matrix with %d rows\n", rows);
        exit(1);
    }

    for (i = 0; i < rows; i++)
    {
        m[i] = malloc(cols * sizeof(float));
        if (!m[i])
        {
            printf
            ("Cannot allocate column %d of matrix with %d rows and %d cols\n",
             i, rows, cols);
            exit(1);
        }
    }
    return m;
}

static void freeIntMatrix(int** m, int rows)
{
    while (rows--)
        free((char*)(m[rows]));
    free((char*)(m));
}

static void freeFloatMatrix(float** m, int rows)
{
    while (rows--)
        free((char*)(m[rows]));
    free((char*)(m));
}

static void fill_char_array(char* buf, int length) {
    for (auto i = 0u; i < length; i++)
    {
        buf[i] = 'X';
    }
    buf[length] = '\0';
}
