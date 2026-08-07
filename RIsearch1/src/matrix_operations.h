#pragma once

#include <stdlib.h>

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

unsigned char * allocate_char_array (int length)
{
    char *string;
    int i;

    string = malloc (length + 1 * sizeof (char));
    for (i = 0; i < length; i++)
    {
        string[i] = 'X';
    }
    string[length] = '\0';
    return string;
}
