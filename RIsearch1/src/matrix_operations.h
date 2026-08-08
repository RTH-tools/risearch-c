#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>

static int** allocIntMatrix(std::uint32_t rows, std::uint32_t cols)
{
    /* function to allocate an Integer Matrix of size rows x cols */
    int** m;

    m = malloc(rows * sizeof(int*));

    if (!m) {
        printf("Cannot allocate integer matrix with %d rows\n", rows);
        exit(1);
    }

    for (auto i = 0u; i < rows; i++) {
        m[i] = malloc(cols * sizeof(int));
        if (!m[i]) {
            printf("Cannot allocate column %d of matrix with %d rows and %d cols\n", i, rows, cols);
            exit(1);
        }
    }
    return m;
}

static float** allocFloatMatrix(std::uint32_t rows, std::uint32_t cols)
{
    /* function to allocate an Integer Matrix of size rows x cols */
    float** m;

    m = malloc(rows * sizeof(float*));

    if (!m) {
        printf("Cannot allocate integer matrix with %d rows\n", rows);
        exit(1);
    }

    for (auto i = 0u; i < rows; i++) {
        m[i] = malloc(cols * sizeof(float));
        if (!m[i]) {
            printf("Cannot allocate column %d of matrix with %d rows and %d cols\n", i, rows, cols);
            exit(1);
        }
    }
    return m;
}

static void freeIntMatrix(int** m, std::uint32_t rows)
{
    while (rows--)
        free((char*)(m[rows]));
    free((char*)(m));
}

static void freeFloatMatrix(float** m, std::uint32_t rows)
{
    while (rows--)
        free((char*)(m[rows]));
    free((char*)(m));
}

static void fill_char_array(char* buf, std::uint32_t length)
{
    for (auto i = 0u; i < length; i++) {
        buf[i] = 'X';
    }
    buf[length] = '\0';
}
