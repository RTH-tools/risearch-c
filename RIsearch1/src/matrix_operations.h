#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>

static int** allocIntMatrix(std::uint32_t rows, std::uint32_t cols)
{
    const auto pointer_bytes = static_cast<std::size_t>(rows) * sizeof(int*);
    const auto data_bytes = static_cast<std::size_t>(rows) * cols * sizeof(int);

    int** m = reinterpret_cast<int**>(malloc(pointer_bytes + data_bytes));
    if (!m) {
        printf("Cannot allocate integer matrix with %u rows and %u cols\n", rows, cols);
        exit(1);
    }

    int* data = reinterpret_cast<int*>(m + rows);
    for (auto i = 0u; i < rows; i++) {
        m[i] = data + static_cast<std::size_t>(i) * cols;
    }
    return m;
}

static float** allocFloatMatrix(std::uint32_t rows, std::uint32_t cols)
{
    const auto pointer_bytes = static_cast<std::size_t>(rows) * sizeof(float*);
    const auto data_bytes = static_cast<std::size_t>(rows) * cols * sizeof(float);

    float** m = reinterpret_cast<float**>(malloc(pointer_bytes + data_bytes));
    if (!m) {
        printf("Cannot allocate float matrix with %u rows and %u cols\n", rows, cols);
        exit(1);
    }

    float* data = reinterpret_cast<float*>(m + rows);
    for (auto i = 0u; i < rows; i++) {
        m[i] = data + static_cast<std::size_t>(i) * cols;
    }
    return m;
}

static void freeIntMatrix(int** m, std::uint32_t)
{
    free(m);
}

static void freeFloatMatrix(float** m, std::uint32_t)
{
    free(m);
}
