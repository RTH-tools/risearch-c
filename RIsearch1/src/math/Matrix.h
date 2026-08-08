#pragma once

#include <cstdint>

class MatrixInt {
public:
    MatrixInt(int rows, int cols);
    ~MatrixInt();
    [[nodiscard]] int** get() const;

private:
    int** m_buffer;
    std::uint32_t m_rows;
    std::uint32_t m_cols;
};

class MatrixFloat {
public:
    MatrixFloat(int rows, int cols);
    ~MatrixFloat();
    [[nodiscard]] float** get() const;

private:
    float** m_buffer;
    std::uint32_t m_rows;
    std::uint32_t m_cols;
};