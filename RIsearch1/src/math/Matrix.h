#pragma once

#include <cstdint>

class MatrixInt {
public:
    MatrixInt(std::uint32_t rows, std::uint32_t cols);
    ~MatrixInt();
    [[nodiscard]] int** get() const;

private:
    int** m_buffer;
    std::uint32_t m_rows;
    std::uint32_t m_cols;
};

class MatrixFloat {
public:
    MatrixFloat(std::uint32_t rows, std::uint32_t cols);
    ~MatrixFloat();
    [[nodiscard]] float** get() const;

private:
    float** m_buffer;
    std::uint32_t m_rows;
    std::uint32_t m_cols;
};