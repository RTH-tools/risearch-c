#include "Matrix.h"
#include "matrix_operations.h"

MatrixInt::MatrixInt(int rows, int cols)
    : m_buffer(allocIntMatrix(rows, cols)), m_rows(rows), m_cols(cols)
{
}

MatrixInt::~MatrixInt()
{
    freeIntMatrix(m_buffer, m_rows);
}

int** MatrixInt::get() const
{
    return m_buffer;
}


MatrixFloat::MatrixFloat(int rows, int cols)
    : m_buffer(allocFloatMatrix(rows, cols)), m_rows(rows), m_cols(cols)
{
}


MatrixFloat::~MatrixFloat()
{
    freeFloatMatrix(m_buffer, m_rows);
}

float** MatrixFloat::get() const
{
    return m_buffer;
}
