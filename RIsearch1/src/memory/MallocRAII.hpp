#pragma once

#include <cstddef>

template<typename T>
class MallocRAII {
public:
    explicit MallocRAII(size_t n) : m_buffer(static_cast<T *>(malloc(n * sizeof(T)))) {}

    ~MallocRAII() { free(m_buffer); }

    MallocRAII(const MallocRAII &) = delete;

    MallocRAII &operator=(const MallocRAII &) = delete;

    T *get() const { return m_buffer; }

private:
    T *m_buffer;
};
