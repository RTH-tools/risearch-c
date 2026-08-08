#pragma once

// ReSharper disable once CppUnusedIncludeDirective
#include <cstddef>
// ReSharper disable once CppUnusedIncludeDirective
#include <cstdlib>

template<typename T>
class MallocRAII {
public:
    MallocRAII() : m_buffer(nullptr)
    {
    }

    explicit MallocRAII(size_t n) : m_buffer(static_cast<T*>(malloc(n * sizeof(T))))
    {
    }

    ~MallocRAII()
    {
        free(m_buffer);
    }

    MallocRAII(const MallocRAII&) = delete;
    MallocRAII& operator=(const MallocRAII&) = delete;
    MallocRAII(MallocRAII&& rhs) noexcept : m_buffer(rhs.m_buffer)
    {
        rhs.m_buffer = nullptr;
    }
    MallocRAII& operator=(MallocRAII&& rhs) noexcept
    {
        if (this != &rhs) {
            free(m_buffer);
            m_buffer = rhs.m_buffer;
            rhs.m_buffer = nullptr;
        }
        return *this;
    }

    T* operator->() const
    {
        return m_buffer;
    }

    T& operator[](size_t i) const
    {
        return m_buffer[i];
    }

    T* get() const
    {
        return m_buffer;
    }

private:
    T* m_buffer;
};
