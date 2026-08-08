#pragma once

#include "memory/MallocRAII.hpp"

#include "fasta/fasta.h"


class FastaRecord {
public:
    FastaRecord() = default;

    bool read(FASTAFILE* file);


    [[nodiscard]] const char* get_sequence() const;
    [[nodiscard]] const char* get_name() const;
    [[nodiscard]] std::uint32_t get_size() const;

    [[nodiscard]] explicit operator bool() const
    {
        return m_sequence.get() != nullptr;
    }

private:
    MallocRAII<char> m_sequence;
    MallocRAII<char> m_name;
    std::uint32_t m_size = 0;
};
