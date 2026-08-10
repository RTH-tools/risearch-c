#pragma once

#include "memory/ByteBuffer.hpp"

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
        return !m_sequence.empty() || !m_name.empty();
    }

private:
    ByteBuffer m_sequence;
    ByteBuffer m_name;
};
