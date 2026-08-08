#pragma once

extern "C" {
#include "fasta.h"
}

class FastaRAII {
public:
    explicit FastaRAII(const char* sequence_file);

    FastaRAII(const FastaRAII&) = delete;
    FastaRAII(FastaRAII&&) = delete;
    FastaRAII& operator=(const FastaRAII&) = delete;
    FastaRAII& operator=(FastaRAII&&) = delete;

    FASTAFILE* handle() const;

    ~FastaRAII();

    const char* m_sequence_file;
    FASTAFILE* m_fasta_file_handle;
};
