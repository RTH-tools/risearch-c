#pragma once

#include "memory/MallocRAII.hpp"

struct IA {
    int qbeg, qend, tbeg, tend;
    int max;

    MallocRAII<char> ali_seq1, ali_seq2, ali_ia;

    explicit IA(int alignment_capacity)
        : qbeg(0), qend(0), tbeg(0), tend(0), max(0), ali_seq1(alignment_capacity),
          ali_seq2(alignment_capacity), ali_ia(alignment_capacity)
    {
    }

    [[nodiscard]] int nucleotide_count() const
    {
        return qend - qbeg + 1 + tend - tbeg + 1;
    }
};
