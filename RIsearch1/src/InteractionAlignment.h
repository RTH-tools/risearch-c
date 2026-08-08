#pragma once

#include "memory/MallocRAII.hpp"

struct IA {
    int qbeg, qend, tbeg, tend;
    int max;

    MallocRAII<char> ali_seq1, ali_seq2, ali_ia;

    [[nodiscard]] int nucleotide_count() const
    {
        return qend - qbeg + 1 + tend - tbeg + 1;
    }

};
