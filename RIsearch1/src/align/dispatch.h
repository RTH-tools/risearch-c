#pragma once

#include <cstdint>

#include "align/force_start.h"
#include "align/linspace.h"
#include "cli.h"

static void run_alignment(unsigned char* qseqIx, unsigned char* tseqIx, std::uint32_t len_seq1,
                          std::uint32_t len_seq2, short (&dsm)[6][6][6][6], const char* nameQ,
                          const char* nameT, const config_st& config)
{
    // can't align empty sequences
    if (len_seq1 == 0 || len_seq2 == 0) {
        return;
    }

    if (uses_force_start(config)) {
        RIs_force_start_end_init(config.force_start_val, config.pos_weights, qseqIx, tseqIx,
                                 len_seq1, len_seq2, dsm, config.mat_name);
    } else {
        RIs_linSpace(qseqIx, tseqIx, len_seq1, len_seq2, dsm, config.extension_penalty,
                     config.min_score, nameQ, nameT, &config);
    }
}