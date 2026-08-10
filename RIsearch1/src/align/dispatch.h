#pragma once

#include <cstdint>

#include "align/force_start.h"
#include "align/linspace.h"
#include "cli.h"
#include "memory/ByteBuffer.hpp"

static void run_alignment(const ByteBuffer& query_seq, const ByteBuffer& target_seq,
                          short (&dsm)[6][6][6][6], const char* nameQ, const char* nameT,
                          const config_st& config)
{
    // can't align empty sequences
    if (query_seq.empty() || target_seq.empty()) {
        return;
    }

    if (uses_force_start(config)) {
        RIs_force_start_end_init(config.force_start_val, config.pos_weights, query_seq, target_seq,
                                 dsm, config.mat_name);
    } else {
        RIs_linSpace(query_seq, target_seq, dsm, config.min_score, nameQ, nameT, config);
    }
}
