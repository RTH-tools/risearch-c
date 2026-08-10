/***********************************************************
  RIsearch v 1.2   --   RNA-RNA interaction search
  Copyright 2012 Anne Wenzel <wenzel@rth.dk> (RIsearch v.1.0 and v.1.1)
  Copyright 2021 Giulia I Corsi <giulia@rth.dk> (Extension of RIsearch v.1.1 in RIsearch v.1.2)

  This file is part of RIsearch.

  RIsearch is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  RIsearch is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with RIsearch, see file COPYING.
  If not, see <http://www.gnu.org/licenses/>.

***********************************************************/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <unistd.h>


#include "nucleotide.h"
#include "dsm.h"
#include "cli.h"

#include "FastaRAII.h"
#include "FastaRecord.h"
#include "align/dispatch.h"


int main(int argc, char* argv[])
{
    /* values filled in by getArgs from the command line */
    static config_st config;

    short dsm[6][6][6][6];

    getArgs(argc, argv, &config);
    if (uses_force_start(config)) {
        validate_force_start_config(config);
    }

    getMat(config.mat_name, &dsm[0][0][0][0], config.extension_penalty,
           config.transpose_matrix_flag);

    if (config.seq2_file_name) {
        /* target given as file - or STDIN */
        FastaRAII target_fasta(config.seq2_file_name);

        if (nullptr == target_fasta.handle()) {
            fprintf(stderr, "Target file %s is not readable\n", config.seq2_file_name);
            return -1;
        }
        FastaRecord target_record;

        ByteBuffer query_seq_indices;
        ByteBuffer target_seq_indices;

        int target_count = 0;
        while (target_record.read(target_fasta.handle())) {
            target_count++;
            /*can be done already when reading in first place */
            if (!seq2ix(target_record.get_size(), target_record.get_sequence(), target_seq_indices,
                        target_record.get_name(), "target")) {
                continue; // non-alpha char in input
            }
            const auto len_seq2 = static_cast<std::uint32_t>(target_seq_indices.size());

            if (config.seq1_file_name) {
                /* query given as file */
                FastaRAII query_fasta(config.seq1_file_name);
                if (nullptr == query_fasta.handle()) {
                    fprintf(stderr, "Query file %s is not readable\n", config.seq1_file_name);
                    return -1;
                }
                FastaRecord query_record;
                auto query_count = 0;
                while (query_record.read(query_fasta.handle())) {
                    query_count++;
                    if (!seq2ix(query_record.get_size(), query_record.get_sequence(), query_seq_indices,
                                query_record.get_name(), "query")) {
                        continue; // non-alpha char in input
                    }
                    const auto len_seq1 = static_cast<std::uint32_t>(query_seq_indices.size());

                    if (!config.all_vs_all && target_count != query_count) {
                        continue;
                    }

                    if (config.printShort < 2)
                        printf("\n\nquery %d: %s (%u nts) vs. target %d: %s (%u nts)\n\n",
                               query_count, query_record.get_name(), len_seq1, target_count,
                               target_record.get_name(), len_seq2);

                    run_alignment(query_seq_indices, target_seq_indices, dsm, query_record.get_name(),
                                  target_record.get_name(), config);
                }
            } else if (config.seq1_cli) {
                /* query given as command line parameter */
                if (!seq2ix(strlen(config.seq1_cli), config.seq1_cli, query_seq_indices, "from command line",
                            "query")) {
                    return -1; /*non-alpha char in input -- break would loop through all query seqs,
                                   no use */
                }
                const auto len_seq1 = static_cast<std::uint32_t>(query_seq_indices.size());

                if (config.printShort < 2)
                    printf("\n\nquery from_cli (%u nts) vs. target %s (%u nts)\n\n", len_seq1,
                           target_record.get_name(), len_seq2);

                run_alignment(query_seq_indices, target_seq_indices, dsm, "from_cli", target_record.get_name(), config);

            } else {
                fprintf(stderr, "No query seq given!");
                /* is caught in getArg already -- alternative run seq against itself!? */
            }
        }
    } else if (config.seq2_cli) {
        /*target given as command line parameter */
        ByteBuffer target_seq_indices;
        if (!seq2ix(strlen(config.seq2_cli), config.seq2_cli, target_seq_indices, "from command line",
                    "target")) {
            return -1; /*non-alpha char in input */
        }
        const auto len_seq2 = static_cast<std::uint32_t>(target_seq_indices.size());

        if (config.seq1_file_name) {
            /* query given as file */

            FastaRAII query_fasta(config.seq1_file_name);
            if (nullptr == query_fasta.handle()) {
                fprintf(stderr, "Query file %s is not readable\n", config.seq1_file_name);
                return -1;
            }
            FastaRecord query_record;
            ByteBuffer qseqIx;
            while (query_record.read(query_fasta.handle())) {
                if (!seq2ix(query_record.get_size(), query_record.get_sequence(), qseqIx,
                            query_record.get_name(), "query")) {
                    continue; /*non-alpha char in input */
                }
                const auto len_seq1 = static_cast<std::uint32_t>(qseqIx.size());

                if (config.printShort < 2) {
                    printf("\n\nquery %s (%u nts) vs. target from_cli (%u nts)\n\n",
                           query_record.get_name(), len_seq1, len_seq2);
                }

                run_alignment(qseqIx, target_seq_indices, dsm, query_record.get_name(), "from_cli", config);
            }
        } else if (config.seq1_cli) {
            /* query given as command line parameter */

            ByteBuffer qseqIx;
            if (!seq2ix(strlen(config.seq1_cli), config.seq1_cli, qseqIx, "from command line",
                        "query")) {
                return -1; /* non-alpha char in input -- break would loop through queries, no use */
            }
            const auto len_seq1 = static_cast<std::uint32_t>(qseqIx.size());

            if (config.printShort < 2)
                printf("\n\nquery from_cli (%u nts) vs. target from_cli (%u nts)\n\n", len_seq1,
                       len_seq2);

            run_alignment(qseqIx, target_seq_indices, dsm, "from_cli", "from_cli", config);
        } else {
            fprintf(stderr, "No query seq given!");
            /* is caught in getArg already -- alternative run seq against itself!? */
        }
    } else {
        fprintf(stderr, "No target seq given!");
        /* is caught in getArg already -- alternative run seq against itself!? */
    }

    return 0;
}
