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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <unistd.h>


#include "fasta.h"
#include "nucleotide.h"
#include "dsm.h"
#include "cli.h"
#include "force_start.h"
#include "linspace.h"

#include "memory/MallocRAII.hpp"
#include "FastaRAII.h"

/* values filled in by getArgs from the command line */
static config_st config;

int main(int argc, char *argv[])
{
    unsigned long len_seq1, len_seq2;
    char *one, *two;

    char *nameQ, *nameT;
    short dsm[6][6][6][6];
    int check;
    int count_q = 0, count_t = 0;

    getArgs(argc, argv, &config);

    getMat(config.mat_name, &dsm[0][0][0][0], config.extension_penalty,
           config.transpose_matrix_flag);

    if (config.seq2_file_name) {
        /* target given as file - or STDIN */
        FastaRAII target_fasta(config.seq2_file_name);

        if (nullptr == target_fasta.handle()) {
            fprintf(stderr, "Target file %s is not readable\n", config.seq2_file_name);
            return -1;
        }
        while (ReadFASTA(target_fasta.handle(), &two, &nameT, &len_seq2)) {
            count_q = 0;
            count_t++;
            MallocRAII<unsigned char> tseqIx(len_seq2);
            /*can be done already when reading in first place */
            check = seq2ix(len_seq2, two, tseqIx.get(), nameT, "target");
            if (check > 0)
                len_seq2 -= check; /*removed gap characters */
            if (check < 0)
                continue; /*non-alpha char in input */
            free(two);    /*free'ing space for full seq, as we have it as ix version */

            if (config.seq1_file_name) {
                /* query given as file */
                FastaRAII query_fasta(config.seq1_file_name);
                if (nullptr == query_fasta.handle()) {
                    fprintf(stderr, "Query file %s is not readable\n", config.seq1_file_name);
                    free(nameT);
                    return -1;
                }
                while (ReadFASTA(query_fasta.handle(), &one, &nameQ, &len_seq1)) {
                    count_q++;
                    MallocRAII<unsigned char> qseqIx(len_seq1);

                    check = seq2ix(len_seq1, one, qseqIx.get(), nameQ, "query");
                    if (check > 0)
                        len_seq1 -= check; /*removed gap characters */
                    if (check < 0)
                        continue; /*non-alpha char in input */
                    free(one);    /*free'ing space for full seq, as we have it as ix version */
                    if (config.printShort < 2 && (config.all_vs_all || count_t == count_q))
                        printf("\n\nquery %d: %s (%lu nts) vs. target %d: %s (%lu nts)\n\n",
                               count_q, nameQ, len_seq1, count_t, nameT, len_seq2);
                    if (config.weighted_positions || (config.force_start_val >= 0)) {
                        if (config.force_start_val < 0) {
                            fprintf(stderr, "Parameter -f must be set when using weights (-w).\n");
                            exit(1);
                        }
                        if (!config.weighted_positions) {
                            fprintf(stderr,
                                    "Parameter -w must be set when using force start (-f). Use "
                                    "array of weights \"noweigths\" to avoid this error.\n");
                            exit(1);
                        }
                        if (config.extension_penalty || config.tblen != 40 || config.doSubopt ||
                            config.filter_e || config.printShort || config.vicinity) {
                            fprintf(stderr, "Options -d -s -n -l -e -p are not available in "
                                            "combination with options -f -w \n");
                            exit(1);
                        }
                        if (config.all_vs_all || count_t == count_q) {
                            RIs_force_start_end_init(config.force_start_val, config.pos_weights,
                                                     qseqIx.get(), tseqIx.get(), len_seq1, len_seq2,
                                                     dsm, config.mat_name);
                        }
                    } else {
                        if (config.all_vs_all || count_t == count_q) {
                            RIs_linSpace(qseqIx.get(), tseqIx.get(), len_seq1, len_seq2, dsm,
                                         config.extension_penalty, config.min_score, nameQ, nameT,
                                         config.mat_name, &config);
                        }
                    }
                    free(nameQ);
                }
            } else if (config.seq1_cli) {
                /* query given as command line parameter */
                len_seq1 = strlen(config.seq1_cli);
                MallocRAII<unsigned char> qseqIx(len_seq1);
                check =
                    seq2ix(len_seq1, config.seq1_cli, qseqIx.get(), "from command line", "query");
                if (check > 0)
                    len_seq1 -= check; /*removed gap characters */
                if (check < 0)
                    return -1; /*non-alpha char in input -- break would loop through all query seqs,
                                  no use */
                if (config.printShort < 2)
                    printf("\n\nquery from_cli (%lu nts) vs. target %s (%lu nts)\n\n", len_seq1,
                           nameT, len_seq2);
                if (config.weighted_positions || (config.force_start_val >= 0)) {
                    if (config.force_start_val < 0) {
                        fprintf(stderr, "Parameter -f must be set when using weights (-w).\n");
                        exit(1);
                    }
                    if (!config.weighted_positions) {
                        fprintf(stderr, "Parameter -w must be set when using force start (-f). Use "
                                        "array of weights \"noweigths\" to avoid this error.\n");
                        exit(1);
                    }
                    if (config.extension_penalty || config.tblen != 40 || config.doSubopt ||
                        config.filter_e || config.printShort || config.vicinity) {
                        fprintf(stderr, "Options -d -s -n -l -e -p are not available in "
                                        "combination with options -f -w \n");
                        exit(1);
                    }
                    RIs_force_start_end_init(config.force_start_val, config.pos_weights,
                                             qseqIx.get(), tseqIx.get(), len_seq1, len_seq2, dsm,
                                             config.mat_name);
                } else {
                    RIs_linSpace(qseqIx.get(), tseqIx.get(), len_seq1, len_seq2, dsm,
                                 config.extension_penalty, config.min_score, "from_cli", nameT,
                                 config.mat_name, &config);
                }
            } else {
                fprintf(stderr, "No query seq given!");
                /* is caught in getArg already -- alternative run seq against itself!? */
            }

            free(nameT);
        }
    } else if (config.seq2_cli) {
        /*target given as command line parameter */
        MallocRAII<unsigned char> tseqIx(len_seq2);
        len_seq2 = strlen(config.seq2_cli);
        check = seq2ix(len_seq2, config.seq2_cli, tseqIx.get(), "from command line", "target");
        if (check > 0)
            len_seq2 -= check; /*removed gap characters */
        if (check < 0)
            return -1; /*non-alpha char in input */

        if (config.seq1_file_name) {
            /* query given as file */

            FastaRAII query_fasta(config.seq1_file_name);
            if (nullptr == query_fasta.handle()) {
                fprintf(stderr, "Query file %s is not readable\n", config.seq1_file_name);
                return -1;
            }
            while (ReadFASTA(query_fasta.handle(), &one, &nameQ, &len_seq1)) {
                MallocRAII<unsigned char> qseqIx(len_seq1);

                check = seq2ix(len_seq1, one, qseqIx.get(), nameQ, "query");
                if (check > 0)
                    len_seq1 -= check; /*removed gap characters */
                if (check < 0)
                    continue; /*non-alpha char in input */
                free(one);

                if (config.printShort < 2)
                    printf("\n\nquery %s (%lu nts) vs. target from_cli (%lu nts)\n\n", nameQ,
                           len_seq1, len_seq2);
                if (config.weighted_positions || (config.force_start_val >= 0)) {
                    if (config.force_start_val < 0) {
                        fprintf(stderr, "Parameter -f must be set when using weights (-w).\n");
                        exit(1);
                    }
                    if (!config.weighted_positions) {
                        fprintf(stderr, "Parameter -w must be set when using force start (-f). Use "
                                        "array of weights \"noweigths\" to avoid this error.\n");
                        exit(1);
                    }
                    if (config.extension_penalty || config.tblen != 40 || config.doSubopt ||
                        config.filter_e || config.printShort || config.vicinity) {
                        fprintf(stderr, "Options -d -s -n -l -e -p are not available in "
                                        "combination with options -f -w \n");
                        exit(1);
                    }
                    RIs_force_start_end_init(config.force_start_val, config.pos_weights,
                                             qseqIx.get(), tseqIx.get(), len_seq1, len_seq2, dsm,
                                             config.mat_name);
                } else {
                    RIs_linSpace(qseqIx.get(), tseqIx.get(), len_seq1, len_seq2, dsm,
                                 config.extension_penalty, config.min_score, nameQ, "from_cli",
                                 config.mat_name, &config);
                }
                free(nameQ);
            }
        } else if (config.seq1_cli) {
            /* query given as command line parameter */
            MallocRAII<unsigned char> qseqIx(len_seq1);

            len_seq1 = strlen(config.seq1_cli);
            check = seq2ix(len_seq1, config.seq1_cli, qseqIx.get(), "from command line", "query");
            if (check > 0)
                len_seq1 -= check; /*removed gap characters */
            if (check < 0)
                return -1; /* non-alpha char in input -- break would loop through queries, no use */
            if (config.printShort < 2)
                printf("\n\nquery from_cli (%lu nts) vs. target from_cli (%lu nts)\n\n", len_seq1,
                       len_seq2);
            if (config.weighted_positions || (config.force_start_val >= 0)) {
                if (config.force_start_val < 0) {
                    fprintf(stderr, "Parameter -f must be set when using weights (-w).\n");
                    exit(1);
                }
                if (!config.weighted_positions) {
                    fprintf(stderr, "Parameter -w must be set when using force start (-f). Use "
                                    "array of weights \"noweigths\" to avoid this error.\n");
                    exit(1);
                }
                if (config.extension_penalty || config.tblen != 40 || config.doSubopt ||
                    config.filter_e || config.printShort || config.vicinity) {
                    fprintf(stderr, "Options -d -s -n -l -e -p are not available in combination "
                                    "with options -f -w \n");
                    exit(1);
                }
                RIs_force_start_end_init(config.force_start_val, config.pos_weights, qseqIx.get(),
                                         tseqIx.get(), len_seq1, len_seq2, dsm, config.mat_name);
            } else {
                RIs_linSpace(qseqIx.get(), tseqIx.get(), len_seq1, len_seq2, dsm,
                             config.extension_penalty, config.min_score, "from_cli", "from_cli",
                             config.mat_name, &config);
            }
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
