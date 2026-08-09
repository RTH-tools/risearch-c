#pragma once

#include "RunningMax.h"


#include <cstdio>
#include <climits>

#include "matrix_operations.h"
#include "force_start.h"
#include "nucleotide.h"

#include "operations.h"
#include "traceback.h"
#include "memory/MallocRAII.hpp"


static void
RIs_linSpace(const unsigned char* query_sequence,  /* query sequence - numeric representation */
             const unsigned char* target_sequence, /* target sequence */
             std::uint32_t m,                      /* query seq length */
             std::uint32_t n,                      /* target seq length */
             short dsm[6][6][6][6],                /* scoring matrix -- TODO variable length!? */
             int extension_penalty,                /* as used in dsm, to calc Score2fakE */
             int threshold,                        /* give out hits higher than that */
             const char* qname,                    /* query name */
             const char* tname,                    /* target name */
             const config_st* config)
{
    // Prepare for RIs calls:
    const auto window = config->tblen + 1;
    MatrixInt M_RI_RAII(window, window);  /* (Mis)Match */
    MatrixInt Ix_RI_RAII(window, window); /* Insertion in x(=query), so x paired to gap (in y) */
    MatrixInt Iy_RI_RAII(window, window); /* Insertion(=bulge) in y(=target) */
    int** M_RI = M_RI_RAII.get();
    int** Ix_RI = Ix_RI_RAII.get();
    int** Iy_RI = Iy_RI_RAII.get();


    const auto reference = reference_from_matrix(config->mat_name);
    RunningMax running_max{};


    MallocRAII<int> hits_score(n);
    MallocRAII<int> hits_pos(n);

    MallocRAII<unsigned char> tmpQseq(config->tblen);
    MallocRAII<unsigned char> tmpTseq(config->tblen);

    const auto alignment_capacity = static_cast<int>(1.5 * config->tblen);
    IA maxHit(alignment_capacity);

    /* matrices for alignment scores ending in different states */
    // Since we only need 2 rows we can optimize the memory layout.
    MallocRAII<int> dp_rows(6 * (m + 1));
    int* const M[2] = {dp_rows.get() + 0 * (m + 1), dp_rows.get() + 1 * (m + 1)};
    int* const Ix[2] = {dp_rows.get() + 2 * (m + 1), dp_rows.get() + 3 * (m + 1)};
    int* const Iy[2] = {dp_rows.get() + 4 * (m + 1), dp_rows.get() + 5 * (m + 1)};


    M[0][0] = Ix[0][0] = Iy[0][0] = 0;

    /**
     * Row j=0 -- the empty prefix before the first target nt processed.
     * (target_sequence[n-1], the 3' end, since the target is walked backwards).
     *
     * M and Iy need a target nt to exist, so neither state is reachable here.
     *
     * Ix is 0 because an alignemnt must open on a base pair, and never on a bulge.
     *
     * This scheme can't score dangling ends.
     *
     */
    for (auto i = 1u; i <= m; ++i) {
        Iy[0][i] = M[0][i] = NEGINF;
        Ix[0][i] = 0;
    }


    /* init j=1 row */
    /*init first col (i=0) */
    Iy[1][0] = 0;
    Ix[1][0] = M[1][0] = NEGINF;

    // n - 1 is the last nt in target
    M[1][1] = dsm[GAP][query_sequence[0]][GAP][target_sequence[n - 1]];

    RunningRowMax running_row_max{};
    running_row_max.set(M[1][1] + dsm[query_sequence[0]][GAP][target_sequence[n - 1]][GAP], 1);

    /* (1,1) cell can not be in Ix or Iy state. */
    Ix[1][1] = Iy[1][1] = NEGINF;

    const auto t_last = target_sequence[n - 1];

    for (auto i = 2u; i <= m; i++) {
        const auto q_prev = query_sequence[i - 2];
        const auto q_cur = query_sequence[i - 1];

        const auto open_score = dsm[GAP][q_cur][GAP][t_last];
        const auto close_score = dsm[q_cur][GAP][t_last][GAP];

        M[1][i] = open_score;
        running_row_max.set_if_better(M[1][i] + close_score, i);

        /* removed one option, namely: start new alignment that starts in gap, reflected by (-, Xi;
         * -, -)  */
        Ix[1][i] =
            max3(0,
                 /* prev. match, now gap (no match possible before!?)  - add (Xi-1, Xi; Y1, -) */
                 M[1][i - 1] != 0 ? M[1][i - 1] + dsm[q_prev][q_cur][t_last][GAP] : -1,
                 /* extending existing gap  - add (Xi-1, Xi; -, -) */
                 Ix[1][i - 1] != 0 ? Ix[1][i - 1] + dsm[q_prev][q_cur][GAP][GAP] : -1);

        // There is no previous row, so there can't be a bulge.
        Iy[1][i] = NEGINF;
    }
    running_max.set(running_row_max.score, running_row_max.pos_i, 1);


    /**
     * The first row can only hold 1nt, so NN terms only now begin to apply.
     */

    // hs[j - 1] and hp[j - 1] hold the best alignment ending at pos j
    int* const hs = hits_score.get();
    int* const hp = hits_pos.get();

    hs[0] = running_row_max.score;
    hp[0] = running_row_max.pos_i;


    for (auto j = 2u; j <= n; j++) {
        /* Begin init of i=1 case */
        const auto currentRow = j % 2;
        const auto lastRow = 1 - currentRow;

        const auto target_current = target_sequence[n - j];
        const auto target_prev = target_sequence[n - j + 1];


        /* Column 1 is the query's first nt, nothing can precede it */
        const auto open_score_1 = dsm[GAP][query_sequence[0]][GAP][target_current];
        const auto close_score_1 = dsm[query_sequence[0]][GAP][target_current][GAP];
        M[currentRow][1] = MAX(0, open_score_1);

        running_row_max.set(M[currentRow][1] + close_score_1, 1);

        // Ix bulges a query nt, impossible in 1st nucleotide
        Ix[currentRow][1] = NEGINF;

        // Iy bulges a target nt, j>=2 so bulge possible.
        // Only M can win row max, so we don't take this as a candidate
        Iy[currentRow][1] = MAX(
            // We are now opening the bulge
            M[lastRow][1] + dsm[query_sequence[0]][GAP][target_prev][target_current],
            // We are extending a bulge
            Iy[lastRow][1] + dsm[GAP][GAP][target_prev][target_current]);

        /* finished init of i=1 col */

        for (auto i = 2u; i <= m; i++) {
            const auto q_prev = query_sequence[i - 2];
            const auto q_cur = query_sequence[i - 1];

            const auto open_score_i = dsm[GAP][q_cur][GAP][target_current];
            const auto close_score_i = dsm[q_cur][GAP][target_current][GAP];

            M[currentRow][i] = max4(
                /* coming from a match */
                M[lastRow][i - 1] != 0
                    ? M[lastRow][i - 1] +
                          dsm[q_prev][q_cur][target_sequence[n - j + 1]][target_sequence[n - j]]
                    : -1,
                /* coming from gap in target */
                Ix[lastRow][i - 1] + dsm[q_prev][q_cur][GAP][target_current],
                /* coming from gap in query */
                Iy[lastRow][i - 1] + dsm[GAP][q_cur][target_prev][target_current],
                /* start fresh */
                open_score_i);

            /* Set the best alignment ending in position i for row j*/
            running_row_max.set_if_better(M[currentRow][i] + close_score_i, i);

            /**
             * Ix: query nt against a gap
             */
            Ix[currentRow][i] = MAX(
                // pair at i - 1, now bulge
                M[currentRow][i - 1] + dsm[q_prev][q_cur][target_current][GAP],
                // already bulging, add one more
                Ix[currentRow][i - 1] + dsm[q_prev][q_cur][GAP][GAP]
            );

            /**
             * Iy: target nt against a gap
             */
            Iy[currentRow][i] = MAX(
                // pair at previous row, now bulge
                M[lastRow][i] + dsm[q_cur][GAP][target_prev][target_current],
                // already bulging, add one more
                Iy[lastRow][i] + dsm[GAP][GAP][target_prev][target_current]
            );
        }

        hs[j - 1] = running_row_max.score;
        hp[j - 1] = running_row_max.pos_i;

        running_max.set_if_better(running_row_max.score, running_row_max.pos_i, j);

    } /*next row j */


    int hitcount = 0;

    /*
      if checking for subopts && set energy threshold, do not guarantee to print best hit first, but
      only once! if checking for subopts and p[2-4], do not check best first.
    */
    if (!(config->doSubopt && (config->filter_e || config->printShort > 1))) {
        /*do backtrack for this one only! by recomputing whole matrix for this subsection */
        /* max going back config->tblen */
        const auto tmpQbeg =
            running_max.pos_i > config->tblen - 1 ? running_max.pos_i - (config->tblen - 1) : 1;
        const auto tmpTend =
            running_max.pos_j > config->tblen - 1 ? running_max.pos_j - (config->tblen - 1) : 1;
        const auto tmpQlen = running_max.pos_i - tmpQbeg + 1;
        const auto tmpTlen = running_max.pos_j - tmpTend + 1;


        for (auto i = 0; i < tmpQlen; i++) {
            tmpQseq[i] = query_sequence[tmpQbeg - 1 + i];
        }

        for (auto i = 0; i < tmpTlen; i++) {
            tmpTseq[i] = target_sequence[n - tmpTend - i];
        }


        /*strncpy will not work, as 0 is needed, no terminate etc. -
         * no need to reset string as we need to give length anyway... otherwise like follows*/
        /*  memset(input_str, '\0', sizeof( input_str )); */

        /*With -249 in the Sugimoto case the results we obtain are the same that we get
         from the scripts used for the off-target paper (without weights). 249 corresponds to the
         value -G-C in the su95 matrix, removing it means we do not want to add the energy of adding
         a firs GC on top of nothing. No initialization contribution is subtracted (which instead
         was the case for the -559 case of the Turner matrices that is -150 for the -G-C cell in the
         matrix, plus -409 as initiation contribution). In the case of the Santa Lucia DNA-DNA
         matrix we applied the same reasoning used for the Sugimoto case, therefore we do
          -363 that is the -G-C case in the matrix. In case another matrix is used options need to
         be added.*/

        RIs(tmpQseq.get(), tmpTseq.get(), tmpQlen, tmpTlen, dsm, &maxHit, config, M_RI, Ix_RI,
            Iy_RI);

        /*number of nt in ia to recalc Score2fakE - only tmp no need to store... */
        const auto energy =
            (maxHit.max + extension_penalty * maxHit.nucleotide_count() - reference) / (-100.0);


        /** TODO :: use maxHit.max OR maxval==hits_score[maxj-1] in output !?!? **/
        if (energy <= config->max_energy) {
            if (config->printShort == 1) {
                printf("%d\t%d\t%d\t%d\t%.2f\n", tmpQbeg + maxHit.qbeg - 1,
                       tmpQbeg + maxHit.qend - 1, n - running_max.pos_j + maxHit.tbeg,
                       n - running_max.pos_j + maxHit.tend, energy);
                /* to be consistent with other output:
                        printf("%d\t%d\t%d\t%d\t%.2f\t%s\n", tmpQbeg+maxHit.qbeg-1,
                   tmpQbeg+maxHit.qend-1, n-maxj+maxHit.tend, n-maxj+maxHit.tbeg, energy,
                   maxHit.ali_ia);
                */
            } else if (config->printShort == 2) {
                printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\n", qname, tmpQbeg + maxHit.qbeg - 1,
                       tmpQbeg + maxHit.qend - 1, tname, n - running_max.pos_j + maxHit.tbeg,
                       n - running_max.pos_j + maxHit.tend, running_max.score, energy);
            } else if (config->printShort == 3) {
                hitcount += 1;
            } else {
                printf("Free energy [kcal/mol]: %.2f (%d)\n", energy, running_max.score);
                /*      printf("no of nucls in ia: %lu + %lu = %lu\n", maxHit.qend-maxHit.qbeg+1 ,
                 * maxHit.tend-maxHit.tbeg+1 ,
                 * maxHit.qend-maxHit.qbeg+1+maxHit.tend-maxHit.tbeg+1); */

                printf("%d - %d\n", tmpQbeg + maxHit.qbeg - 1,
                       tmpQbeg + maxHit.qend - 1); /*alignment in seq1 from to */
                printf("%s\n%s\n%s\n", maxHit.ali_seq1.get(), maxHit.ali_ia.get(),
                       maxHit.ali_seq2.get());
                printf("%d - %d (3' <-- 5')\n", n - running_max.pos_j + maxHit.tend,
                       n - running_max.pos_j + maxHit.tbeg);
            }
        }
    }

    /* break here if -s not set */
    if (!config->doSubopt) {
        if (config->printShort == 3)
            printf("%s\t%s\t%d\n", qname, tname, hitcount);

        return;
    }

    /* handle suboptimals - BEGIN */

    /*MOST naive implementation, puts out one ia for each position in the target,
     given that it is higher than the threshold */

    auto j = n;
    /* highest array index is j-1, lowest 0 */
    /* as opposed to before were it was 1-n; so have to adjust by 1 here! */
    /* OR change indices before and have array size n+1 */
    /* what about i - still from 1-m as before !? */
    while (j--) {
        if (hits_score[j] > threshold) {
#ifdef DEBUG
            printf("j=%d with %d better than threshold %d - testing neighbors\n", j, hits_score[j],
                   threshold);
#endif
            auto tmp = MIN(config->vicinity, j);
            auto tmp_min_j = j - tmp++;
            auto locMax = 0;
            while (--tmp) {
                if (hits_score[j - tmp] > hits_score[j - locMax]) {
#ifdef DEBUG
                    printf("better result %d in distance %d\n", hits_score[j - tmp], tmp);
#endif
                    locMax = tmp;
                }
            }
            j -= locMax;

            /* maxj => j+1 ;; maxi => hits_pos[j] ;; maxval => hits_score[j] */

            /* do backtrack for this hit, recompute whole matrix in subsection config->tblen */
            const auto tmpQbeg =
                hits_pos[j] > config->tblen - 1 ? hits_pos[j] - (config->tblen - 1) : 1;
            const auto tmpTend = j + 1 > config->tblen - 1 ? j + 1 - (config->tblen - 1) : 1;
            const auto tmpQlen = hits_pos[j] - tmpQbeg + 1;
            const auto tmpTlen = j + 1 - tmpTend + 1;


            for (auto i = 0u; i < tmpQlen; i++) {
                tmpQseq[i] = query_sequence[tmpQbeg - 1 + i];
            }
            for (auto i = 0u; i < tmpTlen; i++) {
                tmpTseq[i] = target_sequence[n - tmpTend - i];
            }

            RIs(tmpQseq.get(), tmpTseq.get(), tmpQlen, tmpTlen, dsm, &maxHit, config, M_RI, Ix_RI,
                Iy_RI);

            /*number of nt in ia to recalc Score2fakE - only tmp no need to store... */
            const auto energy =
                (maxHit.max + extension_penalty * maxHit.nucleotide_count() - reference) / (-100.0);


            /* TODO anything about this or ignore !?
                  if (maxHit.max != hits_score[j]) {
                    printf("did some realignment here - probably resulting in a duplicate...\n");
                  }
            */
            /** TODO :: use maxHit.max OR hits_score[j] in output !?!? **/
            if (energy <= config->max_energy) {
                if (config->printShort == 1) {
                    printf("%d\t%d\t%d\t%d\t%.2f\t%s\n", tmpQbeg + maxHit.qbeg - 1,
                           tmpQbeg + maxHit.qend - 1, n - (j + 1) + maxHit.tend,
                           n - (j + 1) + maxHit.tbeg, energy, maxHit.ali_ia.get());
                } else if (config->printShort == 2) {
                    printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\n", qname, tmpQbeg + maxHit.qbeg - 1,
                           tmpQbeg + maxHit.qend - 1, tname, n - (j + 1) + maxHit.tbeg,
                           n - (j + 1) + maxHit.tend, hits_score[j], energy);
                } else if (config->printShort == 3) {
                    hitcount += 1;
                } else {
                    printf("Free energy [kcal/mol]: %.2f (%d)\n", energy, hits_score[j]);
                    /*      printf("no of nucls in ia: %lu + %lu = %lu\n",
                     * maxHit.qend-maxHit.qbeg+1 , maxHit.tend-maxHit.tbeg+1 ,
                     * maxHit.qend-maxHit.qbeg+1+maxHit.tend-maxHit.tbeg+1); */
                    printf("%d - %d\n", tmpQbeg + maxHit.qbeg - 1,
                           tmpQbeg + maxHit.qend - 1); /*alignment in seq1 from to */
                    printf("%s\n%s\n%s\n", maxHit.ali_seq1.get(), maxHit.ali_ia.get(),
                           maxHit.ali_seq2.get());
                    printf("%d - %d (3' <-- 5')\n", n - (j + 1) + maxHit.tend,
                           n - (j + 1) + maxHit.tbeg);
                }
            }
            tmp = j + 1 - maxHit.tend;
            j = tmp_min_j; /*next check will start at first position after(before) the range that
                              was tested here */
                           /* alternative: set to end of that backtraced alignment!? */
#ifdef DEBUG
            printf("set j=%d, equals pos in seq: %d -- what about starting at end pos of ali: %d "
                   "(j=%d) -- next attempt at j-1 resp. (pos in seq)++\n",
                   j, n - j, n - tmp, tmp);
#endif
            /*      printf ("next attempt at %d not %d", j, j-1-locMax);*/
        }
    }

    /* handle suboptimals - END */
    if (config->printShort == 3)
        printf("%s\t%s\t%d\n", qname, tname, hitcount);
}
