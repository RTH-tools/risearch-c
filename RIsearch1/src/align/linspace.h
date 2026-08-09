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
    for (auto i =1u; i <= m; ++i) {
        Iy[0][i] = M[0][i] = NEGINF;
        Ix[0][i] = 0;
    }



    /* init j=1 row */
    /*init first col (i=0) */
    Iy[1][0] = 0; /* MAX(0, dsm[GAP][GAP][GAP][tseq[n-1]]); */ /*n-1 is last nt in target; used to
                                                                  be 0 after reversion */
    Ix[1][0] = M[1][0] = NEGINF;

    M[1][1] = dsm[GAP][query_sequence[0]][GAP][target_sequence[n - 1]];
    /*    maxval = M[1][1] + MAX(0,dsm[qseq[0]][GAP][tseq[n-1]][GAP]); ---- why should we allow a
     * special treatment here!? */
    auto rowMax_score = M[1][1] + dsm[query_sequence[0]][GAP][target_sequence[n - 1]][GAP];
    auto rowMax_pos = 1;

    /* (1,1) cell can not be in Ix or Iy state. */
    Ix[1][1] = Iy[1][1] = NEGINF;

    for (auto i = 2u; i <= m; i++) {
        /* value for M matrix, case we have a pair here (k=0) */
        M[1][i] = dsm[GAP][query_sequence[i - 1]][GAP][target_sequence[n - 1]];
        /*had been 3 possibilities before, removed case
         Ix[0][i-1] != 0 ? Ix[0][i-1] + dsm[qseq[i-2]][qseq[i-1]][GAP][tseq[n-1]] : -1,
         as above all Ix[0][i] are set to 0
        */

        /* max so far? */
        if (const auto testmax =
                M[1][i] + dsm[query_sequence[i - 1]][GAP][target_sequence[n - 1]][GAP];
            (testmax > rowMax_score)) {
            rowMax_score = testmax;
            rowMax_pos = i;
        }

        /* value for Ix matrix, case query sequence paired to gap (k=1) */
        /* removed one option, namely: start new alignment that starts in gap, reflected by (-, Xi;
         * -, -)  */
        Ix[1][i] = max3(
            0,
            /* prev. match, now gap (no match possible before!?)  - add (Xi-1, Xi; Y1, -) */
            M[1][i - 1] != 0
                ? M[1][i - 1] +
                      dsm[query_sequence[i - 2]][query_sequence[i - 1]][target_sequence[n - 1]][GAP]
                : -1,
            /* extending existing gap  - add (Xi-1, Xi; -, -) */
            Ix[1][i - 1] != 0
                ? Ix[1][i - 1] + dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][GAP]
                : -1
            /* OR start new alignment that starts in gap, reflected by (-, Xi; -, -)
               dsm[GAP][qseq[i-1]][GAP][GAP] */
        );
        /* do NOT allow max other than match state!? -- however (Xi, -; -, -) is positive!?
            testmax = Ix[1][i] + dsm[qseq[i-1]][GAP][GAP][GAP];
            if (testmax > maxval) {
              maxval = testmax;
              maxi = i; maxj = 1;
            }
        */
        /* value for Iy matrix, case target sequence paired to gap (k=2) */
        /* not possible in this row */
        Iy[1][i] = NEGINF;
    }

    /* initialization of first rows and columns completed
       recursion to complete alignment with two residues follows
     */

    /*no need to store rowMax_score in the first place, can access hits_score[] just as fast!? (keep
     * pointer to current) */
    /* Raw pointers for the sweep: reading through the owner means the compiler
       reloads its member across the calls in the loop. */
    int* const hs = hits_score.get();
    int* const hp = hits_pos.get();

    hs[0] = rowMax_score;
    hp[0] = rowMax_pos;

    running_max.set(rowMax_score, rowMax_pos, 1);

    for (auto j = 2u; j <= n; j++) { /* new row */
        const auto currentRow = j % 2;
        const auto lastRow = 1 - currentRow; // If currentRow 1 --> 0, if 0 --> 1

        /*changes for pruning:
           switch rows and col / mat[i][j] becomes mat[j][i]
           matIndex  j-1       becomes  lastRow
           matIndex   j        becomes  currentRow
           last nt  tseq[j-2]  becomes  tseq[n-j+1]
           this nt  tseq[j-1]  becomes  tseq[n-j]
         */

        /* TODO better solved with additional running pointer to *(n-j) ?? */

        /* init i=0 column */
        /* no need do do this all over - will not change - no diff if NEGINF or 0 - so already set
           by init of row j=0 and j=1 Ix[currentRow][0] = M[currentRow][0] = NEGINF;
           Iy[currentRow][0] = 0; *//*MAX(0, dsm[GAP][GAP][tseq[n-j+1]][tseq[n-j]]); */

        /* init i=1 column */

        const auto target_current = target_sequence[n - j];
        const auto target_prev = target_sequence[n - j + 1];

        /* value for M matrix, case we have a pair here (k=0) */
        M[currentRow][1] = MAX(0, dsm[GAP][query_sequence[0]][GAP][target_current]);
        /* starting a new alignment with (-, X1; -, Yj)  OR  0 if not possible */
        /* coming from gap in tseq NOT possible as we're looking at first position of query! */
        /*had been 3 possibilities before, removed case 'coming from gap in qseq (Iy-matrix), adding
           (-, X1; Yj-1, Yj)' Iy[lastRow][0] != 0 ? Iy[lastRow][0] +
           dsm[GAP][qseq[0]][tseq[n-j+1]][tseq[n-j]] : -1, as above all Iy[j][0] are set to 0
         */

        rowMax_score =
            M[currentRow][1] +
            dsm[query_sequence[0]][GAP][target_current][GAP]; /*stricter only if not 0 before!? */
        rowMax_pos = 1;
        /* this would be a single pair w/o stacking, could save this check as well!? */

        /* value for Ix matrix, case query sequence paired to gap (k=1) */
        /* not possible in this column */
        Ix[currentRow][1] = NEGINF;

        /* value for Iy matrix, case target sequence paired to gap (k=2) */
        /* removed one option, namely: start new ali starting w/ gap, (-, -; Yj, -)  //
         * dsm[GAP][GAP][GAP][tseq[n-j]] */
        Iy[currentRow][1] = MAX(
            /*coming from match, opening gap - add (X1, -; Yj-1, Yj) */
            M[lastRow][1] + dsm[query_sequence[0]][GAP][target_prev][target_current],
            /* extending an existing gap - add (-, -; Yj-1, Yj) */
            Iy[lastRow][1] + dsm[GAP][GAP][target_prev][target_current]);
        /*  do NOT allow max other than match state!? -- however (-, -; *, -) is positive!?
            testmax = Iy[currentRow][1] + dsm[GAP][GAP][tseq[n-j]][GAP];
            if (testmax > maxval) {
              maxval = testmax;
              maxi = 1; maxj = j;
            }
        */
        /* finished init of i=1 col */

        for (auto i = 2u; i <= m; i++) { /* cols */ /* alt bed: *p1  */

            /* value for M matrix, case we have a pair here (k=0) */
            M[currentRow][i] = max4(
                /* coming from a match, add (Xi-1, Xi; Yi-1, Yi) */
                M[lastRow][i - 1] != 0
                    ? M[lastRow][i - 1] + dsm[query_sequence[i - 2]][query_sequence[i - 1]]
                                             [target_sequence[n - j + 1]][target_sequence[n - j]]
                    : -1,
                /* coming from gap in target, add (Xi-1, Xi; -, Yi) */
                Ix[lastRow][i - 1] +
                    dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][target_current],
                /* coming from gap in query, add (-, Xi; Yi-1, Yi) */
                Iy[lastRow][i - 1] + dsm[GAP][query_sequence[i - 1]][target_prev][target_current],
                /* starting a new alignment with this pair: (-, Xi; -, Yj) */
                dsm[GAP][query_sequence[i - 1]][GAP][target_current]);
            if (const auto testmax =
                    M[currentRow][i] + dsm[query_sequence[i - 1]][GAP][target_current][GAP];
                testmax > rowMax_score) {
                rowMax_score = testmax;
                rowMax_pos = i;
            }
            /* value for Ix matrix, case query paired to gap (k=1) */
            /* removed one option, namely: start new alignment that starts in gap, reflected by (-,
             * Xi; -, -) */
            Ix[currentRow][i] = MAX(
                /*coming from match, add (Xi-1, Xi; Yj, -) */
                M[currentRow][i - 1] +
                    dsm[query_sequence[i - 2]][query_sequence[i - 1]][target_current][GAP],
                /*extend existing gap, add (Xi-1, Xi; -, -) */
                Ix[currentRow][i - 1] + dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][GAP]
                /* start new alignment - starts with gap LEGAL!?
                   dsm[GAP][qseq[i-1]][GAP][GAP] */
            );
            /* only M[][] can be max, point of backtrack...
                  testmax = Ix[currentRow][i] + dsm[qseq[i-1]][GAP][GAP][GAP];
                  if (testmax > maxval) {
                    maxval = testmax;
                    maxi = i; maxj = j;
                  }
            */
            /* value for Iy matrix, case target paired to gap (k=2) */
            /* removed one option, namely: start new ali starting w/ gap, (-, -; Yj, -)  //
             * dsm[GAP][GAP][GAP][tseq[n-j]] */
            Iy[currentRow][i] = MAX(
                /*coming from match, add (Xi, -; Yj-1, Yj) */
                M[lastRow][i] + dsm[query_sequence[i - 1]][GAP][target_prev][target_current],
                /*extend existing gap, add (-, -; Yj-1, Yj) */
                Iy[lastRow][i] + dsm[GAP][GAP][target_prev][target_current]
                /* start new alignment - starts with gap LEGAL!?
                   dsm[GAP][GAP][GAP][tseq[n-j]] */
            );
            /*  only M[][] can be max, point of backtrack...
                  testmax = Iy[currentRow][i] + dsm[GAP][GAP][tseq[n-j]][GAP];
                  if (testmax > maxval) {
                    maxval = testmax;
                    maxi = i; maxj = j;
                  }
            */
        }

        hs[j - 1] = rowMax_score;
        hp[j - 1] = rowMax_pos;

        running_max.set_if_better(rowMax_score, rowMax_pos, j);

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
