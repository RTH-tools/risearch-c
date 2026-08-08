#pragma once


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include "debug_print.h"
#include "nucleotide.h"

#include "InteractionAlignment.h"
#include "alignment.h"
#include "energy.hpp"
#include "matrix_operations.h"
#include "operations.h"
#include "memory/MallocRAII.hpp"

#define GAP 5 /* position of '-' in alphabet, not as define if read from matrix... */



static void fill_char_array(char* buf, std::uint32_t length)
{
    for (auto i = 0u; i < length; i++) {
        buf[i] = 'X';
    }
    buf[length] = '\0';
}


static float find_max_value_f(float** M, float** Ix, float** Iy, int* k, int* i, int j, int n,
                              short dsm[6][6][6][6], unsigned char* qseq, unsigned char* tseq)
{
    float max = 0.0;
    *i = 0; /*row in which maxval is found */
    *k = 0; /*matrix in which the maximum has been found; 0 1 2 if M Ix Iy */
    if (M[n][j] + dsm[qseq[n - 1]][GAP][tseq[j - 1]][GAP] > max) {
        max = M[n][j] + dsm[qseq[n - 1]][GAP][tseq[j - 1]][GAP];
        *i = n;
    }
    /* allow to end with a gap */
    if (Ix[n][j] + dsm[qseq[n - 1]][GAP][GAP][GAP] > max) {
        max = Ix[n][j] + dsm[qseq[n - 1]][GAP][GAP][GAP];
        *k = 1;
        *i = n;
    }
    if (Iy[n][j] + dsm[GAP][GAP][tseq[n - 1]][GAP] > max) {
        max = Iy[n][j] + dsm[GAP][GAP][tseq[n - 1]][GAP];
        *k = 2;
        *i = n;
    }
    return max;
}

static void
RIs_force_start_end_weighted(int force_start_val,
                             unsigned char* qseq,   /* query sequence - numeric representation */
                             unsigned char* tseq,   /* target sequence - reversed */
                             int m,                 /* query seq length */
                             int n,                 /* target seq length */
                             float* weights,        /*array of weights */
                             short dsm[6][6][6][6], /* scoring matrix */
                             [[maybe_unused]] IA*,               /* pointer to struct, fill results */
                             const char* matname)
{
    const auto reference = reference_from_matrix(matname);

    float **M, **Ix, **Iy; /* matrices for alignment scores ending in different states */
    int i, j, k, colj;
    float w1, w2;
    int startpos_found;     /* used in backtrack, true if start pos has been reached */
    float mVal, xVal, yVal; /* values coming from M, Ix, Iy. Not possible to start a new alignment
                               at any position, due to the force start. */

    M = allocFloatMatrix(m + 1, n + 1);  /* (Mis)Match */
    Ix = allocFloatMatrix(m + 1, n + 1); /* Insertion in x(=query), so x paired to gap (in y) */
    Iy = allocFloatMatrix(m + 1, n + 1); /* Insertion(=bulge) in y(=target) */

    Ix[0][0] = Iy[0][0] = 0;
    M[0][0] = force_start_val;
    /*init first row (j=0) -- this is COL - change throughout!? */

    /*filling column zero */
    for (i = 1; i <= m; i++) {
        Iy[i][0] = NEGINF; /* not possible before beginning of target seq */
        Ix[i][0] = 0; /*MAX(0, dsm[qseq[i-2]][qseq[i-1]][GAP][GAP]); */ /* require to always have a
                                                                           match first! */
        M[i][0] = force_start_val; /*To force interaction to start at column 1 (-N in the DNA) */
    }

    /*filling row zero */
    for (j = 1; j <= n; j++) {
        Ix[0][j] = M[0][j] = NEGINF;
        Iy[0][j] = 0; /*MAX(0, dsm[GAP][GAP][tseq[j-2]][tseq[j-1]]); */
    }

    /* The initialization of i=1 column and j=1 row have to be handled explicitly
       since at this point we do not have two residues to use.
       Handle the (1,1) cell explicitly since the boundary recursion includes (i-2) or (j-2) cases.
     */
    M[1][1] = force_start_val +
              (float)dsm[GAP][qseq[0]][GAP][tseq[0]]; /*no weight used, first base pair */
    Ix[1][1] = Iy[1][1] = NEGINF;

    /* Filling first column */
    for (i = 2; i <= m; i++) {
        /* value for M matrix, case we have a pair here (k=0) */
        /* first letter of target sequence, so it can ONLY come from gapped (upstream nucleotides in
         * query have not been used) */
        M[i][1] =
            force_start_val +
            (float)
                dsm[GAP][qseq[i - 1]][GAP][tseq[0]]; /* MAX(0, ); no weight used, first base pair */
        /* value for Ix matrix, case query sequence paired to gap (k=1) */
        mVal =             /* prev. match, now gap - add (Xi-1, Xi; Y1, -) */
            M[i - 1][1] != /*adding a gap on top of a bp, weights used */
                    force_start_val
                ? M[i - 1][1] + dsm[qseq[i - 2]][qseq[i - 1]][tseq[0]][GAP] * weights[i - 2]
                : -1;
        xVal =                  /* extending existing gap  - add (Xi-1, Xi; -, -) */
            Ix[i - 1][1] != 0 ? /*adding a gap on top of a gap, weights used */
                Ix[i - 1][1] + dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] * weights[i - 2]
                              : -1;
        /* OR start new alignment that starts in gap, reflected by (-, Xi; -, -). Forbidden! */
        /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */
        /*Ix[i][1] = max4(mVal,xVal,nVal,0); */
        Ix[i][1] = max3f(mVal, xVal, 0); /*removed nVal option */
        /* value for Iy matrix, case target sequence paired to gap (k=2) not possible in this row */
        Iy[i][1] = NEGINF;
    }

    /* Filling first row */
    for (j = 2; j <= n; j++) {
        /* value for M matrix, case we have a pair here (k=0 */
        /* coming from gap in qseq (Iy-matrix); First letter in query sequence, adding (-, X1; Yj-1,
         * Yj) */
        M[1][j] =
            dsm[GAP][qseq[0]][GAP][tseq[j - 1]]; /* MAX(0, ); no weight used, first base pair */
        /* value for Iy matrix, case target sequence paired to gap (k=2) */
        /*coming from match, opening gap - add (X1, -; Yj-1, Yj) */
        mVal = M[1][j - 1] != /*adding a gap on top of a bp, weights used */
                       0
                   ? M[1][j - 1] + (float)dsm[qseq[0]][GAP][tseq[j - 2]][tseq[j - 1]] * weights[0]
                   : -1;
        /* extending an existing gap - add (-, -; Yj-1, Yj) */
        yVal = Iy[1][j - 1] != 0 ? /*adding a gap on top of a gap, weights used */
                   Iy[1][j - 1] + dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] * weights[0]
                                 : -1;
        /* or start new ali, begin with gap -- forbidden ! */
        /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */
        Iy[1][j] = max3f(mVal, yVal, 0); /*removed nVal option */
        /* value for Ix matrix, case query sequence paired to gap (k=1) */
        /* not possible in this column */
        Ix[1][j] = NEGINF;
    }
    /* initialization of first rows and columns completed
       recursion to complete alignment with two residues follows */
    for (i = 2; i <= m; i++) {
        for (j = 2; j <= n; j++) {
            /* value for M matrix, case we have a pair here (k=0) */
            /* coming from a match, add (Xi-1, Xi; Yi-1, Yi) */
            mVal = M[i - 1][j - 1] != 0
                       ? M[i - 1][j - 1] +
                             ((float)dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] *
                              weights[i - 2])
                       : -1;
            /* coming from gap in target, add (Xi-1, Xi; -, Yi) */
            xVal =
                Ix[i - 1][j - 1] != 0
                    ? Ix[i - 1][j - 1] +
                          ((float)dsm[qseq[i - 2]][qseq[i - 1]][GAP][tseq[j - 1]] * weights[i - 2])
                    : -1;

            /* coming from gap in query, add (-, Xi; Yi-1, Yi) */
            yVal =
                Iy[i - 1][j - 1] != 0
                    ? Iy[i - 1][j - 1] +
                          ((float)dsm[GAP][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] * weights[i - 2])
                    : -1;

            /* starting a new alignment with this pair: (-, Xi; -, Yj) */
            M[i][j] = max4f(mVal, xVal, yVal, 0);
            /* value for Ix matrix, case query paired to gap (k=1) */
            /*coming from match, add (Xi-1, Xi; Yj, -) */
            mVal = /*add a bulge on top of a bp */
                M[i - 1][j] != 0
                    ? M[i - 1][j] +
                          (float)dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 1]][GAP] * weights[i - 2]
                    : -1;
            /*extend existing gap, add (Xi-1, Xi; -, -) */
            xVal =
                Ix[i - 1][j] != 0
                    ? Ix[i - 1][j] + (float)dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] * weights[i - 2]
                    : -1;
            /* start new alignment - starts with gap -> ILLEGAL! */
            /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */

            Ix[i][j] = max3f(mVal, xVal, 0); /*removed option nVal */
            /* value for Iy matrix, case target paired to gap (k=2) */
            /*coming from match, add (Xi, -; Yj-1, Yj) */
            w2 = i - 2 < m - 1 - 1 ? weights[i - 1]
                                   : weights[m - 1 - 1]; /*number of weights=length query -1) */
            mVal = M[i][j - 1] != 0
                       ? M[i][j - 1] + (float)dsm[qseq[i - 1]][GAP][tseq[j - 2]][tseq[j - 1]] *
                                           ((weights[i - 2] + w2) / 2)
                       : -1; /*use average of weights for bulges in target */
            /*extend existing gap, add (-, -; Yj-1, Yj) */
            yVal = Iy[i][j - 1] != 0
                       ? Iy[i][j - 1] + (float)dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] *
                                            ((weights[i - 2] + w2) / 2)
                       : -1; /* No weight used for continuing bulges */
            /* start new alignment - starts with gap LEGAL!? */
            /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */

            Iy[i][j] = max3f(mVal, yVal, 0); /*removed option nVal */
#ifdef DEBUG
            printf("M %f; %i; %i\n", weights[i - 2], i, j);
            printf("Y %f; %i; %i\n", weights[i - 2], i, j);
            printf("X %f; %i; %i\n", weights[i - 2], i, j);
#endif
        }
    }
#ifdef DEBUG
    printf("M matrix:\n");
    printfloatMat(M, m + 1, n + 1, qseq, tseq);
    printf("Bq matrix:\n");
    printfloatMat(Ix, m + 1, n + 1, qseq, tseq);
    printf("Bt matrix:\n");
    printfloatMat(Iy, m + 1, n + 1, qseq, tseq);
#endif
    printf("***Structures and Energies***\n");
    /*backtrack */
    for (colj = n; colj <= n;
         colj++) { /*n is the number of columns, m is the number of rows.  NOTE: here we put colj =
                      n because we only want the interaction that ends at position n (ends at the 5'
                      of the target sequence, which is in 3'-5' direction) */
                   /*printf("%d\n",m);
                      printf("%d\n",n); */

        /*used to store query and target alignment symbols in backtracking */
        MallocRAII<char> query_alignment(m + 1);
        fill_char_array(query_alignment.get(), m);
        MallocRAII<char> target_alignment(n + 1);
        fill_char_array(target_alignment.get(), n);

        /*printf("col : %d\n",colj); */
        const auto max_score = find_max_value_f(
            M, Ix, Iy, &k, &i, colj, m, dsm, qseq,
            tseq); /*given arrays representing columns, we want to find the row position in which
                      the max score is reached. K is 0-1-2 M Ix Iy - should always be 0 to begin
                      with. NOTE: in this new version we search only the maximums in row m (last
                      row) because we want to force the interaction to end at 3' of query */
        /*note that i and k are modified in place in the function, since we are passing the pointers
         */
        j = colj;
        startpos_found = 0;
#ifdef DEBUG
        printf("max score %f\n", max_score);
#endif
        while (i > 0 && j > 0 && !startpos_found) {
            /*printf("%s\n",query_alignment);
               printf("%s\n",target_alignment); */
            /* find which cell gave score */
            if (k == 0) { /* highest score in M matrix, having a match! */
                if (M[i][j] == (dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]]) + force_start_val &&
                    j == 1) { /*      started new alignment */
                    /*printf("means we started alignment at row [%d] col [%d] of matrix %d \n", i,
                     * j, k); */
                    startpos_found = 1;
                } else if (M[i][j] == (dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]]) &&
                           j != 1) { /*      started new alignment */
                    /*printf("means we started alignment at row [%d] col [%d] of matrix %d \n", i,
                     * j, k); */
                    fprintf(stderr,
                            "Force start option did not work: try to increase the number given to "
                            "-f.\nTry with a number > %d\n",
                            max3(2000, m * 500, n * 500));
                    exit(1);
                } else if (M[i][j] ==
                               M[i - 1][j - 1] +
                                   ((float)dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] *
                                    weights[i - 2]) &&
                           i != 1 && j != 1) {
                    /*printf("previous was also match\n"); */ /*      previous was also match */
                } else if (M[i][j] ==
                               Ix[i - 1][j - 1] + dsm[qseq[i - 2]][qseq[i - 1]][GAP][tseq[j - 1]] *
                                                      weights[i - 2] &&
                           i != 1) {
                    /*printf("coming from gap in seq2/target\n"); */ /*      coming from gap in
                                                                        seq2/target */
                    k = 1;
                } else if (M[i][j] ==
                               Iy[i - 1][j - 1] + dsm[GAP][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] *
                                                      weights[i - 2] &&
                           j != 1) {
                    /*printf("coming from gap in seq1/query\n"); */ /*      coming from gap in
                                                                       seq1/query */
                    k = 2;
                } else if (M[i][j] == (dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]]) + force_start_val &&
                           j != 1) { /*      started new alignment */
                    fprintf(stderr,
                            "Force start option did not work: try to increase the number given to "
                            "-f.\nTry with a number > %d\n",
                            max3(2000, m * 500, n * 500));
                    exit(1);
                } else {
                    fprintf(stderr, "ERR: unexpected value in k=0.\n");
                    exit(1);
                }
                set_alignment_symbols(index2nt(qseq[i - 1]), index2nt(tseq[j - 1]),
                                      &query_alignment[i - 1], &target_alignment[j - 1]);
                i = i - 1;
                j = j - 1;
            } else if (k == 1) { /* seq1(query) paired to a gap (in target). Bulge on top of bp */
                if (Ix[i][j] ==
                    M[i - 1][j] +
                        (float)dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 1]][GAP] * weights[i - 2]) {
                    k = 0; /* open a new gap coming from match */
                           /*printf("Gap in Ix coming from match %f\n", Ix[i][j]); */
                } else {
                    if (Ix[i][j] == Ix[i - 1][j] + (float)dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] *
                                                       weights[i - 2]) { /*gap on top of gap */
                        if (i == 1) {
                            k = 3;
                            startpos_found = 1;
                        } else {
                            k = 1; /* extend existing gap */
                        }
                        /*printf("Gap in Ix coming from gap %f\n", Ix[i][j]); */
                    } else {
                        fprintf(stderr, "ERR: unexpected case in k=1 : %f in row %d and col %d\n",
                                Ix[i][j], i, j);
                        exit(1);
                    }
                }
                query_alignment[i - 1] = 'B';
                i = i - 1;
            } else if (k == 2) { /* seq2(target) paired to a gap (in query). bulge on top of bp. */
                if (i == 1) {
                    w1 = weights[0];
                    w2 = weights[0];
                } else {
                    w1 = weights[i - 2];
                    w2 = i - 2 < m - 1 - 1
                             ? weights[i - 1]
                             : weights[m - 1 - 1]; /*number of weights=length query -1) */
                }
                if (Iy[i][j] ==
                    M[i][j - 1] +
                        (float)dsm[qseq[i - 1]][GAP][tseq[j - 2]][tseq[j - 1]] * ((w1 + w2) / 2)) {
                    k = 0; /* open a new gap coming from match */
                           /*printf("Gap in Iy coming from match\n"); */
                } else {
                    if (Iy[i][j] == Iy[i][j - 1] + (float)dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] *
                                                       ((w1 + w2) / 2)) {
                        if (j == 1) {
                            k = 3;
                            startpos_found = 1;
                        } else {
                            k = 2; /* extend existing gap */
                        }
                    } else {
                        fprintf(stderr, "unexpected case in k=2 : %f\n", Iy[i][j]);
                        exit(1);
                    }
                }
                target_alignment[j - 1] = 'B';
                j = j - 1;
            } else {
                fprintf(stderr, "\nThis should really NEVER happen!\n");
                exit(1);
            }
        }
        if (!startpos_found) {
            if (j > 0) {
                fprintf(stderr,
                        "Force start option did not work: try to increase the number given to "
                        "-f.\nTry with a number > %d.\n",
                        (int)max3(2000, m * 200, n * 200));
                exit(1);
            }
        }
        printf("%s\t%s\n", query_alignment.get(), target_alignment.get());

        const auto energy =
            (max_score - force_start_val - reference) / (-100.0);
        printf("Free energy [kcal/mol] (No extension penalty): %.2f (%f)\n", energy,
               max_score - force_start_val);
    }

    freeFloatMatrix(M, m + 1);
    freeFloatMatrix(Ix, m + 1);
    freeFloatMatrix(Iy, m + 1);
}

static void
RIs_force_start_end_init(int force_start_val, const char* pos_weights,
                         unsigned char* qseqIx, /* query sequence - numeric representation */
                         unsigned char* tseqIx, /* target sequence  */
                         std::uint32_t len_seq1,          /* query seq length */
                         std::uint32_t len_seq2,          /* target seq length */
                         short dsm[6][6][6][6], /* scoring matrix */
                         const char* matname    /* name of the scoring matrix */
)
{
    int testmax;

    IA maxHit;
    MallocRAII<unsigned char> tmp(len_seq2);
    testmax = (int)(1.5 * len_seq2);
    maxHit.ali_seq1 = MallocRAII<char>(testmax);
    maxHit.ali_seq2 = MallocRAII<char>(testmax);
    maxHit.ali_ia = MallocRAII<char>(testmax);

    /*reverting seq2 (target) */
    for (auto i = 0u; i < len_seq2; i++) {
        tmp[i] = tseqIx[len_seq2 - 1 - i];
    }
    if (strcmp(pos_weights, "noweights") && force_start_val == 0) {
        fprintf(stderr,
                "ERR: weighted interactions must be forced to start at position 0 in the target. "
                "Please increase -f parameter, or use array \"noweights\" for option -w.\n");
        exit(1);
    }
    if (!strcmp(pos_weights, "CRISPR_20nt_5p_3p")) {
        extern float wC20_5p_3p[19];
        extern int size_wC20_5p_3p;
        if (size_wC20_5p_3p < len_seq1 - 1) {
            fprintf(stderr, "ERR: the array of weights is too short for the given query.\n Weights "
                            "length mush be >= than the length of the query -1.\n");
            exit(1);
        }
        if (size_wC20_5p_3p == len_seq1 - 1) {
            RIs_force_start_end_weighted(force_start_val, qseqIx, tmp.get(), len_seq1, len_seq2,
                                         &wC20_5p_3p[0], dsm, &maxHit, matname);
        } else {
            RIs_force_start_end_weighted(force_start_val, qseqIx, tmp.get(), len_seq1, len_seq2,
                                         &wC20_5p_3p[size_wC20_5p_3p + 1 - len_seq1], dsm, &maxHit,
                                         matname);
        }
    } /*
         else if (!strcmp(pos_weights, "test")){
         extern float test[4];
         extern int size_test;
         if (size_test<len_seq1-1){
         fprintf(stderr,"ERR: the array of weights is too short for the given query.\n Weights
         length mush be >= than the length of the query -1.\n"); exit(1);
         }
         if (size_test==len_seq1-1){
         RIs_force_start_end_weighted(qseqIx,tmp,len_seq1,len_seq2,&test[0],size_test,dsm,maxHit,
         matname);
         }
         else{
         RIs_force_start_end_weighted(qseqIx,tmp,len_seq1,len_seq2,&test[size_test+1-len_seq1],len_seq1-1,dsm,maxHit,
         matname);
         }
         } */
    else if (!strcmp(pos_weights, "noweights")) {
        MallocRAII<float> noweight(len_seq1 - 1);
        for (auto i = 0u; i < (len_seq1 - 1); ++i) {
            noweight[i] = 1.0;
        }
        RIs_force_start_end_weighted(force_start_val, qseqIx, tmp.get(), len_seq1, len_seq2,
                                     &noweight[0], dsm, &maxHit, matname);
    } else {
        fprintf(stderr, "Undefined weights array. Existing weights verctors are CRISPR_20nt_5p_3p "
                        "and noweights. To add a new weights vector, create an array in weights.c "
                        "and an ad-hoc if clause. Use noweights to set all weights to 1.\n");
        exit(1);
    }
}
