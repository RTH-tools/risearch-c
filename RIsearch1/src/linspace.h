#pragma once

#include <cstdio>
#include <cstring>
#include <climits>
#include "debug_print.h"

#include "matrix_operations.h"
#include "force_start.h"
#include "nucleotide.h"
#include "string_util.h"
#include "operations.h"
#include "math/Matrix.h"
#include "memory/MallocRAII.hpp"


static void RIs(const unsigned char* query_seq,  /* query sequence - numeric representation */
                const unsigned char* target_seq, /* target sequence - reversed */
                int m,                     /* query seq length */
                int n,                     /* target seq length */
                short dsm[6][6][6][6],     /* scoring matrix */
                IA* hit,                   /* pointer to struct, fill results */
                const config_st* config)
{
    int maxi, maxj;   /* k is 0,1,2 for M,Ix,Iy !? - max will never be found in gapped anyway!? */
    int maxval, maxk; /* maxk not needed!? k itself could even be char... */
    int mVal, xVal, yVal, nVal; /* values coming from M, Ix, Iy, or starting a NEW alignment */

    MatrixInt M_RAII(m + 1, n + 1);  /* (Mis)Match */
    MatrixInt Ix_RAII(m + 1, n + 1); /* Insertion in x(=query), so x paired to gap (in y) */
    MatrixInt Iy_RAII(m + 1, n + 1); /* Insertion(=bulge) in y(=target) */
    int** M = M_RAII.get();
    int** Ix = Ix_RAII.get();
    int** Iy = Iy_RAII.get();

    maxi = maxj = maxk = 0;

    M[0][0] = Ix[0][0] = Iy[0][0] = 0;

    /*init first row (j=0) -- this is COL - change throughout!? */

    /* explicitly handling of the boundary condition since i-2 is not defined for i = 1
       NOT needed anymore, is now just 0 anyway
       Ix[1][0] = 0; //MAX(0, dsm[GAP][qseq[0]][GAP][GAP]);
       Iy[1][0] = M[1][0] = NEGINF; // not possible to occure
     */
    for (auto i = 1; i <= m; i++) {
        Iy[i][0] = M[i][0] = NEGINF; /* not possible before beginning of target seq */
        Ix[i][0] = 0; /*MAX(0, dsm[qseq[i-2]][qseq[i-1]][GAP][GAP]); */ /* require to always have a
                                                                           match first! */
    }

    /*init first col (i=0) --- this is ROW - change throughout!?
       Iy[0][1] = 0; // MAX(0, dsm[GAP][GAP][GAP][tseq[0]]);
       Ix[0][1] = M[0][1] = NEGINF;
     */
    for (auto j = 1; j <= n; j++) {
        Ix[0][j] = M[0][j] = NEGINF;
        Iy[0][j] = 0; /*MAX(0, dsm[GAP][GAP][tseq[j-2]][tseq[j-1]]); */
    }

    /* The initialization of i=1 column and j=1 row have to be handled explicitly
       since at this point we do not have two residues to use.

       Handle the (1,1) cell explicitly since the boundary recursion includes (i-2) or (j-2) cases.
     */

    M[1][1] =
        dsm[GAP][query_seq[0]][GAP][target_seq[0]]; /*MAX(0,dsm[GAP][qseq[0]][GAP][tseq[0]]); */
    maxval =
        M[1][1] +
        dsm[query_seq[0]][GAP][target_seq[0]][GAP]; /* MAX(0,dsm[qseq[0]][GAP][tseq[0]][GAP]); */
    /* (1,1) cell can not be in Ix or Iy state. */
    Ix[1][1] = Iy[1][1] = NEGINF;

    /* init j=1 row */
    for (auto i = 2; i <= m; i++) {
        /* value for M matrix, case we have a pair here (k=0) */
        /* first letter of target sequence, so it can ONLY come from gapped */
        /* if previous cell has not been 0, than add current pair */
        /* if the case is removed on top, it can be here as well, save checks...

           xVal = Ix[i-1][0] != 0 ? Ix[i-1][0] + dsm[qseq[i-2]][qseq[i-1]][GAP][tseq[0]] : -1; //
           reflects (i-1,i;-,1)  coming from Ix (seq1 paired to gap) otherwise - could start NEW ali
           here OR 0 if nothing else scores positive nVal = dsm[GAP][qseq[i-1]][GAP][tseq[0]];
           M[i][1] = max3(xVal,nVal,0); */
        M[i][1] = dsm[GAP][query_seq[i - 1]][GAP][target_seq[0]]; /* MAX(0, ); */

        const auto tmp = M[i][1] + dsm[query_seq[i - 1]][GAP][target_seq[0]][GAP];
        if (tmp > maxval) {
            maxval = tmp;
            maxi = i;
            maxj = 1;
            maxk = 0;
        }

        /* value for Ix matrix, case query sequence paired to gap (k=1) */
        /* prev. match, now gap - add (Xi-1, Xi; Y1, -) */
        mVal = M[i - 1][1] != 0
                   ? M[i - 1][1] + dsm[query_seq[i - 2]][query_seq[i - 1]][target_seq[0]][GAP]
                   : -1;
        /* extending existing gap  - add (Xi-1, Xi; -, -) */
        xVal = Ix[i - 1][1] != 0 ? Ix[i - 1][1] + dsm[query_seq[i - 2]][query_seq[i - 1]][GAP][GAP]
                                 : -1;
        /* OR start new alignment that starts in gap, reflected by (-, Xi; -, -)  */
        /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */
        /*Ix[i][1] = max4(mVal,xVal,nVal,0); */
        Ix[i][1] = max3(mVal, xVal, 0); /*removed nVal option */
        /* do not allow alignment ending in Xi!
            tmp = Ix[i][1] + dsm[qseq[i-1]][GAP][GAP][GAP];   // (Xi, -; -, -)
            if (tmp > maxval) {
              maxval = tmp;
              maxi = i; maxj = 1; maxk = 1;
            }
        */
        /* value for Iy matrix, case target sequence paired to gap (k=2) */
        /* not possible in this row */
        Iy[i][1] = NEGINF;
    }

    /* init i=1 column */
    for (auto j = 2; j <= n; j++) {
        /* value for M matrix, case we have a pair here (k=0) */
        /* coming from gap in qseq (Iy-matrix), adding (-, X1; Yj-1, Yj) */
        /* not possible, set to 0 before for all of them! */
        /* yVal = Iy[0][j-1] != 0 ? Iy[0][j-1] + dsm[GAP][qseq[0]][tseq[j-2]][tseq[j-1]] : -1; */
        /* OR starting a new alignment with (-, X1; -, Yj) */
        /* nVal = dsm[GAP][qseq[0]][GAP][tseq[j-1]]; */
        /* coming from gap in tseq NOT possible as we're looking at first position of query!? */
        /* M[1][j] = max3(yVal,nVal,0); */
        M[1][j] = dsm[GAP][query_seq[0]][GAP][target_seq[j - 1]]; /* MAX(0, ); */

        const auto tmp = M[1][j] + dsm[query_seq[0]][GAP][target_seq[j - 1]][GAP];
        if (tmp > maxval) {
            maxval = tmp;
            maxi = 1;
            maxj = j;
            maxk = 0;
        }

        /* value for Ix matrix, case query sequence paired to gap (k=1) */
        /* not possible in this column */
        Ix[1][j] = NEGINF;

        /* value for Iy matrix, case target sequence paired to gap (k=2) */
        /*coming from match, opening gap - add (X1, -; Yj-1, Yj) */
        mVal = M[1][j - 1] != 0
                   ? M[1][j - 1] + dsm[query_seq[0]][GAP][target_seq[j - 2]][target_seq[j - 1]]
                   : -1;
        /* extending an existing gap - add (-, -; Yj-1, Yj) */
        yVal = Iy[1][j - 1] != 0
                   ? Iy[1][j - 1] + dsm[GAP][GAP][target_seq[j - 2]][target_seq[j - 1]]
                   : -1;
        /* or start new ali, begin with gap -- forbidden ! */
        /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */
        Iy[1][j] = max3(mVal, yVal, 0); /*removed option nVal */
        /* do not allow alignments ending in gap
            tmp = Iy[1][j] + dsm[GAP][GAP][tseq[j-1]][GAP];
            if (tmp > maxval) {
              maxval = tmp;
              maxi = 1; maxj = j; maxk = 2;
            }
        */
    }
    /* initialization of first rows and columns completed
       recursion to complete alignment with two residues follows
     */
    for (auto i = 2; i <= m; i++) {
        /*alt bed: *p1  */
        for (auto j = 2; j <= n; j++) {
            /* value for M matrix, case we have a pair here (k=0) */
            /* coming from a match, add (Xi-1, Xi; Yi-1, Yi) */
            mVal = M[i - 1][j - 1] != 0
                       ? M[i - 1][j - 1] + dsm[query_seq[i - 2]][query_seq[i - 1]]
                                              [target_seq[j - 2]][target_seq[j - 1]]
                       : -1;
            /* coming from gap in target, add (Xi-1, Xi; -, Yi) */
            xVal = Ix[i - 1][j - 1] != 0
                       ? Ix[i - 1][j - 1] +
                             dsm[query_seq[i - 2]][query_seq[i - 1]][GAP][target_seq[j - 1]]
                       : -1;
            /* coming from gap in query, add (-, Xi; Yi-1, Yi) */
            yVal = Iy[i - 1][j - 1] != 0
                       ? Iy[i - 1][j - 1] +
                             dsm[GAP][query_seq[i - 1]][target_seq[j - 2]][target_seq[j - 1]]
                       : -1;
            /* starting a new alignment with this pair: (-, Xi; -, Yj) */
            nVal = dsm[GAP][query_seq[i - 1]][GAP][target_seq[j - 1]];

            M[i][j] = max4(mVal, xVal, yVal, nVal);
            const auto tmp = M[i][j] + dsm[query_seq[i - 1]][GAP][target_seq[j - 1]][GAP];
            if (tmp > maxval) {
                maxval = tmp;
                maxi = i;
                maxj = j;
                maxk = 0;
            }

            /* value for Ix matrix, case query paired to gap (k=1) */
            /*coming from match, add (Xi-1, Xi; Yj, -) */
            mVal =
                M[i - 1][j] != 0
                    ? M[i - 1][j] + dsm[query_seq[i - 2]][query_seq[i - 1]][target_seq[j - 1]][GAP]
                    : -1;
            /*extend existing gap, add (Xi-1, Xi; -, -) */
            xVal = Ix[i - 1][j] != 0
                       ? Ix[i - 1][j] + dsm[query_seq[i - 2]][query_seq[i - 1]][GAP][GAP]
                       : -1;
            /* start new alignment - starts with gap -> ILLEGAL! */
            /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */

            Ix[i][j] = max3(mVal, xVal, 0); /*removed option nVal */
            /* do not allow alignments ending in gap
                  tmp = Ix[i][j] + dsm[qseq[i-1]][GAP][GAP][GAP];
                  if (tmp > maxval) {
                    maxval = tmp;
                    maxi = i; maxj = j; maxk = 1;
                  }
            */
            /* value for Iy matrix, case target paired to gap (k=2) */
            /*coming from match, add (Xi, -; Yj-1, Yj) */
            mVal =
                M[i][j - 1] != 0
                    ? M[i][j - 1] + dsm[query_seq[i - 1]][GAP][target_seq[j - 2]][target_seq[j - 1]]
                    : -1;
            /*extend existing gap, add (-, -; Yj-1, Yj) */
            yVal = Iy[i][j - 1] != 0
                       ? Iy[i][j - 1] + dsm[GAP][GAP][target_seq[j - 2]][target_seq[j - 1]]
                       : -1;
            /* start new alignment - starts with gap LEGAL!? */
            /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */

            Iy[i][j] = max3(mVal, yVal, 0); /*removed option nVal */
            /* do not allow alignments ending in gap
                  tmp = Iy[i][j] + dsm[GAP][GAP][tseq[j-1]][GAP];
                  if (tmp > maxval) {
                    maxval = tmp;
                    maxi = i; maxj = j; maxk = 2;
                  }
            */
        }
    }
#ifdef VERBOSE
    printf("found maxval %d on pos %d/%d in mat %d\n", maxval, maxi, maxj,
           maxk); /* pos are 1-based */
    printf("print F[0] matrix:\n");
    printMat(M, m + 1, n + 1, query_seq, target_seq);
    printf("print F[1] matrix:\n");
    printMat(Ix, m + 1, n + 1, query_seq, target_seq);
    printf("print F[2] matrix:\n");
    printMat(Iy, m + 1, n + 1, query_seq, target_seq);
#endif

    /*backtrack*/
    auto tmp = (int)(1.5 * config->tblen);

    auto i = maxi;
    auto j = maxj;
    auto k = maxk; /* 0-1-2 M Ix Iy - should always be 0 to begin with */
    auto l = 0;    /* alilen so far -used in backtrack */

    if (k != 0) {
        fprintf(stderr, "\nErr: Found highest value in one of the gap matrices (k=%d)!?\n\n", k);
    }

    tmp -= 2;
    while ((i > 0 && j > 0) && (M[i][j] > 0 || Ix[i][j] > 0 || Iy[i][j] > 0)) {
        if (l > tmp) {
            printf("Interaction longer than max, so the following is only the end of the full "
                   "alignment:\n");
            /*alt: stop here / reallocate? prevent creation of longer alignments in first place? */
            break;
        }
#ifdef DEBUG
#if VERBOSE > 1
        fprintf(stderr, "in backtrack step with i=%d, j=%d, k=%d\n", i, j, k);
#endif
#endif

        /*l++ in end of while instead of every sub? */

        /* find which cell gave score */
        if (k == 0) {
            /* highest score in M matrix, having a match! */
#ifdef DEBUG
#if VERBOSE > 1
            /*printf("check next if M[i][j] (%d) ==  M[i-1][j-1] (%d) +
             * dsm[qseq[i-2]][qseq[i-1]][tseq[j-2]][tseq[j-1]] (%d) \n", M[i][j], M[i-1][j-1],
             * dsm[qseq[i-2]][qseq[i-1]][tseq[j-2]][tseq[j-1]] );  //this line will segfault when at
             * i=1 OR j=1 */
            printf("check next if M[i][j] (%d) \n", M[i][j]);
            printf(" == dsm[GAP][qseq[i-1]][GAP][tseq[j-1]] (%d)\n",
                   dsm[GAP][query_seq[i - 1]][GAP][target_seq[j - 1]]);
#endif
#endif
            if (M[i][j] == dsm[GAP][query_seq[i - 1]][GAP][target_seq[j - 1]]) {
                /*      started new alignment */
#ifdef DEBUG
#if VERBOSE > 1
                printf("means we started alignment at M[%d][%d] \n", i, j);
#endif
#endif
                k = 3;
                /* need to print here already OR later check k and break */
                hit->ali_seq1[l] = index2nt(query_seq[--i]);
                hit->ali_seq2[l] = index2nt(target_seq[--j]);
                if (query_seq[i] + target_seq[j] == 3) {
                    /* AU or CG pair */
                    hit->ali_ia[l] = '|';
                } else if (query_seq[i] + target_seq[j] == 5) {
                    /*GU wobble pair */
                    hit->ali_ia[l] = '.';
                } else {
                    hit->ali_ia[l] = ' ';
                }
                l++;
                break;
            } else if (M[i][j] == M[i - 1][j - 1] + dsm[query_seq[i - 2]][query_seq[i - 1]]
                                                       [target_seq[j - 2]][target_seq[j - 1]]) {
                /* access is save here as matrix always build from subset only, so we need to come
                 * to an end with check before!? */
                /*      previous was also match */
#ifdef DEBUG
#if VERBOSE > 1
                printf("Not the case, so check if it's \n");
                printf(" ==  M[i-1][j-1] (%d)\n", M[i - 1][j - 1]);
                printf(
                    " + dsm[qseq[i-2]][qseq[i-1]][tseq[j-2]][tseq[j-1]] (%d) \n",
                    dsm[query_seq[i - 2]][query_seq[i - 1]][target_seq[j - 2]][target_seq[j - 1]]);
#endif
#endif
                k = 0;
            } else if (M[i][j] ==
                       Ix[i - 1][j - 1] +
                           dsm[query_seq[i - 2]][query_seq[i - 1]][GAP][target_seq[j - 1]]) {
                /*      coming from gap in seq2/target */
                k = 1;
            } else if (M[i][j] ==
                       Iy[i - 1][j - 1] +
                           dsm[GAP][query_seq[i - 1]][target_seq[j - 2]][target_seq[j - 1]]) {
                /*      coming from gap in seq1/query */
                k = 2;
            } else {
                printf("unexpected value in k=0.\n");
            }
            hit->ali_seq1[l] = index2nt(query_seq[--i]);
            hit->ali_seq2[l] = index2nt(target_seq[--j]);
            /*      printf("test to calc %c (%d or %d) + %c (%d or %d) = %d or %d\n",
             * hit.ali_seq1[l], hit.ali_seq1[l],qseq[i], hit.ali_seq2[l], hit.ali_seq2[l],tseq[j],
             * hit.ali_seq1[l]+hit.ali_seq2[l], qseq[i]+tseq[j] ); */
            if (query_seq[i] + target_seq[j] == 3) {
                hit->ali_ia[l] = '|';
            } else if (query_seq[i] + target_seq[j] == 5) {
                hit->ali_ia[l] = '.';
            } else {
                hit->ali_ia[l] = ' ';
            }
            l++;
#ifdef DEBUG
#if VERBOSE > 1
            printf("print pos %d / %d = %c/ %c\n", i, j, hit->ali_seq1[l - 1],
                   hit->ali_seq2[l - 1]);
#endif
#endif
        } else if (k == 1) {
            /* seq1(query) paired to a gap (in target) */
            if (Ix[i][j] ==
                M[i - 1][j] + dsm[query_seq[i - 2]][query_seq[i - 1]][target_seq[j - 1]][GAP]) {
                k = 0; /* open a new gap coming from match */
            } else if (Ix[i][j] ==
                       Ix[i - 1][j] + dsm[query_seq[i - 2]][query_seq[i - 1]][GAP][GAP]) {
                k = 1; /* extend existing gap */
            } else if (Ix[i][j] == dsm[GAP][query_seq[i - 1]][GAP][GAP]) {
                k = 3; /* start new alignment with gap; not possible, prevented by scoring... */
                fprintf(stderr, "\nErr: This alignment starts in a gap - not even an option!?\n");
                hit->ali_seq1[l] = index2nt(query_seq[--i]);
                hit->ali_ia[l] = ' ';
                hit->ali_seq2[l++] = '-';
                break;
            } else {
                printf("unexpected case in k=1 : %d\n", Ix[i][j]);
            }
            hit->ali_seq1[l] = index2nt(query_seq[--i]);
            hit->ali_ia[l] = ' ';
            hit->ali_seq2[l++] = '-';
        } else if (k == 2) {
            /* seq2(target) paired to a gap (in query) */
            if (Iy[i][j] ==
                M[i][j - 1] + dsm[query_seq[i - 1]][GAP][target_seq[j - 2]][target_seq[j - 1]]) {
                k = 0; /* open a new gap coming from match */
            } else if (Iy[i][j] ==
                       Iy[i][j - 1] + dsm[GAP][GAP][target_seq[j - 2]][target_seq[j - 1]]) {
                k = 2; /* extend existing gap */
            } else if (Iy[i][j] == dsm[GAP][GAP][GAP][target_seq[j - 1]]) {
                k = 3; /* start new alignment - with gap LEGAL!? */
                fprintf(stderr, "\nErr: This alignment starts in a gap - not even an option!?\n");
                hit->ali_seq1[l] = '-';
                hit->ali_ia[l] = ' ';
                hit->ali_seq2[l++] = index2nt(target_seq[--j]);
                break;
            } else {
                printf("unexpected case in k=2 : %d\n", Iy[i][j]);
            }
            hit->ali_seq1[l] = '-';
            hit->ali_ia[l] = ' ';
            hit->ali_seq2[l++] = index2nt(target_seq[--j]);
        } else {
            fprintf(stderr, "\nThis should really NEVER happen!\n");
        }
    }
    hit->ali_seq1[l] = '\0';
    hit->ali_ia[l] = '\0';
    hit->ali_seq2[l] = '\0';

    /* reverse sequences in the end*/
    /* printf("my length = %d ; altern. test = %d\n", l, strlen(ali_seq1)); */

    reverse_inplace(hit->ali_seq1.get(), l - 1);
    /*  std::reverse(str, &str[l]); */
    reverse_inplace(hit->ali_ia.get(), l - 1);
    reverse_inplace(hit->ali_seq2.get(), l - 1);

    hit->qbeg = i + 1;
    hit->qend = maxi;
    hit->tbeg = n + 1 - maxj;
    hit->tend = n - j;
    hit->max = maxval;
    /*
      printf("%d - %d\n", i+1, maxi); //alignment in seq1 from to
      printf("%s\n%s\n", ali_seq1, ali_seq2);
      printf("%lu - %lu (3' <-- 5')\n", n-j, n+1-maxj);  //n+1-(j+1)

      printf("Score2fakeE (not considering extpen): %.2f\n", (maxval-559.0)/(-100.0));
      printf("no of nucls in ia: %d + %d = %d\n", maxi-i , maxj-j, maxi-i + maxj-j);
    */
}


static int
RIs_linSpace(unsigned char* query_sequence,  /* query sequence - numeric representation */
             unsigned char* target_sequence, /* target sequence */
             int m,                          /* query seq length */
             int n,                          /* target seq length */
             short dsm[6][6][6][6],          /* scoring matrix -- TODO variable length!? */
             int extensionpenalty, /* as used in dsm, to calc Score2fakE -- now also a global */
             int threshold,        /* give out hits higher than that */
             const char* qname,    /* query name */
             const char* tname,    /* target name */
             const char* matname,  /* name of the scoring matrix */
             const config_st* config)
{
    const auto reference = reference_from_matrix(matname);


    /* create a hit-struct instead and print from main or other sub?*/
    int **M, **Ix, **Iy; /* matrices for alignment scores ending in different states */
    int maxi, maxj;
    int maxval, maxk, testmax; /* maxk not needed!? - max will never be found in gapped anyway! */

    unsigned char currentRow, lastRow; /* alternating 0;1 */
    int rowMax_score;
    int rowMax_pos;
    int tmpQbeg, tmpTend;

    short tmpQlen, tmpTlen;


    int hitcount = 0;

#if VERBOSE > 1
    int tmpi;
#endif

    /* should test if malloc successful! */
    MallocRAII<int> hits_score(n);
    MallocRAII<int> hits_pos(n);

    MallocRAII<unsigned char> tmpQseq(config->tblen);
    MallocRAII<unsigned char> tmpTseq(config->tblen);

    IA maxHit;
    testmax = (int)(1.5 * config->tblen);
    maxHit.ali_seq1 = MallocRAII<char>(testmax);
    maxHit.ali_seq2 = MallocRAII<char>(testmax);
    maxHit.ali_ia = MallocRAII<char>(testmax);


    M = allocIntMatrix(2, m + 1);  /* (Mis)Match */
    Ix = allocIntMatrix(2, m + 1); /*Insertion in x(=query), so x paired to gap (in y) */
    Iy = allocIntMatrix(2, m + 1); /*Insertion(=bulge) in y(=target) */
    maxi = maxj = maxk = 0;
    currentRow = lastRow = 1;

#if VERBOSE > 1
    printf("\t-");
    for (tmpi = 0; tmpi < m; tmpi++)
        printf("\t%c", index2nt(*(query_sequence + tmpi)));
#endif

    M[0][0] = Ix[0][0] = Iy[0][0] = 0;

    /*init first row (j=0) -- refers to "-" before first nt of target */

    /* explicitly handling of the boundary condition since i-2 is not defined for i = 1
       No need anymore!
       Ix[0][1] = 0; // MAX(0, dsm[GAP][qseq[0]][GAP][GAP]);      - do not bother with that
       //S(-,*;-,-) should ALWAYS be negative! OR dangling end which can not be incorp. in scoring
       scheme (not w/o looking further back) Iy[0][1] = M[0][1] = NEGINF; // not possible to occure
       -- TODO: in fact 0 works just as well
     */
    for (auto i = 1; i <= m; i++) {
        Iy[0][i] = M[0][i] = NEGINF; /* not possible before beginning of target seq */
        Ix[0][i] = 0;                /*MAX(0, dsm[qseq[i-2]][qseq[i-1]][GAP][GAP]); */
                                     /* do not bother, require to always have a match first!? */
    }

#if VERBOSE > 1 /* j = 0 */
    printf("\n-\t0");
    for (tmpi = 1; tmpi <= m; tmpi++)
        printf("\t-NI");
/* has been - get back to this if decide not to use NEGINF at all...
    printf("\n-");
    for (tmpi=0; tmpi<=m; tmpi++)
      printf("\t%d", (M[0][tmpi] == NEGINF ? -8 : M[0][tmpi]) ); // -8 as dummy for -inf
  }
 */
#endif

    /* init j=1 row */
    /*init first col (i=0) */
    Iy[1][0] = 0; /* MAX(0, dsm[GAP][GAP][GAP][tseq[n-1]]); */ /*n-1 is last nt in target; used to
                                                                  be 0 after reversion */
    Ix[1][0] = M[1][0] = NEGINF;

    M[1][1] = dsm[GAP][query_sequence[0]][GAP][target_sequence[n - 1]];
    /*    maxval = M[1][1] + MAX(0,dsm[qseq[0]][GAP][tseq[n-1]][GAP]); ---- why should we allow a
     * special treatment here!? */
    rowMax_score = M[1][1] + dsm[query_sequence[0]][GAP][target_sequence[n - 1]][GAP];
    rowMax_pos = 1;

    /* (1,1) cell can not be in Ix or Iy state. */
    Ix[1][1] = Iy[1][1] = NEGINF;

    for (auto i = 2; i <= m; i++) {
        /* value for M matrix, case we have a pair here (k=0) */
        M[1][i] = dsm[GAP][query_sequence[i - 1]][GAP][target_sequence[n - 1]];
        /*had been 3 possibilities before, removed case
         Ix[0][i-1] != 0 ? Ix[0][i-1] + dsm[qseq[i-2]][qseq[i-1]][GAP][tseq[n-1]] : -1,
         as above all Ix[0][i] are set to 0
        */

        /* max so far? */
        if ((testmax = M[1][i] + dsm[query_sequence[i - 1]][GAP][target_sequence[n - 1]][GAP]) >
            rowMax_score) {
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
              maxi = i; maxj = 1; maxk = 1;
            }
        */
        /* value for Iy matrix, case target sequence paired to gap (k=2) */
        /* not possible in this row */
        Iy[1][i] = NEGINF;
    }

    /* initialization of first rows and columns completed
       recursion to complete alignment with two residues follows
     */

#if VERBOSE > 1
    printf("\n%c\t-NI", (index2nt(*(target_sequence + n - 1))));
    for (tmpi = 1; tmpi <= m; tmpi++)
        printf("\t%d", M[1][tmpi]);
#endif

    /*no need to store rowMax_score in the first place, can access hits_score[] just as fast!? (keep
     * pointer to current) */
    /* Raw pointers for the sweep: reading through the owner means the compiler
       reloads its member across the calls in the loop. */
    int* const hs = hits_score.get();
    int* const hp = hits_pos.get();

    hs[0] = rowMax_score;
    hp[0] = rowMax_pos;
    maxval = rowMax_score;
    maxi = rowMax_pos;
    maxj = 1; /*can be set in init */

    for (auto j = 2; j <= n; j++) { /* new row */
        lastRow = currentRow;
        currentRow = j % 2;
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

        /* value for M matrix, case we have a pair here (k=0) */
        M[currentRow][1] = MAX(0, dsm[GAP][query_sequence[0]][GAP][target_sequence[n - j]]);
        /* starting a new alignment with (-, X1; -, Yj)  OR  0 if not possible */
        /* coming from gap in tseq NOT possible as we're looking at first position of query! */
        /*had been 3 possibilities before, removed case 'coming from gap in qseq (Iy-matrix), adding
           (-, X1; Yj-1, Yj)' Iy[lastRow][0] != 0 ? Iy[lastRow][0] +
           dsm[GAP][qseq[0]][tseq[n-j+1]][tseq[n-j]] : -1, as above all Iy[j][0] are set to 0
         */

        rowMax_score = M[currentRow][1] + dsm[query_sequence[0]][GAP][target_sequence[n - j]]
                                             [GAP]; /*stricter only if not 0 before!? */
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
            M[lastRow][1] +
                dsm[query_sequence[0]][GAP][target_sequence[n - j + 1]][target_sequence[n - j]],
            /* extending an existing gap - add (-, -; Yj-1, Yj) */
            Iy[lastRow][1] + dsm[GAP][GAP][target_sequence[n - j + 1]][target_sequence[n - j]]);
        /*  do NOT allow max other than match state!? -- however (-, -; *, -) is positive!?
            testmax = Iy[currentRow][1] + dsm[GAP][GAP][tseq[n-j]][GAP];
            if (testmax > maxval) {
              maxval = testmax;
              maxi = 1; maxj = j; maxk = 2;
            }
        */
        /* finished init of i=1 col */

        for (auto i = 2; i <= m; i++) { /* cols */ /* alt bed: *p1  */

            /* value for M matrix, case we have a pair here (k=0) */
            M[currentRow][i] = max4(
                /* coming from a match, add (Xi-1, Xi; Yi-1, Yi) */
                M[lastRow][i - 1] != 0
                    ? M[lastRow][i - 1] + dsm[query_sequence[i - 2]][query_sequence[i - 1]]
                                             [target_sequence[n - j + 1]][target_sequence[n - j]]
                    : -1,
                /* coming from gap in target, add (Xi-1, Xi; -, Yi) */
                Ix[lastRow][i - 1] +
                    dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][target_sequence[n - j]],
                /* coming from gap in query, add (-, Xi; Yi-1, Yi) */
                Iy[lastRow][i - 1] + dsm[GAP][query_sequence[i - 1]][target_sequence[n - j + 1]]
                                        [target_sequence[n - j]],
                /* starting a new alignment with this pair: (-, Xi; -, Yj) */
                dsm[GAP][query_sequence[i - 1]][GAP][target_sequence[n - j]]);
            if ((testmax = M[currentRow][i] +
                           dsm[query_sequence[i - 1]][GAP][target_sequence[n - j]][GAP]) >
                rowMax_score) {
                rowMax_score = testmax;
                rowMax_pos = i;
            }
            /* value for Ix matrix, case query paired to gap (k=1) */
            /* removed one option, namely: start new alignment that starts in gap, reflected by (-,
             * Xi; -, -) */
            Ix[currentRow][i] = MAX(
                /*coming from match, add (Xi-1, Xi; Yj, -) */
                M[currentRow][i - 1] +
                    dsm[query_sequence[i - 2]][query_sequence[i - 1]][target_sequence[n - j]][GAP],
                /*extend existing gap, add (Xi-1, Xi; -, -) */
                Ix[currentRow][i - 1] + dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][GAP]
                /* start new alignment - starts with gap LEGAL!?
                   dsm[GAP][qseq[i-1]][GAP][GAP] */
            );
            /* only M[][] can be max, point of backtrack...
                  testmax = Ix[currentRow][i] + dsm[qseq[i-1]][GAP][GAP][GAP];
                  if (testmax > maxval) {
                    maxval = testmax;
                    maxi = i; maxj = j; maxk = 1;
                  }
            */
            /* value for Iy matrix, case target paired to gap (k=2) */
            /* removed one option, namely: start new ali starting w/ gap, (-, -; Yj, -)  //
             * dsm[GAP][GAP][GAP][tseq[n-j]] */
            Iy[currentRow][i] = MAX(
                /*coming from match, add (Xi, -; Yj-1, Yj) */
                M[lastRow][i] + dsm[query_sequence[i - 1]][GAP][target_sequence[n - j + 1]]
                                   [target_sequence[n - j]],
                /*extend existing gap, add (-, -; Yj-1, Yj) */
                Iy[lastRow][i] + dsm[GAP][GAP][target_sequence[n - j + 1]][target_sequence[n - j]]
                /* start new alignment - starts with gap LEGAL!?
                   dsm[GAP][GAP][GAP][tseq[n-j]] */
            );
            /*  only M[][] can be max, point of backtrack...
                  testmax = Iy[currentRow][i] + dsm[GAP][GAP][tseq[n-j]][GAP];
                  if (testmax > maxval) {
                    maxval = testmax;
                    maxi = i; maxj = j; maxk = 2;
                  }
            */
        }

        hs[j - 1] = rowMax_score;
        hp[j - 1] = rowMax_pos;

        if (rowMax_score > maxval) {
            maxval = rowMax_score;
            maxi = rowMax_pos;
            maxj = j;
        }
#if VERBOSE > 1
        printf("\n%c\t-NI", (index2nt(*(target_sequence + n - j))));
        for (tmpi = 1; tmpi <= m; tmpi++)
            printf("\t%d", M[currentRow][tmpi]);
#endif

    } /*next row j */

#if VERBOSE > 1
    printf("\n");
#endif

    /*
      if checking for subopts && set energy threshold, do not guarantee to print best hit first, but
      only once! if checking for subopts and p[2-4], do not check best first.
    */
    if (!(config->doSubopt && (config->filter_e || config->printShort > 1))) {
#ifdef VERBOSE
        printf("found maxval %d on pos %d/%d in mat %d\n", maxval, maxi, maxj,
               maxk); /* pos are 1-based */
        printf("equals end pos in query: %d - start pos in target: %d\n", maxi, n + 1 - maxj);
#endif

        /*do backtrack for this one only! by recomputing whole matrix for this subsection */
        /* max going back config->tblen */
        tmpQbeg = maxi > config->tblen - 1 ? maxi - (config->tblen - 1) : 1;
        tmpTend = maxj > config->tblen - 1 ? maxj - (config->tblen - 1) : 1;
        tmpQlen = maxi - tmpQbeg + 1;
        tmpTlen = maxj - tmpTend + 1;

#ifdef VERBOSE
        printf("going to realign query %d - %d (len: %hd) vs. target: %d - %d (len: %hd)\n",
               tmpQbeg, maxi, tmpQlen, n + 1 - tmpTend, n + 1 - maxj, tmpTlen);
#endif

        for (auto i = 0; i < tmpQlen; i++) {
            tmpQseq[i] = query_sequence[tmpQbeg - 1 + i];
#ifdef VERBOSE
            printf("%c", index2nt(tmpQseq[i]));
#endif
        }
#ifdef VERBOSE
        printf("\n");
#endif
        for (auto i = 0; i < tmpTlen; i++) {
            tmpTseq[i] = target_sequence[n - tmpTend - i];
#ifdef VERBOSE
            printf("%c", index2nt(tmpTseq[i]));
#endif
        }
#ifdef VERBOSE
        printf("\n");
#endif

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

        RIs(tmpQseq.get(), tmpTseq.get(), tmpQlen, tmpTlen, dsm, &maxHit, config);

#ifdef VERBOSE
        printf("%d - %d\n", maxHit.qbeg, maxHit.qend); /*alignment in subseq1 from to */
        printf("%s\n%s\n%s\n", maxHit.ali_seq1.get(), maxHit.ali_ia.get(), maxHit.ali_seq2.get());
        printf("%d - %d (3' <-- 5')\n", maxHit.tend, maxHit.tbeg); /* n+1-(j+1) */
#endif
#ifdef VERBOSE
        printf("full final ali:\n");
#endif
        /*number of nt in ia to recalc Score2fakE - only tmp no need to store... */
        const auto energy = (maxHit.max + extensionpenalty * maxHit.nucleotide_count() - reference) / (-100.0);


        /** TODO :: use maxHit.max OR maxval==hits_score[maxj-1] in output !?!? **/
        if (energy <= config->max_energy) {
            if (config->printShort == 1) {
                printf("%d\t%d\t%d\t%d\t%.2f\n", tmpQbeg + maxHit.qbeg - 1,
                       tmpQbeg + maxHit.qend - 1, n - maxj + maxHit.tbeg, n - maxj + maxHit.tend,
                       energy);
                /* to be consistent with other output:
                        printf("%d\t%d\t%d\t%d\t%.2f\t%s\n", tmpQbeg+maxHit.qbeg-1,
                   tmpQbeg+maxHit.qend-1, n-maxj+maxHit.tend, n-maxj+maxHit.tbeg, energy,
                   maxHit.ali_ia);
                */
            } else if (config->printShort == 2) {
                printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\n", qname, tmpQbeg + maxHit.qbeg - 1,
                       tmpQbeg + maxHit.qend - 1, tname, n - maxj + maxHit.tbeg,
                       n - maxj + maxHit.tend, maxval, energy);
            } else if (config->printShort == 3) {
                hitcount += 1;
            } else {
                printf("Free energy [kcal/mol]: %.2f (%d)\n", energy, maxval);
                /*      printf("no of nucls in ia: %lu + %lu = %lu\n", maxHit.qend-maxHit.qbeg+1 ,
                 * maxHit.tend-maxHit.tbeg+1 ,
                 * maxHit.qend-maxHit.qbeg+1+maxHit.tend-maxHit.tbeg+1); */

                printf("%d - %d\n", tmpQbeg + maxHit.qbeg - 1,
                       tmpQbeg + maxHit.qend - 1); /*alignment in seq1 from to */
                printf("%s\n%s\n%s\n", maxHit.ali_seq1.get(), maxHit.ali_ia.get(),
                       maxHit.ali_seq2.get());
                printf("%d - %d (3' <-- 5')\n", n - maxj + maxHit.tend, n - maxj + maxHit.tbeg);
            }
        }
    }

    /* break here if -s not set */
    if (!config->doSubopt) {
        if (config->printShort == 3)
            printf("%s\t%s\t%d\n", qname, tname, hitcount);

        return 0;
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
            tmpQbeg = hits_pos[j] > config->tblen - 1 ? hits_pos[j] - (config->tblen - 1) : 1;
            tmpTend = j + 1 > config->tblen - 1 ? j + 1 - (config->tblen - 1) : 1;
            tmpQlen = hits_pos[j] - tmpQbeg + 1;
            tmpTlen = j + 1 - tmpTend + 1;

#ifdef VERBOSE
            printf(
                "\nfound another high hit %d at: end pos in query: %d , start pos in target: %d\n",
                hits_score[j], hits_pos[j], n - j);
            printf("going to realign query %d - %d (len: %hd) vs. target: %d - %d (len: %hd)\n",
                   tmpQbeg, hits_pos[j], tmpQlen, n + 1 - tmpTend, n - j, tmpTlen);
#endif

            for (auto i = 0; i < tmpQlen; i++) {
                tmpQseq[i] = query_sequence[tmpQbeg - 1 + i];
            }
            for (auto i = 0; i < tmpTlen; i++) {
                tmpTseq[i] = target_sequence[n - tmpTend - i];
            }

            RIs(tmpQseq.get(), tmpTseq.get(), tmpQlen, tmpTlen, dsm, &maxHit, config);

            /*number of nt in ia to recalc Score2fakE - only tmp no need to store... */
            const auto energy = (maxHit.max + extensionpenalty * maxHit.nucleotide_count() - reference) / (-100.0);


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


    return 0;
}
