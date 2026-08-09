#pragma once

#include <cstring>

#include "InteractionAlignment.h"
#include "cli.h"
#include "debug_print.h"
#include "string_util.h"

static void RIs(const unsigned char* query_seq,  /* query sequence - numeric representation */
                const unsigned char* target_seq, /* target sequence - reversed */
                int m,                           /* query seq length */
                int n,                           /* target seq length */
                short dsm[6][6][6][6],           /* scoring matrix */
                IA* hit,                         /* pointer to struct, fill results */
                const config_st* config,
                int** M,
                int** Ix,
                int** Iy
                )
{
    int maxk = 0;               /* maxk not needed!? k itself could even be char... */
    int mVal, xVal, yVal, nVal; /* values coming from M, Ix, Iy, or starting a NEW alignment */

    RunningMax running_max{};

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

    // avoid set_if_better to preserve original behavior
    running_max.score =
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

        running_max.set_if_better(M[i][1] + dsm[query_seq[i - 1]][GAP][target_seq[0]][GAP], i, 1);


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

        running_max.set_if_better(M[1][j] + dsm[query_seq[0]][GAP][target_seq[j - 1]][GAP], 1, j);

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
            running_max.set_if_better(M[i][j] + dsm[query_seq[i - 1]][GAP][target_seq[j - 1]][GAP],
                                      i, j);

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
    printf("found maxval %d on pos %d/%d in mat %d\n", running_max.score, running_max.pos_i, running_max.pos_j,
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

    auto i = running_max.pos_i;
    auto j = running_max.pos_j;
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
    hit->qend = running_max.pos_i;
    hit->tbeg = n + 1 - running_max.pos_j;
    hit->tend = n - j;
    hit->max = running_max.score;
    /*
      printf("%d - %d\n", i+1, maxi); //alignment in seq1 from to
      printf("%s\n%s\n", ali_seq1, ali_seq2);
      printf("%lu - %lu (3' <-- 5')\n", n-j, n+1-maxj);  //n+1-(j+1)

      printf("Score2fakeE (not considering extpen): %.2f\n", (maxval-559.0)/(-100.0));
      printf("no of nucls in ia: %d + %d = %d\n", maxi-i , maxj-j, maxi-i + maxj-j);
    */
}

