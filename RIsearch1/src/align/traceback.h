#pragma once

#include <cstring>

#include "InteractionAlignment.h"
#include "RunningMax.h"
#include "cli.h"
#include "memory/ByteBuffer.hpp"
#include "nucleotide.h" /* GAP, NEGINF */
#include "operations.h"
#include "optimization/QueryProfile.h"
#include "string_util.h"

enum class TraceState { TRACE_M = 0, TRACE_IX = 1, TRACE_IY = 2, TRACE_DONE = 3 };

/* '|' for a Watson-Crick pair, '.' for a G:U wobble, ' ' for a mismatch.
   The indices are chosen so a pair sums to 3 and a wobble to 5. */
static char pair_symbol(unsigned char q, unsigned char t)
{
    if (q + t == 3) {
        return '|';
    }
    if (q + t == 5) {
        return '.';
    }
    return ' ';
}

static void emit_pair(IA* hit, int l, const unsigned char* q, const unsigned char* t, int i, int j)
{
    const auto qn = q[i];
    const auto tn = t[j];
    hit->ali_seq1[l] = index2nt(qn);
    hit->ali_ia[l] = pair_symbol(qn, tn);
    hit->ali_seq2[l] = index2nt(tn);
}

static void emit_query_bulge(IA* hit, int l, const unsigned char* q, int i)
{
    hit->ali_seq1[l] = index2nt(q[i]);
    hit->ali_ia[l] = ' ';
    hit->ali_seq2[l] = '-';
}

static void emit_target_bulge(IA* hit, int l, const unsigned char* t, int j)
{
    hit->ali_seq1[l] = '-';
    hit->ali_ia[l] = ' ';
    hit->ali_seq2[l] = index2nt(t[j]);
}


static void RIs(const unsigned char* query_seq,  /* query sequence - numeric representation */
                const unsigned char* target_seq, /* target sequence - reversed */
                int m,                           /* query seq length */
                int n,                           /* target seq length */
                IA* hit,                         /* pointer to struct, fill results */
                const config_st& config, int** M, int** Ix, int** Iy, const QueryProfile& profile,
                int q_offset)
{
    RunningMax running_max{};

    M[0][0] = Ix[0][0] = Iy[0][0] = 0;

    // First col assign
    for (auto i = 1; i <= m; i++) {
        Iy[i][0] = M[i][0] = NEGINF; /* not possible before beginning of target seq */
        Ix[i][0] = 0;
    }

    // First row assign
    for (auto j = 1; j <= n; j++) {
        Ix[0][j] = M[0][j] = NEGINF;
        Iy[0][j] = 0;
    }

    /*
     * The initialization of i=1 column and j=1 row have to be handled explicitly
     * since at this point we do not have two residues to use.
     *
     * Handle the (1,1) cell explicitly since the boundary recursion includes (i-2) or (j-2) cases.
     */

    // Use this  QueryProfile to fetch terms that ignore the previous target nt (GAP)
    // and relate to first target nt (target_seq[0])
    const auto T = profile.row(QueryProfile::context(GAP, target_seq[0]));
    const int* const ix_from_m_1 = T.ix_from_m;
    const int* const ix_ext = profile.ix_extend();

    M[1][1] = T.m_open[q_offset + 1];
    running_max.set(M[1][1] + T.close[q_offset + 1], 1, 1);

    /* (1,1) cell can not be in Ix or Iy state. */
    Ix[1][1] = Iy[1][1] = NEGINF;


    /* init j=1 col */
    for (auto i = 2; i <= m; i++) {
        const auto qp = q_offset + i;
        M[i][1] = T.m_open[qp];
        running_max.set_if_better(M[i][1] + T.close[qp], i, 1);

        Ix[i][1] = max3(M[i - 1][1] != 0 ? M[i - 1][1] + ix_from_m_1[qp] : -1,
                        Ix[i - 1][1] != 0 ? Ix[i - 1][1] + ix_ext[qp] : -1, 0);

        Iy[i][1] = NEGINF;
    }

    const auto qp_first = q_offset + 1;

    for (auto j = 2; j <= n; j++) {
        const auto context = QueryProfile::context(target_seq[j - 2], target_seq[j - 1]);
        const auto t = profile.row(context);

        M[1][j] = t.m_open[qp_first];
        running_max.set_if_better(M[1][j] + t.close[qp_first], 1, j);

        Ix[1][j] = NEGINF;

        Iy[1][j] = max3(M[1][j - 1] != 0 ? M[1][j - 1] + t.iy_from_m[qp_first] : -1,
                        Iy[1][j - 1] != 0 ? Iy[1][j - 1] + profile.row(context).iy_extend : -1, 0);
    }


    for (auto i = 2; i <= m; i++) {
        const auto qp = q_offset + i;

        for (auto j = 2; j <= n; j++) {
            const auto context = QueryProfile::context(target_seq[j - 2], target_seq[j - 1]);
            const auto t = profile.row(context);
            const int* const IXM = t.ix_from_m;
            const auto iy_ext = t.iy_extend;

            M[i][j] = max4(
                // continue from a pair
                M[i - 1][j - 1] != 0 ? M[i - 1][j - 1] + t.m_from_m[qp] : -1,
                // close a bulge in query
                Ix[i - 1][j - 1] != 0 ? Ix[i - 1][j - 1] + t.m_from_ix[qp] : -1,
                // close a bulge in target
                Iy[i - 1][j - 1] != 0 ? Iy[i - 1][j - 1] + t.m_from_iy[qp] : -1,
                // start fresh
                t.m_open[qp]);

            running_max.set_if_better(M[i][j] + t.close[qp], i, j);

            Ix[i][j] = max3(
                // Open a bulge in query
                M[i - 1][j] != 0 ? M[i - 1][j] + IXM[qp] : -1,
                // Continue a bulge in query
                Ix[i - 1][j] != 0 ? Ix[i - 1][j] + ix_ext[qp] : -1, 0);

            Iy[i][j] = max3(
                // Open a bulge in target
                M[i][j - 1] != 0 ? M[i][j - 1] + t.iy_from_m[qp] : -1,
                // Continue a bulge in target
                Iy[i][j - 1] != 0 ? Iy[i][j - 1] + iy_ext : -1, 0);
        }
    }

    /*backtrack*/
    // we write to capacity and capacity + 1 indices, so we need -2
    // as a buffer zone.
    const auto capacity = static_cast<int>(1.5 * config.tblen) - 2;

    auto i = running_max.pos_i;
    auto j = running_max.pos_j;
    auto k = TraceState::TRACE_M;
    auto l = 0; /* alilen so far -used in backtrack */

    while ((i > 0 && j > 0) && (M[i][j] > 0 || Ix[i][j] > 0 || Iy[i][j] > 0)) {
        if (l > capacity) {
            printf("Interaction longer than max, so the following is only the end of the full "
                   "alignment:\n");
            /*alt: stop here / reallocate? prevent creation of longer alignments in first place? */
            break;
        }


        /*l++ in end of while instead of every sub? */

        const auto qp = q_offset + i;


        switch (k) {
        case TraceState::TRACE_M: {
            if (const auto open_score =
                    profile.row(QueryProfile::context(GAP, target_seq[j - 1])).m_open[qp];
                M[i][j] == open_score) {
                k = TraceState::TRACE_DONE;
            } else {
                const auto& t =
                    profile.row(QueryProfile::context(target_seq[j - 2], target_seq[j - 1]));

                if (M[i][j] == M[i - 1][j - 1] + t.m_from_m[qp]) {
                    k = TraceState::TRACE_M;
                } else if (M[i][j] == Ix[i - 1][j - 1] + t.m_from_ix[qp]) {
                    k = TraceState::TRACE_IX;
                } else if (M[i][j] == Iy[i - 1][j - 1] + t.m_from_iy[qp]) {
                    k = TraceState::TRACE_IY;
                } else {
                    printf("unexpected value in k=0.\n");
                }
            }
            emit_pair(hit, l++, query_seq, target_seq, --i, --j);
            break;
        }
        case TraceState::TRACE_IX: {
            // in this case, t_prev doesn't matter so we put it as GAP
            // avoids bugs in the j == 1 case.
            /* seq1(query) paired to a gap (in target) */
            const auto gap_row = profile.row(QueryProfile::context(GAP, target_seq[j - 1]));
            if (Ix[i][j] == M[i - 1][j] + gap_row.ix_from_m[qp]) {
                k = TraceState::TRACE_M; /* open a new gap coming from match */
            } else if (Ix[i][j] == Ix[i - 1][j] + ix_ext[qp]) {
                k = TraceState::TRACE_IX; /* extend existing gap */
            } else if (Ix[i][j] == profile.row(QueryProfile::context(GAP, GAP)).m_open[qp]) {
                fprintf(stderr, "\nErr: This alignment starts in a gap - not even an option!?\n");
                k = TraceState::TRACE_DONE; /* start new alignment with gap; not possible, prevented
                                               by scoring... */
            } else {
                printf("unexpected case in k=1 : %d\n", Ix[i][j]);
            }
            emit_query_bulge(hit, l++, query_seq, --i);
            break;
        }
        case TraceState::TRACE_IY: {
            const auto context = QueryProfile::context(target_seq[j - 2], target_seq[j - 1]);

            /* seq2(target) paired to a gap (in query) */
            if (Iy[i][j] == M[i][j - 1] + profile.row(context).iy_from_m[qp]) {
                k = TraceState::TRACE_M; /* open a new gap coming from match */
            } else if (Iy[i][j] == Iy[i][j - 1] + profile.row(context).iy_extend) {
                k = TraceState::TRACE_IY; /* extend existing gap */
            } else if (Iy[i][j] ==
                       profile.row(QueryProfile::context(GAP, target_seq[j - 1])).iy_extend) {
                k = TraceState::TRACE_DONE;
                fprintf(stderr, "\nErr: This alignment starts in a gap - not even an option!?\n");
            } else {
                printf("unexpected case in k=2 : %d\n", Iy[i][j]);
            }
            emit_target_bulge(hit, l++, target_seq, --j);
            break;
        }
        default: {
            fprintf(stderr, "\nThis should really NEVER happen!\n");
            break;
        }
        }

        if (TraceState::TRACE_DONE == k) {
            break;
        }
    }
    hit->ali_seq1[l] = '\0';
    hit->ali_ia[l] = '\0';
    hit->ali_seq2[l] = '\0';

    /* reverse sequences in the end*/

    if (l > 0) { // fixes a pre-existing bug, access to invalid memory [-1]
        reverse_inplace(hit->ali_seq1.get(), l - 1);
        reverse_inplace(hit->ali_ia.get(), l - 1);
        reverse_inplace(hit->ali_seq2.get(), l - 1);
    }

    hit->qbeg = i + 1;
    hit->qend = running_max.pos_i;
    hit->tbeg = n + 1 - running_max.pos_j;
    hit->tend = n - j;
    hit->max = running_max.score;
}
