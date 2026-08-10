#pragma once

#include <cstdint>

#include "dsm.h"
#include "memory/MallocRAII.hpp"
#include "nucleotide.h" /* GAP */

/* Scores for one query against every target context, resolved up front.
 *
 * Inside a row both target nucleotides are fixed and the query never changes,
 * so each dsm lookup collapses from four indices to (q_prev, q_cur). There are
 * only DSM_SIDE^2 target contexts, so every term is resolved once per alignment
 * and selected per row with a base pointer -- the hot loop then does linear
 * loads instead of four-level gathers.
 *
 * Each term gets its own run rather than the seven being interleaved per query
 * position. A pass that reads one term walks it without pulling the other six
 * through the cache, and several query positions of one term sit adjacent,
 * which is what a wide load needs. */
class QueryProfile {
public:
    /* One row's terms. The context is resolved once and a query position then
       indexes each run directly: t.m_from_m[i] is the term for column i. */
    struct RowView {
        const int* m_from_m;  /* dsm[q_prev][q_cur][t_prev][t_cur] -- extend a pair        */
        const int* m_from_ix; /* dsm[q_prev][q_cur][GAP][t_cur]    -- close a query bulge  */
        const int* m_from_iy; /* dsm[GAP][q_cur][t_prev][t_cur]    -- close a target bulge */
        const int* m_open;    /* dsm[GAP][q_cur][GAP][t_cur]       -- open on this pair    */
        const int* close;     /* dsm[q_cur][GAP][t_cur][GAP]       -- terminate after it   */
        const int* iy_from_m; /* dsm[q_cur][GAP][t_prev][t_cur]    -- open a target bulge  */
        const int* ix_from_m; /* dsm[q_prev][q_cur][t_cur][GAP]    -- open a query bulge   */
        const int* ix_extend; /* dsm[q_prev][q_cur][GAP][GAP]      -- by query position    */
        /* ix_from_m with ix_prefix taken out, and the running total to put it
           back -- see the constructor. */
        const int* ix_from_m_scan;
        const int* ix_prefix; /* by query position only */
        int iy_extend;        /* dsm[GAP][GAP][t_prev][t_cur]      -- one value per row    */
    };

    QueryProfile(const unsigned char* query_sequence, std::uint32_t m, short dsm[6][6][6][6])
        : m_stride(m + 1), m_m_from_m(kContexts * (m + 1)), m_m_from_ix(kContexts * (m + 1)),
          m_m_from_iy(kContexts * (m + 1)), m_m_open(kContexts * (m + 1)),
          m_close(kContexts * (m + 1)), m_iy_from_m(kContexts * (m + 1)),
          m_ix_from_m(kContexts * (m + 1)), m_ix_extend(m + 1),
          m_ix_from_m_scan(kContexts * (m + 1)), m_ix_prefix(m + 1)
    {
        /* No target dependence: a query bulge over a gap on both sides. */
        for (auto i = 2u; i <= m; i++) {
            m_ix_extend[i] = dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][GAP];
        }

        /* Ix[i] = max(M[i-1] + ix_from_m[i], Ix[i-1] + ix_extend[i]) is a running
           max in which each carried candidate also picks up every ix_extend
           between where it started and i. Subtracting the running total of
           ix_extend from each candidate, and adding it back once at the end,
           leaves a plain running max -- one that several query positions can
           resolve together. ix_prefix is that running total, fixed for the
           whole query. */
        m_ix_prefix[0] = 0;
        m_ix_prefix[1] = 0;
        for (auto i = 2u; i <= m; i++) {
            m_ix_prefix[i] = m_ix_prefix[i - 1] + m_ix_extend[i];
        }

        for (auto t_prev = 0u; t_prev < DSM_SIDE; t_prev++) {
            for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {
                const auto ctx = context(t_prev, t_cur);
                const auto off = ctx * m_stride;

                /* No query dependence: a target bulge over a gap. */
                m_iy_extend[ctx] = dsm[GAP][GAP][t_prev][t_cur];

                for (auto i = 1u; i <= m; i++) {
                    const auto q_cur = query_sequence[i - 1];
                    /* Column 1 has no predecessor, so the q_prev terms are never
                       read there; GAP is a placeholder. */
                    const auto q_prev =
                        i >= 2 ? query_sequence[i - 2] : static_cast<unsigned char>(GAP);

                    m_m_from_m[off + i] = dsm[q_prev][q_cur][t_prev][t_cur];
                    m_m_from_ix[off + i] = dsm[q_prev][q_cur][GAP][t_cur];
                    m_m_from_iy[off + i] = dsm[GAP][q_cur][t_prev][t_cur];
                    m_m_open[off + i] = dsm[GAP][q_cur][GAP][t_cur];
                    m_close[off + i] = dsm[q_cur][GAP][t_cur][GAP];
                    m_iy_from_m[off + i] = dsm[q_cur][GAP][t_prev][t_cur];
                    m_ix_from_m[off + i] = dsm[q_prev][q_cur][t_cur][GAP];
                    m_ix_from_m_scan[off + i] =
                        m_ix_from_m[off + i] - m_ix_prefix[i];
                }
            }
        }
    }

    static unsigned context(unsigned t_prev, unsigned t_cur)
    {
        return t_prev * DSM_SIDE + t_cur;
    }

    std::uint32_t m() const
    {
        return m_stride - 1;
    }

    /* Also carried in RowView; kept here for callers that have no row in hand,
       which is legitimate because this term has no target dependence. */
    const int* ix_extend() const
    {
        return m_ix_extend.get();
    }

    RowView row(unsigned ctx) const
    {
        const auto off = ctx * m_stride;
        return {m_m_from_m.get() + off,  m_m_from_ix.get() + off, m_m_from_iy.get() + off,
                m_m_open.get() + off,    m_close.get() + off,     m_iy_from_m.get() + off,
                m_ix_from_m.get() + off, m_ix_extend.get(),
                m_ix_from_m_scan.get() + off,
                m_ix_prefix.get(),       m_iy_extend[ctx]};
    }

private:
    static constexpr unsigned kContexts = DSM_SIDE * DSM_SIDE;

    std::uint32_t m_stride;
    MallocRAII<int> m_m_from_m;
    MallocRAII<int> m_m_from_ix;
    MallocRAII<int> m_m_from_iy;
    MallocRAII<int> m_m_open;
    MallocRAII<int> m_close;
    MallocRAII<int> m_iy_from_m;
    MallocRAII<int> m_ix_from_m;
    MallocRAII<int> m_ix_extend;
    MallocRAII<int> m_ix_from_m_scan;
    MallocRAII<int> m_ix_prefix;
    int m_iy_extend[kContexts]{};
};
