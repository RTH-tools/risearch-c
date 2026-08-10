#pragma once

#include <cstdint>

#include "RowTerms.h"
#include "dsm.h"
#include "memory/MallocRAII.hpp"
#include "nucleotide.h" /* GAP */

/* Scores for one query against every target context, resolved up front.
 *
 * Inside a row both target nucleotides are fixed and the query never changes,
 * so each dsm lookup collapses from four indices to (q_prev, q_cur). There are
 * only DSM_SIDE^2 target contexts, so every term is resolved once per alignment
 * and selected per row with a base pointer -- the hot loop then does linear
 * loads instead of four-level gathers. */
class QueryProfile {
public:
    QueryProfile(const unsigned char* query_sequence, std::uint32_t m, short dsm[6][6][6][6])
        : m_stride(m + 1), m_terms(DSM_SIDE * DSM_SIDE * (m + 1)),
          m_ix_from_m(DSM_SIDE * DSM_SIDE * (m + 1)), m_ix_extend(m + 1)
    {
        /* No target dependence: a query bulge over a gap on both sides. */
        for (auto i = 2u; i <= m; i++) {
            m_ix_extend[i] = dsm[query_sequence[i - 2]][query_sequence[i - 1]][GAP][GAP];
        }

        for (auto t_prev = 0u; t_prev < DSM_SIDE; t_prev++) {
            for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {
                const auto ctx = context(t_prev, t_cur);

                /* No query dependence: a target bulge over a gap. */
                m_iy_extend[ctx] = dsm[GAP][GAP][t_prev][t_cur];

                RowTerms* const terms = m_terms.get() + ctx * m_stride;
                for (auto i = 1u; i <= m; i++) {
                    const auto q_cur = query_sequence[i - 1];
                    /* Column 1 has no predecessor, so the q_prev terms are never
                       read there; GAP is a placeholder. */
                    const auto q_prev =
                        i >= 2 ? query_sequence[i - 2] : static_cast<unsigned char>(GAP);

                    terms[i].m_from_m = dsm[q_prev][q_cur][t_prev][t_cur];
                    terms[i].m_from_ix = dsm[q_prev][q_cur][GAP][t_cur];
                    terms[i].m_from_iy = dsm[GAP][q_cur][t_prev][t_cur];
                    terms[i].m_open = dsm[GAP][q_cur][GAP][t_cur];
                    terms[i].close = dsm[q_cur][GAP][t_cur][GAP];
                    m_ix_from_m[ctx * m_stride + i] = dsm[q_prev][q_cur][t_cur][GAP];
                    terms[i].iy_from_m = dsm[q_cur][GAP][t_prev][t_cur];
                }
            }
        }
    }

    static unsigned context(unsigned t_prev, unsigned t_cur)
    {
        return t_prev * DSM_SIDE + t_cur;
    }

    const RowTerms* row(unsigned ctx) const
    {
        return m_terms.get() + ctx * m_stride;
    }
    /* Its own contiguous run: the Ix pass reads only this term, so pulling a
       whole RowTerms to use one field of it wastes most of each cache line. */
    const int* ix_from_m(unsigned ctx) const
    {
        return m_ix_from_m.get() + ctx * m_stride;
    }

    const int* ix_extend() const
    {
        return m_ix_extend.get();
    }
    int iy_extend(unsigned ctx) const
    {
        return m_iy_extend[ctx];
    }

private:
    std::uint32_t m_stride;
    MallocRAII<RowTerms> m_terms;
    MallocRAII<int> m_ix_from_m;
    MallocRAII<int> m_ix_extend;
    int m_iy_extend[DSM_SIDE * DSM_SIDE]{};
};