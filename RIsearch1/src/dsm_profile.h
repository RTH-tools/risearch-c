#pragma once

/* dsm specialised to one query. The two query indices of a lookup depend only
on the query position, so they are resolved once instead of per cell. Laid
   out [t_prev][t_cur][i] so fixing a target position gives a contiguous run
   over query positions. */
class QueryProfile {
public:
    QueryProfile(const unsigned char* query, int m, short dsm[6][6][6][6])
        : m_stride(m + 1), m_table(36 * (m + 1))
    {
        for (auto t_prev = 0; t_prev < 6; t_prev++)
            for (auto t_cur = 0; t_cur < 6; t_cur++)
                for (auto i = 2; i <= m; i++)
                    m_table[(t_prev * 6 + t_cur) * m_stride + i] =
                        dsm[query[i - 2]][query[i - 1]][t_prev][t_cur];
    }

    /* Stacking scores at every query position, for one target dinucleotide. */
    [[nodiscard]] const short* stack_row(int t_prev, int t_cur) const
    {
        return m_table.get() + (t_prev * 6 + t_cur) * m_stride;
    }

private:
    int m_stride;
    MallocRAII<short> m_table;
};