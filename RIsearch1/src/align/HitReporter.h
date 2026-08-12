#pragma once

#include <cstdint>
#include <cstdio>

#include "InteractionAlignment.h"
#include "cli.h"
#include "energy.hpp"
#include "math/Matrix.h"
#include "memory/MallocRAII.hpp"
#include "optimization/QueryProfile.h"
#include "traceback.h"


/* Turns a winning DP cell into printed output.
 *
 * RIs_linSpace finds where the best alignments end but keeps only two rows, so
 * it cannot say which nucleotides paired. This clips a window of at most tblen
 * around a hit, re-aligns that window with a full matrix to recover the pairing,
 * converts the score to an energy and prints it.
 *
 * It owns the scratch that re-alignment needs, so those buffers are allocated
 * once per query/target pair and reused for every hit -- there may be hundreds
 * of thousands of them.
 */
class HitReporter {
public:
    HitReporter(const unsigned char* query, const unsigned char* target, std::uint32_t n,
                short dsm[6][6][6][6], const QueryProfile& profile, const config_st& config,
                const char* qname, const char* tname)
        : m_query(query), m_target(target), m_n(n), m_dsm(dsm), m_profile(profile),
          m_config(config), m_qname(qname), m_tname(tname),
          m_reference(reference_from_matrix(config.mat_name)),
          m_M(config.tblen + 1, config.tblen + 1), m_Ix(config.tblen + 1, config.tblen + 1),
          m_Iy(config.tblen + 1, config.tblen + 1), m_qseq(config.tblen), m_tseq(config.tblen),
          m_best(config.tblen + 1), m_hit(static_cast<int>(1.5 * config.tblen))
    {
    }

    /* pos_i, pos_j: the DP cell the hit ends at. score: what the sweep recorded
       for it, which is what gets printed -- see the note in print(). */
    void report(std::uint32_t pos_i, std::uint32_t pos_j, int score, bool is_suboptimal)
    {
        const auto w = extract(pos_i, pos_j);

        RIs(m_qseq.get(), m_tseq.get(), w.qlen, w.tlen, &m_hit, m_config, m_M.get(), m_Ix.get(),
            m_Iy.get(), m_profile, w.qbeg - 1, m_best.get());

        const auto energy =
            (m_hit.max + m_config.extension_penalty * m_hit.nucleotide_count() - m_reference) /
            (-100.0);

        if (energy <= m_config.max_energy) {
            print(w, pos_j, score, energy, is_suboptimal);
        }
    }

    int hitcount() const
    {
        return m_hitcount;
    }

private:
    struct Window {
        std::uint32_t qbeg; /* 1-based query position the window starts at */
        std::uint32_t qlen;
        std::uint32_t tlen;
    };

    /* Clip to at most tblen back from the hit and copy both slices into scratch. */
    Window extract(std::uint32_t pos_i, std::uint32_t pos_j)
    {
        const auto qbeg = pos_i > m_config.tblen - 1 ? pos_i - (m_config.tblen - 1) : 1;
        const auto tbeg = pos_j > m_config.tblen - 1 ? pos_j - (m_config.tblen - 1) : 1;
        const Window w{qbeg, pos_i - qbeg + 1, pos_j - tbeg + 1};

        for (auto i = 0u; i < w.qlen; i++) {
            m_qseq[i] = m_query[qbeg - 1 + i];
        }
        /* Reversed: RIs walks the target forward, and DP row order runs 3'->5'. */
        for (auto i = 0u; i < w.tlen; i++) {
            m_tseq[i] = m_target[m_n - tbeg - i];
        }
        return w;
    }

    void print(const Window& w, int pos_j, int score, double energy, bool is_suboptimal)
    {
        const auto qb = w.qbeg + m_hit.qbeg - 1;
        const auto qe = w.qbeg + m_hit.qend - 1;
        const auto tb = m_n - pos_j + m_hit.tbeg;
        const auto te = m_n - pos_j + m_hit.tend;

        switch (m_config.printShort) {
        case 1:
            /* The two callers have always disagreed here: the best hit prints
               five fields as (tbeg, tend); a suboptimal prints six, as
               (tend, tbeg) plus the interaction string. Preserved so that
               extracting this function cannot move any output. */
            if (is_suboptimal) {
                printf("%d\t%d\t%d\t%d\t%.2f\t%s\n", qb, qe, te, tb, energy, m_hit.ali_ia.get());
            } else {
                printf("%d\t%d\t%d\t%d\t%.2f\n", qb, qe, tb, te, energy);
            }
            break;

        case 2:
            printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\n", m_qname, qb, qe, m_tname, tb, te, score,
                   energy);
            break;

        case 3:
            m_hitcount += 1;
            break;

        default:
            /* score comes from the linear-space sweep, energy from the window
               re-alignment above. They can disagree when the alignment is longer
               than tblen, which is the open question the two TODOs marked. */
            printf("Free energy [kcal/mol]: %.2f (%d)\n", energy, score);
            printf("%d - %d\n", qb, qe); /* alignment in seq1 from to */
            printf("%s\n%s\n%s\n", m_hit.ali_seq1.get(), m_hit.ali_ia.get(), m_hit.ali_seq2.get());
            printf("%d - %d (3' <-- 5')\n", te, tb);
            break;
        }
    }

    const unsigned char* m_query;
    const unsigned char* m_target;
    std::uint32_t m_n;
    short (*m_dsm)[6][6][6];
    const QueryProfile& m_profile;
    const config_st& m_config;
    const char* m_qname;
    const char* m_tname;
    float m_reference;

    /* Scratch, reused across every hit. */
    MatrixInt m_M, m_Ix, m_Iy;
    MallocRAII<unsigned char> m_qseq, m_tseq;
    /* Best M + close per query column; transpose_best_cell reads it. */
    MallocRAII<int> m_best;
    IA m_hit;

    int m_hitcount = 0;
};