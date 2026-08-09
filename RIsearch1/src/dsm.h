#pragma once

#include <cstdio>
#include <cstring>
// ReSharper disable once CppUnusedIncludeDirective
#include <cstdlib>


extern const short dsm_t99[6][6][6][6];
extern const short dsm_t04[6][6][6][6];
extern const short dsm_su95_rev_wGU_pos[6][6][6][6];
extern const short dsm_su95_rev_woGU_pos[6][6][6][6];
extern const short dsm_slh04_woGU_pos[6][6][6][6];

extern const short dsm_extend[6][6][6][6];

constexpr unsigned int DSM_SIDE = 6; /* alphabet size: A C G T/U N - */

struct NamedMatrix {
    const char* name;
    const short* table;
    const char* interaction;
};

static constexpr NamedMatrix DSM_MATRICES[]{
    {"t04", &dsm_t04[0][0][0][0], "RNA-RNA"},
    {"t99", &dsm_t99[0][0][0][0], "RNA-RNA"},
    {"su95", &dsm_su95_rev_wGU_pos[0][0][0][0], "RNA-DNA"},
    {"su95_noGU", &dsm_su95_rev_woGU_pos[0][0][0][0], "RNA-DNA"},
    {"slh04_noGU", &dsm_slh04_woGU_pos[0][0][0][0], "DNA-DNA"},
};

static const short* find_dsm(const char* matname)
{
    for (const auto& matrix : DSM_MATRICES) {
        if (!strcmp(matname, matrix.name)) {
            return matrix.table;
        }
    }
    return nullptr;
}

static void getMat(const char* matname, short* dsm, short extension_penalty, int transpose)
{
    const short* base = find_dsm(matname);
    if (!base) {
        fprintf(stderr,
          "Undefined matrix (%s), -m needs to be set to either t99 or t04 for RNA-RNA "
          "interaction, su95 or su95_noGU for RNA-DNA interaction or slh04_noGU for DNA "
          "interaction\n",
          matname);
        exit(1);
    }

    const short* extend = &dsm_extend[0][0][0][0];
    for (auto q_prev = 0u; q_prev < DSM_SIDE; q_prev++) {
        for (auto q_cur = 0u; q_cur < DSM_SIDE; q_cur++) {
            for (auto t_prev = 0u; t_prev < DSM_SIDE; t_prev++) {
                for (auto t_cur = 0u; t_cur < DSM_SIDE; t_cur++) {

                    const auto src =
                        ((q_prev * DSM_SIDE + q_cur) * DSM_SIDE + t_prev) * DSM_SIDE + t_cur;

                    const auto dst =
                        transpose
                            ? ((t_prev * DSM_SIDE + t_cur) * DSM_SIDE + q_prev) * DSM_SIDE + q_cur
                            : src;

                    dsm[dst] = base[src] - extension_penalty * extend[src];
                }
            }
        }
    }

}