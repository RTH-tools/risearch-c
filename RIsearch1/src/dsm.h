#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int getMat(const char *matname, short *bA_nu, int extPen, int transpose_matrix)
{
    short *bA_bas;
    extern short dsm_extend[6][6][6][6];
    if (!strcmp(matname, "t04")) {
        extern short dsm_t04[6][6][6][6];
        bA_bas = &dsm_t04[0][0][0][0];
    } else if (!strcmp(matname, "t99")) {
        extern short dsm_t99[6][6][6][6];
        bA_bas = &dsm_t99[0][0][0][0];
    } else if (!strcmp(matname, "su95")) {
        extern short dsm_su95_rev_wGU_pos[6][6][6][6];
        bA_bas = &dsm_su95_rev_wGU_pos[0][0][0][0];
    } else if (!strcmp(matname, "slh04_noGU")) {
        extern short dsm_slh04_woGU_pos[6][6][6][6];
        bA_bas = &dsm_slh04_woGU_pos[0][0][0][0];
    } else if (!strcmp(matname, "su95_noGU")) {
        extern short dsm_su95_rev_woGU_pos[6][6][6][6];
        bA_bas = &dsm_su95_rev_woGU_pos[0][0][0][0];
    } else {
        fprintf(stderr,
                "Undefined matrix (%s), -m needs to be set to either t99 or t04 for RNA-RNA "
                "interaction, su95 or su95_noGU for RNA-DNA interaction or slh04_noGU for DNA "
                "interaction\n",
                matname);
        exit(1);
    }
    short *bA_ext = &dsm_extend[0][0][0][0];

    if (transpose_matrix) {
        for (auto i = 0u; i < 6u; ++i) {
            for (auto j = 0u; j < 6u; ++j) {
                for (auto k = 0u; k < 6u; ++k) {
                    for (auto l = 0u; l < 6u; ++l) {
                        // Source: [i][j][k][l]
                        const auto src_idx = (i * 6 * 6 * 6) + (j * 6 * 6) + (k * 6) + l;

                        // Dest: [k][l][i][j] - swapping i,j with k,l
                        const auto dest_idx = (k * 6 * 6 * 6) + (l * 6 * 6) + (i * 6) + j;

                        *(bA_nu + dest_idx) = *(bA_bas + src_idx) - extPen * *(bA_ext + src_idx);
                    }
                }
            }
        }
    } else {
        /* create dsm from   dsm_base - d * dsm_extend   */
        for (auto i = 0u; i < 1296u; i++) {
            *(bA_nu + i) = *(bA_bas + i) - extPen * *(bA_ext + i); /* bA_nu[i] =  ... also works */
        }
    }

    return 0;
}