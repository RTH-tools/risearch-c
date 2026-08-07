#pragma once

static int getMat (const char *matname, short *bA_nu, int extPen, int transpose_matrix)
{
    short *bA_bas, *bA_ext;
    int i, j, k, l;
    int src_idx, dest_idx;
    extern short dsm_extend[6][6][6][6];
    if (!strcmp (matname, "t04"))
    {
        extern short dsm_t04[6][6][6][6];
        bA_bas = &dsm_t04[0][0][0][0];
    }
    else if (!strcmp (matname, "t99"))
    {
        extern short dsm_t99[6][6][6][6];
        bA_bas = &dsm_t99[0][0][0][0];
    }
    else if (!strcmp (matname, "su95"))
    {
        extern short dsm_su95_rev_wGU_pos[6][6][6][6];
        bA_bas = &dsm_su95_rev_wGU_pos[0][0][0][0];
    }
    else if (!strcmp (matname, "slh04_noGU"))
    {
        extern short dsm_slh04_woGU_pos[6][6][6][6];
        bA_bas = &dsm_slh04_woGU_pos[0][0][0][0];
    }
    else if (!strcmp (matname, "su95_noGU"))
    {
        extern short dsm_su95_rev_woGU_pos[6][6][6][6];
        bA_bas = &dsm_su95_rev_woGU_pos[0][0][0][0];
    }
    else
    {
        fprintf(stderr, "Undefined matrix (%s), -m needs to be set to either t99 or t04 for RNA-RNA interaction, su95 or su95_noGU for RNA-DNA interaction or slh04_noGU for DNA interaction\n", matname);
        exit (1);
    }
    bA_ext = &dsm_extend[0][0][0][0];

    if (transpose_matrix) {
        for (i = 0; i < 6; ++i) {
            for (j = 0; j < 6; ++j) {
                for (k = 0; k < 6; ++k) {
                    for (l = 0; l < 6; ++l) {
                        // Source: [i][j][k][l]
                        src_idx = (i * 6 * 6 * 6) + (j * 6 * 6) + (k * 6) + l;

                        // Dest: [k][l][i][j] - swapping i,j with k,l
                        dest_idx = (k * 6 * 6 * 6) + (l * 6 * 6) + (i * 6) + j;

                        *(bA_nu + dest_idx) = *(bA_bas + src_idx) - extPen * *(bA_ext + src_idx);
                    }
                }
            }
        }
    } else {
        /* create dsm from   dsm_base - d * dsm_extend   */
        for (i = 0; i < 1296; i++)
        {
            *(bA_nu + i) = *(bA_bas + i) - extPen * *(bA_ext + i);	/* bA_nu[i] =  ... also works */
        }
    }

    return 0;
}