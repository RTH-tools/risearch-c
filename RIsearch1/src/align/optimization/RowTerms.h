#pragma once

/* Every dsm term the row loop needs, for one (t_prev, t_cur) pair and one query
   position. Both target nucleotides are constant within a row and the query
   never changes, so these are computed once per alignment rather than per cell. */
struct RowTerms {
    int m_from_m;   /* dsm[q_prev][q_cur][t_prev][t_cur] -- extend a pair       */
    int m_from_ix;  /* dsm[q_prev][q_cur][GAP][t_cur]    -- close a query bulge */
    int m_from_iy;  /* dsm[GAP][q_cur][t_prev][t_cur]    -- close a target bulge*/
    int m_open;     /* dsm[GAP][q_cur][GAP][t_cur]       -- open on this pair   */
    int close;      /* dsm[q_cur][GAP][t_cur][GAP]       -- terminate after it  */
    int ix_from_m;  /* dsm[q_prev][q_cur][t_cur][GAP]    -- open a query bulge  */
    int iy_from_m;  /* dsm[q_cur][GAP][t_prev][t_cur]    -- open a target bulge */
};
