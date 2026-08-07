// Tests for getMat: assembling the scoring matrix from a named table, the
// extension penalty, and the -R transpose.

#include <gtest/gtest.h>

#include "dsm.h"

// The raw tables in dsm.c, before any penalty is applied.
extern short dsm_su95_rev_woGU_pos[6][6][6][6];
extern short dsm_extend[6][6][6][6];

TEST(Dsm, WithoutPenaltyReproducesTheRawTable)
{
  short out[6][6][6][6];
  const int extPen = 0;
  const int transpose_matrix = 0;
  char name[] = "su95_noGU";

  ASSERT_EQ(getMat(name, &out[0][0][0][0], extPen, transpose_matrix), 0);

  EXPECT_EQ(out[0][1][2][3], dsm_su95_rev_woGU_pos[0][1][2][3]);
  EXPECT_EQ(out[3][3][0][0], dsm_su95_rev_woGU_pos[3][3][0][0]);
  EXPECT_EQ(out[5][5][5][5], dsm_su95_rev_woGU_pos[5][5][5][5]);
}

TEST(Dsm, SubtractsThePenaltyTimesTheExtensionTable)
{
  short out[6][6][6][6];
  char name[] = "su95_noGU";
  const int extPen = 30;
  const int transpose_matrix = 0;

  ASSERT_EQ(getMat(name, &out[0][0][0][0], extPen, transpose_matrix), 0);

  EXPECT_EQ(out[0][1][2][3],
            dsm_su95_rev_woGU_pos[0][1][2][3] - 30 * dsm_extend[0][1][2][3]);
  EXPECT_EQ(out[3][3][0][0],
            dsm_su95_rev_woGU_pos[3][3][0][0] - 30 * dsm_extend[3][3][0][0]);
}

TEST(Dsm, TransposeSwapsTheQueryAndTargetHalvesOfTheIndex)
{
  // -R scores the duplex with the strands exchanged: the entry for query (i,j)
  // against target (k,l) must land at target (k,l) against query (i,j).
  short out[6][6][6][6];
  const int extPen = 0;
  const int transpose_matrix = 1;
  char name[] = "su95_noGU";

  ASSERT_EQ(getMat(name, &out[0][0][0][0], extPen, transpose_matrix), 0);

  EXPECT_EQ(out[2][3][0][1], dsm_su95_rev_woGU_pos[0][1][2][3]);
  EXPECT_EQ(out[0][0][3][3], dsm_su95_rev_woGU_pos[3][3][0][0]);
}

TEST(Dsm, LoadsEveryDocumentedMatrixName)
{
  short out[6][6][6][6];
  const int extPen = 0;
  const int transpose_matrix = 0;

  char t04[] = "t04";
  char t99[] = "t99";
  char su95[] = "su95";
  char su95_noGU[] = "su95_noGU";
  char slh04_noGU[] = "slh04_noGU";

  EXPECT_EQ(getMat(t04, &out[0][0][0][0], extPen, transpose_matrix), 0);
  EXPECT_EQ(getMat(t99, &out[0][0][0][0], extPen, transpose_matrix), 0);
  EXPECT_EQ(getMat(su95, &out[0][0][0][0], extPen, transpose_matrix), 0);
  EXPECT_EQ(getMat(su95_noGU, &out[0][0][0][0], extPen, transpose_matrix), 0);
  EXPECT_EQ(getMat(slh04_noGU, &out[0][0][0][0], extPen, transpose_matrix), 0);
}
