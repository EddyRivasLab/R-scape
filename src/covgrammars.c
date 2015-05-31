/* covgrammars.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covgrammars.h"

GMX *
GMX_Create(int L)
{
  GMX *gmx = NULL;
  int  pos = 1;
  int  j, d;
  int  status;

  ESL_ALLOC(gmx, sizeof(GMX));
  gmx->L = L;
  
  ESL_ALLOC(gmx->dp,    sizeof(SCVAL *) * (L+1));
  ESL_ALLOC(gmx->dp[0], sizeof(SCVAL  ) * ((L+2)*(L+1))/2);

  for (j = 1; j <= gmx->L; j++)
    {
      gmx->dp[j] = gmx->dp[0] + pos;
      pos += j+1;
    }
  
  /* initialize */
  for (j = 0; j <= gmx->L; j++)
    for (d = 0; d <= j; d++)
      gmx->dp[j][d] = -eslINFINITY;
 
  return gmx;
 ERROR:
  return NULL;
}

void 
GMX_Destroy(GMX *gmx)
{
  if (gmx == NULL) return;
  if (gmx->dp) free(gmx->dp);
  free(gmx);
}

void
GMX_Dump(FILE *fp, GMX *gmx)
{
  int j, d;

  for (j = 0; j <= gmx->L; j++)
    {
      for (d = 0; d <= j; d++)
        fprintf(fp, "%f ", gmx->dp[j][d]);
      fputc('\n', fp);
    }
  fputc('\n', fp);
}


