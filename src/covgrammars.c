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

G6_MX *
G6MX_Create(int L)
{
  G6_MX *g6mx = NULL;
  int    status;

  ESL_ALLOC(g6mx, sizeof(G6_MX));

  g6mx->S = GMX_Create(L);
  g6mx->L = GMX_Create(L);
  g6mx->F = GMX_Create(L);

  return g6mx;
 ERROR:
  return NULL;
}

BGR_MX *
BGRMX_Create(int L)
{
  BGR_MX *bgrmx = NULL;
  int     status;

  ESL_ALLOC(bgrmx, sizeof(BGR_MX));

  bgrmx->S  = GMX_Create(L);
  bgrmx->F0 = GMX_Create(L);
  bgrmx->F5 = GMX_Create(L);
  bgrmx->P  = GMX_Create(L);
  bgrmx->M  = GMX_Create(L);
  bgrmx->R  = GMX_Create(L);
  bgrmx->M1 = GMX_Create(L);

  return bgrmx;
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
G6MX_Destroy(G6_MX *gmx)
{
  if (gmx == NULL) return;
  if (gmx->S) GMX_Destroy(gmx->S);
  if (gmx->L) GMX_Destroy(gmx->L);
  if (gmx->F) GMX_Destroy(gmx->F);
  free(gmx);
}
void 
BGRMX_Destroy(BGR_MX *gmx)
{
  if (gmx == NULL) return;
  if (gmx->S)  GMX_Destroy(gmx->S);
  if (gmx->F0) GMX_Destroy(gmx->F0);
  if (gmx->F5) GMX_Destroy(gmx->F5);
  if (gmx->P)  GMX_Destroy(gmx->P);
  if (gmx->M)  GMX_Destroy(gmx->M);
  if (gmx->R)  GMX_Destroy(gmx->R);
  if (gmx->M1) GMX_Destroy(gmx->M1);
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


