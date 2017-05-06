/* pottsscore.c */

#include "rscape_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

#include "correlators.h"
#include "logsum.h"
#include "pottsscore.h"
#include "pottsbuild.h"

int
potts_UnNormLogp(PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose)
{
  double logp = -eslINFINITY;
  double val;
  int    K = msa->abc->K;
  int    axi, axj;
  int    resi, resj;
  int    i, j;
  int    s;
  
  for (s = 0; s < msa->nseq; s++) {
    
    for (i = 0; i < msa->alen-1; i ++) {
      axi = msa->ax[s][i];
      resi = esl_abc_XIsCanonical(msa->abc, axi);
      val = 0;
      if (resi) val -= pt->h[i][axi];
      
      for (j = 0; j < msa->alen-1; j ++) {
	axj  = msa->ax[s][j];
	resj = esl_abc_XIsCanonical(msa->abc, axj);
	if (j != i && resi && resj) val -= pt->e[i][j][IDX(axi,axj,K)];
      }
    }

    logp = e2_FLogsum(logp, val);
  }
  
  *ret_logp = logp;

  return eslOK;
}


int                 
potts_CalculateCOV(PT *pt, struct data_s *data, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  ESL_DMATRIX *COV = NULL;
  char        *errbuf = data->errbuf;
  double       tol = data->tol;
  double       minCOV = eslINFINITY;
  double       maxCOV = 0;
  double       cov;
  double       eij;
  int          verbose = data->verbose;
  int          K = pt->abc->K;
  int          i, j;
  int          a, b;
  int          status = eslOK;

  // Use the Frobenious norm with zero-sum gauge
  status = potts_GaugeZeroSum(pt, tol, data->errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "zerosum gauge failed");
  
  for (i = 0; i < pt->L-1; i ++) {
    for (j = 0; j < pt->L; j ++) {
      cov = 0;
      
      for (a = 0; a < K; a ++)
	for (b = 0; b < K; b ++) {
	  eij  = pt->e[i][j][IDX(a,b,K)];
	  cov += eij * eij;
	}
      cov = sqrt(cov);
      if (cov > maxCOV) { maxCOV = cov; }
      if (cov < minCOV) { minCOV = cov; }
      COV->mx[i][j]  = cov;
    }
  }

  if (verbose) {
    printf("POTTS[%f,%f]\n", minCOV, maxCOV);
    for (i = 0; i < pt->L-1; i++) 
      for (j = i+1; j < pt->L; j++) {
	printf("POTTS[%d][%d] = %f \n", i, j,COV->mx[i][j]);
      } 
  }

  esl_dmatrix_Destroy(COV);
  return status;

 ERROR:
  if (COV) esl_dmatrix_Destroy(COV);
  return status;
}


