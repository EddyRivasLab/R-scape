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

static double potts_score_oneseq(PT *pt, ESL_DSQ *sq);
static double potts_full_logz(PT *pt);
static double potts_aplm_logz(int i, PT *pt, ESL_DSQ *sq);
static double potts_regularize(PT *pt);


int
potts_SumH(PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose)
{
  double logp = 0.;
  double sqsc;
  int    s;
  
  for (s = 0; s < msa->nseq; s++) {
    sqsc  = potts_score_oneseq(pt, msa->ax[s]);    
    sqsc *= msa->wgt[s];
    logp += sqsc;
  }
  *ret_logp = logp;

  return eslOK;
}


int
potts_FULLLogp(PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose)
{
  double logp;

  potts_SumH(pt, msa, &logp, errbuf, verbose);
  logp /= msa->nseq;
  logp -= potts_full_logz(pt);

  // l2-regularization
  logp += potts_regularize(pt);
  
  *ret_logp = logp;
  return eslOK;
}

int
potts_APLMLogp(PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double   logp = 0.;
  double   sc;
  double   hi, eij;
  int      K = msa->abc->K;
  int      sqi, sqj;
  int      resi, resj;
  int      s;
  int      i, j;
  int      a, b;
  
  for (s = 0; s < msa->nseq; s++) {
    sq = msa->ax[s];
    
    for (i = 0; i < pt->L; i ++) {
      sqi  = sq[i];
      resi = esl_abc_XIsCanonical(msa->abc, sqi);
      
      hi = (resi)? pt->h[i][sqi] : 0;
      sc = hi;
      
      for (j = 0; j < pt->L; j ++) {
	if (j == i) continue;
	
	sqj  = sq[j];
	resj = esl_abc_XIsCanonical(msa->abc, sqj);
	
	eij = (resi && resj)? pt->e[i][j][IDX(sqi,sqj,K)] : 0;
	sc += eij;
      }
      sc -= potts_aplm_logz(i, pt, sq);
      sc *= msa->wgt[s];
      
      logp += sc;
    }
    logp /= msa->nseq;  
  }
  
  // l2-regularization
  logp += potts_regularize(pt);
  
  *ret_logp = logp;

  return eslOK;
}



int                 
potts_CalculateCOV(struct data_s *data, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  PT              *pt = mi->pt;
  char            *errbuf = data->errbuf;
  double           tol = data->tol;
  double           cov;
  double           eij;
  int              verbose = data->verbose;
  int              K = pt->abc->K;
  int              i, j;
  int              a, b;
  int              status = eslOK;

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
      if (cov > mi->maxCOV) { mi->maxCOV = cov; }
      if (cov < mi->minCOV) { mi->minCOV = cov; }
      mi->COV->mx[i][j]  = cov;
    }
  }

  if (verbose) {
    printf("POTTS[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < pt->L-1; i++) 
      for (j = i+1; j < pt->L; j++) {
	printf("POTTS[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }

  return status;

 ERROR:
  return status;
}


static double
potts_score_oneseq(PT *pt, ESL_DSQ *sq)
{
  ESL_ALPHABET *abc = pt->abc;
  double        sc = 0.;
  double        hi, eij;
  int           L = pt->L;
  int           K = abc->K;
  int           sqi, sqj;
  int           resi, resj;
  int           i, j;
  
  for (i = 0; i < L-1; i ++) {
    sqi  = sq[i];
    resi = esl_abc_XIsCanonical(abc, sqi);
    
    hi = (resi)? pt->h[i][sqi] : 0;
    sc += hi;
    
    for (j = i+1; j < L; j ++) {
      sqj  = sq[j];
      resj = esl_abc_XIsCanonical(abc, sqj);
      
      eij = (resi && resj)? pt->e[i][j][IDX(sqi,sqj,K)] : 0;
      sc += eij;
    }
  }
  return sc;
}


// sum_{a,b} exp{ \sim_i hi(a) + \sum_{i<j} eij(a,b) }
//
static double
potts_full_logz(PT *pt)
{
  double logsum = -eslINFINITY;
  double sc;
  int    K = pt->abc->K;
  int    i, j;
  int    a, b;
  
  for (a = 0; a < K; a ++)
    for (b = 0; b < K; b ++) {
      sc = 0.;
      
      for (i = 0; i < pt->L; i ++) {
	sc += pt->h[i][a];	
	for (j = i+1; j < pt->L; j ++) 
	  sc += pt->e[i][j][IDX(a,b,K)];
      }
      
      logsum = e2_FLogsum(logsum, sc);
    }
  
  return logsum;
}

// sum_{a} exp[ hi(a) + \sum_{j\neq i} eij(a,sqj) ]
//
static double
potts_aplm_logz(int i, PT *pt, ESL_DSQ *sq)
{
  double logsum = -eslINFINITY;
  double sc;
  int    K = pt->abc->K;
  int    sqj;
  int    resj;
  int    j;
  int    a;

  for (a = 0; a < K; a ++) {
    sc = pt->h[i][a];
    
    for (j = 0; j < pt->L; j ++) {
      if (j == i) continue;
      sqj  = sq[j];
      resj = esl_abc_XIsCanonical(pt->abc, sqj);
      sc += (resj)? pt->e[i][j][IDX(a,sqj,K)] : 0.;      
    }
    
    logsum = e2_FLogsum(logsum, sc);
  }

  return logsum;
}


static double
potts_regularize(PT *pt)
{
  double reg = 0;
  double hi, eij;
  int    K = pt->abc->K;
  int    i, j;
  int    a, b;
  
  // l2-regularization
  for (a = 0; a < K; a ++) {
    hi = pt->h[i][a];
    reg += hi*hi;
    
    for (j = 0; j < pt->L; j ++) {
      if (j ==i) continue;
      
      for (b = 0; b < K; b ++) {
	eij = pt->e[i][j][IDX(a,b,K)];
	reg += eij*eij;
      }
    }
  }
  reg *= pt->mu;

  return reg;
}
