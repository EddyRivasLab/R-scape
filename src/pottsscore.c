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
static double potts_all_logz(PT *pt);
static double potts_all_regularize(PT *pt);
static double potts_aplm_regularize(int pos, PT *pt);


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
potts_MLLogp(PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose)
{
  double logp;
  double status;

  potts_SumH(pt, msa, &logp, errbuf, verbose);
  logp /= msa->nseq;
  logp -= potts_all_logz(pt);

  // l2-regularization
  logp -= potts_all_regularize(pt);
  
  *ret_logp = logp;
  if (logp > 0) ESL_XFAIL(eslFAIL, errbuf, "potts_MLogp %f is positive!\n", logp); 
   
  return eslOK;

 ERROR:
  return status;
}

int
potts_PLMLogp(PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose)
{
  double   logp = 0.;
  double   sc;
  int      i;
  int      status;
  
  for (i = 0; i < pt->L; i ++) {
    status = potts_APLMLogp(i, pt, msa, &sc, errbuf, verbose);
    if (status != eslOK) return status;
    
    logp += sc;
  } 
  // l2-regularization
  logp -= potts_all_regularize(pt);
  
  *ret_logp = logp;
  if (logp > 0) ESL_XFAIL(eslFAIL, errbuf, "potts_PLMLogp %f is positive!\n", logp); 
  
  return eslOK;

 ERROR:
  return status;
}

int
potts_APLMLogp(int pos, PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double   tol = 1e-3;
  double   logp = 0.;
  double   lognum, logden;
  double   sc;
  int      Kg = msa->abc->K+1;
  int      sqi, sqj;
  int      resi, resj;
  int      s;
  int      status;
  
  for (s = 0; s < msa->nseq; s++) {
    sq = msa->ax[s];
    
    sqi  = sq[pos+1];
    resi = esl_abc_XIsCanonical(msa->abc, sqi);
      
    lognum = potts_APLMLognum(pos, sqi, pt, sq);
    logden = potts_APLMLogz  (pos,      pt, sq);
   
    if (lognum > logden) {
      if (lognum - logden < tol) {
	lognum = logden;
      }
      else {
	printf(" i %d seqi %d num %f > den %f\n", pos, sqi, lognum, logden);
	exit(1);
      }
    }
    sc = lognum - logden;
    sc *= msa->wgt[s];
    
    logp += sc;
  }
  logp /= msa->nseq;  
  
  // l2-regularization
  logp -= potts_aplm_regularize(pos, pt);
  *ret_logp = logp;

  if (logp > 0) ESL_XFAIL(eslFAIL, errbuf, "potts_APLMLogp %f is positive!\n", logp); 
 
  return eslOK;

 ERROR:
  return status;
}



int                 
potts_CalculateCOV(struct data_s *data)
{
  struct mutual_s *mi = data->mi;
  PT              *pt = data->pt;
  int              i, j;
  int              status;

  switch(pt->sctype) {
  case SCNONE:
    break;
  case FROEB:
    status = potts_CalculateCOVFrobenius(data);
    if (status != eslOK) ESL_XFAIL(eslFAIL, data->errbuf, "potts_CalculateCOVFrobenius() error");
    break;
  case AVG:
   status = potts_CalculateCOVAverage(data);
   if (status != eslOK) ESL_XFAIL(eslFAIL, data->errbuf, "potts_CalculateCOVAverage() error");
     break;
  case DI:
    ESL_XFAIL(eslFAIL, data->errbuf, "potts_CalculateCOVDI() not implemented yet");
    break;
  }
  
  if (data->verbose) {
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
int                 
potts_CalculateCOVFrobenius(struct data_s *data)
{
  struct mutual_s *mi = data->mi;
  PT              *pt = data->pt;
  char            *errbuf = data->errbuf;
  double           cov;
  double           eij;
  int              K = pt->abc->K;
  int              Kg = K+1;
  int              i, j;
  int              a, b;
  int              idx;
  int              status = eslOK;

  // Use the Frobenius norm with zero-sum gauge
  status = potts_GaugeZeroSum(pt, data->errbuf, data->verbose);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "zerosum gauge failed");
  
  for (i = 0; i < pt->L; i ++) {
    for (j = i+1; j < pt->L; j ++) {
      cov = 0;
      
      for (a = 0; a < K; a ++)
	for (b = 0; b < K; b ++) {
	  idx  = IDX(a,b,Kg);
	  eij  = pt->e[i][j][idx];
	  cov += eij * eij;
	}
      cov = sqrt(cov);
      if (cov > mi->maxCOV) { mi->maxCOV = cov; }
      if (cov < mi->minCOV) { mi->minCOV = cov; }
      mi->COV->mx[i][j]  = cov;
    }
  }

  return status;

 ERROR:
  return status;
}

int                 
potts_CalculateCOVAverage(struct data_s *data)
{
  struct mutual_s *mi = data->mi;
  PT              *pt = data->pt;
  double           cov;
  double           eij;
  int              K  = pt->abc->K;
  int              Kg = K+1;
  int              i, j;
  int              a, b;
  int              idx;
  int              status = eslOK;

  for (i = 0; i < pt->L; i ++) {
    for (j = i+1; j < pt->L; j ++) {
      cov = 0;
      
      for (a = 0; a < K; a ++)
	for (b = 0; b < K; b ++) {
	  idx  = IDX(a,b,Kg);
	  eij  = pt->e[i][j][idx];
	  cov += eij;
	}
      if (cov > mi->maxCOV) { mi->maxCOV = cov; }
      if (cov < mi->minCOV) { mi->minCOV = cov; }
      mi->COV->mx[i][j] = cov;
    }
  }

  return status;
}

// hi(a) + \sum_{j\neq i} eij(a,sqj)
//
double
potts_APLMLognum(int i, int a, PT *pt, ESL_DSQ *sq)
{
  double lognum = 0.;
  int    L      = pt->L;
  int    Kg     = pt->abc->K+1;
  int    resj;
  int    j;
  
  lognum += pt->h[i][a];
  
  for (j = 0; j < L; j++) {
    resj    = sq[j+1];	  
    lognum += 0.5 * pt->e[i][j][IDX(a,resj,Kg)];
  }
  
  return lognum;
}

// log { sum_{a} exp[ hi(a) + \sum_{j\neq i} eij(a,sqj) ] }
//
double
potts_APLMLogz(int i, PT *pt, ESL_DSQ *sq)
{
  double logz = -eslINFINITY;
  int    Kg     = pt->abc->K+1;
  int    a;
  
  for (a = 0; a < Kg; a++) 
    logz = e2_FLogsum(logz, potts_APLMLognum(i, a, pt, sq));
 
  return logz;
}


static double
potts_score_oneseq(PT *pt, ESL_DSQ *sq)
{
  ESL_ALPHABET *abc = pt->abc;
  double        sc = 0.;
  double        hi, eij;
  int           L = pt->L;
  int           Kg = abc->K+1;
  int           sqi, sqj;
  int           resi, resj;
  int           i, j;
  
  for (i = 0; i < L; i ++) {
    sqi  = sq[i+1];
    resi = esl_abc_XIsCanonical(abc, sqi);
    
    hi = pt->h[i][sqi];
    sc += hi;
    
    for (j = 0; j < L; j ++) {
      sqj  = sq[j+1];
      resj = esl_abc_XIsCanonical(abc, sqj);
      
      eij = pt->e[i][j][IDX(sqi,sqj,Kg)];
      sc += 0.5 * eij;
    }
  }
  return sc;
}


// sum_{a1,..,aL} exp{ \sim_i hi(ai) + \sum_{i<j} eij(ai,aj) }
//
static double
potts_all_logz(PT *pt)
{
  int    *a = NULL;
  double  logsum = -eslINFINITY;
  double  sc;
  int     L  = pt->L;
  int     Kg = pt->abc->K+1;
  int     x = 0;
  int     i, j;
  int     status;
  
  ESL_ALLOC(a, sizeof(int)*L);
  for (i = 0; i < L; i ++) a[i] = 0;

  while(TRUE) {

    sc = 0.;
    for (i = 0; i < L; i ++) {
      sc += pt->h[i][a[i]];	
      for (j = 0; j < L; j ++) 
	sc += pt->e[i][j][IDX(a[i],a[j],Kg)];
    }
    logsum = e2_FLogsum(logsum, sc);
 
    // increment
    a[0] ++;

    // carry
    while (a[x] == Kg) {
      // overflow you are done
      if (x == L-1) { free(a); return logsum; }
      
      a[x++] = 0;
      a[x] ++;
    }
    
    x = 0;
  }

 ERROR:
  return logsum;
}



static double
potts_all_regularize(PT *pt)
{
  double reg = 0;
  double hi, eij;
  int    Kg = pt->abc->K+1;
  int    i, j;
  int    a, b;
  
  // l2-regularization
  for (i = 0; i < pt->L; i ++) 
    for (a = 0; a < Kg; a ++) {
      hi   = pt->h[i][a];
      reg += hi*hi;
      
      for (j = i+1; j < pt->L; j ++) 
	for (b = 0; b < Kg; b ++) {
	  eij  = pt->e[i][j][IDX(a,b,Kg)];
	  reg += eij*eij;
	}
      
    }
  reg *= pt->mu;
  
  return reg;
}

static double
potts_aplm_regularize(int i, PT *pt)
{
  double reg = 0;
  double hi, eij;
  int    Kg = pt->abc->K+1;
  int    j;
  int    a, b;
  
  // l2-regularization
  for (a = 0; a < Kg; a ++) {
    hi   = pt->h[i][a];
    reg += hi*hi;
    
    for (j = i+1; j < pt->L; j ++) 
      for (b = 0; b < Kg; b ++) {
	eij  = pt->e[i][j][IDX(a,b,Kg)];
	reg += eij*eij;
      }
    
  }
  reg *= pt->mu;

  return reg;
}
