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
#include "esl_vectorops.h"

#include "correlators.h"
#include "logsum.h"
#include "pottsscore.h"
#include "pottsbuild.h"

static double potts_seqH   (PT *pt, ESL_DSQ *sq);
static double potts_sumH   (PT *pt, ESL_MSA *msa);
static double potts_ml_logz(PT *pt);
static double potts_plm_regularize_l2        (PT *pt);
static double potts_plm_regularize_l1        (PT *pt);
static double potts_aplm_regularize_l2(int i, PT *pt);
static double potts_aplm_regularize_l1(int i, PT *pt);



int
potts_ML_NLogp(PT *pt, ESL_MSA *msa, double *ret_nlogp, char *errbuf, int verbose)
{
  double nlogp;
  double status;

  nlogp  = -potts_sumH(pt, msa);
  nlogp /= msa->nseq;
  nlogp += potts_ml_logz(pt);

  // l2-regularization
  nlogp += potts_plm_regularize_l2(pt);
  
  *ret_nlogp = nlogp;
  if (nlogp < 0) ESL_XFAIL(eslFAIL, errbuf, "potts_NMLogp %f is negative!\n", nlogp); 
   
  return eslOK;

 ERROR:
  return status;
}

int
potts_PLM_NLogp(PT *pt, ESL_MSA *msa, double *ret_nlogp, double *dnlogp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *H = NULL;
  double   nlogp;
  double   nlogpi;
  double   logzi; 
  double   logdv, logadd;
  double   lwgt;
  int      L   = msa->alen;
  int      Kg  = msa->abc->K+1;
  int      Kg2 = Kg*Kg;
  int      dim = L*Kg + 0.5*L*(L-1)*Kg2;
  int      dofunc;
  int      dodfunc;
  int      resi, resj;
  int      x;
  int      s;
  int      i, j;
  int      a, b;
  int      status;

  dofunc  = (ret_nlogp)? TRUE : FALSE;
  dodfunc = (dnlogp)?    TRUE : FALSE;
  
  /* allocate */
  ESL_ALLOC(H, sizeof(double) * Kg);
  
  // Initialize
  if (dofunc)  nlogp = 0.0;
  if (dodfunc) esl_vec_DSet(dnlogp, dim, 0.);
  
  for (i = 0; i < L; i ++) {
    // Initialize
    if (dofunc)  nlogpi = 0.0;
    
    for (s = 0; s < msa->nseq; s++) {
      // Initialize
      x     = Kg*i + Kg2*i*(L - 0.5*i - 0.5) ; //\sum_{j=0}^{i-1} [ Kg + Kg2*(L-1-j)]
      sq    = msa->ax[s];
      resi  = sq[i+1];
      logzi = -eslINFINITY;
      lwgt  = log(msa->wgt[s]);
      
      for (a = 0; a < Kg; a++) {
	H[a] = potts_aplm_H(i, a, pt, sq);
	logzi = e2_DLogsum(logzi, H[a]); 
      }
      if (dofunc) nlogpi += msa->wgt[s] * (-H[resi] + logzi);
      
      if (dodfunc == FALSE) continue;

      // assign
      // derivative respect to hi(a)
      for (a = 0; a < Kg; a++) {
	logdv       = log(dnlogp[x]);
	logdv       = e2_DLogsum(logdv, H[a]-logzi+lwgt); 
	dnlogp[x++] = exp(logdv);	       
      }
      
      // derivative respect to eij(a,b)
      for (j = i+1; j < L; j ++) {
	resj = sq[j+1];
	for (a = 0; a < Kg; a++)
	  for (b = 0; b < Kg; b++) {
	    logdv       = log(dnlogp[x]);
	    logadd      = (b==resj)? H[a]-logzi+lwgt : -eslINFINITY;
	    logdv       = e2_DLogsum(logdv, logadd); 
	    dnlogp[x++] = exp(logdv);
	  }
      }
    } // for all sequences
    
    if (dofunc) nlogp += nlogpi;
  }
  
  if (dodfunc) {
    for (i = 0; i < L; i ++) {
      for (s = 0; s < msa->nseq; s++) {
	// Initialize
	x    = Kg*i + Kg2*i*(L - 0.5*i - 0.5) ; //\sum_{j=0}^{i-1} [ Kg + Kg2*(L-1-j)]
	sq   = msa->ax[s];
	resi = sq[i+1];
	
	// derivative respect to hi(a)
	for (a = 0; a < Kg; a++)
	  dnlogp[x++] += (a==resi)? -msa->wgt[s]:0.;
	
	// derivative respect to eij(a,b)
	for (j = i+1; j < L; j ++) {
	  resj = sq[j+1];
	  for (a = 0; a < Kg; a++)
	    for (b = 0; b < Kg; b++) 
	      dnlogp[x++] += (a==resi && b==resj)? -msa->wgt[s]:0.;
	}      
      }
    }
  }
  
  // l2-regularization
  if (dofunc) {
    nlogp += potts_plm_regularize_l2(pt);
    if (nlogp < 0) ESL_XFAIL(eslFAIL, errbuf, "potts_PLM_NLogp %f is negative!\n", nlogp);
  }
  
  if (dodfunc) {
    x = 0;
    for (i = 0; i < L; i ++) {     
      for (a = 0; a < Kg; a++) dnlogp[x++] += pt->muh * 2.0 * pt->h[i][a];         
      for (j = i+1; j < L; j ++)
	for (a = 0; a < Kg; a++)
	  for (b = 0; b < Kg; b++) 
	    dnlogp[x++] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
    }
  }
  
  // return
  if (dofunc) *ret_nlogp = nlogp;
  
  free(H);
  return eslOK;
  
 ERROR:
  if (H) free(H);
  return status;
}
  
int
potts_APLM_NLogp(int i, PT *pt, ESL_MSA *msa, double *ret_nlogp, int x0, double *dnlogp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *H = NULL;
  double   logz; 
  double   nlogp;
  double   logdv, logadd;
  double   lwgt;
  int      L  = msa->alen;
  int      Kg = msa->abc->K+1;
  int      dim = Kg + (L-1)*Kg*Kg;
  int      dofunc;
  int      dodfunc;
  int      resi, resj;
  int      x;
  int      s;
  int      j;
  int      a, b;
  int      status;

  dofunc  = (ret_nlogp)? TRUE : FALSE;
  dodfunc = (dnlogp)?    TRUE : FALSE;

  /* allocate */
  ESL_ALLOC(H, sizeof(double) * Kg);

  // Initialize
  if (dofunc)  nlogp = 0.0;
  if (dodfunc) esl_vec_DSet(dnlogp, dim, 0.);
  
  for (s = 0; s < msa->nseq; s++) {

    lwgt = log(msa->wgt[s]);
    
    // Initialize
    x    = x0;
    sq   = msa->ax[s];
    resi = sq[i+1];
    logz = -eslINFINITY;
    
    for (a = 0; a < Kg; a++) {
      H[a] = potts_aplm_H(i, a, pt, sq);
      logz = e2_DLogsum(logz, H[a]); 
    }
    if (dofunc) nlogp += msa->wgt[s] * (-H[resi] + logz);
 
    if (dodfunc == FALSE) continue;
    
    // assign
    // derivative respect to hi(a)
    for (a = 0; a < Kg; a++) {
      logdv       = log(dnlogp[x]);
      logdv       = e2_DLogsum(logdv, H[a]-logz+lwgt); 
      dnlogp[x++] = exp(logdv);	       
    }
    
    // derivative respect to eij(a,b)
    for (j = 0; j < i; j ++) {
      resj = sq[j+1];
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) {
	  logdv       = log(dnlogp[x]);
	  logadd      = (b==resj)? H[a]-logz+lwgt : -eslINFINITY;
	  logdv       = e2_DLogsum(logdv, logadd); 
	  dnlogp[x++] = exp(logdv);	       
	}
    }
    for (j = i+1; j < L; j ++) {
      resj = sq[j+1];
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) {
	  logdv       = log(dnlogp[x]);
	  logadd      = (b==resj)? H[a]-logz+lwgt : -eslINFINITY;
	  logdv       = e2_DLogsum(logdv, logadd); 
	  dnlogp[x++] = exp(logdv);
	}
    }
  } // for all sequences

  if (dodfunc) {
    for (s = 0; s < msa->nseq; s++) {
      // Initialize
      x    = x0;
      sq   = msa->ax[s];
      resi = sq[i+1];
      
      // derivative respect to hi(a)
      for (a = 0; a < Kg; a++)
	dnlogp[x++] += (a==resi)? -msa->wgt[s]:0.;
      
      // derivative respect to eij(a,b)
      for (j = 0; j < i; j ++) {
	resj = sq[j+1];
	for (a = 0; a < Kg; a++)
	  for (b = 0; b < Kg; b++) 
	    dnlogp[x++] += (a==resi && b==resj)? -msa->wgt[s]:0.;
      }
      for (j = i+1; j < L; j ++) {
	resj = sq[j+1];
	for (a = 0; a < Kg; a++)
	  for (b = 0; b < Kg; b++) 
	    dnlogp[x++] += (a==resi && b==resj)? -msa->wgt[s]:0.;
      }      
    }
  }
 
  // l2-regularization
  if (dofunc) {
    nlogp += potts_aplm_regularize_l2(i, pt);
    if (nlogp < 0) ESL_XFAIL(eslFAIL, errbuf, "potts_APLM_NLogp %f is negative!\n", nlogp);
  }

  if (dodfunc) {
    x = x0;
    for (a = 0; a < Kg; a++) dnlogp[x++] += pt->muh * 2.0 * pt->h[i][a];
    
    for (j = 0; j < i; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) 
	  dnlogp[x++] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
    for (j = i+1; j < L; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) 
	  dnlogp[x++] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
  }
  
  // return
  if (dofunc) *ret_nlogp = nlogp;

  free(H);
  return eslOK;

 ERROR:
  if (H) free(H);
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
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = cov;
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
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = cov;
    }
  }

  return status;
}

/* ------ static functions ----- */

static double
potts_seqH(PT *pt, ESL_DSQ *sq)
{
  ESL_ALPHABET *abc = pt->abc;
  double        sc = 0.;
  int           L = pt->L;
  int           Kg = abc->K+1;
  int           resi, resj;
  int           i, j;
  
  for (i = 0; i < L; i ++) {

    resi = sq[i+1];
    sc += pt->h[i][resi];
    
    for (j = 0; j < i; j ++) {
      resj = sq[j+1];
      sc  += pt->e[i][j][IDX(resi,resj,Kg)];
    }
    for (j = i+1; j < L; j ++) {
      resj = sq[j+1];
      sc  += pt->e[i][j][IDX(resi,resj,Kg)];
    }
  }
  return sc;
}

// sumH = sum_s [ \sum_i hi(s_i) + \sum{j!=i} eij(s_i,s_j) ]
static double
potts_sumH(PT *pt, ESL_MSA *msa)
{
  double logp = 0.;
  double sqsc;
  int    s;
    
  for (s = 0; s < msa->nseq; s++) {
    sqsc  = potts_seqH(pt, msa->ax[s]);    
    sqsc *= msa->wgt[s];
    logp += sqsc;
  }

  return logp;
}


// sum_{a1,..,aL} exp{ \sum_i hi(ai) + 1/2 \sum_{i!=j} eij(ai,aj) }
//
static double
potts_ml_logz(PT *pt)
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
      for (j = 0; j < L; j ++)  {
	if (j==i) continue;
	sc += 0.5 * pt->e[i][j][IDX(a[i],a[j],Kg)];
      }
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

// hi(a) + \sum_{j\neq i} eij(a,sqj)
//
double
potts_aplm_H(int i, int a, PT *pt, ESL_DSQ *sq)
{
  double lognum = 0.;
  int    L      = pt->L;
  int    Kg     = pt->abc->K+1;
  int    resj;
  int    j;
  
  lognum += pt->h[i][a];

  for (j = 0; j < i; j++) {
    resj    = sq[j+1];	  
    lognum += pt->e[i][j][IDX(a,resj,Kg)];
  }
  for (j = i+1; j < L; j++) {
    resj    = sq[j+1];	  
    lognum += pt->e[i][j][IDX(a,resj,Kg)];
  }

  return lognum;
}

// log { sum_{a} exp[ hi(a) + \sum_{j\neq i} eij(a,sqj) ] }
//
double
potts_aplm_logz(int i, PT *pt, ESL_DSQ *sq)
{
  double logz = -eslINFINITY;
  int    Kg     = pt->abc->K+1;
  int    a;
  
  for (a = 0; a < Kg; a++) 
    logz = e2_DLogsum(logz, potts_aplm_H(i, a, pt, sq));
 
  return logz;
}




static double
potts_plm_regularize_l2(PT *pt)
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
      reg += pt->muh * hi * hi;
      
      for (j = i+1; j < pt->L; j ++) {
	if (j==i) continue;
	for (b = 0; b < Kg; b ++) {
	  eij  = pt->e[i][j][IDX(a,b,Kg)];
	  reg += pt->mue * eij * eij;
	}
      }
    }
  
  return reg;
}

static double
potts_plm_regularize_l1(PT *pt)
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
      reg += pt->muh * fabs(hi);
      
      for (j = i+1; j < pt->L; j ++) {
	if (j==i) continue;
	for (b = 0; b < Kg; b ++) {
	  eij  = pt->e[i][j][IDX(a,b,Kg)];
	  reg += pt->mue * fabs(eij);
	}
      }
    }
  
  return reg;
}

static double
potts_aplm_regularize_l2(int i, PT *pt)
{
  double reg = 0;
  double hi, eij;
  int    Kg = pt->abc->K+1;
  int    j;
  int    a, b;
  
  // l2-regularization
  for (a = 0; a < Kg; a ++) {
    hi   = pt->h[i][a];
    reg += pt->muh * hi * hi;
      
    for (j = 0; j < pt->L; j ++) {
      if (j==i) continue;
      for (b = 0; b < Kg; b ++) {
	eij  = pt->e[i][j][IDX(a,b,Kg)];
	reg += pt->mue * eij * eij;
      }
    }
    
  }

  return reg;
}
static double
potts_aplm_regularize_l1(int i, PT *pt)
{
  double reg = 0;
  double hi, eij;
  int    Kg = pt->abc->K+1;
  int    j;
  int    a, b;
  
  // l2-regularization
  for (a = 0; a < Kg; a ++) {
    hi   = pt->h[i][a];
    reg += pt->muh * fabs(hi);
      
    for (j = 0; j < pt->L; j ++) {
      if (j==i) continue;
      for (b = 0; b < Kg; b ++) {
	eij  = pt->e[i][j][IDX(a,b,Kg)];
	reg += pt->mue * fabs(eij);
      }
    }
    
  }

  return reg;
}
