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

static double potts_H      (PT *pt, ESL_DSQ *sq);
static double potts_sumH   (PT *pt, ESL_MSA *msa);
static double potts_ml_logz(PT *pt);
static double potts_plm_regularize_l2        (PT *pt);
static double potts_plm_regularize_l1        (PT *pt);
static double potts_aplm_regularize_l2(int i, PT *pt);
static double potts_aplm_regularize_l1(int i, PT *pt);
static double potts_plm_regularize_l2_packed        (double *p, PT *pt);
static double potts_aplm_regularize_l2_packed(int i, double *p, PT *pt);
static double potts_logzi_packed_plm (int i, double *p, int L, int Kg, ESL_DSQ *sq, double *Hi);
static double potts_logzi_packed_aplm(int i, double *p, int L, int Kg, ESL_DSQ *sq, double *Hi);

int
potts_NLogp_ML(PT *pt, ESL_MSA *msa, double *ret_nlogp, char *errbuf, int verbose)
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
potts_NLogp_PLM(PT *pt, ESL_MSA *msa, double *ret_nlogp, double *dnlogp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *Hi = NULL;
  double   nlogp;
  double   nlogpi;
  double   logzi; 
  double   logadd;
  double   wgt, lwgt;
  int      L   = msa->alen;
  int      Kg  = msa->abc->K+1;
  int      dim = PLMDIM(L,L,Kg);
  int      dofunc;
  int      dodfunc;
  int      resi, resj;
  int      x0i, x;
  int      s;
  int      i, j;
  int      a, b;
  int      status;

  dofunc  = (ret_nlogp)? TRUE : FALSE;
  dodfunc = (dnlogp)?    TRUE : FALSE;
  
  // Initialize
  if (dofunc)  nlogp = 0.0;
  if (dodfunc) esl_vec_DSet(dnlogp, dim, -eslINFINITY); // working in log space
  ESL_ALLOC(Hi, sizeof(double)*Kg);
  
  for (i = 0; i < L; i ++) {
    // Initialize
    x0i = PLMDIM(i,L,Kg); //\sum_{j=0}^{i-1} [ Kg + Kg2*(L-1-j)]
    if (dofunc)  nlogpi = 0.0;
    
    for (s = 0; s < msa->nseq; s++) {
      // Initialize
      sq    = msa->ax[s];
      resi  = sq[i+1];
      lwgt  = log(msa->wgt[s]);

      // the hamiltonian
      //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
      // and logzi    = log \sum_a exp H^s_i[a]
      logzi = potts_Logzi(i, pt, sq, Hi);
       
      // the function 
      if (dofunc) nlogpi += msa->wgt[s] * (-Hi[resi] + logzi);
      
      // the gradient
      if (dodfunc == FALSE) continue;
      for (a = 0; a < Kg; a++) {
	logadd    = Hi[a]-logzi+lwgt;

	// derivative respect to hi(a)
	x         = x0i + a;
  	dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
    
	// derivative respect to eij(a,b) 
	for (j = i+1; j < L; j ++) {
	  resj      = sq[j+1];
	  x         = PLMIDX(i,j,L,Kg) + IDX(a,resj,Kg);
	  dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
	}
      
	for (j = 0; j < i; j ++) {
	  resj      = sq[j+1];
	  x         = PLMIDX(j,i,L,Kg) + IDX(resj,a,Kg);
	  dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
	}
      }
      
    } // for all sequences
    
    if (dofunc) nlogp += nlogpi;
  } // for all positions i

  
  // one more term to the gradient
  if (dodfunc) {
    esl_vec_DExp(dnlogp, dim);
    
    for (i = 0; i < L; i ++) {
      x0i = PLMDIM(i,L,Kg); //\sum_{j=0}^{i-1} [ Kg + Kg2*(L-1-j)]
      
      for (s = 0; s < msa->nseq; s++) {
	// Initialize
	sq   = msa->ax[s];
	resi = sq[i+1];
	wgt  = msa->wgt[s];
	
	// derivative respect to hi(a)
	x          = x0i + resi;
	dnlogp[x] -= wgt;
	
	// derivative respect to eij(a,b)
	for (j = i+1; j < L; j ++) {
	  resj       = sq[j+1];
	  x          = PLMIDX(i,j,L,Kg) + IDX(resi,resj,Kg);
	  dnlogp[x] -= wgt;	    
	}      
	for (j = 0; j < i; j ++) {
	  resj       = sq[j+1];
	  x          = PLMIDX(j,i,L,Kg) + IDX(resj,resi,Kg);
	  dnlogp[x] -= wgt;
	}      
      }
    }
  }
  
  // l2-regularization
  if (dofunc) 
    nlogp += potts_plm_regularize_l2(pt);
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
  
  free(Hi);
  return eslOK;
  
 ERROR:
  if (Hi) free(Hi);
  return status;
}

int
potts_NLogp_APLM(int i, PT *pt, ESL_MSA *msa, double *ret_nlogp, double *dnlogp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *Hi = NULL;
  double   logzi; 
  double   nlogp;
  double   logadd;
  double   wgt, lwgt;
  int      L   = msa->alen;
  int      Kg  = msa->abc->K+1;
  int      dim = APLMDIM(L,Kg);
  int      dofunc;
  int      dodfunc;
  int      resi, resj;
  int      x;
  int      s;
  int      j;
  int      a, b, ab;
  int      status;

  dofunc  = (ret_nlogp)? TRUE : FALSE;
  dodfunc = (dnlogp)?    TRUE : FALSE;

  // Initialize
  if (dofunc)  nlogp = 0.0;
  if (dodfunc) esl_vec_DSet(dnlogp, dim, -eslINFINITY); // working in logspace
  ESL_ALLOC(Hi, sizeof(double)*Kg);

  for (s = 0; s < msa->nseq; s++) {

    // Initialize
    sq   = msa->ax[s];
    resi = sq[i+1];
    lwgt = log(msa->wgt[s]);
    
    // the hamiltonian
    //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
    // and logzi    = log \sum_a exp H^s_i[a]
    logzi = potts_Logzi(i, pt, sq, Hi);   

    // the function
    if (dofunc)
      nlogp += msa->wgt[s] * (-Hi[resi] + logzi);
    
    // the gradient
    if (dodfunc == FALSE) continue;
    for (a = 0; a < Kg; a++) {
      logadd = Hi[a]-logzi+lwgt;
      
      // derivative respect to hi(a)
      dnlogp[a] = e2_DLogsum(dnlogp[a], logadd);
      
      // derivative respect to eij(a,b)
      for (j = 0; j < i; j ++) {
	resj      = sq[j+1];
	x         = APLMIDXL(j,Kg) + IDX(a,resj,Kg);
  	dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
      }
      for (j = i+1; j < L; j ++) {
	resj      = sq[j+1];
	x         = APLMIDXG(j,Kg) + IDX(a,resj,Kg);
  	dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
      }
    }
  } // for all sequences
  
  // one more term to the gradient
  if (dodfunc) {
    esl_vec_DExp(dnlogp, dim);
    
    for (s = 0; s < msa->nseq; s++) {
      // Initialize
      sq   = msa->ax[s];
      resi = sq[i+1];
      wgt  = msa->wgt[s];
      
      // derivative respect to hi(a)
      dnlogp[resi] -= wgt;
      
      // derivative respect to eij(a,b)
      for (j = 0; j < i; j ++) {
	resj       = sq[j+1];
	x          = APLMIDXL(j,Kg) + IDX(resi,resj,Kg);
	dnlogp[x] -= wgt;
      }
      for (j = i+1; j < L; j ++) {
	resj       = sq[j+1];
	x          = APLMIDXG(j,Kg) + IDX(resi,resj,Kg);
	dnlogp[x] -= wgt;
      }     
    }
  }
 
  // l2-regularization
  if (dofunc) 
    nlogp += potts_aplm_regularize_l2(i, pt);
  if (dodfunc) {
    for (a = 0; a < Kg; a++) dnlogp[a] += pt->muh * 2.0 * pt->h[i][a];
    
    for (j = 0; j < i; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) {
	  x          = APLMIDXL(j,Kg) + IDX(a,b,Kg);
	  dnlogp[x] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
	}
    for (j = i+1; j < L; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) {
	  x          = APLMIDXG(j,Kg) + IDX(a,b,Kg);
	  dnlogp[x] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
	}
  }
  
  // return
  if (dofunc) *ret_nlogp = nlogp;

  free(Hi);
  return eslOK;

 ERROR:
  if (Hi) free(Hi);
  return status;
}

int
potts_NLogp_PLM_Packed(int npt, double *p, PT *pt, ESL_MSA *msa, double *ret_nlogp, double *dnlogp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *Hi = NULL;
  double   nlogp;
  double   nlogpi;
  double   logzi; 
  double   logadd;
  double   wgt, lwgt;
  int      L   = msa->alen;
  int      Kg  = msa->abc->K+1;
  int      dim = PLMDIM(L,L,Kg);
  int      dofunc;
  int      dodfunc;
  int      resi, resj;
  int      x, x0i;
  int      s;
  int      i, j;
  int      a, b;
  int      status;
  
  dofunc  = (ret_nlogp)? TRUE : FALSE;
  dodfunc = (dnlogp)?    TRUE : FALSE;
  
  // Initialize
  if (dofunc)  nlogp = 0.0;
  if (dodfunc) esl_vec_DSet(dnlogp, dim, -eslINFINITY); // working in log space
  
  for (i = 0; i < L; i ++) {
    // Initialize
    x0i = PLMDIM(i,L,Kg); //\sum_{j=0}^{i-1} [ Kg + Kg2*(L-1-j)] 
    if (dofunc) nlogpi = 0.0;
    
    for (s = 0; s < msa->nseq; s++) {
      // Initialize
      sq   = msa->ax[s];
      resi = sq[i+1];
      lwgt = log(msa->wgt[s]);
      
      // the hamiltonian
      //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
      // and logzi    = log \sum_a exp H^s_i[a]
      logzi = potts_logzi_packed_plm(i, p, L, Kg, sq, Hi);   

      // the function
      if (dofunc) nlogpi += msa->wgt[s] * (-Hi[resi] + logzi);
      
      // the gradient
      if (dodfunc == FALSE) continue;
      for (a = 0; a < Kg; a++) {
	logadd = Hi[a]-logzi+lwgt;

	// derivative respect to hi(a)
	x         = x0i + a;
	dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
	
	// derivative respect to eij(a,b)
	for (j = i+1; j < L; j ++) {
	  resj      = sq[j+1];
	  x         = PLMIDX(i,j,L,Kg) + IDX(a,resj,Kg);
	  dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
	}
	for (j = 0; j < i; j ++) {
	  resj      = sq[j+1];
	  x         = PLMIDX(j,i,L,Kg) + IDX(resj,a,Kg);
	  dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
	}
      }
    } // for all sequences
    
    if (dofunc) nlogp += nlogpi;
  
  } // for all positions i

  // one more term to the gradient
  if (dodfunc) {
    esl_vec_DExp(dnlogp, dim);
      
    for (i = 0; i < L; i ++) {
      x0i = PLMDIM(i,L,Kg); //\sum_{j=0}^{i-1} [ Kg + Kg2*(L-1-j)]
      
      for (s = 0; s < msa->nseq; s++) {
	// Initialize
	sq   = msa->ax[s];
	resi = sq[i+1];
	wgt  = msa->wgt[s];
	      
	// derivative respect to hi(a)
	x          = x0i + resi;
	dnlogp[x] -= wgt;
	
	// derivative respect to eij(a,b)
	for (j = i+1; j < L; j ++) {
	  resj       = sq[j+1];
	  x          = PLMIDX(i,j,L,Kg) + IDX(resi,resj,Kg);
	  dnlogp[x] -= wgt;
	}      
	for (j = 0; j < i; j ++) {
	  resj       = sq[j+1];
	  x          = PLMIDX(j,i,L,Kg) + IDX(resj,resi,Kg);
	  dnlogp[x] -= wgt;
	}      
      }
    }
  }

  // l2-regularization
  if (dofunc) 
    nlogp += potts_plm_regularize_l2_packed(p, pt);
  if (dodfunc) {
    x = 0;
    for (i = 0; i < L; i ++) {     
      for (a = 0; a < Kg; a++)     { dnlogp[x] += pt->muh * 2.0 * p[x]; x ++; }        
      for (j = i+1; j < L; j ++)
	for (a = 0; a < Kg; a++)
	  for (b = 0; b < Kg; b++) { dnlogp[x] += pt->mue * 2.0 * p[x]; x ++; }
    }
  }
  
  // return
  if (dofunc) *ret_nlogp = nlogp;
  
  free(Hi);
  return eslOK;
  
 ERROR:
  if (Hi) free(Hi);
  return status;
}
  
int
potts_NLogp_APLM_Packed(int i, int np, double *p, PT *pt, ESL_MSA *msa, double *ret_nlogp, double *dnlogp, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *Hi = NULL;
  double   logzi; 
  double   nlogp;
  double   logadd;
  double   wgt, lwgt;
  int      L   = msa->alen;
  int      Kg  = msa->abc->K+1;
  int      dim = APLMDIM(L,Kg);
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

  // Initialize 
  if (dofunc)  nlogp = 0.;
  if (dodfunc) esl_vec_DSet(dnlogp, dim, -eslINFINITY); // working in log space
  
  for (s = 0; s < msa->nseq; s++) {
    
    // Initialize
    sq   = msa->ax[s];
    resi = sq[i+1];
    lwgt = log(msa->wgt[s]);
   
    // the hamiltonian
    //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
    // and logzi    = log \sum_a exp H^s_i[a]
    logzi = potts_logzi_packed_aplm(i, p, L, Kg, sq, Hi);   
    
    // the function
    if (dofunc)
      nlogp += msa->wgt[s] * (-Hi[resi] + logzi);
    
     // the gradient
    if (dodfunc == FALSE) continue;
    for (a = 0; a < Kg; a++) {
      logadd = Hi[a]-logzi+lwgt;
      
      // derivative respect to hi(a)
      dnlogp[a] = e2_DLogsum(dnlogp[a], logadd);
     
      // derivative respect to eij(a,b)
      for (j = 0; j < i; j ++) {
	resj      = sq[j+1];
 	x         = APLMIDXL(j,Kg) + IDX(a,resj,Kg);
	dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
      }
      for (j = i+1; j < L; j ++) {
	resj      = sq[j+1];
	x         = APLMIDXG(j,Kg) + IDX(a,resj,Kg);
	dnlogp[x] = e2_DLogsum(dnlogp[x], logadd);
      }
    }
    
  } // for all sequences
  
  // one more term to the gradient
  if (dodfunc) {
    esl_vec_DExp(dnlogp, dim); // first exponentiate
    
    for (s = 0; s < msa->nseq; s++) {
      
      // Initialize
      sq   = msa->ax[s];
      resi = sq[i+1];
      wgt  = msa->wgt[s];
      
      // derivative respect to hi(a)
      dnlogp[resi] -= wgt;
      
      // derivative respect to eij(a,b)
      for (j = 0; j < i; j ++) {
	resj       = sq[j+1];
	x          = APLMIDXL(j,Kg) + IDX(resi,resj,Kg);
	dnlogp[x] -= wgt;
      }
      for (j = i+1; j < L; j ++) {
	resj       = sq[j+1];
	x          = APLMIDXG(j,Kg) + IDX(resi,resj,Kg);
	dnlogp[x] -= wgt;
      }      
    }
  }
 
  // l2-regularization
  if (dofunc) 
    nlogp += potts_aplm_regularize_l2_packed(i, p, pt);
  if (dodfunc) {
    x = 0;
    for (a = 0; a < Kg; a++) { dnlogp[x] += pt->muh * 2.0 * p[x]; x++; }
    
    for (j = 0; j < i; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) {
	  dnlogp[x] += pt->mue * 2.0 * p[x]; x ++;
	}
    for (j = i+1; j < L; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) {
	  dnlogp[x] += pt->mue * 2.0 * p[x]; x ++;
	}
  }
  
  // return
  if (dofunc) *ret_nlogp = nlogp;

  free(Hi);
  return eslOK;

 ERROR:
  if (Hi) free(Hi);
  return status;
}

// hi(a) + \sum_{j\neq i} eij(a,sqj)
//
double
potts_Hi(int i, int a, PT *pt, ESL_DSQ *sq)
{
  double val = 0.;
  int    L   = pt->L;
  int    Kg  = pt->abc->K+1;
  int    resj;
  int    j;
  
  val += pt->h[i][a];

  for (j = 0; j < i; j++) {
    resj  = sq[j+1];	  
    val  += pt->e[i][j][IDX(a,resj,Kg)];
  }
  for (j = i+1; j < L; j++) {
    resj  = sq[j+1];	  
    val  += pt->e[i][j][IDX(a,resj,Kg)];
  }

  return val;
}
double
potts_Hi_APLM_Packed(int i, int a, double *p, int L, int Kg, ESL_DSQ *sq)
{
  double val = 0.;
  int    resj;
  int    j;
  int    x;
  
  val += p[a];
 
  for (j = 0; j < i; j++) {
    resj  = sq[j+1];	  
    x     = APLMIDXL(j,Kg) + IDX(a,resj,Kg);
    val  += p[x];
  }
  for (j = i+1; j < L; j++) {
    resj  = sq[j+1];	  
    x     = APLMIDXG(j,Kg) + IDX(a,resj,Kg);
    val  += p[x];
  }
  return val;
}

//exp[ hi(a) + \sum_{j\neq i} eij(a,sqj)
double
potts_Hi_PLM_Packed(int i, int a, double *p, int L, int Kg, ESL_DSQ *sq)
{
  double val = 0.;
  int    resj;
  int    j;
  int    x;
  
  x    = PLMDIM(i,L,Kg) + a;
  val += p[x];
  
  for (j = i+1; j < L; j++) {
    resj  = sq[j+1];
    x     = PLMIDX(i,j,L,Kg) + IDX(a,resj,Kg);
    val  += p[x];
  }
  for (j = 0; j < i; j++) {
    resj  = sq[j+1];
    x     = PLMIDX(j,i,L,Kg) + IDX(resj,a,Kg);
    val  += p[x];
  }
  return val;
}

// log { sum_{a} exp[ hi(a) + \sum_{j\neq i} eij(a,sqj) ] }
//
double
potts_Logzi(int i, PT *pt, ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  logzi = -eslINFINITY;
  int     Kg    = pt->abc->K+1;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi(i, a, pt, sq);
    if (Hi) Hi[a] = Ha;
     
    logzi = e2_DLogsum(logzi, Ha);
  }
  return logzi;

 ERROR:
  return status;
}

static double
potts_logzi_packed_plm(int i, double *p, int L, int Kg, ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  logzi = -eslINFINITY;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi_PLM_Packed(i, a, p, L, Kg, sq);
    if (Hi) Hi[a] = Ha;
     
    logzi = e2_DLogsum(logzi, Ha);
  }
  return logzi;

 ERROR:
  return status;
}

static double
potts_logzi_packed_aplm(int i, double *p, int L, int Kg, ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  logzi = -eslINFINITY;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi_APLM_Packed(i, a, p, L, Kg, sq);
    if (Hi) Hi[a] = Ha;
     
    logzi = e2_DLogsum(logzi, Ha);
  }
  return logzi;

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

      // only for residues -- no gaps
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
      
      // only for residues -- no gaps
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

// H = \sum_i hi(s_i) + \sum{j!=i} eij(s_i,s_j) 
static double
potts_H(PT *pt, ESL_DSQ *sq)
{
  ESL_ALPHABET *abc = pt->abc;
  double        H = 0.;
  int           L = pt->L;
  int           Kg = abc->K+1;
  int           resi, resj;
  int           i, j;
  
  for (i = 0; i < L; i ++) 
    H += potts_Hi(i, sq[i+1], pt, sq);
  
  return H;
}

// sumH = sum_s [ \sum_i hi(s_i) + \sum{j!=i} eij(s_i,s_j) ]
static double
potts_sumH(PT *pt, ESL_MSA *msa)
{
  double logp = 0.;
  double sumH;
  int    s;
    
  for (s = 0; s < msa->nseq; s++) {
    sumH  = potts_H(pt, msa->ax[s]);    
    sumH *= msa->wgt[s];
    logp += sumH;
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



static double
potts_plm_regularize_l2(PT *pt)
{
  double reg = 0.;
  double hi, eij;
  int    L  = pt->L;
  int    Kg = pt->abc->K+1;
  int    i, j;
  int    a, b;
  
  // l2-regularization
  for (i = 0; i < L; i ++) 
    for (a = 0; a < Kg; a ++) {
      hi   = pt->h[i][a];
      reg += pt->muh * hi * hi;
      
      for (j = i+1; j < L; j ++) 
	for (b = 0; b < Kg; b ++) {
	  eij  = pt->e[i][j][IDX(a,b,Kg)];
	  reg += pt->mue * eij * eij;
	}
    }
  
  return reg;
}

static double
potts_plm_regularize_l2_packed(double *p, PT *pt)
{
  double reg = 0.;
  double hi, eij;
  int    L  = pt->L;
  int    Kg = pt->abc->K+1;
  int    x0i, x;
  int    i, j;
  int    a, b;
  
  // l2-regularization
  for (i = 0; i < L; i ++) {
    x0i = PLMDIM(i,L,Kg);
    for (a = 0; a < Kg; a ++) {
      hi   = p[x0i+a];
      reg += pt->muh * hi * hi;
  
      for (j = i+1; j < L; j ++) 
	for (b = 0; b < Kg; b ++) {
	  x    = PLMIDX(i,j,L,Kg) + IDX(a,b,Kg);
	  eij  = p[x];
	  reg += pt->mue * eij * eij;
	}
    }
  }
    
  return reg;
}

static double
potts_plm_regularize_l1(PT *pt)
{
  double reg = 0.;
  double hi, eij;
  int    Kg = pt->abc->K+1;
  int    i, j;
  int    a, b;
  
  // l2-regularization
  for (i = 0; i < pt->L; i ++) 
    for (a = 0; a < Kg; a ++) {
      hi   = pt->h[i][a];
      reg += pt->muh * fabs(hi);
      
      for (j = i+1; j < pt->L; j ++) 
	for (b = 0; b < Kg; b ++) {
	  eij  = pt->e[i][j][IDX(a,b,Kg)];
	  reg += pt->mue * fabs(eij);
	}
    }

  return reg;
}

static double
potts_aplm_regularize_l2(int i, PT *pt)
{
  double reg = 0.;
  double hi, eij;
  int    Kg = pt->abc->K+1;
  int    j;
  int    a, b;
  
  // l2-regularization
  for (a = 0; a < Kg; a ++) {
    hi   = pt->h[i][a];
    reg += pt->muh * hi * hi;
    
    for (j = 0; j < i; j ++) 
      for (b = 0; b < Kg; b ++) {
	eij  = pt->e[i][j][IDX(a,b,Kg)];
	reg += pt->mue * eij * eij;
      }
    for (j = 0; j < pt->L; j ++) 
      for (b = 0; b < Kg; b ++) {
	eij  = pt->e[i][j][IDX(a,b,Kg)];
	reg += pt->mue * eij * eij;
      }
  }
  return reg;
}

static double
potts_aplm_regularize_l2_packed(int i, double *p, PT *pt)
{
  double reg = 0.;
  double hi, eij;
  int    L = pt->L;
  int    Kg = pt->abc->K+1;
  int    x;
  int    j;
  int    a, b;
  
  // l2-regularization
  for (a = 0; a < Kg; a ++) {
    hi   = p[a];
    reg += pt->muh * hi * hi;
    
    for (j = 0; j < i; j ++) 
      for (b = 0; b < Kg; b ++) {
	x    = APLMIDXL(j,Kg) + IDX(a,b,Kg);
	eij  = p[x];
	reg += pt->mue * eij * eij;
      }  
    for (j = i+1; j < L; j ++) 
      for (b = 0; b < Kg; b ++) {
	x    = APLMIDXG(j,Kg) + IDX(a,b,Kg);
	eij  = p[x];
	reg += pt->mue * eij * eij;
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
      
    for (j = 0; j < i; j ++) 
      for (b = 0; b < Kg; b ++) {
	eij  = pt->e[i][j][IDX(a,b,Kg)];
	reg += pt->mue * fabs(eij);
      }
    for (j = i+1; j < pt->L; j ++) 
      for (b = 0; b < Kg; b ++) {
	eij  = pt->e[i][j][IDX(a,b,Kg)];
	reg += pt->mue * fabs(eij);
      }    
  }
  return reg;
}
