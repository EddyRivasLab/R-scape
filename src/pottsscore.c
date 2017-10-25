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
static double potts_logzi_packed_plm (int i, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq, double *Hi);
static double potts_logzi_packed_aplm(int i, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq, double *Hi);
static double potts_Zi_packed_plm (int i, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq, double *Hi);
static double potts_Zi_packed_aplm(int i, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq, double *Hi);

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
potts_NLogp_PLM(PT *pt, ESL_MSA *msa, double *ret_nlogp, PT *gr, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *Hi = NULL;
  double   nlogp;
  double   nlogpi;
  double   zi; 
  double   add;
  double   wgt;
  int      L   = msa->alen;
  int      Kg  = pt->Kg;
  int      dofunc;
  int      dodfunc;
  int      resi, resj;
  int      s;
  int      i, j;
  int      a, b;
  int      status;

  dofunc  = (ret_nlogp)? TRUE : FALSE;
  dodfunc = (gr)?        TRUE : FALSE;
  
  // Initialize
  if (dofunc)  nlogp = 0.0;
  if (dodfunc) potts_AssignZero(gr, errbuf, verbose);
  
  ESL_ALLOC(Hi, sizeof(double)*Kg);
  
  for (i = 0; i < L; i ++) {
    // Initialize
    if (dofunc)  nlogpi = 0.0;
    
    for (s = 0; s < msa->nseq; s++) {
      // Initialize
      sq   = msa->ax[s];
      resi = sq[i+1];
      wgt  = msa->wgt[s];

      // the hamiltonian
      //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
      // and logzi = log \sum_a exp H^s_i[a]
      // and zi    = \sum_a exp H^s_i[a]
      //logzi = potts_Logzi(i, pt, sq, Hi);
      zi = potts_Zi(i, pt, sq, Hi);
       
      // the function 
      if (dofunc) nlogpi += wgt * (-Hi[resi] + log(zi));
      
      // the gradient
      if (dodfunc == FALSE) continue;
      for (a = 0; a < Kg; a++) {
	//add = wgt * exp(Hi[a]-logzi);
	add = wgt * exp(Hi[a]) / zi;

	// derivative respect to hi(a)
  	gr->h[i][a] += add;
    
	// derivative respect to eij(a,b) 
	for (j = i+1; j < L; j ++) {
	  resj = sq[j+1];
	  gr->e[i][j][IDX(a,resj,Kg)] += add;
	}
      
	for (j = 0; j < i; j ++) {
	  resj = sq[j+1];
	  gr->e[j][i][IDX(resj,a,Kg)] += add;
	}
      }
      
      // one more term to the gradient
      // derivative respect to hi(resi)
      gr->h[i][resi] -= wgt;

      // derivative respect to eij(resi,resj)
      for (j = i+1; j < L; j ++) {
	resj = sq[j+1];	    
	gr->e[i][j][IDX(resi,resj,Kg)] -= wgt;
      }      
      for (j = 0; j < i; j ++) {
	resj       = sq[j+1];
	gr->e[j][i][IDX(resj,resi,Kg)] -= wgt;
      }      
      
    } // for all sequences
    
    if (dofunc) nlogp += nlogpi;
  } // for all positions i

  
  // l2-regularization
  if (dofunc) 
    nlogp += potts_plm_regularize_l2(pt);
  
  if (dodfunc) {
    for (i = 0; i < L; i ++) {     
      for (a = 0; a < Kg; a++) gr->h[i][a] += pt->muh * 2.0 * pt->h[i][a];
      
      for (j = 0; j < i; j ++)
	for (a = 0; a < Kg; a++)
	  for (b = 0; b < Kg; b++) 
	    gr->e[i][j][IDX(a,b,Kg)] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
      for (j = i+1; j < L; j ++)
	for (a = 0; a < Kg; a++)
	  for (b = 0; b < Kg; b++) 
	    gr->e[i][j][IDX(a,b,Kg)] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
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
potts_NLogp_APLM(int i, PT *pt, ESL_MSA *msa, double *ret_nlogp, PT *gr, char *errbuf, int verbose)
{
  ESL_DSQ *sq;
  double  *Hi = NULL;
  double   zi; 
  double   nlogp;
  double   add;
  double   wgt;
  int      L   = msa->alen;
  int      Kg  = pt->Kg;
  int      dofunc;
  int      dodfunc;
  int      resi, resj;
  int      s;
  int      j;
  int      a, b, ab;
  int      status;

  dofunc  = (ret_nlogp)? TRUE : FALSE;
  dodfunc = (gr)?        TRUE : FALSE;

  // Initialize
  if (dofunc)  nlogp = 0.0;
  if (dodfunc) potts_AssignZero(gr, errbuf, verbose);
  ESL_ALLOC(Hi, sizeof(double)*Kg);

  for (s = 0; s < msa->nseq; s++) {

    // Initialize
    sq   = msa->ax[s];
    resi = sq[i+1];
    wgt = msa->wgt[s];
    
    // the hamiltonian
    //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
    // and logzi    = log \sum_a exp H^s_i[a]
    // and zi       = \sum_a exp H^s_i[a]
    //logzi = potts_Logzi(i, pt, sq, Hi);   
    zi = potts_Zi(i, pt, sq, Hi);   

    // the function
    if (dofunc)
      nlogp += wgt * (-Hi[resi] + log(zi));
    
    // the gradient
    if (dodfunc == FALSE) continue;
    for (a = 0; a < Kg; a++) {
      add = wgt * exp(Hi[a]) / zi;
      
      // derivative respect to hi(a)
      gr->h[i][a] += add;
      
      // derivative respect to eij(a,b)
      for (j = 0; j < i; j ++) {
	resj = sq[j+1];
  	gr->e[i][j][IDX(a,resj,Kg)] += add;
      }
      for (j = i+1; j < L; j ++) {
	resj = sq[j+1];
  	gr->e[j][i][IDX(resj,a,Kg)] += add;
      }
    }

    // one more term to the gradient
    // derivative respect to hi(resi)
    gr->h[i][resi] -= wgt;
    
    // derivative respect to eij(resi,resj)
    for (j = 0; j < i; j ++) {
      resj       = sq[j+1];
      gr->e[i][j][IDX(resi,resj,Kg)] -= wgt;
    }
    for (j = i+1; j < L; j ++) {
      resj       = sq[j+1];
      gr->e[j][i][IDX(resj,resi,Kg)] -= wgt;
     }
    
  } // for all sequences
  
  // l2-regularization
  if (dofunc) 
    nlogp += potts_aplm_regularize_l2(i, pt);
  if (dodfunc) {
    for (a = 0; a < Kg; a++) gr->h[i][a] += pt->muh * 2.0 * pt->h[i][a];
    
    for (j = 0; j < i; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) 
	  gr->e[i][j][IDX(a,b,Kg)] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
    
    for (j = i+1; j < L; j ++)
      for (a = 0; a < Kg; a++)
	for (b = 0; b < Kg; b++) 
	  pt->e[i][j][IDX(a,b,Kg)] += pt->mue * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
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
  double   zi; 
  double   add;
  double   wgt;
  int      L   = msa->alen;
  int      Kg  = pt->Kg;
  int      Kg2 = pt->Kg2;
  int      dim = PLMDIM(L,L,Kg,Kg2);
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
  if (dodfunc) esl_vec_DSet(dnlogp, dim, 0.0); 
  ESL_ALLOC(Hi, sizeof(double)*Kg);
 
  for (i = 0; i < L; i ++) {
    // Initialize
    x0i = PLMDIM(i,L,Kg,Kg2); //\sum_{j=0}^{i-1} [ Kg + Kg2*(L-1-j)] 
    if (dofunc) nlogpi = 0.0;
    
    for (s = 0; s < msa->nseq; s++) {
      // Initialize
      sq   = msa->ax[s];
      resi = sq[i+1];
      wgt  = msa->wgt[s];
      
      // the hamiltonian
      //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
      // and logzi    = log \sum_a exp H^s_i[a]
      // and zi       = \sum_a exp H^s_i[a]
      //logzi = potts_logzi_packed_plm(i, p, L, Kg, sq, Hi);   
      zi = potts_Zi_packed_plm(i, p, L, Kg, Kg2, sq, Hi);   

      // the function
      if (dofunc) nlogpi += wgt * (-Hi[resi] + log(zi));
      
      // the gradient
      if (dodfunc == FALSE) continue;
      for (a = 0; a < Kg; a++) {
	add = wgt * exp(Hi[a]) / zi;

	// derivative respect to hi(a)
	x          = x0i + a;
	dnlogp[x] += add;
	
	// derivative respect to eij(a,b)
	for (j = i+1; j < L; j ++) {
	  resj       = sq[j+1];
	  x          = x0i + PLMIDXR(i,j,L,Kg,Kg2) + IDX(a,resj,Kg);
	  dnlogp[x] += add;
	}
	for (j = 0; j < i; j ++) {
	  resj       = sq[j+1];
	  x          = PLMIDX(j,i,L,Kg,Kg2) + IDX(resj,a,Kg);
	  dnlogp[x] += add;
	}
      }

      // one more term to the gradient
      // derivative respect to hi(resi)
      x          = x0i + resi;
      dnlogp[x] -= wgt;
      
      // derivative respect to eij(resi,resj)
      for (j = i+1; j < L; j ++) {
	resj       = sq[j+1];
	x          = x0i + PLMIDXR(i,j,L,Kg,Kg2) + IDX(resi,resj,Kg);
	dnlogp[x] -= wgt;
      }      
      for (j = 0; j < i; j ++) {
	resj       = sq[j+1];
	x          = PLMIDX(j,i,L,Kg,Kg2) + IDX(resj,resi,Kg);
	dnlogp[x] -= wgt;
      }      

    } // for all sequences
    
    if (dofunc) nlogp += nlogpi;
  
  } // for all positions i

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
  double   zi; 
  double   nlogp;
  double   add;
  double   wgt;
  int      L   = msa->alen;
  int      Kg  = pt->Kg;
  int      Kg2 = pt->Kg2;
  int      dim = APLMDIM(L,Kg,Kg2);
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
  if (dodfunc) esl_vec_DSet(dnlogp, dim, 0.0); 
  ESL_ALLOC(Hi, sizeof(double)*Kg);
 
  for (s = 0; s < msa->nseq; s++) {
    
    // Initialize
    sq   = msa->ax[s];
    resi = sq[i+1];
    wgt  = msa->wgt[s];
   
    // the hamiltonian
    //     H^s_i[a] = hi[a] + \sum_j(\neq i) eij(a,sj)
    // and logzi    = log \sum_a exp H^s_i[a]
    // and zi       = \sum_a exp H^s_i[a]
    //logzi = potts_logzi_packed_aplm(i, p, L, Kg, sq, Hi);   
    zi = potts_Zi_packed_aplm(i, p, L, Kg, Kg2, sq, Hi);   
    
    // the function
    if (dofunc)
      nlogp += wgt * (-Hi[resi] + log(zi));
    
     // the gradient
    if (dodfunc == FALSE) continue;
    for (a = 0; a < Kg; a++) {
      add = wgt * exp(Hi[a]) / zi;
      
      // derivative respect to hi(a)
      dnlogp[a] += add;
     
      // derivative respect to eij(a,b)
      for (j = 0; j < i; j ++) {
	resj       = sq[j+1];
 	x          = APLMIDXL(j,Kg,Kg2) + IDX(a,resj,Kg);
	dnlogp[x] += add;
      }
      for (j = i+1; j < L; j ++) {
	resj       = sq[j+1];
	x          = APLMIDXG(j,Kg,Kg2) + IDX(a,resj,Kg);
	dnlogp[x] += add;
      }
    }

    // one more term to the gradient
    // derivative respect to hi(resi)
    dnlogp[resi] -= wgt;
    
    // derivative respect to eij(resi,resj)
    for (j = 0; j < i; j ++) {
      resj       = sq[j+1];
      x          = APLMIDXL(j,Kg,Kg2) + IDX(resi,resj,Kg);
      dnlogp[x] -= wgt;
    }
    for (j = i+1; j < L; j ++) {
      resj       = sq[j+1];
      x          = APLMIDXG(j,Kg,Kg2) + IDX(resi,resj,Kg);
      dnlogp[x] -= wgt;
    }      
    
  } // for all sequences
  
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
  int    Kg  = pt->Kg;
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
potts_Hi_APLM_Packed(int i, int a, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq)
{
  double val = 0.;
  int    resj;
  int    j;
  int    x;
  
  val += p[a];
 
  for (j = 0; j < i; j++) {
    resj  = sq[j+1];	  
    x     = APLMIDXL(j,Kg,Kg2) + IDX(a,resj,Kg);
    val  += p[x];
  }
  for (j = i+1; j < L; j++) {
    resj  = sq[j+1];	  
    x     = APLMIDXG(j,Kg,Kg2) + IDX(a,resj,Kg);
    val  += p[x];
  }
  return val;
}

//exp[ hi(a) + \sum_{j\neq i} eij(a,sqj)
double
potts_Hi_PLM_Packed(int i, int a, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq)
{
  double val = 0.;
  int    x0i = PLMDIM(i,L,Kg,Kg2);
  int    resj;
  int    j;
  int    x;
  
  x    = x0i + a;
  val += p[x];
  
  for (j = i+1; j < L; j++) {
    resj  = sq[j+1];
    x     = x0i + PLMIDXR(i,j,L,Kg,Kg2) + IDX(a,resj,Kg);
    val  += p[x];
  }
  for (j = 0; j < i; j++) {
    resj  = sq[j+1];
    x     = PLMIDX(j,i,L,Kg,Kg2) + IDX(resj,a,Kg);
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
  int     Kg    = pt->Kg;
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
double
potts_Zi(int i, PT *pt, ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  Zi = 0.;
  int     Kg = pt->Kg;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi(i, a, pt, sq);
    if (Hi) Hi[a] = Ha;
     
    Zi += exp(Ha);
  }
  return Zi;

 ERROR:
  return status;
}

static double
potts_logzi_packed_plm(int i, double *p, int L, int Kg, int Kg2,  ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  logzi = -eslINFINITY;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi_PLM_Packed(i, a, p, L, Kg, Kg2, sq);
    if (Hi) Hi[a] = Ha;
     
    logzi = e2_DLogsum(logzi, Ha);
  }
  return logzi;

 ERROR:
  return status;
}
static double
potts_Zi_packed_plm(int i, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  Zi = 0.0;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi_PLM_Packed(i, a, p, L, Kg, Kg2, sq);
    if (Hi) Hi[a] = Ha;
     
    Zi += exp(Ha);
  }
  return Zi;

 ERROR:
  return status;
}

static double
potts_logzi_packed_aplm(int i, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  logzi = -eslINFINITY;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi_APLM_Packed(i, a, p, L, Kg, Kg2, sq);
    if (Hi) Hi[a] = Ha;
     
    logzi = e2_DLogsum(logzi, Ha);
  }
  return logzi;

 ERROR:
  return status;
}
static double
potts_Zi_packed_aplm(int i, double *p, int L, int Kg, int Kg2, ESL_DSQ *sq, double *Hi)
{
  double  Ha;
  double  Zi = 0.0;
  int     a;
  int     status;

  for (a = 0; a < Kg; a++) {
    Ha = potts_Hi_APLM_Packed(i, a, p, L, Kg, Kg2, sq);
    if (Hi) Hi[a] = Ha;
     
    Zi += exp(Ha);
  }
  return Zi;

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
  int              K  = pt->abc->K;
  int              Kg = pt->Kg;
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
  int              Kg = pt->Kg;
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
  double        H = 0.;
  int           L = pt->L;
  int           Kg = pt->Kg;
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
  int     Kg = pt->Kg;
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
  int    Kg = pt->Kg;
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
  int    Kg  = pt->Kg;
  int    Kg2 = pt->Kg2;
  int    x0i, x;
  int    i, j;
  int    a, b;
  
  // l2-regularization
  for (i = 0; i < L; i ++) {
    x0i = PLMDIM(i,L,Kg,Kg2);
    for (a = 0; a < Kg; a ++) {
      hi   = p[x0i+a];
      reg += pt->muh * hi * hi;
  
      for (j = i+1; j < L; j ++) 
	for (b = 0; b < Kg; b ++) {
	  x    = x0i + PLMIDXR(i,j,L,Kg,Kg2) + IDX(a,b,Kg);
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
  int    Kg = pt->Kg;
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
  int    Kg = pt->Kg;
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
  int    Kg  = pt->Kg;
  int    Kg2 = pt->Kg2;
  int    x;
  int    j;
  int    a, b;
  
  // l2-regularization
  for (a = 0; a < Kg; a ++) {
    hi   = p[a];
    reg += pt->muh * hi * hi;
    
    for (j = 0; j < i; j ++) 
      for (b = 0; b < Kg; b ++) {
	x    = APLMIDXL(j,Kg,Kg2) + IDX(a,b,Kg);
	eij  = p[x];
	reg += pt->mue * eij * eij;
      }  
    for (j = i+1; j < L; j ++) 
      for (b = 0; b < Kg; b ++) {
	x    = APLMIDXG(j,Kg,Kg2) + IDX(a,b,Kg);
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
  int    Kg = pt->Kg;
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
