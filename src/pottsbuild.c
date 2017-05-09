/* pottsbuild.c */

#include "rscape_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "minimize.h"
#include "pottsbuild.h"
#include "pottsscore.h"

static int    optimize_pack_paramvector   (double *p, long np, struct optimize_data *data);
static int    optimize_unpack_paramvector (double *p, long np, struct optimize_data *data);
static void   optimize_bracket_define_direction(double *p, long np, struct optimize_data *data);
static double optimize_potts_func              (double *p, long np, void *dptr);
static double func_potts_full                  (PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
static double func_potts_aplm                  (PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);


PT *
potts_Create(int64_t L, ESL_ALPHABET *abc, double mu, POTTSTRAIN pottstrain)
{
  PT  *pt = NULL;
  int  K  = abc->K;
  int  K2 = K * K;
  int  i, j;
  int  status;
  
  ESL_ALLOC(pt, sizeof(PT));
  pt->L     = L;
  pt->abc   = abc;
  pt->mu    = mu;
  pt->train = pottstrain;

  ESL_ALLOC(pt->e,            sizeof(double **) * L);
  ESL_ALLOC(pt->h,            sizeof(double  *) * L);
  for (i = 0; i < L; i++) {
    ESL_ALLOC(pt->e[i],       sizeof(double  *) * L);
    ESL_ALLOC(pt->h[i],       sizeof(double   ) * K2);
    for (j = 0; j < L; j++) {
       ESL_ALLOC(pt->e[i][j], sizeof(double   ) * K2);
    }
  }
   
  /* initialize for adding counts */
  for (i = 0; i < L; i++) {
    esl_vec_DSet(pt->h[i], K, 0.0); 
 
    for (j = 0; j < L; j++) 
      esl_vec_DSet(pt->e[i][j], K2, 0.0); 
  }

  return pt;
  
 ERROR:
  return NULL;
}

int
potts_Assign(PT *pt, double val)
{
  int L = pt->L;
  int K = pt->abc->K;
  int i, j;
  int a, b;
  
  for (i = 0; i < L; i++) {
    for (a = 0; a < K; a ++) pt->h[i][a] = val;
    for (j = i; j < L; j++)
      for (a = 0; a < K; a ++) 
	for (b = 0; b < K; b ++) {
	  pt->e[i][j][IDX(a,b,K)] = ((i==j||j>i+1)&&a+b==K)? 0.:val;
	  pt->e[j][i][IDX(b,a,K)] = pt->e[i][j][IDX(a,b,K)];
      }
  }
  return eslOK;
}

void                
potts_Destroy(PT *pt)
{
  int i, j;

  if (pt) {
    for (i = 0; i < pt->L; i++) {
      for (j = 0; j < pt->L; j++) 
	free(pt->e[i][j]);
     free(pt->e[i]);
     free(pt->h[i]);
    }
    free(pt->e);
    free(pt->h);
    free(pt);
  }
}

int
potts_Build(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  float    firststep;
  int      status;

  /* init */
  firststep = 1.0;

  if (pt == NULL) ESL_XFAIL(eslFAIL, errbuf, "error in potts_Build()");
  potts_Assign(pt, 0.1);

  status = potts_OptimizeGD(pt, msa, firststep, tol, errbuf, verbose);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "error optimizing potts");
    
 return eslOK;

 ERROR:
 return status;
}

int
potts_OptimizeGD(PT *pt, ESL_MSA *msa, float firststep, float tol, char *errbuf, int verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 logp;
  int                    L = msa->alen;
  int                    K = msa->abc->K;
  int                    np;
  int                    status;

  np = L*K*(1+L*K);     /* the variables hi eij */
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);
  
 /* Copy shared info into the "data" structure
   */
  data.pt         = pt;
  data.msa        = msa;
  data.firststep  = firststep;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.verbose    = verbose;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, (long)np, &data);
 
  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (long)np, &data);
  status = min_Bracket(p, u, np, data.firststep,
		       &optimize_potts_func,
		       (void *) (&data), 
		       tol, wrk, &logp);
  if (status != eslOK) 
    esl_fatal("optimize_potts(): bad bracket minimization");	
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  if (verbose) printf("END POTTS OPTIMIZATION\n");
  status = potts_GaugeZeroSum(pt, tol, errbuf);
  if (status != eslOK) esl_fatal("%s.optimize_potts(): bad zero gauge sum", errbuf);
    
  /* clean up */
  if (u   != NULL) free(u);
  if (p   != NULL) free(p);
  if (wrk != NULL) free(wrk);
  return eslOK;

 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;
}


// sum_b eij(a,b) = sum_a eij(a,b) = sum_a hi[a] = 0
//
// e(a,b) → e(a,b) − e(·,b) − e(a,·) + e(·,·)
// h(a)   → h(a)   - e(a,·)
int
potts_GaugeZeroSum(PT *pt, double tol, char *errbuf)
{
  int      K = pt->abc->K;
  double   sumi[K];
  double   sumj[K];
  double **sumh;
  double   sum;
  int      i, j;
  int      a, b;
  int      status;

  ESL_ALLOC(sumh,    sizeof(double *) * pt->L);
  ESL_ALLOC(sumh[0], sizeof(double)   * pt->L*K);
  for (i = 1; i < pt->L; i++) sumh[i] = sumh[i-1] + K;
  
  for (i = 0; i < pt->L; i++) {
    
    for (a = 0; a < K; a ++) sumh[i][a] = 0;
    
    for (j = 0; j < pt->L; j++) {
      
      sum = 0;     
      for (a = 0; a < K; a ++) {
	sumi[a] = 0;
	sumj[a] = 0; 
	for (b = 0; b < K; b ++) {
	  sum        += pt->e[i][j][IDX(a,b,K)];
	  sumi[a]    += pt->e[i][j][IDX(a,b,K)];
	  sumj[a]    += pt->e[i][j][IDX(b,a,K)];
	  sumh[i][a] += pt->e[i][j][IDX(a,b,K)];
	}
      }
      
      for (a = 0; a < K; a ++) 
	for (b = 0; b < K; b ++)
	  pt->e[i][j][IDX(a,b,K)] -= sumi[b]/K + sumj[a]/K - sum/K;
    }
    for (a = 0; a < K; a ++) pt->h[i][a] -= sumh[i][a]; 
  }
  
#if 0
  // check it is a zero sum
  for (i = 0; i < pt->L; i++) {
    sum = 0;
    for (a = 0; a < K; a ++) sum += pt->h[i][a];
    if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for h failed. sum = %f", sum);
    
    for (j = i+1; j < pt->L; j++) {
      sum = 0;
      for (a = 0; a < K; a ++)
	for (b = 0; b < K; b ++) 
	  sum += pt->e[i][j][IDX(a,b,K)];
      if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for eij failed");

      sum = 0;
      for (a = 0; a < K; a ++)
	for (b = 0; b < K; b ++) 
	  sum += pt->e[i][j][IDX(b,a,K)];
      if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for eij failed");
    }
  }
  #endif
      
  free(sumh[0]);
  free(sumh);
  return eslOK;

 ERROR:
  if (sumh[0]) free(sumh[0]);
  if (sumh)    free(sumh);
  return status;
}


/*------------------------------- internal functions ----------------------------------*/

static int
optimize_pack_paramvector(double *p, long np, struct optimize_data *data)
{
  int   L = data->msa->alen;
  int   K = data->msa->abc->K;
  int   x = 0;
  int   i, j;
  int   a;

  for (i = 0; i < L; i++) {
    for (a = 0; a < K; a++)     p[x++] = data->pt->h[i][a];
    for (j = 0; j < L; j++) 
      for (a = 0; a < K*K; a++) p[x++] = data->pt->e[i][j][a];
  }
  return eslOK;  
}

static int
optimize_unpack_paramvector(double *p, long np, struct optimize_data *data)
{
  int   L = data->msa->alen;
  int   K = data->msa->abc->K;
  int   x = 0;
  int   i, j;
  int   a;
  
  for (i = 0; i < L; i++) {
    for (a = 0; a < K; a++)     data->pt->h[i][a]    = p[x++];
    for (j = 0; j < L; j++) 
      for (a = 0; a < K*K; a++) data->pt->e[i][j][a] = p[x++];
  }

 return eslOK;
}


static void
optimize_bracket_define_direction(double *u, long np, struct optimize_data *data)
{
  int x;
  for (x = 0; x < np; x++) u[x] = 0.1;
  u[np] = 0.5;
}

static double
optimize_potts_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  
  optimize_unpack_paramvector(p, np, data);
  
  switch(data->pt->train) {
  case FULL:
    data->logp = func_potts_full(data->pt, data->msa, data->tol, data->errbuf, data->verbose);
    break;
  case APLM:
    data->logp = func_potts_aplm(data->pt, data->msa, data->tol, data->errbuf, data->verbose);
    break;
  case GINV:
  case ACE:
  case BML:
  case NONE:
    data->logp = -eslINFINITY;
    break;
  }
  
  return (double)-data->logp;
}

static double
func_potts_full(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_FULLLogp(pt, msa, &logp, errbuf, verbose);
  if (status != eslOK) exit(1);
  
#if 1
  printf("logp FULL %f\n", logp);
#endif
  
  return -logp;
}
 
static double
func_potts_aplm(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_APLMLogp(pt, msa, &logp, errbuf, verbose);

  if (status != eslOK) exit(1);
  
#if 1
  printf("logp APLM %f\n", logp);
#endif
  
  return -logp;
}
