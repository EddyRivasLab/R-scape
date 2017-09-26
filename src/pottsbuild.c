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
#include "esl_fileparser.h"
#include "esl_minimizer.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include <lbfgs.h>

#include "logsum.h"
#include "minimize.h"
#include "pottsbuild.h"
#include "pottsscore.h"
#include "correlators.h"

static int             optimize_all_pack_paramvector         (double *p,          int np, struct optimize_data *data);
static int             optimize_all_unpack_paramvector       (double *p,          int np, struct optimize_data *data);
static int             optimize_aplm_pack_paramvector        (double *p,          int np, struct optimize_data *data);
static int             optimize_aplm_unpack_paramvector      (double *p,          int np, struct optimize_data *data);
static int             optimize_aplm_lbfgs_pack_paramvector  (lbfgsfloatval_t *p, int np, struct optimize_data *data);
static int             optimize_aplm_lbfgs_unpack_paramvector(lbfgsfloatval_t *p, int np, struct optimize_data *data);
static void            optimize_bracket_define_direction     (double *p,          int np, struct optimize_data *data);
static double          optimize_potts_func                   (double *p,          int np, void *dptr);
static void            optimize_potts_dfunc_aplm             (double *p,          int np, void *dptr, double *dx);
static double          logp_potts_ml                         (PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
static double          logp_potts_plm                        (PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
static double          logp_potts_aplm              (int pos, PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
static int             symmetrize                            (PT *pt);
static int             progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
				const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
				int n, int k, int ls);
static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

PT *
potts_Build(ESL_RANDOMNESS *r, ESL_MSA *msa, double ptmu, PTTRAIN pttrain, PTSCTYPE ptsctype, FILE *pottsfp, float tol, char *errbuf, int verbose)
{
  PT              *pt = NULL;
  int              status;

  tol = 0.1;

  pt = potts_Create(msa->alen, msa->abc->K+1, msa->abc, ptmu, pttrain, ptsctype);
  if (pt == NULL) ESL_XFAIL(eslFAIL, errbuf, "error creating potts");

  // Initialize 
  //status = potts_AssignGaussian(r, pt, 0., 1.0);
  status = potts_AssignGT(r, msa, pt, tol, errbuf, verbose);
  if (status != eslOK) { printf("%s\n", errbuf); goto ERROR; }

  status = potts_GaugeZeroSum(pt, errbuf, verbose);
  if (status != eslOK) { printf("%s\n", errbuf); goto ERROR; }
  
  /* init */
  switch (pt->train) {
  case NONE:
   ESL_XFAIL(eslFAIL, errbuf, "error, you should not be here");
     break;
  case ML:
  case PLM:
    status = potts_OptimizeCGD_ALL(pt, msa, tol, errbuf, verbose);
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "error all optimizing potts");
    break;
  case APLM:
    status = potts_OptimizeCGD_APLM(pt, msa, tol, errbuf, verbose);
    //status = potts_OptimizeLBFGS_APLM(pt, msa, tol, errbuf, verbose);
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "error aplm optimizing potts");
    break;
  case GINV:
  case ACE:
  case BML:
    ESL_XFAIL(eslFAIL, errbuf, "error optimization method not implemented");
    break;
  }

  if (pottsfp) potts_Write(pottsfp, pt);
 return pt;

 ERROR:
 return NULL;
}

int
potts_OptimizeCGD_ALL(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 logp;
  int                    L = msa->alen;
  int                    Kg = msa->abc->K+1;
  int                    np;
  int                    status;

  np = L*Kg*(1+L*Kg);     /* the variables hi eij */
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);
  
 /* Copy shared info into the "data" structure
   */
  data.pt         = pt;
  data.msa        = msa;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.verbose    = verbose;
 
  /* Create the parameter vector.
   */
  optimize_all_pack_paramvector(p, (int)np, &data);
 
  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (int)np, &data);
  status = esl_min_ConjugateGradientDescent(p, u, np, 
					    &optimize_potts_func, NULL,
					    (void *) (&data), 
					    tol, wrk, &logp);
  if (status != eslOK) 
    esl_fatal("optimize_potts(): esl_min_ConjugateGradientDescent failed");	
  
  /* unpack the final parameter vector */
  optimize_all_unpack_paramvector(p, (int)np, &data);
  if (verbose) printf("END POTTS ALL OPTIMIZATION\n");

  if (verbose) potts_Write(stdout, pt);

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


int
potts_OptimizeCGD_APLM(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 logp;
  int                    L = msa->alen;
  int                    Kg = msa->abc->K+1;
  int                    np;
  int                    i;
  int                    status;

  np = Kg*(1+L*Kg);     /* the variables hi eij */
  
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);

  for (i = 0; i < L; i ++) {
    /* Copy shared info into the "data" structure
     */
    data.pt         = pt;
    data.msa        = msa;
    data.pos        = i;
    data.tol        = tol;
    data.errbuf     = errbuf;
    data.verbose    = verbose;
    
    /* Create the parameter vector.
     */
    optimize_aplm_pack_paramvector(p, (int)np, &data);
    
    /* pass problem to the optimizer
     */
    optimize_bracket_define_direction(u, (int)np, &data);
    status = esl_min_ConjugateGradientDescent(p, u, np, 
					      &optimize_potts_func, &optimize_potts_dfunc_aplm,
					      (void *) (&data), 
					      tol, wrk, &logp);
    if (status != eslOK) 
      esl_fatal("optimize_potts(): esl_min_ConjugateGradientDescent failed");	
    
    /* unpack the final parameter vector */
    optimize_aplm_unpack_paramvector(p, (int)np, &data);
    if (1||verbose) printf("END POTTS CGD APLM OPTIMIZATION for position %d\n", i);
  }
  if (verbose) printf("END POTTS CGD APLM OPTIMIZATION\n");

  // symmetrize
  symmetrize(pt);
 
  if (verbose) potts_Write(stdout, pt);

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

int
potts_OptimizeLBFGS_APLM(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  struct optimize_data   data;
  int                    L  = msa->alen;
  int                    Kg = msa->abc->K+1;
  int                    np = Kg*(1+L*Kg);
  lbfgsfloatval_t       *x = lbfgs_malloc(np+1);
  lbfgsfloatval_t        fx;
  lbfgs_parameter_t      param;
  int                    i;
  int                    ret = 0;
  int                    status;
  
  if (x == NULL) { printf("ERROR: Failed to allocate a memory block for variables.\n"); return eslFAIL; }
  
  /* Copy shared info into the "data" structure
   */
  data.pt         = pt;
  data.msa        = msa;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.verbose    = verbose;
  
  for (i = 0; i < L; i ++) {

    data.pos = i;
    
    // Initialize the variables.
    optimize_aplm_lbfgs_pack_paramvector(x, (int)np, &data);

    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    param.m              = 3;
    param.max_iterations = 50;
    param.epsilon        = 0.01;

    // The L-BFGS optimization
    ret = lbfgs(np, x, &fx, evaluate, progress, (void *)(&data), &param);
    if (ret == LBFGS_ALREADY_MINIMIZED) { printf("LBFGS_ALREADY_MINIMIZED \n"); }
    if (ret < 0) { printf("LBFGS failed with code %d\n", ret); exit(1); }
    
    /* Bail out if the function is now +/-inf: this can happen if the caller
     * has screwed something up.
     */
    if (fx == eslINFINITY || fx == -eslINFINITY) ESL_EXCEPTION(eslERANGE, "minimum not finite");

    // recover the variable into data->pt
    optimize_aplm_lbfgs_unpack_paramvector(x, (int)np, &data);
    
    if (1||verbose) printf("END POTTS LBFGS APLM OPTIMIZATION for position %d\n", i);
  }
  if (verbose) printf("END POTTS LBFGS APLM OPTIMIZATION\n");

  // symmetrize
  symmetrize(pt);
 
  if (verbose) potts_Write(stdout, pt);
  
  lbfgs_free(x); x = NULL;
  return eslOK;

 ERROR:
  if (x) lbfgs_free(x);
  return status;
}

PT *
potts_Create(int64_t L, int Kg, ESL_ALPHABET *abc, double mu, PTTRAIN pttrain, PTSCTYPE  ptsctype)
{
  PT  *pt  = NULL;
  int  Kg2 = Kg * Kg;
  int  i, j;
  int  status;

  if (abc && Kg != abc->K+1) return NULL;
  
  ESL_ALLOC(pt, sizeof(PT));
  pt->L      = L;
  pt->abc    = abc;
  pt->mu     = mu;
  pt->train  = pttrain;
  pt->sctype = ptsctype;

  ESL_ALLOC(pt->e,            sizeof(double **) * L);
  ESL_ALLOC(pt->h,            sizeof(double  *) * L);
  for (i = 0; i < L; i++) {
    ESL_ALLOC(pt->e[i],       sizeof(double  *) * L);
    ESL_ALLOC(pt->h[i],       sizeof(double   ) * Kg2);
    for (j = 0; j < L; j++) {
       ESL_ALLOC(pt->e[i][j], sizeof(double   ) * Kg2);
    }
  }
   
  /* initialize for adding counts */
  for (i = 0; i < L; i++) {
    esl_vec_DSet(pt->h[i], Kg, 0.0); 
 
    for (j = 0; j < L; j++) 
      esl_vec_DSet(pt->e[i][j], Kg2, 0.0); 
  }

  return pt;
  
 ERROR:
  return NULL;
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
potts_AssignGaussian(ESL_RANDOMNESS *r, PT *pt, double mu, double sigma)
{
  int L  = pt->L;
  int Kg = pt->abc->K+1;
  int i, j;
  int a, b;
  
  for (i = 0; i < L; i++) {
    for (a = 0; a < Kg; a ++)
      pt->h[i][a] = esl_rnd_Gaussian(r, mu, sigma);
    
    for (j = i; j < L; j++)
      for (a = 0; a < Kg; a ++) 
	for (b = 0; b < Kg; b ++) {
	  pt->e[i][j][IDX(a,b,Kg)] = (i==j)? 0 : esl_rnd_Gaussian(r, mu, sigma);
	  pt->e[j][i][IDX(b,a,Kg)] = pt->e[i][j][IDX(a,b,Kg)];
      }
  }

  return eslOK;
}


int
potts_AssignGT(ESL_RANDOMNESS *r, ESL_MSA *msa, PT *pt, float tol, char *errbuf, int verbose)
{
  struct mutual_s *mi = NULL;
  double           exp;
  double           obs;
  double           gt;
  double           pseudoc = 0.1;
  double           hi;
  int              L  = pt->L;
  int              K  = pt->abc->K;
  int              Kg = K+1;
  int              i, j;
  int              a, b;
  int              status;

  mi = corr_Create(L, msa->nseq, FALSE, 0, 0, pt->abc, C16);
  if (mi == NULL) ESL_XFAIL(eslFAIL, errbuf, "could not create mi");

  corr_ReuseCOV(mi, GT, C16);
  status = corr_Probs(r, msa, NULL, NULL, mi, POTTS, FALSE, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  for (i = 0; i < L; i++) {
    for (j = i; j < L; j++)
      for (a = 0; a < K; a ++) 
	for (b = 0; b < K; b ++) {

	  exp = mi->nseff[i][j] * mi->pm[i][a] * mi->pm[j][b];
	  obs = mi->nseff[i][j] * mi->pp[i][j][IDX(a,b,Kg)];
	  gt  = (exp > 0. && obs > 0.) ? obs * log (obs / exp) : pseudoc;
	  
	  pt->e[i][j][IDX(a,b,Kg)] = (i==j)? 0 : gt;
	  pt->e[j][i][IDX(b,a,Kg)] = pt->e[i][j][IDX(a,b,Kg)];
      }
  }
  
  for (i = 0; i < L; i++) 
    for (a = 0; a < Kg; a ++) 
     pt->h[i][a] = 0.;

  //potts_Write(stdout, pt);

  corr_Destroy(mi);
  return eslOK;

 ERROR:
  if (mi) corr_Destroy(mi);
  return status;
}

// sum_b eij(a,b) = sum_a eij(a,b) = sum_a hi[a] = 0
//
// eij(a,b) → eij(a,b) − 1/K * \sum_c eij(c,b) − 1/k * \sum_d eij(a,d) + 1/K^2 * \sum_c \sum_d eij(c,d)
// hi(a)    → hi(a)    + \sum_{j\neq i} 1/K * \sum_d eij(a,d) - 1/K \sum_c hi(c) - 1/K^2 * \sum_c \sum_d eij(c,d)
int
potts_GaugeZeroSum(PT *pt, char *errbuf, int verbose)
{
  int      Kg = pt->abc->K+1;
  int      Kg2 = Kg*Kg;
  double   suma;
  double   sumi[Kg];
  double   sumj[Kg];
  double   sumhh;
  double   sumh;
  double   sum;
  int      i, j;
  int      a, b;
  int      status;

  for (i = 0; i < pt->L; i++) {
    
    suma = 0;
    for (a = 0; a < Kg; a ++)
      suma += pt->h[i][a];
    
    sumhh = 0;
    for (a = 0; a < Kg; a ++) {
      
      sumh = 0.0;      
      for (j = 0; j < pt->L; j++) {
	if (i==j) continue;
	for (b = 0; b < Kg; b ++) {
	  sumh  += pt->e[i][j][IDX(a,b,Kg)];
	  sumhh += pt->e[i][j][IDX(a,b,Kg)];
	}
      }
      
      pt->h[i][a] += sumh/Kg - suma/Kg;    
    }
    
    for (a = 0; a < Kg; a ++)
      pt->h[i][a] -= sumhh/Kg2;  
  }
  
  for (i = 0; i < pt->L; i++)    
    for (j = 0; j < pt->L; j++) {
      
      sum = 0;     
      for (a = 0; a < Kg; a ++) {
	sumi[a] = 0;
	sumj[a] = 0; 
	for (b = 0; b < Kg; b ++) {
	  sum        += pt->e[i][j][IDX(a,b,Kg)];
	  sumi[a]    += pt->e[i][j][IDX(b,a,Kg)];
	  sumj[a]    += pt->e[i][j][IDX(a,b,Kg)];
	}
      }
      
      for (a = 0; a < Kg; a ++) 
	for (b = 0; b < Kg; b ++) {
	  pt->e[i][j][IDX(a,b,Kg)] -= sumi[b]/Kg + sumj[a]/Kg - sum/Kg2;
	  pt->e[j][i][IDX(a,b,Kg)]  = pt->e[i][j][IDX(a,b,Kg)];
	}
    } 
  
  
#if 0
  double tol = 1e-5;
  // check it is a zero sum
  for (i = 0; i < pt->L; i++) 
    for (j = 0; j < pt->L; j++) {
      
      for (a = 0; a < Kg; a ++) {
	sum = 0;	
	for (b = 0; b < Kg; b ++) 
	  sum += pt->e[i][j][IDX(a,b,Kg)];
	if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for eij failed sumj %f | %d %d %d", sum, i, j, a);
      }
      
      for (a = 0; a < Kg; a ++) {
	sum = 0;
	for (b = 0; b < Kg; b ++) 
	  sum += pt->e[i][j][IDX(b,a,Kg)];
	if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for eij failed sumi %f | %d %d %d", sum, i, j, a);
      }
      
    } // for all i,j
  
  for (i = 0; i < pt->L; i++) {
    sum = 0;
    for (a = 0; a < Kg; a ++) sum += pt->h[i][a];
    if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for h failed. sum = %f", sum);
  }
  
#endif
  
  return eslOK;

 ERROR:
  return status;
}


PT *
potts_Read(char *paramfile, ESL_ALPHABET *abc, char *errbuf)
{
  ESL_FILEPARSER  *efp   = NULL;
  char            *tok;
  PT              *pt = NULL;
  double           hi;
  int              ln = 0;
  int              Kg = 0;
  int              L  = 0;
  int              t  = 0;
  int              i, j;
  int              status;

  // One pass to figure out L
  if (esl_fileparser_Open(paramfile, NULL, &efp) != eslOK)  ESL_XFAIL(eslFAIL, errbuf, "file open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      t = 0;
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", paramfile);
      i = atoi(tok);
	
      while (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) == eslOK) {
	hi = atof(tok);
	t ++;
      }
      if (ln == 0) Kg = t;
      
      if (t  == Kg) ln ++;
      else { L = ln; break; }
    }
  esl_fileparser_Close(efp);
  printf("L %d K %d\n", L, Kg);

  pt = potts_Create(L, Kg, NULL, 0.0, NONE, SCNONE);
  if (abc) {
    if (Kg != abc->K+1) ESL_XFAIL(eslFAIL, errbuf, "wrong alphabet for file %s", paramfile);
    pt->abc = abc;
  }
  
  // now the real parsing
  if (esl_fileparser_Open(paramfile, NULL, &efp) != eslOK)  ESL_XFAIL(eslFAIL, errbuf, "file open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  
  ln = 0;
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      t = 0;
      if (ln < L) {
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", paramfile);
	i = atoi(tok);
	
	while (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) == eslOK) {
	  pt->h[ln][t] = atof(tok);
	  t ++;
	}
	ln ++;
      }
      else {
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", paramfile);
	i = atoi(tok);
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", paramfile);
	j = atoi(tok);
	
	while (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) == eslOK) {
	  pt->e[i][j][t] = atof(tok);
	  t ++;
	}
      }
    }
  esl_fileparser_Close(efp);
  
  return pt;
  
 ERROR:
  return NULL;
}

void
potts_Write(FILE *fp, PT *pt)
{
  int L = pt->L;
  int Kg = pt->abc->K+1;
  int i, j;
  int a, b;
  
  for (i = 0; i < L; i++) {
    fprintf(fp, "%d ", i);
    for (a = 0; a < Kg; a ++) fprintf(fp, "%f ", pt->h[i][a]);
    fprintf(fp, "\n");
  }
  
  for (i = 0; i < L; i++) 
    for (j = i+1; j < L; j++) {
      fprintf(fp, "%d %d ", i, j);
      for (a = 0; a < Kg; a ++) 
	for (b = 0; b < Kg; b ++) {
	  fprintf(fp, "%f ", pt->e[i][j][IDX(a,b,Kg)]);
	}
      fprintf(fp, "\n");
    }
}


/*------------------------------- internal functions ----------------------------------*/

static int
optimize_all_pack_paramvector(double *p, int np, struct optimize_data *data)
{
  int   L   = data->msa->alen;
  int   Kg  = data->msa->abc->K+1;
  int   Kg2 = Kg*Kg;
  int   x   = 0;
  int   i, j;
  int   a;

  for (i = 0; i < L; i++) {
    for (a = 0; a < Kg; a++)      p[x++] = data->pt->h[i][a];
    for (j = 0; j < L; j++) 
      for (a = 0; a < Kg2; a++) { p[x++] = data->pt->e[i][j][a]; p[x++] = data->pt->e[j][i][a]; }
  }
  return eslOK;  
}

static int
optimize_aplm_pack_paramvector(double *p, int np, struct optimize_data *data)
{
  int   L   = data->msa->alen;
  int   Kg  = data->msa->abc->K+1;
  int   Kg2 = Kg*Kg;
  int   x   = 0;
  int   i   = data->pos;
  int   j;
  int   a;

  for (a = 0; a < Kg; a++)    p[x++] = data->pt->h[i][a];
  for (j = 0; j < L; j++) {
    for (a = 0; a < Kg2; a++) p[x++] = data->pt->e[i][j][a];
  }
  return eslOK;  
}
static int
optimize_aplm_lbfgs_pack_paramvector(lbfgsfloatval_t *p, int np, struct optimize_data *data)
{
  int   L   = data->msa->alen;
  int   Kg  = data->msa->abc->K+1;
  int   Kg2 = Kg*Kg;
  int   x   = 0;
  int   i   = data->pos;
  int   j;
  int   a;

  for (a = 0; a < Kg; a++)    p[x++] = (lbfgsfloatval_t)data->pt->h[i][a];
  for (j = 0; j < L; j++) {
    for (a = 0; a < Kg2; a++) p[x++] = (lbfgsfloatval_t)data->pt->e[i][j][a];
  }
  return eslOK;  
}

static int
optimize_all_unpack_paramvector(double *p, int np, struct optimize_data *data)
{
  double eij, eji;
  int    L   = data->msa->alen;
  int    Kg  = data->msa->abc->K+1;
  int    Kg2 = Kg*Kg;
  int    x   = 0;
  int    i, j;
  int    a;
  
  for (i = 0; i < L; i++) {
    for (a = 0; a < Kg; a++)   data->pt->h[i][a]    = p[x++];
    for (j = 0; j < L; j++) 
      for (a = 0; a < Kg2; a++) {
	eij = p[x++];
	eji = p[x++];
	data->pt->e[i][j][a] = data->pt->e[j][i][a] = 0.5*(eij+eji);
      }
  }

 return eslOK;
}

static int
optimize_aplm_unpack_paramvector(double *p, int np, struct optimize_data *data)
{
  int    L   = data->msa->alen;
  int    Kg  = data->msa->abc->K+1;
  int    Kg2 = Kg*Kg;
  int    x   = 0;
  int    i   = data->pos;
  int    j;
  int    a;
  
  for (a = 0; a < Kg; a++)      data->pt->h[i][a]    = p[x++];
    for (j = 0; j < L; j++) {
      for (a = 0; a < Kg2; a++) data->pt->e[i][j][a] = p[x++];
  }
  
  return eslOK;
}

static int
optimize_aplm_lbfgs_unpack_paramvector(lbfgsfloatval_t *p, int np, struct optimize_data *data)
{
  int    L   = data->msa->alen;
  int    Kg  = data->msa->abc->K+1;
  int    Kg2 = Kg*Kg;
  int    x   = 0;
  int    i   = data->pos;
  int    j;
  int    a;
  
  for (a = 0; a < Kg; a++)      data->pt->h[i][a]    = (double)p[x++];
    for (j = 0; j < L; j++) {
      for (a = 0; a < Kg2; a++) data->pt->e[i][j][a] = (double)p[x++];
  }
  
  return eslOK;
}


static void
optimize_bracket_define_direction(double *u, int np, struct optimize_data *data)
{
  int x;
  for (x = 0; x < np; x++) u[x] = 10.1;
  u[np] = 0.1;
}

static double
optimize_potts_func(double *p, int np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  
  
  switch(data->pt->train) {
  case ML:
    optimize_all_unpack_paramvector(p, np, data);
    data->minfunc = -logp_potts_ml(data->pt, data->msa, data->tol, data->errbuf, data->verbose);
    break;
  case PLM:
    optimize_all_unpack_paramvector(p, np, data);
    data->minfunc = -logp_potts_plm(data->pt, data->msa, data->tol, data->errbuf, data->verbose);
    break;
  case APLM:
    optimize_aplm_unpack_paramvector(p, np, data);
    data->minfunc = -logp_potts_aplm(data->pos, data->pt, data->msa, data->tol, data->errbuf, data->verbose);
    break;
  case GINV:
  case ACE:
  case BML:
  case NONE:
    data->minfunc = -eslINFINITY;
    break;
  }
  
  return (double)data->minfunc;
}

static void
optimize_potts_dfunc_aplm(double *p, int np, void *dptr, double *dx)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *sq;
  double                tmp;
  double                logadd;
  double                val;
  PT                   *pt  = data->pt;
  int                   nsq = data->msa->nseq;
  int                   L   = data->msa->alen;
  int                   Kg  = data->msa->abc->K+1;
  int                   Kg2 = Kg*Kg;
  int                   i   = data->pos;
  int                   x   = 0;
  int                   j;
  int                   resi, resj;
  int                   a, b;
  int                   s;

  optimize_aplm_unpack_paramvector(p, np, data);
  
  for (a = 0; a < Kg; a++) {
    tmp    = 0.;
    logadd = -eslINFINITY;
    
    for (s = 0; s < nsq; s++) {
      sq   = data->msa->ax[s];
      resi = sq[i+1];
      if (resi == a) tmp -= 1.;
      
      val = potts_APLMLognum(i, a, pt, sq) - potts_APLMLogz(i, pt, sq);
      logadd = e2_FLogsum(logadd, val);
    }
    tmp += exp(logadd);
    tmp /= (double)nsq;
    
    // add regularization term
    tmp += pt->mu * 2.0 * pt->h[i][a];
    
    dx[x++] = tmp;
  }
  
  for (j = 0; j < L; j++) 
    for (a = 0; a < Kg; a++) 
      for (b = 0; b < Kg; b++) {
	tmp    = 0.;
	logadd = -eslINFINITY;
	
	for (s = 0; s < nsq; s++) {
	  sq   = data->msa->ax[s];
	  resi = sq[i+1];
	  resj = sq[j+1];

	  if (resi == a && resj == b) tmp -= 1.;

	  if (resj == b) {
	    val    = 0.5 * (potts_APLMLognum(i, a, pt, sq) - potts_APLMLogz(i, pt, sq));
	    logadd = e2_FLogsum(logadd, val);
	  }	  
	}
	tmp += exp(logadd);
	tmp /= (double)nsq;
	
	// add regularization term
	tmp += pt->mu * 2.0 * pt->e[i][j][IDX(a,b,Kg)];
	
	dx[x++] = tmp;
  }
 
}


static double
logp_potts_ml(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_MLLogp(pt, msa, &logp, errbuf, verbose);
  if (status != eslOK) { printf("BAD: logp ML %f\n", logp); exit(1); }
  
#if 0
  printf("logp ML %f\n", logp);
#endif
  
  return logp;
}
 
static double
logp_potts_plm(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_PLMLogp(pt, msa, &logp, errbuf, verbose);
  if (status != eslOK) { printf("logp_potts_plm() failed: %s\n", errbuf); exit(1); }
  
#if 0
  printf("logp PLM %f\n", logp);
#endif
  
  return logp;
}
 
static double
logp_potts_aplm(int pos, PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_APLMLogp(pos, pt, msa, &logp, errbuf, verbose);
  if (status != eslOK) { printf("logp_potts_aplm() failed: pos %d. %s\n", pos, errbuf); exit(1); }
  
#if 1
  printf("pos %d logp APLM %f\n", pos, logp);
#endif
  
  return logp;
}


static int
symmetrize(PT *pt)
{
  double eij, eji;
  int    L  = pt->L;
  int    Kg = pt->abc->K+1;
  int    i, j;
  int    a, b;

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) 
      for (a = 0; a < Kg; a++)     
	for (b = 0; b < Kg; b++) {
	  eij = pt->e[i][j][IDX(a,b,Kg)];
	  eji = pt->e[j][i][IDX(b,a,Kg)];
	  pt->e[i][j][IDX(a,b,Kg)] = pt->e[j][i][IDX(b,a,Kg)] = 0.5 * (eij+eji); 
	}
  
  return eslOK;  
}

static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
		    const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls)
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f\n", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
  lbfgsfloatval_t  fx;
  int              i;
  int              status;

  fx = optimize_potts_func ((double *)x, n, instance);
  optimize_potts_dfunc_aplm((double *)x, n, instance, (double *)g);

  return fx;

 ERROR:
  return -eslINFINITY;
}
