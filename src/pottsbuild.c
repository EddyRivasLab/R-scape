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
#include "esl_random.h"
#include "esl_vectorops.h"

#include "minimize.h"
#include "pottsbuild.h"
#include "pottsscore.h"
#include "correlators.h"

static int    optimize_all_pack_paramvector   (double *p, long np, struct optimize_data *data);
static int    optimize_all_unpack_paramvector (double *p, long np, struct optimize_data *data);
static int    optimize_aplm_pack_paramvector   (double *p, long np, struct optimize_data *data);
static int    optimize_aplm_unpack_paramvector (double *p, long np, struct optimize_data *data);
static void   optimize_bracket_define_direction(double *p, long np, struct optimize_data *data);
static double optimize_potts_func              (double *p, long np, void *dptr);
static double func_potts_ml                             (PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
static double func_potts_plm                            (PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
static double func_potts_aplm                  (int pos, PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
static int    symmetrize(PT *pt);

PT *
potts_Build(ESL_RANDOMNESS *r, ESL_MSA *msa, double ptmu, PTTRAIN pttrain, PTSCTYPE ptsctype, FILE *pottsfp, float tol, char *errbuf, int verbose)
{
  PT              *pt = NULL;
  float            firststep;
  int              status;
  
  pt = potts_Create(msa->alen, msa->abc->K, msa->abc, ptmu, pttrain, ptsctype);
  if (pt == NULL) ESL_XFAIL(eslFAIL, errbuf, "error creating potts");
  
  // Initialize 
  //status = potts_AssignGaussian(r, pt, 0., 1.0);
  status = potts_AssignGT(r, msa, pt, tol, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  
  /* init */
  firststep = 1.0;
  switch (pt->train) {
  case NONE:
   ESL_XFAIL(eslFAIL, errbuf, "errork, you should not be here");
     break;
  case ML:
  case PLM:
    status = potts_OptimizeGDALL(pt, msa, firststep, tol, errbuf, verbose);
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "error all optimizing potts");
    break;
  case APLM:
    status = potts_OptimizeGDAPLM(pt, msa, firststep, tol, errbuf, verbose);
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
potts_OptimizeGDALL(PT *pt, ESL_MSA *msa, float firststep, float tol, char *errbuf, int verbose)
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
  optimize_all_pack_paramvector(p, (long)np, &data);
 
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
  optimize_all_unpack_paramvector(p, (long)np, &data);
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
potts_OptimizeGDAPLM(PT *pt, ESL_MSA *msa, float firststep, float tol, char *errbuf, int verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 logp;
  int                    L = msa->alen;
  int                    K = msa->abc->K;
  int                    np;
  int                    i;
  int                    status;

  np = K*(1+L*K);     /* the variables hi eij */
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
    data.firststep  = firststep;
    data.tol        = tol;
    data.errbuf     = errbuf;
    data.verbose    = verbose;
    
    /* Create the parameter vector.
     */
    optimize_aplm_pack_paramvector(p, (long)np, &data);
    
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
    optimize_aplm_unpack_paramvector(p, (long)np, &data);
    if (verbose) printf("END POTTS APLM OPTIMIZATION for position %d\n", i);
  }
  if (1||verbose) printf("END POTTS APLM OPTIMIZATION\n");

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


PT *
potts_Create(int64_t L, int K, ESL_ALPHABET *abc, double mu, PTTRAIN pttrain, PTSCTYPE  ptsctype)
{
  PT  *pt = NULL;
  int  K2 = K * K;
  int  i, j;
  int  status;
  
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
  int L = pt->L;
  int K = pt->abc->K;
  int i, j;
  int a, b;
  
  for (i = 0; i < L; i++) {
    for (a = 0; a < K; a ++)
      pt->h[i][a] = esl_rnd_Gaussian(r, mu, sigma);
    
    for (j = i; j < L; j++)
      for (a = 0; a < K; a ++) 
	for (b = 0; b < K; b ++) {
	  pt->e[i][j][IDX(a,b,K)] = (i==j)? 0 : esl_rnd_Gaussian(r, mu, sigma);
	  pt->e[j][i][IDX(b,a,K)] = pt->e[i][j][IDX(a,b,K)];
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
  int              L = pt->L;
  int              K = pt->abc->K;
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
	  obs = mi->nseff[i][j] * mi->pp[i][j][IDX(a,b,K)];
	  gt  = (exp > 0. && obs > 0.) ? obs * log (obs / exp) : pseudoc;
	  
	  pt->e[i][j][IDX(a,b,K)] = (i==j)? 0 : gt;
	  pt->e[j][i][IDX(b,a,K)] = pt->e[i][j][IDX(a,b,K)];
      }
  }
  
  for (i = 0; i < L; i++) 
    for (a = 0; a < K; a ++) {
      hi = 0;
      
      for (j = 0; j < L; j++)
	for (b = 0; b < K; b ++) 
	  hi += pt->e[i][j][IDX(a,b,K)];
      pt->h[i][a] = hi/K;
    }

  //potts_Write(stdout, pt);

  corr_Destroy(mi);
  return eslOK;

 ERROR:
  if (mi) corr_Destroy(mi);
  return status;
}

// sum_b eij(a,b) = sum_a eij(a,b) = sum_a hi[a] = 0
//
// eij(a,b) → eij(a,b) − 1/K * \sum_c eij(c,b) − 1/k * \sum_d eij(a,d) + \1/K^2 * \sum_c \sum_d eij(c,d)
// hi(a)    → hi(a)    - \sum_{j\neq i} 1/K * \sum_d eij(a,d)
int
potts_GaugeZeroSum(PT *pt, char *errbuf, int verbose)
{
  int      K = pt->abc->K;
  int      K2 = K*K; 
  double   sumi[K];
  double   sumj[K];
  double   sumh;
  double   sum;
  double   tol = 1e-5;
  int      i, j;
  int      a, b;
  int      status;

  for (i = 0; i < pt->L; i++) 
    for (a = 0; a < K; a ++) {
      
      sumh = 0.0;  
      for (j = 0; j < pt->L; j++) {
	if (i==j) continue;
	for (b = 0; b < K; b ++) 
	  sumh += pt->e[i][j][IDX(a,b,K)];
      }
      
      pt->h[i][a] -= sumh/K;    
    }
  
  for (i = 0; i < pt->L; i++) 
    for (j = i+1; j < pt->L; j++) {    
      sum = 0;     
      for (a = 0; a < K; a ++) {
	sumi[a] = 0;
	sumj[a] = 0; 
	for (b = 0; b < K; b ++) {
	  sum        += pt->e[i][j][IDX(a,b,K)];
	  sumi[a]    += pt->e[i][j][IDX(b,a,K)];
	  sumj[a]    += pt->e[i][j][IDX(a,b,K)];
	}
      }
      
      for (a = 0; a < K; a ++) 
	for (b = 0; b < K; b ++) {
	  pt->e[i][j][IDX(a,b,K)] -= sumi[b]/K + sumj[a]/K - sum/K2;
	  pt->e[j][i][IDX(a,b,K)]  = pt->e[i][j][IDX(a,b,K)];
	}
    }
  
#if 0
  // check it is a zero sum
  for (i = 0; i < pt->L; i++) {
    
    sum = 0;
    for (a = 0; a < K; a ++) sum += pt->h[i][a];
    if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for h failed. sum = %f", sum);
    
    for (j = i+1; j < pt->L; j++) {
      
      for (a = 0; a < K; a ++) {
	sum = 0;	
	for (b = 0; b < K; b ++) 
	  sum += pt->e[i][j][IDX(a,b,K)];
	if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for eij failed sumj %f | %d %d %d", sum, i, j, a);
      }
      
      for (a = 0; a < K; a ++) {
	sum = 0;
	for (b = 0; b < K; b ++) 
	  sum += pt->e[i][j][IDX(b,a,K)];
	if (fabs(sum) > tol) ESL_XFAIL(eslFAIL, errbuf, "zerosum for eij failed sumi %f | %d %d %d", sum, i, j, a);
      }
      
    } // for all i,j   
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
  double           hi, eij;
  int              ln = 0;
  int              K  = 0;
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
      if (ln == 0) K = t;
      
      if (t  == K) ln ++;
      else { L = ln; break; }
    }
  esl_fileparser_Close(efp);
  printf("L %d K %d\n", L, K);

  pt = potts_Create(L, K, NULL, 0.0, NONE, SCNONE);
  if (abc) {
    if (K != abc->K) ESL_XFAIL(eslFAIL, errbuf, "wrong alphabet for file %s", paramfile);
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
  int K = pt->abc->K;
  int i, j;
  int a, b;
  
  for (i = 0; i < L; i++) {
    fprintf(fp, "%d ", i);
    for (a = 0; a < K; a ++) fprintf(fp, "%f ", pt->h[i][a]);
    fprintf(fp, "\n");
  }
  
  for (i = 0; i < L; i++) 
    for (j = i+1; j < L; j++) {
      fprintf(fp, "%d %d ", i, j);
      for (a = 0; a < K; a ++) 
	for (b = 0; b < K; b ++) {
	  fprintf(fp, "%f ", pt->e[i][j][IDX(a,b,K)]);
	}
      fprintf(fp, "\n");
    }
}


/*------------------------------- internal functions ----------------------------------*/

static int
optimize_all_pack_paramvector(double *p, long np, struct optimize_data *data)
{
  int   L = data->msa->alen;
  int   K = data->msa->abc->K;
  int   x = 0;
  int   i, j;
  int   a;

  for (i = 0; i < L; i++) {
    for (a = 0; a < K; a++)       p[x++] = data->pt->h[i][a];
    for (j = i+1; j < L; j++) 
      for (a = 0; a < K*K; a++) { p[x++] = data->pt->e[i][j][a]; p[x++] = data->pt->e[j][i][a]; }
  }
  return eslOK;  
}

static int
optimize_aplm_pack_paramvector(double *p, long np, struct optimize_data *data)
{
  int   L = data->msa->alen;
  int   K = data->msa->abc->K;
  int   x = 0;
  int   i = data->pos;
  int   j;
  int   a;

  for (a = 0; a < K; a++)     p[x++] = data->pt->h[i][a];
  for (j = 0; j < L; j++) {
    if (j==i) continue;
    for (a = 0; a < K*K; a++) p[x++] = data->pt->e[i][j][a];
  }
  return eslOK;  
}

static int
optimize_all_unpack_paramvector(double *p, long np, struct optimize_data *data)
{
  double eij, eji;
  int    L = data->msa->alen;
  int    K = data->msa->abc->K;
  int    x = 0;
  int    i, j;
  int    a;
  
  for (i = 0; i < L; i++) {
    for (a = 0; a < K; a++)    data->pt->h[i][a]    = p[x++];
    for (j = i+1; j < L; j++) 
      for (a = 0; a < K*K; a++) {
	eij = p[x++];
	eji = p[x++];
	data->pt->e[i][j][a] = data->pt->e[j][i][a] = 0.5*(eij+eji);
      }
  }

 return eslOK;
}

static int
optimize_aplm_unpack_paramvector(double *p, long np, struct optimize_data *data)
{
  int    L = data->msa->alen;
  int    K = data->msa->abc->K;
  int    x = 0;
  int    i = data->pos;
  int    j;
  int    a;
  
  for (a = 0; a < K; a++)     data->pt->h[i][a]    = p[x++];
    for (j = 0; j < L; j++) {
    if (j==i) continue;
    for (a = 0; a < K*K; a++) data->pt->e[i][j][a] = p[x++];
  }
  
  return eslOK;
}


static void
optimize_bracket_define_direction(double *u, long np, struct optimize_data *data)
{
  int x;
  for (x = 0; x < np; x++) u[x] = 0.01;
  u[np] = 0.5;
}

static double
optimize_potts_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  
  
  switch(data->pt->train) {
  case ML:
    optimize_all_unpack_paramvector(p, np, data);
    data->logp = func_potts_ml(data->pt, data->msa, data->tol, data->errbuf, data->verbose);
    break;
  case PLM:
    optimize_all_unpack_paramvector(p, np, data);
    data->logp = func_potts_plm(data->pt, data->msa, data->tol, data->errbuf, data->verbose);
    break;
  case APLM:
    optimize_aplm_unpack_paramvector(p, np, data);
    data->logp = func_potts_aplm(data->pos, data->pt, data->msa, data->tol, data->errbuf, data->verbose);
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
func_potts_ml(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_MLLogp(pt, msa, &logp, errbuf, verbose);
  if (status != eslOK) exit(1);
  
#if 1
  printf("logp ML %f\n", logp);
#endif
  
  return -logp;
}
 
static double
func_potts_plm(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_PLMLogp(pt, msa, &logp, errbuf, verbose);
  if (status != eslOK) exit(1);
  
#if 1
  printf("logp PLM %f\n", logp);
#endif
  
  return -logp;
}
 
static double
func_potts_aplm(int pos, PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  double  logp;
  int     status;
  
  status = potts_APLMLogp(pos, pt, msa, &logp, errbuf, verbose);

  if (status != eslOK) exit(1);
  
#if 0
  printf("pos %d logp APLM %f\n", pos, logp);
#endif
  
  return -logp;
}

static int
symmetrize(PT *pt)
{
  double eij, eji;
  int    L = pt->L;
  int    K = pt->abc->K;
  int    x = 0;
  int    i, j;
  int    a, b;

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) {
      if (i==j) continue;
      
      for (a = 0; a < K; a++)     
	for (b = 0; b < K; b++) {
	  eij = pt->e[i][j][IDX(a,b,K)];
	  eji = pt->e[j][i][IDX(b,a,K)];
	  pt->e[i][j][IDX(a,b,K)] = pt->e[j][i][IDX(b,a,K)] = 0.5 * (eij+eji); 
	}
    }
  
  return eslOK;  
}
