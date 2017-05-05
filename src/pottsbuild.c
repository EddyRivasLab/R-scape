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

static int    optimize_pack_paramvector        (double *p, long np, struct optimize_data *data);
static int    optimize_unpack_paramvector      (double *p, long np, struct optimize_data *data);
static void   optimize_bracket_define_direction(double *p, long np, struct optimize_data *data);
static double optimize_potts_func              (double *p, long np, void *dptr);
static double func_potts                       (PT *pt, ESL_MSA *msa, float  time,     float tol, char *errbuf, int verbose);
static int    optimize_potts                   (PT *pt, ESL_MSA *msa, float *ret_time, float tol, char *errbuf, int verbose);


int
potts_Build(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  int status;

  status = optimize_potts(pt, msa, tol, errbuf, verbose);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "error optimizing potts");
  
  return eslOK;

 ERROR:
  return status;
}

PT *
potts_Create(int64_t L, ESL_ALPHABET *abc)
{
  PT  *pt = NULL;
  int  K  = abc->K;
  int  K2 = K * K;
  int  i, j;
  int  status;
  
  ESL_ALLOC(pt, sizeof(PT));
  pt->L   = L;
  pt->abc = abc;

  ESL_ALLOC(pt->e,                   sizeof(double **) * L);
  ESL_ALLOC(pt->h,                   sizeof(double  *) * L);
  for (i = 0; i < L; i++) {
    ESL_ALLOC(pt->e[i],              sizeof(double  *) * L);
    for (j = 0; j < L; j++) {
       ESL_ALLOC(pt->e[i][j],        sizeof(double  ) * K2);
    }
  }
   
  /* initialize for adding counts */
  for (i = 0; i < L; i++) {
    esl_vec_DSet(pt->h[i], K, 0.0); 
 
    for (j = 0; j < L; j++) {
      esl_vec_DSet(pt->e[i][j], K2, 0.0); 
    }
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
      for (j = 0; j < pt->L; j++) {
	free(pt->e[i][j]);
      }
     free(pt->h[i]);
    }
    free(pt->e);
    free(pt->h);
    free(pt);
  }
}

/*------------------------------- internal functions ----------------------------------*/

static int
optimize_pack_paramvector(double *p, long np, struct optimize_data *data)
{
  int   x = 0;
  
  p[x] = (data->time < 1.0)? log(data->time) : data->time - 1.0;

  return eslOK;  
}


static int
optimize_unpack_paramvector(double *p, long np, struct optimize_data *data)
{
  float time;
  float tmax = 10.0;
  int   x = 0;
  
  time = (p[x] < 0.0)? exp(p[x]) : p[x] + 1.0; 
  if (time > tmax) time = tmax;
  
  data->time = time;
  return eslOK;
}

static void
optimize_bracket_define_direction(double *u, long np, struct optimize_data *data)
{
  int x;
  for (x = 0; x < np; x++) u[x] = 0.25;
  u[np] = 0.25;
}


static int
optimize_potts(PT *pt, ESL_MSA *msa, double *ret_time, float tol, char *errbuf, int verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 sc;
  double                 sc_init;
  float                  time;
  int                    np;
  int                    status;
  
  time = *ret_time;
  sc_init = func_potts(fx, time, tol, errbuf, be_verbose);
  if (sc_init == eslINFINITY) {
    *ret_sc  = sc_init;
    *ret_time = time;
    return eslOK;
  }

  np = 1;     /* variable: time */
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);
  
 /* Copy shared info into the "data" structure
   */
  data.time       = time;
  data.firststep  = firststep;
  data.fx         = fx;
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
		       tol, wrk, &sc);
  if (status != eslOK) 
    esl_fatal("optimize_potts(): bad bracket minimization");	
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  data.sc = -sc;
  if (be_verbose) printf("END MSV OPTIMIZATION: time %f usc %f --> %f\n", data.time, sc_init, data.sc);
  
  *ret_sc   = data.sc;
  *ret_time = data.time;
  
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

static double
optimize_potts_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, np, data);
  
  data->sc = func_potts(ddata->tol, data->errbuf, data->be_verbose);
  
  if (data->usc == eslINFINITY) data->sc = 1000.;
  return -(double)data->sc;
}

static double
func_potts(PT *pt, ESL_MSA *msa, float time, float tol, char *errbuf, int verbose)
{
  float  sc;
  
  potts_Score(pt, msa, fx, &(sc));
  
#if 0
  printf("time %f sc %f\n", time, sc);
#endif

  return (double)sc;
 }
