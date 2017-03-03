/* The rates for the evolutionary model
 * 
 * Contents:
 *   1. The E1_RATE object: allocation, initialization, destruction.
 *   2. Convenience routines for setting fields in an E1.
 *   3. Renormalization and rescaling counts in E1.
 *   4. Debugging and development code.
 *   5. Other routines in the API.
 *   6. Unit tests.
 *   7. Test driver. 
 *   8. Copyright and license.
 * 
 */
#include "p7_config.h"

#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_random.h"
#include "esl_rootfinder.h"
#include "esl_dirichlet.h"
#include "esl_dmatrix.h"

#include "hmmer.h"

#include "e2.h"
#include "tkf_rate.h"
#include "ratematrix.h"


/*****************************************************************
 *# 1. The TKF_MODEL object: allocation, initialization, destruction.
 *****************************************************************/
TKF_RATE *
tkf_rate_Create(const ESL_ALPHABET *abc)
{
  TKF_RATE *R = NULL;
  int      status;

   ESL_ALLOC(R, sizeof(TKF_RATE));

  /* initialize all transition rates to zero */
  R->mu = 0.0;
  R->lambda = 0.0;
  R->etaz = 0.0;
 
  /* add the substitution rates */
 if (abc != NULL) {
     R->em = ratematrix_emrate_Create(abc, 1);
  }

  return R;

 ERROR:
  return NULL;
}

TKF_RATE *
tkf_rate_CreateWithValues(const ESL_ALPHABET *abc, double mu, double ld, double etaz,  
			  char *subsrate, ESL_DMATRIX *rate, int subsratescale, double tol, char *errbuf, int verbose)
{
  TKF_RATE *R = NULL;
  int      status;

  if ((R = tkf_rate_Create(abc)) == NULL) ESL_XFAIL(eslFAIL, errbuf, "tkf_rate_Create failed");
  
  R->mu     = mu;
  R->lambda = ld;
  R->etaz   = etaz;

  R->em = NULL;
  if (ratematrix_emrate_LoadRate(R->em, subsrate, rate, NULL, subsratescale, tol, errbuf, verbose) != eslOK) goto ERROR;
   
  return R;

 ERROR:
  return NULL;
}



void
tkf_rate_Destroy(TKF_RATE *R)
{
  if (R == NULL) return;
  if (R->em != NULL) ratematrix_emrate_Destroy(R->em, 1);
  free(R);
  return;
}


