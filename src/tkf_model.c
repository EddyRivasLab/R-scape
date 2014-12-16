/* The evolutionary model.
 * 
 * Contents:
 *   1. The TKF_MODEL object: allocation, initialization, destruction.
 *   2. Convenience routines for setting fields in an TKF.
 *   3. Renormalization and rescaling counts in TKF.
 *   4. Debugging and development code.
 *   5. Other routines in the API.
 *   6. Unit tests.
 *   7. Test driver. 
 *   8. Copyright and license.
 * 
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_random.h"
#include "esl_dirichlet.h"
#include "esl_dmatrix.h"
#include "esl_stats.h"

#include "hmmer.h"

#include "e2.h"
#include "tkf_rate.h"
#include "tkf_model.h"
#include "ratematrix.h"

/*****************************************************************
 *# 1. The TKF_MODEL object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  tkf_model_Create()
 * Synopsis:  Allocate a new <TKF_MODEL>.
 *
 * Purpose:   Allocate a <TKF_MODEL> for symbol
 *            alphabet <abc>, and return a pointer to it.
 *            
 *            The TKF_MODEL only keeps a copy of the <abc> alphabet
 *            pointer. The caller is responsible for providing the
 *            alphabet, keeping it around while the  TKF_MODEL is in use,
 *            and (eventually) free'ing the alphabet when it's
 *            not needed any more. (Basically, just a step removed
 *            from keeping the alphabet as a global.)
 *
 * Throws:    <NULL> on allocation failure.
 */
TKF_MODEL *
tkf_model_Create(TKF_RATE *R, double time, double alpha, const double *fins, int L, const ESL_ALPHABET *abc, double tol, int verbose) 
{
  TKF_MODEL *tkf = NULL;
  double     tsubs;
  int        status;
 
  ESL_ALLOC(tkf, sizeof(TKF_MODEL));
  tkf->sub    = NULL;
  tkf->abc    = abc;
  tkf->R      = R;
  tkf->time   = time;

  tkf->gammat = 1.0 - exp(-R->mu*time);
  status = tkf_model_BetaFunc(R, time, &(tkf->betat)); if (tkf->betat < 0. || tkf->betat > 1.0 || isnan(tkf->betat)) { 
    printf("tkfbeta failed %f\n", tkf->betat); goto ERROR;
  }

  if (abc != NULL) {
    tsubs = time * alpha;
    tkf->sub = ratematrix_ConditionalsFromRate(tsubs, R->em->Qstar, tol, NULL, verbose);
  }
  if (fins != NULL) esl_vec_DCopy(fins, e2_MAXABET, tkf->ins);
  
 return tkf;
  
 ERROR:
  if (tkf != NULL) tkf_model_Destroy(tkf);
  return NULL;
}  


/* Function:  tkf_model_Destroy()
 * Synopsis:  Free a <TKF_MODEL>.
 *
 * Purpose:   Frees both the shell and body of an <tkf>.
 *            Works even if the <tkf> is damaged (incompletely allocated)
 *            or even <NULL>.
 *
 * Note:      Remember, leave reference pointers like abc, gm, and
 *            bg alone. These are under the application's control not ours.
 *
 * Returns:   (void).
 */
void
tkf_model_Destroy(TKF_MODEL *tkf)
{
  if (tkf == NULL) return;
  if (tkf->sub)  esl_dmatrix_Destroy(tkf->sub);

  free(tkf);
  return;
}

/* 
 */
int
tkf_model_BetaFunc(TKF_RATE *R, double time, double *ret_beta)
{
  double a;                /* a = ldI-muI                                    */
  double E;                /* E = exp(a*t)                                   */
  double A;                /* A = lambda * exp(a*t) = exp(a*t + log(lambda)) */
  double num, den;
  double beta;
  int    status;

  if (R->lambda >= R->mu) { printf("tkf condition: ld < mu is violated \n"); status = eslFAIL; goto ERROR;} 

  /* some asignments */
  a = R->lambda - R->mu;
  E = exp(time * a);
  A = exp(time * a + log(R->lambda));

  num = R->lambda * (1.0 - E);
  den = R->mu - A;
  
  beta = (fabs(den) > 0.)? num/den : 1e-8;

#if 0
  printf("\neta[t=%.3f] = %.8f || ld=%.8f\tmu=%.8f\tetaz=%.8f num %f den %f\n", time, beta, R->lambda, R->mu, R->etaz, num, den);
#endif
if (beta < 0.) { printf("tkfbeta is negative \n"); status = eslFAIL; goto ERROR; }
  *ret_beta = beta;

  return eslOK;

 ERROR:
  return status;
}

