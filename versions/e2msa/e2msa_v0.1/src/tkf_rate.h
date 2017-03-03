/* tkf_rate - 
 *
 */
#ifndef TKFRATE_INCLUDED
#define TKFRATE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"

#include "ratematrix.h"

/*****************************************************************
 * 1. TKF_RATE: the rate parameters.
 *****************************************************************/

typedef struct tkf_rate_s {
  double  mu;      /* deletion  rate of deletions   */
  double  lambda;  /* insertion rate of insertions  */

  double  etaz;                  /* instantaneous insertion at time zero eta_0  */      

  double  tsat;                  /* time of saturation */
  EMRATE *em;                    /* substitution rates */
} TKF_RATE;

extern TKF_RATE *tkf_rate_Create(const ESL_ALPHABET *abc);
extern TKF_RATE *tkf_rate_CreateWithValues(const ESL_ALPHABET *abc, double mu, double ld, double etaz, 
					   char *subsrate, ESL_DMATRIX *rate,  int subsratescale, double tol, char *errbuf, int verbose);
extern void      tkf_rate_Destroy(TKF_RATE *R);


#endif /*TKFRATE_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
