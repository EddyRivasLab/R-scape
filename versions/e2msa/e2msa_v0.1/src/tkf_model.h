/* tkf_model
 *
 */
#ifndef TKFMODEL_INCLUDED
#define TKFMODEL_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"

#include "e2.h"
#include "e2_config.h"
#include "tkf_rate.h"

/*****************************************************************
 * 2. TKF_MODEL: a model for one evolved sequence.
 *****************************************************************/
typedef struct tkf_model_s {
  double       betat;              /* lambda * beta(TKF)                            */
  double       gammat;             /* 1-exp(mu*t)                                   */
  ESL_DMATRIX *sub;                /* subsitution emissions.  sub[0..K-1][0..K-1]   */ 
  double       ins[e2_MAXABET];    /* insert emissions.       ins[0..K-1]           */
  
  double       time;
    
  const TKF_RATE      *R;        /* ptr to TKF_RATE used to build the model                    */
  const ESL_ALPHABET *abc;       /* ptr to alphabet info (tkf->abc->K is alphabet size)        */
} TKF_MODEL;


/* tkf_model.c */
extern TKF_MODEL *tkf_model_Create(TKF_RATE *R, double time, double alpha, const double *fins, int L, const ESL_ALPHABET *abc, double tol, int verbose);
extern void       tkf_model_Destroy(TKF_MODEL *evom);
extern int        tkf_model_BetaFunc(TKF_RATE *R, double time, double *ret_func);


#endif /*TKFMODEL_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
