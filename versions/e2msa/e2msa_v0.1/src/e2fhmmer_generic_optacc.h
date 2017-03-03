/* e2fhmmer_generic_optacc
 *
 *   
*/
#ifndef E2FHMMER_GENERIC_OPTACC_INCLUDED
#define E2FHMMER_GENERIC_OPTACC_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_random.h"

#include "e1_bg.h"
#include "e1_model.h"
#include "e2_gmx.h"
#include "e2_profile.h"
#include "e2_profilesq.h"
#include "e2hmmer_profile.h"

#include "hmmer.h"
#include "evohmmer.h"

/* e2fhmmer_generic_optacc.c */
extern int e2fhmmer_GOptimalAccuracy(const E2HMMER_PROFILE *gm, const E2_GMX *pp, E2_GMX *gx, 
				     const PSQ *psql, const PSQ *psqr, float *ret_e, int *ret_k);
extern int e2fhmmer_GOATrace        (ESL_RANDOMNESS *r, P7_RATE *R7, float timel, float timer, 
				     P7_BG *bg7, const E2HMMER_PROFILE *gm, 
				     int kmax, const E2_GMX *pp, const E2_GMX *gx, 
				     E2_TRACE *tr, PSQ *psq1, PSQ *psq2, PSQ **ret_psqA);

#endif
