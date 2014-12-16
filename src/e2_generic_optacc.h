/* e2_generic_optacc
 *
 *   
*/
#ifndef E2_GENERIC_OPTACC_INCLUDED
#define E2_GENERIC_OPTACC_INCLUDED

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

#include "hmmer.h"

/* e2f_generic_optacc.c */
extern int e2_GOptimalAccuracy(const E2_PROFILE *gm, const E2_GMX *pp,       E2_GMX *gx, float *ret_e);
extern int e2_GOATrace        (ESL_RANDOMNESS *r, const E1_MODEL *evo1, const E1_MODEL *evo2, const float *frq, const E2_PROFILE *gm, 
			       const E2_GMX *pp, const E2_GMX *gx, E2_TRACE *tr, 
			       const PSQ *psq1, const PSQ *psq2, PSQ **ret_psqA);


#endif
