/* e2fhmmer_generic_fwdback
 *
 *   
*/
#ifndef E2FHMMER_GENERIC_FWDBACK_INCLUDED
#define E2FHMMER_GENERIC_FWDBACK_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"

#include "e2_gmx.h"
#include "e2_profile.h"
#include "e2_profilesq.h"
#include "e2hmmer_profile.h"
#include "hmmer.h"

/* e2fhmmer_generic_fwdback.c */
extern int e2fhmmer_GForward (const PSQ *psq1, const PSQ *psq2, const E2HMMER_PROFILE *gm, E2_GMX *gx, float *opt_sc);
extern int e2fhmmer_GBackward(const PSQ *psq1, const PSQ *psq2, const E2HMMER_PROFILE *gm, E2_GMX *gx, float *opt_sc);
extern int isgap(float *p, int K);
#endif
