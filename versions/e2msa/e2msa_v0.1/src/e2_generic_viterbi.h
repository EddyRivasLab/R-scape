/* viterbi
 *
 *   
*/
#ifndef E2_GENERIC_VITERBI_INCLUDED
#define E2_GENERIC_VITERBI_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"

#include "e2_gmx.h"
#include "e2_profile.h"
#include "e2_profilesq.h"

#include "hmmer.h"

/* e2_generic_viterbi.c */
extern int e2_GViterbi (const PSQ *psq1, const PSQ *psq2, const E2_PROFILE *gm, E2_GMX *gx, float *opt_sc);


#endif
