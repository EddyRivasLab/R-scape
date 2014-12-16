/* e2fhmmer_deneric_decoding
 *
 *   
*/
#ifndef E2FHMMER_GENERIC_DECODING_INCLUDED
#define E2FHMMER_GENERIC_DECODING_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"

#include "e2_profile.h"
#include "e2_gmx.h"
#include "e2hmmer_profile.h"

#include "hmmer.h"

/* e2fhmmer_generic_decoding.c */
extern int e2fhmmer_GDecoding (const E2HMMER_PROFILE *gm, const E2_GMX *fwd, E2_GMX *bck, float overal_sc, E2_GMX *pp);

#endif
