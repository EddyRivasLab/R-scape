/* e2_gmx
 *
 *   
*/
#ifndef E2F_GENERIC_DECODING_INCLUDED
#define E2F_GENERIC_DECODING_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"

#include "e2_gmx.h"
#include "e2_profile.h"

#include "hmmer.h"

/* e2f_generic_decoding.c */
extern int e2f_GDecoding (const E2_PROFILE *gm, const E2_GMX *fwd, E2_GMX *bck, E2_GMX *pp, char *errbuf);


#endif
