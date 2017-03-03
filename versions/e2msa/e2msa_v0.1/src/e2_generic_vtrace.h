/* viterbi traceback
 *
 *   
*/
#ifndef E2_GENERIC_VTRACE_INCLUDED
#define E2_GENERIC_VTRACE_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"

#include "e2_gmx.h"
#include "e2_profile.h"
#include "e2_profilesq.h"
#include "e2_trace.h"

#include "hmmer.h"

/* e2_generic_vtrace.c */
extern int e2_GTrace(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, E2_GMX *gx, E2_TRACE *tr);
#endif
