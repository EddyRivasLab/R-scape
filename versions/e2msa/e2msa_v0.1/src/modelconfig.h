/* modelconfig
 *
 *   
*/
#ifndef MODELCONFIG_INCLUDED
#define MODELCONFIG_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"

#include "hmmer.h"

#include "e1_model.h"
#include "e2_gmx.h"
#include "e2_profile.h"
#include "e2hmmer_profile.h"
#include "e2_profilesq.h"

#include "evohmmer.h"

extern int e2_ProfileConfig     (const E1_MODEL *evoml, const E1_MODEL *evomr, float *frq, E2_PROFILE *gm, 
				 float L, E2_ALI e2ali, int verbose);
extern int e2hmmer_ProfileConfig(const P7_RATE *R, float tl, float tr, const P7_HMM *evo7l, const P7_HMM *evo7r, P7_BG *bg7, 
				 E2HMMER_PROFILE *gm7, float L, E2_ALI e2ali, int mode, int verbose);
extern int select_discrete_time(const P7_RATE *R, float time);
#endif
