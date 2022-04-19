/* r3d.c */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covgrammars.h"
#include "cacofold.h"
#include "correlators.h"
#include "contactmap.h"
#include "e2_profilesq.h"
#include "logsum.h"
#include "r3d.h"
#include "structure.h"


extern int
R3D_GetParam (R3Dparam  **ret_r3dp, char *errbuf, int verbose)
{
  R3Dparam  *r3dp;
  int        status;

   ESL_ALLOC(r3dp,      sizeof(R3Dparam));
   ESL_ALLOC(r3dp->HLp, sizeof(R3D_HLparam));
   ESL_ALLOC(r3dp->BLp, sizeof(R3D_BLparam));
   ESL_ALLOC(r3dp->ILp, sizeof(R3D_ILparam));

   *ret_r3dp = r3dp;
   return eslOK;

 ERROR:
   if (r3dp) {
     if (r3dp->HLp) free(r3dp->HLp);
     if (r3dp->BLp) free(r3dp->BLp);
     if (r3dp->ILp) free(r3dp->ILp);
     free(r3dp);
   }
  return status;
}
