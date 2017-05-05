/* pottsscore.c */

#include "rscape_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_minimize.h"
#include "esl_vectorops.h"

#include "pottsscore.h"
#include "pottsbuild.h"


int
potts_Score(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose)
{
  int status;

 if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "error scoring with potts");
  
  return eslOK;

 ERROR:
  return status;
}

