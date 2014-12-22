/* ribosum_matrix.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_msa.h"
#include "esl_vectorops.h"


int
RibosumMatrix(ESL_MSA *msa, float thresh1, float thresh2, char *errbuf)
{
  int status;

  /* calculate the weight BLOSUM-style */
  status = esl_msaweight_BLOSUM(msa, thresh1);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed msaweight_BLOSUM");
  
  return eslOK;
}



