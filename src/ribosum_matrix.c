/* ribosum_matrix.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_msa.h"
#include "esl_msaweights.h"
#include "esl_vectorops.h"

#include "ribosum_matrix.h"


int
Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, char *errbuf)
{
  int status;
  
  /* calculate the weight BLOSUM-style */
  status = esl_msaweight_BLOSUM(msa, thresh1);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed msaweight_BLOSUM");
  
  /* add the counts */
  status = Ribosum_matrix_AddCounts(msa, ribosum, thresh2);

  return eslOK;
}


int 
Ribosum_matrix_AddCounts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float *thresh)
{
  int status;

  return eslOK;
}

struct ribomatrix_s *
Ribosum_matrix_Create()
{
  struct ribomatrix_s *ribosum = NULL;

  return ribosum;
}

void           
Ribosum_matrix_Destroy(struct ribomatrix_s *ribosum);
{
  if (ribosum) {
  }
}
