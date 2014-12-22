/* ribosum_matrix.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_msa.h"
#include "esl_msaweights.h"
#include "esl_vectorops.h"

#include "ribosum_matrix.h"

static int ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float *thresh2);

int
Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, char *errbuf)
{
  int status;
  
  /* calculate the weight BLOSUM-style */
  status = esl_msaweight_BLOSUM(msa, thresh1);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed msaweight_BLOSUM");
  
  /* add the counts */
  status = ribosum_matrix_add_counts(msa, ribosum, thresh2);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed ribosum_matrix_add_counts()");

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
    if (ribosum->prnaP) esl_dmatrix_Destroy(ribosum->prnaP);
    if (ribosum->prnaC) esl_dmatrix_Destroy(ribosum->prnaC);
    if (ribosum->prnaQ) esl_dmatrix_Destroy(ribosum->prnaQ);
    if (ribosum->urnaP) esl_dmatrix_Destroy(ribosum->urnaP);
    if (ribosum->urnaC) esl_dmatrix_Destroy(ribosum->urnaC);
    if (ribosum->urnaQ) esl_dmatrix_Destroy(ribosum->urnaQ);
    if (ribosum->bg)    free(ribosum->bg);
    if (ribosum->name)  free(ribosum->name);
  }
}


static int                  
ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float *thresh2)
{

  return eslOK;
}

