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
static int prna_add_counts(char *seq1, char *seq2, ELS_DMATRIX *prnaP);
static int urna_add_counts(char *seq1, char *seq2, ELS_DMATRIX *urnaP);
static int bg_add_counts  (char *seq1, char *seq2, double *bg, int dim);

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
Ribosum_matrix_Create(ESL_ALPHABET *abc)
{
  struct ribomatrix_s *ribosum = NULL;
  int    udim = abc->K;
  int    pdim = abc->K * abc->K;

  ESL_ALLOC(ribosum, sizeof(struct ribomatrix_s));
  ribosum->abc = abc;

  ribosum->prnaP = esl_dmatrix_Create(pdim,pdim);
  ribosum->prnaC = esl_dmatrix_Create(pdim,pdim);
  ribosum->prnaQ = esl_dmatrix_Create(pdim,pdim);
  ribosum->urnaP = esl_dmatrix_Create(pdim,pdim);
  ribosum->urnaC = esl_dmatrix_Create(pdim,pdim);
  ribosum->urnaQ = esl_dmatrix_Create(pdim,pdim);

  ESL_ALLOC(ribosum->bg, sizeof(double)*udim);

  /* set joints to zero, ready to add counts */
  esl_dmatrix_Set(ribosum->prnaP, 0.0);
  esl_dmatrix_Set(ribosum->urnaP, 0.0);
  esl_vec_DSet(ribosum->bg, udim, 0.0):
  
  esl_dmatrix_Set(ribosum->prnaC, -1.0);
  esl_dmatrix_Set(ribosum->urnaQ, -1.0);
  esl_dmatrix_Set(ribosum->prnaC, -1.0);
  esl_dmatrix_Set(ribosum->urnaQ, -1.0);

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
    free(ribosum);
  }
}


static int                  
ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float *threshid)
{
  double wgt;
  double pid;
  int    i, j;

  for (i = 0; i < msa->nseq; i++) {
    for (j = 0; j < i; j++) {
      
      esl_dst_CPairId(msa->aseq[i], msa->aseq[j], &pid, NULL, NULL);
      if (pid < thresid) continue; // not above cutoff

      wgt = msa->wgt[i] * msa->wgt[j]; 
      
      prna_add_counts(msa->aseq[i], msa->aseq[j], ribosum->prnaC);
      urna_add_counts(msa->aseq[i], msa->aseq[j], ribosum->urnaC);
      bg_add_counts  (msa->aseq[i], msa->aseq[j], ribosum->bg);

    }
  }

  return eslOK;
}


static int
prna_add_counts(char *seq1, char *seq2, ELS_DMATRIX *prnaP)
{
  return eslOK;
}

static int
urna_add_counts(char *seq1, char *seq2, ELS_DMATRIX *urnaP)
{
  return eslOK;
}
static int
bg_add_counts(char *seq1, char *seq2, double *bg, int dim)
{
  return eslOK;
}

