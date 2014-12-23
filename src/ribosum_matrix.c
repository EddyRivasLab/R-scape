/* ribosum_matrix.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_msa.h"
#include "esl_msaweights.h"
#include "esl_vectorops.h"

#include "e2_gmx.h"
#include "ribosum_matrix.h"

static int ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float *thresh2, char errbuf);
static int prna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ELS_DMATRIX *prnaP);
static int urna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ELS_DMATRIX *urnaP);
static int bg_add_counts  (ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, int alen, double *bg, int dim);

int
Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, char *errbuf)
{
  int status;
  
  /* calculate the weight BLOSUM-style */
  status = esl_msaweight_BLOSUM(msa, thresh1);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed msaweight_BLOSUM");
  
  /* add the counts */
  status = ribosum_matrix_add_counts(msa, ribosum, thresh2, errbuf);
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
ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float *threshid, char errbuf)
{
  double   wgt;
  double   pid;
  int    **ct = NULL;   // ct array for each seq
  int      i, j;
  int      status;
  
  ESL_ALLOC(ct, sizeof(int *)*msa->nseq);
  for (i = 0; i < nsma->nseq; i ++) {
    ESL_ALLOC(ct[i], sizeof(int)*msa->alen);
    if      (msa->ss_cons) esl_wuss2ct(msa->ss_cons, msa->alen, ct[i]);
    else if (msa->ss[i])   esl_wuss2ct(msa->ss[i],   msa->alen, ct[i]);
    else                   ESL_XFAIL(eslFAIL, "no ss for sequence %d\n", errbuf);
  }
  

  for (i = 0; i < msa->nseq; i++) {
    for (j = 0; j < i; j++) {
      esl_dst_CPairId(msa->aseq[i], msa->aseq[j], &pid, NULL, NULL);
      if (pid < thresid) continue; // not above cutoff

      wgt = msa->wgt[i] * msa->wgt[j]; 
      
      prna_add_counts(ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], msa->ss[i], msa->ss[j], ct[i], ct[j], ribosum->prnaC);
      urna_add_counts(ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], msa->ss[i], msa->ss[j], ct[i], ct[j], ribosum->urnaC);
      bg_add_counts  (ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], ribosum->bg);
    }
  }

  for (i = 0; i < nsma->nseq; i ++) free(ct[i]);
  free(ct);

  return eslOK;
}


static int
prna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ELS_DMATRIX *prnaP)
{
  int  col;
  int  cpos1, cpos2;
  char c1, c2;
  char cc1, cc2;
  int  left, right;

  for (col = 0; col < alen; col++) {
    c1 = seq1[col];
    c2 = seq2[col];
    
    if (esl_abc_CIsResidue(abc, c1) && esl_abc_CIsResidue(abc, c2)) 
      {
	cpos1 = ct1[pos+1];
	cpos2 = ct2[pos+1];
	if (cpos1 != cpos2 || cpos1 < col) continue;

	cc1 = ss1[cpos1-1];
	cc2 = ss2[cpos2-1];
	if (esl_abc_CIsResidue(abc, cc1) && esl_abc_CIsResidue(abc, cc2)) 
	  {
	    left  = ID(esl_abc_DigitizeSymbol(abc, c1),  esl_abc_DigitizeSymbol(abc, c2),  abc->K);
	    right = ID(esl_abc_DigitizeSymbol(abc, cc1), esl_abc_DigitizeSymbol(abc, cc2), abc->K);

	    prnaP->mx[left][right] += wgt;
	    prnaP->mx[right][left] += wgt;
	  }
      }
  }

  return eslOK;
}

static int
urna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ELS_DMATRIX *urnaP)
{
  int  col;
  char c1, c2;
  
  for (col = 0; col < alen; col++) {
    c1 = seq1[col];
    c2 = seq2[col];
    if (ct1[pos+1] != 0 || ct2[pos+1] != 0) continue;
    if (esl_abc_CIsResidue(abc, c1) && esl_abc_CIsResidue(abc, c2)) {
      urnaP->mx[esl_abc_DigitizeSymbol(abc, c1)[esl_abc_DigitizeSymbol(abc, c2)] += wgt;
      urnaP->mx[esl_abc_DigitizeSymbol(abc, c2)[esl_abc_DigitizeSymbol(abc, c2)] += wgt;
    }
  }

  return eslOK;
}

static int
bg_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, int alen, double *bg, int dim)
{
  int col;
  int c1, c2;
  
  for (col = 0; col < alen; col++) {
    
    c1 = seq1[col];
    c2 = seq2[col];
    
    if (esl_abc_CIsResidue(abc, c1) && esl_abc_CIsResidue(abc, c2)) 
      {
	bg[esl_abc_DigitizeSymbol(abc, c1)] += wgt;
	bg[esl_abc_DigitizeSymbol(abc, c2)] += wgt;
      }
  }

  return eslOK;
}

