/* ribosum_matrix.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_distance.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "ribosum_matrix.h"

static int ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh, char *errbuf);
static int prna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *prnaP);
static int urna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *urnaP);
static int bg_add_counts  (ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, int alen, double *bg);

int
Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, FILE *fp, int verbose, char *errbuf)
{						
  float sum;
  int   status;
  
  /* calculate the weight BLOSUM-style */
  status = esl_msaweight_BLOSUM(msa, thresh1);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed msaweight_BLOSUM");
  
  /* add the counts */
  status = ribosum_matrix_add_counts(msa, ribosum, thresh2, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed ribosum_matrix_add_counts()");

  /* normalize */
  sum = esl_dmx_Sum(ribosum->prnaP);
  if (sum > 0.) esl_dmx_Scale(ribosum->prnaP, 1.0/sum);
  sum = esl_dmx_Sum(ribosum->urnaP);
  if (sum > 0.) esl_dmx_Scale(ribosum->urnaP, 1.0/sum);
  sum = esl_vec_DSum(ribosum->bg, ribosum->abc->K);
  if (sum > 0.) esl_vec_DScale(ribosum->bg, ribosum->abc->K, 1.0/sum);

  /* calculate conditionals and marginals */

  /* calculate rates */

  if (verbose) Ribosum_matrix_Write(stdout, ribosum); 
  Ribosum_matrix_Write(fp, ribosum);
  return eslOK;

 ERROR:
  return status;
}

struct ribomatrix_s *
Ribosum_matrix_Create(ESL_ALPHABET *abc, char *name)
{
  struct ribomatrix_s *ribosum = NULL;
  int    udim = abc->K;
  int    pdim = abc->K * abc->K;
  int    status;
  
  ESL_ALLOC(ribosum, sizeof(struct ribomatrix_s));
  ribosum->abc = abc;
  ribosum->name = NULL;
  if (name) ribosum->name = name;

  ribosum->prnaP = esl_dmatrix_Create(pdim, pdim);
  ribosum->prnaC = esl_dmatrix_Create(pdim, pdim);
  ribosum->prnaQ = esl_dmatrix_Create(pdim, pdim);
  ribosum->urnaP = esl_dmatrix_Create(udim, udim);
  ribosum->urnaC = esl_dmatrix_Create(udim, udim);
  ribosum->urnaQ = esl_dmatrix_Create(udim, udim);

  ESL_ALLOC(ribosum->prna, sizeof(double) * pdim);
  ESL_ALLOC(ribosum->urna, sizeof(double) * udim);
  ESL_ALLOC(ribosum->bg,   sizeof(double) * udim);

  /* set joints to zero, ready to add counts */
  esl_dmatrix_Set(ribosum->prnaP, 0.0);
  esl_dmatrix_Set(ribosum->urnaP, 0.0);
  esl_vec_DSet(ribosum->prna, pdim, 0.0);
  esl_vec_DSet(ribosum->urna, udim, 0.0);
  esl_vec_DSet(ribosum->bg,   udim, 0.0);
  
  esl_dmatrix_Set(ribosum->prnaC, -1.0);
  esl_dmatrix_Set(ribosum->urnaQ, -1.0);
  esl_dmatrix_Set(ribosum->prnaC, -1.0);
  esl_dmatrix_Set(ribosum->urnaQ, -1.0);

  return ribosum;

 ERROR:
  return NULL;
}

void           
Ribosum_matrix_Destroy(struct ribomatrix_s *ribosum)
{
  if (ribosum) {
    if (ribosum->prnaP) esl_dmatrix_Destroy(ribosum->prnaP);
    if (ribosum->prnaC) esl_dmatrix_Destroy(ribosum->prnaC);
    if (ribosum->prnaQ) esl_dmatrix_Destroy(ribosum->prnaQ);
    if (ribosum->urnaP) esl_dmatrix_Destroy(ribosum->urnaP);
    if (ribosum->urnaC) esl_dmatrix_Destroy(ribosum->urnaC);
    if (ribosum->urnaQ) esl_dmatrix_Destroy(ribosum->urnaQ);
    if (ribosum->prna)  free(ribosum->prna);
    if (ribosum->urna)  free(ribosum->urna);
    if (ribosum->bg)    free(ribosum->bg);
    if (ribosum->name)  free(ribosum->name);
    free(ribosum);
  }
}

int 
Ribosum_matrix_Write(FILE *fp, struct ribomatrix_s *ribosum)
{
  fprintf(fp, "%s\n", ribosum->name);
  if (ribosum->prnaP) esl_dmatrix_Dump(fp, ribosum->prnaP, NULL, NULL);
  if (ribosum->urnaP) esl_dmatrix_Dump(fp, ribosum->urnaP, "ACGU", "ACGU");
  fprintf(fp, "Background frequencies\n");
  if (ribosum->bg)    esl_vec_DDump   (fp, ribosum->bg, ribosum->abc->K, "ACGU");

  fprintf(fp, "Conditionals and Marginals\n");
  if (ribosum->prnaC) esl_dmatrix_Dump(fp, ribosum->prnaC, NULL, NULL);
  if (ribosum->prna)  esl_vec_DDump   (fp, ribosum->prna, ribosum->prnaC->m, NULL);

  if (ribosum->urnaC) esl_dmatrix_Dump(fp, ribosum->urnaC, "ACGU", "ACGU");
  if (ribosum->urna)  esl_vec_DDump   (fp, ribosum->urna, ribosum->urnaC->m, NULL);

  fprintf(fp, "Rates\n");
  if (ribosum->prnaQ) esl_dmatrix_Dump(fp, ribosum->prnaQ, NULL, NULL);
  if (ribosum->urnaQ) esl_dmatrix_Dump(fp, ribosum->urnaQ, "ACGU", "ACGU");

  return eslOK;
}

static int                  
ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float threshid, char *errbuf)
{
  double   wgt;
  double   pid;
  char    *ssi, *ssj;
  int    **ct = NULL;   // ct array for each seq
  int      i, j;
  int      status;
  
  ESL_ALLOC(ct, sizeof(int *) * msa->nseq);
  for (i = 0; i < msa->nseq; i ++) {
    ESL_ALLOC(ct[i], sizeof(int)*msa->alen);
    if      (msa->ss_cons) esl_wuss2ct(msa->ss_cons, msa->alen, ct[i]);
    else if (msa->ss[i])   esl_wuss2ct(msa->ss[i],   msa->alen, ct[i]);
    else                   ESL_XFAIL(eslFAIL, "no ss for sequence %d\n", errbuf);
  }  

  for (i = 0; i < msa->nseq; i++) {
    for (j = 0; j < i; j++) {
      esl_dst_CPairId(msa->aseq[i], msa->aseq[j], &pid, NULL, NULL);
      if (pid < threshid) continue; // not above cutoff

      wgt = msa->wgt[i] * msa->wgt[j]; 
      ssi = (msa->ss_cons)? msa->ss_cons: msa->ss[i];
      ssj = (msa->ss_cons)? msa->ss_cons: msa->ss[j];
      
      prna_add_counts(ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], ssi, ssj, ct[i], ct[j], msa->alen, ribosum->prnaP);
      urna_add_counts(ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], ssi, ssj, ct[i], ct[j], msa->alen, ribosum->urnaP);
      bg_add_counts  (ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], msa->alen, ribosum->bg);
    }
  }

  for (i = 0; i < msa->nseq; i ++) free(ct[i]);
  free(ct);

  return eslOK;

 ERROR:
  return status;
}


static int
prna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *prnaP)
{
  int  pos;
  int  cpos1, cpos2;
  char c1, c2;
  char cc1, cc2;
  int  left, right;

  for (pos = 0; pos < alen; pos++) {
    c1 = seq1[pos];
    c2 = seq2[pos];
    
    if (esl_abc_CIsCanonical(abc, c1) && esl_abc_CIsCanonical(abc, c2)) 
      {
	cpos1 = ct1[pos+1] - 1;
	cpos2 = ct2[pos+1] - 1;
	if (cpos1 != cpos2 || cpos1 < pos || cpos1 == 0) continue;
	    
	cc1 = seq1[cpos1];
	cc2 = seq2[cpos2];
	if (esl_abc_CIsCanonical(abc, cc1) && esl_abc_CIsCanonical(abc, cc2)) 
	  {

	    left  = IDX(esl_abc_DigitizeSymbol(abc, c1),  esl_abc_DigitizeSymbol(abc, c2),  abc->K);
	    right = IDX(esl_abc_DigitizeSymbol(abc, cc1), esl_abc_DigitizeSymbol(abc, cc2), abc->K);

	    prnaP->mx[left][right] += wgt;
	    prnaP->mx[right][left] += wgt;
	  }
      }
  }

  return eslOK;
}

static int
urna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *urnaP)
{
  int  pos;
  char c1, c2;
  
  for (pos = 0; pos < alen; pos++) {
    c1 = seq1[pos];
    c2 = seq2[pos];
    if (ct1[pos+1] != 0 || ct2[pos+1] != 0) continue;
    if (esl_abc_CIsCanonical(abc, c1) && esl_abc_CIsCanonical(abc, c2)) {
      urnaP->mx[esl_abc_DigitizeSymbol(abc, c1)][esl_abc_DigitizeSymbol(abc, c2)] += wgt;
      urnaP->mx[esl_abc_DigitizeSymbol(abc, c2)][esl_abc_DigitizeSymbol(abc, c1)] += wgt;
    }
  }

  return eslOK;
}

static int
bg_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, int alen, double *bg)
{
  int pos;
  int c1, c2;
  
  for (pos = 0; pos < alen; pos++) {
    
    c1 = seq1[pos];
    c2 = seq2[pos];
    
    if (esl_abc_CIsCanonical(abc, c1) && esl_abc_CIsCanonical(abc, c2)) 
      {
	bg[esl_abc_DigitizeSymbol(abc, c1)] += wgt;
	bg[esl_abc_DigitizeSymbol(abc, c2)] += wgt;
      }
  }

  return eslOK;
}

