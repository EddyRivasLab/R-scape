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

#include "ratematrix.h"
#include "ribosum_matrix.h"

static int ribosum_matrix_add_counts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh, char *errbuf);
static int prna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *prnaP);
static int urna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *urnaP);
static int bg_add_counts  (ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, int alen, double *bg);
static int conditionals_from_joints(ESL_DMATRIX *J, ESL_DMATRIX *C, double *marginals, double tol, int verbose, char *errbuf);
static int rate_from_conditionals(ESL_DMATRIX *C, ESL_DMATRIX *Q, double tol, int verbose, char *errbuf);

int
Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, FILE *fp, double tol, int verbose, char *errbuf)
{	
  ESL_DMATRIX *pC = NULL;
  ESL_DMATRIX *uC = NULL;
  int          status;
  
  /* calculate Joints */
  status = Ribosum_matrix_JointsFromMSA(msa, ribosum, thresh1, thresh2, tol, verbose, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed Ribosum_matrix_Joints()");

  /* calculate conditionals and marginals */
  status = Ribosum_matrix_ConditionalsFromJoint(ribosum, tol, verbose, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed Ribosum_matrix_ConditionalsFromJoint()");

  /* calculate rates */
  status = Ribosum_matrix_RateFromConditionals(ribosum, tol, verbose, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed Ribosum_matrix_RateFromConditionals()");

  pC = ratematrix_ConditionalsFromRate(1.0, ribosum->prnaQ, tol, errbuf, verbose);
  uC = ratematrix_ConditionalsFromRate(1.0, ribosum->urnaQ, tol, errbuf, verbose);
  ratematrix_specialDump(pC);
  ratematrix_specialDump(uC);

  ratematrix_Rescale(ribosum->prnaQ, NULL, ribosum->prna);
  ratematrix_Rescale(ribosum->urnaQ, NULL, ribosum->urna);
  if (verbose) Ribosum_matrix_Write(stdout, ribosum); 
  Ribosum_matrix_Write(fp, ribosum);

  if (pC) esl_dmatrix_Destroy(pC);
  if (uC) esl_dmatrix_Destroy(uC);
  return eslOK;
  
 ERROR:
  if (pC) esl_dmatrix_Destroy(pC);
  if (uC) esl_dmatrix_Destroy(uC);
  return status;
}

int
Ribosum_matrix_JointsFromMSA(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, double tol, int verbose, char *errbuf)
{	
  float        sum;
  int          status;
  
  /* calculate the weight BLOSUM-style */
  status = esl_msaweight_BLOSUM(msa, thresh1);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed msaweight_BLOSUM");
  
  /* add the counts */
  status = ribosum_matrix_add_counts(msa, ribosum, thresh2, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed ribosum_matrix_add_counts()");

  /* normalize */
  sum = esl_dmx_Sum(ribosum->prnaJ);
  if (sum > 0.) esl_dmx_Scale(ribosum->prnaJ, 1.0/sum);
  sum = esl_dmx_Sum(ribosum->urnaJ);
  if (sum > 0.) esl_dmx_Scale(ribosum->urnaJ, 1.0/sum);
  sum = esl_vec_DSum(ribosum->bg, ribosum->abc->K);
  if (sum > 0.) esl_vec_DScale(ribosum->bg, ribosum->abc->K, 1.0/sum);

  return eslOK;
  
 ERROR:
  return status;
}

int 
Ribosum_matrix_ConditionalsFromJoint(struct ribomatrix_s *ribosum, double tol, int verbose, char *errbuf)
{
  double sum;
  int    i, j;
  int    status;
  
  for (i = 0; i < ribosum->bprsJ->m; i ++) {
    sum = 0.;
    for (j = 0; j < ribosum->bprsJ->n; j ++) {
      ribosum->bprsC->mx[i][j] = ribosum->bprsJ->mx[i][j];
      sum += ribosum->bprsJ->mx[i][j];
    }
    if (sum > 0) {
      ribosum->bprsM[i] = sum;
      for (j = 0; j < ribosum->bprsJ->n; j ++) ribosum->bprsC->mx[i][j] *= 1.0/sum;
    }
  }
  status = ratematrix_ValidateP(ribosum->bprsC, tol, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "bprsC did not validate");

  for (i = 0; i < ribosum->urnaJ->m; i ++) {
    sum = 0.;
    for (j = 0; j < ribosum->urnaJ->n; j ++) {
      ribosum->urnaC->mx[i][j] = ribosum->urnaJ->mx[i][j];
      sum += ribosum->urnaJ->mx[i][j];
    }
    if (sum > 0) {
      ribosum->urnaM[i] = sum;
      for (j = 0; j < ribosum->urnaJ->n; j ++) ribosum->urnaC->mx[i][j] *= 1.0/sum;
    }
  }
  status = ratematrix_ValidateP(ribosum->urnaC, tol, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "urnaC did not validate");
  
  return eslOK;

 ERROR:
  return status;
}

struct ribomatrix_s *
Ribosum_matrix_Create(ESL_ALPHABET *abc, char *name)
{
  struct ribomatrix_s *ribosum = NULL;
  int    dim = abc->K;
  int    dim2 = abc->K * abc->K;
  int    status;
  
  ESL_ALLOC(ribosum, sizeof(struct ribomatrix_s));
  ribosum->abc = abc;
  ribosum->name = NULL;
  if (name) ribosum->name = name;
  
  ribosum->bprsJ = esl_dmatrix_Create(dim2, dim2);
  ribosum->bprsC = esl_dmatrix_Create(dim2, dim2);
  ribosum->bprsQ = esl_dmatrix_Create(dim2, dim2);
  ESL_ALLOC(ribosum->bprsM, sizeof(double) * dim2);

  ribosum->prnaJ = esl_dmatrix_Create(dim, dim);
  ribosum->prnaC = esl_dmatrix_Create(dim, dim);
  ribosum->prnaQ = esl_dmatrix_Create(dim, dim);
  ESL_ALLOC(ribosum->prnaM, sizeof(double) * dim);

  ribosum->urnaJ = esl_dmatrix_Create(dim, dim);
  ribosum->urnaC = esl_dmatrix_Create(dim, dim);
  ribosum->urnaQ = esl_dmatrix_Create(dim, dim);
  ESL_ALLOC(ribosum->urnaM, sizeof(double) * dim);

  ribosum->xrnaJ = esl_dmatrix_Create(dim, dim);
  ribosum->xrnaC = esl_dmatrix_Create(dim, dim);
  ribosum->xrnaQ = esl_dmatrix_Create(dim, dim);
  ESL_ALLOC(ribosum->xrnaM, sizeof(double) * dim);

  ESL_ALLOC(ribosum->bg, sizeof(double) * dim);
  
  /* set joints to zero, ready to add counts */
  esl_dmatrix_Set(ribosum->bprsJ, 0.0);
  esl_dmatrix_Set(ribosum->prnaJ, 0.0);
  esl_dmatrix_Set(ribosum->urnaJ, 0.0);
  esl_dmatrix_Set(ribosum->xrnaJ, 0.0);
  esl_vec_DSet(ribosum->bg, dim,  0.0);
  
  /* conditionals and marginals */
  esl_dmatrix_Set(ribosum->bprsC, 0.0);
  esl_dmatrix_Set(ribosum->prnaC, 0.0);
  esl_dmatrix_Set(ribosum->urnaC, 0.0);
  esl_dmatrix_Set(ribosum->xrnaC, 0.0);
  esl_vec_DSet(ribosum->bprsM, dim2, 0.0);
  esl_vec_DSet(ribosum->prnaM, dim,  0.0);
  esl_vec_DSet(ribosum->urnaM, dim,  0.0);
  esl_vec_DSet(ribosum->xrnaM, dim,  0.0);
  
  /* rates */
  esl_dmatrix_Set(ribosum->bprsQ, 0.0);
  esl_dmatrix_Set(ribosum->prnaQ, 0.0);
  esl_dmatrix_Set(ribosum->urnaQ, 0.0);
  esl_dmatrix_Set(ribosum->xrnaQ, 0.0);
  
  return ribosum;
  
 ERROR:
  return NULL;
}

void           
Ribosum_matrix_Destroy(struct ribomatrix_s *ribosum)
{
  if (ribosum) {
    if (ribosum->bprsJ) esl_dmatrix_Destroy(ribosum->bprsJ);
    if (ribosum->bprsC) esl_dmatrix_Destroy(ribosum->bprsC);
    if (ribosum->bprsQ) esl_dmatrix_Destroy(ribosum->bprsQ);

    if (ribosum->prnaJ) esl_dmatrix_Destroy(ribosum->prnaJ);
    if (ribosum->prnaC) esl_dmatrix_Destroy(ribosum->prnaC);
    if (ribosum->prnaQ) esl_dmatrix_Destroy(ribosum->prnaQ);
    if (ribosum->urnaJ) esl_dmatrix_Destroy(ribosum->urnaJ);
    if (ribosum->urnaC) esl_dmatrix_Destroy(ribosum->urnaC);
    if (ribosum->urnaQ) esl_dmatrix_Destroy(ribosum->urnaQ);
    if (ribosum->xrnaJ) esl_dmatrix_Destroy(ribosum->xrnaJ);
    if (ribosum->xrnaC) esl_dmatrix_Destroy(ribosum->xrnaC);
    if (ribosum->xrnaQ) esl_dmatrix_Destroy(ribosum->xrnaQ);
    if (ribosum->bprsM)  free(ribosum->bprsM);
    if (ribosum->prnaM)  free(ribosum->prnaM);
    if (ribosum->urnaM)  free(ribosum->urnaM);
    if (ribosum->xrnaM)  free(ribosum->xrnaM);
    if (ribosum->bg)     free(ribosum->bg);
    if (ribosum->name)   free(ribosum->name);
    free(ribosum);
  }
}

int 
Ribosum_matrix_RateFromConditionals(struct ribomatrix_s *ribosum, double tol, int verbose, char *errbuf)
{
  int status;

  status = rate_from_conditionals(ribosum->bprsC, ribosum->bprsQ, tol, verbose, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed bprs RateFromConditionals()");

  status = rate_from_conditionals(ribosum->prnaC, ribosum->prnaQ, tol, verbose, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed prna RateFromConditionals()");

  status = rate_from_conditionals(ribosum->urnaC, ribosum->urnaQ, tol, verbose, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed urna RateFromConditionals()");

  status = rate_from_conditionals(ribosum->xrnaC, ribosum->xrnaQ, tol, verbose, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed xrna RateFromConditionals()");

  return eslOK;

 ERROR:
  return status;
}

int
rate_from_conditionals(ESL_DMATRIX *C, ESL_DMATRIX *Q, double tol, int verbose, char *errbuf)
{
  int i;
  int status;

  if (C == NULL) return eslOK;

  status = ratematrix_ValidatePLog(C, tol, errbuf); /* Make sure that we can take the log of P */
  if (status != eslOK) { printf("failed to validate PLog\n%s\n", errbuf); goto ERROR; }
    
  dmx_Log(C, Q, tol);                      /* take the log */
  for (i = 0; i < Q->n; i++)               /* regularize in case some entries are not good */
    ratematrix_QOGRegularization(Q->mx[i], Q->m, i, tol, errbuf);
  status = ratematrix_ValidateQ(Q, tol, errbuf);    /* Make sure Q is a rate matrix */
  if (status != eslOK) { printf("failed to validate rate Q\n%s\n", errbuf); goto ERROR; }

  return eslOK;

 ERROR:
  return status;
}

int
conditionals_from_joints(ESL_DMATRIX *J, ESL_DMATRIX *C, double *marginals, double tol, int verbose, char *errbuf)
{
  double sum;
  int    i, j;
  int    status;
  
  for (i = 0; i < J->m; i ++) {
    sum = 0.;
    for (j = 0; j < J->n; j ++) {
      C->mx[i][j] = J->mx[i][j];
      sum += J->mx[i][j];
    }
    if (sum > 0) {
      marginals[i] = sum;
      for (j = 0; j < J->n; j ++) C->mx[i][j] *= 1.0/sum;
    }
  }
  status = ratematrix_ValidateP(C, tol, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "Conditionals did not validate");

  return eslOK;
}

int 
Ribosum_matrix_Write(FILE *fp, struct ribomatrix_s *ribosum)
{
  fprintf(fp, "%s\n", ribosum->name);
  if (ribosum->bprsJ) esl_dmatrix_Dump(fp, ribosum->bprsJ, NULL, NULL);
  if (ribosum->prnaJ) esl_dmatrix_Dump(fp, ribosum->prnaJ, "ACGU", "ACGU");
  if (ribosum->urnaJ) esl_dmatrix_Dump(fp, ribosum->urnaJ, "ACGU", "ACGU");
  if (ribosum->xrnaJ) esl_dmatrix_Dump(fp, ribosum->xrnaJ, "ACGU", "ACGU");
  fprintf(fp, "Background frequencies\n");
  if (ribosum->bg)    esl_vec_DDump   (fp, ribosum->bg, ribosum->abc->K, "ACGU");

  fprintf(fp, "Conditionals and Marginals\n");
  if (ribosum->bprsC) esl_dmatrix_Dump(fp, ribosum->bprsC, NULL, NULL);
  if (ribosum->bprsM) esl_vec_DDump   (fp, ribosum->bprsM, ribosum->bprsC->m, NULL);
  if (ribosum->prnaC) esl_dmatrix_Dump(fp, ribosum->prnaC, "ACGU", "ACGU");
  if (ribosum->prnaM) esl_vec_DDump   (fp, ribosum->prnaM, ribosum->prnaC->m, NULL);
  if (ribosum->urnaC) esl_dmatrix_Dump(fp, ribosum->urnaC, "ACGU", "ACGU");
  if (ribosum->urnaM) esl_vec_DDump   (fp, ribosum->urnaM, ribosum->urnaC->m, NULL);
  if (ribosum->xrnaC) esl_dmatrix_Dump(fp, ribosum->xrnaC, "ACGU", "ACGU");
  if (ribosum->xrnaM) esl_vec_DDump   (fp, ribosum->xrnaM, ribosum->xrnaC->m, NULL);

  fprintf(fp, "Rates\n");
  if (ribosum->bprsQ) esl_dmatrix_Dump(fp, ribosum->bprsQ, NULL, NULL);
  if (ribosum->prnaQ) esl_dmatrix_Dump(fp, ribosum->prnaQ, "ACGU", "ACGU");
  if (ribosum->urnaQ) esl_dmatrix_Dump(fp, ribosum->urnaQ, "ACGU", "ACGU");
  if (ribosum->xrnaQ) esl_dmatrix_Dump(fp, ribosum->xrnaQ, "ACGU", "ACGU");

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
      
      prna_add_counts(ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], ssi, ssj, ct[i], ct[j], msa->alen, ribosum->prnaJ);
      urna_add_counts(ribosum->abc, wgt, msa->aseq[i], msa->aseq[j], ssi, ssj, ct[i], ct[j], msa->alen, ribosum->urnaJ);
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
prna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *prnaJ)
{
  int  pos;
  int  cpos1, cpos2;
  char c1, c2;
  char cc1, cc2;
  int  left, right;
  int  pair1, pair2;

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

	    pair1 = IDX(esl_abc_DigitizeSymbol(abc, c1), esl_abc_DigitizeSymbol(abc, cc1), abc->K);
	    pair2 = IDX(esl_abc_DigitizeSymbol(abc, c2), esl_abc_DigitizeSymbol(abc, cc2), abc->K);

	    //prnaJ->mx[left][right] += wgt;
	    prnaJ->mx[pair1][pair2] += wgt;
	    prnaJ->mx[pair2][pair1] += wgt;
	  }
      }
  }

  return eslOK;
}

static int
urna_add_counts(ESL_ALPHABET *abc, double wgt, char *seq1, char *seq2, char *ss1, char *ss2, int *ct1, int *ct2, int alen, ESL_DMATRIX *urnaJ)
{
  int  pos;
  char c1, c2;
  
  for (pos = 0; pos < alen; pos++) {
    c1 = seq1[pos];
    c2 = seq2[pos];
    if (ct1[pos+1] != 0 || ct2[pos+1] != 0) continue;
    if (esl_abc_CIsCanonical(abc, c1) && esl_abc_CIsCanonical(abc, c2)) {
      urnaJ->mx[esl_abc_DigitizeSymbol(abc, c1)][esl_abc_DigitizeSymbol(abc, c2)] += wgt;
      urnaJ->mx[esl_abc_DigitizeSymbol(abc, c2)][esl_abc_DigitizeSymbol(abc, c1)] += wgt;
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

