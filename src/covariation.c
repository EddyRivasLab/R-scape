/* mutualinfo.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covariation.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"

static int is_wc(int x, int y);
static int mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				 ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);

int                 
Mutual_Calculate(ESL_MSA *msa, int *msamap, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, COVTYPE covtype, COVCLASS covclass, 
		 int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, double tol, int verbose, char *errbuf)
{
  HIT  *hitlist = NULL;
  int   status;

  status = Mutual_Probs(msa, T, ribosum, mi, method, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
 
   switch(covtype) {
   case CHIa: 
     status = Mutual_CalculateCHI         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case CHIp:
     status = Mutual_CalculateCHI         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;  
     break;
   case CHI: 
     status = Mutual_CalculateCHI         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;     
     break;
   case GTa: 
     status = Mutual_CalculateGT          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC,  TRUE, &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
      break;
   case GTp: 
     status = Mutual_CalculateGT          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC,  TRUE, &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case GT: 
     status = Mutual_CalculateGT          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,       TRUE, &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case MIa: 
     status = Mutual_CalculateMI          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MIp: 
     status = Mutual_CalculateMI          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MI: 
     status = Mutual_CalculateMI          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case MIra: 
     status = Mutual_CalculateMIr         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(         mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC,  TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MIrp:
     status = Mutual_CalculateMIr         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;  
     break;
   case MIr: 
     status = Mutual_CalculateMIr         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE, &hitlist,  tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case MIga: 
     status = Mutual_CalculateMIg         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MIgp:
     status = Mutual_CalculateMIg         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;  
     break;
   case MIg: 
     status = Mutual_CalculateMIg          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,     TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case OMESa: 
     status = Mutual_CalculateOMES        (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case OMESp: 
     status = Mutual_CalculateOMES        (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case OMES: 
     status = Mutual_CalculateOMES        (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
      break;
   case COVALL: 
     status = Mutual_CalculateCHI         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = Mutual_CalculateCHI         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     
     status = Mutual_CalculateOMES        (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = Mutual_CalculateOMES        (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     
     status = Mutual_CalculateGT          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = Mutual_CalculateGT          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     
     status = Mutual_CalculateMI          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = Mutual_CalculateMI          (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     
     status = Mutual_CalculateMIr         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = Mutual_CalculateMIr         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     
     status = Mutual_CalculateMIg         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      TRUE,  FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, APC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = Mutual_CalculateMIg         (covclass, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs,      FALSE, FALSE,    tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = Mutual_CalculateCOVCorrected(          mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, ASC, TRUE,  &hitlist, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   default:
     ESL_XFAIL(eslFAIL, errbuf, "wrong covariation type\n");
     break;
   }
   fprintf(sumfp, "\n");   
      
   status = Mutual_FormatR2R(msa, hitlist, verbose, errbuf);
   if  (status != eslOK) goto ERROR;

   free(hitlist);
  return eslOK;
  
 ERROR:
  if (hitlist) free(hitlist);
  return status;
}


int                 
Mutual_Probs(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf)
{
  int i, j;
  int x, y;
  int K = msa->abc->K;
  int status;

  switch(method) {
  case NAIVE:
    status = Mutual_NaivePP(msa, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    break;
  case PHYLO:
    status = Mutual_PostOrderPP(msa, T, ribosum, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    break;
  case DCA:
    break;
  case AKMAEV:
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "bad method option");
  }
  
  /* pm are the marginals */
  for (i = 0; i < mi->alen; i ++) {
    esl_vec_DSet(mi->pm[i], K, 0.0);

    for (j = 0; j < mi->alen; j ++)     
      for (x = 0; x < K; x ++) 
	for (y = 0; y < K; y ++) 
	  mi->pm[i][x] += mi->pp[i][j][IDX(x,y,K)];
    esl_vec_DNorm(mi->pm[i], K);
    status = esl_vec_DValidate(mi->pm[i], K, tol, errbuf);
    if (status != eslOK) {
      printf("pm[%d]\n", i);
      esl_vec_DDump(stdout, mi->pm[i], K, "ACGU");
      goto ERROR;
    }
  }

  status = Mutual_ValidateProbs(mi, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  if (verbose) {
    for (i = 0; i < mi->alen-1; i ++) {
      for (j = i+1; j < mi->alen; j ++) {
	if (i==5&&j==118) {
	  printf("pp[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->pp[i][j][IDX(x,y,K)]);
	    }
	  printf("\n");
	  printf("pm*pm[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->pm[i][x]*mi->pm[j][y]);
	    }
	  printf("\n");
	  printf("Cp[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ",  mi->nseff[i][j]*mi->pp[i][j][IDX(x,y,K)]);
	    }
	  printf("\n");
	  printf("\n");
	  printf("Ex[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->nseff[i][j]*mi->pm[i][x]*mi->pm[j][y]);
	    }
	  printf("\n");
	}
      }
      if (i==5) esl_vec_DDump(stdout, mi->pm[i], K, NULL);
    }
  }
  
  return status;

 ERROR:
  return status;

}


int 
Mutual_ValidateProbs(struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int    i, j;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  /* pp validation */
  for (i = 0; i < mi->alen-1; i ++) 
    for (j = i+1; j < mi->alen; j ++) {
      status = esl_vec_DValidate(mi->pp[i][j], K*K, tol, errbuf);
      if (status != eslOK) {
	printf("pp[%d][%d]\n", i, j);
	esl_vec_DDump(stdout, mi->pp[i][j], K*K, NULL);
	ESL_XFAIL(eslFAIL, errbuf, "pp validation failed");
      }
    }
  
  /* pm validation */
  for (i = 0; i < mi->alen; i ++) {
    status = esl_vec_DValidate(mi->pm[i], K, tol, errbuf);
    if (status != eslOK) {
      printf("pm[%d]\n", i);
      esl_vec_DDump(stdout, mi->pm[i], K, NULL);
      ESL_XFAIL(eslFAIL, errbuf, "pm validation failed");
    }
  }
  
  return eslOK;
  
 ERROR:
  return status;
}

int                 
Mutual_CalculateCHI(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		    int analyze, HIT **ret_hitlist, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = Mutual_CalculateCHI_C16  (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case C2:
    status = Mutual_CalculateCHI_C2   (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = Mutual_CalculateCHI_C2 (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    else
      status = Mutual_CalculateCHI_C16(mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
}
  
  if (verbose) {
    printf("CHI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("CHI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }

  if (analyze) status = Mutual_SignificantPairs_Ranking(ret_hitlist, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}


int                 
Mutual_CalculateCHI_C16(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
			double tol, int verbose, char *errbuf)
{
  double chi;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, CHI, C16);
  
  // CHI
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      chi  = 0.0;

      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  exp = (double)mi->nseff[i][j] * mi->pm[i][x] * mi->pm[j][y];
	  obs = (double)mi->nseff[i][j] * mi->pp[i][j][IDX(x,y,K)];
	  chi += (exp > 0.)? (obs-exp) * (obs-exp) / exp : 0.0 ;
	}	  
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = chi;
      if (chi < mi->minCOV) mi->minCOV = chi;
      if (chi > mi->maxCOV) mi->maxCOV = chi;
    }
  
  return status;
}

int                 
Mutual_CalculateCHI_C2(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		       double tol, int verbose, char *errbuf)
{
  double chi;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  double obs_wc, obs_nwc;
  double exp_wc, exp_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, CHI, C2);
  
   // CHI
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      chi = 0.0;
      pij_wc = pij_nwc = 0.;
      qij_wc = qij_nwc = 0.;
      
      for (x = 0; x < K; x ++) 
	for (y = 0; y < K; y ++) {
	  if (is_wc(x,y)) {
	    pij_wc  += mi->pp[i][j][IDX(x,y,K)];
	    qij_wc  += mi->pm[i][x] * mi->pm[j][y];
	  }
	  else {
	    pij_nwc += mi->pp[i][j][IDX(x,y,K)];
	    qij_nwc += mi->pm[i][x] * mi->pm[j][y];
	  }
	}	  
      
      exp_wc  = (double)mi->nseff[i][j] * qij_wc;
      exp_nwc = (double)mi->nseff[i][j] * qij_nwc;
      obs_wc  = (double)mi->nseff[i][j] * pij_wc;
      obs_nwc = (double)mi->nseff[i][j] * pij_nwc;

      chi += (exp_wc  > 0.)? (obs_wc -exp_wc)  * (obs_wc -exp_wc)  / exp_wc  : 0.0 ;
      chi += (exp_nwc > 0.)? (obs_nwc-exp_nwc) * (obs_nwc-exp_nwc) / exp_nwc : 0.0 ;
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = chi;
      if (chi < mi->minCOV) mi->minCOV = chi;
      if (chi > mi->maxCOV) mi->maxCOV = chi;
    }

  return status;
}


int                 
Mutual_CalculateOMES(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		     int analyze, HIT **ret_hitlist, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status;
  
  switch (covclass) {
  case C16:
    status = Mutual_CalculateOMES_C16  (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case C2:
    status = Mutual_CalculateOMES_C2   (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = Mutual_CalculateOMES_C2 (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    else
      status = Mutual_CalculateOMES_C16(mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("OMES[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("OMES[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(ret_hitlist, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;
  
 ERROR:
  return status;
}

int                 
Mutual_CalculateOMES_C16(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
			 double tol, int verbose, char *errbuf)
{
  double omes;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, OMES, C16);
  
  // OMES
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      omes  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  exp = (double)mi->nseff[i][j] * mi->pm[i][x] * mi->pm[j][y];
	  obs = (double)mi->nseff[i][j] * mi->pp[i][j][IDX(x,y,K)];
	  omes += (exp > 0.)? (obs-exp) * (obs-exp) / (double)mi->nseff[i][j] : 0.0;
	}	  
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = omes;
      if (omes < mi->minCOV) mi->minCOV = omes;
      if (omes > mi->maxCOV) mi->maxCOV = omes;
    }
  
  return status;
}

int                 
Mutual_CalculateOMES_C2(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
			double tol, int verbose, char *errbuf)
{
  double omes;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  double obs_wc, obs_nwc;
  double exp_wc, exp_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, OMES, C2);
  
  // OMES
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      omes  = 0.0;
      pij_wc = pij_nwc = 0.;
      qij_wc = qij_nwc = 0.;
      
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  if (is_wc(x,y)) {
	    pij_wc  += mi->pp[i][j][IDX(x,y,K)];
	    qij_wc  += mi->pm[i][x] * mi->pm[j][y];
	  }
	  else {
	    pij_nwc += mi->pp[i][j][IDX(x,y,K)];
	    qij_nwc += mi->pm[i][x] * mi->pm[j][y];
	  }
	}	  
      
      exp_wc  = (double)mi->nseff[i][j] * qij_wc;
      exp_nwc = (double)mi->nseff[i][j] * qij_nwc;
      obs_wc  = (double)mi->nseff[i][j] * pij_wc;
      obs_nwc = (double)mi->nseff[i][j] * pij_nwc;

      omes += (exp_wc  > 0.)? (obs_wc -exp_wc)  * (obs_wc -exp_wc)  / (double)mi->nseff[i][j] : 0.0;
      omes += (exp_nwc > 0.)? (obs_nwc-exp_nwc) * (obs_nwc-exp_nwc) / (double)mi->nseff[i][j] : 0.0;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = omes;
      if (omes < mi->minCOV) mi->minCOV = omes;
      if (omes > mi->maxCOV) mi->maxCOV = omes;
    }


  return status;
}

int                 
Mutual_CalculateGT(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		   int analyze, HIT **ret_hitlist, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = Mutual_CalculateGT_C16  (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case C2:
    status = Mutual_CalculateGT_C2   (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh) {
      status = Mutual_CalculateGT_C2 (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    }
    else {
      status = Mutual_CalculateGT_C16(mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    }
    break;
  }
  
  if (verbose) {
    printf("GT[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("GT[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(ret_hitlist, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}


int                 
Mutual_CalculateGT_C16(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		       double tol, int verbose, char *errbuf)
{
  double gt;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, GT, C16);
  
  // GT
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      gt  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  exp = (double)mi->nseff[i][j] * mi->pm[i][x] * mi->pm[j][y];
	  obs = (double)mi->nseff[i][j] * mi->pp[i][j][IDX(x,y,K)];
	  gt += (exp > 0.) ? obs * log (obs / exp) : 0.0;
	}	  
      gt *= 2.0;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = gt;
      if (gt < mi->minCOV) mi->minCOV = gt;
      if (gt > mi->maxCOV) mi->maxCOV = gt;
    }
  
   return status;
}

int                 
Mutual_CalculateGT_C2(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		      double tol, int verbose, char *errbuf)
{
  double gt;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  double obs_wc, obs_nwc;
  double exp_wc, exp_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, GT, C2);
  
  // GT
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      gt = 0.0;
      pij_wc = pij_nwc = 0.;
      qij_wc = qij_nwc = 0.;
 
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  if (is_wc(x,y)) {
	    pij_wc  += mi->pp[i][j][IDX(x,y,K)];
	    qij_wc  += mi->pm[i][x] * mi->pm[j][y];
	  }
	  else {
	    pij_nwc += mi->pp[i][j][IDX(x,y,K)];
	    qij_nwc += mi->pm[i][x] * mi->pm[j][y];
	  }
	}

      exp_wc  = (double)mi->nseff[i][j] * qij_wc;
      exp_nwc = (double)mi->nseff[i][j] * qij_nwc;
      obs_wc  = (double)mi->nseff[i][j] * pij_wc;
      obs_nwc = (double)mi->nseff[i][j] * pij_nwc;
      
      gt += (exp_wc  > 0.) ? obs_wc  * log (obs_wc  / exp_wc)  : 0.0;
      gt += (exp_nwc > 0.) ? obs_nwc * log (obs_nwc / exp_nwc) : 0.0;
      gt *= 2.0;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = gt;
      if (gt < mi->minCOV) mi->minCOV = gt;
      if (gt > mi->maxCOV) mi->maxCOV = gt;
    }

  return status;
}


int                 
Mutual_CalculateMI(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		   int analyze, HIT **ret_hitlist, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = Mutual_CalculateMI_C16  (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case C2:
    status = Mutual_CalculateMI_C2   (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = Mutual_CalculateMI_C2 (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    else
      status = Mutual_CalculateMI_C16(mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("MI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(ret_hitlist, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}

int                 
Mutual_CalculateMI_C16(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		       double tol, int verbose, char *errbuf)
{
  double mutinf;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MI, C16);
  
  // MI
  for (i = 0; i < mi->alen-1; i++) {
    for (j = i+1; j < mi->alen; j++) {
      mutinf  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  mutinf += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->pm[i][x]) - log(mi->pm[j][y]) );
	}	  
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = mutinf;
      if (mutinf < mi->minCOV) mi->minCOV = mutinf;
      if (mutinf > mi->maxCOV) mi->maxCOV = mutinf;
    }
  }
  
  return status;
}

int                 
Mutual_CalculateMI_C2(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		      double tol, int verbose, char *errbuf)
{
  double mutinf;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MI, C2);
  
  // MI
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      mutinf = 0.0;
      pij_wc = pij_nwc = 0.;
      qij_wc = qij_nwc = 0.;
      
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  if (is_wc(x,y)) {
	    pij_wc  += mi->pp[i][j][IDX(x,y,K)];
	    qij_wc  += mi->pm[i][x] * mi->pm[j][y];
	  }
	  else {
	    pij_nwc += mi->pp[i][j][IDX(x,y,K)];
	    qij_nwc += mi->pm[i][x] * mi->pm[j][y];
	  }
	}
      
      mutinf += pij_wc  * ( log(pij_wc)  - log(qij_wc)  );
      mutinf += pij_nwc * ( log(pij_nwc) - log(qij_nwc) );
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = mutinf;
      if (mutinf < mi->minCOV) mi->minCOV = mutinf;
      if (mutinf > mi->maxCOV) mi->maxCOV = mutinf; 
    }

  return status;
}


int                 
Mutual_CalculateMIr(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		    int analyze, HIT **ret_hitlist, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = Mutual_CalculateMIr_C16  (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case C2:
    status = Mutual_CalculateMIr_C2   (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = Mutual_CalculateMIr_C2 (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    else
      status = Mutual_CalculateMIr_C16(mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("MIr[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MIr[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(ret_hitlist, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}

int                 
Mutual_CalculateMIr_C16(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs,
			double tol, int verbose, char *errbuf)
{
  double mutinf, HH;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MIr, C16);
  
  //MIr
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      HH  = 0.0;
      mutinf  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  HH -= mi->pp[i][j][IDX(x,y,K)] * log(mi->pp[i][j][IDX(x,y,K)]);
	  mutinf += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->pm[i][x]) - log(mi->pm[j][y]) );
	}	  
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = (HH > 0.0)? mutinf/HH : 0.0;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
    }
  
  return status;
}

int                 
Mutual_CalculateMIr_C2(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		       double tol, int verbose, char *errbuf)
{
  double mutinf, HH;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MIr, C2);
  
  // MIr
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      HH = mutinf = 0.0;
      pij_wc = pij_nwc = 0.;
      qij_wc = qij_nwc = 0.;
      
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  if (is_wc(x,y)) {
	    pij_wc  += mi->pp[i][j][IDX(x,y,K)];
	    qij_wc  += mi->pm[i][x] * mi->pm[j][y];
	  }
	  else {
	    pij_nwc += mi->pp[i][j][IDX(x,y,K)];
	    qij_nwc += mi->pm[i][x] * mi->pm[j][y];
	  }
	}
      
      HH -= pij_wc  * log(pij_wc);
      HH -= pij_nwc * log(pij_nwc);

      mutinf += pij_wc  * ( log(pij_wc)  - log(qij_wc)  );
      mutinf += pij_nwc * ( log(pij_nwc) - log(qij_nwc) );
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = (HH > 0.0)? mutinf/HH : 0.0;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j]; 
    }

  return status;
}


int                 
Mutual_CalculateMIg(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		    int analyze, HIT **ret_hitlist, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = Mutual_CalculateMIg_C16  (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case C2:
    status = Mutual_CalculateMIg_C2   (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = Mutual_CalculateMIg_C2 (mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    else
      status = Mutual_CalculateMIg_C16(mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, tol, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("MIg[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MIg[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(ret_hitlist, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}

int                 
Mutual_CalculateMIg_C16(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs,
			double tol, int verbose, char *errbuf)
{
  double mutinf;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MIg, C16);
  
  //MIg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      mutinf  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  mutinf += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->pm[i][x]) - log(mi->pm[j][y]) );
	}
	  
      /* the negative correction of gaps */
      mutinf -= (mi->nseff[i][j] > 0)? (double)mi->ngap[i][j] / (double)mi->nseff[i][j] : 0.0;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = mutinf;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
    }
  
  return status;
}

int                 
Mutual_CalculateMIg_C2(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
		       double tol, int verbose, char *errbuf)
{
  double mutinf;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MIg, C2);
  
  // MIg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      mutinf = 0.0;
      pij_wc = pij_nwc = 0.;
      qij_wc = qij_nwc = 0.;
      
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  if (is_wc(x,y)) {
	    pij_wc  += mi->pp[i][j][IDX(x,y,K)];
	    qij_wc  += mi->pm[i][x] * mi->pm[j][y];
	  }
	  else {
	    pij_nwc += mi->pp[i][j][IDX(x,y,K)];
	    qij_nwc += mi->pm[i][x] * mi->pm[j][y];
	  }
	}
      
      mutinf += pij_wc  * ( log(pij_wc)  - log(qij_wc)  );
      mutinf += pij_nwc * ( log(pij_nwc) - log(qij_nwc) );
      
      /* the negative correction of gaps */
      mutinf -= (mi->nseff[i][j] > 0)? (double)mi->ngap[i][j] / (double)mi->nseff[i][j] : 0.0;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = mutinf;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j]; 
    }

  return status;
}


int                 
Mutual_CalculateCOVCorrected(struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
			     CORRTYPE corrtype, int analyze, HIT **ret_hitlist, double tol, int verbose, char *errbuf)
{
  char        *covtype = NULL;
  ESL_DMATRIX *COV  = NULL;
  double      *COVx = NULL;
  double       COVavg = 0.0;
  int          i, j;
  int          status = eslOK;
  
  Mutual_COVTYPEString(&covtype, mi->type, errbuf);

  switch(corrtype) {
  case APC: esl_sprintf(&covtype, "%sp", covtype); break;
  case ASC: esl_sprintf(&covtype, "%sa", covtype); break;
  default:  
    ESL_XFAIL(eslFAIL, errbuf, "wrong correction type\n");
    break;
  }
  COV = esl_dmatrix_Clone(mi->COV);
  
  Mutual_String2COVTYPE(covtype, &mi->type, errbuf);
  Mutual_ReuseCOV(mi, mi->type, mi->class);  
 
  // COVavg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) 
      COVavg += COV->mx[i][j];
  if (mi->alen > 1) COVavg /= (double)mi->alen * ((double)mi->alen-1.);
  COVavg *= 2.;

  //COVx
  ESL_ALLOC(COVx, sizeof(double) * mi->alen);
  for (i = 0; i < mi->alen; i++) {
    COVx[i] = 0.0;
    for (j = 0; j < mi->alen; j++) {
      if (j != i) COVx[i] += COV->mx[i][j];
    }
    if (mi->alen > 1) COVx[i] /= (double)mi->alen-1.;
  }

  //COVp
  for (i = 0; i < mi->alen; i++) 
    for (j = 0; j < mi->alen; j++) {

      if (corrtype == APC) 
	mi->COV->mx[i][j] = (COVavg != 0.0)? COV->mx[i][j] - COVx[i] * COVx[j] / COVavg : 0.0;
      else if (corrtype == ASC) 
	mi->COV->mx[i][j] = COV->mx[i][j] - (COVx[i] + COVx[j] - COVavg); 
      else 
	ESL_XFAIL(eslFAIL, errbuf, "wrong correction type\n");

      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
    }

  if (verbose) {
    printf("%s-[%f,%f] \n", covtype,  mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("%s-[%d][%d] = %f | MI %f | COVx %f COVy %f | COVavg %f\n", 
				 covtype, i, j, mi->COV->mx[i][j], COV->mx[i][j], COVx[i], COVx[j], COVavg);
      } 
  }

  if (analyze) status = Mutual_SignificantPairs_Ranking(ret_hitlist, mi, msamap, ct, rocfp, sumfp, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  if (covtype) free(covtype);
  esl_dmatrix_Destroy(COV);
  free(COVx);
  return status;

 ERROR:
  if (covtype) free(covtype);
  if (COV)     esl_dmatrix_Destroy(COV);
  if (COVx)    free(COVx);
  return status;
}

int 
Mutual_COVTYPEString(char **ret_covtype, COVTYPE type, char *errbuf)
{
  int status;

  switch(type) {
  case CHI:   esl_sprintf(ret_covtype, "CHI");   break;
  case GT:    esl_sprintf(ret_covtype, "GT");    break;
  case OMES:  esl_sprintf(ret_covtype, "OMES");  break;
  case MI:    esl_sprintf(ret_covtype, "MI");    break;
  case MIr:   esl_sprintf(ret_covtype, "MIr");   break; 
  case MIg:   esl_sprintf(ret_covtype, "MIg");   break; 

  case CHIp:  esl_sprintf(ret_covtype, "CHIp");  break;
  case GTp:   esl_sprintf(ret_covtype, "GTp");   break;
  case OMESp: esl_sprintf(ret_covtype, "OMESp"); break;
  case MIp:   esl_sprintf(ret_covtype, "MIp");   break;
  case MIrp:  esl_sprintf(ret_covtype, "MIrp");  break; 
  case MIgp:  esl_sprintf(ret_covtype, "MIgp");  break; 

  case CHIa:  esl_sprintf(ret_covtype, "CHIa");  break;
  case GTa:   esl_sprintf(ret_covtype, "GTa");   break;
  case OMESa: esl_sprintf(ret_covtype, "OMESa"); break;
  case MIa:   esl_sprintf(ret_covtype, "MIa");   break;
  case MIra:  esl_sprintf(ret_covtype, "MIra");  break; 
  case MIga:  esl_sprintf(ret_covtype, "MIga");  break; 

  default: ESL_XFAIL(eslFAIL, errbuf, "wrong COVTYPE");
  }

  return eslOK;
  
 ERROR:
  return status;
}

int 
Mutual_String2COVTYPE(char *covtype, COVTYPE *ret_type, char *errbuf)
{
  COVTYPE type;
  int     status;

  if      (!esl_strcmp(covtype, "CHI"))    type = CHI;
  else if (!esl_strcmp(covtype, "GT"))     type = GT;
  else if (!esl_strcmp(covtype, "OMES"))   type = OMES;
  else if (!esl_strcmp(covtype, "MI"))     type = MI;
  else if (!esl_strcmp(covtype, "MIr"))    type = MIr;
  else if (!esl_strcmp(covtype, "MIg"))    type = MIg;
  else if (!esl_strcmp(covtype, "CHIp"))   type = CHIp;
  else if (!esl_strcmp(covtype, "GTp"))    type = GTp;
  else if (!esl_strcmp(covtype, "OMESp"))  type = OMESp;
  else if (!esl_strcmp(covtype, "MIp"))    type = MIp;
  else if (!esl_strcmp(covtype, "MIrp"))   type = MIrp;
  else if (!esl_strcmp(covtype, "MIgp"))   type = MIgp;
  else if (!esl_strcmp(covtype, "CHIa"))   type = CHIa;
  else if (!esl_strcmp(covtype, "GTa"))    type = GTa;
  else if (!esl_strcmp(covtype, "OMESa"))  type = OMESa;
  else if (!esl_strcmp(covtype, "MIa"))    type = MIa;
  else if (!esl_strcmp(covtype, "MIra"))   type = MIra;
  else if (!esl_strcmp(covtype, "MIga"))   type = MIga;
  else
    ESL_XFAIL(eslFAIL, errbuf, "wrong COVTYPE %s", covtype);

  *ret_type = type;
  return eslOK;
  
 ERROR:
  return status;
}

struct mutual_s *
Mutual_Create(int64_t alen, int64_t nseq, int ishuffled, int nseqthresh, ESL_ALPHABET *abc)
{
  struct mutual_s *mi = NULL;
  int              K  = abc->K;
  int              K2 = K * K;
  int              i, j;
  int              status;
  
  ESL_ALLOC(mi, sizeof(struct mutual_s));
  mi->alen       = alen;
  mi->nseq       = nseq;
  mi->nseqthresh = nseqthresh;
  mi->ishuffled  = ishuffled;
  mi->abc        = abc;

  ESL_ALLOC(mi->pp,           sizeof(double **) * alen);
  ESL_ALLOC(mi->nseff,        sizeof(int     *) * alen);
  ESL_ALLOC(mi->ngap,         sizeof(int     *) * alen);
  ESL_ALLOC(mi->pm,           sizeof(double  *) * alen);
  for (i = 0; i < alen; i++) {
    ESL_ALLOC(mi->pp[i],      sizeof(double  *) * alen);
    ESL_ALLOC(mi->nseff[i],   sizeof(int      ) * alen);
    ESL_ALLOC(mi->ngap[i],    sizeof(int      ) * alen);
    ESL_ALLOC(mi->pm[i],      sizeof(double   ) * K);
    for (j = 0; j < alen; j++) {
       ESL_ALLOC(mi->pp[i][j], sizeof(double  ) * K2);
    }
  }
   
  mi->COV  = esl_dmatrix_Create(alen, alen);
 
  /* initialize for adding counts */
  for (i = 0; i < alen; i++) {
    esl_vec_DSet(mi->pm[i], K, 0.0); 
 
    for (j = 0; j < alen; j++) {
      mi->nseff[i][j] = 0;
      mi->ngap[i][j]  = 0;
      esl_vec_DSet(mi->pp[i][j], K2, 0.0); 
    }
  }

  /* inititalize to zero the COV matrix */
  Mutual_ReuseCOV(mi, COVNONE, C16);
  
  return mi;
  
 ERROR:
  return NULL;
}

int
Mutual_ReuseCOV(struct mutual_s *mi, COVTYPE mitype, COVCLASS miclass)
{
  int  i, j;
  
  mi->type  = mitype;
  mi->class = miclass;

  for (i = 0; i < mi->alen; i++) 
    for (j = 0; j < mi->alen; j++) 
      mi->COV->mx[i][j]  = 0.0;

  mi->besthreshCOV = -eslINFINITY;
  mi->minCOV       =  eslINFINITY;
  mi->maxCOV       = -eslINFINITY;

  return eslOK;
}


void                
Mutual_Destroy(struct mutual_s *mi)
{
  int i, j;

  if (mi) {
    for (i = 0; i < mi->alen; i++) {
      for (j = 0; j < mi->alen; j++) {
	free(mi->pp[i][j]);
      }
      free(mi->nseff[i]);
      free(mi->ngap[i]);
      free(mi->pp[i]);
      free(mi->pm[i]);
    }
    esl_dmatrix_Destroy(mi->COV);
    free(mi->nseff);
    free(mi->ngap);
    free(mi->pp);
    free(mi->pm);
    free(mi);
  }
}


int 
Mutual_NaivePP(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i, j;
  int     status;

  for (i = 0; i < alen-1; i ++)
    for (j = i+1; j < alen; j ++) {
      status = mutual_naive_ppij(i, j, msa, mi, tol, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
  
  return eslOK;

 ERROR:
  return status;
}


int 
Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  ESL_DMATRIX **CL = NULL;
  ESL_DMATRIX **CR = NULL;
  int64_t       alen = msa->alen;
  int           nnodes;
  int           v;
  int           i, j;
  int           status;

  nnodes = (T->N > 1)? T->N - 1 : 1;
  ESL_ALLOC(CL, sizeof(ESL_DMATRIX *) * nnodes);
  ESL_ALLOC(CR, sizeof(ESL_DMATRIX *) * nnodes);
  for (v = 0; v < nnodes; v++) {
    CL[v] = ratematrix_ConditionalsFromRate(T->ld[v], ribosum->bprsQ, tol, errbuf, verbose);
    if (CL[v] == NULL) goto ERROR;
    CR[v] = ratematrix_ConditionalsFromRate(T->rd[v], ribosum->bprsQ, tol, errbuf, verbose);
    if (CR[v] == NULL) goto ERROR;
  }
 
  for (i = 0; i < alen-1; i ++) 
    for (j = i+1; j < alen; j ++) {
      status = mutual_postorder_ppij(i, j, msa, T, ribosum, mi, CL, CR, tol, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
  
  for (v = 0; v < nnodes; v++) {
    esl_dmatrix_Destroy(CL[v]);
    esl_dmatrix_Destroy(CR[v]);
  }
  free(CL);
  free(CR);
  return eslOK;

 ERROR:
  for (v = 0; v < nnodes; v++) {
    if (CL[v]) esl_dmatrix_Destroy(CL[v]);
    if (CR[v]) esl_dmatrix_Destroy(CR[v]);
  }
  if (CL) free(CL);
  if (CR) free(CR);
  return status;
}


int
Mutual_SignificantPairs_Ranking(HIT **ret_hitlist, struct mutual_s *mi, int *msamap, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP,
				int nbpairs, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx = mi->COV;
  char        *covtype = NULL;
  double      *list_exp = NULL;
  double      *list_sc = NULL;
  HIT         *hitlist = NULL;
  double       pval;
  double       delta = 200.;
  double       inc;
  double       min = mi->minCOV;
  double       max = mi->maxCOV;
  double       sen;
  double       ppv;
  double       F;
  double       expect;
  double       bestF = 0.0;
  double       bestsen;
  double       bestppv;
  double       thresh;
  double       besthresh = mi->besthreshCOV;
  double       maxFPF;
  double       maxFPsen;
  double       maxFPppv;
  double       maxFPthresh;
  double       expectFPF;
  double       expectFPsen;
  double       expectFPppv;
  double       expectFPthresh;
  double       oneFPF;
  double       oneFPsen;
  double       oneFPppv;
  double       oneFPthresh;
  double       expectTF_frac_total;   // fraction of covarying basepairs relative to the total number of basepairs
  double       expectTF_frac_surv;    // fraction of covarying basepairs relative to the basepairs that survive the gapthresh
  int          fp, tf, t, f, neg;
  int          oneFP_tf, oneFP_t, oneFP_f, oneFP_fp;
  int          best_tf, best_t, best_f, best_fp;
  int          maxFP_tf, maxFP_t, maxFP_f, maxFP_fp;
  int          expectFP_tf, expectFP_t, expectFP_f, expectFP_fp;
  int          N;
  int          n = 0;
  int          i, j;
  int          status;

  Mutual_COVTYPEString(&covtype, mi->type, errbuf);

  fprintf(rocfp, "\n# %s ", covtype);  
  if (mi->ishuffled) fprintf(rocfp, "shuffled thresh fp tf found true negatives sen ppv F\n"); 
  else               fprintf(rocfp, "thresh fp tf found true negatives sen ppv F\n"); 
  
  inc = (max - min) / delta;
  N = (int)delta;
  ESL_ALLOC(list_sc,  sizeof(double *) * N);
  ESL_ALLOC(list_exp, sizeof(double *) * N);

  for (thresh = max; thresh > min-inc; thresh -= inc) {

    f = t = tf = 0;
    for (i = 0; i < mi->alen-1; i ++) 
      for (j = i+1; j < mi->alen; j ++) {
	if (mtx->mx[i][j] > thresh)   f  ++;
	if (ct[i+1] == j+1) {         t  ++;
	  if (mtx->mx[i][j] > thresh) tf ++;
	}
      }
    
    fp     = f - tf;
    sen    = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv    = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F      = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;
    expect = (mi->alen > 0)? (double)fp/(double)mi->alen : 0.0;
    
    neg = mi->alen * (mi->alen-1) / 2 - t;
    fprintf(rocfp, "%.5f %d %d %d %d %d %.2f %.2f %.2f\n", thresh, fp, tf, f, t, neg, sen, ppv, F);
     
    if (n < N) {
      list_exp[n] = expect;
      list_sc[n]  = thresh;
      n++;
    }
    
    if (expect <= expectFP) {
      expectFPF      = F;
      expectFPsen    = sen;
      expectFPppv    = ppv;
      expectFP_tf    = tf;
      expectFP_f     = f;
      expectFP_t     = t;
      expectFP_fp    = fp;
      expectFPthresh = thresh;
     }
    
    if (maxFP >= 0) {
      if (fp < 1) {
	oneFPF      = F;
	oneFPsen    = sen;
	oneFPppv    = ppv;
	oneFP_tf    = tf;
	oneFP_f     = f;
	oneFP_t     = t;
	oneFP_fp    = fp;
	oneFPthresh = thresh;
      }
      if (fp <= maxFP) {
	maxFPF      = F;
	maxFPsen    = sen;
	maxFPppv    = ppv;
	maxFP_tf    = tf;
	maxFP_f     = f;
	maxFP_t     = t;
	maxFP_fp    = fp;
	maxFPthresh = thresh;
      }
      if (F > bestF) { 
	bestF     = F; 
	bestsen   = sen;
	bestppv   = ppv;
	best_tf   = tf;
	best_t    = t;
	best_f    = f;
	best_fp   = fp;
	besthresh = thresh; 
      }  
    }
  }
      
  if (maxFP >= 0) {
    if (best_fp < maxFP_fp) {
      if (best_fp == 0 && oneFP_tf > best_tf) {
	printf("# %s before 1FP %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, oneFPthresh, min, max,
	       oneFP_fp, oneFP_tf, oneFP_t, oneFP_f, oneFPsen, oneFPppv, oneFPF);
	Mutual_CreateHitList(&hitlist, oneFPthresh, mi, msamap, ct, N, list_sc, list_exp, verbose, errbuf);
      }
      else {
	printf("# %s optimalF %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, besthresh, min, max,
	       best_fp, best_tf, best_t, best_f, bestsen, bestppv, bestF);
	Mutual_CreateHitList(&hitlist, besthresh, mi, msamap, ct, N, list_sc, list_exp, verbose, errbuf);
      } 
    }
    
    printf("# %s maxFP=%d %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, maxFP, maxFPthresh, min, max,
	   maxFP_fp, maxFP_tf, maxFP_t, maxFP_f, maxFPsen, maxFPppv, maxFPF);
    Mutual_CreateHitList(&hitlist, maxFPthresh, mi, msamap, ct, N, list_sc, list_exp, verbose, errbuf);
  }
  else {
    expectTF_frac_total = (nbpairs    > 0)? 100.*(double)expectFP_tf/(double)nbpairs    : 0.0;
    expectTF_frac_surv  = (expectFP_t > 0)? 100.*(double)expectFP_tf/(double)expectFP_t : 0.0;
    fprintf(sumfp, "%s\t%d\t%d\t%d\t%.2f\t%.2f\t", covtype, expectFP_tf, expectFP_t, nbpairs, expectTF_frac_surv, expectTF_frac_total);
    
    Mutual_FisherExactTest(&pval, expectFP_tf, expectFP_fp, expectFP_t, mi->alen);
    
    printf("# %s expectFP=%f (%d cov nonBPs) pval = %f | cov_BP %d/%d covariations %d | sen %f ppv %f F %f] \n", covtype, expectFP, expectFP_fp, pval,
	   expectFP_tf, expectFP_t, expectFP_f, expectFPsen, expectFPppv, expectFPF);
    Mutual_CreateHitList(&hitlist, expectFPthresh, mi, msamap, ct, N, list_sc, list_exp, verbose, errbuf); 
  }

  

  if (ret_hitlist) *ret_hitlist = hitlist; else free(hitlist);
  free(list_sc);
  free(list_exp);
  if (covtype) free(covtype); 
  return eslOK;

 ERROR:
  if (hitlist)  free(hitlist);
  if (list_sc)  free(list_sc);
  if (list_exp) free(list_exp);
  if (covtype)  free(covtype); 
  return status;
}

int 
Mutual_CreateHitList(HIT **ret_hitlist, double threshsc, struct mutual_s *mi, int *msamap, int *ct, int N, double *list_sc, double *list_exp, int verbose, char *errbuf)
{
  HIT    *hit = NULL;
  int     alloc_nhit = 5;
  int     nhit;
  int     h = 0;
  int     i, j;
  int     ih, jh;
  int     n;
  int     status;
  
  nhit = alloc_nhit;
  ESL_ALLOC(hit, sizeof(HIT) * nhit);

 for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      if (mi->COV->mx[i][j] > threshsc) {

	if (h == nhit - 1) {
	  nhit += alloc_nhit;
	  ESL_REALLOC(hit, sizeof(HIT) * nhit);
	}
	for (n = 0; n < N; n ++) {
	  if (mi->COV->mx[i][j] <= list_sc[n]) hit[h].exp = list_exp[n];
	  else break;
	}
	hit[h].i = i;
	hit[h].j = j;
	hit[h].sc = mi->COV->mx[i][j];
	if (ct[i+1] == j+1) { hit[h].is_bpair = TRUE;  }
	else                { 
	  hit[h].is_bpair = FALSE; 
	  if (ct[i+1] == 0 && ct[j+1] == 0) hit[h].is_compatible = TRUE;
	  else                              hit[h].is_compatible = FALSE;
	} 
 
	h ++;
      }
    }
 nhit = h;

 //sort by sc? should go here
  
  for (h = 0; h < nhit; h ++) {
    ih = hit[h].i;
    jh = hit[h].j;
    if      (hit[h].is_bpair)      { printf("* %d %d\t%.4f\t%.2f\n", msamap[ih]+1, msamap[jh]+1, hit[h].exp, hit[h].sc); }
    else if (hit[h].is_compatible) { printf("~ %d %d\t%.4f\t%.2f\n", msamap[ih]+1, msamap[jh]+1, hit[h].exp, hit[h].sc); }
    else                           { printf("  %d %d\t%.4f\t%.2f\n", msamap[ih]+1, msamap[jh]+1, hit[h].exp, hit[h].sc); } 
    
  }

  if (ret_hitlist) *ret_hitlist = hit; else free(hit);
  return eslOK;

 ERROR:
  if (hit) free(hit);
  return status;
}


int
Mutual_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int *ct, int verbose, char *errbuf)
{
  double       avgi, avgj;
  double       stdi, stdj;
  double       zscorei, zscorej;
  double       zscore;
  int          i, j, k;
  int          ipair;

  for (i = 0; i < mi->alen; i ++) {
    if (ct[i+1] > 0 && ct[i+1] > i+1) {
      ipair = ct[i+1]-1;

      avgi = 0.0;
      stdi = 0.0;
      
      for (j = 0; j < mi->alen; j ++) {
	if (j != ipair) {   
	  if (j != i) {
	    avgi += mi->COV->mx[i][j];
	    stdi += mi->COV->mx[i][j] * mi->COV->mx[i][j];
	  }
	}
	else {
	  avgj = 0.0;
	  stdj = 0.0;
	  for (k = 0; k < mi->alen; k ++) {
	    if (k != j && k != i) {       
	      avgj += mi->COV->mx[j][k];
	      stdj += mi->COV->mx[j][k] * mi->COV->mx[j][k];
	    }
	  }
	  avgj /= (mi->alen-2);
	  stdj /= (mi->alen-2);
	  stdj -= avgj*avgj;
	  stdj = sqrt(stdj);
	}
      }
      avgi /= (mi->alen-2);
      stdi /= (mi->alen-2);
      stdi -= avgi*avgi;
      stdi = sqrt(stdi);

      zscorei = (mi->COV->mx[i][ipair] - avgi) / stdi;
      zscorej = (mi->COV->mx[i][ipair] - avgj) / stdj;
      zscore  = ESL_MIN(zscorej, zscorej);
      printf("[%d][%d] %f | %f | %f %f | %f %f\n", msamap[i], msamap[ipair], zscore, mi->COV->mx[i][ipair], avgi, stdi, avgj, stdj);
    }
  }  
  return eslOK;
}

int
Mutual_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen)
{
  double pval = 0.0;
  double add, add_fixed;
  double factorial_total, factorial_BP, factorial_NBP;
  double factorial_cnbp,  factorial_ncnbp;
  double factorial_cbp,   factorial_ncbp;
  double factorial_cov,   factorial_non;
  double tol = 1e-3;
  int    NBP;
  int    total;
  int    cov;
  int    non;
  int    cbp, ncbp;
  int    cnbp, ncnbp;
 
  NBP   = alen * (alen - 1) / 2;
  total = BP    + NBP;
  cov   = cBP   + cNBP;
  non   = total - cov;

  esl_stats_LogGamma(total+1, &factorial_total);
  esl_stats_LogGamma(BP+1,    &factorial_BP);
  esl_stats_LogGamma(NBP+1,   &factorial_NBP);
  esl_stats_LogGamma(cov+1,   &factorial_cov);
  esl_stats_LogGamma(non+1,   &factorial_non);

  add_fixed  = factorial_BP + factorial_NBP + factorial_cov + factorial_non - factorial_total;

  for (cnbp = cNBP; cnbp <= NBP; cnbp ++) {
    cbp   = cov - cnbp; if (cbp < 0) break;
    ncbp  = BP  - cbp;
    ncnbp = non - ncbp;
    
    esl_stats_LogGamma(cbp+1,   &factorial_cbp);
    esl_stats_LogGamma(ncbp+1,  &factorial_ncbp);
    esl_stats_LogGamma(cnbp+1,  &factorial_cnbp);
    esl_stats_LogGamma(ncnbp+1, &factorial_ncnbp);

    add = add_fixed - factorial_cbp - factorial_ncbp - factorial_cnbp - factorial_ncnbp;

    pval += exp(add);
  }
  if (pval > 1.0 && pval < 1.0 + tol) pval = 1.0;

  *ret_pval = pval;

  return eslOK;
}


int
Mutual_FormatR2R(ESL_MSA *msa, HIT *hitlist, int verbose, char *errbuf)
{

  return eslOK;
}


/*---------------- internal functions --------------------- */

static int
is_wc(int x, int y) 
{
  if (x+y == 3 || x+y == 5) return TRUE;
  
  return FALSE;
}

static int    
mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  double *pp = mi->pp[i][j];
  int          K = mi->abc->K;
  int          K2 = K*K;
  int          s;
  int          resi, resj;
  int          x, y;
  
  esl_vec_DSet(pp, K2, 1.0/(double)K2); // laplace prior
  mi->nseff[i][j] = 1;
  
  for (s = 0; s < msa->nseq; s ++) {
    resi = msa->ax[s][i+1];
    resj = msa->ax[s][j+1];
    
    if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) { mi->nseff[i][j] ++; pp[IDX(resi,resj,K)] += 1.0; }
    else if (esl_abc_XIsCanonical(msa->abc, resi)) { mi->nseff[i][j] ++; mi->ngap[i][j] ++; for (y = 0; y < K; y ++) pp[IDX(resi,y,   K)] += 1./(double)K; }
    else if (esl_abc_XIsCanonical(msa->abc, resj)) { mi->nseff[i][j] ++; mi->ngap[i][j] ++; for (x = 0; x < K; x ++) pp[IDX(x,   resj,K)] += 1./(double)K; }
#if 0
    else { 
      mi->nseff[i][j] ++; 
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) 
	  pp[IDX(x,y,K)] += 1./(double)(K*K);
    }
#endif
  }

  /* the probabilities */
  esl_vec_DNorm(pp, K2);                // normalize

  /* symmetrize */
  esl_vec_DCopy(pp, K2, mi->pp[j][i]);
  mi->nseff[j][i] = mi->nseff[i][j];

  return eslOK;
}



int 
mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, ESL_DMATRIX **CL, ESL_DMATRIX **CR, 
		      double tol, int verbose, char *errbuf)
{
  ESL_STACK     *vs = NULL;   /* node index stack */
  ESL_DMATRIX  **lk = NULL;
  ESL_DMATRIX   *lkl, *lkr;
  ESL_DMATRIX   *cl, *cr;
  double         sc;
  int            dim;
  int            K = mi->abc->K;
  int            K2 = K*K;
  int            nnodes;
  int            v;
  int            which;
  int            idx;
  int            resi, resj;
  int            x, y;
  int            xl, yl;
  int            xr, yr;
  int            status;
  
  p7_FLogsumInit();    

  /* allocate the single and pair probs for theinternal nodes */
  nnodes = (T->N > 1)? T->N-1 : T->N;
  dim    = nnodes + T->N;
  ESL_ALLOC(lk, sizeof(ESL_DMATRIX *) * dim);
  for (v = 0; v < dim; v ++)  lk[v] = NULL;
 
  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL)   { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      //printf("%d %d lk-pop v %d | l %d r %d\n", i,  j, v, T->left[v], T->right[v]);

      which = (T->left[v] <= 0)? -T->left[v] : T->left[v];
      idx   = (T->left[v] <= 0)? nnodes + which : which;
      if (T->left[v] <= 0) {
	lk[idx] = esl_dmatrix_Create(K, K);
	esl_dmatrix_Set(lk[idx], -eslINFINITY);
	resi = msa->ax[which][i+1];
	resj = msa->ax[which][j+1];
#if 1
	if (i==5&&j==118) {
	  printf("v=%d parent %d lk l time %f %d which %d i %d %d (%d) j %d %d (%d)\n", 
		 v, T->parent[v], T->ld[v], idx, which, 
		 i, resi, esl_abc_XIsCanonical(msa->abc, resi), 
		 j, resj, esl_abc_XIsCanonical(msa->abc, resj));
	}
#endif

	if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) {
	  lk[idx]->mx[resi][resj] = 0.0;
	}
	else if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  for (y = 0; y < K; y ++) lk[idx]->mx[resi][y] = -log((double)K);
	}
	else if (esl_abc_XIsCanonical(msa->abc, resj)) {
	  for (x = 0; x < K; x ++) lk[idx]->mx[x][resj] = -log((double)K);
	}
	else {
	  esl_dmatrix_Set(lk[idx], -log((double)K2));
	}
      }
      lkl = lk[idx];
 
      which = (T->right[v] <= 0)? -T->right[v] : T->right[v];
      idx   = (T->right[v] <= 0)? nnodes + which : which;
      if (T->right[v] <= 0) {
	lk[idx] = esl_dmatrix_Create(K, K);
	esl_dmatrix_Set(lk[idx], -eslINFINITY); 
	resi = msa->ax[which][i+1];
	resj = msa->ax[which][j+1];
#if 1
	if (i==5&&j==118) {
	  printf("v=%d parent %d lk r time %f %d which %d i %d %d (%d) j %d %d (%d)\n", 
		 v, T->parent[v], T->rd[v], idx, which, 
		 i, resi, esl_abc_XIsCanonical(msa->abc, resi), 
		 j, resj, esl_abc_XIsCanonical(msa->abc, resj));
	}
#endif

	if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) {
	  lk[idx]->mx[resi][resj] = 0.0;
	}
	else if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  for (y = 0; y < K; y ++) lk[idx]->mx[resi][y] = -log((double)K);
	}
	else if (esl_abc_XIsCanonical(msa->abc, resj)) {
	  for (x = 0; x < K; x ++) lk[idx]->mx[x][resj] = -log((double)K);
	}
	else {
	  esl_dmatrix_Set(lk[idx], -log((double)K2));
	}
      }
      lkr = lk[idx];

      if (lkl != NULL && lkr != NULL) { /* ready to go: calculate ps and lk at the parent node */
	cl = CL[v];
	cr = CR[v];
     
	lk[v] = esl_dmatrix_Create(K, K);
	for (x = 0; x < K; x ++) 
	  for (y = 0; y < K; y ++) {
	    sc = -eslINFINITY;
	    for (xl = 0; xl < K; xl ++) 
	      for (yl = 0; yl < K; yl ++) 
		for (xr = 0; xr < K; xr ++) 
		  for (yr = 0; yr < K; yr ++) 
		    sc = p7_FLogsum(sc, lkl->mx[xl][yl] + log(cl->mx[IDX(x,y,K)][IDX(xl,yl,K)]) + lkr->mx[xr][yr] + log(cr->mx[IDX(x,y,K)][IDX(xr,yr,K)]));
	    lk[v]->mx[x][y] = sc; 
	  }
	
#if 0
	double sum;
	sum = -eslINFINITY;
	for (x = 0; x < K; x ++) 
	  for (y = 0; y < K; y ++) 
	    sum = p7_FLogsum(sum, lk[v]->mx[x][y]);
	for (x = 0; x < K; x ++) 
	  for (y = 0; y < K; y ++) 
	    lk[v]->mx[x][y] -= sum;
#endif

#if 1
	if (i==5&&j==118) {
	  printf("l %d r %d v %d\n", T->left[v], T->right[v], v);
	  esl_dmatrix_Dump(stdout, cl,    NULL,   NULL);
	  esl_dmatrix_Dump(stdout, cr,    NULL,   NULL);
	  esl_dmatrix_Dump(stdout, lkl,   "ACGU", "ACGU");
	  esl_dmatrix_Dump(stdout, lkr,   "ACGU", "ACGU");
	  esl_dmatrix_Dump(stdout, lk[v], "ACGU", "ACGU");
	}
#endif

	/* push parent into stack unless already at the root */
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslFAIL; goto ERROR; }; 
      }
      else if (lkl == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])  != eslOK) { status = eslFAIL; goto ERROR; };
      }
      else if (lkr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v]) != eslOK) { status = eslFAIL; goto ERROR; }
      }
    }
  if (v != 0) ESL_XFAIL(eslFAIL, errbuf, "lk did not transverse tree to the root");
  
  for (x = 0; x < K; x ++) 
    for (y = 0; y < K; y ++) {
      //mi->pp[i][j][IDX(x,y,K)] = exp(lk[v]->mx[x][y]) * ribosum->bprsM[IDX(x,y,K)];
      mi->pp[i][j][IDX(x,y,K)] = exp(lk[v]->mx[x][y]);
    }
  esl_vec_DNorm(mi->pp[i][j], K*K);

  for (v = 0; v < dim; v ++) esl_dmatrix_Destroy(lk[v]);
  free(lk);
  esl_stack_Destroy(vs);
  return eslOK;

 ERROR:
  if (vs) esl_stack_Destroy(vs);
  for (v = 0; v < dim; v ++) if (lk[v]) free(lk[v]);
  if (lk) free(lk);
  return status;
}

