/* covariation.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_histogram.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covariation.h"
#include "cococyk.h"
#include "covgrammars.h"
#include "cykcov.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"

static int is_wc(int x, int y);
static int is_stacked_pair(int i, int j, int L, int *ct);
static int number_pairs(int L, int *ct);
static int is_cannonical_pair(char nti, char ntj);
static int mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, ESL_DMATRIX **CL, ESL_DMATRIX **CR, 
				 double tol, int verbose, char *errbuf);
static int cykcov_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop);
static int cov_histogram_plotsurvival(FILE *pipe, ESL_HISTOGRAM *h, int logscale, int style1, int style2);

int                 
cov_Calculate(ESL_MSA **omsa, int *msamap, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
	      RANKLIST *ranklist_null, RANKLIST **ret_ranklist, HITLIST **ret_hitlist,
	      METHOD method, COVTYPE covtype, COVCLASS covclass, int *ct, 
	      FILE *outfp, FILE *rocfp, FILE *sumfp, char *gnuplot, char *dplotfile, char *r2rfile, char *r2rversion, int r2rall, 
	      THRESH *thresh, MODE mode, int nbpairs, double tol, int verbose, char *errbuf)
{
  ESL_MSA   *msa = *omsa;
  RANKLIST  *ranklist = NULL;
  HITLIST   *hitlist = NULL;
  int        status;

  status = cov_Probs(msa, T, ribosum, mi, method, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  switch(covtype) {
   case CHIa: 
     status = cov_CalculateCHI         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case CHIp:
     status = cov_CalculateCHI         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;  
     break;
   case CHI: 
     status = cov_CalculateCHI         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;     
     break;
   case GTa: 
     status = cov_CalculateGT          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC,  TRUE, &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
      break;
   case GTp: 
     status = cov_CalculateGT          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC,  TRUE, &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case GT: 
     status = cov_CalculateGT          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,       TRUE, &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case MIa: 
     status = cov_CalculateMI          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MIp: 
     status = cov_CalculateMI          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MI: 
     status = cov_CalculateMI          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case MIra: 
     status = cov_CalculateMIr         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(         mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC,  TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MIrp:
     status = cov_CalculateMIr         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;  
     break;
   case MIr: 
     status = cov_CalculateMIr         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE, &ranklist, &hitlist,  ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case MIga: 
     status = cov_CalculateMIg         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case MIgp:
     status = cov_CalculateMIg         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;  
     break;
   case MIg: 
     status = cov_CalculateMIg          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,     TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   case OMESa: 
     status = cov_CalculateOMES        (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case OMESp: 
     status = cov_CalculateOMES        (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     break;
   case OMES: 
     status = cov_CalculateOMES        (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
      break;
   case COVALL: 
     status = cov_CalculateCHI         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = cov_CalculateCHI         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     
     status = cov_CalculateOMES        (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = cov_CalculateOMES        (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     
     status = cov_CalculateGT          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = cov_CalculateGT          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     
     status = cov_CalculateMI          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = cov_CalculateMI          (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     
     status = cov_CalculateMIr         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = cov_CalculateMIr         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     
     status = cov_CalculateMIg         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      TRUE,  NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, APC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR; 
     status = cov_CalculateMIg         (covclass, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs,      FALSE, NULL,      NULL,     ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     status = cov_CalculateCOVCorrected(          mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, ASC, TRUE,  &ranklist, &hitlist, ranklist_null, tol, verbose, errbuf);
     if (status != eslOK) goto ERROR;
     break;
   default:
     ESL_XFAIL(eslFAIL, errbuf, "wrong covariation type\n");
     break;
   }
   fprintf(sumfp, "\n");   
      
   if (mode == GIVSS || mode == CYKSS) { // do the plots only GIVSS or CYKSS
     status = cov_DotPlot(gnuplot, dplotfile, msa, ct, mi, msamap, hitlist, TRUE, verbose, errbuf);
     if  (status != eslOK) goto ERROR;
     status = cov_DotPlot(gnuplot, dplotfile, msa, ct, mi, msamap, hitlist, FALSE, verbose, errbuf);
     if  (status != eslOK) goto ERROR;
     
     status = cov_R2R(r2rfile, r2rversion, r2rall, &msa, ct, msamap, hitlist, TRUE, TRUE, verbose, errbuf);
     if  (status != eslOK) goto ERROR;
   }

   *omsa = msa;
   if (ret_ranklist) *ret_ranklist = ranklist; else if (ranklist) cov_FreeRankList(ranklist);
   if (ret_hitlist)  *ret_hitlist = hitlist;   else if (hitlist)  cov_FreeHitList(hitlist);
   return eslOK;
   
 ERROR:
   if (ranklist) cov_FreeRankList(ranklist);
   if (hitlist)  cov_FreeHitList(hitlist);
  return status;
}

int                 
cov_Probs(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf)
{
  int i, j;
  int x, y;
  int K = msa->abc->K;
  int status;

  switch(method) {
  case NAIVE:
    status = cov_NaivePP(msa, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    break;
  case PHYLO:
    status = cov_PostOrderPP(msa, T, ribosum, mi, tol, verbose, errbuf);
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

  status = cov_ValidateProbs(mi, tol, verbose, errbuf);
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
cov_ValidateProbs(struct mutual_s *mi, double tol, int verbose, char *errbuf)
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
cov_CalculateCHI(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
		    int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateCHI_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateCHI_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = cov_CalculateCHI_C2 (mi, verbose, errbuf);
    else
      status = cov_CalculateCHI_C16(mi, verbose, errbuf);
    break;
}
  
  if (verbose) {
    printf("CHI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("CHI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }

  if (analyze) 
    status = cov_SignificantPairs_Ranking(ranklist_null, ret_ranklist, ret_hitlist, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}


int                 
cov_CalculateCHI_C16(struct mutual_s *mi, int verbose, char *errbuf)
{
  double chi;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, CHI, C16);
  
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
cov_CalculateCHI_C2(struct mutual_s *mi, int verbose, char *errbuf)
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
  
  cov_ReuseCOV(mi, CHI, C2);
  
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
cov_CalculateOMES(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
		     int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateOMES_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateOMES_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = cov_CalculateOMES_C2 (mi, verbose, errbuf);
    else
      status = cov_CalculateOMES_C16(mi, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("OMES[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("OMES[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) 
    status = cov_SignificantPairs_Ranking(ranklist_null, ret_ranklist, ret_hitlist, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;
  
 ERROR:
  return status;
}

int                 
cov_CalculateOMES_C16(struct mutual_s *mi, int verbose, char *errbuf)
{
  double omes;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, OMES, C16);
  
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
cov_CalculateOMES_C2(struct mutual_s *mi, int verbose, char *errbuf)
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
  
  cov_ReuseCOV(mi, OMES, C2);
  
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
cov_CalculateGT(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
		   int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateGT_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateGT_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh) {
      status = cov_CalculateGT_C2 (mi, verbose, errbuf);
    }
    else {
      status = cov_CalculateGT_C16(mi, verbose, errbuf);
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
  
  if (analyze) 
    status = cov_SignificantPairs_Ranking(ranklist_null, ret_ranklist, ret_hitlist, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}


int                 
cov_CalculateGT_C16(struct mutual_s *mi, int verbose, char *errbuf)
{
  double gt;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, GT, C16);
  
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
cov_CalculateGT_C2(struct mutual_s *mi, int verbose, char *errbuf)
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
  
  cov_ReuseCOV(mi, GT, C2);
  
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
cov_CalculateMI(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
		   int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateMI_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateMI_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = cov_CalculateMI_C2 (mi, verbose, errbuf);
    else
      status = cov_CalculateMI_C16(mi, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("MI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) 
    status = cov_SignificantPairs_Ranking(ranklist_null, ret_ranklist, ret_hitlist, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}

int                 
cov_CalculateMI_C16(struct mutual_s *mi, int verbose, char *errbuf)
{
  double mutinf;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, MI, C16);
  
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
cov_CalculateMI_C2(struct mutual_s *mi, int verbose, char *errbuf)
{
  double mutinf;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, MI, C2);
  
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
cov_CalculateMIr(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
		    int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateMIr_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateMIr_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = cov_CalculateMIr_C2 (mi, verbose, errbuf);
    else
      status = cov_CalculateMIr_C16(mi, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("MIr[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MIr[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) 
    status = cov_SignificantPairs_Ranking(ranklist_null, ret_ranklist, ret_hitlist, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}

int                 
cov_CalculateMIr_C16(struct mutual_s *mi, int verbose, char *errbuf)
{
  double mutinf, HH;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, MIr, C16);
  
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
cov_CalculateMIr_C2(struct mutual_s *mi, int verbose, char *errbuf)
{
  double mutinf, HH;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, MIr, C2);
  
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
cov_CalculateMIg(COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
		    int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, double tol, int verbose, char *errbuf)
{
  int i,j;
  int status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateMIg_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateMIg_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh)
      status = cov_CalculateMIg_C2 (mi, verbose, errbuf);
    else
      status = cov_CalculateMIg_C16(mi, verbose, errbuf);
    break;
  }
  
  if (verbose) {
    printf("MIg[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MIg[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) 
    status = cov_SignificantPairs_Ranking(ranklist_null, ret_ranklist, ret_hitlist, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return status;
  
 ERROR:
  return status;
}

int                 
cov_CalculateMIg_C16(struct mutual_s *mi, int verbose, char *errbuf)
{
  double mutinf;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, MIg, C16);
  
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
cov_CalculateMIg_C2(struct mutual_s *mi, int verbose, char *errbuf)
{
  double mutinf;
  double pij_wc, pij_nwc;
  double qij_wc, qij_nwc;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  cov_ReuseCOV(mi, MIg, C2);
  
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
cov_CalculateCOVCorrected(struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
			     CORRTYPE corrtype, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, double tol, int verbose, char *errbuf)
{
  char        *covtype = NULL;
  ESL_DMATRIX *COV  = NULL;
  double      *COVx = NULL;
  double       COVavg = 0.0;
  int          i, j;
  int          status = eslOK;
  
  cov_COVTYPEString(&covtype, mi->type, errbuf);

  switch(corrtype) {
  case APC: esl_sprintf(&covtype, "%sp", covtype); break;
  case ASC: esl_sprintf(&covtype, "%sa", covtype); break;
  default:  
    ESL_XFAIL(eslFAIL, errbuf, "wrong correction type\n");
    break;
  }
  COV = esl_dmatrix_Clone(mi->COV);
  
  cov_String2COVTYPE(covtype, &mi->type, errbuf);
  cov_ReuseCOV(mi, mi->type, mi->class);  
 
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

  if (analyze) 
    status = cov_SignificantPairs_Ranking(ranklist_null, ret_ranklist, ret_hitlist, mi, msamap, ct, outfp, rocfp, sumfp, thresh, mode, nbpairs, verbose, errbuf);
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
cov_THRESHTYPEString(char **ret_threshtype, THRESHTYPE type, char *errbuf)
{
  int status;

  switch(type) {
  case covNBP:   esl_sprintf(ret_threshtype, "covNBP");  break;
  case covNBPu:  esl_sprintf(ret_threshtype, "covNBPu"); break;
  case covNBPf:  esl_sprintf(ret_threshtype, "covNBPf"); break;
  case covRBP:   esl_sprintf(ret_threshtype, "covRBP");  break;
  case covRBPu:  esl_sprintf(ret_threshtype, "covRBPu"); break;
  case covRBPf:  esl_sprintf(ret_threshtype, "covRBPu"); break;
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong THRESHTYPE");
  }

  return eslOK;
  
 ERROR:
  return status;
}
int 
cov_COVTYPEString(char **ret_covtype, COVTYPE type, char *errbuf)
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
cov_String2COVTYPE(char *covtype, COVTYPE *ret_type, char *errbuf)
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
cov_Create(int64_t alen, int64_t nseq, int ishuffled, int nseqthresh, ESL_ALPHABET *abc)
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
  cov_ReuseCOV(mi, COVNONE, C16);
  
  return mi;
  
 ERROR:
  return NULL;
}

int
cov_ReuseCOV(struct mutual_s *mi, COVTYPE mitype, COVCLASS miclass)
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
cov_Destroy(struct mutual_s *mi)
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
cov_NaivePP(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
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
cov_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
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
cov_SignificantPairs_Ranking(RANKLIST *ranklist_null, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, struct mutual_s *mi, int *msamap, int *ct, 
			     FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx = mi->COV;
  char        *threshtype = NULL;
  char        *covtype = NULL;
  RANKLIST    *ranklist = NULL;
  HITLIST     *hitlist = NULL;
  double       sen;
  double       ppv;
  double       F;
  double       cvBP, cvNBP, cvNBPu, cvNBPf;
  double       cvRBP, cvRBPu, cvRBPf;
  double       cov;
  double       val;
  double       bmax, bmin;
  int          fp, tf, t, f, neg;
  int          i, j;
  int          x, nullx;
  int          status;

  cov_COVTYPEString(&covtype, mi->type, errbuf);
  cov_THRESHTYPEString(&threshtype, thresh->type, NULL);

  if (rocfp) {
    fprintf(rocfp, "\n# %s ", covtype);  
    if (mi->ishuffled) fprintf(rocfp, "shuffled thresh fp tf found true negatives sen ppv F\n"); 
    else               fprintf(rocfp, "thresh fp tf found true negatives sen ppv F\n"); 
  }

  bmax = mi->maxCOV+W;
  bmin = mi->minCOV-W;
  ranklist = cov_CreateRankList(bmax, bmin, W);
  ranklist->scmax = bmax-W;
  ranklist->scmin = bmin+W;

  for (i = 0; i < mi->alen-1; i ++) 
    for (j = i+1; j < mi->alen; j ++) {
      /* add to the histogram if not a real basepair */
      if (mode == GIVSS && ct[i+1] != j+1) 
	esl_histogram_Add(ranklist->h, mtx->mx[i][j]);
      else if (mode == RANSS) 
	esl_histogram_Add(ranklist->h, mtx->mx[i][j]);
    }
  
  for (x = ranklist->nb-1; x >= 0; x --) {
    cov = cov_ranklist_Bin2Mid(ranklist,x);
     
    f = t = tf = 0;
    for (i = 0; i < mi->alen-1; i ++) 
      for (j = i+1; j < mi->alen; j ++) {
	if (mtx->mx[i][j] > cov)   f  ++;
	if (ct[i+1] == j+1) {      t  ++;
	  if (mtx->mx[i][j] > cov) tf ++;
	}	
     }
    
    fp  = f - tf;
    sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;   
    neg = mi->alen * (mi->alen-1) / 2 - t;
    if (rocfp) fprintf(rocfp, "%.5f %d %d %d %d %d %.2f %.2f %.2f\n", cov, fp, tf, f, t, neg, sen, ppv, F);
    
    cvBP   = (double)tf;
    cvNBP  = (double)fp;
    cvNBPu = (mi->alen > 0)? cvNBP/(double)mi->alen : 0.0;
    cvNBPf = (neg > 0)?      cvNBP/(double)neg      : 0.0;
    ranklist->covBP[x]  = cvBP;
    ranklist->covNBP[x] = cvNBP;
    
    if (mode == GIVSS) {
      if (ranklist_null) {
	if      (cov > ranklist_null->bmax) cvRBP = ranklist_null->covBP[ranklist_null->nb-1];
	else if (cov < ranklist_null->bmin) cvRBP = ranklist_null->covBP[0];
	else {
	  cov_ranklist_Score2Bin(ranklist_null, cov, &nullx);
	  cvRBP = ranklist_null->covBP[nullx];
	}
      }
      cvRBPu = (mi->alen > 0)? cvRBP/(double)mi->alen : 0.0;
      cvRBPf = (neg > 0)?      cvRBP/(double)neg      : 0.0;
      
      switch(thresh->type) {
      case covNBP:  val = cvNBP;  break;
      case covNBPu: val = cvNBPu; break;
      case covNBPf: val = cvNBPf; break;
      case covRBP:  val = cvRBP;  break;
      case covRBPu: val = cvRBPu; break;
      case covRBPf: val = cvRBPf; break;
      }
      if (val <= thresh->val) { ranklist->scthresh = cov + ranklist->w; thresh->sc = ranklist->scthresh; }
    }
    if (mode == CYKSS) {
      if (cov < thresh->sc) { ranklist->scthresh = cov + ranklist->w; break; }
    }
  }
  
  if (mode == GIVSS || mode == CYKSS) {
    status = cov_CreateHitList(outfp, &hitlist, thresh, mi, msamap, ct, ranklist, ranklist_null, covtype, threshtype, verbose, errbuf);
    if (status != eslOK) goto ERROR;
  }    
  
  if (ret_ranklist) *ret_ranklist = ranklist; else if (ranklist) cov_FreeRankList(ranklist);
  if (ret_hitlist)  *ret_hitlist  = hitlist;  else if (hitlist)  cov_FreeHitList(hitlist);
  
  if (threshtype) free(threshtype); 
  if (covtype)    free(covtype); 
  return eslOK;
  
 ERROR:
  if (ranklist)   cov_FreeRankList(ranklist);
  if (hitlist)    cov_FreeHitList(hitlist);
  if (threshtype) free(threshtype); 
  if (covtype)    free(covtype); 
  return status;
}

RANKLIST *
cov_CreateRankList(double bmax, double bmin, double w)
{
  RANKLIST *ranklist = NULL;
  int       status;
  
  ESL_ALLOC(ranklist, sizeof(RANKLIST));

  ranklist->bmax  = bmax;
  ranklist->bmin  = bmin;
  ranklist->scmax = bmax-w;
  ranklist->scmin = bmin+w;
  ranklist->w     = w;
  ranklist->nb    = (int)((bmax-bmin)/w);
  ranklist->h     = esl_histogram_Create(bmin, bmax, w);
 
  ESL_ALLOC(ranklist->covBP,  sizeof(double) * ranklist->nb);
  ESL_ALLOC(ranklist->covNBP, sizeof(double) * ranklist->nb);
 
  esl_vec_DSet(ranklist->covBP,  ranklist->nb, 0.);
  esl_vec_DSet(ranklist->covNBP, ranklist->nb, 0.);

  ranklist->scthresh = ranklist->scmax;
  return ranklist;

 ERROR:
  return NULL;
}

int
cov_GrowRankList(RANKLIST **oranklist, double bmax, double bmin)
{
  RANKLIST *ranklist = *oranklist;
  RANKLIST *new = NULL;
  int       x, newx;
  int       status;

  /* for the low bound bmin, we stay at the highest value, not need
   * to take them all */
  new = cov_CreateRankList(ESL_MAX(bmax, ranklist->bmax), ESL_MAX(bmin, ranklist->bmin), ranklist->w);
  if (new == NULL) goto ERROR;
  new->h->bmin = new->bmin;
  new->h->bmax = new->bmax;
  new->h->xmin = ranklist->h->xmin;
  new->h->xmax = ranklist->h->xmax;
  new->h->imin = ranklist->h->imin;
  new->h->imax = ranklist->h->imax;
  new->h->Nc   = ranklist->h->Nc;
  new->h->No   = ranklist->h->No;

  for (x = 0; x < ranklist->nb; x ++) {
    cov_ranklist_Score2Bin(new, cov_ranklist_Bin2Mid(ranklist,x), &newx);
    if (newx >= 0 && newx < new->nb) {
      new->covBP[newx]   = ranklist->covBP[x];
      new->covNBP[newx]  = ranklist->covNBP[x];
      new->h->obs[newx]  = ranklist->h->obs[x];
    }
  }
    
  cov_FreeRankList(*oranklist);
  *oranklist = new;

  return eslOK;

 ERROR:
  if (new) cov_FreeRankList(new);
  return status;
}

int              
cov_DumpRankList(FILE *fp, RANKLIST *ranklist)
{
  int x;

  printf("bmin %f bmax %f covmin %f covmax %f\n", ranklist->bmin, ranklist->bmax, ranklist->scmin, ranklist->scmax);
  for (x = ranklist->nb-1; x >= 0; x --) 
    printf("cov %f covBP %f covNBP %f\n", (double)x*ranklist->w + ranklist->bmin, ranklist->covBP[x],ranklist->covNBP[x]); 
  
  return eslOK;
}

int 
cov_CreateHitList(FILE *outfp, HITLIST **ret_hitlist, THRESH *thresh, struct mutual_s *mi, int *msamap, int *ct, RANKLIST *ranklist, RANKLIST *ranklist_null, 
		  char *covtype, char *threshtype, int verbose, char *errbuf)
{
  HITLIST *hitlist = NULL;
  double   sen, ppv, F;
  int      alloc_nhit = 5;
  int      tf = 0;
  int      f  = 0;
  int      t, fp;
  int      BP, NBP;
  int      nhit;
  int      h = 0;
  int      i, j;
  int      ih, jh;
  int      x;
  int      status;
  
  ESL_ALLOC(hitlist, sizeof(HITLIST));

  nhit = alloc_nhit;
  ESL_ALLOC(hitlist->hit, sizeof(HIT) * nhit);
  hitlist->covthresh = ranklist->scthresh;

  BP  = number_pairs(mi->alen, ct);
  NBP = mi->alen * (mi->alen-1) / 2 - BP;

 for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      if (mi->COV->mx[i][j] > ranklist->scthresh) {

	if (h == nhit - 1) {
	  nhit += alloc_nhit;
	  ESL_REALLOC(hitlist->hit, sizeof(HIT) * nhit);
	}
	/* initialize */
	 hitlist->hit[h].sc            = -eslINFINITY;
	 hitlist->hit[h].Eval          = 1000.;
	 hitlist->hit[h].covNBP        = 0.0;
	 hitlist->hit[h].covNBPu       = 0.0;
	 hitlist->hit[h].covNBPf       = 0.0;
	 hitlist->hit[h].covRBP        = 0.0;
	 hitlist->hit[h].covRBPu       = 0.0;
	 hitlist->hit[h].covRBPf       = 0.0;
	 hitlist->hit[h].is_bpair      = FALSE;
	 hitlist->hit[h].is_compatible = FALSE;
	 
	 for (x = ranklist->nb-1; x >= 0; x --) {
	   if (mi->COV->mx[i][j] <= ranklist->bmin+(double)x*ranklist->w) {
	     hitlist->hit[h].covNBP  = ranklist->covNBP[x];
	     hitlist->hit[h].covNBPu = (mi->alen > 0)? ranklist->covNBP[x]/(double)mi->alen:0.;
	     hitlist->hit[h].covNBPf = ranklist->covNBP[x]/(double)NBP;
	   }
	   else break;
	 }
	 if (ranklist_null) {
	   for (x = ranklist_null->nb-1; x >= 0; x --) {
	     if (mi->COV->mx[i][j] <= ranklist_null->bmin+(double)x*ranklist_null->w) {
	       hitlist->hit[h].covRBP  = ranklist_null->covBP[x];
	       hitlist->hit[h].covRBPu = (mi->alen > 0)? ranklist_null->covBP[x]/(double)mi->alen:0.;
	       hitlist->hit[h].covRBPf = ranklist_null->covNBP[x]/(double)BP;
	     }
	     else break;
	   }
	 }
	 
	 hitlist->hit[h].i = i;
	 hitlist->hit[h].j = j;
	 hitlist->hit[h].sc = mi->COV->mx[i][j];
	 if (ct[i+1] == j+1) { hitlist->hit[h].is_bpair = TRUE;  }
	 else                { 
	   hitlist->hit[h].is_bpair = FALSE; 
	   if (ct[i+1] == 0 && ct[j+1] == 0) hitlist->hit[h].is_compatible = TRUE;
	 } 	 
	 h ++;
      }
    }
 nhit = h;
 hitlist->nhit = nhit;

 t = BP;
 for (h = 0; h < nhit; h ++) {
   f ++;
   if (hitlist->hit[h].is_bpair) tf ++;
 } 
 fp = f - tf;
 sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
 ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
 F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;   
 
 if (outfp) fprintf(outfp, "# %s thresh %s %f cov=%f [%f,%f] [%d | %d %d %d | %f %f %f] \n", 
		    covtype, threshtype, thresh->val, ranklist->scthresh, ranklist->scmin, ranklist->scmax, fp, tf, t, f, sen, ppv, F);
 for (h = 0; h < nhit; h ++) {
   ih = hitlist->hit[h].i;
   jh = hitlist->hit[h].j;
   if (outfp) {
     if (hitlist->hit[h].is_bpair)      { 
       fprintf(outfp, "*\t%10d\t%10d\t%.2f\t%.4f\n", 
	       msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc, hitlist->hit[h].Eval); 
     }
     else if (hitlist->hit[h].is_compatible) { 
       fprintf(outfp, "~\t%10d\t%10d\t%.2f\t%.4f\n", 
	       msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc, hitlist->hit[h].Eval); 
    }
     else { 
       fprintf(outfp, " \t%10d\t%10d\t%.2f\t%.4f\n",
	       msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc, hitlist->hit[h].Eval); 
    } 
   }
 }
 
 if (ret_hitlist) *ret_hitlist = hitlist; else cov_FreeHitList(hitlist);
 return eslOK;
 
 ERROR:
 if (hitlist) cov_FreeHitList(hitlist);
 return status;
}

void
cov_FreeRankList(RANKLIST *ranklist)
{
  if (ranklist == NULL) return;

  if (ranklist->h)      esl_histogram_Destroy(ranklist->h);
  if (ranklist->covBP)  free(ranklist->covBP);
  if (ranklist->covNBP) free(ranklist->covNBP);
  free(ranklist);
}

void
cov_FreeHitList(HITLIST *hitlist)
{
  if (hitlist == NULL) return;

  if (hitlist->hit) free(hitlist->hit);
  free(hitlist);
}

int
cov_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int *ct, int verbose, char *errbuf)
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
cov_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen)
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
cov_CYKCOVCT(FILE *outfp, char *gnuplot, char *dplotfile, char *R2Rcykfile, char *R2Rversion, int R2Rall,  ESL_RANDOMNESS *r, ESL_MSA **omsa, struct mutual_s *mi, 
	     int *msamap, int minloop, enum grammar_e G, THRESH *thresh, double covthresh, int nbpairs, char *errbuf, int verbose)
{
  ESL_MSA *msa = *omsa;
  HITLIST *hitlist = NULL;
  int     *cykct = NULL;
  char    *ss = NULL;	
  SCVAL    sc;
  int      status;
            
  /* calculate the cykcov ct vector */
  status = CYKCOV(r, mi, &cykct, &sc, minloop, covthresh, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  if (verbose) printf("cykcov score = %f\n", sc);

  /* impose the ct on the msa GC line 'cons_ss' */
  ESL_ALLOC(ss, sizeof(char) * (msa->alen+1));
  esl_ct2simplewuss(cykct, msa->alen, ss);
  /* replace the 'SS_cons' GC line with the new ss */
  esl_sprintf(&(msa->ss_cons), "%s", ss);  
  
  /* redo the hitlist since the ct has now changed */
  status = cov_SignificantPairs_Ranking(NULL, NULL, &hitlist, mi, msamap, cykct, outfp, NULL, NULL, thresh, CYKSS, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  /* R2R */
  status = cov_R2R(NULL, R2Rversion, R2Rall, &msa, cykct, msamap, hitlist, FALSE, FALSE, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  esl_msa_Digitize(mi->abc, msa, errbuf);

  esl_ct2simplewuss(cykct, msa->alen, ss);
  //printf("ss:%s\n", ss);
  
  /* expand the CT with compatible/stacked A:U C:G G:U pairs */
  status = cov_ExpandCT(R2Rcykfile, R2Rall, r, msa, &cykct, minloop, G, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  if (verbose) eslx_msafile_Write(stdout, msa, eslMSAFILE_PFAM);
  
  /* R2Rpdf */
  status = cov_R2Rpdf(R2Rcykfile, R2Rversion, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  /* R2Rsvg */
  status = cov_R2Rsvg(R2Rcykfile, R2Rversion, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  /* DotPlots (pdf,svg) */
  status = cov_DotPlot(gnuplot, dplotfile, msa, cykct, mi, msamap, hitlist, TRUE, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = cov_DotPlot(gnuplot, dplotfile, msa, cykct, mi, msamap, hitlist, FALSE, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  *omsa = msa;

  cov_FreeHitList(hitlist);
  free(cykct);
  free(ss);
  return eslOK;
  
 ERROR:
  if (hitlist) cov_FreeHitList(hitlist);
  if (cykct)   free(cykct);
  if (ss)      free(ss);
  return status;
}

int 
cov_WriteHistogram(char *gnuplot, char *covhisfile, char *nullcovhisfile, RANKLIST *ranklist, RANKLIST *ranklist_null, int dosvg, char *errbuf)
{
  FILE    *fp = NULL;
  int      status;

  if ((fp = fopen(covhisfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "could not open file %s\n", covhisfile);
  esl_histogram_PlotSurvival(fp, ranklist->h);
  fclose(fp);
  if (ranklist_null) {
    if ((fp = fopen(nullcovhisfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "could not open file %s\n", nullcovhisfile);
    esl_histogram_PlotSurvival(fp, ranklist_null->h);
    fclose(fp);
  }

  status = cov_PlotHistogramSurvival(gnuplot, covhisfile, ranklist, ranklist_null, dosvg, errbuf);
  if (status != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;

}

int 
cov_PlotHistogramSurvival(char *gnuplot, char *covhisfile, RANKLIST *ranklist, RANKLIST *ranklist_null, int dosvg, char *errbuf)
{
  FILE    *pipe;
  char    *filename = NULL;
  char    *outplot = NULL;
  int      status;

  if (!gnuplot) return eslOK;

  esl_FileTail(covhisfile, FALSE, &filename);

  pipe = popen(gnuplot, "w");

  if (dosvg) {
    esl_sprintf(&outplot, "%s.svg", covhisfile);
    fprintf(pipe, "set terminal svg fname 'Verdana' fsize 10 \n");
  }
  else {
    esl_sprintf(&outplot, "%s.ps", covhisfile);
    fprintf(pipe, "set terminal postscript color 14\n");
  }

  fprintf(pipe, "set output '%s'\n", outplot);
  fprintf(pipe, "unset key\n");
  //fprintf(pipe, "set pointsize %f\n", pointsize);
  fprintf(pipe, "set title '%s' \n", filename);

  fprintf(pipe, "set style line 1   lt 1 lc rgb 'grey' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 2   lt 1 lc rgb 'brown' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 3   lt 1 lc rgb 'cyan' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 4   lt 1 lc rgb 'red' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 5   lt 1 lc rgb 'orange' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 6   lt 1 lc rgb 'turquoise' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 7   lt 1 lc rgb 'black' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 8   lt 1 lc rgb 'green' pt 7 lw 2 ps 0.5\n");
  fprintf(pipe, "set style line 9   lt 1 lc rgb 'blue' pt 7 lw 2 ps 0.5\n");

  // survival plot for ranklist and ranklist_null
  fprintf(pipe, "set multiplot\n");
  fprintf(pipe, "set xlabel 'covariation score'\n");
  if (ranklist_null)
    fprintf(pipe, "set xrange [%f:%f]\n", ESL_MIN(ranklist->h->xmin,ranklist_null->h->xmin), ESL_MAX(ranklist->h->xmax,ranklist_null->h->xmax));
  else 
    fprintf(pipe, "set xrange [%f:%f]\n", ranklist->h->xmin, ranklist->h->xmax);
  fprintf(pipe, "set ylabel 'P(x > score)'\n");
  fprintf(pipe, "set yrange [0:1.0]\n");
     status = cov_histogram_plotsurvival  (pipe, ranklist->h,      FALSE, 4, 2);
  if (status != eslOK) goto ERROR;
  if (ranklist_null) {
    status = cov_histogram_plotsurvival(pipe, ranklist_null->h, FALSE, 1, 7);
    if (status != eslOK) goto ERROR;
  }

  // log survival plot for ranklist and ranklist_null
  fprintf(pipe, "set multiplot\n");
  fprintf(pipe, "set xlabel 'covariation score'\n");
  if (ranklist_null)
    fprintf(pipe, "set xrange [%f:%f]\n", ESL_MIN(ranklist->h->xmin,ranklist_null->h->xmin), ESL_MAX(ranklist->h->xmax,ranklist_null->h->xmax));
  else 
    fprintf(pipe, "set xrange [%f:%f]\n", ranklist->h->xmin, ranklist->h->xmax);
  fprintf(pipe, "set ylabel 'logP(x > score)'\n");
  fprintf(pipe, "set yrange [%f:%f]\n", -log(ranklist->h->No), log(0.01));
    status = cov_histogram_plotsurvival  (pipe, ranklist->h,      TRUE, 4, 2);
  if (status != eslOK) goto ERROR;
  if (ranklist_null) {
    status = cov_histogram_plotsurvival(pipe, ranklist_null->h, TRUE, 1, 7);
    if (status != eslOK) goto ERROR;
  }

  pclose(pipe);
  
  free(outplot);
  free(filename);

  return eslOK;

 ERROR:
  return status;
}


int 
cov_CreateNullCov(char *gnuplot, char *nullcovfile, int L, int *ct, RANKLIST *ranklist, RANKLIST *ranklist_null, int dosvg, char *errbuf)
{
  FILE    *fp = NULL;
  double   cov;
  double   covBP, covBP_prv;
  double   covNBP, covNBP_prv;
  double   covRBP, covRBP_prv;
  int      nullx;
  int      BP;
  int      x;
  int      status;

  if (ranklist_null == NULL) return eslOK;

  BP = number_pairs(L, ct); if (BP <= 0.) return eslOK;
  
  if ((fp = fopen(nullcovfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "could not open file %s\n", nullcovfile);

  covRBP_prv = -1.0;
  covBP_prv  = 0.0;
  covNBP_prv = 0.0;
  for (x = ranklist->nb-1; x >= 0; x--) {
    cov = cov_ranklist_Bin2Mid(ranklist,x);
    covBP  = ranklist->covBP[x];
    covNBP = ranklist->covNBP[x];
    
    if      (cov > ranklist_null->bmax) covRBP = 0;
    else if (cov < ranklist_null->bmin) covRBP = ranklist_null->covBP[0];
    else {
      cov_ranklist_Score2Bin(ranklist_null, cov, &nullx);
      covRBP = ranklist_null->covBP[nullx];
    }
    
    if (covRBP > covRBP_prv) {
      fprintf(fp, "%f %f %f %f %f %f\n", cov, covRBP, 100.*covRBP/(double)BP, covBP, 100.0*covBP/(double)BP, covNBP);
    }
    
    covBP_prv  = covBP;
    covNBP_prv = covNBP;
    covRBP_prv = covRBP;
  } 
  fclose(fp);

  if (gnuplot) cov_PlotNullCov(gnuplot, nullcovfile, 2.*(double)BP, 0.2, 0.2*(double)BP, dosvg);
  return eslOK;

 ERROR:
  return status;
}

int 
cov_PlotNullCov(char *gnuplot, char *nullcovfile, double maxBP, double maxcovRBP, double maxcovRBPf, int dosvg)
{
  FILE    *pipe;
  char    *filename = NULL;
  char    *outplot = NULL;
 
  esl_FileTail(nullcovfile, FALSE, &filename);

  pipe = popen(gnuplot, "w");

  if (dosvg) {
    esl_sprintf(&outplot, "%s.svg", nullcovfile);
    fprintf(pipe, "set terminal svg fname 'Verdana' fsize 10 \n");
  }
  else {
    esl_sprintf(&outplot, "%s.ps", nullcovfile);
    fprintf(pipe, "set terminal postscript color 14\n");
  }

  fprintf(pipe, "set output '%s'\n", outplot);
  fprintf(pipe, "unset key\n");
  //fprintf(pipe, "set pointsize %f\n", pointsize);
  fprintf(pipe, "set title '%s' \n", filename);

  fprintf(pipe, "set style line 1   lt 1 lc rgb 'grey' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 2   lt 1 lc rgb 'brown' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 3   lt 1 lc rgb 'cyan' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 4   lt 1 lc rgb 'red' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 5   lt 1 lc rgb 'orange' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 6   lt 1 lc rgb 'turquoise' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 7   lt 1 lc rgb 'black' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 8   lt 1 lc rgb 'green' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 9   lt 1 lc rgb 'blue' pt 7 lw 2 ps variable\n");

  // covarying BP and covaryingNBPs (covBP & covNBP) / null convarying bpairs (covRBP)
  fprintf(pipe, "set ylabel '# covarying pairs'\n");
  fprintf(pipe, "set xlabel '# null covarying basepairs'\n");
  fprintf(pipe, "set yrange [0:%f]\n", maxBP);
  fprintf(pipe, "set xrange [0:%f]\n", maxcovRBP);
  fprintf(pipe, "plot '%s' u 2:4 with linespoints ls 8, ", nullcovfile);
  fprintf(pipe, "'%s' u 2:6 with linespoints ls 7\n", nullcovfile);
  
  // % covarying bpairs (covBPf) / % null convarying bpairs (covRBPf)
  fprintf(pipe, "set size ratio -1\n");
  fprintf(pipe, "set ylabel '%% covarying basepairs'\n");
  fprintf(pipe, "set xlabel '%% null covarying basepairs'\n");
  fprintf(pipe, "set yrange [0:100]\n");
  fprintf(pipe, "set xrange [0:100]\n");
  fprintf(pipe, "plot '%s' u 3:5 with linespoints ls 8\n", nullcovfile);
  
  pclose(pipe);
  
  free(outplot);
  free(filename);

  return eslOK;
}

int              
cov_DotPlot(char *gnuplot, char *dplotfile, ESL_MSA *msa, int *ct, struct mutual_s *mi, int *msamap, HITLIST *hitlist, int dosvg, int verbose, char *errbuf)
{
  FILE    *pipe;
  char    *filename = NULL;
  char    *outplot = NULL;
  double   pointsize;
  double   ps_max = 0.40;
  double   ps_min = 0.0003;
  int      L = msamap[msa->alen-1]+1;
  int      h;           /* index for hitlist */
  int      i, ipair;
  int      ih, jh;
  

  esl_FileTail(dplotfile, FALSE, &filename);

  pipe = popen(gnuplot, "w");
  
  if (dosvg) {
    ps_max = 0.40;
    ps_min = 0.0003;
    esl_sprintf(&outplot, "%s.svg", dplotfile);
    fprintf(pipe, "set terminal svg fname 'Verdana' fsize 10 \n");
  }
  else {
    ps_max = 1.40;
    ps_min = 0.30;
    esl_sprintf(&outplot, "%s.ps", dplotfile);
    fprintf(pipe, "set terminal postscript color 14\n");
  }
  fprintf(pipe, "set output '%s'\n", outplot);
  pointsize = (mi->maxCOV > 0.)? ps_max/mi->maxCOV : ps_min;

  /* matlab's 'jet' colormap scale */
  fprintf(pipe, "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0)\n");

  fprintf(pipe, "unset key\n");
  fprintf(pipe, "set size ratio -1\n");
  fprintf(pipe, "set pointsize %f\n", pointsize);
  fprintf(pipe, "set title '%s' \n", filename);
  fprintf(pipe, "set ylabel 'alignment position'\n");
  fprintf(pipe, "set xlabel 'alignment position'\n");

  fprintf(pipe, "set style line 1   lt 1 lc rgb 'grey' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 2   lt 1 lc rgb 'brown' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 3   lt 1 lc rgb 'cyan' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 4   lt 1 lc rgb 'red' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 5   lt 1 lc rgb 'orange' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 6   lt 1 lc rgb 'turquoise' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 7   lt 1 lc rgb 'black' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 8   lt 1 lc rgb 'green' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 9   lt 1 lc rgb 'blue' pt 7 lw 2 ps variable\n");

  fprintf(pipe, "set yrange [1:%d]\n", L);
  fprintf(pipe, "set xrange [1:%d]\n", L);

  fprintf(pipe, "set multiplot\n");

  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "fun(x)=x\n");
  fprintf(pipe, "plot fun(x) with lines ls 7\n");

  // the actual ct
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 9\n");
  for (i = 1; i <= msa->alen; i ++) {
    ipair = ct[i];
    if (ipair > 0) {
      fprintf(pipe, "%d %d %f\n", msamap[i-1]+1,     msamap[ipair-1]+1, (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
      fprintf(pipe, "%d %d %f\n", msamap[ipair-1]+1, msamap[i-1]+1,     (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
    }	
  } 
  fprintf(pipe, "e\n");

  // the covarying basepairs
  //fprintf(pipe, "plot '-' u 1:2:3 with image \n");
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 8 \n");
  for (h = 0; h < hitlist->nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    if (hitlist->hit[h].is_bpair) {
      fprintf(pipe, "%d %d %f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc);
      fprintf(pipe, "%d %d %f\n", msamap[jh]+1, msamap[ih]+1, hitlist->hit[h].sc);
    }	
  } 
  fprintf(pipe, "e\n");

  // covarying pairs compatible with the given structure
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 5\n");
  for (h = 0; h < hitlist->nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    if (hitlist->hit[h].is_compatible) {
      fprintf(pipe, "%d %d %f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc);	
      fprintf(pipe, "%d %d %f\n", msamap[jh]+1, msamap[ih]+1, hitlist->hit[h].sc);	
    }
  } 
  fprintf(pipe, "e\n");
  
  // covarying pairs incompatible with the given structure
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 7\n");
  for (h = 0; h < hitlist->nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    if (!hitlist->hit[h].is_bpair && !hitlist->hit[h].is_compatible) {
      fprintf(pipe, "%d %d %f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc);	
      fprintf(pipe, "%d %d %f\n", msamap[jh]+1, msamap[ih]+1, hitlist->hit[h].sc);	
    }
  } 
  fprintf(pipe, "e\n");
  
  pclose(pipe);
  
  free(outplot);
  free(filename);
  return eslOK;
}

int
cov_R2R(char *r2rfile, char *r2rversion, int r2rall, ESL_MSA **ret_msa, int *ct, int *msamap, HITLIST *hitlist, int makepdf, int makesvg, int verbose, char *errbuf)
 {
  ESLX_MSAFILE *afp = NULL;
  FILE         *fp = NULL;
  ESL_MSA      *msa = *ret_msa;
  ESL_MSA      *r2rmsa = NULL;
  char          tmpinfile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char          covtag[12] = "cov_SS_cons";
  char         *args = NULL;
  char         *s = NULL;
  char         *ssstr = NULL;
  char         *covstr = NULL;
  char         *tok;
  int           found;
  int           i;
  int           h;
  int           ih, jh;
  int           tagidx;
  int           status;
 
  /* first modify the ss to a simple <> format. R2R cannot deal with fullwuss 
   */
  ESL_ALLOC(ssstr, sizeof(char) * (msa->alen+1));
  esl_ct2simplewuss(ct, msa->alen, ssstr);
  
  /* replace the 'SS_cons' GC line with the new ss */
  esl_sprintf(&(msa->ss_cons), "%s", ssstr);  
  
  /* R2R input and output in PFAM format (STOCKHOLM in one single block) */
  if ((status = esl_tmpfile_named(tmpinfile,  &fp))                      != eslOK) ESL_XFAIL(status, errbuf, "failed to create input file");
  if ((status = eslx_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_PFAM)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write PFAM file\n");
  fclose(fp);
  
  /* run R2R */
  if ("R2RDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("R2RDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpoutfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  esl_sprintf(&args, "%s/%s/src/r2r --GSC-weighted-consensus %s %s 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1", s, r2rversion, tmpinfile, tmpoutfile);
  system(args);
  fclose(fp);
 
  /* convert output to r2rmsa */
  if (eslx_msafile_Open(NULL, tmpoutfile, NULL, eslMSAFILE_PFAM, NULL, &afp) != eslOK) eslx_msafile_OpenFailure(afp, status);
  afp->format = eslMSAFILE_PFAM;
  if (eslx_msafile_Read(afp, &r2rmsa) != eslOK) eslx_msafile_ReadFailure(afp, status);
  eslx_msafile_Close(afp);

  /* modify the cov_cons_ss line acording to our hitlist */
  if (msa->alen != r2rmsa->alen) ESL_XFAIL(eslFAIL, errbuf, "r2r has modified the alignment\n");
  for (i = 1; i <= msa->alen; i ++) {
    found = FALSE;
    for (h = 0; h < hitlist->nhit; h ++) {
      ih = hitlist->hit[h].i+1;
      jh = hitlist->hit[h].j+1;

      if ((i == ih || i == jh) && hitlist->hit[h].is_bpair) { 
	esl_sprintf(&tok, "2"); 
	found = TRUE; 
      }
      
      if (found) break;
    }
    if (!found) esl_sprintf(&tok, ".");  

    if (i == 1) esl_sprintf(&covstr, "%s", tok);
    else        esl_sprintf(&covstr, "%s%s", covstr, tok);
  }
 
  /* add line #=GF R2R keep allpairs 
   * so that it does not truncate ss.
   * cannot use the standard esl_msa_addGF:
   *             esl_msa_AddGF(msa, "R2R", -1, " keep allpairs", -1);
   * since it does not parse with r2r
   *
   * turns out the above solution can only deal with the  <> annotation
   */
  if (r2rall) 
    esl_msa_AddGF(r2rmsa, "R2R keep all", -1, "", -1);
  else
    esl_msa_AddGF(r2rmsa, "R2R keep allpairs", -1, "", -1);
  
  /* replace the r2r 'cov_SS_cons' GC line with our own */
  for (tagidx = 0; tagidx < r2rmsa->ngc; tagidx++)
    if (strcmp(r2rmsa->gc_tag[tagidx], covtag) == 0) break;
  if (tagidx == r2rmsa->ngc) {
    ESL_REALLOC(r2rmsa->gc_tag, (r2rmsa->ngc+1) * sizeof(char **));
    ESL_REALLOC(r2rmsa->gc,     (r2rmsa->ngc+1) * sizeof(char **));
    r2rmsa->gc[r2rmsa->ngc] = NULL;
    r2rmsa->ngc++;
  }
  if ((status = esl_strdup(covtag, -1, &(r2rmsa->gc_tag[tagidx]))) != eslOK) goto ERROR;
  esl_sprintf(&(r2rmsa->gc[tagidx]), "%s", covstr);

  if (verbose) eslx_msafile_Write(stdout, r2rmsa, eslMSAFILE_PFAM);
  
  /* write the R2R annotated to PFAM format */
  if (r2rfile) {
    if ((fp = fopen(r2rfile, "w")) == NULL) esl_fatal("Failed to open output file %s", r2rfile);
    eslx_msafile_Write(fp, r2rmsa, eslMSAFILE_PFAM);
    fclose(fp);
    
    /* produce the R2R pdf */
    if (makepdf) {
      status = cov_R2Rpdf(r2rfile, r2rversion, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
    if (makesvg) {
      status = cov_R2Rsvg(r2rfile, r2rversion, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
  }
  
  remove(tmpinfile);
  remove(tmpoutfile);
  
  *ret_msa = r2rmsa;
  esl_msa_Destroy(msa);
  if (args) free(args);
  if (ssstr)  free(ssstr);
  if (covstr) free(covstr);
  return eslOK;

 ERROR:
  remove(tmpinfile);
  remove(tmpoutfile);
  
  if (msa)    esl_msa_Destroy(msa);
  if (r2rmsa) esl_msa_Destroy(r2rmsa);
  if (args)   free(args);
  if (ssstr)  free(ssstr);
  if (covstr) free(covstr);
  return status;
}

int
cov_R2Rpdf(char *r2rfile, char *r2rversion, int verbose, char *errbuf)
{
  char *r2rpdf = NULL;
  char *args = NULL;
  char *s = NULL;

  /* produce the R2R pdf */
  if ("R2RDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("R2RDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&r2rpdf, "%s.pdf", r2rfile);
  esl_sprintf(&args, "%s/%s/src/r2r %s %s >/dev/null", s, r2rversion, r2rfile, r2rpdf);
  //esl_sprintf(&args, "%s/%s/src/r2r %s %s ", s, r2rversion, r2rfile, r2rpdf);
  system(args);
  
  free(args);
  free(r2rpdf);
  
  return eslOK;
 }

int
cov_R2Rsvg(char *r2rfile, char *r2rversion, int verbose, char *errbuf)
{
  char *r2rsvg = NULL;
  char *args = NULL;
  char *s = NULL;

  /* produce the R2R svg */
  if ("R2RDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("R2RDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&r2rsvg, "%s.svg", r2rfile);
  esl_sprintf(&args, "%s/%s/src/r2r %s %s >/dev/null", s, r2rversion, r2rfile, r2rsvg);
  //esl_sprintf(&args, "%s/%s/src/r2r %s %s ", s, r2rversion, r2rfile, r2rsvg);
  system(args);
  
  free(args);
  free(r2rsvg);
  
  return eslOK;
 }

int
cov_ExpandCT(char *r2rfile, int r2rall, ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, int minloop, enum grammar_e G, int verbose, char *errbuf)
{
  FILE *fp = NULL;
  char *ss = NULL;
  int   tagidx;
  int   L = msa->alen;
  int   status;
  
  /* replace any R2R line with a new one not tab delimited 
   * R2R chockes on tab delimited lines
   */
  for (tagidx = 0; tagidx < msa->ngf; tagidx++) 
    if (strcmp(msa->gf_tag[tagidx], "R2R") == 0) {
      esl_sprintf(&(msa->gf_tag[tagidx]), "R2R %s", msa->gf[tagidx]);
      esl_sprintf(&(msa->gf[tagidx]), "");
  }

#if 0 // naive method
  status = cov_ExpandCT_Naive(msa, *ret_ct, minloop, verbose, errbuf);
#else // covariance-constrain CYK using a probabilistic grammar
  status = cov_ExpandCT_CCCYK(r, msa, ret_ct, G, minloop, verbose, errbuf);
#endif
  if (status != eslOK) goto ERROR;

  /* replace the 'SS_cons' GC line with the new ss */
  ESL_ALLOC(ss, sizeof(char) * (L+1));
  esl_ct2simplewuss(*ret_ct, L, ss);
  esl_sprintf(&(msa->ss_cons), "%s", ss);  
  //printf("ss:%s\n", msa->ss_cons);

  if ((fp = fopen(r2rfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to open output file %s", r2rfile);
  eslx_msafile_Write(fp, msa, eslMSAFILE_PFAM);
  fclose(fp);
  
  free(ss);
  return eslOK;

 ERROR:
  if (ss) free(ss);
  return status;
}

int
cov_ExpandCT_Naive(ESL_MSA *msa, int *ct, int minloop, int verbose, char *errbuf)
{
  char  tag[10] = "cons";
  char *cons;
  int   tagidx;
  int   L = msa->alen;
  int   i, j, d;
  
  /* get the line #=GC cons */
  for (tagidx = 0; tagidx < msa->ngc; tagidx++)
    if (strcmp(msa->gc_tag[tagidx], tag) == 0) break;
  if (tagidx == msa->ngc) return eslOK; // no cons line to expand the CT
  cons =  msa->gc[tagidx];

  for (j = 0; j < L; j++)
    for (d = 0; d <= j; d++)
      {
	i = j-d+1;

	if ( d >= minloop && is_stacked_pair(i+1, j+1, L, ct) && is_cannonical_pair(cons[i], cons[j]) )
	  { // add this pair
	    //printf("%c %c | %d %d %d\n", cons[i], cons[j], i, j, d);
	    ct[i+1] = j+1;
	    ct[j+1] = i+1;
	  }	
      }

return eslOK;
}

int
cov_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct,  enum grammar_e G, int minloop, int verbose, char *errbuf)
{
  char  tag[10] = "RF_cons";
  char    *rfline = NULL;
  ESL_SQ  *sq = NULL;
  int     *ct = *ret_ct;
  int     *cct = NULL;
  SCVAL    sc;
  float    idthresh = 0.3;
  int      tagidx;
  int      status;
  
  /* create an RF sequence */
  ESL_ALLOC(rfline, sizeof(char) * (msa->alen+1));
  esl_msa_ReasonableRF(msa, idthresh, TRUE, rfline);
  sq = esl_sq_CreateFrom(msa->name,  rfline, msa->desc, msa->acc, msa->ss_cons); 
  
  //printf("sq:%s\n", sq->seq);
  esl_sq_Digitize((const ESL_ALPHABET *)msa->abc, sq);
 
  cykcov_remove_inconsistencies(sq, ct, minloop);
 
 /* calculate the convariance-constrain CYK structure using a probabilistic grammar */
  status = COCOCYK(r, G, sq, ct, &cct, &sc, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  if (verbose) printf("coco-cyk score = %f\n", sc);

  if (cct) {
    free(ct); ct = NULL;
    *ret_ct = cct;
  }
  else *ret_ct = ct;

  if (rfline) free(rfline);
  esl_sq_Destroy(sq);
  return eslOK;

 ERROR:
  if (sq) esl_sq_Destroy(sq);
  if (rfline) free(rfline);
  if (cct) free(cct);
  return status;
}

int
cov_ranklist_Score2Bin(RANKLIST *ranklist, double x, int *ret_b)
{
  int status;

  if (! isfinite(x)) ESL_XEXCEPTION(eslERANGE, "value added to histogram is not finite");

  x = ceil( ((x - ranklist->bmin) / ranklist->w) - 1.); 
  
  /* x is now the bin number as a double, which we will convert to
   * int. Because x is a double (64-bit), we know all ints are exactly
   * represented.  Check for under/overflow before conversion.
   */
  if (x < (double) INT_MIN || x > (double) INT_MAX) 
    ESL_XEXCEPTION(eslERANGE, "value %f isn't going to fit in histogram", x);

  *ret_b = (int) x;
  return eslOK;

 ERROR:
  *ret_b = 0;
  return status;
}



/*---------------- internal functions --------------------- */

static int
is_wc(int x, int y) 
{
  if (x+y == 3 || x+y == 5) return TRUE;
  
  return FALSE;
}


static int
is_stacked_pair(int i, int j, int L, int *ct) 
{
  int is_stacked = FALSE;

  if (ct[i] > 0 || ct[j] > 0) return FALSE; // both have to be unpaired

  if (ct[i+1] == j-1 && ct[j-1] == i+1                  ) is_stacked = TRUE;
  if (ct[i-1] == j+1 && ct[j+1] == i-1 && i > 1 && j < L) is_stacked = TRUE;
	
  return is_stacked;
}

static int
number_pairs(int L, int *ct) 
{
  int nbp = 0;
  int i;

  for (i = 1; i <= L; i ++) 
    if (ct[i] > 0 && i < ct[i]) nbp ++;
 	
  return nbp;
}

static int
is_cannonical_pair(char nti, char ntj) 
{
  int is_cannonical = FALSE;

  if (nti == 'A' && ntj == 'U') is_cannonical = TRUE;
  if (nti == 'C' && ntj == 'G') is_cannonical = TRUE;
  if (nti == 'G' && ntj == 'C') is_cannonical = TRUE;
  if (nti == 'G' && ntj == 'U') is_cannonical = TRUE;
  if (nti == 'G' && ntj == 'Y') is_cannonical = TRUE;
  if (nti == 'U' && ntj == 'A') is_cannonical = TRUE;
  if (nti == 'U' && ntj == 'G') is_cannonical = TRUE;
  if (nti == 'U' && ntj == 'R') is_cannonical = TRUE;
  if (nti == 'R' && ntj == 'Y') is_cannonical = TRUE;
  if (nti == 'Y' && ntj == 'R') is_cannonical = TRUE;

  return is_cannonical;
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

static int
cykcov_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop)
{
  int L = sq->n;
  int n;
  int ipair;
  int i;
  int x;

  /* remove covariation that correspond to gap-gap in
   * the RF sequence */
  for (i = 1; i <= L; i++) {
    ipair = ct[i];
    if (ipair > 0 && sq->dsq[i] >= NB && sq->dsq[ipair] >= NB) { // remove this covariation
      ct[i]     = 0;
      ct[ipair] = 0;
    }
  }

  /* remove covariation that are closer than minloop in RF 
   */
  for (i = 1; i <= L; i++) {
    ipair = ct[i];
    if (ipair == 0 || ipair < i) continue;

    n = 0;
    for (x = i+1; x < ipair; x ++) 
      if (sq->dsq[x] < NB) n ++;

    if (n < minloop-2) { // remove this covariation
      ct[i]     = 0;
      ct[ipair] = 0;
    }
  }

  return eslOK;
}

static int
cov_histogram_plotsurvival(FILE *pipe, ESL_HISTOGRAM *h, int logscale, int style1, int style2)
{
  int      i;
  uint64_t c = 0;
  double   esum;
  double   ai;
 
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "plot '-' using 1:2 with points ls %d\n", style1);
  if (h->obs[h->imax] > 1) 
    if (fprintf(pipe, "%f\t%f\n", 
		h->xmax, (logscale)? -log((double)h->Nc) : 1.0/(double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
  for (i = h->imax; i >= h->imin; i--)
    {
      if (h->obs[i] > 0) {
	c   += h->obs[i];
	ai = esl_histogram_Bin2LBound(h, i);
 	if (fprintf(pipe, "%f\t%f\n", 
		    ai, (logscale)? log((double)c)-log((double)h->Nc) : (double)c/(double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  fprintf(pipe, "e\n");

  /* The expected binned counts:
   */
  if (h->expect != NULL) 
    {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");
      fprintf(pipe, "plot '-' using 1:2 with points ls %d\n", style2);
      
      esum = 0.;
      for (i = h->nb-1; i >= 0; i--)
	{
	  if (h->expect[i] > 0.) { 
	    esum += h->expect[i];        /* some worry about 1+eps=1 problem here */
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(pipe, "%f\t%f\n", 
			ai, (logscale)? log(esum)-log((double)h->Nc) : esum/(double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      fprintf(pipe, "e\n"); 
    }

  return eslOK;
}
