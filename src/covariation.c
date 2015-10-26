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
#include "esl_exponential.h"
#include "esl_gamma.h"
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

static int    is_wc(int x, int y);
static int    is_stacked_pair(int i, int j, int L, int *ct);
static int    number_pairs(int L, int *ct);
static int    is_cannonical_pair(char nti, char ntj);
static int    mutual_naive_ppij(ESL_RANDOMNESS *r, int i, int j, ESL_MSA *msa, struct mutual_s *mi,
				int donull2b, double tol, int verbose, char *errbuf);
static int    shuffle_null2b_col(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, int nseq, int *col, int *paircol, int **ret_shcol, char *errbuf);
static int    shuffle_col(ESL_RANDOMNESS *r, int nseq, int *useme, int *col, int **ret_shcol, char *errbuf);
static int    mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi,
				    ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);
static int    cykcov_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop);
static double cov2evalue(struct data_s *data, double cov, int Nc, ESL_HISTOGRAM *h);
static double evalue2cov(struct data_s *data, double eval, int Nc, ESL_HISTOGRAM *h);
static int    cov_histogram_plotdensity(FILE *pipe, ESL_HISTOGRAM *h, char *key, double posx, double posy, int logval, int subsample, int style1, int style2);
static int    cov_histogram_plotsurvival(FILE *pipe, ESL_HISTOGRAM *h, char *key, double posx, double posy, int logval, int subsample, int style1, int style2);
static int    cov_histogram_plotexpectsurv(FILE *pipe, int Nc, ESL_HISTOGRAM *h, char *key, double posx, double posy, int logval, int subsample, 
					   int linespoints, int style1, int style2);
static int    cov_plot_lineatexpcov(FILE *pipe, struct data_s *data, double expsurv, int Nc, ESL_HISTOGRAM *h,
				    double ymin, double ymax, char *key, double offx, double offy, int style);
static int    cov_histogram_bin2expectsurv(int i, ESL_HISTOGRAM *h, double *ret_expsurv);

int                 
cov_Calculate(struct data_s *data, ESL_MSA *msa, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, int analize)
{
  RANKLIST      *ranklist = NULL;
  HITLIST       *hitlist = NULL;
  COVCLASS       covclass = data->mi->class;
  int            status;
  
  /* Calculate the covariation matrix */
  if ( !(data->covtype == RAF  || data->covtype == RAFp  || data->covtype == RAFa ||
	 data->covtype == RAFS || data->covtype == RAFSp || data->covtype == RAFSa ) ) {
    status = cov_Probs(data->r, msa, data->T, data->ribosum, data->mi, data->method, data->donull2b, data->tol, data->verbose, data->errbuf);
    if (status != eslOK) goto ERROR;
  }

  switch(data->covtype) {
  case CHIa: 
    status = cov_CalculateCHI         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case CHIp:
    status = cov_CalculateCHI         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;  
    break;
  case CHI: 
    status = cov_CalculateCHI         (covclass, data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;     
    break;
  case GTa: 
    status = cov_CalculateGT          (covclass, data, FALSE,  NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case GTp: 
    status = cov_CalculateGT          (covclass, data, FALSE,  NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
     if (status != eslOK) goto ERROR; 
     break;
  case GT: 
     status = cov_CalculateGT          (covclass, data, analize, &ranklist, &hitlist);
     if (status != eslOK) goto ERROR;
     break;
  case MIa: 
    status = cov_CalculateMI          (covclass, data, FALSE,  NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case MIp: 
    status = cov_CalculateMI          (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case MI: 
    status = cov_CalculateMI          (covclass, data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case MIra: 
    status = cov_CalculateMIr         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case MIrp:
    status = cov_CalculateMIr         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;  
    break;
  case MIr: 
    status = cov_CalculateMIr         (covclass, data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case MIga: 
    status = cov_CalculateMIg         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case MIgp:
    status = cov_CalculateMIg         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;  
    break;
  case MIg: 
    status = cov_CalculateMIg          (covclass, data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case OMESa: 
    status = cov_CalculateOMES        (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case OMESp: 
    status = cov_CalculateOMES        (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case OMES: 
    status = cov_CalculateOMES        (covclass, data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
 case RAFa: 
   status = cov_CalculateRAF         (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case RAFp: 
    status = cov_CalculateRAF         (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case RAF: 
    status = cov_CalculateRAF         (covclass, data, msa, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case RAFSa: 
    status = cov_CalculateRAFS        (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(ASC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case RAFSp: 
    status = cov_CalculateRAFS        (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = cov_CalculateCOVCorrected(APC,      data, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR; 
    break;
  case RAFS: 
    status = cov_CalculateRAFS        (covclass, data, msa, analize, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  default:
    ESL_XFAIL(eslFAIL, data->errbuf, "wrong covariation type\n");
    break;
  }
  if (data->mode != RANSS) fprintf(data->sumfp, "\n");   
  
  if (data->mode == GIVSS) { // do the plots only for GIVSS
    status = cov_DotPlot(data->gnuplot, data->dplotfile, msa, data->ct, data->mi, data->msamap, data->firstpos, hitlist, TRUE, data->verbose, data->errbuf);
    if  (status != eslOK) goto ERROR;
    status = cov_DotPlot(data->gnuplot, data->dplotfile, msa, data->ct, data->mi, data->msamap, data->firstpos, hitlist, FALSE, data->verbose, data->errbuf);
    if  (status != eslOK) goto ERROR;

    status = cov_R2R(data->R2Rfile, data->R2Rversion, data->R2Rall, msa, data->ct, hitlist, TRUE, TRUE, data->verbose, data->errbuf);
    if  (status != eslOK) goto ERROR;
  }
  
  if (ret_ranklist) *ret_ranklist = ranklist; else if (ranklist) cov_FreeRankList(ranklist);
  if (ret_hitlist)  *ret_hitlist = hitlist;   else if (hitlist)  cov_FreeHitList(hitlist);
  return eslOK;
  
 ERROR:
  if (ranklist) cov_FreeRankList(ranklist);
  if (hitlist)  cov_FreeHitList(hitlist);
  return status;
}

int                 
cov_Probs(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
	  METHOD method, int donull2b, double tol, int verbose, char *errbuf)
{
  int i, j;
  int x, y;
  int K = msa->abc->K;
  int status;

  switch(method) {
  case NAIVE:
    status = cov_NaivePP(r, msa, mi, donull2b, tol, verbose, errbuf);
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

  status = cov_Marginals(mi, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;    

  status = cov_ValidateProbs(mi, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  if (verbose) {
    for (i = 0; i < mi->alen-1; i ++) {
      for (j = i+1; j < mi->alen; j ++) {
	if ((i==5&&j==118)||verbose) {
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
cov_CalculateCHI(COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  char            *errbuf = data->errbuf;
  int              verbose = data->verbose;
  int              i, j;
  int              status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateCHI_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateCHI_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh || mi->alen <= mi->alenthresh)
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

  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }

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
cov_CalculateOMES(COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  char            *errbuf = data->errbuf;
  int              verbose = data->verbose;
  int              i, j;
  int              status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateOMES_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateOMES_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
   if (mi->nseq <= mi->nseqthresh || mi->alen <= mi->alenthresh)
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
  
  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }

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
cov_CalculateGT(COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  char            *errbuf = data->errbuf;
  int              verbose = data->verbose;
  int              i, j;
  int              status = eslOK;
 
  switch (covclass) {
  case C16:
    status = cov_CalculateGT_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateGT_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh || mi->alen <= mi->alenthresh) {
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
  
  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }
  
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
	  gt += (exp > 0. && obs > 0.) ? obs * log (obs / exp) : 0.0;
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
      
      gt += (exp_wc  > 0. && obs_wc  > 0.) ? obs_wc  * log (obs_wc  / exp_wc)  : 0.0;
      gt += (exp_nwc > 0. && obs_nwc > 0.) ? obs_nwc * log (obs_nwc / exp_nwc) : 0.0;
      gt *= 2.0;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = gt;
      if (gt < mi->minCOV) mi->minCOV = gt;
      if (gt > mi->maxCOV) mi->maxCOV = gt;
    }

  return status;
}


int                 
cov_CalculateMI(COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  char            *errbuf = data->errbuf;
  int              verbose = data->verbose;
  int              i, j;
  int              status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateMI_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateMI_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh || mi->alen <= mi->alenthresh)
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
  
  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }
  
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
	  mutinf += (mi->pp[i][j][IDX(x,y,K)] > 0.0 && mi->pm[i][x] > 0.0 && mi->pm[j][y] > 0.0)? 
	    mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->pm[i][x]) - log(mi->pm[j][y]) ) : 0.0;
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
      
      mutinf += (pij_wc  > 0.)? pij_wc  * ( log(pij_wc)  - log(qij_wc)  ) : 0.0;
      mutinf += (pij_nwc > 0.)? pij_nwc * ( log(pij_nwc) - log(qij_nwc) ) : 0.0;
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = mutinf;
      if (mutinf < mi->minCOV) mi->minCOV = mutinf;
      if (mutinf > mi->maxCOV) mi->maxCOV = mutinf; 
    }

  return status;
}


int                 
cov_CalculateMIr(COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  char            *errbuf = data->errbuf;
  int              verbose = data->verbose;
  int              i, j;
  int              status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateMIr_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateMIr_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh || mi->alen <= mi->alenthresh)
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
  
  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }
  
  return status;
  
 ERROR:
  return status;
}

int                 
cov_CalculateMIr_C16(struct mutual_s *mi, int verbose, char *errbuf)
{
  double mutinf, HH;
  double tol = 1e-2;
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
	  HH -= (mi->pp[i][j][IDX(x,y,K)] > 0.0)? mi->pp[i][j][IDX(x,y,K)] * log(mi->pp[i][j][IDX(x,y,K)]) : 0.0;
	  mutinf += (mi->pp[i][j][IDX(x,y,K)] > 0.0  && mi->pm[i][x] > 0.0 && mi->pm[j][y] > 0.0)? 
	    mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->pm[i][x]) - log(mi->pm[j][y]) ) : 0.0;
	}	  
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = (HH > tol)? mutinf/HH : 0.0;

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
  double tol = 1e-2;
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
      
      HH -= (pij_wc  > 0.)? pij_wc  * log(pij_wc)  : 0.0;
      HH -= (pij_nwc > 0.)? pij_nwc * log(pij_nwc) : 0.0;

      mutinf += (pij_wc  > 0. && qij_wc  > 0.0)? pij_wc  * ( log(pij_wc)  - log(qij_wc)  ) : 0.0;
      mutinf += (pij_nwc > 0. && qij_nwc > 0.0)? pij_nwc * ( log(pij_nwc) - log(qij_nwc) ) : 0.0;
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = (HH > tol)? mutinf/HH : 0.0;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j]; 
    }

  return status;
}


int                 
cov_CalculateMIg(COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  char            *errbuf = data->errbuf;
  int              verbose = data->verbose;
  int              i, j;
  int              status = eslOK;
  
  switch (covclass) {
  case C16:
    status = cov_CalculateMIg_C16  (mi, verbose, errbuf);
    break;
  case C2:
    status = cov_CalculateMIg_C2   (mi, verbose, errbuf);
    break;
  case CSELECT:
    if (mi->nseq <= mi->nseqthresh || mi->alen <= mi->alenthresh)
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
  
  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }
  
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
	  mutinf += (mi->pp[i][j][IDX(x,y,K)] > 0.0 && mi->pm[i][x] > 0.0  && mi->pm[j][y] > 0.0)? 
	    mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->pm[i][x]) - log(mi->pm[j][y]) ) : 0.0;
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
      
      mutinf += (pij_wc  > 0.0)? pij_wc  * ( log(pij_wc)  - log(qij_wc)  ) : 0.0;
      mutinf += (pij_nwc > 0.0)? pij_nwc * ( log(pij_nwc) - log(qij_nwc) ) : 0.0;
      
      /* the negative correction of gaps */
      mutinf -= (mi->nseff[i][j] > 0)? (double)mi->ngap[i][j] / (double)mi->nseff[i][j] : 0.0;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = mutinf;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j]; 
   }

  return status;
}

int                 
cov_CalculateRAF(COVCLASS covclass, struct data_s *data, ESL_MSA *msa,  int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  int              verbose = data->verbose;
  double           psi = 1.0;
  double           cij, qij;
  int              i, j;
  int              s1, s2;
  int              ai, aj, bi, bj;
  int              status = eslOK;
  
  cov_ReuseCOV(mi, RAF, C2);

  // RAF
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      cij = 0.0;
      qij = 0.0;
      
      for (s1 = 0; s1 < msa->nseq; s1 ++) {
	ai = msa->ax[s1][i+1];
	aj = msa->ax[s1][j+1];

	if ( !is_wc(ai,aj) ) qij += 1.0;
	
	for (s2 = s1+1; s2 < msa->nseq; s2 ++) {
	  bi = msa->ax[s2][i+1];
	  bj = msa->ax[s2][j+1];

	  if ( is_wc(ai,aj) && is_wc(bi,bj) ) 
	    {
	      if      (ai != bi && aj != bj) cij += 2.0;
	      else if (ai != bi || aj != bj) cij += 1.0;
	    }
	}
      }
      qij /= msa->nseq;
      
      cij /= (msa->nseq > 1)? (double)msa->nseq * ((double)msa->nseq-1.0): 1.0;
      cij *= 2.0;
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = cij - psi * qij;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j]; 

    }

  if (verbose) {
    printf("RAF[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("RAF[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }

  return status;
  
 ERROR:
  return status;
}

int                 
cov_CalculateRAFS(COVCLASS covclass, struct data_s *data, ESL_MSA *msa, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  int              verbose = data->verbose;
  ESL_DMATRIX     *bij = NULL;
  double           bijs;
  int              i, j;
  int              status = eslOK;

  cov_ReuseCOV(mi, RAF, C2);
  cov_CalculateRAF(RAF, data, msa, FALSE, NULL, NULL);

  // RAFS
  bij = esl_dmatrix_Clone(mi->COV);
  cov_ReuseCOV(mi, RAFS, C2);
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      
      bijs = 2.0 * bij->mx[i][j];
      if (i > 0 && j < mi->alen-1) bijs += bij->mx[i-1][j+1];
      if (j > 0 && i < mi->alen-1) bijs += bij->mx[i+1][j-1];
      bijs *= 0.25;

      mi->COV->mx[i][j] = mi->COV->mx[j][i] = bijs;     
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j]; 

    }

  if (verbose) {
    printf("RAFS[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("RAF[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) {
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
    if (status != eslOK) goto ERROR;
  }

  esl_dmatrix_Destroy(bij);
  return status;
  
 ERROR:
  if (bij) esl_dmatrix_Destroy(bij);
  return status;
}


int                 
cov_CalculateCOVCorrected(CORRTYPE corrtype, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  char            *errbuf = data->errbuf;
  int              verbose = data->verbose;
  char            *type = NULL;
  char            *covtype = NULL;
  ESL_DMATRIX     *COV  = NULL;
  double          *COVx = NULL;
  double           COVavg = 0.0;
  int              i, j;
  int              status = eslOK;
  
  cov_COVTYPEString(&type, mi->type, errbuf);

  switch(corrtype) {
  case APC: esl_sprintf(&covtype, "%sp", type); break;
  case ASC: esl_sprintf(&covtype, "%sa", type); break;
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

      if (isnan(mi->COV->mx[i][j])) ESL_XFAIL(eslFAIL, errbuf, "bad covariation\n");

      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
    }

  if (verbose) {
    printf("%s-[%f,%f] \n", covtype,  mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if ((i==0&&j==44)||(i==2&&j==43)||(i==7&&j==34)||(i==9&&j==32)) printf("%s-[%d][%d] = %f | COV %f | COVx %f COVy %f | COVavg %f\n", 
				 covtype, i, j, mi->COV->mx[i][j], COV->mx[i][j], COVx[i], COVx[j], COVavg);
      } 
  }

  if (analyze) 
    status = cov_SignificantPairs_Ranking(data, ret_ranklist, ret_hitlist);
  if (status != eslOK) goto ERROR;
  
  free(type);
  free(covtype);
  esl_dmatrix_Destroy(COV);
  free(COVx);
  return status;

 ERROR:
  if (type)    free(type);
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
  case Eval:     esl_sprintf(ret_threshtype, "Eval");    break;
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
  case RAF:   esl_sprintf(ret_covtype, "RAF");   break; 
  case RAFS:  esl_sprintf(ret_covtype, "RAFS");  break; 

  case CHIp:  esl_sprintf(ret_covtype, "CHIp");  break;
  case GTp:   esl_sprintf(ret_covtype, "GTp");   break;
  case OMESp: esl_sprintf(ret_covtype, "OMESp"); break;
  case MIp:   esl_sprintf(ret_covtype, "MIp");   break;
  case MIrp:  esl_sprintf(ret_covtype, "MIrp");  break; 
  case MIgp:  esl_sprintf(ret_covtype, "MIgp");  break; 
  case RAFp:  esl_sprintf(ret_covtype, "RAFp");  break; 
  case RAFSp: esl_sprintf(ret_covtype, "RAFSp"); break; 

  case CHIa:  esl_sprintf(ret_covtype, "CHIa");  break;
  case GTa:   esl_sprintf(ret_covtype, "GTa");   break;
  case OMESa: esl_sprintf(ret_covtype, "OMESa"); break;
  case MIa:   esl_sprintf(ret_covtype, "MIa");   break;
  case MIra:  esl_sprintf(ret_covtype, "MIra");  break; 
  case MIga:  esl_sprintf(ret_covtype, "MIga");  break; 
  case RAFa:  esl_sprintf(ret_covtype, "RAFa");  break; 
  case RAFSa: esl_sprintf(ret_covtype, "RAFSa"); break; 

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
  else if (!esl_strcmp(covtype, "RAF"))    type = RAF;
  else if (!esl_strcmp(covtype, "RAFS"))   type = RAFS;

  else if (!esl_strcmp(covtype, "CHIp"))   type = CHIp;
  else if (!esl_strcmp(covtype, "GTp"))    type = GTp;
  else if (!esl_strcmp(covtype, "OMESp"))  type = OMESp;
  else if (!esl_strcmp(covtype, "MIp"))    type = MIp;
  else if (!esl_strcmp(covtype, "MIrp"))   type = MIrp;
  else if (!esl_strcmp(covtype, "MIgp"))   type = MIgp;
  else if (!esl_strcmp(covtype, "RAFp"))   type = RAFp;
  else if (!esl_strcmp(covtype, "RAFSp"))  type = RAFSp;

  else if (!esl_strcmp(covtype, "CHIa"))   type = CHIa;
  else if (!esl_strcmp(covtype, "GTa"))    type = GTa;
  else if (!esl_strcmp(covtype, "OMESa"))  type = OMESa;
  else if (!esl_strcmp(covtype, "MIa"))    type = MIa;
  else if (!esl_strcmp(covtype, "MIra"))   type = MIra;
  else if (!esl_strcmp(covtype, "MIga"))   type = MIga;
  else if (!esl_strcmp(covtype, "RAFa"))   type = RAFa;
  else if (!esl_strcmp(covtype, "RAFSa"))  type = RAFSa;

  else
    ESL_XFAIL(eslFAIL, errbuf, "wrong COVTYPE %s", covtype);

  *ret_type = type;
  return eslOK;
  
 ERROR:
  return status;
}

struct mutual_s *
cov_Create(int64_t alen, int64_t nseq, int ishuffled, int nseqthresh, int alenthresh, ESL_ALPHABET *abc, COVCLASS covclass)
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
  mi->alenthresh = alenthresh;
  mi->ishuffled  = ishuffled;
  mi->abc        = abc;

  ESL_ALLOC(mi->pp,                  sizeof(double **) * alen);
  ESL_ALLOC(mi->nseff,               sizeof(int     *) * alen);
  ESL_ALLOC(mi->ngap,                sizeof(int     *) * alen);
  ESL_ALLOC(mi->pm,                  sizeof(double  *) * alen);
  for (i = 0; i < alen; i++) {
    ESL_ALLOC(mi->pp[i],             sizeof(double  *) * alen);
    ESL_ALLOC(mi->nseff[i],          sizeof(int      ) * alen);
    ESL_ALLOC(mi->ngap[i],           sizeof(int      ) * alen);
    ESL_ALLOC(mi->pm[i],             sizeof(double   ) * K);
    for (j = 0; j < alen; j++) {
       ESL_ALLOC(mi->pp[i][j],       sizeof(double  ) * K2);
    }
  }
   
  mi->COV  = esl_dmatrix_Create(alen, alen);
 
  /* initialize for adding counts */
  for (i = 0; i < alen; i++) {
    esl_vec_DSet(mi->pm[i], K, 0.0); 
 
    for (j = 0; j < alen; j++) {
      mi->nseff[i][j] = 0;
      mi->ngap[i][j]  = 0;
      esl_vec_DSet(mi->pp[i][j],       K2, 0.0); 
    }
  }

  /* inititalize to zero the COV matrix */
  cov_ReuseCOV(mi, COVNONE, covclass);
  
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
cov_NaivePP(ESL_RANDOMNESS *r, ESL_MSA *msa, struct mutual_s *mi, int donull2b, double tol, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i, j;
  int     status;

  for (i = 0; i < alen-1; i ++)
    for (j = i+1; j < alen; j ++) {
      status = mutual_naive_ppij(r, i, j, msa, mi, donull2b, tol, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
  
  return eslOK;

 ERROR:
  return status;
}

int
cov_Marginals(struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int     K = mi->abc->K;
  int     i, j;
  int     x, y;
  int     status;

  /* pm are the marginals */
  for (i = 0; i < mi->alen; i ++) {
    esl_vec_DSet(mi->pm[i], K, 0.0);

    for (j = 0; j < mi->alen; j ++)     
      for (x = 0; x < K; x ++) 
	for (y = 0; y < K; y ++) 
	  mi->pm[i][x] += mi->pp[i][j][IDX(x,y,K)];

    // it should be normalized, but just in case
    esl_vec_DNorm(mi->pm[i], K);              
  
    status = esl_vec_DValidate(mi->pm[i], K, tol, errbuf);
    if (status != eslOK) {
      printf("pm[%d]\n", i);
      esl_vec_DDump(stdout, mi->pm[i], K, "ACGU");
      goto ERROR;
    }
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
cov_SignificantPairs_Ranking(struct data_s *data, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi = data->mi;
  ESL_DMATRIX     *mtx = mi->COV;
  char            *threshtype = NULL;
  char            *covtype = NULL;
  RANKLIST        *ranklist = NULL;
  HITLIST         *hitlist = NULL;
  double           pmass, newmass;
  double           sen;
  double           ppv;
  double           F;
  double           cov;
  double           eval;
  double           val;
  double           bmax;
  double           add;
  double           tol = 1e-2;
  int              usenull = TRUE; // otherwise use ranklist->ht (the non_SS covariations)
  int              fp, tf, t, f, neg;
  int              i, j;
  int              b;
  int              status;

  if (usenull && (data->mode == GIVSS || data->mode == CYKSS) && data->ranklist_null == NULL) usenull = FALSE;

  cov_COVTYPEString(&covtype, mi->type, data->errbuf);
  cov_THRESHTYPEString(&threshtype, data->thresh->type, NULL);

  if (data->rocfp && (data->mode == GIVSS || data->mode == CYKSS)) {
    fprintf(data->rocfp, "\n# %s ", covtype);  
    if (mi->ishuffled) fprintf(data->rocfp, "shuffled thresh fp tf found true negatives sen ppv F evalue\n"); 
    else               fprintf(data->rocfp, "thresh fp tf found true negatives sen ppv F evalue\n"); 
  }
  
  bmax = mi->maxCOV+data->w;
  while (fabs(bmax-data->bmin) < tol) bmax += data->w;

  ranklist = cov_CreateRankList(bmax, data->bmin, data->w);
  for (i = 0; i < mi->alen-1; i ++) 
    for (j = i+1; j < mi->alen; j ++) {

      /* add to the ha histogram  */
      add = ESL_MAX(mtx->mx[i][j], data->bmin+data->w);
      esl_histogram_Add(ranklist->ha, add);

      /* add to the ht histogram if not a real basepair */
      if (data->ct[i+1] != j+1) 
  	esl_histogram_Add(ranklist->ht, add);
    }
  /* initialize to something impossible*/
  ranklist->scthresh = ranklist->ha->xmax + ranklist->ha->w;

  /* histogram and exponential fit */
  if (data->mode == GIVSS || data->mode == CYKSS) {
    if (data->ranklist_null && usenull) {
       if (data->ranklist_null->ha->nb < ranklist->ha->nb) {
	ESL_REALLOC(data->ranklist_null->ha->obs, sizeof(uint64_t) * ranklist->ha->nb);
	for (i = data->ranklist_null->ha->nb; i < ranklist->ha->nb; i++) data->ranklist_null->ha->obs[i] = 0;
	data->ranklist_null->ha->nb = ranklist->ha->nb;
      }
    }
    
    /* censor the histogram and do an exponential fit to the tail */
    if (!usenull) pmass = (data->Nfit < ranklist->ht->Nc)?            (double)data->Nfit/(double)data->ranklist_null->ha->Nc : data->pmass;
    else          pmass = (data->Nfit < data->ranklist_null->ha->Nc)? (double)data->Nfit/(double)data->ranklist_null->ha->Nc : data->pmass;
    if (data->doexpfit) {
      if (!usenull) status = cov_NullFitExponential(ranklist->ht,            pmass, &newmass, &data->mu, &data->lambda, data->verbose, data->errbuf);
      else          status = cov_NullFitExponential(data->ranklist_null->ha, pmass, &newmass, &data->mu, &data->lambda, data->verbose, data->errbuf);
      if (status != eslOK) goto ERROR;
      if (1||data->verbose) {
	if (!usenull) 
	  printf("ExpFIT: pmass %f mu %f lambda %f\n", newmass, data->mu, data->lambda);
	else          
	  printf("ExpFIT: pmass %f mu %f lambda %f\n", newmass, data->mu, data->lambda);
      }
    }
    else { // a gamma fit
      if (!usenull) status = cov_NullFitGamma(ranklist->ht,            pmass, &newmass, &data->mu, &data->lambda, &data->tau, data->verbose, data->errbuf);
      else          status = cov_NullFitGamma(data->ranklist_null->ha, pmass, &newmass, &data->mu, &data->lambda, &data->tau, data->verbose, data->errbuf);
      if (status != eslOK) goto ERROR;
     if (1||data->verbose) {
	if (!usenull) 
	  printf("GammaFIT: pmass %f mu %f lambda %f tau %f\n", newmass, data->mu, data->lambda, data->tau);
	else          
	  printf("GammaFIT: pmass %f mu %f lambda %f tau %f\n", newmass, data->mu, data->lambda, data->tau);
      }
    }
  }

  /* ranklist from histogram ha */
  for (b = ranklist->ha->imax; b >= ranklist->ha->imin; b --) {
    cov = esl_histogram_Bin2LBound(ranklist->ha, b);
    
    f = t = tf = 0;
    for (i = 0; i < mi->alen-1; i ++) 
      for (j = i+1; j < mi->alen; j ++) {
	if (mtx->mx[i][j] > cov)    f  ++;
	if (data->ct[i+1] == j+1) { t  ++;
	  if (mtx->mx[i][j] > cov)  tf ++;
	}	
     }
    
    fp  = f - tf;
    sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;   
    neg  = mi->alen * (mi->alen-1) / 2 - t;
    if (data->ranklist_null) {
      eval = cov2evalue(data, cov, ranklist->ha->Nc, data->ranklist_null->ha);
    }
    else {
      eval = eslINFINITY;
    }
    if (data->rocfp && (data->mode == GIVSS || data->mode == CYKSS)) 
      fprintf(data->rocfp, "%.5f %d %d %d %d %d %.2f %.2f %.2f %g\n", cov, fp, tf, f, t, neg, sen, ppv, F, eval);
  
    if (data->mode == GIVSS || data->mode == CYKSS) { 
      ranklist->eval[b] = eval; // evalues
      
      switch(data->thresh->type) {
      case Eval:    val = ranklist->eval[b]; break;
      default: 
	ESL_XFAIL(eslFAIL, data->errbuf, "do not recognize thresholding type\n");
	break;
      } 
      
      if (val > 0.0 && val <= data->thresh->val) { 
	ranklist->scthresh = cov; 
	data->thresh->sc   = cov; 
      }
    }
  }

  if (data->mode == GIVSS || data->mode == CYKSS) {
    status = cov_CreateHitList(data, mi, ranklist, &hitlist, covtype, threshtype, usenull);
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

  ranklist->ha = NULL;
  ranklist->ht = NULL;
  ranklist->ha = esl_histogram_CreateFull(bmin, bmax, w);
  ranklist->ht = esl_histogram_CreateFull(bmin, bmax, w);
  if (ranklist->ha == NULL) goto ERROR;
  if (ranklist->ht == NULL) goto ERROR;

  ranklist->eval = NULL;
  ESL_ALLOC(ranklist->eval, sizeof(double) * ranklist->ha->nb);
  esl_vec_DSet(ranklist->eval, ranklist->ha->nb, eslINFINITY);

  ranklist->scthresh = ranklist->ha->xmax;
  return ranklist;

 ERROR:
  return NULL;
}

int
cov_GrowRankList(RANKLIST **oranklist, double bmax, double bmin)
{
  RANKLIST *ranklist = *oranklist;
  RANKLIST *new = NULL;
  double    new_bmin;
  int       b, newb;
  int       status;

  /* bmin has to be a w-multiple of ranklist->bin */
  new_bmin = ranklist->ha->bmin;
  if (bmin < ranklist->ha->bmin) new_bmin -= fabs(bmin) * 2. * ranklist->ha->w;
  
  new = cov_CreateRankList(ESL_MAX(bmax, ranklist->ha->bmax), new_bmin, ranklist->ha->w);
  if (new == NULL) goto ERROR;

  new->ha->n    = ranklist->ha->n;
  new->ha->xmin = ranklist->ha->xmin;
  new->ha->xmax = ranklist->ha->xmax;
  new->ha->imin = ranklist->ha->imin;
  new->ha->imax = ranklist->ha->imax;
  new->ha->Nc   = ranklist->ha->Nc;
  new->ha->No   = ranklist->ha->No;

  new->ht->n    = ranklist->ht->n;
  new->ht->xmin = ranklist->ht->xmin;
  new->ht->xmax = ranklist->ht->xmax;
  new->ht->imin = ranklist->ht->imin;
  new->ht->imax = ranklist->ht->imax;
  new->ht->Nc   = ranklist->ht->Nc;
  new->ht->No   = ranklist->ht->No;

  for (b = ranklist->ha->imin; b <= ranklist->ha->imax; b ++) {
    cov_ranklist_Bin2Bin(b, ranklist->ha, new->ha, &newb);
    if (newb < new->ha->nb) {
      new->eval[newb]    = ranklist->eval[b];
      new->ha->obs[newb] = ranklist->ha->obs[b];
    }
  }
    
  cov_FreeRankList(ranklist);
  *oranklist = new;
  return eslOK;

 ERROR:
  if (new) cov_FreeRankList(new);
  return status;
}

int              
cov_DumpRankList(FILE *fp, RANKLIST *ranklist)
{
  int b;

  printf("imin %d imax %d covmin %f covmax %f\n", ranklist->ha->imin, ranklist->ha->imax, ranklist->ha->xmin, ranklist->ha->xmax);
  for (b = ranklist->ha->imax; b >= ranklist->ha->imin; b --) 
    printf("cov %f eval %f\n",  esl_histogram_Bin2LBound(ranklist->ha,b), ranklist->eval[b]);
 
  return eslOK;
}
int              
cov_DumpHistogram(FILE *fp, ESL_HISTOGRAM *h)
{
  double nobs = 0.;
  double cum_obs = 0.;
  double cum_exp = 0.;
  int    b;

  for (b = h->imax; b >= h->imin; b --) nobs += (double)h->obs[b];

  printf("imin %d imax %d covmin %f covmax %f nobs %f\n", h->imin, h->imax, h->xmin, h->xmax, nobs);
  for (b = h->imax; b >= h->imin; b --) {
    cum_obs += (double)h->obs[b];
    if (h->expect) cum_exp += (double)h->expect[b];
    if (h->expect) printf("b %d val %f obs %f cum_obs %f expect %f cum_exp %f\n", 
			  b, esl_histogram_Bin2LBound(h, b), (double)h->obs[b], cum_obs, h->expect[b], cum_exp); 
    else           printf("b %d val %f obs %f cum_obs %f\n", b, esl_histogram_Bin2LBound(h, b), (double)h->obs[b], cum_obs); 
  }
  
  return eslOK;
}

int 
cov_CreateHitList(struct data_s *data, struct mutual_s *mi, RANKLIST *ranklist, HITLIST **ret_hitlist,
		  char *covtype, char *threshtype, int usenull)
{
  HITLIST  *hitlist = NULL;
  double    sen, ppv, F;
  int       alloc_nhit = 5;
  int       bin;
  int       tf = 0;
  int       f  = 0;
  int       t, fp;
  int       BP, NBP;
  int       nhit;
  int       h = 0;
  int       i, j;
  int       b;
  int       status;
  
  ESL_ALLOC(hitlist, sizeof(HITLIST));
  hitlist->hit    = NULL;
  hitlist->srthit = NULL;

  nhit = alloc_nhit;
  ESL_ALLOC(hitlist->hit,    sizeof(HIT)   * nhit);
  ESL_ALLOC(hitlist->srthit, sizeof(HIT *) * nhit);
  hitlist->srthit[0] = hitlist->hit;

  hitlist->covthresh = ranklist->scthresh;

  BP  = number_pairs(mi->alen, data->ct);
  NBP = mi->alen * (mi->alen-1) / 2 - BP;
  
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      if (mi->COV->mx[i][j] >= hitlist->covthresh) {
	
	if (usenull) esl_histogram_Score2Bin(ranklist->ht,            mi->COV->mx[i][j], &bin);
	else         esl_histogram_Score2Bin(data->ranklist_null->ha, mi->COV->mx[i][j], &bin);
	
	if (h == nhit - 1) {
 	  nhit += alloc_nhit;
	  
	  ESL_REALLOC(hitlist->hit,    sizeof(HIT)   * nhit);
	  ESL_REALLOC(hitlist->srthit, sizeof(HIT *) * nhit);
	}
	/* initialize */
	hitlist->hit[h].sc            = -eslINFINITY;
	hitlist->hit[h].Eval          = +eslINFINITY;
	hitlist->hit[h].is_bpair      = FALSE;
	hitlist->hit[h].is_compatible = FALSE;
	
	for (b = ranklist->ha->imax; b >= ranklist->ha->imin; b --) {
	  if (mi->COV->mx[i][j] <= ranklist->ha->bmin+(double)(b+1)*ranklist->ha->w) {
	    hitlist->hit[h].Eval    = ranklist->eval[b];
	  }
	  else break;
	}
	
	hitlist->hit[h].i    = i;
	hitlist->hit[h].j    = j;
	hitlist->hit[h].sc   = mi->COV->mx[i][j];
	
	if (data->ct[i+1] == j+1) { hitlist->hit[h].is_bpair = TRUE;  }
	else                { 
	  hitlist->hit[h].is_bpair = FALSE; 
	  if (data->ct[i+1] == 0 && data->ct[j+1] == 0) hitlist->hit[h].is_compatible = TRUE;
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
  
  if (data->mode != RANSS && data->sumfp) {
    fprintf(data->sumfp, " %s %d %d %d %f %f ", 
	    covtype, tf, t, data->onbpairs, (t > 0)? 100.*(double)tf/(double)t:0.0, (data->onbpairs>0)? 100.*(double)tf/(double)data->onbpairs:0.0);
  }
  if (data->outfp) {
    if (data->mode == CYKSS) fprintf(data->outfp, "# cyk-cov structure\n");
    fprintf(data->outfp,    "# %s thresh %s %f cov=%f [%f,%f] [%d | %d %d %d | %f %f %f] \n", 
	    covtype, threshtype, data->thresh->val, ranklist->scthresh, ranklist->ha->xmin, ranklist->ha->xmax, fp, tf, t, f, sen, ppv, F);
    cov_WriteHitList(data->outfp,    nhit, hitlist, data->msamap, data->firstpos);
  }
  
  if (data->outsrtfp) {
    if (data->mode == CYKSS) fprintf(data->outsrtfp, "# cyk-cov structure\n");
    fprintf(data->outsrtfp, "# %s thresh %s %f cov=%f [%f,%f] [%d | %d %d %d | %f %f %f] \n", 
	    covtype, threshtype, data->thresh->val, ranklist->scthresh, ranklist->ha->xmin, ranklist->ha->xmax, fp, tf, t, f, sen, ppv, F);
    cov_WriteRankedHitList(data->outsrtfp, nhit, hitlist, data->msamap, data->firstpos);
  }
  
  if (ret_hitlist) *ret_hitlist = hitlist; else cov_FreeHitList(hitlist);
  return eslOK;
  
 ERROR:
  if (hitlist) cov_FreeHitList(hitlist);
  return status;
}

int 
cov_WriteHitList(FILE *fp, int nhit, HITLIST *hitlist, int *msamap, int firstpos)
{
  int h;
  int ih, jh;

  if (fp == NULL) return eslOK;

  for (h = 0; h < nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    
    if (hitlist->hit[h].is_bpair)      { 
      fprintf(fp, "*\t%10d\t%10d\t%.2f\t%g\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval); 
    }
    else if (hitlist->hit[h].is_compatible) { 
      fprintf(fp, "~\t%10d\t%10d\t%.2f\t%g\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval); 
    }
    else { 
      fprintf(fp, " \t%10d\t%10d\t%.2f\t%g\n",
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval); 
    }  
  }

  return eslOK;
}

static int
hit_sorted_by_eval(const void *vh1, const void *vh2)
{
  HIT *h1 = *((HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  HIT *h2 = *((HIT **) vh2);

  if      (h1->sc < h2->sc) return  1;
  else if (h1->sc > h2->sc) return -1;
  else {
 
    /* report first pair first */
    int dir1 = (h1->i < h2->i ? 1 : -1);
    int dir2 = (h1->j < h2->j ? 1 : -1);
    if (dir1 != dir2) return dir2; // so if dir1 is pos (1), and dir2 is neg (-1), this will return -1, placing h1 before h2;  otherwise, vice versa
    else              return dir1;

  }
}

int 
cov_WriteRankedHitList(FILE *fp, int nhit, HITLIST *hitlist, int *msamap, int firstpos)
{
  int h;
  int ih, jh;

  if (fp == NULL) return eslOK;

  for (h = 0; h < nhit; h++) hitlist->srthit[h] = hitlist->hit + h;
  if (nhit > 1) qsort(hitlist->srthit, nhit, sizeof(HIT *), hit_sorted_by_eval);

  for (h = 0; h < nhit; h ++) {
    ih = hitlist->srthit[h]->i;
    jh = hitlist->srthit[h]->j;
    
    if (hitlist->srthit[h]->is_bpair)      { 
      fprintf(fp, "*\t%10d\t%10d\t%.2f\t%g\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval); 
    }
    else if (hitlist->srthit[h]->is_compatible) { 
      fprintf(fp, "~\t%10d\t%10d\t%.2f\t%g\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval); 
    }
    else { 
      fprintf(fp, " \t%10d\t%10d\t%.2f\t%g\n",
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval); 
    }  
  }

  return eslOK;
}


void
cov_FreeRankList(RANKLIST *ranklist)
{
  if (ranklist == NULL) return;

  if (ranklist->ha)   esl_histogram_Destroy(ranklist->ha);
  if (ranklist->ht)   esl_histogram_Destroy(ranklist->ht);
  if (ranklist->eval) free(ranklist->eval);
  free(ranklist);
}

void
cov_FreeHitList(HITLIST *hitlist)
{
  if (hitlist == NULL) return;

  if (hitlist->srthit) free(hitlist->srthit);
  if (hitlist->hit)    free(hitlist->hit);
  free(hitlist);
}

int
cov_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int firstpos, int *ct, int verbose, char *errbuf)
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
      printf("[%d][%d] %f | %f | %f %f | %f %f\n", msamap[i]+firstpos, msamap[ipair]+firstpos, zscore, mi->COV->mx[i][ipair], avgi, stdi, avgj, stdj);
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
cov_CYKCOVCT(struct data_s *data, ESL_MSA *msa, RANKLIST **ret_ranklist, int minloop, enum grammar_e G, double covthresh)
{
  RANKLIST      *ranklist = NULL;
  HITLIST       *hitlist = NULL;
  int           *cykct = NULL;
  char          *ss = NULL;	
  SCVAL          sc;
  int            status;
            
   /* calculate the cykcov ct vector */
  status = CYKCOV(data->r, data->mi, &cykct, &sc, minloop, covthresh, data->errbuf, data->verbose);
  if (status != eslOK) goto ERROR;

  /* impose the ct on the msa GC line 'cons_ss' */
  ESL_ALLOC(ss, sizeof(char) * (msa->alen+1));
  esl_ct2simplewuss(cykct, msa->alen, ss);
  /* replace the 'SS_cons' GC line with the new ss */
  esl_sprintf(&(msa->ss_cons), "%s", ss);
  if (!msa->ax) esl_msa_Digitize(data->mi->abc, msa, data->errbuf);
  if (data->verbose) {
    printf("cykcov score = %f minloop %d covthresh %f\n", sc, minloop, covthresh);
    printf("ss:%s\n", ss);
  }

  /* expand the CT with compatible/stacked A:U C:G G:U pairs */
  status = cov_ExpandCT(data->R2Rcykfile, data->R2Rall, data->r, msa, &cykct, minloop, G, data->verbose, data->errbuf);
  if (status != eslOK) goto ERROR;
  data->ct = cykct;
  
  /* redo the hitlist since the ct has now changed */
  status = cov_SignificantPairs_Ranking(data, &ranklist, &hitlist);
  if (status != eslOK) goto ERROR;

  /* R2R */
  status = cov_R2R(data->R2Rcykfile, data->R2Rversion, data->R2Rall, msa, cykct, hitlist, TRUE, TRUE, data->verbose, data->errbuf);
  if (status != eslOK) goto ERROR;

  /* DotPlots (pdf,svg) */
  status = cov_DotPlot(data->gnuplot, data->dplotfile, msa, cykct, data->mi, data->msamap, data->firstpos, hitlist, TRUE,  data->verbose, data->errbuf);
  if (status != eslOK) goto ERROR;
  status = cov_DotPlot(data->gnuplot, data->dplotfile, msa, cykct, data->mi, data->msamap, data->firstpos, hitlist, FALSE, data->verbose, data->errbuf);
  if (status != eslOK) goto ERROR;

  *ret_ranklist = ranklist;

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
cov_WriteHistogram(struct data_s *data, char *gnuplot, char *covhisfile, char *nullcovhisfile, RANKLIST *ranklist, char *title)
{
  FILE     *fp = NULL;
  RANKLIST *ranklist_null = data->ranklist_null;
  char     *errbuf = data->errbuf;
  int       status;

  if (ranklist == NULL) return eslOK;

  if (covhisfile) {
    if ((fp = fopen(covhisfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "could not open covhisfile %s\n", covhisfile);
    esl_histogram_PlotSurvival(fp, ranklist->ht);
    fclose(fp);
  }
  if (ranklist_null && nullcovhisfile) {
    if ((fp = fopen(nullcovhisfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "could not open nullcovhisfile %s\n", nullcovhisfile);
    esl_histogram_PlotSurvival(fp, ranklist_null->ha);
    fclose(fp);
  }

  status = cov_PlotHistogramSurvival(data, gnuplot, covhisfile, ranklist, title, FALSE);
  status = cov_PlotHistogramSurvival(data, gnuplot, covhisfile, ranklist, title, TRUE);
  if (status != eslOK) goto ERROR;

  return eslOK;

 ERROR:
  return status;

}

int 
cov_NullFitExponential(ESL_HISTOGRAM *h, double pmass, double *ret_newmass, double *ret_mu, double *ret_lambda, int verbose, char *errbuf)
{
  double ep[2];  	/* estimated mu, lambda  */
  double newmass;
  int    status;

  /* set the tail by mass */
  status = esl_histogram_SetTailByMass(h, pmass, &newmass);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not set TailByMass");

  /* exponential fit to tail */
  status = esl_exp_FitCompleteBinned(h, &ep[0], &ep[1]);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not do exponential fit");

  while (ep[1] == eslINFINITY && pmass < 0.5) {
    pmass *= 1.5;
    status = esl_histogram_SetTailByMass(h, pmass, &newmass);
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not set TailByMass");
    status = esl_exp_FitCompleteBinned(h, &ep[0], &ep[1]);
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not do exponential fit");
  }

  /* add the expected data to the histogram */
  if (ep[1] < eslINFINITY) {  
    status = esl_histogram_SetExpectedTail(h, ep[0], newmass, &esl_exp_generic_cdf, ep);
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not set expected tail");
  }

  *ret_mu      = ep[0];
  *ret_lambda  = ep[1];
  *ret_newmass = newmass;
  return eslOK;

 ERROR:
  return status;
}

int 
cov_NullFitGamma(ESL_HISTOGRAM *h, double pmass, double *ret_newmass, double *ret_mu, double *ret_lambda, double *ret_k, int verbose, char *errbuf)
{
  double ep[3];  	/* ep[0] = mu; ep[1]=lambda=1/2; ep[2] = tau = k/2 */
  double newmass;
  int    status;

  /* set the tail by mass */
  status = esl_histogram_SetTailByMass(h, pmass, &newmass);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not set TailByMass");

  /* chi-square fit to tail */
  status = esl_gam_FitCompleteBinned(h, &ep[0], &ep[1], &ep[2]);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not do chi-square fit");
  if (verbose) printf("GammaFit mu = %f lambda = %f tau = %f\n", ep[0], ep[1], ep[2]);

  /* add the expected data to the histogram */
  if (ep[1] < eslINFINITY) {  
    status = esl_histogram_SetExpectedTail(h, ep[0], newmass, &esl_gam_generic_cdf, ep);
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not set expected tail");
  }

  *ret_mu     = ep[0];
  *ret_lambda = ep[1];
  *ret_k      = ep[2];
  *ret_newmass = newmass;
  return eslOK;

 ERROR:
  return status;
}


int 
cov_PlotHistogramSurvival(struct data_s *data, char *gnuplot, char *covhisfile, RANKLIST *ranklist, char *title, int dosvg)
{
  FILE     *pipe;
  RANKLIST *ranklist_null = data->ranklist_null;
  RANKLIST *ranklist_aux  = data->ranklist_null;
  char     *filename = NULL;
  char     *outplot = NULL;
  char     *key1 = NULL;
  char     *key2 = NULL;
  char     *key3 = NULL;
  char     *key4 = NULL;
  double    minphi;
  double    minmass = 0.005;
  int       pointype;
  double    pointintbox;
  int       linew;
  double    pointsize;
  double    xmin, xmax;
  double    ymin, ymax;
  double    posx, posy;
  double    incx, incy;
  double    expsurv;
  double    offx, offy;
  int       subsample;
  int       linespoints;
  int       i;
  int       status;

  if (gnuplot    == NULL) return eslOK;
  if (covhisfile == NULL) return eslOK;
  if (ranklist   == NULL) return eslOK;

  esl_FileTail(covhisfile, FALSE, &filename);

  esl_sprintf(&key1, "all pairs");
  esl_sprintf(&key2, "not proposed pairs");
  esl_sprintf(&key3, "null distribution");
  if (ranklist_aux) esl_sprintf(&key4, "null-null distribution");

  /* for plotting */
  esl_histogram_SetTailByMass(ranklist->ha, minmass, NULL);
  minphi = ranklist->ha->phi; 
 
  /* do the plotting */
  pipe = popen(gnuplot, "w");
  if (dosvg) {
    esl_sprintf(&outplot, "%s.svg", covhisfile);
    fprintf(pipe, "set terminal svg dynamic fname 'Arial' fsize 12 \n");
    pointype    = 71;
    pointsize   = 0.6;
    pointintbox = 0.4;
    linew       = 1;
  }
  else {
    esl_sprintf(&outplot, "%s.ps", covhisfile);
    fprintf(pipe, "set terminal postscript color 14\n");
    pointype    = 65;
    pointsize   = 0.7;
    pointintbox = 1.0;
    linew       = 2;
  }

  fprintf(pipe, "set output '%s'\n", outplot);
  fprintf(pipe, "set title '%s' \n", title);

  fprintf(pipe, "set style line 1   lt 1 lc rgb 'grey'      pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 2   lt 1 lc rgb 'brown'     pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 3   lt 1 lc rgb 'cyan'      pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 4   lt 1 lc rgb 'red'       pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 5   lt 1 lc rgb 'orange'    pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 6   lt 1 lc rgb 'turquoise' pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 7   lt 1 lc rgb 'black'     pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 8   lt 1 lc rgb 'green'     pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 9   lt 1 lc rgb 'blue'      pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 11  lt 1 lc rgb 'grey'      pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 22  lt 1 lc rgb 'brown'     pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 33  lt 1 lc rgb 'cyan'      pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 44  lt 1 lc rgb 'red'       pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 55  lt 1 lc rgb 'orange'    pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 66  lt 1 lc rgb 'turquoise' pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 77  lt 1 lc rgb 'black'     pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 88  lt 1 lc rgb 'green'     pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 99  lt 1 lc rgb 'blue'      pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);

  // plot evalue
  fprintf(pipe, "set multiplot\n");  
  fprintf(pipe, "set ylabel 'Expected or Observed #pairs(x > score)'\n");
  ymax = 100.;
  ymin = 0.1*ranklist->ha->Nc/ranklist_null->ha->Nc;
  incy = (ymax-ymin)/26.;
  posy = ymax - 8*incy;
  fprintf(pipe, "set logscale y\n");
  fprintf(pipe, "set yrange [%g:%f]\n", ymin, ymax);

  fprintf(pipe, "set xlabel 'covariation score'\n");
  xmin = (ranklist_null)? ESL_MIN(ranklist->ht->phi, ESL_MIN(minphi, ranklist_null->ha->phi)) : ESL_MIN(minphi, ranklist->ht->phi);
  xmax = (ranklist_null)? ESL_MAX(ranklist->ha->xmax,ranklist_null->ha->xmax) : ranklist->ha->xmax;
  xmin = ESL_MIN(evalue2cov(data, ymax, ranklist->ha->Nc, ranklist->ha), evalue2cov(data, ymax, ranklist->ha->Nc, ranklist_null->ha));
  incx = (xmax-xmin)/12.;
  xmax += incx;
  posx = xmin + 11.*incx;
  fprintf(pipe, "set xrange [%f:%f]\n", xmin, xmax);

  subsample = (int)(0.5*ranklist_null->ha->Nc/ranklist->ha->Nc);
  if (subsample < 1) subsample = 1;
  subsample = 1;

  offx = incx * 1/2;
  offy = incy * 16;
  expsurv = 0.00001;
  cov_plot_lineatexpcov(pipe, data, expsurv, ranklist->ha->Nc, ranklist_null->ha, ymin, ymax, "E 1e-5", offx, offy, 1);
  
  expsurv = 1.0;
  offy = incy * 11;
  cov_plot_lineatexpcov(pipe, data, expsurv, ranklist->ha->Nc, ranklist_null->ha, ymin, ymax, "E 1.00", offx, offy, 1);
  
  expsurv = 10.0;
  offy = incy * 6;
  cov_plot_lineatexpcov(pipe, data, expsurv, ranklist->ha->Nc, ranklist_null->ha, ymin, ymax, "E 10.0", offx, offy, 1);

  if (ranklist_null) {
    linespoints = FALSE;
    status = cov_histogram_plotexpectsurv(pipe, ranklist->ha->Nc, ranklist_null->ha, key3, posx, posy-12.*incy, FALSE, subsample, linespoints, 77, 7);
    if (status != eslOK) goto ERROR;
  }
  if (ranklist_aux) {
    linespoints = FALSE;
    status = cov_histogram_plotexpectsurv(pipe, ranklist->ha->Nc, ranklist_aux->ha,  key4, posx, posy-16.*incy, FALSE, subsample, linespoints, 11, 7);
    if (status != eslOK) goto ERROR;
  }
  linespoints = TRUE; 
#if 0
  status = cov_histogram_plotexpectsurv  (pipe, ranklist->ha->Nc, ranklist->ht, key2, posx, posy-8*incy,        FALSE, 1, linespoints, 44, 2);
  if (status != eslOK) goto ERROR;
#endif
  status = cov_histogram_plotexpectsurv  (pipe, ranklist->ha->Nc, ranklist->ha, key1, posx, posy,               FALSE, 1, linespoints, 99, 2);
  if (status != eslOK) goto ERROR;
  
  if (!dosvg) {
    
    // log survival plot for ranklist and ranklist_null
    fprintf(pipe, "unset logscale y\n");
    fprintf(pipe, "set multiplot\n");
    fprintf(pipe, "set xlabel 'covariation score'\n");
    xmin = (ranklist_null)? ESL_MIN(ranklist->ht->phi,ESL_MIN(minphi,ranklist_null->ha->phi)) : ESL_MIN(minphi,ranklist->ht->phi);
    xmax = (ranklist_null)? ESL_MAX(ranklist->ha->xmax,ranklist_null->ha->xmax) : ranklist->ha->xmax;
    incx = (xmax-xmin)/12.;
    xmax += incx;
    posx = xmin + 6.*incx;
    fprintf(pipe, "set xrange [%f:%f]\n", xmin, xmax);
    
    ymin = -log(ranklist->ht->Nc);
    ymax = log(ESL_MAX(minmass, data->pmass));
    incy = (ymax-ymin)/12.;
    ymin -= incy;
    posy = ymax - incy;
    fprintf(pipe, "set yrange [%f:%f]\n", ymin, ymax);
    fprintf(pipe, "set ylabel 'lnP(x > score)'\n");
    if (ranklist_null) {
      status = cov_histogram_plotsurvival(pipe, ranklist_null->ha, key3, posx, posy-2.*incy, TRUE, subsample, 77, 7);
      if (status != eslOK) goto ERROR;
    }
    if (ranklist_aux) {
      status = cov_histogram_plotsurvival(pipe, ranklist_aux->ha, key4, posx, posy-3.*incy, TRUE, subsample, 11, 7);
      if (status != eslOK) goto ERROR;
    }

    status = cov_histogram_plotsurvival  (pipe, ranklist->ht, key2, posx, posy-incy, TRUE, 1, 44, 2);
    if (status != eslOK) goto ERROR;
    status = cov_histogram_plotsurvival  (pipe, ranklist->ha, key1, posx, posy,      TRUE, 1, 99, 2);
    if (status != eslOK) goto ERROR;

    // survival plot for ranklist and ranklist_null
    fprintf(pipe, "set multiplot\n");
    xmin = (ranklist_null)? ESL_MIN(ranklist_null->ha->xmin, ranklist->ha->xmin) : ranklist->ha->xmin;
    xmax = (ranklist_null)? ESL_MAX(ranklist_null->ha->xmax, ranklist->ha->xmax) : ranklist->ha->xmax;
    ymin = 0.0;
    ymax = 1.0;
    incx = (xmax-xmin)/12.;
    incy = (ymax-ymin)/12.;
    xmin -= incx;
    posx = xmin + 4.*incx;
    posy = ymax - incy;
    fprintf(pipe, "set xrange [%f:%f]\n", xmin, xmax);
    fprintf(pipe, "set yrange [%f:%f]\n", ymin, ymax);
    fprintf(pipe, "set ylabel 'P(x > score)'\n");

    status = cov_histogram_plotsurvival  (pipe, ranklist->ha, key1, posx, posy,      FALSE, 1, 9, 2);
    if (status != eslOK) goto ERROR;
    status = cov_histogram_plotsurvival  (pipe, ranklist->ht, key2, posx, posy-incy, FALSE, 1, 4, 2);
    if (status != eslOK) goto ERROR;
    if (ranklist_null) {
      status = cov_histogram_plotsurvival(pipe, ranklist_null->ha, key3, posx, posy-2.*incy, FALSE, subsample, 7, 7);
      if (status != eslOK) goto ERROR;
    }
    if (ranklist_aux) {
      status = cov_histogram_plotsurvival(pipe, ranklist_null->ha, key4, posx, posy-3.*incy, FALSE, subsample, 1, 7);
      if (status != eslOK) goto ERROR;
    }

    // plot the density distribution
    fprintf(pipe, "set multiplot\n");
    ymin = 0.0;
    ymax = 0.0;
    for (i = ranklist->ha->imin; i < ranklist->ha->imax; i ++) 
      if ((double)ranklist->ha->obs[i]/(double)ranklist->ha->Nc > ymax) ymax = (double)ranklist->ha->obs[i]/(double)ranklist->ha->Nc;

    fprintf(pipe, "set yrange [%f:%f]\n", ymin, ymax);
    fprintf(pipe, "set ylabel 'P(score)'\n");

    status = cov_histogram_plotdensity  (pipe, ranklist->ha, key1, posx, posy,      FALSE, 1, 9, 2);
    if (status != eslOK) goto ERROR;
    status = cov_histogram_plotdensity  (pipe, ranklist->ht, key2, posx, posy-incy, FALSE, 1, 4, 2);
    if (status != eslOK) goto ERROR;
    if (ranklist_null) {
      status = cov_histogram_plotdensity(pipe, ranklist_null->ha, key3, posx, posy-2.*incy, FALSE, subsample, 7, 7);
      if (status != eslOK) goto ERROR;
    }
    if (ranklist_aux) {
      status = cov_histogram_plotdensity(pipe, ranklist_null->ha, key4, posx, posy-3.*incy, FALSE, subsample, 1, 7);
      if (status != eslOK) goto ERROR;
    }

  }
  
  pclose(pipe);
  
  free(key1);
  free(key2);
  free(key3);
  free(key4);
  free(outplot);
  free(filename);
  return eslOK;

 ERROR:
  if (key1) free(key1);
  if (key2) free(key2);
  if (key3) free(key3);
  if (key4) free(key4);
  if (outplot) free(outplot);
  if (filename) free(filename);
  return status;
}




int              
cov_DotPlot(char *gnuplot, char *dplotfile, ESL_MSA *msa, int *ct, struct mutual_s *mi, int *msamap, int firstpos, HITLIST *hitlist, 
	    int dosvg, int verbose, char *errbuf)
{
  FILE    *pipe;
  char    *filename = NULL;
  char    *outplot = NULL;
  double   pointsize;
  double   ps_max = 0.40;
  double   ps_min = 0.0003;
  int      ileft, iright;
  int      h;           /* index for hitlist */
  int      i, ipair;
  int      ih, jh;
  int      status;

  if (dplotfile == NULL) return eslOK;
  if (hitlist == NULL) return eslOK;
  
  esl_FileTail(dplotfile, FALSE, &filename);

  pipe = popen(gnuplot, "w");
  
  if (dosvg) {
    ps_max = 1.00;
    ps_min = 0.3;
    esl_sprintf(&outplot, "%s.svg", dplotfile);
    fprintf(pipe, "set terminal svg fname 'Arial' fsize 12 \n");
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

  ileft  = msa->alen;
  iright = 1;
  for (i = 1; i <= msa->alen; i ++) {
    if (ct[i] > 0) {
      if (ct[i] > i && i < ileft)  ileft  = i;
      if (ct[i] < i && i > iright) iright = i;
    }
  }
  if (ileft == msa->alen && iright == 1) {//no structure
    ileft = 1;
    iright = msa->alen;
  }
  if (ileft > iright)     ESL_XFAIL(eslFAIL, errbuf, "error in cov_DotPlot()");
  if (iright > msa->alen) ESL_XFAIL(eslFAIL, errbuf, "error in cov_DotPlot()");

  ileft  = (hitlist->nhit > 0)? ESL_MIN(ileft,  hitlist->hit[0].i+1) : ileft;
  iright = (hitlist->nhit > 0)? ESL_MAX(iright, hitlist->hit[hitlist->nhit-1].j+1) : iright;

  fprintf(pipe, "set yrange [%d:%d]\n", msamap[ileft-1]+firstpos, msamap[iright-1]+firstpos);
  fprintf(pipe, "set xrange [%d:%d]\n", msamap[ileft-1]+firstpos, msamap[iright-1]+firstpos);

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
      fprintf(pipe, "%d %d %f\n", msamap[i-1]+firstpos,     msamap[ipair-1]+firstpos, 
	      (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
      fprintf(pipe, "%d %d %f\n", msamap[ipair-1]+firstpos, msamap[i-1]+firstpos,     
	      (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
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
      fprintf(pipe, "%d %d %f\n", msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc);
      fprintf(pipe, "%d %d %f\n", msamap[jh]+firstpos, msamap[ih]+firstpos, hitlist->hit[h].sc);
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
      fprintf(pipe, "%d %d %f\n", msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc);	
      fprintf(pipe, "%d %d %f\n", msamap[jh]+firstpos, msamap[ih]+firstpos, hitlist->hit[h].sc);	
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
      fprintf(pipe, "%d %d %f\n", msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc);	
      fprintf(pipe, "%d %d %f\n", msamap[jh]+firstpos, msamap[ih]+firstpos, hitlist->hit[h].sc);	
    }
  } 
  fprintf(pipe, "e\n");
  
  pclose(pipe);
  
  free(outplot);
  free(filename);
  return eslOK;

 ERROR:
  if (outplot) free(outplot);
  if (filename) free(filename);
  return status;
}

int
cov_R2R(char *r2rfile, char *r2rversion, int r2rall, ESL_MSA *msa, int *ct,
	HITLIST *hitlist, int makepdf, int makesvg, int verbose, char *errbuf)
 {
  ESLX_MSAFILE *afp = NULL;
  FILE         *fp = NULL;
  ESL_MSA      *r2rmsa = NULL;
  char          tmpinfile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char          covtag[12] = "cov_SS_cons";
  char         *args = NULL;
  char         *s = NULL;
  char         *ssstr = NULL;
  char         *covstr = NULL;
  char         *prv_covstr = NULL;
  char         *tok;
  int           found;
  int           i;
  int           h;
  int           ih, jh;
  int           tagidx;
  int           idx;
  int           do_r2rcovmarkup = FALSE;
  int           status;
 
  if (r2rfile == NULL) return eslOK;
  if (hitlist == NULL) return eslOK;

  /* first modify the ss to a simple <> format. R2R cannot deal with fullwuss 
   */
  ESL_ALLOC(ssstr, sizeof(char) * (msa->alen+1));
  esl_ct2simplewuss(ct, msa->alen, ssstr);

  /* replace the 'SS_cons' GC line with the new ss */
  if (msa->ss_cons) free(msa->ss_cons); msa->ss_cons = NULL;
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

  /* modify the cov_cons_ss line according to our hitlist */
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
    else        esl_sprintf(&covstr, "%s%s", prv_covstr, tok);

    free(prv_covstr); prv_covstr = NULL;
    esl_sprintf(&prv_covstr, "%s", covstr);
    free(tok); tok = NULL;
    free(covstr); covstr = NULL;
  }
  
  /* add line #=GF R2R keep allpairs 
   * so that it does not truncate ss.
   * cannot use the standard esl_msa_addGF:
   *             esl_msa_AddGF(msa, "R2R", -1, " keep allpairs", -1);
   * since it does not parse with r2r
   *
   * turns out the above solution can only deal with the  <> annotation
   */
  if (r2rall) {
    for (tagidx = 0; tagidx < r2rmsa->ngf; tagidx++) {
      esl_strchop(r2rmsa->gf[tagidx], -1);
      if (strcmp(r2rmsa->gf[tagidx], "keep all") == 0) break;
    }
    
    if (tagidx < r2rmsa->ngf) { //remove 
      for (idx = tagidx; idx < r2rmsa->ngf-1; idx++) {
	esl_sprintf(&r2rmsa->gf_tag[idx], r2rmsa->gf_tag[idx+1]);
	esl_sprintf(&r2rmsa->gf[idx],     r2rmsa->gf[idx+1]);
      }
      r2rmsa->ngf --;
      }
    
    esl_msa_AddGF(r2rmsa, "R2R keep all", -1, "", -1);
  }
  else {
    for (tagidx = 0; tagidx < r2rmsa->ngf; tagidx++) {
      esl_strchop(r2rmsa->gf[tagidx], -1);
      if (strcmp(r2rmsa->gf[tagidx], "keep allpairs") == 0) break;
    }
    if (tagidx < r2rmsa->ngf) { //remove 
      for (idx = tagidx; idx < r2rmsa->ngf-1; idx++) {
	free(r2rmsa->gf_tag[idx]); r2rmsa->gf_tag[idx] = NULL; 
	free(r2rmsa->gf[idx]);      r2rmsa->gf[idx] = NULL; 
	esl_sprintf(&r2rmsa->gf_tag[idx], r2rmsa->gf_tag[idx+1]);
	esl_sprintf(&r2rmsa->gf[idx],     r2rmsa->gf[idx+1]);
      }
      r2rmsa->ngf --;
    }
    esl_msa_AddGF(r2rmsa, "R2R keep allpairs", -1, "", -1);
  }
  
  /* replace the r2r 'cov_SS_cons' GC line with our own */
  if (!do_r2rcovmarkup) {
    for (tagidx = 0; tagidx < r2rmsa->ngc; tagidx++)
      if (strcmp(r2rmsa->gc_tag[tagidx], covtag) == 0) break;
    if (tagidx == r2rmsa->ngc) {
      ESL_REALLOC(r2rmsa->gc_tag, (r2rmsa->ngc+1) * sizeof(char **));
      ESL_REALLOC(r2rmsa->gc,     (r2rmsa->ngc+1) * sizeof(char **));
      r2rmsa->gc[r2rmsa->ngc] = NULL;
      r2rmsa->ngc++;
    }
    if (r2rmsa->gc_tag[tagidx]) free(r2rmsa->gc_tag[tagidx]); r2rmsa->gc_tag[tagidx] = NULL;
    if ((status = esl_strdup(covtag, -1, &(r2rmsa->gc_tag[tagidx]))) != eslOK) goto ERROR;
    free(r2rmsa->gc[tagidx]); r2rmsa->gc[tagidx] = NULL; 
    esl_sprintf(&(r2rmsa->gc[tagidx]), "%s", prv_covstr);
    if (verbose) eslx_msafile_Write(stdout, r2rmsa, eslMSAFILE_PFAM);
  }
  
  /* write the R2R annotated to PFAM format */
  if ((fp = fopen(r2rfile, "w")) == NULL) esl_fatal("Failed to open r2rfile %s", r2rfile);
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
  
  esl_msa_Destroy(r2rmsa);

  remove(tmpinfile);
  remove(tmpoutfile);
  
  free(tok);
  if (args) free(args);
  if (ssstr) free(ssstr);
  if (covstr) free(covstr);
  if (prv_covstr) free(prv_covstr);
  return eslOK;

 ERROR:
  remove(tmpinfile);
  remove(tmpoutfile);
  
  if (msa)    esl_msa_Destroy(msa);
  if (r2rmsa) esl_msa_Destroy(r2rmsa);
  if (tok)    free(tok);
  if (args)   free(args);
  if (ssstr)  free(ssstr);
  if (covstr) free(covstr);
  if (prv_covstr) free(prv_covstr);
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
  
  if (r2rfile == NULL) return eslOK;

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

  if ((fp = fopen(r2rfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to open r2rfile %s", r2rfile);
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
  char    *rfline = NULL;
  ESL_SQ  *sq = NULL;
  int     *ct = *ret_ct;
  int     *cct = NULL;
  SCVAL    sc;
  float    idthresh = 0.3;
  int      status;
  
  /* create an RF sequence */
  ESL_ALLOC(rfline, sizeof(char) * (msa->alen+1));
  esl_msa_ReasonableRF(msa, idthresh, TRUE, rfline);
  sq = esl_sq_CreateFrom(msa->name, rfline, msa->desc, msa->acc, msa->ss_cons); 
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
cov_ranklist_Bin2Bin(int b, ESL_HISTOGRAM *h, ESL_HISTOGRAM *new, int *ret_newb)
{
  double x;
  int    newb;
  int    status;
  
  x = esl_histogram_Bin2UBound(h, b);
  if (! isfinite(x)) ESL_XEXCEPTION(eslERANGE, "value added to histogram is not finite");

  x = round( ((x - new->bmin) / new->w) - 1.); 

  /* x is now the bin number as a double, which we will convert to
   * int. Because x is a double (64-bit), we know all ints are exactly
   * represented.  Check for under/overflow before conversion.
   */
  if (x < (double) INT_MIN || x > (double) INT_MAX) 
    ESL_XEXCEPTION(eslERANGE, "value %f isn't going to fit in histogram", x);
  
  newb = (int) x;
  if (newb > new->nb) 
    ESL_XEXCEPTION(eslERANGE, "bin value %d isn't going to fit in histogram bin %d ", newb, new->nb);
  
  *ret_newb = newb;
  return eslOK;

 ERROR:
  *ret_newb = -1;
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
mutual_naive_ppij(ESL_RANDOMNESS *r, int i, int j, ESL_MSA *msa, struct mutual_s *mi,
		  int donull2b, double tol, int verbose, char *errbuf)
{
  int    *coli = NULL;
  int    *colj = NULL;
  int    *shcoli = NULL;
  int    *shcolj = NULL;
  int     K = mi->abc->K;
  int     K2 = K*K;
  int     s;
  int     resi, resj;
  int     x, y;
  int     status;

  esl_vec_DSet(mi->pp[i][j], K2, 1e-5); //some prior to avoid zeros
  mi->nseff[i][j] = 0;

  ESL_ALLOC(coli, sizeof(int)*msa->nseq);
  ESL_ALLOC(colj, sizeof(int)*msa->nseq);
  for (s = 0; s < msa->nseq; s ++) {
    coli[s] = msa->ax[s][i+1];
    colj[s] = msa->ax[s][j+1];
  }

  /* shuffle in each column residues that appear to be canonical (i,j) pairs */
  if (donull2b) {
    status = shuffle_null2b_col(r, msa->abc, msa->nseq, coli, colj, &shcoli, errbuf);
    if (status != eslOK) goto ERROR;
    status = shuffle_null2b_col(r, msa->abc, msa->nseq, colj, coli, &shcolj, errbuf);
    if (status != eslOK) goto ERROR;
  }

  for (s = 0; s < msa->nseq; s ++) {
    resi = (donull2b)? shcoli[s] : coli[s];
    resj = (donull2b)? shcolj[s] : colj[s];
    
    if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) { 
      mi->nseff[i][j] ++; 
      mi->pp[i][j][IDX(resi,resj,K)] += msa->wgt[s]; 
    }
    else if (esl_abc_XIsCanonical(msa->abc, resi)) { 
      mi->nseff[i][j] ++; 
      mi->ngap[i][j]  ++; 
      for (y = 0; y < K; y ++) mi->pp[i][j][IDX(resi,y,K)] += msa->wgt[s]/(double)K; 
    }
    else if (esl_abc_XIsCanonical(msa->abc, resj)) { 
      mi->nseff[i][j] ++; 
      mi->ngap[i][j]  ++; 
      for (x = 0; x < K; x ++) mi->pp[i][j][IDX(x,resj,K)] += msa->wgt[s]/(double)K; 
    }
#if 0
    else { 
      mi->nseff[i][j] ++; 
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) 
	  mi->pp[i][j][IDX(x,y,K)] += msa->wgt[s]/(double)K2;
    }
#endif
  }

  // normalize
  esl_vec_DNorm(mi->pp[i][j], K2);             

  /* symmetrize */
  for (x = 0; x < K; x ++)
    for (y = 0; y < K; y ++) 
      mi->pp[j][i][IDX(y,x,K)] = mi->pp[i][j][IDX(x,y,K)];
  mi->nseff[j][i] = mi->nseff[i][j];
  
  free(coli);
  free(colj);
  if (shcoli) free(shcoli);
  if (shcolj) free(shcolj);
  return eslOK;

 ERROR:
  if (coli)   free(coli);
  if (colj)   free(colj);
  if (shcoli) free(shcoli);
  if (shcolj) free(shcolj);

  return status;
}

int 
shuffle_null2b_col(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, int nseq, int *col, int *paircol, int **ret_shcol, char *errbuf)
{
  int *useme = NULL;
  int *shcol = NULL;
  int  s;
  int  status;
  
  /* allocation */
  ESL_ALLOC(useme, sizeof(int) * nseq); // vecto to mark residues in the column
 
  /* shuffle only positions with residues and with canonical pair in the other column */
  esl_vec_ISet(useme, nseq, FALSE);
  for (s = 0; s < nseq; s ++) 
    if ( esl_abc_XIsResidue(abc,col[s]) && is_wc(col[s], paircol[s]) ) useme[s] = TRUE;
 
  /* within colum permutation */
  status = shuffle_col(r, nseq, useme, col, &shcol, errbuf);
  if (status != eslOK) goto ERROR;

  *ret_shcol = shcol;

  free(useme);
  return eslOK;
  
 ERROR:
  if (shcol)  free(shcol);
  if (useme)  free(useme);
  return status;
}

int 
shuffle_col(ESL_RANDOMNESS *r, int nseq, int *useme, int *col, int **ret_shcol, char *errbuf)
{
  int *shcol  = NULL;
  int *seqidx = NULL;
  int *perm   = NULL;
  int  nuse;
  int  s;
  int  u;
  int  status;
 
  /* allocation */
  ESL_ALLOC(shcol, sizeof(int) * nseq);
  esl_vec_ICopy(col, nseq, shcol);
 
 /* within colum permutation */
  nuse = nseq;
  for (s = 0; s < nseq; s ++) if (useme[s] == FALSE) nuse --;
  if (nuse == 0) {
    *ret_shcol = shcol;
    return eslOK;
  }

  ESL_ALLOC(seqidx, sizeof(int) * nuse);
  u = 0;
  for (s = 0; s < nseq; s++) if (useme[s] == TRUE) { seqidx[u] = s; u++; }
  ESL_ALLOC(perm,  sizeof(int) * nuse);
  for (u = 0; u < nuse; u ++) perm[u] = u;
  if ((status = esl_vec_IShuffle(r, perm, nuse)) != eslOK) ESL_XFAIL(status, errbuf, "failed to randomize perm");

  u = 0;
  for (s = 0; s < nseq; s++) {
    if (useme[s] == TRUE) {
      shcol[s] = col[seqidx[perm[u]]];
      u ++;
    }
  }

  *ret_shcol = shcol;

  free(perm);
  free(seqidx);
  return eslOK;

 ERROR:
  if (shcol)  free(shcol);
  if (perm)   free(perm);
  if (seqidx) free(seqidx);
  return status;
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

static double
cov2evalue(struct data_s *data, double cov, int Nc, ESL_HISTOGRAM *h)
{
  double eval;
  int    c = 0;
  int    icov;
  int    i;
  
  esl_histogram_Score2Bin(h, cov, &icov);

  /* use the sampled distribution if possible */
  if (icov <= h->imax) {

    if (icov < h->imin) icov = h->imin;
    if (icov == h->imax && h->obs[h->imax] > 1) eval = (double)Nc / (double)h->Nc;
    
    for (i = h->imax; i >= icov; i--) c += h->obs[i];
    eval = (double)c * (double)Nc / (double)h->Nc;
    return eval;
  }

  /* otherwise use the fit */
  if (data->doexpfit) {
    eval = (data->lambda < eslINFINITY)? data->pmass * esl_exp_surv(cov, data->mu, data->lambda) * (double)Nc : eslINFINITY;
  }
  else {
    eval = (data->lambda < eslINFINITY)? data->pmass * esl_gam_surv(cov, data->mu, data->lambda, data->tau) * (double)Nc : eslINFINITY;
  }
  
  return eval;
}
 
static double
evalue2cov(struct data_s *data, double eval, int Nc, ESL_HISTOGRAM *h)
{
  double cov;
  double p  = eval/(double)Nc/data->pmass;
  int    c = 0;
  int    i;
  
  /* use the sampled distribution if possible */
  if (eval >= (double)Nc / (double)h->Nc) {

    for (i = h->imax; i >= h->imin; i--) {
      c += h->obs[i];
      if ((double)c * (double)Nc / (double)h->Nc > eval) break;
    }    
    cov = esl_histogram_Bin2UBound(h, i+1); 

    return cov;
  }

  /* otherwise use the fit */
  if (data->doexpfit) cov = esl_exp_invsurv(p, data->mu, data->lambda);
  else                cov = -eslINFINITY;
  
  return cov;
}


static int
cov_histogram_plotdensity(FILE *pipe, ESL_HISTOGRAM *h, char *key, double posx, double posy, int logval, int subsample, int style1, int style2)
{
  int       i;
  uint64_t  obs;
  double    exp;
  double    ai;
 
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %f,%f '%s' center tc ls %d\n", posx, posy, key, style1);
  fprintf(pipe, "plot '-' using 1:2 with linespoints ls %d \n", style1);

  if (h->obs[h->imax] > 1) 
    if (fprintf(pipe, "%f\t%f\n", 
		h->xmax, (logval)? -log((double)h->Nc) : 1.0/(double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
  for (i = h->imax; i >= h->imin; i--)
    {
      obs = h->obs[i];

      if (obs > 0 && (i-h->imin)%subsample == 0) {
	ai = esl_histogram_Bin2LBound(h, i);
	if (fprintf(pipe, "%f\t%f\n", 
		    ai, (logval)? log((double)obs)-log((double)h->Nc) : (double)obs/(double) h->Nc) < 0) 
	  ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  fprintf(pipe, "e\n");
  
  /* The expected binned counts:
   */
  if (h->expect != NULL) 
    {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");
      fprintf(pipe, "set key off\n");
      fprintf(pipe, "plot '-' using 1:2 with lines ls %d \n", style2);
      
      for (i = h->nb-1; i >= 0; i--)
	{
	  exp = h->expect[i];        /* some worry about 1+eps=1 problem here */

	  if (exp > 0. && (i-h->imin)%subsample == 0) { 
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(pipe, "%f\t%f\n", 
			ai, (logval)? log(exp)-log((double)h->Nc) : exp/(double) h->Nc) < 0) 
	      ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      fprintf(pipe, "e\n"); 
    }
  
  return eslOK;
}
static int
cov_histogram_plotsurvival(FILE *pipe, ESL_HISTOGRAM *h, char *key, double posx, double posy, int logval, int subsample, int style1, int style2)
{
  int       i;
  uint64_t  c = 0;
  double    esum;
  double    ai;
 
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %f,%f '%s' center tc ls %d\n", posx, posy, key, style1);
  fprintf(pipe, "plot '-' using 1:2 with linespoints ls %d \n", style1);

  if (h->obs[h->imax] > 1) 
    if (fprintf(pipe, "%f\t%f\n", 
		h->xmax, (logval)? -log((double)h->Nc) : 1.0/(double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
  for (i = h->imax; i >= h->imin; i--)
    {
      c += h->obs[i];

      if (h->obs[i] > 0 && (i-h->imin)%subsample == 0) {
	ai = esl_histogram_Bin2LBound(h, i);
	if (fprintf(pipe, "%f\t%f\n", 
		    ai, (logval)? log((double)c)-log((double)h->Nc) : (double)c/(double) h->Nc) < 0) 
	  ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  fprintf(pipe, "e\n");
  
  /* The expected binned counts:
   */
  if (h->expect != NULL) 
    {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");
      fprintf(pipe, "set key off\n");
      fprintf(pipe, "plot '-' using 1:2 with lines ls %d \n", style2);
      
      esum = 0.;
      for (i = h->nb-1; i >= 0; i--)
	{
	  esum += h->expect[i];        /* some worry about 1+eps=1 problem here */

	  if (h->expect[i] > 0. && (i-h->imin)%subsample == 0) { 
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(pipe, "%f\t%f\n", 
			ai, (logval)? log(esum)-log((double)h->Nc) : esum/(double) h->Nc) < 0) 
	      ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      fprintf(pipe, "e\n"); 
    }
  
  return eslOK;
}

static int
cov_histogram_plotexpectsurv(FILE *pipe, int Nc, ESL_HISTOGRAM *h, char *key, double posx, double posy, int logval, 
			     int subsample, int linespoints, int style1, int style2)
{
  int       i;
  uint64_t  c = 0;
  double    esum;
  double    ai;
 
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %f,%f '%s' center tc ls %d\n", posx, posy, key, style1);
  if (linespoints) fprintf(pipe, "plot '-' using 1:2 with linespoints ls %d \n", style1);
  else             fprintf(pipe, "plot '-' using 1:2 with points ls %d \n", style1);

  if (h->obs[h->imax] > 1) 
    if (fprintf(pipe, "%f\t%f\n", 
		h->xmax, (logval)? log((double)Nc) - log((double)h->Nc) : (double)Nc / (double)h->Nc) < 0) 
      ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");

  for (i = h->imax; i >= h->imin; i--)
    {
      c += h->obs[i];
      if (h->obs[i] > 0 && (i-h->imin)%subsample == 0) {
	ai = esl_histogram_Bin2LBound(h, i);
	if (fprintf(pipe, "%f\t%f\n", 
		    ai, (logval)? log((double)c) + log((double)Nc) - log((double)h->Nc) : (double)c * (double)Nc / (double)h->Nc) < 0) 
	  ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  fprintf(pipe, "e\n");
  
  /* The expected binned counts:
   */
  if (h->expect != NULL) 
    {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");
      fprintf(pipe, "set key off\n");
      fprintf(pipe, "plot '-' using 1:2 with lines ls %d \n", style2);
      
      esum = 0.;
      for (i = h->nb-1; i >= 0; i--)
	{
	    esum += h->expect[i];        /* some worry about 1+eps=1 problem here */

	  if (h->expect[i] > 0. && (i-h->imin)%subsample == 0) { 
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(pipe, "%f\t%f\n", 
			ai, (logval)? log(esum) + log((double)Nc) - log((double)h->Nc) : esum * (double)Nc / (double)h->Nc) < 0) 
	      ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      fprintf(pipe, "e\n"); 
    }
  
  return eslOK;
}

static int
cov_plot_lineatexpcov(FILE *pipe, struct data_s *data, double expsurv, int Nc, ESL_HISTOGRAM *h, double ymin, double ymax, char *key, 
		      double offx, double offy, int style)
{
  double cov;
  double posx, posy;
 
  cov = evalue2cov(data, expsurv, Nc, h);
  if (cov <= -eslINFINITY) return eslOK;

  posx = cov + offx;
  posy = ymax - offy;

  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %f,%f '%s' center tc ls %d\n", posx, posy, key, style);
  fprintf(pipe, "plot '-' using 1:2 with lines ls %d \n", style);

  fprintf(pipe, "%f\t%f\n", cov, ymin);
  fprintf(pipe, "%f\t%f\n", cov, ymax);
  fprintf(pipe, "e\n");
   
  return eslOK;
}

int 
cov_histogram_bin2expectsurv(int i, ESL_HISTOGRAM *h, double *ret_expsurv)
{
  double expsurv = 0.;
  int    j;

  if (i > h->imax) 
    expsurv = (double)h->Nc;

  for (j = h->imax; j >= i; j--)
    expsurv += h->obs[j];

  *ret_expsurv = expsurv;

  return eslOK;
}

