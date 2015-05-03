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

#include "mutualinfo.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"

static int mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				 ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);

int                 
Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, int *ct, FILE *rocfp, FILE *sumfp, 
		 int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs, int ishuffled, double tol, int verbose, char *errbuf)
{
   int         status;
  
  status = Mutual_Probs(msa, T, ribosum, mi, method, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  status = Mutual_CalculateH(mi, tol, FALSE, errbuf);
  if (status != eslOK) goto ERROR;
    
  status = Mutual_CalculateCHI         (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, APC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 
  status = Mutual_CalculateCHI         (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      FALSE, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, ASC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 

  status = Mutual_CalculateOMES        (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, APC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 
  status = Mutual_CalculateOMES        (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      FALSE, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, ASC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 
  
  status = Mutual_CalculateGT          (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, APC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 
  status = Mutual_CalculateGT          (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      FALSE, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, ASC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 

  status = Mutual_CalculateMI          (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, APC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 
  status = Mutual_CalculateMI          (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      FALSE, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, ASC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  status = Mutual_CalculateMIr         (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, APC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR; 
  status = Mutual_CalculateMIr         (mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled,      FALSE, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateCOVCorrected(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, ASC, TRUE,  tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  fprintf(sumfp, "\n");
  return eslOK;
  
 ERROR:
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
Mutual_CalculateH(struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  double H;
  int    i;
  int    x;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  // H
  for (i = 0; i < mi->alen; i++) {
    H = 0.0;
    for (x = 0; x < K; x ++)
      H -= log(mi->pm[i][x]) * mi->pm[i][x];
    mi->H[i] = H;
  }
  
  if (verbose) {
    printf("H\n");
    for (i = 0; i < mi->alen; i++) 
      printf("H[%d] = %f \n", i, mi->H[i]);
  } 
  
  return status;
}


int                 
Mutual_CalculateCHI(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs, 
		    int ishuffled, int analyze, double tol, int verbose, char *errbuf)
{
  double chi;
  double chip;
  double val;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, CHI);
  
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
      
      /* chi is distributed approximately chi^2. */
      if (chi == 0.) 
	chip = 1.0;
      else if (chi != eslINFINITY) {
	if ((status = esl_stats_ChiSquaredTest(mi->nseff[i][j], chi, &chip)) != eslOK) goto ERROR;
      }
      else 
	chip = 0.;

      val = chi;
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = val;
      if (val < mi->minCOV) mi->minCOV = val;
      if (val > mi->maxCOV) mi->maxCOV = val;
    }
  
  if (verbose) {
    printf("CHI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("CHI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }

  if (analyze) status = Mutual_SignificantPairs_Ranking(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}

int                 
Mutual_CalculateOMES(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs, 
		     int ishuffled, int analyze, double tol, int verbose, char *errbuf)
{
  double omes;
  double omesp;
  double val;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, OMES);
  
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
      
      /* omes is distributed approximately omes^2. */
      if (omes == 0.) 
	omesp = 1.0;
      else if (omes != eslINFINITY) {
	if ((status = esl_stats_ChiSquaredTest(mi->nseff[i][j], omes, &omesp)) != eslOK) goto ERROR;
      }
      else 
	omesp = 0.;

      val = -omes;
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = val;
      if (val < mi->minCOV) mi->minCOV = val;
      if (val > mi->maxCOV) mi->maxCOV = val;
    }
  
  if (verbose) {
    printf("OMES[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("OMES[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}


int                 
Mutual_CalculateGT(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs, 
		   int ishuffled, int analyze, double tol, int verbose, char *errbuf)
{
  double gt;
  double gtp;
  double val;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, GT);
  
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

      /* GT is distributed approximately chi^2. */
      if (gt == 0.) 
	gtp = 1.0;
      else if (gt != eslINFINITY) {
	if ((status = esl_stats_ChiSquaredTest(mi->nseff[i][j], gt, &gtp)) != eslOK) goto ERROR;
      }
      else 
	gtp = 0.;

      val = gt;
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = val;
      if (val < mi->minCOV) mi->minCOV = val;
      if (val > mi->maxCOV) mi->maxCOV = val;
    }
  
  if (verbose) {
    printf("GT[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("GT[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}

int                 
Mutual_CalculateMI(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs, 
		   int ishuffled, int analyze, double tol, int verbose, char *errbuf)
{
  double mutinf;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MI);
  
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
  
  if (verbose) {
    printf("MI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==188) printf("MI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}



int                 
Mutual_CalculateMIr(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs,
		    int ishuffled, int analyze, double tol, int verbose, char *errbuf)
{
  double mutinf, HH;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi, MIr);
  
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
  
  if (verbose) {
    printf("MIr[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MIr[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  if (analyze) status = Mutual_SignificantPairs_Ranking(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}



int                 
Mutual_CalculateCOVCorrected(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs, 
			     int ishuffled, CORRTYPE corrtype, int analyze, double tol, int verbose, char *errbuf)
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
  Mutual_ReuseCOV(mi, mi->type);  
 
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

  if (analyze) status = Mutual_SignificantPairs_Ranking(mi, ct, rocfp, sumfp, maxFP, maxDecoy, expectFP, ratioFP, nbpairs, ishuffled, verbose, errbuf);
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

  case CHIp:  esl_sprintf(ret_covtype, "CHIp");  break;
  case GTp:   esl_sprintf(ret_covtype, "GTp");   break;
  case OMESp: esl_sprintf(ret_covtype, "OMESp"); break;
  case MIp:   esl_sprintf(ret_covtype, "MIp");   break;
  case MIrp:  esl_sprintf(ret_covtype, "MIrp");  break; 

  case CHIa:  esl_sprintf(ret_covtype, "CHIa");  break;
  case GTa:   esl_sprintf(ret_covtype, "GTa");   break;
  case OMESa: esl_sprintf(ret_covtype, "OMESa"); break;
  case MIa:   esl_sprintf(ret_covtype, "MIa");   break;
  case MIra:  esl_sprintf(ret_covtype, "MIra");  break; 

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
  else if (!esl_strcmp(covtype, "CHIp"))   type = CHIp;
  else if (!esl_strcmp(covtype, "GTp"))    type = GTp;
  else if (!esl_strcmp(covtype, "OMESp"))  type = OMESp;
  else if (!esl_strcmp(covtype, "MIp"))    type = MIp;
  else if (!esl_strcmp(covtype, "MIrp"))   type = MIrp;
  else if (!esl_strcmp(covtype, "CHIa"))   type = CHIa;
  else if (!esl_strcmp(covtype, "GTa"))    type = GTa;
  else if (!esl_strcmp(covtype, "OMESa"))  type = OMESa;
  else if (!esl_strcmp(covtype, "MIa"))    type = MIa;
  else if (!esl_strcmp(covtype, "MIra"))   type = MIra;
  else
    ESL_XFAIL(eslFAIL, errbuf, "wrong COVTYPE %s", covtype);

  *ret_type = type;
  return eslOK;
  
 ERROR:
  return status;
}

struct mutual_s *
Mutual_Create(int64_t alen, int64_t nseq, ESL_ALPHABET *abc)
{
  struct mutual_s *mi = NULL;
  int              K  = abc->K;
  int              K2 = K * K;
  int              i, j;
  int              status;
  
  ESL_ALLOC(mi, sizeof(struct mutual_s));
  mi->alen = alen;
  mi->nseq = nseq;
  mi->abc  = abc;

  ESL_ALLOC(mi->pp,           sizeof(double **) * alen);
  ESL_ALLOC(mi->nseff,        sizeof(int     *) * alen);
  ESL_ALLOC(mi->pm,           sizeof(double  *) * alen);
  for (i = 0; i < alen; i++) {
    ESL_ALLOC(mi->pp[i],      sizeof(double  *) * alen);
    ESL_ALLOC(mi->nseff[i],   sizeof(int      ) * alen);
    ESL_ALLOC(mi->pm[i],      sizeof(double   ) * K);
    for (j = 0; j < alen; j++) {
       ESL_ALLOC(mi->pp[i][j], sizeof(double  ) * K2);
    }
  }
   
  mi->COV  = esl_dmatrix_Create(alen, alen);
  ESL_ALLOC(mi->H, sizeof(double) * alen);
 
  /* initialize for adding counts */
  for (i = 0; i < alen; i++) {
    esl_vec_DSet(mi->pm[i], K, 0.0); 
    mi->H[i]  = 0.0;

    for (j = 0; j < alen; j++) {
      mi->nseff[i][j] = 0;
      esl_vec_DSet(mi->pp[i][j], K2, 0.0); 
    }
  }

  /* inititalize to zero the COV matrix */
  Mutual_ReuseCOV(mi, COVNONE);
  
  return mi;
  
 ERROR:
  return NULL;
}

int
Mutual_ReuseCOV(struct mutual_s *mi, COVTYPE mitype)
{
  int  i, j;
  
  mi->type = mitype;

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
      free(mi->pp[i]);
      free(mi->pm[i]);
    }
    esl_dmatrix_Destroy(mi->COV);
    free(mi->nseff);
    free(mi->pp);
    free(mi->pm);
    free(mi->H);
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
Mutual_SignificantPairs_Ranking(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, int maxDecoy, double expectFP, double ratioFP, int nbpairs, 
				int ishuffled, int verbose, char *errbuf)
{
  int status;

  if (!ishuffled) {
    status = Mutual_SignificantPairs(mi, ct, rocfp, sumfp, maxFP, expectFP, ratioFP, nbpairs, verbose, errbuf);
  }
  else {
    status = Mutual_SignificantPairs_Shuffled(mi, ct, rocfp, sumfp, maxDecoy, expectFP, ratioFP, nbpairs, verbose, errbuf);
  }
  
  return status;
}

int
Mutual_SignificantPairs(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, double ratioFP,
			int nbpairs, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx = mi->COV;
  char        *covtype = NULL;
  double       delta = 200.;
  double       inc;
  double       min = mi->minCOV;
  double       max = mi->maxCOV;
  double       ratio;
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
  double       ratioFPF;
  double       ratioFPsen;
  double       ratioFPppv;
  double       ratioFPthresh;
  double       expectFPF;
  double       expectFPsen;
  double       expectFPppv;
  double       expectFPthresh;
  double       oneFPF;
  double       oneFPsen;
  double       oneFPppv;
  double       oneFPthresh;
  double       ratioTF_frac_total;    // fraction of covarying basepairs relative to the total number of basepairs
  double       ratioTF_frac_surv;     // fraction of covarying basepairs relative to the basepairs that survive the gapthresh
  double       expectTF_frac_total;   // fraction of covarying basepairs relative to the total number of basepairs
  double       expectTF_frac_surv;    // fraction of covarying basepairs relative to the basepairs that survive the gapthresh
  int          fp, tf, t, f, neg;
  int          oneFP_tf, oneFP_t, oneFP_f, oneFP_fp;
  int          best_tf, best_t, best_f, best_fp;
  int          maxFP_tf, maxFP_t, maxFP_f, maxFP_fp;
  int          ratioFP_tf, ratioFP_t, ratioFP_f, ratioFP_fp;
  int          expectFP_tf, expectFP_t, expectFP_f, expectFP_fp;
  int          nt = 0;
  int          nf = 0;
  int          i, j;

  Mutual_COVTYPEString(&covtype, mi->type, errbuf);

  fprintf(rocfp, "\n# %s ", covtype);  
  fprintf(rocfp, "thresh fp tf found true negatives sen ppv F\n"); 
  
  inc = (max - min) / delta;
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
    ratio  = (t > 0)? (double)fp/(double)t : 0.0;
    sen    = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv    = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F      = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;
    expect = (mi->alen > 0)? (double)fp/(double)mi->alen : 0.0;
    
    neg = mi->alen * (mi->alen-1) / 2 - t;
    fprintf(rocfp, "%.5f %d %d %d %d %d %.2f %.2f %.2f\n", thresh, fp, tf, f, t, neg, sen, ppv, F);
    
    if (ratioFP >= 0.) {
      if (ratio <= ratioFP) {
	ratioFPF      = F;
	ratioFPsen    = sen;
	ratioFPppv    = ppv;
	ratioFP_tf    = tf;
	ratioFP_f     = f;
	ratioFP_t     = t;
	ratioFP_fp    = fp;
	ratioFPthresh = thresh;
      }
    }
    if (expectFP >= 0.) {
      if (ratio <= expectFP) {
	expectFPF      = F;
	expectFPsen    = sen;
	expectFPppv    = ppv;
	expectFP_tf    = tf;
	expectFP_f     = f;
	expectFP_t     = t;
	expectFP_fp    = fp;
	expectFPthresh = thresh;
      }
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
      
  if (expectFP >= 0.0) {
    expectTF_frac_total = (nbpairs    > 0)? 100.*(double)expectFP_tf/(double)nbpairs    : 0.0;
    expectTF_frac_surv  = (expectFP_t > 0)? 100.*(double)expectFP_tf/(double)expectFP_t : 0.0;
    fprintf(sumfp, "%s\t%.2f\t%.2f\t", covtype, expectTF_frac_surv, expectTF_frac_total);

    printf("%s expectFP=%f %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, expectFP, expectFPthresh, min, max,
	   expectFP_fp, expectFP_tf, expectFP_t, expectFP_f, expectFPsen, expectFPppv, expectFPF);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (mtx->mx[i][j] > expectFPthresh) {
	  if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	  else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); } 
	}
      }
  }

  if (ratioFP >= 0.0) {
    ratioTF_frac_total = (nbpairs   > 0)? 100.*(double)ratioFP_tf/(double)nbpairs   : 0.0;
    ratioTF_frac_surv  = (ratioFP_t > 0)? 100.*(double)ratioFP_tf/(double)ratioFP_t : 0.0;
    fprintf(sumfp, "%s\t%.2f\t%.2f\t", covtype, ratioTF_frac_surv, ratioTF_frac_total);

    printf("%s ratioFP=%f %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, ratioFP, ratioFPthresh, min, max,
	   ratioFP_fp, ratioFP_tf, ratioFP_t, ratioFP_f, ratioFPsen, ratioFPppv, ratioFPF);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (mtx->mx[i][j] > ratioFPthresh) {
	  if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	  else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); } 
	}
      }

    ratioTF_frac_total = (nbpairs   > 0)? 100.*(double)ratioFP_tf/(double)nbpairs   : 0.0;
    ratioTF_frac_surv  = (ratioFP_t > 0)? 100.*(double)ratioFP_tf/(double)ratioFP_t : 0.0;
    fprintf(sumfp, "%s\t%.2f\t%.2f\t", covtype, ratioTF_frac_surv, ratioTF_frac_total);
  }    
  
  if (maxFP >= 0) {
    if (best_fp < maxFP_fp) {
      if (best_fp == 0 && oneFP_tf > best_tf) {
	printf("%s before 1FP %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, oneFPthresh, min, max,
	       oneFP_fp, oneFP_tf, oneFP_t, oneFP_f, oneFPsen, oneFPppv, oneFPF);
	for (i = 0; i < mi->alen-1; i++) 
	  for (j = i+1; j < mi->alen; j++) {
	    if (mtx->mx[i][j] > oneFPthresh) {
	      if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	      else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); }
	    }
	  }
      }
      else {
	printf("%s optimalF %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, besthresh, min, max,
	       best_fp, best_tf, best_t, best_f, bestsen, bestppv, bestF);
	for (i = 0; i < mi->alen-1; i++) 
	  for (j = i+1; j < mi->alen; j++) {
	    if (mtx->mx[i][j] > besthresh) {
	      if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	      else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); }
	    }
	  }
      }
    }
    else {
      printf("%s maxFP=%d %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, maxFP, maxFPthresh, min, max,
	     maxFP_fp, maxFP_tf, maxFP_t, maxFP_f, maxFPsen, maxFPppv, maxFPF);
      for (i = 0; i < mi->alen-1; i++) 
	for (j = i+1; j < mi->alen; j++) {
	  if (mtx->mx[i][j] > maxFPthresh) {
	    if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	    else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); } 
	  }
	}
    }
  }

  if (covtype) free(covtype); 
  return eslOK;
}

int
Mutual_SignificantPairs_Shuffled(struct mutual_s *mi, int *ct, FILE *rocfp, FILE *sumfp, int maxDecoy, double expectFP, double ratioFP, 
				 int nbpairs, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx = mi->COV;
  char        *covtype = NULL;
  double       delta = 200.;
  double       inc;
  double       min = mi->minCOV;
  double       max = mi->maxCOV;
  double       ratio;
  double       sen;
  double       ppv;
  double       F;
  double       expect;
  double       thresh;
  double       besthresh = mi->besthreshCOV;
  double       ratioFPF;
  double       ratioFPsen;
  double       ratioFPppv;
  double       ratioFPthresh;
  double       expectFPF;
  double       expectFPsen;
  double       expectFPppv;
  double       expectFPthresh;
  double       decoyF;
  double       decoysen;
  double       decoyppv;
  double       decoythresh;
  double       ratioTF_frac_total;    // fraction of covarying basepairs relative to the total number of basepairs
  double       ratioTF_frac_surv;     // fraction of covarying basepairs relative to the basepairs that survive the gapthresh
  double       expectTF_frac_total;   // fraction of covarying basepairs relative to the total number of basepairs
  double       expectTF_frac_surv;    // fraction of covarying basepairs relative to the basepairs that survive the gapthresh
  int          fp, tf, t, f, neg;
  int          ratioFP_tf, ratioFP_t, ratioFP_f, ratioFP_fp;
  int          expectFP_tf, expectFP_t, expectFP_f, expectFP_fp;
  int          decoy_tf, decoy_t, decoy_f, decoy_fp;
  int          nt, nf;
  int          i, j;

  Mutual_COVTYPEString(&covtype, mi->type, errbuf);

  fprintf(rocfp, "\n# %s ", covtype);  
  fprintf(rocfp, " shuffled ");
  fprintf(rocfp, "thresh fp tp true found negatives sen ppv F\n"); 
  
  inc = (max - min) / delta;
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
    ratio  = (t > 0)? (double)fp/(double)t : 0.0;
    sen    = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv    = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F      = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;
    expect = (mi->alen > 0)? (double)fp/(double)mi->alen : 0.0;
    
    neg = mi->alen * (mi->alen-1) / 2 - t;
    fprintf(rocfp, "%.5f %d %d %d %d %d %.2f %.2f %.2f\n", thresh, fp, tf, f, t, neg, sen, ppv, F);
    
    if (ratioFP >= 0.) {
      if (ratio <= ratioFP) {
	ratioFPF      = F;
	ratioFPsen    = sen;
	ratioFPppv    = ppv;
	ratioFP_tf    = tf;
	ratioFP_f     = f;
	ratioFP_t     = t;
	ratioFP_fp    = fp;
	ratioFPthresh = thresh;
      }
    }
    if (expectFP >= 0.) {
      if (ratio <= expectFP) {
	expectFPF      = F;
	expectFPsen    = sen;
	expectFPppv    = ppv;
	expectFP_tf    = tf;
	expectFP_f     = f;
	expectFP_t     = t;
	expectFP_fp    = fp;
	expectFPthresh = thresh;
      }
    }

    if (tf <= maxDecoy) {
	decoyF      = F;
	decoysen    = sen;
	decoyppv    = ppv;
	decoy_tf    = tf;
	decoy_f     = f;
	decoy_t     = t;
	decoy_fp    = fp;
	decoythresh = thresh;
    }
  }
      
  if (expectFP >= 0.0) {
    expectTF_frac_total = (nbpairs    > 0)? 100.*(double)expectFP_tf/(double)nbpairs    : 0.0;
    expectTF_frac_surv  = (expectFP_t > 0)? 100.*(double)expectFP_tf/(double)expectFP_t : 0.0;
    fprintf(sumfp, "%s\t%.2f\t%.2f\t", covtype, expectTF_frac_surv, expectTF_frac_total);

    printf("%s expectFP=%f %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, expectFP, expectFPthresh, min, max,
	   expectFP_fp, expectFP_tf, expectFP_t, expectFP_f, expectFPsen, expectFPppv, expectFPF);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (mtx->mx[i][j] > expectFPthresh) {
	  if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	  else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); } 
	}
      }
  }

  if (ratioFP >= 0.0) {
    ratioTF_frac_total = (nbpairs   > 0)? 100.*(double)ratioFP_tf/(double)nbpairs   : 0.0;
    ratioTF_frac_surv  = (ratioFP_t > 0)? 100.*(double)ratioFP_tf/(double)ratioFP_t : 0.0;
    fprintf(sumfp, "%s\t%.2f\t%.2f\t", covtype, ratioTF_frac_surv, ratioTF_frac_total);

 #if 0
    printf("%s ratioFP=%f %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, ratioFP, ratioFPthresh, min, max,
	   ratioFP_fp, ratioFP_tf, ratioFP_t, ratioFP_f, ratioFPsen, ratioFPppv, ratioFPF);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (mtx->mx[i][j] > ratioFPthresh) {
	  if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	  else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); } 
	}
      }
#endif
  }
    
  if (maxDecoy >= 0) {
    printf("%s maxDecoy=%d %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, maxDecoy, decoythresh, min, max,
	   decoy_fp, decoy_tf, decoy_t, decoy_f, decoysen, decoyppv, decoyF);
#if 0
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (mtx->mx[i][j] > decoythresh) {
	  if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, covtype, i, j, mtx->mx[i][j]); }
	  else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, covtype, i, j, mtx->mx[i][j]); } 
	}
      }
#endif

  }    

  if (covtype) free(covtype); 
  return eslOK;
}


int
Mutual_SignificantPairs_ZScore(struct mutual_s *mi, int *ct, int verbose, char *errbuf)
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
      printf("[%d][%d] %f | %f | %f %f | %f %f\n", i, ipair, zscore, mi->COV->mx[i][ipair], avgi, stdi, avgj, stdj);
    }
  }  
  return eslOK;
}

/*---------------- internal functions --------------------- */


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
    else if (esl_abc_XIsCanonical(msa->abc, resi)) { mi->nseff[i][j] ++; for (y = 0; y < K; y ++) pp[IDX(resi,y,   K)] += 1./(double)K; }
    else if (esl_abc_XIsCanonical(msa->abc, resj)) { mi->nseff[i][j] ++; for (x = 0; x < K; x ++) pp[IDX(x,   resj,K)] += 1./(double)K; }
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

