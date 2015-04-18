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
static int mutual_naive_psi(int i, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);


int                 
Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, int *ct, FILE *rocfp, int maxFP, int ishuffled,
		 double tol, int verbose, char *errbuf)
{
   int         status;
  
  status = Mutual_Probs(msa, T, ribosum, mi, method, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  status = Mutual_CalculateH  (mi, tol, FALSE, errbuf);
  if (status != eslOK) goto ERROR;

  status = Mutual_CalculateCHI (mi, ct, rocfp, maxFP, ishuffled, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateGTST(mi, ct, rocfp, maxFP, ishuffled, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateMI  (mi, ct, rocfp, maxFP, ishuffled, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateMIp (mi, ct, rocfp, maxFP, ishuffled, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;  
  status = Mutual_CalculateMIa (mi, ct, rocfp, maxFP, ishuffled, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;

#if 0
  status = Mutual_CalculateMIr(mi, ct, rocfp, maxFP, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
#endif

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
  
  /* ps are the marginals */
  for (i = 0; i < mi->alen; i ++) {
    esl_vec_DSet(mi->ps[i], K, 0.0);

    for (j = 0; j < mi->alen; j ++)     
      for (x = 0; x < K; x ++) 
	for (y = 0; y < K; y ++) 
	  mi->ps[i][x] += mi->pp[i][j][IDX(x,y,K)];
    esl_vec_DNorm(mi->ps[i], K);
    status = esl_vec_DValidate(mi->ps[i], K, tol, errbuf);
    if (status != eslOK) {
      printf("ps[%d]\n", i);
      esl_vec_DDump(stdout, mi->ps[i], K, "ACGU");
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
	  printf("ps*ps[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->ps[i][x]*mi->ps[j][y]);
	    }
	  printf("\n");
	  printf("Cp[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->cp[i][j][IDX(x,y,K)]);
	    }
	  printf("\n");
	  printf("\n");
	  printf("Ex[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->nseq*mi->ps[i][x]*mi->ps[j][y]);
	    }
	  printf("\n");
	}
      }
      if (i==5) esl_vec_DDump(stdout, mi->ps[i], K, NULL);
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
  
  /* ps validation */
  for (i = 0; i < mi->alen; i ++) {
    status = esl_vec_DValidate(mi->ps[i], K, tol, errbuf);
    if (status != eslOK) {
      printf("ps[%d]\n", i);
      esl_vec_DDump(stdout, mi->ps[i], K, NULL);
      ESL_XFAIL(eslFAIL, errbuf, "ps validation failed");
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
      H -= log(mi->ps[i][x]) * mi->ps[i][x];
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
Mutual_CalculateCHI(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, double tol, int verbose, char *errbuf)
{
  double chi;
  double chip;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi);
  
  // CHI
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      chi  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  exp = mi->nseq * mi->ps[i][x] * mi->ps[j][y];
	  obs = mi->cp[i][j][IDX(x,y,K)];
	  //printf("%d %d obs  %f exp %f\n", i, j, obs, exp);
	  chi += (obs-exp) * (obs-exp) / exp;
	}	  
      
      /* chi is distributed approximately chi^2. */
      if (chi == 0.) 
	chip = 1.0;
      else if (chi != eslINFINITY) {
	if ((status = esl_stats_ChiSquaredTest(mi->nseq, chi, &chip)) != eslOK) goto ERROR;
      }
      else 
	chip = 0.;

      chip = log(1. - chip);
      chip = chi;
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = chip;
      if (chip < mi->minCOV) mi->minCOV = chip;
      if (chip > mi->maxCOV) mi->maxCOV = chip;
    }
  
  if (verbose) {
    printf("CHI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("CHI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  status = Mutual_SignificantPairs_Ranking(mi, ct, CHI, rocfp, maxFP, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}

int                 
Mutual_CalculateGTST(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, double tol, int verbose, char *errbuf)
{
  double G;
  double Gp;
  double obs;
  double exp;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi);
  
  // G
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      G  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  exp = mi->nseq * mi->ps[i][x] * mi->ps[j][y];
	  obs = mi->cp[i][j][IDX(x,y,K)];
	  G += obs * log (obs / exp);
	}	  
      G *= 2.0;

      /* G is distributed approximately chi^2. */
      if (G == 0.) 
	Gp = 1.0;
      else if (G != eslINFINITY) {
	if ((status = esl_stats_ChiSquaredTest(mi->nseq, G, &Gp)) != eslOK) goto ERROR;
      }
      else 
	Gp = 0.;

      Gp = log(1. - Gp);
      Gp = G;
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = Gp;
      if (Gp < mi->minCOV) mi->minCOV = Gp;
      if (Gp > mi->maxCOV) mi->maxCOV = Gp;
    }
  
  if (verbose) {
    printf("GTST[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("GTST[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  status = Mutual_SignificantPairs_Ranking(mi, ct, GTST, rocfp, maxFP, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}

int                 
Mutual_CalculateMI(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, double tol, int verbose, char *errbuf)
{
  double mutinf;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi);
  
  // MI
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      mutinf  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  mutinf += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i]  = mutinf;
      if (mutinf < mi->minCOV) mi->minCOV = mutinf;
      if (mutinf > mi->maxCOV) mi->maxCOV = mutinf;
    }
  
  if (verbose) {
    printf("MI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==188) printf("MI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  status = Mutual_SignificantPairs_Ranking(mi, ct, MI, rocfp, maxFP, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}


int                 
Mutual_CalculateMIa(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, double tol, int verbose, char *errbuf)
{
  ESL_DMATRIX *mutinf = NULL;
  double      *MIx = NULL;
  double       MIval;
  double       MIavg = 0.0;
  int          i, j;
  int          x, y;
  int          K = mi->abc->K;
  int          status = eslOK;

  Mutual_ReuseCOV(mi);
  
  // MI
  mutinf = esl_dmatrix_Create(mi->alen, mi->alen);
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
       MIval  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  MIval += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
 
      mutinf->mx[i][j] = mutinf->mx[j][i] = MIval;
     }

  // MIavg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) 
      MIavg += mutinf->mx[i][j];
  if (mi->alen > 1) MIavg /= (double)mi->alen * ((double)mi->alen-1.);

  //MIx
  ESL_ALLOC(MIx, sizeof(double) * mi->alen);
  for (i = 0; i < mi->alen; i++) {
    MIx[i] = 0.0;
    for (j = 0; j < mi->alen; j++) {
      if (j != i) MIx[i] += mutinf->mx[i][j];
    }
    if (mi->alen > 1) MIx[i] /= (double)mi->alen - 1.;
  }

  //MIa
  for (i = 0; i < mi->alen; i++) 
    for (j = 0; j < mi->alen; j++) {
      mi->COV->mx[i][j] = mutinf->mx[i][j] - (MIx[i] + MIx[j] - MIavg);
      
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
   }
  
  if (verbose) {
    printf("MIa[%f,%f] \n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MIa[%d][%d] = %f | MI %f | MIx %f MIy %f | MIavg %f\n", 
				i, j, mi->COV->mx[i][j], mutinf->mx[i][j], MIx[i], MIx[j], MIavg);
      } 
  }
  
  status = Mutual_SignificantPairs_Ranking(mi, ct, MIa, rocfp, maxFP, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  free(MIx);
  esl_dmatrix_Destroy(mutinf);
  return status;

 ERROR:
  if (MIx) free(MIx);
  if (mutinf)  esl_dmatrix_Destroy(mutinf);
  return status;
}

int                 
Mutual_CalculateMIp(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, double tol, int verbose, char *errbuf)
{
  ESL_DMATRIX *mutinf = NULL;
  double      *MIx = NULL;
  double       MIval;
  double       MIavg = 0.0;
  int          i, j;
  int          x, y;
  int          K = mi->abc->K;
  int          status = eslOK;
  
  Mutual_ReuseCOV(mi);
  
  // MI
  mutinf = esl_dmatrix_Create(mi->alen, mi->alen);
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
       MIval  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  MIval += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
 
      mutinf->mx[i][j] = mutinf->mx[j][i]  = MIval;
     }

  // MIavg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) 
      MIavg += mutinf->mx[i][j];
  if (mi->alen > 1) MIavg /= (mi->alen * (mi->alen-1));

  //MIx
  ESL_ALLOC(MIx, sizeof(double) * mi->alen);
  for (i = 0; i < mi->alen; i++) {
    MIx[i] = 0.0;
    for (j = 0; j < mi->alen; j++) {
      if (j != i) MIx[i] += mutinf->mx[i][j];
    }
    if (mi->alen > 1) MIx[i] /= mi->alen -1;
  }

  //MIp
  for (i = 0; i < mi->alen; i++) 
    for (j = 0; j < mi->alen; j++) {
      mi->COV->mx[i][j] = (MIavg != 0.0)? mutinf->mx[i][j] - MIx[i] * MIx[j] / MIavg : 0.0;

      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
   }

  if (verbose) {
    printf("MIp[%f,%f] \n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==5&&j==118) printf("MIp[%d][%d] = %f | MI %f | MIx %f MIy %f | MIavg %f\n", 
				i, j, mi->COV->mx[i][j], mutinf->mx[i][j], MIx[i], MIx[j], MIavg);
      } 
  }

  status = Mutual_SignificantPairs_Ranking(mi, ct, MIp, rocfp, maxFP, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  free(MIx);
  esl_dmatrix_Destroy(mutinf);
  return status;

 ERROR:
  if (MIx) free(MIx);
  if (mutinf)  esl_dmatrix_Destroy(mutinf);
  return status;
}

int                 
Mutual_CalculateMIr(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, double tol, int verbose, char *errbuf)
{
  double mutinf, HH;
  int    i, j;
  int    x, y;
  int    K = mi->abc->K;
  int    status = eslOK;
  
  Mutual_ReuseCOV(mi);
  
  //MIr
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      HH  = 0.0;
      mutinf  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  HH -= mi->pp[i][j][IDX(x,y,K)] * log(mi->pp[i][j][IDX(x,y,K)]);
	  mutinf += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
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
  
  status = Mutual_SignificantPairs_Ranking(mi, ct, MIr, rocfp, maxFP, ishuffled, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  return status;

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
  
  ESL_ALLOC(mi->cp,           sizeof(double **) * alen);
  ESL_ALLOC(mi->pp,           sizeof(double **) * alen);
  ESL_ALLOC(mi->ps,           sizeof(double  *) * alen);
  for (i = 0; i < alen; i++) {
    ESL_ALLOC(mi->cp[i],      sizeof(double  *) * alen);
    ESL_ALLOC(mi->pp[i],      sizeof(double  *) * alen);
    ESL_ALLOC(mi->ps[i],      sizeof(double   ) * K);
    for (j = 0; j < alen; j++) {
      ESL_ALLOC(mi->cp[i][j], sizeof(double   ) * K2);
      ESL_ALLOC(mi->pp[i][j], sizeof(double   ) * K2);
    }
  }
   
  mi->COV  = esl_dmatrix_Create(alen, alen);
  ESL_ALLOC(mi->H, sizeof(double) * alen);
 
  /* initialize for adding counts */
  for (i = 0; i < alen; i++) {
    esl_vec_DSet(mi->ps[i], K, 0.0); 
    mi->H[i]  = 0.0;

    for (j = 0; j < alen; j++) {
      esl_vec_DSet(mi->cp[i][j], K2, 0.0); 
      esl_vec_DSet(mi->pp[i][j], K2, 0.0); 
    }
  }

  /* inititalize to zero the COV matrix */
  Mutual_ReuseCOV(mi);
  
  return mi;
  
 ERROR:
  return NULL;
}

int
Mutual_ReuseCOV(struct mutual_s *mi)
{
  int  i, j;
  
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
	free(mi->cp[i][j]);
	free(mi->pp[i][j]);
      }
      free(mi->cp[i]);
      free(mi->pp[i]);
      free(mi->ps[i]);
    }
    esl_dmatrix_Destroy(mi->COV);
    free(mi->cp);
    free(mi->pp);
    free(mi->ps);
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
Mutual_NaivePS(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i;
  int     status;

  for (i = 0; i < alen; i ++) {
    status = mutual_naive_psi(i, msa, mi, tol, verbose, errbuf);
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
Mutual_PostOrderPS(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  ESL_DMATRIX **CL = NULL;
  ESL_DMATRIX **CR = NULL;
  int64_t       alen = msa->alen;
  int           nnodes;
  int           v;
  int           i;
  int           status;

  nnodes = (T->N > 1)? T->N - 1 : 1;
  ESL_ALLOC(CL, sizeof(ESL_DMATRIX *) * nnodes);
  ESL_ALLOC(CR, sizeof(ESL_DMATRIX *) * nnodes);
  for (v = 0; v < nnodes; v++) {
    CL[v] = ratematrix_ConditionalsFromRate(T->ld[v], ribosum->prnaQ, tol, errbuf, verbose);
    if (CL[v] == NULL) goto ERROR;
    CR[v] = ratematrix_ConditionalsFromRate(T->rd[v], ribosum->prnaQ, tol, errbuf, verbose);
    if (CR[v] == NULL) goto ERROR;
  }
 
  for (i = 0; i < alen; i ++) {
    status = mutual_postorder_psi(i, msa, T, ribosum, mi, CL, CR, tol, verbose, errbuf);
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

int
Mutual_SignificantPairs_Ranking(struct mutual_s *mi, int *ct, MITYPE whichmi, FILE *rocfp, int maxFP, int ishuffled, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx = mi->COV;
  char        *mitype;
  double       delta = 200.;
  double       inc;
  double       min = mi->minCOV;
  double       max = mi->maxCOV;
  double       sen;
  double       ppv;
  double       F;
  double       bestF = 0.0;
  double       bestsen;
  double       bestppv;
  double       thresh;
  double       besthresh = mi->besthreshCOV;
  double       maxFPF;
  double       maxFPsen;
  double       maxFPppv;
  double       maxFPthresh;
  double       oneFPF;
  double       oneFPsen;
  double       oneFPppv;
  double       oneFPthresh;
  int          fp, tf, t, f;
  int          oneFP_tf, oneFP_t, oneFP_f, oneFP_fp;
  int          best_tf, best_t, best_f, best_fp;
  int          maxFP_tf, maxFP_t, maxFP_f, maxFP_fp;
  int          nt = 0;
  int          nf = 0;
  int          i, j;
  int          status;

  switch(whichmi) {
  case CHI:  fprintf(rocfp, "\n# CHI "); esl_sprintf(&mitype, "CHI");  break;
  case GTST: fprintf(rocfp, "\n# GTST"); esl_sprintf(&mitype, "GTST"); break;
  case MI:   fprintf(rocfp, "\n# MI  "); esl_sprintf(&mitype, "MI");   break;
  case MIa:  fprintf(rocfp, "\n# MIa "); esl_sprintf(&mitype, "MIa");  break; 
  case MIp:  fprintf(rocfp, "\n# MIp "); esl_sprintf(&mitype, "MIp");  break; 
  case MIr:  fprintf(rocfp, "\n# MIr "); esl_sprintf(&mitype, "MIr");  break; 
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MITYPE");
  }
  if (ishuffled) fprintf(rocfp, " shuffled ");

  fprintf(rocfp, "thresh fp tp true found sen ppv F\n"); 
  
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
    
    fp  = f - tf;
    sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;
    
    fprintf(rocfp, "%.5f %d %d %d %d %.2f %.2f %.2f\n", thresh, fp, tf, t, f, sen, ppv, F);
    
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
 
  if (best_fp < maxFP_fp) {
    if (best_fp == 0 && oneFP_tf > best_tf) {
      printf("%s before 1FP %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", mitype, oneFPthresh, min, max,
	     oneFP_fp, oneFP_tf, oneFP_t, oneFP_f, oneFPsen, oneFPppv, oneFPF);
      for (i = 0; i < mi->alen-1; i++) 
	for (j = i+1; j < mi->alen; j++) {
	  if (mtx->mx[i][j] > oneFPthresh) {
	    if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, mitype, i, j, mtx->mx[i][j]); }
	    else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, mitype, i, j, mtx->mx[i][j]); }
	  }
	}
    }
    else {
      printf("%s optimalF %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", mitype, besthresh, min, max,
	     best_fp, best_tf, best_t, best_f, bestsen, bestppv, bestF);
      for (i = 0; i < mi->alen-1; i++) 
	for (j = i+1; j < mi->alen; j++) {
	  if (mtx->mx[i][j] > besthresh) {
	    if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, mitype, i, j, mtx->mx[i][j]); }
	    else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, mitype, i, j, mtx->mx[i][j]); }
	  }
	}
    }
  }
  else {
    printf("%s maxFP=%d %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", mitype, maxFP, maxFPthresh, min, max,
	   maxFP_fp, maxFP_tf, maxFP_t, maxFP_f, maxFPsen, maxFPppv, maxFPF);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (mtx->mx[i][j] > maxFPthresh) {
	  if (ct[i+1] == j+1) { nt ++; printf("*[%d] %s[%d][%d] = %f\n", nt, mitype, i, j, mtx->mx[i][j]); }
	  else                { nf ++; printf("[%d]  %s[%d][%d] = %f\n", nf, mitype, i, j, mtx->mx[i][j]); } 
	}
      }
  }
  
  return eslOK;

 ERROR:
  return status;
}


/*---------------- internal functions --------------------- */


static int    
mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  double *cp = mi->cp[i][j];
  double *pp = mi->pp[i][j];
  int          K = mi->abc->K;
  int          s;
  int          resi, resj;
  int          x, y;

  esl_vec_DSet(cp, K*K, 1.0); // laplace prior

  for (s = 0; s < msa->nseq; s ++) {
    resi = msa->ax[s][i+1];
    resj = msa->ax[s][j+1];
    if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) cp[IDX(resi,resj,K)] += 1.0;
    else if (esl_abc_XIsCanonical(msa->abc, resi)) for (y = 0; y < K; y ++)           cp[IDX(resi,y,   K)] += 1./(double)K;
    else if (esl_abc_XIsCanonical(msa->abc, resj)) for (x = 0; x < K; x ++)           cp[IDX(x,   resj,K)] += 1./(double)K;
    else {
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) 
	  cp[IDX(x,y,K)] += 1./(double)(K*K);
    }
  }
  esl_vec_DCopy(cp, K*K, mi->cp[j][i]);  // symetrize counts

  /* the probabilities */
  esl_vec_DCopy(cp, K*K, pp);            // get the counts
  esl_vec_DNorm(pp, K*K);                // normalize
  esl_vec_DCopy(pp, K*K, mi->pp[j][i]);  // symetrize

  return eslOK;
}



static int    
mutual_naive_psi(int i, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  double *ps = mi->ps[i];
  int     K = mi->abc->K;
  int     s;
  int     resi;
  int     x;

  esl_vec_DSet(ps, K, 1.0); // laplace prior

  for (s = 0; s < msa->nseq; s ++) {
    resi = msa->ax[s][i+1];
    if (esl_abc_XIsCanonical(msa->abc, resi)) ps[resi] += 1.0;
    else {             
      for (x = 0; x < K; x ++) ps[x] += 1./(double)K;
    }
  }

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
mutual_postorder_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, ESL_DMATRIX **CL, ESL_DMATRIX **CR, 
		     double tol, int verbose, char *errbuf)
{
  ESL_STACK     *vs = NULL;   /* node index stack */
  double       **lk = NULL;
  double        *lkl, *lkr;
  ESL_DMATRIX   *cl, *cr;
  double         sc;
  int            dim;
  int            K = mi->abc->K;
  int            nnodes;
  int            v;
  int            idx;
  int            which;
  int            x;
  int            yl, yr;
  int            resi;
  int            status;
  
  p7_FLogsumInit();    

  /* allocate the single and pair probs for theinternal nodes */
  nnodes = (T->N > 1)? T->N-1 : T->N;
  dim    = nnodes + T->N;
  ESL_ALLOC(lk, sizeof(double *) * dim);
  for (v = 0; v < dim; v ++) lk[v] = NULL;
 
  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      which = (T->left[v] <= 0)? -T->left[v] : T->left[v];
      idx   = (T->left[v] <= 0)? nnodes + which : which;
      if (T->left[v] <= 0) {
	resi = msa->ax[which][i+1];
	ESL_ALLOC(lk[idx], sizeof(double)*K);
	esl_vec_DSet(lk[idx], K, -eslINFINITY);
	if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  lk[idx][resi] = 0.0;
	}
	else {
	  esl_vec_DSet(lk[idx], K, -log((double)K));
	}

      }
      lkl = lk[idx];
      
      which = (T->right[v] <= 0)? -T->right[v] : T->right[v];
      idx   = (T->right[v] <= 0)? nnodes + which : which;
      if (T->right[v] <= 0) {
	resi = msa->ax[which][i+1];
	ESL_ALLOC(lk[idx], sizeof(double)*K);
	esl_vec_DSet(lk[idx], K, -eslINFINITY);
	if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  lk[idx][resi] = 0.0;
	}
	else {
	  esl_vec_DSet(lk[idx], K, -log((double)K));
	}
      }
      lkr = lk[idx];
      
      if (lkl != NULL && lkr != NULL) { /* ready to go: calculate lk at the parent node */
	cl = CL[v];
	cr = CR[v];

	ESL_ALLOC(lk[v], sizeof(double)*K);
	for (x = 0; x < K; x ++) {
	  sc = -eslINFINITY;
	  for (yl = 0; yl < K; yl ++)
	    for (yr = 0; yr < K; yr ++)
	      sc = p7_FLogsum(sc, lkl[yl] + log(cl->mx[x][yl]) + lkr[yr] + log(cr->mx[x][yr]));
	}
	
	/* push node into stack unless already at the root */
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
  for (x = 0; x < K; x ++) mi->ps[i][x] = exp(lk[v][x]) * ribosum->xrnaM[x];
  esl_vec_DNorm(mi->ps[i], K);

  for (v = 0; v < dim; v ++) free(lk[v]);
  free(lk);
  esl_stack_Destroy(vs);
  return eslOK;

 ERROR:
  if (vs) esl_stack_Destroy(vs);
  for (v = 0; v < dim; v ++) if (lk[v]) free(lk[v]);
  if (lk) free(lk);
  return status;
}


