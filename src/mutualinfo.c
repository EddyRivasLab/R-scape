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
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "mutualinfo.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"

static int mutual_naive_ppij_counts(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				 ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);
static int mutual_naive_psi_counts(int i, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);


int                 
Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf)
{
  double     *MIx = NULL;
  double      H, HH, MI;
  double      MIavg = 0.0;
  int         i, j;
  int         x, y;
  int         K = msa->abc->K;
  int         status;
  
  status = Mutual_Counts(msa, T, ribosum, mi, method, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  //status = Mutual_CalculateCHI(mi, ct, plotroc, maxFP, tol, verbose, errbuf);
  //if (status != eslOK) goto ERROR;
  
  status = Mutual_Counts2Probs(mi, method, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  status = Mutual_CalculateH  (mi, ct, plotroc, maxFP, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateMI (mi, ct, plotroc, maxFP, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateMIa(mi, ct, plotroc, maxFP, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateMIp(mi, ct, plotroc, maxFP, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  status = Mutual_CalculateMIr(mi, ct, plotroc, maxFP, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  return eslOK;
  
 ERROR:
  return status;
}


int                 
Mutual_Counts(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf)
{
  switch(method) {
  case NAIVE:
    status = Mutual_NaivePPCounts(msa, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    status = Mutual_NaivePSCounts(msa, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    break;
  case PHYLO:
    status = Mutual_PostOrderPP(msa, T, ribosum, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    status = Mutual_PostOrderPD(msa, T, ribosum, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    break;
  case DCA:
    break;
  case AKMAEV:
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "bad method option");
  }
  
  if (verbose) {
    for (i = 0; i < mi->alen-1; i ++) {
      for (j = i+1; j < mi->alen; j ++) {
	if (i==3&&j==28) {
	  printf("pp[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->pp[i][j][IDX(x,y,K)]);
	    }
	  printf("\n");
	}
      }
      if (i==3) esl_vec_DDump(stdout, mi->ps[i], K, NULL);
    }
  }
  
  return status;
}

int 
Mutual_Counts2Probs(struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf)
{
  int    i, j;
  int    x, y;
  int    K = msa->abc->K;
  int    status = eslOK;
  
  for (i = 0; i < mi->alen; i ++) 
    esl_vec_DNorm(mi->ps[i], K);
  
  for (i = 0; i < mi->alen; i ++) 
    for (j = 0; j < mi->alen; j ++) 
      esl_vec_DNorm(mi->pp[i][j], K*K);
  
  status = Mutual_ValidateProbs(mi, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  if (verbose) {
    for (i = 0; i < mi->alen-1; i ++) {
      for (j = i+1; j < mi->alen; j ++) {
	if (i==3&&j==28) {
	  printf("pp[%d][%d] = ", i, j);
	  for (x = 0; x < K; x ++) 
	    for (y = 0; y < K; y ++) {
	      printf(" %f ", mi->pp[i][j][IDX(x,y,K)]);
	    }
	  printf("\n");
	}
      }
      if (i==3) esl_vec_DDump(stdout, mi->ps[i], K, NULL);
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
  int    K = msa->abc->K;
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
      printf("H[%d] = %f \n", i, mi->H->mx[i]);
  } 
  
  return status;
}

int                 
Mutual_CalculateMI(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf)
{
  double MI;
  int    i, j;
  int    x, y;
  int    K = msa->abc->K;
  int    status = eslOK;
  
  Mutual_ResetCOV(mi);
  
  // MI
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      MI  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  MI += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
      
      mi->COV->mx[i][j]  = mi->COV->mx[j][i]  = MI;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
    }
  
  if (verbose) {
    printf("MI[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==3&&j==28) printf("MI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  status = Mutual_SignificantPairs_Ranking(ct, mi, MI, plotroc, maxFP, verbose, cerrbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}


int                 
Mutual_CalculateMIa(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf)
{
  ESL_DMATRIX *MI = NULL;
  double      *MIx = NULL;
  double       MIval;
  double       MIavg = 0.0;
  int          i, j;
  int          x, y;
  int          K = msa->abc->K;
  int          status = eslOK;

  Mutual_ResetCOV(mi);
  
  // MI
  MI = esl_dmatrix_Create(mi->alen, mi->alen);
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
       MIval  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  MIval += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
 
      MI->mx[i][j]  = MI->mx[j][i]  = MIval;
     }

  // MIavg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) 
      MIavg += MI->mx[i][j];
  if (mi->alen > 1) MIavg /= (mi->alen * (mi->alen-1));

  //MIx
  ESL_ALLOC(MIx, sizeof(double) * mi->alen);
  for (i = 0; i < mi->alen; i++) {
    MIx[i] = 0.0;
    for (j = 0; j < mi->alen; j++) {
      if (j != i) MIx[i] += MI->mx[i][j];
    }
    if (mi->alen > 1) MIx[i] /= mi->alen -1;
  }

  //MIa
  for (i = 0; i < mi->alen; i++) 
    for (j = 0; j < mi->alen; j++) {
      mi->COV->mx[i][j] = MI->mx[i][j] - (MIx[i] + MIx[j] - MIavg);
      
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
   }
  
  if (verbose) {
    printf("MIp[%f,%f] \n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==3&&j==28) printf("MIa[%d][%d] = %f | MI %f | MIx %f MIy %f | MIavg %f\n", 
				i, j, mi->COV->mx[i][j], MI->mx[i][j], MIx[i], MIx[j], MIavg);
      } 
  }
  
  status = Mutual_SignificantPairs_Ranking(ct, mi, MI, plotroc, maxFP, verbose, cerrbuf);
  if (status != eslOK) goto ERROR;

  free(MIx);
  esl_dmatrix_Destroy(MI);
  return status;

 ERROR:
  if (MIx) free(MIx);
  if (MI)  esl_dmatrix_Destroy(MI);
  return status;
}

int                 
Mutual_CalculateMIp(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf)
{
  ESL_DMATRIX *MI = NULL;
  double      *MIx = NULL;
  double       MIval;
  double       MIavg = 0.0;
  int          i, j;
  int          x, y;
  int          K = msa->abc->K;
  int          status = eslOK;
  
  Mutual_ResetCOV(mi);
  
  // MI
  MI = esl_dmatrix_Create(mi->alen, mi->alen);
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
       MIval  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  MIval += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
 
      MI->mx[i][j]  = MI->mx[j][i]  = MIval;
     }

  // MIavg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) 
      MIavg += MI->mx[i][j];
  if (mi->alen > 1) MIavg /= (mi->alen * (mi->alen-1));

  //MIx
  ESL_ALLOC(MIx, sizeof(double) * mi->alen);
  for (i = 0; i < mi->alen; i++) {
    MIx[i] = 0.0;
    for (j = 0; j < mi->alen; j++) {
      if (j != i) MIx[i] += MI->mx[i][j];
    }
    if (mi->alen > 1) MIx[i] /= mi->alen -1;
  }

  //MIp
  for (i = 0; i < mi->alen; i++) 
    for (j = 0; j < mi->alen; j++) {
      mi->COV->mx[i][j] = (MIavg != 0.0)? MI->mx[i][j] - MIx[i] * MIx[j] / MIavg : 0.0;

      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
   }

  if (verbose) {
    printf("MIp[%f,%f] \n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==3&&j==28) printf("MIp[%d][%d] = %f | MI %f | MIx %f MIy %f | MIavg %f\n", 
				i, j, mi->COV->mx[i][j], MI->mx[i][j], MIx[i], MIx[j], MIavg);
      } 
  }

  status = Mutual_SignificantPairs_Ranking(ct, mi, MI, plotroc, maxFP, verbose, cerrbuf);
  if (status != eslOK) goto ERROR;

  free(MIx);
  esl_dmatrix_Destroy(MI);
  return status;

 ERROR:
  if (MIx) free(MIx);
  if (MI)  esl_dmatrix_Destroy(MI);
  return status;
}

int                 
Mutual_CalculateMIr(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf)
{
  double MI, HH;
  int    i, j;
  int    x, y;
  int    K = msa->abc->K;
  int    status = eslOK;
  
  Mutual_ResetCOV(mi);
  
  //MIr
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      HH  = 0.0;
      MI  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  HH -= mi->pp[i][j][IDX(x,y,K)] * log(mi->pp[i][j][IDX(x,y,K)]);
	  MI += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
      
      mi->COV->mx[i][j] = mi->COV->mx[j][i] = (HH > 0.0)? MI/HH : 0.0;
      if (mi->COV->mx[i][j] < mi->minCOV) mi->minCOV = mi->COV->mx[i][j];
      if (mi->COV->mx[i][j] > mi->maxCOV) mi->maxCOV = mi->COV->mx[i][j];
    }
  
  if (verbose) {
    printf("MIr[%f,%f]\n", mi->minCOV, mi->maxCOV);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==3&&j==28) printf("MI[%d][%d] = %f \n", i, j, mi->COV->mx[i][j]);
      } 
  }
  
  status = Mutual_SignificantPairs_Ranking(ct, mi, MI, plotroc, maxFP, verbose, cerrbuf);
  if (status != eslOK) goto ERROR;

  return status;

 ERROR:
  return status;
}


struct mutual_s *
Mutual_Create(int64_t alen, int64_t nseq, int K)
{
  struct mutual_s *mi = NULL;
  int              K2 = K * K;
  int              i, j;
  int              status;
  
  ESL_ALLOC(mi, sizeof(struct mutual_s));
  mi->alen = alen;
  mi->nseq = nseq;
  
  ESL_ALLOC(mi->pp,           sizeof(double **) * alen);
  ESL_ALLOC(mi->ps,           sizeof(double  *) * alen);
  for (i = 0; i < alen; i++) {
    ESL_ALLOC(mi->pp[i],      sizeof(double  *) * alen);
    ESL_ALLOC(mi->ps[i],      sizeof(double   ) * K);
    for (j = 0; j < alen; j++) 
      ESL_ALLOC(mi->pp[i][j], sizeof(double   ) * K2);
  }
   
  mi->COV  = esl_dmatrix_Create(alen, alen);
  ESL_ALLOC(mi->H, sizeof(double) * alen);
 
  /* initialize for adding counts */
  for (i = 0; i < alen; i++) {
    esl_vec_DSet(mi->ps[i], K, 0.0); 
    mi->H[i]  = 0.0;

    for (j = 0; j < alen; j++) 
      esl_vec_DSet(mi->pp[i][j], K2, 0.0); 
  }

  /* inititalize to zero the COV matrix */
  Mutual_ResetCOV(mi);
  
  return mi;
  
 ERROR:
  return NULL;
}

int
Mutual_ResetCOV(struct mutual_s *mi)
{
  int  i, j;
  
  for (i = 0; i < mi->alen; i++) {
    mi->H[i]  = 0.0;
    for (j = 0; j < mi->alen; j++) 
      mi->COV->mx[i][j]  = 0.0;
  }

  mi->besthresCOV = -eslINFINITY;
  mi->minCOV      =  eslINFINITY;
  mi->maxCOV      = -eslINFINITY;

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
      free(mi->pp[i]);
      free(mi->ps[i]);
    }
    esl_dmatrix_Destroy(mi->COV);
    free(mi->pp);
    free(mi->ps);
    free(mi->H);
    free(mi);
  }
}


int 
Mutual_NaivePPCounts(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i, j;
  int     status;

  for (i = 0; i < alen-1; i ++)
    for (j = i+1; j < alen; j ++) {
      status = mutual_naive_ppij_counts(i, j, msa, mi, tol, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
  
  return eslOK;

 ERROR:
  return status;
}

int 
Mutual_NaivePSCounts(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i;
  int     status;

  for (i = 0; i < alen; i ++) {
    status = mutual_naive_psi_counts(i, msa, mi, tol, verbose, errbuf);
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


static int
Mutual_SignificantPairs_ZScore(int *ct, struct mutual_s *mi, int verbose, char *errbuf)
{
  double       avgi, avgj;
  double       stdi, stdj;
  double       zscorei, zscorej;
  double       zscore;
  int          i, j, k;
  int          ipair;
  int          status;

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
  
 ERROR:
  return status;
}

static int
Mutual_SignificantPairs_Ranking(int *ct, struct mutual_s *mi, MITYPE whichmi, int *ct, int plotroc, int maxFP, int verbose, char *errbuf)
{
  MITYPE       which;
  ESL_DMATRIX *mtx = mi->COV;
  char        *mitype;
  double       delta = 500;
  double       inc;
  double       min = mi->maxCOV;
  double       max = mi->maxCOV;
  double       sen;
  double       ppv;
  double       F;
  double       bestF = 0.0;
  double       bestsen;
  double       bestppv;
  double       thresh;
  double       besthresh = mi->bestreshCOV;
  double       maxFPF;
  double       maxFPsen;
  double       maxFPppv;
  double       maxFPthresh;
  int          tf, t, f;
  int          best_tf, best_t, best_f, best_fp;
  int          maxFP_tf, maxFP_t, maxFP_f, maxFP_fp;
  int          nt = 0;
  int          nf = 0;
  int          i, j;
  int          status;

  switch(whichmi) {
  case CHI: if (plotroc) printf("\n# CHI "); break;
  case MI:  if (plotroc) printf("\n# MI ");  break;
  case MIa: if (plotroc) printf("\n# MIa "); break; 
  case MIp: if (plotroc) printf("\n# MIp "); break; 
  case MIr: if (plotroc) printf("\n# MIr "); break; 
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MITYPE");
  }
  if (plotroc) printf("thresh tp true found sen ppv F\n"); 
  
  inc = (max - min) / delta;
  for (thresh = max; thresh > min-inc; thresh -= inc) {

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
    
    if (plotroc) printf("%.5f %d %d %d %d %.2f %.2f %.2f\n", thresh, fp, tf, t, f, sen, ppv, F);
    
    if (f-tf <= maxFP) {
      maxFPF      = F;
      maxFPsen    = sen;
      maxFPppv    = ppv;
      maxFP_tf    = tf;
      maxFP_f     = f;
      maxFP_t     = t;
      maxFP_fp    = f - tf;
      maxFPthresh = thresh;
    }
    if (F > bestF) { 
      bestF     = F; 
      bestsen   = sen;
      bestppv   = ppv;
      best_tf   = tf;
      best_t    = t;
      best_f    = f;
      best_fp   = f - tf;
      besthresh = thresh; 
    }
  }
 
  if (best_fp < maxFP_fp) {
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



static int    
mutual_naive_ppij_counts(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  double *pp = mi->pp[i][j];
  int          K = msa->abc->K;
  int          s;
  int          resi, resj;
  int          x, y;

  esl_vec_DSet(pp, K*K, 1.0); // laplace prior

  for (s = 0; s < msa->nseq; s ++) {
    resi = msa->ax[s][i+1];
    resj = msa->ax[s][j+1];
    if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) pp[IDX(resi,resj,K)] += 1.0;
    else if (esl_abc_XIsCanonical(msa->abc, resi)) for (y = 0; y < K; y ++)           pp[IDX(resi,y,   K)] += 1.0;
    else if (esl_abc_XIsCanonical(msa->abc, resj)) for (x = 0; x < K; x ++)           pp[IDX(x,   resj,K)] += 1.0;
    else {
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) 
	  pp[IDX(x,y,K)] += 1.0;
    }
  }
  esl_vec_DCopy(pp, K*K, mi->pp[j][i]);
  
  return eslOK;
}


/*---------------- internal functions --------------------- */

static int    
mutual_naive_psi_counts(int i, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  double *ps = mi->ps[i];
  int     K = msa->abc->K;
  int     s;
  int     resi;
  int     x;

  esl_vec_DSet(ps, K, 1.0); // laplace prior

  for (s = 0; s < msa->nseq; s ++) {
    resi = msa->ax[s][i+1];
    if (esl_abc_XIsCanonical(msa->abc, resi)) ps[resi] += 1.0;
    else {             
      for (x = 0; x < K; x ++) ps[x] += 1.0;
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
  int            K = msa->abc->K;
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
	if (i==3&&j==28) {
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
	if (i==3&&j==28) {
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
	if (i==3&&j==28) {
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
  int            K = msa->abc->K;
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


