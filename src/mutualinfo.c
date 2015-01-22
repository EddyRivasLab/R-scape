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

static int mutual_analyze_ranking(int *ct, struct mutual_s *mi, MITYPE whichmi, int plotroc, int maxFP, int verbose, char *errbuf);
static int mutual_analyze_ranking_thresh(int *ct, struct mutual_s *mi, MITYPE whichmi, double thresh, 
					 int *ret_tf, int *ret_t, int *ret_f, double *ret_sen, double *ret_ppv, double *ret_F, 
					 int plotroc, int maxFP, int verbose, char *errbuf);
static int mutual_analyze_significant_pairs(int *ct, struct mutual_s *mi, MITYPE whichmi, int verbose, char *errbuf);
static int mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				 ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);
static int mutual_naive_psi(int i, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				ESL_DMATRIX **CL, ESL_DMATRIX **CR, double tol, int verbose, char *errbuf);

int                 
Mutual_Analyze(int *ct, struct mutual_s *mi, int plotroc, int maxFP, int verbose, char *errbuf)
{
  int status;

  status = Mutual_AnalyzeRanking(ct, mi, plotroc, maxFP, verbose, errbuf);
  // status = Mutual_AnalyzeSignificantPairs(ct, mi, verbose, errbuf);

  return eslOK;
}

int                 
Mutual_AnalyzeSignificantPairs(int *ct, struct mutual_s *mi, int verbose, char *errbuf)
{
  int status;
  
  status = mutual_analyze_significant_pairs(ct, mi, MI, verbose, errbuf);
  status = mutual_analyze_significant_pairs(ct, mi, MIa, verbose, errbuf);
  status = mutual_analyze_significant_pairs(ct, mi, MIp, verbose, errbuf);
  //status = mutual_analyze_significant_pairs(ct, mi, MIr, verbose, errbuf);

  return eslOK;
}

int                 
Mutual_AnalyzeRanking(int *ct, struct mutual_s *mi, int plotroc, int maxFP, int verbose, char *errbuf)
{
  int    status;
  
  status = mutual_analyze_ranking(ct, mi, MI,  plotroc, maxFP, verbose, errbuf);
  status = mutual_analyze_ranking(ct, mi, MIa, plotroc, maxFP, verbose, errbuf);
  status = mutual_analyze_ranking(ct, mi, MIp, plotroc, maxFP, verbose, errbuf);
  //status = mutual_analyze_ranking(ct, mi, MIr, plotroc, maxFP, verbose, errbuf);
 
  return eslOK;
}


int                 
Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf)
{
  double     *MIx = NULL;
  double      H, HH, MI;
  double      MIavg = 0.0;
  int         i, j;
  int         x, y;
  int         K = msa->abc->K;
  int         status;
  
  switch(method) {
  case NAIVE:
    status = Mutual_NaivePP(msa, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    //status = Mutual_NaivePS(msa, mi, tol, verbose, errbuf);
    //if (status != eslOK) goto ERROR;
    break;
  case PHYLO:
    status = Mutual_PostOrderPP(msa, T, ribosum, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;    
    //status = Mutual_PostOrderPS(msa, T, ribosum, mi, tol, verbose, errbuf);
    //if (status != eslOK) goto ERROR;
    break;
  case DCA:
    break;
  case AKMAEV:
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "bad method option");
  }
  
  /* pp validation */
  for (i = 0; i < mi->alen-1; i ++) 
    for (j = i+1; j < mi->alen; j ++) {
      status = esl_vec_DValidate(mi->pp[i][j], K*K, tol, errbuf);
      if (status != eslOK) {
	printf("pp[%d][%d]\n", i, j);
	esl_vec_DDump(stdout, mi->pp[i][j], K*K, NULL);
	goto ERROR;
      }
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

  // H
  for (i = 0; i < mi->alen; i++) {
    H = 0.0;
    for (x = 0; x < K; x ++)
      H -= log(mi->ps[i][x]) * mi->ps[i][x];
    mi->H[i] = H;
  }

  // MI, MIr
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      HH  = 0.0;
      MI  = 0.0;
      for (x = 0; x < K; x ++)
	for (y = 0; y < K; y ++) {
	  HH -= mi->pp[i][j][IDX(x,y,K)] * log(mi->pp[i][j][IDX(x,y,K)]);
	  MI += mi->pp[i][j][IDX(x,y,K)] * ( log(mi->pp[i][j][IDX(x,y,K)]) - log(mi->ps[i][x]) - log(mi->ps[j][y]) );
	}	  
 
      mi->MI->mx[i][j]  = mi->MI->mx[j][i]  = MI;
      mi->MIr->mx[i][j] = mi->MIr->mx[j][i] = (HH > 0.0)? MI/HH : 0.0;
      if (mi->MI->mx[i][j] < mi->minMI) mi->minMI = mi->MI->mx[i][j];
      if (mi->MI->mx[i][j] > mi->maxMI) mi->maxMI = mi->MI->mx[i][j];
    }

  // MIavg
  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) 
      MIavg += mi->MI->mx[i][j];
  if (mi->alen > 1) MIavg /= (mi->alen * (mi->alen-1));

  //MIx
  ESL_ALLOC(MIx, sizeof(double) * mi->alen);
  for (i = 0; i < mi->alen; i++) {
    MIx[i] = 0.0;
    for (j = 0; j < mi->alen; j++) {
      if (j != i) MIx[i] += mi->MI->mx[i][j];
    }
    if (mi->alen > 1) MIx[i] /= mi->alen -1;
  }

  //MIa, MIp
  for (i = 0; i < mi->alen; i++) 
    for (j = 0; j < mi->alen; j++) {
      mi->MIp->mx[i][j] = (MIavg != 0.0)? mi->MI->mx[i][j] - MIx[i] * MIx[j] / MIavg : 0.0;
      mi->MIa->mx[i][j] = mi->MI->mx[i][j] - (MIx[i] + MIx[j] - MIavg);
      if (mi->MIa->mx[i][j] < mi->minMIa) mi->minMIa = mi->MIa->mx[i][j];
      if (mi->MIa->mx[i][j] > mi->maxMIa) mi->maxMIa = mi->MIa->mx[i][j];
      if (mi->MIp->mx[i][j] < mi->minMIp) mi->minMIp = mi->MIp->mx[i][j];
      if (mi->MIp->mx[i][j] > mi->maxMIp) mi->maxMIp = mi->MIp->mx[i][j];
      if (mi->MIr->mx[i][j] < mi->minMIr) mi->minMIr = mi->MIr->mx[i][j];
      if (mi->MIr->mx[i][j] > mi->maxMIr) mi->maxMIr = mi->MIr->mx[i][j];
 
   }

  if (verbose) {
    printf("MI[%f,%f] MIa[%f,%f] MIp[%f,%f] MIr[%f,%f] \n", 
	   mi->minMI, mi->maxMI, mi->minMIa, mi->maxMIa, mi->minMIp, mi->maxMIp, mi->minMIr, mi->maxMIr);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (i==3&&j==28) printf("MI[%d][%d] = %f | %f %f %f | %f %f | %f %f | %f\n", 
					   i, j, mi->MI->mx[i][j], mi->MIa->mx[i][j], mi->MIp->mx[i][j], mi->MIr->mx[i][j],
					   mi->H[i], mi->H[j], MIx[i], MIx[j], MIavg);
      } 
  }

  free(MIx);
  return eslOK;

 ERROR:
  if (MIx) free(MIx);
  return status;
}

struct mutual_s *
Mutual_Create(int64_t alen, int K)
{
  struct mutual_s *mi = NULL;
  int              K2 = K * K;
  int              i, j;
  int              status;

  ESL_ALLOC(mi, sizeof(struct mutual_s));
  mi->alen = alen;
  
  ESL_ALLOC(mi->pp,           sizeof(double **) * alen);
  ESL_ALLOC(mi->ps,           sizeof(double  *) * alen);
  for (i = 0; i < alen; i++) {
    ESL_ALLOC(mi->pp[i],      sizeof(double  *) * alen);
    ESL_ALLOC(mi->ps[i],      sizeof(double   ) * K);
    for (j = 0; j < alen; j++) 
      ESL_ALLOC(mi->pp[i][j], sizeof(double   ) * K2);
  }
   
  mi->MI  = esl_dmatrix_Create(alen, alen);
  mi->MIa = esl_dmatrix_Create(alen, alen);
  mi->MIp = esl_dmatrix_Create(alen, alen);
  mi->MIr = esl_dmatrix_Create(alen, alen);
  ESL_ALLOC(mi->H, sizeof(double) * alen);
 
  for (i = 0; i < alen; i++) {
    esl_vec_DSet(mi->ps[i], K, 0.0); 
    mi->H[i]  = 0.0;
    for (j = 0; j < alen; j++) {
      esl_vec_DSet(mi->pp[i][j], K2, 0.0); 
       
      mi->MI->mx[i][j]  = 0.0;
      mi->MIa->mx[i][j] = 0.0;
      mi->MIp->mx[i][j] = 0.0;
      mi->MIr->mx[i][j] = 0.0;
    }
  }

  mi->besthreshMI  = 0.0;
  mi->besthreshMIa = 0.0;
  mi->besthreshMIp = 0.0;
  mi->besthreshMIr = 0.0;

  mi->minMI  = eslINFINITY;
  mi->minMIa = eslINFINITY;
  mi->minMIp = eslINFINITY;
  mi->minMIr = eslINFINITY;

  mi->maxMI  = -eslINFINITY;
  mi->maxMIa = -eslINFINITY;
  mi->maxMIp = -eslINFINITY;
  mi->maxMIr = -eslINFINITY;

  return mi;

 ERROR:
  return NULL;
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
    esl_dmatrix_Destroy(mi->MI);
    esl_dmatrix_Destroy(mi->MIa);
    esl_dmatrix_Destroy(mi->MIp);
    esl_dmatrix_Destroy(mi->MIr);
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
    //CL[v] = ratematrix_ConditionalsFromRate(T->ld[v], ribosum->bprsQ, tol, errbuf, verbose);
     CL[v] = ratematrix_ConditionalsFromRate(0.1, ribosum->bprsQ, tol, errbuf, verbose);
    if (CL[v] == NULL) goto ERROR;
    //CR[v] = ratematrix_ConditionalsFromRate(T->rd[v], ribosum->bprsQ, tol, errbuf, verbose);
    CR[v] = ratematrix_ConditionalsFromRate(0.1, ribosum->bprsQ, tol, errbuf, verbose);
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

/*---------------- internal functions --------------------- */

static int
mutual_analyze_significant_pairs(int *ct, struct mutual_s *mi, MITYPE whichmi, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx;
  double       avgi, avgj;
  double       stdi, stdj;
  double       zscorei, zscorej;
  double       zscore;
  int          i, j, k;
  int          ipair;
  int          status;

  switch(whichmi) {
  case MI:  mtx = mi->MI;  printf("\n# MI ");  break;
  case MIa: mtx = mi->MIa; printf("\n# MIa "); break; 
  case MIp: mtx = mi->MIp; printf("\n# MIp "); break; 
  case MIr: mtx = mi->MIr; printf("\n# MIr "); break; 
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MITYPE");
  }
  printf("zscore mij avg std\n");

  for (i = 0; i < mi->alen; i ++) {
    if (ct[i+1] > 0 && ct[i+1] > i+1) {
      ipair = ct[i+1]-1;

      avgi = 0.0;
      stdi = 0.0;
      
      for (j = 0; j < mi->alen; j ++) {
	if (j != ipair) {   
	  if (j != i) {
	    avgi += mtx->mx[i][j];
	    stdi += mtx->mx[i][j] * mtx->mx[i][j];
	  }
	}
	else {
	  avgj = 0.0;
	  stdj = 0.0;
	  for (k = 0; k < mi->alen; k ++) {
	    if (k != j && k != i) {       
	      avgj += mtx->mx[j][k];
	      stdj += mtx->mx[j][k] * mtx->mx[j][k];
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

      zscorei = (mtx->mx[i][ipair] - avgi) / stdi;
      zscorej = (mtx->mx[i][ipair] - avgj) / stdj;
      zscore  = ESL_MIN(zscorej, zscorej);
      printf("[%d][%d] %f | %f | %f %f | %f %f\n", i, ipair, zscore, mtx->mx[i][ipair], avgi, stdi, avgj, stdj);
    }
  }  
  return eslOK;
  
 ERROR:
  return status;
}

static int
mutual_analyze_ranking(int *ct, struct mutual_s *mi, MITYPE whichmi, int plotroc, int maxFP, int verbose, char *errbuf)
{
  MITYPE       which;
  ESL_DMATRIX *mtx;
  char        *mitype;
  double       delta = 500;
  double       inc;
  double       min;
  double       max;
  double       sen;
  double       ppv;
  double       F;
  double       bestF = 0.0;
  double       bestsen;
  double       bestppv;
  double       thresh;
  double       besthresh;
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
  case MI:  if (plotroc) printf("\n# MI ");  max = mi->maxMI;  min = mi->minMI;  which = MI;  break;
  case MIa: if (plotroc) printf("\n# MIa "); max = mi->maxMIa; min = mi->minMIa; which = MIa; break; 
  case MIp: if (plotroc) printf("\n# MIp "); max = mi->maxMIp; min = mi->minMIp; which = MIp; break; 
  case MIr: if (plotroc) printf("\n# MIr "); max = mi->maxMIr; min = mi->minMIr; which = MIr; break; 
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MITYPE");
  }
  if (plotroc) printf("thresh tp true found sen ppv F\n"); 
  
  inc = (max - min) / delta;
  for (thresh = max; thresh > min-inc; thresh -= inc) {
    mutual_analyze_ranking_thresh(ct, mi, which, thresh, &tf, &t, &f, &sen, &ppv, &F, plotroc, maxFP, verbose, errbuf);
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

  switch(whichmi) {
  case MI:  mi->besthreshMI  = besthresh;  break;
  case MIa: mi->besthreshMIa = besthresh;  break;
  case MIp: mi->besthreshMIp = besthresh;  break;
  case MIr: mi->besthreshMIr = besthresh;  break;
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MITYPE");
  }
  
  switch(whichmi) {
  case MI:  mtx = mi->MI;  esl_sprintf(&mitype, "MI");  break;
  case MIa: mtx = mi->MIa; esl_sprintf(&mitype, "MIa"); break;
  case MIp: mtx = mi->MIp; esl_sprintf(&mitype, "MIp"); break; 
  case MIr: mtx = mi->MIr; esl_sprintf(&mitype, "MIr"); break;
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MITYPE");
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
mutual_analyze_ranking_thresh(int *ct, struct mutual_s *mi, MITYPE whichmi, double thresh, 
			      int *ret_tf, int *ret_t, int *ret_f, double *ret_sen, double *ret_ppv, double *ret_F, 
			      int plotroc, int maxFP, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx;
  double       sen;
  double       ppv;
  double       F;
  int          tf = 0;
  int          f  = 0;
  int          t  = 0;
  int          fp;
  int          i, j;
  int          status;

  switch(whichmi) {
  case MI:  mtx = mi->MI;  break;
  case MIa: mtx = mi->MIa; break; 
  case MIp: mtx = mi->MIp; break; 
  case MIr: mtx = mi->MIr; break; 
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MITYPE");
  }
  
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

  if (plotroc && fp <= maxFP) printf("%.5f %d %d %d %d %.2f %.2f %.2f\n", thresh, fp, tf, t, f, sen, ppv, F);

  if (ret_tf)  *ret_tf  = tf;
  if (ret_t)   *ret_t   = t;
  if (ret_f)   *ret_f   = f;
  if (ret_sen) *ret_sen = sen;
  if (ret_ppv) *ret_ppv = ppv;
  if (ret_F)   *ret_F   = F;

  return eslOK;

 ERROR:
  return status;
}

static int    
mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
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
  esl_vec_DNorm(pp, K*K);
  esl_vec_DCopy(pp, K*K, mi->pp[j][i]);
  
  return eslOK;
}


static int    
mutual_naive_psi(int i, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
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

  esl_vec_DNorm(ps, K);

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
	  esl_dmatrix_Dump(stdout, cl, NULL, NULL);
	  esl_dmatrix_Dump(stdout, cr, NULL, NULL);
	  esl_dmatrix_Dump(stdout, lkl, "ACGU", "ACGU");
	  esl_dmatrix_Dump(stdout, lkr, "ACGU", "ACGU");
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
      mi->pp[i][j][IDX(x,y,K)] = exp(lk[v]->mx[x][y]) * ribosum->bprsM[IDX(x,y,K)];
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


