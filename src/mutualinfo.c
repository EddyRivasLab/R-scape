/* mutualinfo.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

static int mutual_analyze_ranking(int *ct, struct mutual_s *mi, MItype whichmi, double thresh, int verbose, char *errbuf);
static int mutual_analyze_significant_pairs(int *ct, struct mutual_s *mi, MItype whichmi, int verbose, char *errbuf);
static int mutual_naive_ppij(int i, int j, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
static int mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, 
				  int verbose, char *errbuf);
static int mutual_naive_psi(int i, ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);

static int mutual_postorder_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				 double tol, int verbose, char *errbuf);

int                 
Mutual_Analyze(int *ct, struct mutual_s *mi, int verbose, char *errbuf)
{
  int status;

  status = Mutual_AnalyzeRanking(ct, mi, verbose, errbuf);
  status = Mutual_AnalyzeSignificantPairs(ct, mi, verbose, errbuf);

  return eslOK;
}

int                 
Mutual_AnalyzeSignificantPairs(int *ct, struct mutual_s *mi, int verbose, char *errbuf)
{
  int status;
  
  status = mutual_analyze_significant_pairs(ct, mi, MI, verbose, errbuf);
  status = mutual_analyze_significant_pairs(ct, mi, MIa, verbose, errbuf);
  status = mutual_analyze_significant_pairs(ct, mi, MIp, verbose, errbuf);
  status = mutual_analyze_significant_pairs(ct, mi, MIr, verbose, errbuf);

  return eslOK;
}

int                 
Mutual_AnalyzeRanking(int *ct, struct mutual_s *mi, int verbose, char *errbuf)
{
  double thresh;
  double inc;
  double delta = 50;
  int    status;
  
  inc = (mi->maxMI - mi->minMI)/ delta;
  for (thresh = mi->minMI; thresh < mi->maxMI+inc; thresh += inc) {
    status = mutual_analyze_ranking(ct, mi, MI, thresh, verbose, errbuf);
  }

  inc = (mi->maxMIa - mi->minMIa)/ delta;
  for (thresh = mi->minMIa; thresh < mi->maxMIa+inc; thresh += inc) {
    status = mutual_analyze_ranking(ct, mi, MIa, thresh, verbose, errbuf);
  }

  inc = (mi->maxMIp - mi->minMIp)/ delta;
  for (thresh = mi->minMIp; thresh < mi->maxMIp+inc; thresh += inc) {
    status = mutual_analyze_ranking(ct, mi, MIp, thresh, verbose, errbuf);
  }

  inc = (mi->maxMIr - mi->minMIr)/ delta;
  for (thresh = mi->minMIr; thresh < mi->maxMIr+inc; thresh += inc) {
    status = mutual_analyze_ranking(ct, mi, MIr, thresh, verbose, errbuf);
  }

  return eslOK;
}


int                 
Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int naive, double tol, int verbose, char *errbuf)
{
  double     *MIx = NULL;
  double      H, HH, MI;
  double      MIavg = 0.0;
  int         i, j;
  int         x, y;
  int         K = msa->abc->K;
  int         status;
  
  if (naive) {
    status = Mutual_NaivePP(msa, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;
    
    //status = Mutual_NaivePS(msa, mi, tol, verbose, errbuf);
    //if (status != eslOK) goto ERROR;
  }
  else {
    status = Mutual_PostOrderPP(msa, T, ribosum, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;
    
    //status = Mutual_PostOrderPS(msa, T, ribosum, mi, tol, verbose, errbuf);
    //if (status != eslOK) goto ERROR;
  }
  
  // pp validation
  for (i = 0; i < mi->alen-1; i ++) 
    for (j = i+1; j < mi->alen; j ++) {
      status = esl_vec_DValidate(mi->pp[i][j], K*K, tol, errbuf);
      if (status != eslOK) {
	printf("pp[%d][%d]\n", i, j);
	esl_vec_DDump(stdout, mi->pp[i][j], K*K, NULL);
	goto ERROR;
      }
    }

  // ps are the marginals
  for (i = 0; i < mi->alen; i ++) {
    esl_vec_DSet(mi->ps[i], K, 0.0);

    for (j = 0; j < mi->alen; j ++)     
      for (x = 0; x < K; x ++) 
	for (y = 0; y < K; y ++) 
	  if (j != i) mi->ps[i][x] += mi->pp[i][j][IDX(x,y,K)];
      
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
	printf("pp[%d][%d] = ", i, j);
	for (x = 0; x < K; x ++) 
	  for (y = 0; y < K; y ++) {
	    printf(" %f ", mi->pp[i][j][IDX(x,y,K)]);
	  }
	printf("\n");
      }
      esl_vec_DDump(stdout, mi->ps[i], K, NULL);
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

  if (1||verbose) {
    printf("MI[%f,%f] MIa[%f,%f] MIp[%f,%f] MIr[%f,%f] \n", 
	   mi->minMI, mi->maxMI, mi->minMIa, mi->maxMIa, mi->minMIp, mi->maxMIp, mi->minMIr, mi->maxMIr);
    for (i = 0; i < mi->alen-1; i++) 
      for (j = i+1; j < mi->alen; j++) {
	if (mi->MIp->mx[i][j] > 0.4) printf("MI[%d][%d] = %f | %f %f %f | %f %f | %f %f | %f\n", 
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

  mi->threshMI  = 0.0;
  mi->threshMIa = 0.0;
  mi->threshMIp = 0.0;
  mi->threshMIr = 0.0;

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
  int64_t alen = msa->alen;
  int     i, j;
  int     status;

  for (i = 0; i < alen-1; i ++) 
    for (j = i+1; j < alen; j ++) {
      status = mutual_postorder_ppij(i, j, msa, T, ribosum, mi, tol, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
  
  return eslOK;

 ERROR:
  return status;
}

int 
Mutual_PostOrderPS(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i;
  int     status;

  for (i = 0; i < alen; i ++) {
    status = mutual_postorder_psi(i, msa, T, ribosum, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;
  }

  return eslOK;

 ERROR:
  return status;
}

/*---------------- internal functions --------------------- */

static int
mutual_analyze_significant_pairs(int *ct, struct mutual_s *mi, MItype whichmi, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx;
  double       avg;
  double       std;
  double       zscore;
  int          i, j;
  int          status;

  switch(whichmi) {
  case MI:  mtx = mi->MI; break;
  case MIa: mtx = mi->MIa; break; 
  case MIp: mtx = mi->MIp; break; 
  case MIr: mtx = mi->MIr; break; 
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MItype");
  }
  
  for (i = 0; i < mi->alen; i ++) {
    if (ct[i+1] > 0 && ct[i+1] > i+1) {
      avg    = 0.0;
      std    = 0.0;
      zscore = 0.0;

      for (j = 0; j < mi->alen; j ++) {
	if (ct[i+1] != j+1) {       
	  avg += mtx->mx[i][j];
	  std += mtx->mx[i][j] * mtx->mx[i][j];
	}
      }
      avg /= mi->alen-1;
      std /= mi->alen-1;
      std -= avg*avg;
      std = sqrt(std);

      zscore = (mtx->mx[i][ct[i+1]-1] - avg) / std;
      printf("[%d][%d] %f | %f %f %f\n", i, ct[i+1]-1, mtx->mx[i][ct[i+1]-1], avg, std, zscore);
    }
  }  
  return eslOK;
  
 ERROR:
  return status;
}

static int
mutual_analyze_ranking(int *ct, struct mutual_s *mi, MItype whichmi, double thresh, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx;
  double       sen;
  double       ppv;
  double       F;
  int          tf = 0;
  int          f  = 0;
  int          t  = 0;
  int          i, j;
  int          status;

  switch(whichmi) {
  case MI:  mtx = mi->MI; break;
  case MIa: mtx = mi->MIa; break; 
  case MIp: mtx = mi->MIp; break; 
  case MIr: mtx = mi->MIr; break; 
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong MItype");
  }
  
  for (i = 0; i < mi->alen-1; i ++) 
    for (j = i+1; j < mi->alen; j ++) {
      if (mtx->mx[i][j] > thresh)   f  ++;
      if (ct[i+1] == j+1) {         t  ++;
	if (mtx->mx[i][j] > thresh) tf ++;
      }
    }

  sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
  ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
  F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;

  printf("%.5f %d %d %d %.2f %.2f %.2f\n", thresh, tf, t, f, sen, ppv, F);
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
mutual_postorder_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  ESL_STACK     *vs = NULL;   /* node index stack */
  ESL_DMATRIX  **pp = NULL;
  ESL_DMATRIX   *ppl, *ppr;
  ESL_DMATRIX   *CL = NULL;
  ESL_DMATRIX   *CR = NULL;
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
  
  /* allocate the single and pair probs for theinternal nodes */
  nnodes = (T->N > 1)? T->N-1 : T->N;
  dim    = nnodes + T->N;
  ESL_ALLOC(pp, sizeof(ESL_DMATRIX *) * dim);
  for (v = 0; v < dim; v ++)  pp[v] = NULL;
  CL = esl_dmatrix_Create(K2, K2);
  CR = esl_dmatrix_Create(K2, K2);
 
  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL)   { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      //printf("%d %d pp-pop v %d | l %d r %d\n", i,  j, v, T->left[v], T->right[v]);

      which = (T->left[v] <= 0)? -T->left[v] : T->left[v];
      idx   = (T->left[v] <= 0)? nnodes + which : which;
      if (T->left[v] <= 0) {
	pp[idx] = esl_dmatrix_Create(K, K);
	esl_dmatrix_Set(pp[idx], 0.0);
	resi = msa->ax[which][i+1];
	resj = msa->ax[which][j+1];
	//printf("v=%d parent %d pp l time %f %d which %d i %d %d (%d) j %d %d (%d)\n", 
	//       v, T->parent[v], T->ld[v], idx, which, 
	//       i, resi, esl_abc_XIsCanonical(msa->abc, resi), 
	//       j, resj, esl_abc_XIsCanonical(msa->abc, resj));

	if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) {
	  pp[idx]->mx[resi][resj] = 1.0;
	}
	else if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  for (y = 0; y < K; y ++) pp[idx]->mx[resi][y] = 1.0/(double)K;
	}
	else if (esl_abc_XIsCanonical(msa->abc, resj)) {
	  for (x = 0; x < K; x ++) pp[idx]->mx[x][resj] = 1.0/(double)K;
	}
	else {
	  esl_dmatrix_Set(pp[idx], 1.0/(double)K2);
	}
      }
      ppl = pp[idx];

      which = (T->right[v] <= 0)? -T->right[v] : T->right[v];
      idx   = (T->right[v] <= 0)? nnodes + which : which;
      if (T->right[v] <= 0) {
	pp[idx] = esl_dmatrix_Create(K, K);
	esl_dmatrix_Set(pp[idx], 0.0); 
	resi = msa->ax[which][i+1];
	resj = msa->ax[which][j+1];
	//printf("v=%d parent %d pp r time %f %d which %d i %d %d (%d) j %d %d (%d)\n", 
	//       v, T->parent[v], T->rd[v], idx, which, 
	//       i, resi, esl_abc_XIsCanonical(msa->abc, resi), 
	//       j, resj, esl_abc_XIsCanonical(msa->abc, resj));


	if (esl_abc_XIsCanonical(msa->abc, resi) && esl_abc_XIsCanonical(msa->abc, resj)) {
	  pp[idx]->mx[resi][resj] = 1.0;
	}
	else if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  for (y = 0; y < K; y ++) pp[idx]->mx[resi][y] = 1.0/(double)K;
	}
	else if (esl_abc_XIsCanonical(msa->abc, resj)) {
	  for (x = 0; x < K; x ++) pp[idx]->mx[x][resj] = 1.0/(double)K;
	}
	else {
	  esl_dmatrix_Set(pp[idx], 1.0/(double)K2);
	}
      }
      ppr = pp[idx];

     if (ppl != NULL && ppr != NULL) { /* ready to go: calculate ps and pp at the parent node */
       status = ratematrix_CalculateConditionalsFromRate(T->ld[v], ribosum->bprsQ, CL, tol, errbuf, verbose);
       if (status != eslOK) goto ERROR;
       status = ratematrix_CalculateConditionalsFromRate(T->rd[v], ribosum->bprsQ, CR, tol, errbuf, verbose);
       if (status != eslOK) goto ERROR;

       pp[v] = esl_dmatrix_Create(K, K);
	esl_dmatrix_Set(pp[v], 0.0);
	for (x = 0; x < K; x ++) 
	  for (y = 0; y < K; y ++) {
	    for (xl = 0; xl < K; xl ++) 
	      for (yl = 0; yl < K; yl ++) 
		for (xr = 0; xr < K; xr ++) 
		  for (yr = 0; yr < K; yr ++) 
		    pp[v]->mx[x][y] += ppl->mx[xl][yl] * CL->mx[IDX(x,y,K)][IDX(xl,yl,K)] * ppr->mx[xr][yr] * CR->mx[IDX(x,y,K)][IDX(xr,yr,K)];
	  }
	
	/* push parent into stack unless already at the root */
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslFAIL; goto ERROR; }; 
      }
      else if (ppl == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])  != eslOK) { status = eslFAIL; goto ERROR; };
      }
      else if (ppr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v]) != eslOK) { status = eslFAIL; goto ERROR; }
      }
    }
  if (v != 0) ESL_XFAIL(eslFAIL, errbuf, "pp did not transverse tree to the root");
  
  for (x = 0; x < K; x ++) 
    for (y = 0; y < K; y ++) 
      mi->pp[i][j][IDX(x,y,K)] = pp[v]->mx[x][y] * ribosum->bprsM[IDX(x,y,K)];
  esl_vec_DNorm(mi->pp[i][j], K*K);

  for (v = 0; v < dim; v ++) esl_dmatrix_Destroy(pp[v]);
  free(pp);
  esl_stack_Destroy(vs);
  esl_dmatrix_Destroy(CL);
  esl_dmatrix_Destroy(CR);
  return eslOK;

 ERROR:
  if (vs) esl_stack_Destroy(vs);
  for (v = 0; v < dim; v ++) if (pp[v]) free(pp[v]);
  if (pp) free(pp);
  if (CL) esl_dmatrix_Destroy(CL);
  if (CR) esl_dmatrix_Destroy(CR);
  return status;
}

static int    
mutual_postorder_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  ESL_STACK     *vs = NULL;   /* node index stack */
  double       **ps = NULL;
  double        *psl, *psr;
  ESL_DMATRIX   *CL = NULL;
  ESL_DMATRIX   *CR = NULL;
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
  
  /* allocate the single and pair probs for theinternal nodes */
  nnodes = (T->N > 1)? T->N-1 : T->N;
  dim    = nnodes + T->N;
  ESL_ALLOC(ps, sizeof(double *) * dim);
  for (v = 0; v < dim; v ++) ps[v] = NULL;
  CL = esl_dmatrix_Create(K, K);
  CR = esl_dmatrix_Create(K, K);
 
  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      which = (T->left[v] <= 0)? -T->left[v] : T->left[v];
      idx   = (T->left[v] <= 0)? nnodes + which : which;
      if (T->left[v] <= 0) {
	resi = msa->ax[which][i+1];
	ESL_ALLOC(ps[idx], sizeof(double)*K);
	esl_vec_DSet(ps[idx], K, 0.0);
	if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  ps[idx][resi] = 1.0;
	}
	else {
	  esl_vec_DSet(ps[idx], K, 1.0/(double)K);
	}

      }
      psl = ps[idx];
      
      which = (T->right[v] <= 0)? -T->right[v] : T->right[v];
      idx   = (T->right[v] <= 0)? nnodes + which : which;
      if (T->right[v] <= 0) {
	resi = msa->ax[which][i+1];
	ESL_ALLOC(ps[idx], sizeof(double)*K);
	esl_vec_DSet(ps[idx], K, 0.0);
	esl_vec_DSet(ps[idx], K, 0.0);
	if (esl_abc_XIsCanonical(msa->abc, resi)) {
	  ps[idx][resi] = 1.0;
	}
	else {
	  esl_vec_DSet(ps[idx], K, 1.0/(double)K);
	}
      }
      psr = ps[idx];
      
      if (psl != NULL && psr != NULL) { /* ready to go: calculate ps at the parent node */
	status = ratematrix_CalculateConditionalsFromRate(T->ld[v], ribosum->xrnaQ, CL, tol, errbuf, verbose);
	if (status != eslOK) goto ERROR;
	status = ratematrix_CalculateConditionalsFromRate(T->rd[v], ribosum->xrnaQ, CR, tol, errbuf, verbose);
	if (status != eslOK) goto ERROR;
	
	ESL_ALLOC(ps[v], sizeof(double)*K);
	esl_vec_DSet(ps[v], K, 0.0);
	for (x = 0; x < K; x ++) {
	  for (yl = 0; yl < K; yl ++)
	    for (yr = 0; yr < K; yr ++)
	      ps[v][x] += psl[yl] * CL->mx[x][yl] * psr[yr] * CR->mx[x][yr];
	}
	
	/* push node into stack unless already at the root */
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslFAIL; goto ERROR; }; 
      }
      else if (psl == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])  != eslOK) { status = eslFAIL; goto ERROR; };
      }
      else if (psr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v]) != eslOK) { status = eslFAIL; goto ERROR; }
      }
    }
  if (v != 0) ESL_XFAIL(eslFAIL, errbuf, "ps did not transverse tree to the root");
  for (x = 0; x < K; x ++) mi->ps[i][x] = ps[v][x] * ribosum->xrnaM[x];
  esl_vec_DNorm(mi->ps[i], K);

  for (v = 0; v < dim; v ++) free(ps[v]);
  free(ps);
  if (CL) esl_dmatrix_Destroy(CL);
  if (CR) esl_dmatrix_Destroy(CR);
  esl_stack_Destroy(vs);
  return eslOK;

 ERROR:
  if (vs) esl_stack_Destroy(vs);
  for (v = 0; v < dim; v ++) if (ps[v]) free(ps[v]);
  if (ps) free(ps);
  if (CL) esl_dmatrix_Destroy(CL);
  if (CR) esl_dmatrix_Destroy(CR);
  return status;
}


