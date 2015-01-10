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

static int mutual_post_order_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, 
				  int verbose, char *errbuf);
static int mutual_post_order_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				 double tol, int verbose, char *errbuf);

int                 
Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int         status;
  
  status = Mutual_PostOrderPP(msa, T, ribosum, mi, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  status = Mutual_PostOrderPS(msa, T, ribosum, mi, tol, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  // MI

  // MIa

  // MIp

  // MIr

  return eslOK;

 ERROR:
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
      ESL_ALLOC(mi->pp[i][j], sizeof(double ) * K2);
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
Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i, j;
  int     status;

  for (i = 0; i < alen; i ++)
    for (j = i; j < alen; j ++) {
      status = mutual_post_order_ppij(i, j, msa, T, ribosum, mi, tol, verbose, errbuf);
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
    status = mutual_post_order_psi(i, msa, T, ribosum, mi, tol, verbose, errbuf);
    if (status != eslOK) goto ERROR;
  }

  return eslOK;

 ERROR:
  return status;
}

/*---------------- internal functions --------------------- */


int 
mutual_post_order_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
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
  int            xx, yy;
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
      printf("%d %d pp-pop v %d | l %d r %d\n", i,  j, v, T->left[v], T->right[v]);

      which = (T->left[v] <= 0)? -T->left[v] : T->left[v];
      idx   = (T->left[v] <= 0)? nnodes + which : which;
      if (T->left[v] <= 0) {
	pp[idx] = esl_dmatrix_Create(K, K);
	esl_dmatrix_Set(pp[idx], 0.0);
	resi = msa->ax[which][i+1];
	resj = msa->ax[which][j+1];
	printf("v=%d parent %d pp l time %f %d which %d i %d %d (%d) j %d %d (%d)\n", 
	       v, T->parent[v], T->ld[v], idx, which, 
	       i, resi, esl_abc_XIsResidue(msa->abc, resi), 
	       j, resj, esl_abc_XIsResidue(msa->abc, resj));

	if (esl_abc_XIsResidue(msa->abc, resi) && esl_abc_XIsResidue(msa->abc, resj)) {
	  pp[idx]->mx[resi][resj] = 1.0;
	}
	else if (esl_abc_XIsResidue(msa->abc, resi)) {
	  for (y = 0; y < K; y ++) pp[idx]->mx[resi][y] = 1.0/(double)K;
	}
	else if (esl_abc_XIsResidue(msa->abc, resj)) {
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
	printf("v=%d parent %d pp r time %f %d which %d i %d %d (%d) j %d %d (%d)\n", 
	       v, T->parent[v], T->rd[v], idx, which, 
	       i, resi, esl_abc_XIsResidue(msa->abc, resi), 
	       j, resj, esl_abc_XIsResidue(msa->abc, resj));


	if (esl_abc_XIsResidue(msa->abc, resi) && esl_abc_XIsResidue(msa->abc, resj)) {
	  pp[idx]->mx[resi][resj] = 1.0;
	}
	else if (esl_abc_XIsResidue(msa->abc, resi)) {
	  for (y = 0; y < K; y ++) pp[idx]->mx[resi][y] = 1.0/(double)K;
	}
	else if (esl_abc_XIsResidue(msa->abc, resj)) {
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

	    for (xx = 0; xx < K; xx ++) 
	      for (yy = 0; yy < K; yy ++) 
		pp[v]->mx[x][y] += ppl->mx[xx][yy] * CL->mx[IDX(x,y,K)][IDX(xx,yy,K)];
	    for (xx = 0; xx < K; xx ++) 
	      for (yy = 0; yy < K; yy ++) 
		pp[v]->mx[x][y] += ppr->mx[xx][yy] * CR->mx[IDX(x,y,K)][IDX(xx,yy,K)];
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
  for (xx = 0; xx < K; xx ++) 
    for (yy = 0; yy < K; yy ++) 
      mi->pp[i][j][IDX(xx,yy,K)] = pp[v]->mx[xx][yy];

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
mutual_post_order_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, double tol, int verbose, char *errbuf)
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
  int            x, y;
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
	if (esl_abc_XIsResidue(msa->abc, resi)) {
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
	if (esl_abc_XIsResidue(msa->abc, resi)) {
	  ps[idx][resi] = 1.0;
	}
	else {
	  esl_vec_DSet(ps[idx], K, 1.0/(double)K);
	}
      }
      psr = ps[idx];
      
      if (psl != NULL && psr != NULL) { /* ready to go: calculate ps at the parent node */
	status = ratematrix_CalculateConditionalsFromRate(T->ld[v], ribosum->prnaQ, CL, tol, errbuf, verbose);
	if (status != eslOK) goto ERROR;
	status = ratematrix_CalculateConditionalsFromRate(T->rd[v], ribosum->prnaQ, CR, tol, errbuf, verbose);
	if (status != eslOK) goto ERROR;
	
	ESL_ALLOC(ps[v], sizeof(double)*K);
	esl_vec_DSet(ps[v], K, 0.0);
	for (x = 0; x < K; x ++) {
	  for (y = 0; y < K; y ++)
	    ps[v][x] += psl[y] * CL->mx[x][y];
	  for (y = 0; y < K; y ++)
	    ps[v][x] += psr[y] * CR->mx[x][y];
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
  esl_vec_DCopy(ps[v], mi->pp[i], K);

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


