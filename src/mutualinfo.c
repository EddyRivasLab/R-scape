/* mutualinfo.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "mutualinfo.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"

static int mutual_post_order_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf);
static int mutual_post_order_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf);

int                 
Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf)
{
  int         status;
  
  status = Mutual_PostOrderPP(msa, T, ribosum, mi, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  
  status = Mutual_PostOrderPS(msa, T, ribosum, mi, verbose, errbuf);
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
Mutual_Create(int64_t alen)
{
  struct mutual_s *mi = NULL;
  int              i, j;
  int              status;

  ESL_ALLOC(mi, sizeof(struct mutual_s));
  mi->alen = alen;
  
  mi->pp  = esl_dmatrix_Create(alen, alen);
  ESL_ALLOC(mi->ps, sizeof(double) * alen);

  mi->MI  = esl_dmatrix_Create(alen, alen);
  mi->MIa = esl_dmatrix_Create(alen, alen);
  mi->MIp = esl_dmatrix_Create(alen, alen);
  mi->MIr = esl_dmatrix_Create(alen, alen);
  ESL_ALLOC(mi->H, sizeof(double) * alen);
 
  for (i = 0; i < alen; i++) {
    mi->ps[i] = 0.0;
    mi->H[i]  = 0.0;
    for (j = 0; j < alen; j++) {
      mi->pp->mx[i][j]  = 0.0;
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
  if (mi) {
    esl_dmatrix_Destroy(mi->pp);
    esl_dmatrix_Destroy(mi->MI);
    esl_dmatrix_Destroy(mi->MIa);
    esl_dmatrix_Destroy(mi->MIp);
    esl_dmatrix_Destroy(mi->MIr);
    free(mi->ps);
    free(mi->H);
    free(mi);
  }
}


int 
Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i, j;
  int     status;

  for (i = 0; i < alen; i ++)
    for (j = i; j < alen; j ++) {
      status = mutual_post_order_ppij(i, j, msa, T, ribosum, mi, verbose, errbuf);
    }
  
  return eslOK;

 ERROR:
  return status;
}

int 
Mutual_PostOrderPS(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf)
{
  int64_t alen = msa->alen;
  int     i;
  int     status;

  for (i = 0; i < alen; i ++)
    status = mutual_post_order_psi(i, msa, T, ribosum, mi, verbose, errbuf);
  
  return eslOK;

 ERROR:
  return status;
}

/*---------------- internal functions --------------------- */


int 
mutual_post_order_ppij(int i, int j, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf)
{
  ESL_STACK     *vs = NULL;   /* node index stack */
  ESL_DMATRIX  **pp = NULL;
  ESL_DMATRIX   *ppl, *ppr;
  int            dim;
  int            K = msa->abc->K;
  int            nnodes;
  int            v;
  int            idx;
  int            status;
  
  /* allocate the single and pair probs for theinternal nodes */
  nnodes = (T->N > 1)? T->N-1 : T->N;
  dim    = nnodes + T->N;
  ESL_ALLOC(pp, sizeof(ESL_DMATRIX *) * dim);
  for (v = 0; v < dim; v ++)  pp[v] = NULL;
 
  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      if (T->left[v] <= 0)  idx = nnodes - T->left[v];
      else                  idx = T->left[v];  
      pp[idx] = esl_dmatrix_Create(K, K);
      ppl = pp[idx];

      if (T->right[v] <= 0) idx = nnodes - T->right[v];
      else                  idx = T->right[v];  
      pp[idx] = esl_dmatrix_Create(K, K);
      ppr = pp[idx];

     if (ppl != NULL && ppr != NULL) { /* ready to go: calculate ps and pp at the parent node */
	
	/* push parent into stack unless already at the root */
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslFAIL; goto ERROR; }; 
      }
      else if (ppl == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])   != eslOK) { status = eslFAIL; goto ERROR; };
      }
      else if (ppr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v])  != eslOK) { status = eslFAIL; goto ERROR; }
      }
    }

  mi->pp->mx[i][j] = pp[0]->mx[i][j];

  for (v = 0; v < dim; v ++) esl_dmatrix_Destroy(pp[v]);
  free(pp);
  esl_stack_Destroy(vs);

  return eslOK;

 ERROR:
  if (vs) esl_stack_Destroy(vs);
  for (v = 0; v < dim; v ++) if (pp[v]) free(pp[v]);
  if (pp)    free(pp);
  return status;
}

static int    
mutual_post_order_psi(int i, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf)
{
  ESL_STACK  *vs = NULL;   /* node index stack */
  double    **ps = NULL;
  double     *psl, *psr;
  int         dim;
  int         K = msa->abc->K;
  int         nnodes;
  int         v;
  int         idx;
  int         status;
  
  /* allocate the single and pair probs for theinternal nodes */
  nnodes = (T->N > 1)? T->N-1 : T->N;
  dim    = nnodes + T->N;
  ESL_ALLOC(ps, sizeof(double *) * dim);
  for (v = 0; v < dim; v ++) ps[v] = NULL;
 
  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      if (T->left[v] <= 0) idx = nnodes - T->left[v];
      else                 idx = T->left[v];  
      ESL_ALLOC(ps[idx], sizeof(double) * K);
      psl = ps[idx];

      if (T->right[v] <= 0) idx = nnodes - T->right[v];
      else                  idx = T->right[v];  
      ESL_ALLOC(ps[idx], sizeof(double) * K);
      psr = ps[idx];
      
      if (psl != NULL && psr != NULL) { /* ready to go: calculate ps at the parent node */
	
	/* push parent into stack unless already at the root */
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslFAIL; goto ERROR; }; 
      }
      else if (psl == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])   != eslOK) { status = eslFAIL; goto ERROR; };
      }
      else if (psr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v])  != eslOK) { status = eslFAIL; goto ERROR; }
      }
    }

  mi->ps[i] = ps[0][i];

  for (v = 0; v < dim; v ++) free(ps[v]);
  free(ps);
  esl_stack_Destroy(vs);

  return eslOK;

 ERROR:
  if (vs) esl_stack_Destroy(vs);
  for (v = 0; v < dim; v ++) if (ps[v]) free(ps[v]);
  if (ps) free(ps);
  return status;
}


