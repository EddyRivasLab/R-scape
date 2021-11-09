/* cacofold.c */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covgrammars.h"
#include "cacofold.h"
#include "correlators.h"
#include "contactmap.h"
#include "e2_profilesq.h"
#include "logsum.h"
#include "structure.h"


// Index for (i,j,L)
//     for i=0 < L-1; j=i+1 < L
//
//     IDX(i,j,L) = \sum_{k=0}^{i-1}\sum_{l=k+1}^{L-1}          + j - i - 1
//
//                = \sum_{k=0}^{i-1} (L - 1 - k)                + j - i - 1
//
//                = (L-1)*\sum_{k=0}^{i-1} - \sum_{k=0}^{i-1} k + j - i - 1
//
//                = (L-1)*i                - i*(i-1)/2          + j - i - 1
//
#define INDEX(i, j, L) ( ((L) - 1)*(i) - (i)*((i)-1)/2 + (j) - (i) - 1 )

static int   dp_recursion_mea_cyk              (FOLDPARAM *foldparam, G6Xparam *p,  POST *post, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose);
static int   dp_recursion_g6x_cyk              (FOLDPARAM *foldparam, G6Xparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose);
static int   dp_recursion_g6x_inside           (FOLDPARAM *foldparam, G6Xparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_g6x_outside          (FOLDPARAM *foldparam, G6Xparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_g6x_posterior_single (FOLDPARAM *foldparam, G6Xparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx,
						int w, int j, POST *post, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_g6x_posterior_pair   (FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx,
						int w, int j, int d, POST *post, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_g6xs_cyk             (FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX  *cyk,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose);
static int   dp_recursion_g6xs_inside          (FOLDPARAM *foldparam, G6XSparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_g6xs_outside         (FOLDPARAM *foldparam, G6XSparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_g6xs_posterior_single(FOLDPARAM *foldparam, G6XSparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx,
						int w, int j, POST *post, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_g6xs_posterior_pair  (FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx,
						int w, int j, int d, POST *post, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_rbg_cyk              (FOLDPARAM *foldparam, RBGparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *cyk,
					        int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose);
static int   dp_recursion_rbg_inside           (FOLDPARAM *foldparam, RBGparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_rbg_outside          (FOLDPARAM *foldparam, RBGparam *p,  PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *omx, RBG_MX *imx,
						int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose);
static int   dp_recursion_rbg_posterior_pair   (FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx, RBG_MX *omx,
						int w, int j, int d, POST *post, int *ret_nneg, char *errbuf, int verbose);
static int   allow_bpair(double power_thresh,  double neg_eval_thresh, int hloop_min, int i, int j, int L, int *covct,
			 COVLIST *exclude, SPAIR *spair, int *ret_nneg);
static int   force_bpair(int i, int j, int *covct);
static int   allow_hairpin(int hloop_min, int i, int j, int L, int *covct);
static int   allow_loop(int i, int j, int *covct);
static int   allow_single(int i, int *covct);
static SCVAL emitsc_stck     (int i, int j, int L, ESL_DSQ *dsq, SCVAL e_pair[NP], SCVAL e_stck[NP][NP]);
static SCVAL emitsc_stck_prof(int i, int j, int L, PSQ     *psq, SCVAL e_pair[NP], SCVAL e_stck[NP][NP]);
static SCVAL emitsc_pair     (int i, int j,        ESL_DSQ *dsq, SCVAL e_pair[NP]);
static SCVAL emitsc_pair_prof(int i, int j,        PSQ     *psq, SCVAL e_pair[NP]);
static SCVAL emitsc_sing     (int i,               ESL_DSQ *dsq, SCVAL e_sing[NB]);
static SCVAL emitsc_sing_prof(int i,               PSQ     *psq, SCVAL e_sing[NB]);
static SCVAL score_loop_hairpin      (int i, int j, RBGparam *p, ESL_DSQ *dsq);
static SCVAL score_loop_hairpin_prof (int i, int j, RBGparam *p, PSQ *psq);
static SCVAL score_loop_bulge        (int i, int j, RBGparam *p, ESL_DSQ *dsq);
static SCVAL score_loop_bulge_prof   (int i, int j, RBGparam *p, PSQ *psq);
static SCVAL score_loop_intloop      (int i, int j, RBGparam *p, ESL_DSQ *dsq);
static SCVAL score_loop_intloop_prof (int i, int j, RBGparam *p, PSQ *psq);
static int   segment_remove_gaps     (int i, int j, ESL_DSQ *dsq);
static int   segment_remove_gaps_prof(int i, int j, PSQ *psq);
static int   vec_SCVAL_LogNorm(SCVAL *vec, int n);
static int   dvec_SCVAL_LogNorm(int n1, int n2, SCVAL dvec[n1][n2]);

/* G6X/G6XS
 *----------------------------------------------------------
 *   S -> LS   | L   | epsilon
 *   L -> aFa' | aa' | a
 *   F -> aFa' | aa' | LS
 *
 *
 * Basic Grammar (RBG)
 *-----------------------------------------------------------
 *  S  -> a S     | F0 S    | e
 *  F0 -> a F5 a' | a P a'  | aa'
 *  F5 -> a F5 a' | a P a'  | aa'
 *  P  -> m..m    | m..m F0 | F0 m..m | m..m F0 m..m | M1 M
 *  M  -> M1 M    | R
 *  R  ->    R a  | M1
 *  M1 -> a M1    | F0
 *
 */

/* Inputs are:
 *
 *        covct[L+1]  - covct[i]  >  0 a covarying pair forced to  basepair in structure s (the covariation skeleton of s)
 *                      covct[i]  =  0 unrestricted
 *
 *        exclude     - a CLIST with those covarying pairs forced to remain unpaired in this strcture.
 *
 * Output: fills  ct[L+1]  - A complete structure in ct format
 *
 */
int
CACO_CYK(ESL_RANDOMNESS *r, enum grammar_e G, FOLDPARAM *foldparam, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6Xparam  *g6p  = NULL;
  G6XSparam *g6sp = NULL;
  RBGparam *rbgp = NULL;
  int       status;

  /* get the grammar parameters and run the corresponding CYK */
  switch(G) {
  case G6X:
    /* Transfer scores from static built-in storage */
    status = CACO_G6X_GetParam(&g6p, errbuf, verbose);
    if (status != eslOK) goto ERROR;    
    status = CACO_G6X_CYK(r, foldparam, g6p, psq, spair, covct,  exclude, ct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case G6XS:
    status = CACO_G6XS_GetParam(&g6sp, errbuf, verbose);
    if (status != eslOK) goto ERROR; 
    status = CACO_G6XS_CYK(r, foldparam, g6sp, psq, spair, covct, exclude, ct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case RBG:
    status = CACO_RBG_GetParam(&rbgp, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    status = CACO_RBG_CYK(r, foldparam, rbgp, psq, spair, covct, exclude, ct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  }

  if (g6p)  free(g6p);
  if (g6sp) free(g6sp);
  if (rbgp) free(rbgp);
  return eslOK;

 ERROR:
  if (g6p)  free(g6p);
  if (g6sp) free(g6sp);
  if (rbgp) free(rbgp);
  return status;
}

int
CACO_DECODING(ESL_RANDOMNESS *r, enum grammar_e G, FOLDPARAM *foldparam, PSQ *psq, SPAIR *spair, int *covct,  COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  POST      *post = NULL;   /* the posterior probabilities */
  G6Xparam  *g6p  = NULL;
  G6XSparam *g6sp = NULL;
  RBGparam  *rbgp = NULL;
  int        status;

  post = POST_Create(psq->n);
  if (post == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_G6X_Decoding() allocation error\n");

  /* get the grammar parameters and run the corresponding DECODING */
  switch(G) {
  case G6X:
    /* Transfer scores from static built-in storage */
    status = CACO_G6X_GetParam(&g6p, errbuf, verbose);
    if (status != eslOK) goto ERROR;    
    status = CACO_G6X_DECODING(r, foldparam, g6p, psq, spair, covct, exclude, post, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case G6XS:
    status = CACO_G6XS_GetParam(&g6sp, errbuf, verbose);
    if (status != eslOK) goto ERROR; 
    status = CACO_G6XS_DECODING(r,foldparam, g6sp, psq, spair, covct, exclude, post, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case RBG:
    status = CACO_RBG_GetParam(&rbgp, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    status = CACO_RBG_DECODING(r, foldparam, rbgp, psq, spair, covct, exclude, post, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  }

  /* MEA (Maximal Expected Accuracy */
  if ((status = CACO_MEA(r, foldparam, post, spair, covct, exclude, ct, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;

  if (g6p)  free(g6p);
  if (g6sp) free(g6sp);
  if (rbgp) free(rbgp);
  if (post) POST_Destroy(post);
  return eslOK;

 ERROR:
  if (g6p)  free(g6p);
  if (g6sp) free(g6sp);
  if (rbgp) free(rbgp);
  if (post) POST_Destroy(post);
  return status;
}


int
CACO_MEA(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, POST *post, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6Xparam  *meap = NULL;
  G6X_MX    *gmx  = NULL;
  int        status;

  gmx = G6XMX_Create(post->L);
  
  CACO_G6X_MEA_GetParam(&meap, foldparam->gamma, errbuf, verbose);

  /* Fill the cyk matrix */
  if ((status = CACO_MEA_Fill_CYK        (foldparam, meap, post, spair, covct, exclude, gmx, ret_sc, errbuf, verbose))  != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = CACO_MEA_Traceback_CYK(r, foldparam, meap, post, spair, covct, exclude, gmx, ct,      errbuf, verbose)) != eslOK) goto ERROR;
  
  free(meap);
  G6XMX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (meap) free(meap);
  if (gmx)  G6XMX_Destroy(gmx);
  return status;
}


int
CACO_G6X_GetParam(G6Xparam **ret_p, char *errbuf, int verbose)
{
  G6Xparam *p = NULL;
  int       x;
  int       status;

 ESL_ALLOC(p, sizeof(G6Xparam));

  p->t1[0] = G6X_PRELOADS_TrATrBTrB.t1[0];
  p->t1[1] = G6X_PRELOADS_TrATrBTrB.t1[1];
  p->t1[2] = G6X_PRELOADS_TrATrBTrB.t1[2];
  p->t2[0] = G6X_PRELOADS_TrATrBTrB.t2[0];
  p->t2[1] = G6X_PRELOADS_TrATrBTrB.t2[1];
  p->t2[2] = G6X_PRELOADS_TrATrBTrB.t2[2];
  p->t3[0] = G6X_PRELOADS_TrATrBTrB.t3[0];
  p->t3[1] = G6X_PRELOADS_TrATrBTrB.t3[1];
  p->t3[2] = G6X_PRELOADS_TrATrBTrB.t3[2];

  for (x = 0; x < NB; x ++) p->e_sing[x] = G6X_PRELOADS_TrATrBTrB.e_sing[x];
  for (x = 0; x < NP; x ++) p->e_pair[x] = G6X_PRELOADS_TrATrBTrB.e_pair[x];

  // renormalize, just in case
  vec_SCVAL_LogNorm(p->t1, 3);
  vec_SCVAL_LogNorm(p->t2, 3);
  vec_SCVAL_LogNorm(p->t3, 3);
  vec_SCVAL_LogNorm(p->e_sing, NB);
  vec_SCVAL_LogNorm(p->e_pair, NP);

  *ret_p = p;
  return eslOK;

 ERROR:
  if (p) free(p);
  return status;
}


int
CACO_G6XS_GetParam(G6XSparam **ret_p, char *errbuf, int verbose)
{
 G6XSparam *p = NULL;
 int        x, y;
 int        status;

 ESL_ALLOC(p, sizeof(G6XSparam));

  p->t1[0] = G6XS_PRELOADS_TrATrBTrB.t1[0];
  p->t1[1] = G6XS_PRELOADS_TrATrBTrB.t1[1];
  p->t1[2] = G6XS_PRELOADS_TrATrBTrB.t1[2];
  p->t2[0] = G6XS_PRELOADS_TrATrBTrB.t2[0];
  p->t2[1] = G6XS_PRELOADS_TrATrBTrB.t2[1];
  p->t2[2] = G6XS_PRELOADS_TrATrBTrB.t2[2];
  p->t3[0] = G6XS_PRELOADS_TrATrBTrB.t3[0];
  p->t3[1] = G6XS_PRELOADS_TrATrBTrB.t3[1];
  p->t3[2] = G6XS_PRELOADS_TrATrBTrB.t3[2];

  for (x = 0; x < NB; x ++)   p->e_sing[x]    = G6XS_PRELOADS_TrATrBTrB.e_sing[x];
  for (x = 0; x < NP; x ++)   p->e_pair[x]    = G6XS_PRELOADS_TrATrBTrB.e_pair[x];
  for (x = 0; x < NP; x ++) 
    for (y = 0; y < NP; y ++) p->e_stck[x][y] = G6XS_PRELOADS_TrATrBTrB.e_stck[x][y];

  // renormalize, just in case
  vec_SCVAL_LogNorm(p->t1, 3);
  vec_SCVAL_LogNorm(p->t2, 3);
  vec_SCVAL_LogNorm(p->t3, 3);
  vec_SCVAL_LogNorm(p->e_sing, NB);
  vec_SCVAL_LogNorm(p->e_pair, NP);
  for (x = 0; x < NP; x ++)
    vec_SCVAL_LogNorm(p->e_stck[x], NP);

  *ret_p = p;
  return eslOK;

 ERROR:
  if (p) free(p);
  return status;

  return eslOK;
}

int
CACO_RBG_GetParam(RBGparam **ret_p, char *errbuf, int verbose)
{
  RBGparam *p = NULL;
  int       x, y;
  int       l, l1, l2;
  int       status;
  
  ESL_ALLOC(p, sizeof(RBGparam));

  p->tS[0] = RBG_PRELOADS_TrATrBTrB.tS[0];
  p->tS[1] = RBG_PRELOADS_TrATrBTrB.tS[1];
  p->tS[2] = RBG_PRELOADS_TrATrBTrB.tS[2];
  
  p->tF0[0] = RBG_PRELOADS_TrATrBTrB.tF0[0];
  p->tF0[1] = RBG_PRELOADS_TrATrBTrB.tF0[1];
  p->tF0[2] = RBG_PRELOADS_TrATrBTrB.tF0[2];
  
  p->tF5[0] = RBG_PRELOADS_TrATrBTrB.tF5[0];
  p->tF5[1] = RBG_PRELOADS_TrATrBTrB.tF5[1];
  p->tF5[2] = RBG_PRELOADS_TrATrBTrB.tF5[2];
  
  p->tP[0] = RBG_PRELOADS_TrATrBTrB.tP[0];
  p->tP[1] = RBG_PRELOADS_TrATrBTrB.tP[1];
  p->tP[2] = RBG_PRELOADS_TrATrBTrB.tP[2];
  p->tP[3] = RBG_PRELOADS_TrATrBTrB.tP[3];
  p->tP[4] = RBG_PRELOADS_TrATrBTrB.tP[4];
  
  p->tM[0]  = RBG_PRELOADS_TrATrBTrB.tM[0];
  p->tM[1]  = RBG_PRELOADS_TrATrBTrB.tM[1];

  p->tR[0]  = RBG_PRELOADS_TrATrBTrB.tR[0];
  p->tR[1]  = RBG_PRELOADS_TrATrBTrB.tR[1];

  p->tM1[0] = RBG_PRELOADS_TrATrBTrB.tM1[0];
  p->tM1[1] = RBG_PRELOADS_TrATrBTrB.tM1[1];
  
  for (x = 0; x < NB;  x ++) {
    p->e_sing[x]    = RBG_PRELOADS_TrATrBTrB.e_sing[x];
    p->e_sing_l1[x] = RBG_PRELOADS_TrATrBTrB.e_sing_l1[x];
    p->e_sing_l2[x] = RBG_PRELOADS_TrATrBTrB.e_sing_l2[x];
    p->e_sing_l3[x] = RBG_PRELOADS_TrATrBTrB.e_sing_l3[x];
  }
  
  for (x = 0; x < NP; x ++) {
    p->e_pair1[x] = RBG_PRELOADS_TrATrBTrB.e_pair1[x];
    p->e_pair2[x] = RBG_PRELOADS_TrATrBTrB.e_pair2[x];
    
    for (y = 0; y < NP; y ++) {
      p->e_stck1[x][y] = RBG_PRELOADS_TrATrBTrB.e_stck1[x][y];
      p->e_stck2[x][y] = RBG_PRELOADS_TrATrBTrB.e_stck2[x][y];
    }
  }

  for (l    = 0; l  < MAXLOOP_H; l ++) p->l1[l] = RBG_PRELOADS_TrATrBTrB.l1[l];
  for (l    = 0; l  < MAXLOOP_B; l ++) p->l2[l] = RBG_PRELOADS_TrATrBTrB.l2[l];
  for (l1   = 0; l1 < MAXLOOP_I; l1 ++) 
    for (l2 = 0; l2 < MAXLOOP_I; l2 ++) 
      p->l3[l1][l2] = RBG_PRELOADS_TrATrBTrB.l3[l1][l2];

  // renormalize, just in case
  vec_SCVAL_LogNorm(p->tS,  3);
  vec_SCVAL_LogNorm(p->tF0, 3);

  vec_SCVAL_LogNorm(p->tF5, 3);
  vec_SCVAL_LogNorm(p->tP,  5);
  vec_SCVAL_LogNorm(p->tM,  2);
  vec_SCVAL_LogNorm(p->tR,  2);
  vec_SCVAL_LogNorm(p->tM1, 2);
  vec_SCVAL_LogNorm(p->e_sing,    NB);
  vec_SCVAL_LogNorm(p->e_sing_l1, NB);
  vec_SCVAL_LogNorm(p->e_sing_l2, NB);
  vec_SCVAL_LogNorm(p->e_sing_l3, NB);
  vec_SCVAL_LogNorm(p->e_pair1,   NP);
  vec_SCVAL_LogNorm(p->e_pair2,   NP);
  for (x = 0; x < NP; x ++)  {
    vec_SCVAL_LogNorm(p->e_stck1[x], NP);
    vec_SCVAL_LogNorm(p->e_stck2[x], NP);
  }
  vec_SCVAL_LogNorm(p->l1, MAXLOOP_H);
  vec_SCVAL_LogNorm(p->l2, MAXLOOP_B);
  dvec_SCVAL_LogNorm(MAXLOOP_I, MAXLOOP_I, p->l3);
  
  *ret_p = p;
  return eslOK;

 ERROR:
  if (p) free(p);
  return status;

  return eslOK;
}


int
CACO_G6X_CYK(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, G6Xparam  *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6X_MX *gmx = NULL;
  int    status;

  gmx = G6XMX_Create(psq->n);

  /* Fill the cyk matrix */
  if ((status = CACO_G6X_Fill_CYK        (foldparam, p, psq, spair, covct, exclude, gmx, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = CACO_G6X_Traceback_CYK(r, foldparam, p, psq, spair, covct, exclude, gmx, ct,     errbuf, verbose)) != eslOK) goto ERROR;
  
  G6XMX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) G6XMX_Destroy(gmx);
  return status;
}

int
CACO_G6X_DECODING(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, POST *post, char *errbuf, int verbose) 
{
  G6X_MX *imx  = NULL;
  G6X_MX *omx  = NULL;
  SCVAL   isc;
  SCVAL   osc;
  int     status;

  imx = G6XMX_Create(psq->n);
  if (imx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_G6X_Inside()  allocation error\n");
  omx = G6XMX_Create(psq->n);
  if (omx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_G6X_Outside() allocation error\n");

  /* Inside algorithm */  
  if ((status = CACO_G6X_Inside   (foldparam, p, psq, spair, covct, exclude, imx,      &isc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Outside algorithm */
  if ((status = CACO_G6X_Outside  (foldparam, p, psq, spair, covct, exclude, omx, imx, &osc, errbuf, verbose)) != eslOK) goto ERROR;
  /* Posterior decoding */
  if ((status = CACO_G6X_Posterior(foldparam, p, psq, spair, covct, exclude, imx, omx, post, errbuf, verbose)) != eslOK) goto ERROR;
  
  G6XMX_Destroy(imx);
  G6XMX_Destroy(omx);
  return eslOK;

 ERROR:
  if (imx)  G6XMX_Destroy(imx);
  if (omx)  G6XMX_Destroy(omx);
  return status;
}

int
CACO_G6XS_CYK(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, G6XSparam  *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6X_MX *gmx = NULL;
  int    status;

  gmx = G6XMX_Create(psq->n);
  if (gmx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_G6X_CYK() allocation error\n");

  /* Fill the cyk matrix */
  if ((status = CACO_G6XS_Fill_CYK        (foldparam, p, psq, spair, covct, exclude, gmx, ret_sc, errbuf, verbose))  != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = CACO_G6XS_Traceback_CYK(r, foldparam, p, psq, spair, covct, exclude, gmx, ct,     errbuf, verbose)) != eslOK) goto ERROR;
  
  G6XMX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) G6XMX_Destroy(gmx);
  return status;
}

int
CACO_G6XS_DECODING(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, POST *post, char *errbuf, int verbose) 
{
  G6X_MX *imx  = NULL;
  G6X_MX *omx  = NULL;
  SCVAL   isc;
  SCVAL   osc;
  int     status;

  imx = G6XMX_Create(psq->n);
  if (imx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_G6XS_Inside()  allocation error\n");
  omx = G6XMX_Create(psq->n);
  if (omx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_G6XS_Outside() allocation error\n");

  /* Inside algorithm */  
  if ((status = CACO_G6XS_Inside   (foldparam, p, psq, spair, covct, exclude, imx,      &isc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Outside algorithm */
  if ((status = CACO_G6XS_Outside  (foldparam, p, psq, spair, covct, exclude, omx, imx, &osc, errbuf, verbose)) != eslOK) goto ERROR;
  /* Posterior decoding */
  if ((status = CACO_G6XS_Posterior(foldparam, p, psq, spair, covct, exclude, imx, omx, post, errbuf, verbose)) != eslOK) goto ERROR;
  
  G6XMX_Destroy(imx);
  G6XMX_Destroy(omx);
  return eslOK;

 ERROR:
  if (imx)  G6XMX_Destroy(imx);
  if (omx)  G6XMX_Destroy(omx);
  return status;
}

int
CACO_RBG_CYK(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  RBG_MX *gmx = NULL;
  int     status;

  gmx = RBGMX_Create(psq->n);
  if (gmx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_RBG_CYK() allocation error\n");

  /* Fill the cyk matrix */
  if ((status = CACO_RBG_Fill_CYK        (foldparam, p, psq, spair, covct, exclude, gmx, ret_sc, errbuf, verbose))  != eslOK) goto ERROR;
  /* Report a traceback */
  if ((status = CACO_RBG_Traceback_CYK(r, foldparam, p, psq, spair, covct, exclude, gmx, ct,     errbuf, verbose)) != eslOK) goto ERROR;

  RBGMX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) RBGMX_Destroy(gmx);
  return status;
}

int
CACO_RBG_DECODING(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, POST *post, char *errbuf, int verbose) 
{
  RBG_MX *imx  = NULL;
  RBG_MX *omx  = NULL;
  SCVAL   isc;
  SCVAL   osc;
  int     status;

  imx = RBGMX_Create(psq->n);
  if (imx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_RBG_Inside()  allocation error\n");
  omx = RBGMX_Create(psq->n);
  if (omx == NULL) ESL_XFAIL(eslFAIL, errbuf, "CACO_RBG_Outside() allocation error\n");

  /* Inside algorithm */  
  if ((status = CACO_RBG_Inside   (foldparam, p, psq, spair, covct, exclude, imx,      &isc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Outside algorithm */
  if ((status = CACO_RBG_Outside  (foldparam, p, psq, spair, covct, exclude, omx, imx, &osc, errbuf, verbose)) != eslOK) goto ERROR;
  /* Posterior decoding */
  if ((status = CACO_RBG_Posterior(foldparam, p, psq, spair, covct, exclude, imx, omx, post, errbuf, verbose)) != eslOK) goto ERROR;
  
  RBGMX_Destroy(imx);
  RBGMX_Destroy(omx);
  return eslOK;

 ERROR:
  if (imx)  RBGMX_Destroy(imx);
  if (omx)  RBGMX_Destroy(omx);
  return status;
}


int
CACO_G6X_Fill_CYK(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk, SCVAL *ret_sc,
		  char *errbuf, int verbose) 
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;

  /* G6X grammar
  */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {  // order is: L, F, S
	status = dp_recursion_g6x_cyk(foldparam, p, psq, spair, covct, exclude, cyk, G6X_L, j, d, &(cyk->L->dp[j][d]), &nneg, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L caco failed");
	status = dp_recursion_g6x_cyk(foldparam, p, psq, spair, covct, exclude, cyk, G6X_F, j, d, &(cyk->F->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F caco failed");
	status = dp_recursion_g6x_cyk(foldparam, p, psq, spair, covct, exclude, cyk, G6X_S, j, d, &(cyk->S->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S caco failed");
	if (verbose)
	  printf("\nG6X caco S=%f L=%f F=%f | i=%d j=%d d=%d L=%d | ct %d %d\n",
		 cyk->S->dp[j][d], cyk->L->dp[j][d], cyk->F->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]);
      } 
  sc = cyk->S->dp[L][L];
  if (verbose) printf("G6X CYK-score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
CACO_G6X_Inside(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;
  
  e2_FLogsumInit();

  /* G6X grammar
   */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {  // order is: L, F, S
	status = dp_recursion_g6x_inside(foldparam, p, psq, spair, covct, exclude, imx, G6X_L, j, d, &(imx->L->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L caco failed");
	status = dp_recursion_g6x_inside(foldparam, p, psq, spair, covct, exclude, imx, G6X_F, j, d, &(imx->F->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F caco failed");
	status = dp_recursion_g6x_inside(foldparam, p, psq, spair, covct, exclude, imx, G6X_S, j, d, &(imx->S->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S caco failed");
	if (verbose)
	  printf("\nG6X Inside S=%f L=%f F=%f | i=%d j=%d d=%d L=%d | ct %d %d\n",
		 imx->S->dp[j][d], imx->L->dp[j][d], imx->F->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]);
      }
  
  sc = imx->S->dp[L][L];
  if (verbose) printf("G6X Inside score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}


int
CACO_G6X_Outside(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;

  e2_FLogsumInit();
  
  /* G6X grammar
  */
  /* Outside fills j,d in the reverse order than Inside */
  for (j = L; j >= 0; j--) {
    for (d = j; d >= 0; d--)
      { // order is: S, F, L
	status = dp_recursion_g6x_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, G6X_S, j, d, &(omx->S->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S caco failed");
	status = dp_recursion_g6x_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, G6X_F, j, d, &(omx->F->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F caco failed");
	status = dp_recursion_g6x_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, G6X_L, j, d, &(omx->L->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L caco failed");
	if (verbose)
	  printf("\nG6X Outside S=%f L=%f F=%f | i=%d j=%d d=%d L=%d | ct %d %d\n",
		 omx->S->dp[j][d], omx->L->dp[j][d], omx->F->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]);
      }
  }

  sc = omx->S->dp[L][1];
  if (verbose) printf("G6X Outside score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;
  
 ERROR:
  return status;
}


int
CACO_G6X_Posterior(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, POST *post, char *errbuf, int verbose)
{
  SCVAL  sc = -eslINFINITY;
  SCVAL  sum;
  double tol = TOLVAL;
  int    nneg = 0;
  int    L;
  int    j, d;
  int    i;
  int    status;

  L = psq->n;
  
  e2_FLogsumInit();
  
  for (j = 0; j <= L; j++)
    {
      status = dp_recursion_g6x_posterior_single(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_S, j, post, NULL, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S posterior failed");
      status = dp_recursion_g6x_posterior_single(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_F, j, post, NULL, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F posterior failed");
      status = dp_recursion_g6x_posterior_single(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_L, j, post, &nneg, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L posterior failed");

      post->ps[j] -= imx->S->dp[L][L];
      if (post->ps[j] > tol) ESL_XFAIL(eslFAIL, errbuf, "G6X posterior failed ps[%d] = %f\n", j, exp(post->ps[j]));
      if (verbose) printf("\nG6X posterior ps = %f | j=%d | ct %d\n", post->ps[j], j, covct[j]);
      
      for (d = 1; d <= j; d++)
	{
	  i = j - d + 1;
	  status = dp_recursion_g6x_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_S, j, d, post, NULL, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S posterior failed");
	  status = dp_recursion_g6x_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_F, j, d, post, NULL, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F posterior failed");
	  status = dp_recursion_g6x_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_L, j, d, post, &nneg, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L posterior failed");

	  post->pp[i][j] -= imx->S->dp[L][L];
	  if (post->pp[i][j] > tol) ESL_XFAIL(eslFAIL, errbuf, "G6X posterior failed pp[%d][%d] = %f\n", i, j, exp(post->pp[i][j]));

	  post->pp[j][i]  = post->pp[i][j];
 	  if (verbose) printf("G6X posterior pp = %f | i=%d j=%d d=%d | ct %d %d\n", post->pp[i][j], i, j, d, covct[i], covct[j]);
	}
    }
  
  // normalize
  for (i = 1; i <= L; i++) {
    sum = post->ps[i];
    for (j = 1; j <= L; j++) sum = e2_FLogsum(sum, post->pp[i][j]);
    
    post->ps[i] -= sum;
    if (verbose) printf("\nG6X posterior ps = %f | j=%d | ct %d\n", post->ps[i], i, covct[i]);
    
    for (j = 1; j <= L; j++) {
      post->pp[i][j] -= sum;
      if (verbose) printf("G6X posterior pp[%d][%d] = %f | ct %d %d\n", i, j, post->pp[i][j], covct[i], covct[j]);
    }
  }
  
  return eslOK;
  
 ERROR:
  return status;
}


int
CACO_G6XS_Fill_CYK(FOLDPARAM *foldparam, G6XSparam  *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk, SCVAL *ret_sc,
		   char *errbuf, int verbose) 
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;
  
  /* G6XS grammar
  */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      { // order is: L, F, S
	status = dp_recursion_g6xs_cyk(foldparam, p, psq, spair, covct, exclude, cyk, G6X_L, j, d, &(cyk->L->dp[j][d]), &nneg, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L caco failed");
	status = dp_recursion_g6xs_cyk(foldparam, p, psq, spair, covct, exclude, cyk, G6X_F, j, d, &(cyk->F->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F caco failed");
	status = dp_recursion_g6xs_cyk(foldparam, p, psq, spair, covct, exclude, cyk, G6X_S, j, d, &(cyk->S->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S caco failed");
	if (verbose) printf("\nG6XS caco S=%f L=%f F=%f j=%d d=%d L=%d\n", cyk->S->dp[j][d], cyk->L->dp[j][d], cyk->F->dp[j][d], j, d, L); 
     } 
  sc = cyk->S->dp[L][L];
  if (verbose) printf("G6XS CYK-score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
CACO_G6XS_Inside(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;
  
  e2_FLogsumInit();

  /* G6X grammar
   */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {  // order is: L, F, S
	status = dp_recursion_g6xs_inside(foldparam, p, psq, spair, covct, exclude, imx, G6X_L, j, d, &(imx->L->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS L caco failed");
	status = dp_recursion_g6xs_inside(foldparam, p, psq, spair, covct, exclude, imx, G6X_F, j, d, &(imx->F->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS F caco failed");
	status = dp_recursion_g6xs_inside(foldparam, p, psq, spair, covct, exclude, imx, G6X_S, j, d, &(imx->S->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS S caco failed");
	if (verbose)
	  printf("\nG6XS Inside S=%f L=%f F=%f | i=%d j=%d d=%d L=%d | ct %d %d\n",
		 imx->S->dp[j][d], imx->L->dp[j][d], imx->F->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]);
      }
  
  sc = imx->S->dp[L][L];
  if (verbose) printf("G6XS Inside score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
CACO_G6XS_Outside(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;

  e2_FLogsumInit();
  
  /* G6X grammar
  */
  /* Outside fills j,d in the reverse order than Inside */
  for (j = L; j >= 0; j--) {
    for (d = j; d >= 0; d--)
      { // order is: S, F, L
	status = dp_recursion_g6xs_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, G6X_S, j, d, &(omx->S->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS S caco failed");
	status = dp_recursion_g6xs_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, G6X_F, j, d, &(omx->F->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS F caco failed");
	status = dp_recursion_g6xs_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, G6X_L, j, d, &(omx->L->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS L caco failed");
	if (verbose)
	  printf("\nG6XS Outside S=%f L=%f F=%f | i=%d j=%d d=%d L=%d | ct %d %d\n",
		 omx->S->dp[j][d], omx->L->dp[j][d], omx->F->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]);

      }
  }

  sc = omx->S->dp[L][1];
  if (verbose) printf("G6XS Outside score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;
  
 ERROR:
  return status;
}

int
CACO_G6XS_Posterior(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, POST *post, char *errbuf, int verbose)
{
  SCVAL  sc = -eslINFINITY;
  SCVAL  sum;
  double tol = TOLVAL;
  int    nneg = 0;
  int    L;
  int    j, d;
  int    i;
  int    status;

  L = psq->n;
  
  e2_FLogsumInit();
  
  for (j = 0; j <= L; j++)
    {
      status = dp_recursion_g6xs_posterior_single(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_S, j, post, NULL, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS S posterior failed");
      status = dp_recursion_g6xs_posterior_single(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_F, j, post, NULL, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS F posterior failed");
      status = dp_recursion_g6xs_posterior_single(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_L, j, post, &nneg, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS L posterior failed");

      post->ps[j] -= imx->S->dp[L][L];
      if (post->ps[j] > tol) ESL_XFAIL(eslFAIL, errbuf, "G6XS posterior failed ps[%d] = %f\n", j, exp(post->ps[j]));
      if (verbose) printf("\nG6XS posterior ps = %f | j=%d | ct %d\n", post->ps[j], j, covct[j]);
      
      for (d = 1; d <= j; d++)
	{
	  i = j - d + 1;
	  status = dp_recursion_g6xs_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_S, j, d, post, NULL, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S posterior failed");
	  status = dp_recursion_g6xs_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_F, j, d, post, NULL, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F posterior failed");
	  status = dp_recursion_g6xs_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, G6X_L, j, d, post, &nneg, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L posterior failed");

	  post->pp[i][j] -= imx->S->dp[L][L];
	  if (post->pp[i][j] > tol) ESL_XFAIL(eslFAIL, errbuf, "G6XS posterior failed pp[%d][%d] = %f\n", i, j, exp(post->pp[i][j]));

	  post->pp[j][i]  = post->pp[i][j];
 	  if (verbose) printf("G6XS posterior pp = %f | i=%d j=%d d=%d | ct %d %d\n", post->pp[i][j], i, j, d, covct[i], covct[j]);
	}
    }
  
  // normalize
  for (i = 1; i <= L; i++) {
    sum = post->ps[i];
    for (j = 1; j <= L; j++) sum = e2_FLogsum(sum, post->pp[i][j]);
    
    post->ps[i] -= sum;
    if (verbose) printf("\nG6XS posterior ps = %f | j=%d | ct %d\n", post->ps[i], i, covct[i]);
    
    for (j = 1; j <= L; j++) {
      post->pp[i][j] -= sum;
      if (verbose) printf("G6XS posterior pp[%d][%d] = %f | ct %d %d\n", i, j, post->pp[i][j], covct[i], covct[j]);
    }
  }
  
  return eslOK;
  
 ERROR:
  return status;
}

int
CACO_RBG_Fill_CYK(FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *cyk, SCVAL *ret_sc,
		  char *errbuf, int verbose) 
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;

  /* RBG grammar
  * order: (F0 before M1) AND (M1 before R) AND (S last)
  */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {
	status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, RBG_P,  j, d, &(cyk->P->dp[j][d]), &nneg, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG P caco failed");
	status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, RBG_F5, j, d, &(cyk->F5->dp[j][d]), NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG F5 caco failed");
	status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, RBG_F0, j, d, &(cyk->F0->dp[j][d]), NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG F0 caco failed");
	status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, RBG_M1, j, d, &(cyk->M1->dp[j][d]), NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG M1 caco failed");
	status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, RBG_R,  j, d, &(cyk->R->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG R caco failed");
	status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, RBG_M,  j, d, &(cyk->M->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG M caco failed");
	status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, RBG_S,  j, d, &(cyk->S->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG S caco failed");
	if (verbose) 
	  printf("\nRBG CYK P=%f M=%f M1=%f R=%f F5=%f F0=%f S=%f | i=%d j=%d d=%d L=%d | covct %d %d\n", 
		 cyk->P->dp[j][d], cyk->M->dp[j][d], cyk->M1->dp[j][d], cyk->R->dp[j][d],
		 cyk->F5->dp[j][d], cyk->F0->dp[j][d], cyk->S->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]); 
      } 
  sc = cyk->S->dp[L][L];
  if (verbose) printf("RBG CYK-score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
CACO_RBG_Inside(FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;
  
  e2_FLogsumInit();

  /* RBG grammar
   */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      { // (F0 before M1) AND (M1 before R) AND (S last)
	// order is: P F5 F0 M1 R M S
	status = dp_recursion_rbg_inside(foldparam, p, psq, spair, covct, exclude, imx, RBG_P,  j, d, &(imx->P->dp[j][d]),  &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Inside failed at P");
	status = dp_recursion_rbg_inside(foldparam, p, psq, spair, covct, exclude, imx, RBG_F5, j, d, &(imx->F5->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Inside failed at F5");
	status = dp_recursion_rbg_inside(foldparam, p, psq, spair, covct, exclude, imx, RBG_F0, j, d, &(imx->F0->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Inside failed at F0");
	status = dp_recursion_rbg_inside(foldparam, p, psq, spair, covct, exclude, imx, RBG_M1, j, d, &(imx->M1->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Inside failed at M1");
	status = dp_recursion_rbg_inside(foldparam, p, psq, spair, covct, exclude, imx, RBG_R,   j, d, &(imx->R->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Inside failed at R");
	status = dp_recursion_rbg_inside(foldparam, p, psq, spair, covct, exclude, imx, RBG_M,   j, d, &(imx->M->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Inside failed at M");
	status = dp_recursion_rbg_inside(foldparam, p, psq, spair, covct, exclude, imx, RBG_S,   j, d, &(imx->S->dp[j][d]), &nneg, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Inside failed at S");
	
	if (verbose)
	  printf("\nRBG Inside P=%f M=%f M1=%f R=%f F5=%f F0=%f S=%f | i=%d j=%d d=%d L=%d | covct %d %d\n", 
		 imx->P->dp[j][d],  imx->M->dp[j][d],  imx->M1->dp[j][d], imx->R->dp[j][d],
		 imx->F5->dp[j][d], imx->F0->dp[j][d], imx->S->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]); 
      }
  
  sc = imx->S->dp[L][L];
  if (verbose) printf("RBG Inside score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
CACO_RBG_Outside(FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *omx, RBG_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = psq->n;

  e2_FLogsumInit();
  
  /* RBG grammar
  */
  /* Outside fills j,d in the reverse order than Inside */
  for (j = L; j >= 0; j--) {
    for (d = j; d >= 0; d--)
      { // order is: S M R M1 F0 F5 P
	status = dp_recursion_rbg_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, RBG_S,  j, d, &(omx->S->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Outsice faile at S");
	status = dp_recursion_rbg_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, RBG_M,  j, d, &(omx->M->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Outsice faile at M");
	status = dp_recursion_rbg_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, RBG_R,  j, d, &(omx->R->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Outsice faile at R");
	status = dp_recursion_rbg_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, RBG_M1, j, d, &(omx->M1->dp[j][d]), NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Outsice faile at M1");
	status = dp_recursion_rbg_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, RBG_F0, j, d, &(omx->F0->dp[j][d]), NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Outsice faile at F0");
	status = dp_recursion_rbg_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, RBG_F5, j, d, &(omx->F5->dp[j][d]), NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Outsice faile at F5");
	status = dp_recursion_rbg_outside(foldparam, p, psq, spair, covct, exclude, omx, imx, RBG_P,  j, d, &(omx->P->dp[j][d]),  NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG Outsice faile at S");
	if (verbose)
	  printf("\nRBG Outside S=%f F0 %f F5 %f P %f M %f M1 %f R %f| i=%d j=%d d=%d L=%d | ct %d %d\n",
		 omx->S->dp[j][d], omx->F0->dp[j][d], omx->F5->dp[j][d], omx->P->dp[j][d], omx->M->dp[j][d], omx->M1->dp[j][d], omx->R->dp[j][d],
		 j-d+1, j, d, L, covct[j-d+1], covct[j]);
      }
  }

  sc = omx->S->dp[L][1];
  if (verbose) printf("RBG Outside score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;
  
 ERROR:
  return status;
}

int
CACO_RBG_Posterior(FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx, RBG_MX *omx, POST *post, char *errbuf, int verbose)
{
  SCVAL  sc = -eslINFINITY;
  SCVAL  sum;
  double tol = TOLVAL;
  int    nneg = 0;
  int    L;
  int    j, d;
  int    i;
  int    status;

  L = psq->n;
  
  e2_FLogsumInit();
  
  for (j = 0; j <= L; j++)
    {
     for (d = 1; d <= j; d++)
	{
	  i = j - d + 1;
	  status = dp_recursion_rbg_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, RBG_F0, j, d, post, NULL, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG F0 pp posterior failed");
	  status = dp_recursion_rbg_posterior_pair(foldparam, p, psq, spair, covct, exclude, imx, omx, RBG_F5, j, d, post, NULL, errbuf, verbose);
	  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "RBG F5 pp posterior failed");

	  post->pp[i][j] -= imx->S->dp[L][L];
	  if (verbose) printf("RBG posterior pp = %f | i=%d j=%d d=%d | ct %d %d | %f\n", post->pp[i][j], i, j, d, covct[i], covct[j], imx->S->dp[L][L]);
		  
	  if (post->pp[i][j] > tol) ESL_XFAIL(eslFAIL, errbuf, "RBG posterior failed pp[%d][%d] = %f | ct %d %d\n", i, j, exp(post->pp[i][j]), covct[i], covct[j]);

	  post->pp[j][i]  = post->pp[i][j];
 	}
    }
  
  // obtain the ps probabilities from the pp
  for (i = 1; i <= L; i++) {
    sum = -eslINFINITY;
    for (j = 1; j <= L; j++) sum = e2_FLogsum(sum, post->pp[i][j]);
    
    post->ps[i] = log(1-exp(sum));
    if (verbose) printf("\nRBG posterior ps = %f | j=%d | ct %d\n", post->ps[i], i, covct[i]);
  }
  
  return eslOK;
  
 ERROR:
  return status;
}

int
CACO_G6X_Traceback_CYK(ESL_RANDOMNESS *rng, FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude,
		       G6X_MX *cyk, int *ct, char *errbuf, int verbose) 
{
  ESL_STACK      *ns = NULL;             /* integer pushdown stack for traceback */
  ESL_STACK      *alts = NULL;           /* stack of alternate equal-scoring tracebacks */
  SCVAL           bestsc;                /* max score over possible rules */
  SCVAL           fillsc;                /* max score in fill */
  int             L ;
  int             nequiv;                /* number of equivalent alternatives for a traceback */
  int             x;                     /* a random choice from nequiv */
  int             w;                     /* index of a non terminal S (w=0) L (w=1) F (w=2) */
  int             r;                     /* index of a rule */
  int             d, d1;                 /* optimum values of d1 iterator */
  int             i,j,k;                 /* seq coords */
  float           tol = TOLVAL;
  int             status;

  L = psq->n;

  /* is sq score is -infty, nothing to traceback */
  if (cyk->S->dp[L][L] == -eslINFINITY) {
    if (1||verbose) printf("G6X no traceback.\n");
    return eslOK;  
  }

  /* We're going to do a simple traceback that only
   * remembers who was a base pair, and keeps a ct[]
   * array. 
   */
  esl_vec_ISet(ct, L+1, 0);
  
  /* We implement a "stochastic" traceback, which chooses randomly
   * amongst equal-scoring alternative parse trees. This is particularly
   * essential for working with ambiguous grammars, for which 
   * choosing an arbitrary optimal parse tree by order of evaluation
   * can easily result in infinite loops. To do this, we will keep
   * a stack of the alternate solutions.
   */
  alts = esl_stack_ICreate();

  /* Start an integer stack for traversing the traceback.
   * push w,i,j = G->ntS_idx,1,L to init. 
   */
  w = G6X_S;
  ns = esl_stack_ICreate();
  esl_stack_IPush(ns, w);
  esl_stack_IPush(ns, 1);
  esl_stack_IPush(ns, L);
  
  while (esl_stack_ObjectCount(ns) != 0)
    {
      esl_stack_IPop(ns, &j);
      esl_stack_IPop(ns, &i);
      esl_stack_IPop(ns, &w);
      d = j-i+1;
      
      status = dp_recursion_g6x_cyk(foldparam, p, psq, spair, covct, exclude, cyk, w, j, d, &bestsc, NULL, alts, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X cyk failed");
      
      /* Some assertions.
       */
      switch(w) {
      case G6X_S: fillsc = cyk->S->dp[j][d]; break;
      case G6X_L: fillsc = cyk->L->dp[j][d]; break;
      case G6X_F: fillsc = cyk->F->dp[j][d]; break;
      }
      if (fabs(bestsc - fillsc) > tol) 
	ESL_XFAIL(eslFAIL, errbuf, "CACO_G6X_Traceback(): that can't happen either. i=%d j=%d d=%d bestsc %f cyk %f", 
		  j-d+1, j, d, bestsc, fillsc); 
      
      /* Now we know one or more equiv solutions, and they're in
       * the stack <alts>, which keeps 2 numbers (r, d1) for each
       * solution. Choose one of them at random.
       */
      nequiv = esl_stack_ObjectCount(alts) / 2; /* how many solutions? */
      x = esl_rnd_Roll(rng, nequiv);            /* uniformly, 0.nequiv-1 */
      esl_stack_DiscardTopN(alts, x*2);         /* dig down to choice */
      esl_stack_IPop(alts, &d1);
      esl_stack_IPop(alts, &r);
      
      /* Now we know a best rule; figure out where we came from,
       * and push that info onto the <ns> stack.
       */
      if (verbose) {
        printf("-----------------------------------\n"); 
        printf("w=%d i=%d j=%d d=%d d1=%d\n", w, j-d+1, j, d, d1);
	printf("tracing %f\n", bestsc);
        printf("       rule(%d)\n", r);
      }
 
      if (w == G6X_S && r != G6X_S_1 && r != G6X_S_2 && r != G6X_S_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with S", r);
      if (w == G6X_L && r != G6X_L_1 && r != G6X_L_2 && r != G6X_L_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with L", r);
      if (w == G6X_F && r != G6X_F_1 && r != G6X_F_2 && r != G6X_F_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with F", r);
      
      i = j - d  + 1;
      k = i + d1 - 1;
      
      switch(r) {
      case G6X_S_1: // S -> LS
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	
	esl_stack_IPush(ns, G6X_S);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      case G6X_S_2: // S -> L
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, j);
	break;
      case G6X_S_3: // S -> epsilon
	break;
      case G6X_L_1: // L -> a F a'
	esl_stack_IPush(ns, G6X_F);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_L_2: // L -> a a'
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_L_3: // L -> a
	break;
      case G6X_F_1: // F -> a F a'
	esl_stack_IPush(ns, G6X_F);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_F_2: // F -> a a'
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_F_3: // F -> LS
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	
	esl_stack_IPush(ns, G6X_S);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      default: ESL_XFAIL(eslFAIL, errbuf, "rule %d disallowed. Max number is %d", r, G6X_NR);
      }
    }

  esl_stack_Destroy(ns);
  esl_stack_Destroy(alts);
  return eslOK;
  
 ERROR:
  if (ns)   esl_stack_Destroy(ns); ns = NULL;
  if (alts) esl_stack_Destroy(alts); alts = NULL;
  return status;
}



int
CACO_G6XS_Traceback_CYK(ESL_RANDOMNESS *rng, FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude,
			G6X_MX *cyk, int *ct, char *errbuf, int verbose) 
{
  ESL_STACK      *ns = NULL;             /* integer pushdown stack for traceback */
  ESL_STACK      *alts = NULL;           /* stack of alternate equal-scoring tracebacks */
  SCVAL           bestsc;                /* max score over possible rules */
  SCVAL           fillsc;                /* max score in fill */
  int             L;
  int             nequiv;                /* number of equivalent alternatives for a traceback */
  int             x;                     /* a random choice from nequiv */
  int             w;                     /* index of a non terminal S (w=0) L (w=1) F (w=2) */
  int             r;                     /* index of a rule */
  int             d, d1;                 /* optimum values of d1 iterator */
  int             i,j,k;                 /* seq coords */
  float           tol = TOLVAL;
  int             status;

  L = psq->n;

  /* is sq score is -infty, nothing to traceback */
  if (cyk->S->dp[L][L] == -eslINFINITY) {
    if (1||verbose) printf("G6XS no traceback.\n");
    return eslOK;
  }

  /* We're going to do a simple traceback that only
   * remembers who was a base pair, and keeps a ct[]
   * array. 
   */
  esl_vec_ISet(ct, L+1, 0);
  
  /* We implement a "stochastic" traceback, which chooses randomly
   * amongst equal-scoring alternative parse trees. This is particularly
   * essential for working with ambiguous grammars, for which 
   * choosing an arbitrary optimal parse tree by order of evaluation
   * can easily result in infinite loops. To do this, we will keep
   * a stack of the alternate solutions.
   */
  alts = esl_stack_ICreate();

  /* Start an integer stack for traversing the traceback.
   * push w,i,j = G->ntS_idx,1,L to init. 
   */
  w = G6X_S;
  ns = esl_stack_ICreate();
  esl_stack_IPush(ns, w);
  esl_stack_IPush(ns, 1);
  esl_stack_IPush(ns, L);
  
  while (esl_stack_ObjectCount(ns) != 0)
    {
      esl_stack_IPop(ns, &j);
      esl_stack_IPop(ns, &i);
      esl_stack_IPop(ns, &w);
      d = j-i+1;
     
      status = dp_recursion_g6xs_cyk(foldparam, p, psq, spair, covct, exclude, cyk, w, j, d, &bestsc, NULL, alts, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6XS cyk failed");
      
      /* Some assertions.
       */
      switch(w) {
      case G6X_S: fillsc = cyk->S->dp[j][d]; break;
      case G6X_L: fillsc = cyk->L->dp[j][d]; break;
      case G6X_F: fillsc = cyk->F->dp[j][d]; break;
      }
      if (fabs(bestsc - fillsc) > tol) 
	ESL_XFAIL(eslFAIL, errbuf, "CACO_G6X_Traceback(): that can't happen either. i=%d j=%d d=%d bestsc %f cyk %f", 
		  j-d+1, j, d, bestsc, fillsc); 
      
      /* Now we know one or more equiv solutions, and they're in
       * the stack <alts>, which keeps 2 numbers (r, d1) for each
       * solution. Choose one of them at random.
       */
      nequiv = esl_stack_ObjectCount(alts) / 2; /* how many solutions? */
      x = esl_rnd_Roll(rng, nequiv);            /* uniformly, 0.nequiv-1 */
      esl_stack_DiscardTopN(alts, x*2);         /* dig down to choice */
      esl_stack_IPop(alts, &d1);
      esl_stack_IPop(alts, &r);
      
      if (w == G6X_S && r != G6X_S_1 && r != G6X_S_2 && r != G6X_S_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with S", r);
      if (w == G6X_L && r != G6X_L_1 && r != G6X_L_2 && r != G6X_L_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with L", r);
      if (w == G6X_F && r != G6X_F_1 && r != G6X_F_2 && r != G6X_F_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with F", r);
      
      /* Now we know a best rule; figure out where we came from,
       * and push that info onto the <ns> stack.
       */
      if (verbose) {
        printf("-----------------------------------\n"); 
        printf("i=%d j=%d d=%d d1=%d\n", j-d+1, j, d, d1);
	printf("tracing %f\n", bestsc);
        printf("   rule(%d)\n", r);
      }
      
      i = j - d  + 1;
      k = i + d1 - 1;
      
      switch(r) {
      case G6X_S_1: // S -> LS
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	
	esl_stack_IPush(ns, G6X_S);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      case G6X_S_2: // S -> L
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, j);
	break;
     case G6X_S_3: // S -> epsilon
	break;
      case G6X_L_1: // L -> a F a'
	esl_stack_IPush(ns, G6X_F);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_L_2: // L -> a a'
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_L_3: // L -> a
	break;
      case G6X_F_1: // F -> a F a'
	esl_stack_IPush(ns, G6X_F);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_F_2: // F -> a a'
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_F_3: // F -> LS
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	
	esl_stack_IPush(ns, G6X_S);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      default: ESL_XFAIL(eslFAIL, errbuf, "rule %d disallowed. Max number is %d", r, G6X_NR);
      }
    }

  esl_stack_Destroy(ns);
  esl_stack_Destroy(alts);
  return eslOK;
  
 ERROR:
  if (ns)   esl_stack_Destroy(ns); ns = NULL;
  if (alts) esl_stack_Destroy(alts); alts = NULL;
  return status;
}

int
CACO_RBG_Traceback_CYK(ESL_RANDOMNESS *rng, FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude,
		       RBG_MX *cyk, int *ct, char *errbuf, int verbose) 
{
  ESL_STACK      *ns   = NULL;           /* integer pushdown stack for traceback */
  ESL_STACK      *alts = NULL;           /* stack of alternate equal-scoring tracebacks */
  SCVAL           bestsc;                /* max score over possible rules */
  SCVAL           fillsc;                /* max score in fill */
  int             L;
  int             nequiv;                /* number of equivalent alternatives for a traceback */
  int             x;                     /* a random choice from nequiv */
  int             w;                     /* index of a non terminal S (w=0) L (w=1) F (w=2) */
  int             r;                     /* index of a rule */
  int             d, d1, d2;             /* optimum values of d1 iterator */
  int             i,j,k,l;               /* seq coords */
  float           tol = TOLVAL;
  int             status;

  L = psq->n;

  /* is sq score is -infty, nothing to traceback */
  if (cyk->S->dp[L][L] == -eslINFINITY) {
    if (1||verbose) printf("RBG no traceback.\n");
    return eslOK;
  }

   /* We're going to do a simple traceback that only
   * remembers who was a base pair, and keeps a ct[]
   * array. 
   */
  esl_vec_ISet(ct, L+1, 0);
  
  /* We implement a "stochastic" traceback, which chooses randomly
   * amongst equal-scoring alternative parse trees. This is particularly
   * essential for working with ambiguous grammars, for which 
   * choosing an arbitrary optimal parse tree by order of evaluation
   * can easily result in infinite loops. To do this, we will keep
   * a stack of the alternate solutions.
   */
  alts = esl_stack_ICreate();

  /* Start an integer stack for traversing the traceback.
   * push w,i,j = G->ntS_idx,1,L to init. 
   */
  w = RBG_S;
  ns = esl_stack_ICreate();
  esl_stack_IPush(ns, w);
  esl_stack_IPush(ns, 1);
  esl_stack_IPush(ns, L);
  
  while (esl_stack_ObjectCount(ns) != 0)
    {
      esl_stack_IPop(ns, &j);
      esl_stack_IPop(ns, &i);
      esl_stack_IPop(ns, &w);
      d = j-i+1;

      status = dp_recursion_rbg_cyk(foldparam, p, psq, spair, covct, exclude, cyk, w, j, d, &bestsc, NULL, alts, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "CYK failed");
      
      /* Some assertions.
       */
      switch(w) {
      case RBG_S:  fillsc = cyk->S->dp[j][d];  break;
      case RBG_F0: fillsc = cyk->F0->dp[j][d]; break;
      case RBG_F5: fillsc = cyk->F5->dp[j][d]; break;
      case RBG_P:  fillsc = cyk->P->dp[j][d];  break;
      case RBG_M:  fillsc = cyk->M->dp[j][d];  break;
      case RBG_R:  fillsc = cyk->R->dp[j][d];  break;
      case RBG_M1: fillsc = cyk->M1->dp[j][d]; break;
      }
 
      if (fabs(bestsc - fillsc) > tol) 
	ESL_XFAIL(eslFAIL, errbuf, "CACO_RBG_Traceback(): that can't happen either. i=%d j=%d d=%d bestsc %f cyk %f", 
		  j-d+1, j, d, bestsc, fillsc); 
      
      /* Now we know one or more equiv solutions, and they're in
       * the stack <alts>, which keeps 3 numbers (r, d1, d2) for each
       * solution. Choose one of them at random.
       */
      nequiv = esl_stack_ObjectCount(alts) / 3; /* how many solutions? */
      x = esl_rnd_Roll(rng, nequiv);            /* uniformly, 0.nequiv-1 */
      esl_stack_DiscardTopN(alts, x*3);         /* dig down to choice */
      esl_stack_IPop(alts, &d2);
      esl_stack_IPop(alts, &d1);
      esl_stack_IPop(alts, &r);
      
      if (w == RBG_S  && r != RBG_S_1  && r != RBG_S_2  && r != RBG_S_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with S",  r);
      if (w == RBG_F0 && r != RBG_F0_1 && r != RBG_F0_2 && r != RBG_F0_3) ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with F0", r);
      if (w == RBG_F5 && r != RBG_F5_1 && r != RBG_F5_2 && r != RBG_F5_3) ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with F5", r);
      if (w == RBG_M  && r != RBG_M_1  && r != RBG_M_2)                   ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with M",  r);
      if (w == RBG_R  && r != RBG_R_1  && r != RBG_R_2)                   ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with R",  r);
      if (w == RBG_M1 && r != RBG_M1_1 && r != RBG_M1_2)                  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with M1", r);
      if (w == RBG_P  && r != RBG_P_1  && r != RBG_P_2 && r != RBG_P_3 && r != RBG_P_4 && r != RBG_P_5)               
	ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with P", r);
     
      /* Now we know a best rule; figure out where we came from,
       * and push that info onto the <ns> stack.
       */
      if (verbose) {
        printf("-----------------------------------\n"); 
        printf("i=%d j=%d d=%d d1=%d d2=%d\n", j-d+1, j, d, d1, d2);
	printf("tracing %f\n", bestsc);
        printf("   w=%d rule(%d)\n", w, r);
      }
      
      i = j - d  + 1;
      k = i + d1 - 1;
      l = j - d2 + 1;
      
      switch(r) {
      case RBG_S_1: // S -> a S
	esl_stack_IPush(ns, RBG_S);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j);
	break;
      case RBG_S_2: // S -> F0 S
	esl_stack_IPush(ns, RBG_F0);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	esl_stack_IPush(ns, RBG_S);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      case RBG_S_3: // S -> epsilon
	break;

      case RBG_F0_1: // F0 -> a F5 a'
	esl_stack_IPush(ns, RBG_F5);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case RBG_F0_2: // F0 -> a P a'
	esl_stack_IPush(ns, RBG_P);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case RBG_F0_3: // F0 -> a a'
	ct[i] = j;
	ct[j] = i;
	break;
 
      case RBG_F5_1: // F5 -> a F5 a'
	esl_stack_IPush(ns, RBG_F5);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case RBG_F5_2: // F5 -> a P a'
	esl_stack_IPush(ns, RBG_P);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case RBG_F5_3: // F5 -> a a'
	ct[i] = j;
	ct[j] = i;
	break;

      case RBG_P_1: // P -> m..m
	break;
      case RBG_P_2: // P -> m..m F0
	esl_stack_IPush(ns, RBG_F0);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      case RBG_P_3: // P -> F0 m..m
	esl_stack_IPush(ns, RBG_F0);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, l-1);
	break;
      case RBG_P_4: // P -> m..m F0 m..m
	esl_stack_IPush(ns, RBG_F0);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, l-1);
	break;
      case RBG_P_5: // P -> M1 M
	esl_stack_IPush(ns, RBG_M1);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	esl_stack_IPush(ns, RBG_M);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
	
      case RBG_M_1: // M -> M1 M
	esl_stack_IPush(ns, RBG_M1);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	esl_stack_IPush(ns, RBG_M);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      case RBG_M_2: // M -> R
	esl_stack_IPush(ns, RBG_R);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, j);
	break;

      case RBG_R_1: // R -> R a
	esl_stack_IPush(ns, RBG_R);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, j-1);
	break;
      case RBG_R_2: // R -> M1
	esl_stack_IPush(ns, RBG_M1);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, j);
	break;
 
     case RBG_M1_1: // M1 -> a M1
	esl_stack_IPush(ns, RBG_M1);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j);
	break;
      case RBG_M1_2: // M1 -> F0
	esl_stack_IPush(ns, RBG_F0);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, j);
	break;

     default: 
       printf("rule %d disallowed. Max number is %d", r, RBG_NR);
       ESL_XFAIL(eslFAIL, errbuf, "rule %d disallowed. Max number is %d", r, RBG_NR);
       break;
      }
    }

  esl_stack_Destroy(ns);
  esl_stack_Destroy(alts);
  return eslOK;
  
 ERROR:
  if (ns)   esl_stack_Destroy(ns); ns = NULL;
  if (alts) esl_stack_Destroy(alts); alts = NULL;
  return status;
}

int
CACO_G6X_MEA_GetParam(G6Xparam **ret_p, double gamma, char *errbuf, int verbose)
{
  G6Xparam *p = NULL;
  double   lg = (gamma > 1)? log(gamma) : 0;
  int       x;
  int       status;

 ESL_ALLOC(p, sizeof(G6Xparam));

  p->t1[0] = 0;
  p->t1[1] = 0;
  p->t1[2] = 0;
  p->t2[0] = lg;
  p->t2[1] = lg;
  p->t2[2] = 0;
  p->t3[0] = lg;
  p->t3[1] = lg;
  p->t3[2] = 0;

  *ret_p = p;
  return eslOK;

 ERROR:
  if (p) free(p);
  return status;
}

int
CACO_MEA_Fill_CYK(FOLDPARAM *foldparam, G6Xparam *meap, POST *post, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *gmx, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   nneg = 0;
  int   L;
  int   j, d;
  int   status;

  L = post->L;

  /* G6X grammar
  */
  for (j = 0; j <= L; j++)
    for (d = 1; d <= j; d++)
      {  // order is: L, F, S
	status = dp_recursion_mea_cyk(foldparam, meap, post, spair, covct, exclude, gmx, G6X_L, j, d, &(gmx->L->dp[j][d]), &nneg, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X L caco failed");
	status = dp_recursion_mea_cyk(foldparam, meap, post, spair, covct, exclude, gmx, G6X_F, j, d, &(gmx->F->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X F caco failed");
	status = dp_recursion_mea_cyk(foldparam, meap, post, spair, covct, exclude, gmx, G6X_S, j, d, &(gmx->S->dp[j][d]),  NULL, NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6X S caco failed");
	if (verbose)
	  printf("\nG6X MEA S=%f L=%f F=%f | i=%d j=%d d=%d L=%d | ct %d %d\n",
		 gmx->S->dp[j][d], gmx->L->dp[j][d], gmx->F->dp[j][d], j-d+1, j, d, L, covct[j-d+1], covct[j]);
      } 
  sc = gmx->S->dp[L][L];
  if (verbose) printf("MEA cyk score = %f\n# negatives = %d\n", sc, nneg);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
CACO_MEA_Traceback_CYK(ESL_RANDOMNESS *rng, FOLDPARAM *foldparam, G6Xparam *meap, POST *post, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *gmx, int *ct,
		       char *errbuf, int verbose)
{
  ESL_STACK      *ns = NULL;             /* integer pushdown stack for traceback */
  ESL_STACK      *alts = NULL;           /* stack of alternate equal-scoring tracebacks */
  SCVAL           bestsc;                /* max score over possible rules */
  SCVAL           fillsc;                /* max score in fill */
  int             L ;
  int             nequiv;                /* number of equivalent alternatives for a traceback */
  int             x;                     /* a random choice from nequiv */
  int             w;                     /* index of a non terminal S (w=0) L (w=1) F (w=2) */
  int             r;                     /* index of a rule */
  int             d, d1;                 /* optimum values of d1 iterator */
  int             i,j,k;                 /* seq coords */
  float           tol = TOLVAL;
  int             status;

  L = post->L;

  /* is sq score is -infty, nothing to traceback */
  if (gmx->S->dp[L][L] == -eslINFINITY) {
    if (1||verbose) printf("MEA no traceback.\n");
    return eslOK;  
  }

  /* We're going to do a simple traceback that only
   * remembers who was a base pair, and keeps a ct[]
   * array. 
   */
  esl_vec_ISet(ct, L+1, 0);
  
  /* We implement a "stochastic" traceback, which chooses randomly
   * amongst equal-scoring alternative parse trees. This is particularly
   * essential for working with ambiguous grammars, for which 
   * choosing an arbitrary optimal parse tree by order of evaluation
   * can easily result in infinite loops. To do this, we will keep
   * a stack of the alternate solutions.
   */
  alts = esl_stack_ICreate();

  /* Start an integer stack for traversing the traceback.
   * push w,i,j = G->ntS_idx,1,L to init. 
   */
  w = G6X_S;
  ns = esl_stack_ICreate();
  esl_stack_IPush(ns, w);
  esl_stack_IPush(ns, 1);
  esl_stack_IPush(ns, L);
  
  while (esl_stack_ObjectCount(ns) != 0)
    {
      esl_stack_IPop(ns, &j);
      esl_stack_IPop(ns, &i);
      esl_stack_IPop(ns, &w);
      d = j-i+1;
      
      status = dp_recursion_mea_cyk(foldparam, meap, post, spair, covct, exclude, gmx, w, j, d, &bestsc, NULL, alts, errbuf, verbose);
      if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "MEA cyk failed");
      
      /* Some assertions.
       */
      switch(w) {
      case G6X_S: fillsc = gmx->S->dp[j][d]; break;
      case G6X_L: fillsc = gmx->L->dp[j][d]; break;
      case G6X_F: fillsc = gmx->F->dp[j][d]; break;
      }
      if (fabs(bestsc - fillsc) > tol) 
	ESL_XFAIL(eslFAIL, errbuf, "CACO_MEA_Traceback(): that can't happen either. i=%d j=%d d=%d bestsc %f gmx %f", 
		  j-d+1, j, d, bestsc, fillsc); 
      
      /* Now we know one or more equiv solutions, and they're in
       * the stack <alts>, which keeps 2 numbers (r, d1) for each
       * solution. Choose one of them at random.
       */
      nequiv = esl_stack_ObjectCount(alts) / 2; /* how many solutions? */
      x = esl_rnd_Roll(rng, nequiv);            /* uniformly, 0.nequiv-1 */
      esl_stack_DiscardTopN(alts, x*2);         /* dig down to choice */
      esl_stack_IPop(alts, &d1);
      esl_stack_IPop(alts, &r);
      
      /* Now we know a best rule; figure out where we came from,
       * and push that info onto the <ns> stack.
       */
      if (verbose) {
        printf("-----------------------------------\n"); 
        printf("w=%d i=%d j=%d d=%d d1=%d\n", w, j-d+1, j, d, d1);
	printf("tracing %f\n", bestsc);
        printf("       rule(%d)\n", r);
      }
 
      if (w == G6X_S && r != G6X_S_1 && r != G6X_S_2 && r != G6X_S_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with S", r);
      if (w == G6X_L && r != G6X_L_1 && r != G6X_L_2 && r != G6X_L_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with L", r);
      if (w == G6X_F && r != G6X_F_1 && r != G6X_F_2 && r != G6X_F_3)  ESL_XFAIL(eslFAIL, errbuf, "rule %d cannot appear with F", r);
      
      i = j - d  + 1;
      k = i + d1 - 1;
      
      switch(r) {
      case G6X_S_1: // S -> LS
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	
	esl_stack_IPush(ns, G6X_S);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      case G6X_S_2: // S -> L
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, j);
	break;
      case G6X_S_3: // S -> epsilon
	break;
      case G6X_L_1: // L -> a F a'
	esl_stack_IPush(ns, G6X_F);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_L_2: // L -> a a'
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_L_3: // L -> a
	break;
      case G6X_F_1: // F -> a F a'
	esl_stack_IPush(ns, G6X_F);
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j-1);
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_F_2: // F -> a a'
	ct[i] = j;
	ct[j] = i;
	break;
      case G6X_F_3: // F -> LS
	esl_stack_IPush(ns, G6X_L);
	esl_stack_IPush(ns, i);
	esl_stack_IPush(ns, k);
	
	esl_stack_IPush(ns, G6X_S);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	break;
      default: ESL_XFAIL(eslFAIL, errbuf, "rule %d disallowed. Max number is %d", r, G6X_NR);
      }
    }

  esl_stack_Destroy(ns);
  esl_stack_Destroy(alts);
  return eslOK;
  
 ERROR:
  if (ns)   esl_stack_Destroy(ns); ns = NULL;
  if (alts) esl_stack_Destroy(alts); alts = NULL;
  return status;
}



/*----------------------------------- internal functions --------------------------------------------------*/

static int
dp_recursion_mea_cyk(FOLDPARAM *foldparam, G6Xparam *p, POST *post, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk,
		     int w, int j, int d, SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose)
{
  SCVAL    bestsc = -eslINFINITY;
  SCVAL    sc;
  double   emitsc_singi;
  double   emitsc_pairij;
  int      allow_bp, force_bp;
  int      allow_si;
  int      L = post->L;
  int      d1;
  int      i, k;
  int      status;

  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;

  if (d < 1 && w == G6X_L)  { *ret_sc = -eslINFINITY; return eslOK; }  // L  has at least 1 residues
  if (d < 2 && w == G6X_F)  { *ret_sc = -eslINFINITY; return eslOK; }  // F  has at least 2 residues

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  allow_si = allow_single(i, covct);
  
  // emission scores
  emitsc_singi  = post->ps[i];
  emitsc_pairij = post->pp[i][j];
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    /* rule0: S -> LS */
    if (d > 0) {
      for (d1 = 0; d1 <= d; d1++) {
	k = i + d1 - 1;
	
	sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1] + p->t1[0];      
	if (sc >= bestsc) {
	  if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	    if (alts) esl_stack_Reuse(alts);
	    bestsc = sc;
	  }     
	  if (alts) {
	    esl_stack_IPush(alts, G6X_S_1);
	    esl_stack_IPush(alts, d1);
	  }
	}
      }
    }
    
    /* rule1: S -> L */
    d1 = 0;
    sc = cyk->L->dp[j][d] + p->t1[1];
    if (sc >= bestsc) {
      if (sc > bestsc) {   /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_S_2);
	esl_stack_IPush(alts, d1);
      }
    }
    
    /* rule2: S -> epsilon */
    d1 = 0;
    if (d == 0) {
      sc = p->t1[2];
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_S_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;     

  case G6X_L:
    /* rule3: L -> a F a' */
    
    if (force_bp) {
      sc = cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij;
    }
    else 
      sc = (allow_bp)?
	cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij : -eslINFINITY;
    
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_L_1);
	esl_stack_IPush(alts, d1);
      }
    }
     /* rule4: L -> a a' */
    d1 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->t2[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t2[1] + emitsc_pairij : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_L_2);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    
    /* rule5: L -> a */
    d1 = 0;
    if (d == 1) {
      sc = (allow_si)? p->t2[2] + emitsc_singi : -eslINFINITY;
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_L_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;

  case G6X_F:
    /* rule6: F -> a F a' */
    d1 = 0;
    if (force_bp) 
      sc = cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij;
    else 
      sc = (allow_bp)?
	cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij : -eslINFINITY;
      
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_F_1);
	esl_stack_IPush(alts, d1);
      }
    }
    
     /* rule7: F -> a a' */
    d1 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->t3[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t3[1] + emitsc_pairij : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_F_2);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    
    /* rule8: F -> LS */
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;

      sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1] + p->t3[2];

      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_F_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  *ret_sc = bestsc;
  return eslOK;

 ERROR:
  return status;
}


static int 
dp_recursion_g6x_cyk(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk, int w, int j, int d,
		     SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose)
{
  SCVAL    bestsc = -eslINFINITY;
  SCVAL    sc;
  double   emitsc_singi;
  double   emitsc_pairij;
  int      allow_bp, force_bp;
  int      allow_si;
  int      L = psq->n;
  int      d1;
  int      i, k;
  int      status;

  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;

  if (d < 1 && w == G6X_L)  { *ret_sc = -eslINFINITY; return eslOK; }  // L  has at least 1 residues
  if (d < 2 && w == G6X_F)  { *ret_sc = -eslINFINITY; return eslOK; }  // F  has at least 2 residues

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  allow_si = allow_single(i, covct);
  
  // emission scores
  emitsc_singi  = emitsc_sing_prof(i,    psq, p->e_sing);
  emitsc_pairij = emitsc_pair_prof(i, j, psq, p->e_pair);
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    /* rule0: S -> LS */
    if (d > 0) {
      for (d1 = 0; d1 <= d; d1++) {
	k = i + d1 - 1;
	
	sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1] + p->t1[0];      
	if (sc >= bestsc) {
	  if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	    if (alts) esl_stack_Reuse(alts);
	    bestsc = sc;
	  }     
	  if (alts) {
	    esl_stack_IPush(alts, G6X_S_1);
	    esl_stack_IPush(alts, d1);
	  }
	}
      }
    }
    
    /* rule1: S -> L */
    d1 = 0;
    sc = cyk->L->dp[j][d] + p->t1[1];
    if (sc >= bestsc) {
      if (sc > bestsc) {   /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_S_2);
	esl_stack_IPush(alts, d1);
      }
    }
    
    /* rule2: S -> epsilon */
    d1 = 0;
    if (d == 0) {
      sc = p->t1[2];
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_S_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;     

  case G6X_L:
    /* rule3: L -> a F a' */
    d1 = 0;
    if (force_bp) {
      sc = cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij;
    }
    else 
      sc = (allow_bp)?
	cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij : -eslINFINITY;
    
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_L_1);
	esl_stack_IPush(alts, d1);
      }
    }
     /* rule4: L -> a a' */
    d1 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->t2[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t2[1] + emitsc_pairij : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_L_2);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    
    /* rule5: L -> a */
    d1 = 0;
    if (d == 1) {
      sc = (allow_si)? p->t2[2] + emitsc_singi : -eslINFINITY;
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_L_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;

  case G6X_F:
    /* rule6: F -> a F a' */
    d1 = 0;
    if (force_bp) 
      sc = cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij;
    else 
      sc = (allow_bp)?
	cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij : -eslINFINITY;
      
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_F_1);
	esl_stack_IPush(alts, d1);
      }
    }
    
     /* rule7: F -> a a' */
    d1 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->t3[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t3[1] + emitsc_pairij : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_F_2);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    
    /* rule8: F -> LS */
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;

      sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1] + p->t3[2];

      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_F_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  *ret_sc = bestsc;
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6x_inside(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, int w, int j, int d,
			SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    sumsc = -eslINFINITY;	/* sum score so far */
  SCVAL    sc;	        	        /* score for a rule */
  double   emitsc_singi;
  double   emitsc_pairij;
  int      allow_bp, force_bp;
  int      allow_si;
  int      L = psq->n;
  int      d1;
  int      i, k;
  int      status;

  i = j - d + 1;

  if (d < 1 && w == G6X_L)  { *ret_sc = -eslINFINITY; return eslOK; }  // L  has at least 1 residues
  if (d < 1 && w == G6X_F)  { *ret_sc = -eslINFINITY; return eslOK; }  // F  has at least 1 residues

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  allow_si = allow_single(i, covct);
  
  // emission scores
  emitsc_singi  = emitsc_sing_prof(i,    psq, p->e_sing);
  emitsc_pairij = emitsc_pair_prof(i, j, psq, p->e_pair);
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    /* rule0: S -> LS */
    //     L        S
    //  i_____k k+1____j
    //  i______________j
    //          S
   
    if (d > 0) {
      for (d1 = 0; d1 <= d; d1++) {
	k = i + d1 - 1;
	
	sc    = imx->L->dp[k][d1] + imx->S->dp[j][d-d1] + p->t1[0];
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
    
    /* rule1: S -> L */
    sc    = imx->L->dp[j][d] + p->t1[1];
    sumsc = e2_FLogsum(sumsc, sc);
    
    /* rule2: S -> epsilon */
    if (d == 0) {
      sc    = p->t1[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;     

  case G6X_L:
    /* rule3: L -> a F a' */
    if (d > 2) {
      if (force_bp) 
	sc = imx->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij;
      else 
	sc = (allow_bp)? imx->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule4: L -> a a' */
    else if (d == 2) {
      if (force_bp) 
	sc = p->t2[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t2[1] + emitsc_pairij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule5: L -> a */
    if (d == 1) {
      sc    = (allow_si)? p->t2[2] + emitsc_singi : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;

  case G6X_F:
    /* rule6: F -> a F a' */
    if (d > 2) {
      if (force_bp) 
	sc = imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
        
    /* rule7: F -> a a' */
    else if (d == 2) {
      if (force_bp) 
	sc = p->t3[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t3[1] + emitsc_pairij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule8: F -> LS */
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;
      
      sc    = imx->L->dp[k][d1] + imx->S->dp[j][d-d1] + p->t3[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  *ret_sc = sumsc;
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6x_outside(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx, int w, int j, int d,
			 SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    sumsc = -eslINFINITY;	/* sum score so far */
  SCVAL    sc;	        	        /* score for a rule */
  double   emitsc_pairij;
  int      allow_bp, force_bp;
  int      L = psq->n;
  int      im, jp;
  int      d1;
  int      i, k;
  int      status;

  if (d == L && w == G6X_S)  { *ret_sc = 0; return eslOK; }  // Initialization

  i  = j - d + 1;
  im = i - 1;
  jp = j + 1;
   
  // decide on constraints
  force_bp = (i > 1 && j < L)? force_bpair(im, jp, covct) : FALSE;
  allow_bp = (i > 1 && j < L)? allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, im, jp, L, covct, exclude, spair, ret_nneg) : FALSE;

  // emission scores
  emitsc_pairij = (i > 1 && j < L)? emitsc_pair_prof(im, jp, psq, p->e_pair) : -eslINFINITY;

  // Follow the grammar
  switch(w) {
  case G6X_S:
    // rule0: S -> LS 
    //
    //     L         S
    //  k_____i-1 i____j
    //  k______________j
    //          S
    //
    for (d1 = 0; d1 <= i; d1++) {
      k = i - d1;
      
      sc    = omx->S->dp[j][d1+d] + imx->L->dp[im][d1] + p->t1[0];
      sumsc = e2_FLogsum(sumsc, sc);
    }

    // rule8: F -> LS 
    for (d1 = 0; d1 <= i; d1++) {
      k = i - d1;
      
      sc    = omx->F->dp[j][d1+d] + imx->L->dp[im][d1] + p->t3[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;     

  case G6X_L:
    // rule0: S -> LS 
    //
    //     L        S
    //  i_____j j+1____k
    //  i______________k
    //          S
    //
     for (d1 = 0; d1 <= L-j; d1++) {
      k = j + d1;
      sc    = omx->S->dp[k][d1+d] + imx->S->dp[k][d1] + p->t1[0];
      sumsc = e2_FLogsum(sumsc, sc);

     }
    
    // rule1: S -> L
    sc    = omx->S->dp[j][d] + p->t1[1];
    sumsc = e2_FLogsum(sumsc, sc);
    
    // rule8: F -> LS 
    for (d1 = 0; d1 <= L-j; d1++) {
      k = j + d1;
      sc    = omx->F->dp[k][d1+d] + imx->S->dp[k][d1] + p->t3[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    break;
    
  case G6X_F:
    // rule3: L -> a F a'
    if (force_bp || allow_bp) {
      sc    = omx->L->dp[jp][d+2] + p->t2[0] + emitsc_pairij;
      sumsc = e2_FLogsum(sumsc, sc);
    }
   
    // rule6: F -> a F a'       
    if (force_bp || allow_bp) {
      sc    = omx->F->dp[jp][d+2] + p->t3[0] + emitsc_pairij;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  *ret_sc = sumsc;
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6x_posterior_single(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, int w, int j, 
				  POST *post, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    thisps;	  // add to the ps posterior
  SCVAL    ps;
  double   emitsc_singj; 
  int      allow_sj;
  int      L = psq->n;
  int      status;

  ps = post->ps[j];
  
 // decide on constraints
  allow_sj = allow_single(j, covct);

  // emission pair scores
  emitsc_singj  = emitsc_sing_prof(j, psq, p->e_sing);
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    break;     

  case G6X_L:
    /* rule5: L -> a */
    if (j > 0) {
      thisps = (allow_sj)? omx->L->dp[j][1] + p->t2[2] + emitsc_singj : -eslINFINITY;
      ps     = e2_FLogsum(ps, thisps);
    }
    break;
    
  case G6X_F:
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  post->ps[j] = ps;
  
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6x_posterior_pair(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, int w, int j, int d,
				POST *post, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    thispp;        // add to the pp posterior
  SCVAL    pp;
  double   emitsc_pairij;
  int      allow_bp, force_bp;
  int      L = psq->n;
  int      d1;
  int      i, k;
  int      status;

  i = j - d + 1;

  pp = post->pp[i][j];
  
 // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);

  // emission pair scores
  emitsc_pairij = emitsc_pair_prof(i, j, psq, p->e_pair);
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    break;     

  case G6X_L:
    /* rule3: L -> a F a' */
    if (force_bp || allow_bp) {
      thispp = omx->L->dp[j][d] + imx->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij;
      pp     = e2_FLogsum(pp, thispp);
    }
    
    /* rule4: L -> a a' */
    if (d == 2) {
      if (force_bp) 
	thispp =
	  omx->L->dp[j][d] + p->t2[1] + emitsc_pairij;
      else 
	thispp = (allow_bp)?
	  omx->L->dp[j][d] + p->t2[1] + emitsc_pairij : -eslINFINITY;

      pp = e2_FLogsum(pp, thispp);
    }
    break;
    
  case G6X_F:
    /* rule6: F -> a F a' */
    if (force_bp) 
      thispp =
	omx->F->dp[j][d] + imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij;
    else 
      thispp = (allow_bp)?
	omx->F->dp[j][d] + imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij : -eslINFINITY;
    pp = e2_FLogsum(pp, thispp);
        
    /* rule7: F -> a a' */
    if (d == 2) {
      if (force_bp) 
	thispp =
	  omx->F->dp[j][d] + p->t3[1] + emitsc_pairij;
      else 
	thispp = (allow_bp)?
	  omx->F->dp[j][d] + p->t3[1] + emitsc_pairij : -eslINFINITY;
      pp = e2_FLogsum(pp, thispp);
    }
    
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  post->pp[i][j] = pp;
  post->pp[j][i] = pp;										
  
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6xs_cyk(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX  *cyk, int w, int j, int d,
		      SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose)
{
  SCVAL    bestsc = -eslINFINITY;
  SCVAL    sc;
  double   emitsc_singi;
  double   emitsc_pairij;
  double   emitsc_stckij;
  int      allow_bp, force_bp;
  int      allow_si;
  int      L = psq->n;
  int      d1;
  int      i, k;
  int      status;
  
  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;

  if (d < 1 && w == G6X_L)  { *ret_sc = -eslINFINITY; return eslOK; }  // L  has at least 1 residues
  if (d < 1 && w == G6X_F)  { *ret_sc = -eslINFINITY; return eslOK; }  // F  has at least 1 residues

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  allow_si = allow_single(i, covct);

    // emission scores
  emitsc_singi  = emitsc_sing_prof(i,       psq, p->e_sing);
  emitsc_pairij = emitsc_pair_prof(i, j,    psq, p->e_pair);
  emitsc_stckij = emitsc_stck_prof(i, j, L, psq, p->e_pair, p->e_stck);
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    /* rule0: S -> LS */
    if (d > 0) {
      for (d1 = 0; d1 <= d; d1++) {
	k = i + d1 - 1;
	
	sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1] + p->t1[0];
	
	if (sc >= bestsc) {
	  if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	    if (alts) esl_stack_Reuse(alts);
	    bestsc = sc;
	  }     
	  if (alts) {
	    esl_stack_IPush(alts, G6X_S_1);
	    esl_stack_IPush(alts, d1);
	  }
	}
      }
    }

    /* rule1: S -> L */
    d1 = 0;
    sc = cyk->L->dp[j][d] + p->t1[1];
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_S_2);
	esl_stack_IPush(alts, d1);
      }
    }
   
    /* rule2: S -> epsilon */
    d1 = 0;
    if (d == 0) {
      sc = p->t1[2];
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_S_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    
    break;
    
  case G6X_L:
    /* rule3: L -> a F a' */
    d1 = 0;
    if (force_bp) 
      sc = cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij;
    else 
      sc = (allow_bp)?
	cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij : -eslINFINITY;
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_L_1);
	esl_stack_IPush(alts, d1);
      }
    }
    
    /* rule4: L -> a a' */
    d1 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->t2[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t2[1] + emitsc_pairij : -eslINFINITY;
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_L_2);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    
    /* rule5: L -> a */
    d1 = 0;
    if (d == 1) {
      sc = (allow_si)? p->t2[2] + emitsc_singi : -eslINFINITY;
	
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_L_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;

  case G6X_F:
    /* rule6: F -> a F a' */
    d1 = 0;
    if (force_bp) 
      sc = cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_pairij;     
    else 
      sc = (allow_bp)?
	cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_stckij : -eslINFINITY;
      
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, G6X_F_1);
	esl_stack_IPush(alts, d1);
      }
    }
    /* rule7: F -> a a' */
    d1 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->t3[1] + emitsc_stckij;     
      else 
	sc = (allow_bp)?
	  p->t3[1] + emitsc_stckij : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_F_2);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    /* rule8: F -> LS */
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;
      
      sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1] + p->t3[2];
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, G6X_F_3);
	  esl_stack_IPush(alts, d1);
	}
      }
    }
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }
 
  *ret_sc = bestsc;
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6xs_inside(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, int w, int j, int d,
			 SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    sumsc = -eslINFINITY;	/* sum score so far */
  SCVAL    sc;	        	        /* score for a rule */
  double   emitsc_singi;
  double   emitsc_pairij;
  double   emitsc_stckij;
  int      allow_bp, force_bp;
  int      allow_si;
  int      L = psq->n;
  int      d1;
  int      i, k;
  int      status;

  i = j - d + 1;

  if (d < 1 && w == G6X_L)  { *ret_sc = -eslINFINITY; return eslOK; }  // L  has at least 1 residues
  if (d < 1 && w == G6X_F)  { *ret_sc = -eslINFINITY; return eslOK; }  // F  has at least 1 residues

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  allow_si = allow_single(i, covct);
  
  // emission scores
  emitsc_singi  = emitsc_sing_prof(i,       psq, p->e_sing);
  emitsc_pairij = emitsc_pair_prof(i, j,    psq, p->e_pair);
  emitsc_stckij = emitsc_stck_prof(i, j, L, psq, p->e_pair, p->e_stck);
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    /* rule0: S -> LS */
    //     L        S
    //  i_____k k+1____j
    //  i______________j
    //          S
   
    if (d > 0) {
      for (d1 = 0; d1 <= d; d1++) {
	k = i + d1 - 1;
	
	sc    = imx->L->dp[k][d1] + imx->S->dp[j][d-d1] + p->t1[0];
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
    
    /* rule1: S -> L */
    sc    = imx->L->dp[j][d] + p->t1[1];
    sumsc = e2_FLogsum(sumsc, sc);
    
    /* rule2: S -> epsilon */
    if (d == 0) {
      sc    = p->t1[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;     

  case G6X_L:
    /* rule3: L -> a F a' */
    if (d > 2) {
      if (force_bp) 
	sc = imx->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij;
      else 
	sc = (allow_bp)? imx->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule4: L -> a a' */
    else if (d == 2) {
      if (force_bp) 
	sc = p->t2[1] + emitsc_pairij;
      else 
	sc = (allow_bp)?
	  p->t2[1] + emitsc_pairij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule5: L -> a */
    if (d == 1) {
      sc    = (allow_si)? p->t2[2] + emitsc_singi : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;

  case G6X_F:
    /* rule6: F -> a F a' */
    if (d > 2) {
      if (force_bp) 
	sc = imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_stckij;
      else 
	sc = (allow_bp)?
	  imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_stckij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
        
    /* rule7: F -> a a' */
    else if (d == 2) {
      if (force_bp) 
	sc = p->t3[1] + emitsc_stckij;
      else 
	sc = (allow_bp)?
	  p->t3[1] + emitsc_stckij : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule8: F -> LS */
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;
      
      sc    = imx->L->dp[k][d1] + imx->S->dp[j][d-d1] + p->t3[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  *ret_sc = sumsc;
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6xs_outside(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx, int w, int j, int d,
			  SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    sumsc = -eslINFINITY;	/* sum score so far */
  SCVAL    sc;	        	        /* score for a rule */
  double   emitsc_pairij;
  double   emitsc_stckij;
  int      allow_bp, force_bp;
  int      L = psq->n;
  int      im, jp;
  int      d1;
  int      i, k;
  int      status;

  if (d == L && w == G6X_S)  { *ret_sc = 0; return eslOK; }  // Initialization

  i  = j - d + 1;
  im = i - 1;
  jp = j + 1;
   
   // decide on constraints
  force_bp = (i > 1 && j < L)? force_bpair(im, jp, covct) : FALSE;
  allow_bp = (i > 1 && j < L)? allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, im, jp, L, covct, exclude, spair, ret_nneg) : FALSE;

  // emission scores
  emitsc_pairij = (i > 1 && j < L)? emitsc_pair_prof(im, jp,    psq, p->e_pair)            : -eslINFINITY;
  emitsc_stckij = (i > 1 && j < L)? emitsc_stck_prof(im, jp, L, psq, p->e_pair, p->e_stck) : -eslINFINITY;

  // Follow the grammar
  switch(w) {
  case G6X_S:
    // rule0: S -> LS 
    //
    //     L         S
    //  k_____i-1 i____j
    //  k______________j
    //          S
    //
    for (d1 = 0; d1 <= i; d1++) {
      k = i - d1;
      
      sc    = omx->S->dp[j][d1+d] + imx->L->dp[im][d1] + p->t1[0];
      sumsc = e2_FLogsum(sumsc, sc);
    }

    // rule8: F -> LS 
    for (d1 = 0; d1 <= i; d1++) {
      k = i - d1;
      
      sc    = omx->F->dp[j][d1+d] + imx->L->dp[im][d1] + p->t3[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;     

  case G6X_L:
    // rule0: S -> LS 
    //
    //     L        S
    //  i_____j j+1____k
    //  i______________k
    //          S
    //
     for (d1 = 0; d1 <= L-j; d1++) {
      k = j + d1;
      sc    = omx->S->dp[k][d1+d] + imx->S->dp[k][d1] + p->t1[0];
      sumsc = e2_FLogsum(sumsc, sc);

     }
    
    // rule1: S -> L
    sc    = omx->S->dp[j][d] + p->t1[1];
    sumsc = e2_FLogsum(sumsc, sc);
    
    // rule8: F -> LS 
    for (d1 = 0; d1 <= L-j; d1++) {
      k = j + d1;
      sc    = omx->F->dp[k][d1+d] + imx->S->dp[k][d1] + p->t3[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    break;
    
  case G6X_F:
    // rule3: L -> a F a'
    if (force_bp || allow_bp) {
      sc    = omx->L->dp[jp][d+2] + p->t2[0] + emitsc_pairij;
      sumsc = e2_FLogsum(sumsc, sc);
    }
   
    // rule6: F -> a F a'       
    if (force_bp || allow_bp) {
      sc    = omx->F->dp[jp][d+2] + p->t3[0] + emitsc_stckij;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  *ret_sc = sumsc;
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6xs_posterior_single(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, int w, int j, 
				   POST *post, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    thisps;	  // add to the ps posterior
  SCVAL    ps;
  double   emitsc_singj; 
  int      allow_sj;
  int      L = psq->n;
  int      status;

  ps = post->ps[j];
  
 // decide on constraints
  allow_sj = allow_single(j, covct);

  // emission pair scores
  emitsc_singj  = emitsc_sing_prof(j, psq, p->e_sing);
  
  // Follow the grammar
  switch(w) {
  case G6X_S:
    break;     

  case G6X_L:
    /* rule5: L -> a */
    if (j > 0) {
      thisps = (allow_sj)? omx->L->dp[j][1] + p->t2[2] + emitsc_singj : -eslINFINITY;
      ps     = e2_FLogsum(ps, thisps);
    }
    break;
    
  case G6X_F:
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  post->ps[j] = ps;
  
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_g6xs_posterior_pair(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, int w, int j, int d,
				 POST *post, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    thispp;        // add to the pp posterior
  SCVAL    pp;
  double   emitsc_pairij;
  double   emitsc_stckij;
  int      allow_bp, force_bp;
  int      L = psq->n;
  int      d1;
  int      i, k;
  int      status;

  i = j - d + 1;

  pp = post->pp[i][j];
  
 // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);

  // emission pair scores
  emitsc_pairij = emitsc_pair_prof(i, j,    psq, p->e_pair);
  emitsc_stckij = emitsc_stck_prof(i, j, L, psq, p->e_pair, p->e_stck);

  // Follow the grammar
  switch(w) {
  case G6X_S:
    break;     

  case G6X_L:
    /* rule3: L -> a F a' */
    if (force_bp || allow_bp) {
      thispp = omx->L->dp[j][d] + imx->F->dp[j-1][d-2] + p->t2[0] + emitsc_pairij;
      pp     = e2_FLogsum(pp, thispp);
    }
    
    /* rule4: L -> a a' */
    if (d == 2) {
      if (force_bp) 
	thispp =
	  omx->L->dp[j][d] + p->t2[1] + emitsc_pairij;
      else 
	thispp = (allow_bp)?
	  omx->L->dp[j][d] + p->t2[1] + emitsc_pairij : -eslINFINITY;

      pp = e2_FLogsum(pp, thispp);
    }
    break;
    
  case G6X_F:
    /* rule6: F -> a F a' */
    if (force_bp) 
      thispp =
	omx->F->dp[j][d] + imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_stckij;
    else 
      thispp = (allow_bp)?
	omx->F->dp[j][d] + imx->F->dp[j-1][d-2] + p->t3[0] + emitsc_stckij : -eslINFINITY;
    pp = e2_FLogsum(pp, thispp);
        
    /* rule7: F -> a a' */
    if (d == 2) {
      if (force_bp) 
	thispp =
	  omx->F->dp[j][d] + p->t3[1] + emitsc_stckij;
      else 
	thispp = (allow_bp)?
	  omx->F->dp[j][d] + p->t3[1] + emitsc_pairij : -eslINFINITY;
      pp = e2_FLogsum(pp, thispp);
    }
    
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize G6X nt %d\n", w);
  }

  post->pp[i][j] = pp;
  post->pp[j][i] = pp;										
  
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_rbg_cyk(FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *cyk, int w, int j, int d,
		     SCVAL *ret_sc, int *ret_nneg, ESL_STACK *alts, char *errbuf, int verbose)
{
  SCVAL    bestsc = -eslINFINITY;
  SCVAL    sc;
  double   emitsc_singi, emitsc_singj;
  double   emitsc_pair1, emitsc_pair2;
  double   emitsc_stck1, emitsc_stck2;
  int      allow_si, allow_sj;
  int      allow_bp, force_bp;
  int      allow_hp;
  int      L = psq->n;
  int      d1, d2;
  int      i, k, l;
  int      d_ng, d1_ng, d2_ng;
  int      status;

  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;

  if (d < 2 && w == RBG_M)  { *ret_sc = -eslINFINITY; return eslOK; }  // M  has at least 2 residues
  if (d < 2 && w == RBG_F0) { *ret_sc = -eslINFINITY; return eslOK; }  // F0 has at least 2 residues
  if (d < 2 && w == RBG_F5) { *ret_sc = -eslINFINITY; return eslOK; }  // F5 has at least 2 residues
  if (d < 2 && w == RBG_R)  { *ret_sc = -eslINFINITY; return eslOK; }  // R  has at least 2 residues
  if (d < 2 && w == RBG_M1) { *ret_sc = -eslINFINITY; return eslOK; }  // M1 has at least 2 residues

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  allow_hp = allow_hairpin(foldparam->hloop_min, i, j, L, covct);
  allow_si = allow_single(i, covct);
  allow_sj = allow_single(j, covct);
  
  // special case
  if (force_bp && w == RBG_P) { *ret_sc = -eslINFINITY; return eslOK; }  // P does not allow i-j pairing
			
  // emission scores
  emitsc_singi  = emitsc_sing_prof(i,       psq, p->e_sing);
  emitsc_singj  = emitsc_sing_prof(j,       psq, p->e_sing);
  emitsc_pair1  = emitsc_pair_prof(i, j,    psq, p->e_pair1);
  emitsc_pair2  = emitsc_pair_prof(i, j,    psq, p->e_pair2);
  emitsc_stck1  = emitsc_stck_prof(i, j, L, psq, p->e_pair1, p->e_stck1);
  emitsc_stck2  = emitsc_stck_prof(i, j, L, psq, p->e_pair2, p->e_stck2);
  
  // Follow the grammar
  // order: S F0 F5 P M1 R M
  switch(w) {
  case RBG_S:
    /* rule0: S -> a S */
    d1 = d2 = 0;
    sc = -eslINFINITY;
    
    if (d > 0) {
      sc = (allow_si)? cyk->S->dp[j][d-1] + p->tS[0] + emitsc_singi : -eslINFINITY;
	
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_S_1);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
    
    /* rule1: S -> F0 S */
    d2 = 0;
    for (d1 = 0; d1 <= d; d1++) {
      
      k = i + d1 - 1;
      
      sc = cyk->F0->dp[k][d1] + cyk->S->dp[j][d-d1] + p->tS[1];
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_S_2);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
    /* rule2: S -> epsilon */
    d1 = d2 = 0;
    if (d == 0) {
      sc = p->tS[2];
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}
	
	if (alts) {
	  esl_stack_IPush(alts, RBG_S_3);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
    break;
    
  case RBG_F0:
    /* rule3: F0 -> a F5 a' */
    d1 = d2 = 0;
   if (force_bp) 
      sc = cyk->F5->dp[j-1][d-2] + p->tF0[0] + emitsc_pair1;
    else 
      sc = (allow_bp)?
	cyk->F5->dp[j-1][d-2] + p->tF0[0] + emitsc_pair1 : -eslINFINITY;
   
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     	
      if (alts) {
	esl_stack_IPush(alts, RBG_F0_1);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }
  
    /* rule4: F0 -> a P a' */
    d1 = d2 = 0;
    if (force_bp) 
      sc = cyk->P->dp[j-1][d-2] + p->tF0[1] + emitsc_pair2;     
     else 
       sc = (allow_bp)?
	 cyk->P->dp[j-1][d-2] + p->tF0[1] + emitsc_pair2 : -eslINFINITY;
    
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     	
      if (alts) {
	esl_stack_IPush(alts, RBG_F0_2);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }
  
    /* rule5: F0 -> a a' */
    d1 = d2 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->tF0[2] + emitsc_pair2;     
      else 
	sc = (allow_bp)? p->tF0[2] + emitsc_pair2 : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     	
	if (alts) {
	  esl_stack_IPush(alts, RBG_F0_3);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
  
    break;

  case RBG_F5:
    /* rule6: F5 -> a F5^{bb'} a' */
    d1 = d2 = 0;
    if (force_bp) 
      sc = cyk->F5->dp[j-1][d-2] + p->tF5[0] + emitsc_stck1;
    else 
      sc = (allow_bp)?
	cyk->F5->dp[j-1][d-2] + p->tF5[0] + emitsc_stck1 : -eslINFINITY;
    
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     	
      if (alts) {
	esl_stack_IPush(alts, RBG_F5_1);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }
    
    /* rule7: F5 -> a P^{bb'} a' */
    d1 = d2 = 0;
    if (force_bp) 
      sc = cyk->P->dp[j-1][d-2] + p->tF5[1] + emitsc_stck2;
    else 
      sc = (allow_bp)?
	cyk->P->dp[j-1][d-2] + p->tF5[1] + emitsc_stck2 : -eslINFINITY;
  
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     	
      if (alts) {
	esl_stack_IPush(alts, RBG_F5_2);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }
    break;
    
    /* rule8: F5 -> a a' */
    d1 = d2 = 0;
    if (d == 2) {
      if (force_bp) 
	sc = p->tF5[2] + emitsc_stck2;
      else 
	sc = (allow_bp)?
	  p->tF5[2] + emitsc_stck2 : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     	
	if (alts) {
	  esl_stack_IPush(alts, RBG_F5_3);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
    break;

  case RBG_P:
    /* rule9: P -> m..m */
    d1 = d2 = 0;
    if (d >= MAXLOOP_H) {
      sc = -eslINFINITY;
    }
    else {
      d_ng = segment_remove_gaps_prof(i,j,psq); if (d_ng == 0) d_ng = d;
      sc = (allow_hp)? p->tP[0] + p->l1[d_ng-1] + score_loop_hairpin_prof(i, j, p, psq) : -eslINFINITY;
    }
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, RBG_P_1);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }
    
    /* rule10: P -> m..m F0 */
    d2 = 0;
    for (d1 = 1; d1 <= ESL_MIN(d,MAXLOOP_B); d1++) {
      
      k = i + d1 - 1;
      
      d1_ng = segment_remove_gaps_prof(i,k,psq); if (d1_ng == 0) d1_ng = d1;
      sc = allow_loop(i, k, covct)? cyk->F0->dp[j][d-d1] + p->tP[1] + p->l2[d1_ng-1] + score_loop_bulge_prof(i, k, p, psq) : -eslINFINITY;
	
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_P_2);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}	
      }
    }
    
    /* rule11: P -> F0 m..m */
    d1 = 0;
    for (d2 = 1; d2 <= ESL_MIN(d,MAXLOOP_B); d2++) {
      
      l = j - d2 + 1;
      
      d2_ng = segment_remove_gaps_prof(l,j,psq); if (d2_ng == 0) d2_ng = d2;
      sc = allow_loop(l, j, covct)? cyk->F0->dp[l-1][d-d2] + p->tP[2] + p->l2[d2_ng-1] + score_loop_bulge_prof(l, j, p, psq) : -eslINFINITY;

      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_P_3);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
    
    /* rule12: P -> m..m F0 m..m */
    for (d1 = 1; d1 <= ESL_MIN(d,MAXLOOP_I); d1++) {
      for (d2 = 1; d2 <= ESL_MIN(d-d1,MAXLOOP_I); d2++) {
	
	if (d1 + d2 > MAXLOOP_I) break;

	k = i + d1 - 1;
	l = j - d2 + 1;

	d1_ng = segment_remove_gaps_prof(i,k,psq); if (d1_ng == 0) d1_ng = d1;
	d2_ng = segment_remove_gaps_prof(l,j,psq); if (d2_ng == 0) d2_ng = d2;
	sc = (l > 0 && allow_loop(i,k,covct) && allow_loop(l,j,covct))?
	    cyk->F0->dp[l-1][d-d1-d2] + p->tP[3] + p->l3[d1_ng-1][d2_ng-1] + score_loop_intloop_prof(i, k, p, psq) + score_loop_intloop_prof(l, j, p, psq) : -eslINFINITY;
	  
	if (sc >= bestsc) {
	  if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	    if (alts) esl_stack_Reuse(alts);
	    bestsc = sc;
	  }     
	  if (alts) {
	    esl_stack_IPush(alts, RBG_P_4);
	    esl_stack_IPush(alts, d1);
	    esl_stack_IPush(alts, d2);
	  }
	}
      }
    }

    /* rule13: P -> M1 M */
    d2 = 0;
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;
      
      sc = cyk->M1->dp[k][d1] + cyk->M->dp[j][d-d1] + p->tP[4];
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_P_5);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
    break;
    
  case RBG_M:
    d2 = 0;
    /* rule14: M -> M1 M */
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;
      
      sc = cyk->M1->dp[k][d1] + cyk->M->dp[j][d-d1] + p->tM[0];
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_M_1);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
  
    /* rule15: M -> R */
    d1 = d2 = 0;
    sc = cyk->R->dp[j][d] + p->tM[1];
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, RBG_M_2);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }
    break;
    
  case RBG_R:
    /* rule16: R -> R a */
    d1 = d2 = 0;
    if (d > 0) {
      sc = (allow_sj)? cyk->R->dp[j-1][d-1] + p->tR[0] + emitsc_singj : -eslINFINITY;
	
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_R_1);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }

    /* rule17: R -> M1 */
    d1 = d2 = 0;
    sc = cyk->M1->dp[j][d] + p->tR[1];
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, RBG_R_2);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }

    break;
    
  case RBG_M1:
    /* rule18: M1 -> a M1 */
    d1 = d2 = 0;
    if (d > 0) {
      sc = (allow_si)? cyk->M1->dp[j][d-1] + p->tM1[0] + emitsc_singi : -eslINFINITY;
      
      if (sc >= bestsc) {
	if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	  if (alts) esl_stack_Reuse(alts);
	  bestsc = sc;
	}     
	if (alts) {
	  esl_stack_IPush(alts, RBG_M1_1);
	  esl_stack_IPush(alts, d1);
	  esl_stack_IPush(alts, d2);
	}
      }
    }
 
    /* rule19: M1 -> F0 */
    d1 = d2 = 0;
    sc = cyk->F0->dp[j][d] + p->tM1[1];
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, RBG_M1_2);
	esl_stack_IPush(alts, d1);
	esl_stack_IPush(alts, d2);
      }
    }
  
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize RBG nt %d\n", w);
    
  }
  
  *ret_sc = bestsc;

  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_rbg_inside(FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx, int w, int j, int d,
			 SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    sumsc = -eslINFINITY;	/* sum score so far */
  SCVAL    sc;	        	        /* score for a rule */
  double   emitsc_singi, emitsc_singj;
  double   emitsc_pair1, emitsc_pair2;
  double   emitsc_stck1, emitsc_stck2;
  int      allow_bp, force_bp;
  int      allow_hp;
  int      allow_si, allow_sj;
  int      L = psq->n;
  int      d1, d2;
  int      i, k, l;
  int      d_ng, d1_ng, d2_ng;
  int      status;

  i = j - d + 1;

  if (d < 2 && w == RBG_M)  { *ret_sc = -eslINFINITY; return eslOK; }  // M  has at least 2 residues
  if (d < 2 && w == RBG_F0) { *ret_sc = -eslINFINITY; return eslOK; }  // F0 has at least 2 residues
  if (d < 2 && w == RBG_F5) { *ret_sc = -eslINFINITY; return eslOK; }  // F5 has at least 2 residues
  if (d < 2 && w == RBG_R)  { *ret_sc = -eslINFINITY; return eslOK; }  // R  has at least 2 residues
  if (d < 2 && w == RBG_M1) { *ret_sc = -eslINFINITY; return eslOK; }  // M1 has at least 2 residues

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  allow_hp = allow_hairpin(foldparam->hloop_min, i, j, L, covct);
  allow_si = allow_single(i, covct);
  allow_sj = allow_single(j, covct);

  // special case
  if (force_bp && w == RBG_P) { *ret_sc = -eslINFINITY; return eslOK; }  // P does not allow i-j pairing

  // emission scores
  emitsc_singi  = emitsc_sing_prof(i,       psq, p->e_sing);
  emitsc_singj  = emitsc_sing_prof(j,       psq, p->e_sing);
  emitsc_pair1  = emitsc_pair_prof(i, j,    psq, p->e_pair1);
  emitsc_pair2  = emitsc_pair_prof(i, j,    psq, p->e_pair2);
  emitsc_stck1  = emitsc_stck_prof(i, j, L, psq, p->e_pair1, p->e_stck1);
  emitsc_stck2  = emitsc_stck_prof(i, j, L, psq, p->e_pair2, p->e_stck2);
  
  // Follow the grammar
  switch(w) {
  case RBG_S:
    /* rule0: S -> a S */ 
    //  a      S
    //  i i+1_____j
    //  i_________j
    //       S
     if (d > 0) {
      sc = (allow_si)? imx->S->dp[j][d-1] + p->tS[0] + emitsc_singi : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule1: S -> F0 S */
    //    F0        S
    //  i_____k k+1____j
    //  i______________j
    //          S
    //
    for (d1 = 0; d1 <= d; d1++) {
      
      k = i + d1 - 1;
      
      sc = imx->F0->dp[k][d1] + imx->S->dp[j][d-d1] + p->tS[1];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    /* rule2: S -> epsilon */
    if (d == 0) {
      sc    = p->tS[2];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;
    
  case RBG_F0:
    /* rule3: F0 -> a F5 a' */
    //             
    //      a      F5      a'
    //      i i+1______j-1 j
    //      i______________j
    //             F0
    //
    if (d > 2) {
      if (force_bp || allow_bp) {
	sc = imx->F5->dp[j-1][d-2] + p->tF0[0] + emitsc_pair1; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
	
    /* rule4: F0 -> a P a' */
    //             
    //      a      P       a'
    //      i i+1______j-1 j
    //      i______________j
    //             F0
    //
    if (d > 2) {
      if (force_bp || allow_bp) {
	sc = imx->P->dp[j-1][d-2] + p->tF0[1] + emitsc_pair2; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
      
    /* rule5: F0 -> a a' */
    if (d == 2) {
      if (force_bp || allow_bp) {
	sc =  p->tF0[2] + emitsc_pair2; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
    
    break;
    
  case RBG_F5:
    /* rule6: F5 -> a F5^{bb'} a' */
    //             
    //      a      F5       a'
    //      i i+1______j-1 j
    //      i______________j
    //             F5
    //
    if (d > 2) {
      if (force_bp || allow_bp) {
	sc = imx->F5->dp[j-1][d-2] + p->tF5[0] + emitsc_stck1; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
        
    /* rule7: F5 -> a P^{bb'} a' */
    //             
    //      a      P       a'
    //      i i+1______j-1 j
    //      i______________j
    //             F5
    //
    if (d > 2) {
      if (force_bp || allow_bp) {
	sc = imx->P->dp[j-1][d-2] + p->tF5[1] + emitsc_stck2; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
   break;
    
    /* rule8: F5 -> a a' */
    if (d == 2) {
      if (force_bp || allow_bp) {
	sc = p->tF5[2] + emitsc_stck2; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
     }
    break;

  case RBG_P:
    /* rule9: P -> m..m */
    if (d < MAXLOOP_H) {
      d_ng = segment_remove_gaps_prof(i,j,psq); if (d_ng == 0) d_ng = d;
      sc = (allow_hp)? p->tP[0] + p->l1[d_ng-1] + score_loop_hairpin_prof(i, j, p, psq) : -eslINFINITY;
    }
    sumsc = e2_FLogsum(sumsc, sc); 

    /* rule10: P -> m..m F0 */
    //
    //  m......m      F0
    //  i______k k+1_____j
    //  i________________j
    //          P
    //
    d2 = 0;
    for (d1 = 1; d1 <= ESL_MIN(d,MAXLOOP_B); d1++) {
      
      k = i + d1 - 1;
      
      d1_ng = segment_remove_gaps_prof(i,k,psq); if (d1_ng == 0) d1_ng = d1;
      sc = allow_loop(i, k, covct)? imx->F0->dp[j][d-d1] + p->tP[1] + p->l2[d1_ng-1] + score_loop_bulge_prof(i, k, p, psq) : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc); 
    }
    
    /* rule11: P -> F0 m..m */
    //
    //     F0     m......m
    //  i_____l-1 l______j
    //  i________________j
    //          P
    //
    d1 = 0;
    for (d2 = 1; d2 <= ESL_MIN(d,MAXLOOP_B); d2++) {
      
      l = j - d2 + 1;
      
      d2_ng = segment_remove_gaps_prof(l,j,psq); if (d2_ng == 0) d2_ng = d2;
      sc = allow_loop(l, j, covct)? imx->F0->dp[l-1][d-d2] + p->tP[2] + p->l2[d2_ng-1] + score_loop_bulge_prof(l, j, p, psq) : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc); 
    }
    
    /* rule12: P -> m..m F0 m..m */
    //
    //  m.....m      F0      m.....m
    //  i_____k k+1______l-1 l_____j
    //  i__________________________j
    //                P
    //
    for (d1 = 1; d1 <= ESL_MIN(d,MAXLOOP_I); d1++) {
      for (d2 = 1; d2 <= ESL_MIN(d-d1,MAXLOOP_I); d2++) {
	
	if (d1 + d2 > MAXLOOP_I) break;

	k = i + d1 - 1;
	l = j - d2 + 1;

	d1_ng = segment_remove_gaps_prof(i,k,psq); if (d1_ng == 0) d1_ng = d1;
	d2_ng = segment_remove_gaps_prof(l,j,psq); if (d2_ng == 0) d2_ng = d2;
	sc = (l > 0 && allow_loop(i,k,covct) && allow_loop(l,j,covct))?
	    imx->F0->dp[l-1][d-d1-d2] + p->tP[3] + p->l3[d1_ng-1][d2_ng-1] + score_loop_intloop_prof(i, k, p, psq) + score_loop_intloop_prof(l, j, p, psq) : -eslINFINITY;
	sumsc = e2_FLogsum(sumsc, sc); 
      }
    }

    /* rule13: P -> M1 M */
    //
    //     M1         M
    //  i_____k k+1______j
    //  i________________j
    //          P
    //
    d2 = 0;
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;
      
      sc = imx->M1->dp[k][d1] + imx->M->dp[j][d-d1] + p->tP[4];
      sumsc = e2_FLogsum(sumsc, sc);       
    }
    break;
    
  case RBG_M:
    d2 = 0;
    /* rule14: M -> M1 M */
    //
    //     M1         M
    //  i_____k k+1______j
    //  i________________j
    //          M
    //
    for (d1 = 0; d1 <= d; d1++) {
      k = i + d1 - 1;
      
      sc = imx->M1->dp[k][d1] + imx->M->dp[j][d-d1] + p->tM[0];
      sumsc = e2_FLogsum(sumsc, sc);             
    }
  
    /* rule15: M -> R */
    sc    = imx->R->dp[j][d] + p->tM[1];
    sumsc = e2_FLogsum(sumsc, sc);                 
    break;
    
  case RBG_R:
    /* rule16: R -> R a */
    //
    //      R     a
    //  i_____j-1 j
    //  i_________j
    //        R
    //
    if (d > 0) {
      sc    = (allow_sj)? imx->R->dp[j-1][d-1] + p->tR[0] + emitsc_singj : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);                 	
    }

    /* rule17: R -> M1 */
    sc    = imx->M1->dp[j][d] + p->tR[1];
    sumsc = e2_FLogsum(sumsc, sc);                 	
    break;
    
  case RBG_M1:
    /* rule18: M1 -> a M1 */
    //
    //  a     M1
    //  i i+1_____j
    //  i_________j
    //       M1
    //
    if (d > 0) {
      sc    = (allow_si)? imx->M1->dp[j][d-1] + p->tM1[0] + emitsc_singi : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc);                 	    
    }
 
    /* rule19: M1 -> F0 */
    sc    = imx->F0->dp[j][d] + p->tM1[1];
    sumsc = e2_FLogsum(sumsc, sc);                 	        
  
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize RBG nt %d\n", w);
    
  }

  *ret_sc = sumsc;
  return eslOK;

 ERROR:
  return status;
}

static int 
dp_recursion_rbg_outside(FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *omx, RBG_MX *imx, int w, int j, int d,
			  SCVAL *ret_sc, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    sumsc = -eslINFINITY;	/* sum score so far */
  SCVAL    sc;	        	        /* score for a rule */
  double   emitsc_pair1, emitsc_pair2;
  double   emitsc_stck1, emitsc_stck2;
  double   emitsc_singi, emitsc_singj;
  int      allow_bp, force_bp;
  int      allow_si, allow_sj;
  int      L = psq->n;
  int      im, jp;
  int      d1, d2;
  int      i, k, l;
  int      d_ng, d1_ng, d2_ng;
  int      status;

  if (d == L && w == RBG_S)  { *ret_sc = 0;            return eslOK; }  // Initialization
  if (d == L && w == RBG_P)  { *ret_sc = -eslINFINITY; return eslOK; }  // 
  if (d == L && w == RBG_M1) { *ret_sc = -eslINFINITY; return eslOK; }  // 
  if (d == L && w == RBG_R)  { *ret_sc = -eslINFINITY; return eslOK; }  // 
  if (d == L && w == RBG_M)  { *ret_sc = -eslINFINITY; return eslOK; }  // 
  if (d == L && w == RBG_F5) { *ret_sc = -eslINFINITY; return eslOK; }  // 

  i  = j - d + 1;
  im = i - 1;
  jp = j + 1;

  // special case
  if (force_bpair(i, j, covct) && w == RBG_P) { *ret_sc = -eslINFINITY; return eslOK; }  // P does not allow i-j pairing
  
  // decide on constraints
  force_bp = (i > 1 && j < L)? force_bpair(im, jp, covct) : -eslINFINITY;
  allow_bp = (i > 1 && j < L)? allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, im, jp, L, covct, exclude, spair, ret_nneg) : -eslINFINITY;
  allow_si = (i > 1)? allow_single(im, covct) : -eslINFINITY;
  allow_sj = (j < L)? allow_single(jp, covct) : -eslINFINITY;
  
  // emission scores
  emitsc_singi  = (i > 1)?          emitsc_sing_prof(im,        psq, p->e_sing)  : -eslINFINITY;
  emitsc_singj  = (j < L)?          emitsc_sing_prof(jp,        psq, p->e_sing)  : -eslINFINITY;
  emitsc_pair1  = (i > 1 && j < L)? emitsc_pair_prof(im, jp,    psq, p->e_pair1) : -eslINFINITY;
  emitsc_pair2  = (i > 1 && j < L)? emitsc_pair_prof(im, jp,    psq, p->e_pair2) : -eslINFINITY;
  emitsc_stck1  = (i > 1 && j < L)? emitsc_stck_prof(im, jp, L, psq, p->e_pair1, p->e_stck1) : -eslINFINITY;
  emitsc_stck2  = (i > 1 && j < L)? emitsc_stck_prof(im, jp, L, psq, p->e_pair2, p->e_stck2) : -eslINFINITY;
 
  // Follow the grammar
  switch(w) { // order: 
  case RBG_S:
    // rule0: S -> aS 
    //
    //   a    S
    //  i-1 i____j
    //  i-1______j
    //       S
    if (d < L && allow_si) {
      sc    = omx->S->dp[j][d+1] + emitsc_singi + p->tS[0];
      sumsc = e2_FLogsum(sumsc, sc);
    }

    // rule1: S -> F0 S 
    //
    //    F0        S
    //  k_____i-1 i____j
    //  k______________j
    //          S
    //
    for (d1 = 1; d1 <= i; d1++) {
      k = i - d1;
      
      sc    = omx->S->dp[j][d1+d] + imx->F0->dp[im][d1] + p->tS[1];
      sumsc = e2_FLogsum(sumsc, sc);
    }
   
    break;
    
  case RBG_F0:
    // rule1: S -> F0 S
    //
    //     F0       S
    //  i_____j j+1____l
    //  i______________l
    //          S
    for (d2 = 0; d2 <= L-j; d2++) {
      l = j + d2;
      sc    = omx->S->dp[l][d2+d] + imx->S->dp[l][d2] + p->tS[1];
      sumsc = e2_FLogsum(sumsc, sc);
     }
     
    // rule10: P -> m..m F0 
    //
    //  m.....m     F0
    //  k_____i-1 i____j
    //  k______________j
    //          P
    //
    for (d1 = 1; d1 <= ESL_MIN(i-1,MAXLOOP_B); d1++) {
      
      k = i - d1;
      
      d1_ng = segment_remove_gaps_prof(k,im,psq); if (d1_ng == 0) d1_ng = d1;
      sc = allow_loop(k, im, covct)? omx->P->dp[j][d+d1] + p->tP[1] + p->l2[d1_ng-1] + score_loop_bulge_prof(k, im, p, psq) : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc); 
    }

    // rule11: P -> F0 m..m
    //
    //     F0   m......m
    //  i_____j j+1____l
    //  i______________l
    //          P
    //
    for (d2 = 1; d2 <= ESL_MIN(L-j,MAXLOOP_B); d2++) {
      
      l = j + d2;
      
      d2_ng = segment_remove_gaps_prof(jp,l,psq); if (d2_ng == 0) d2_ng = d2;
      sc = allow_loop(jp, l, covct)? omx->P->dp[l][d+d2] + p->tP[2] + p->l2[d2_ng-1] + score_loop_bulge_prof(jp, l, p, psq) : -eslINFINITY;
      sumsc = e2_FLogsum(sumsc, sc); 
    }

    // rule12: P -> m..m F0 m..m
    //
    //  m....m     F0    m.....m
    //  k___i-1 i_____j j+1____l
    //  k______________________l
    //             P
     for (d1 = 1; d1 <= ESL_MIN(i-1,MAXLOOP_I); d1++) {
      for (d2 = 1; d2 <= ESL_MIN(L-j,MAXLOOP_I); d2++) {
	
	if (d1 + d2 > MAXLOOP_I) break;

	k = i - d1;
	l = j + d2;

	d1_ng = segment_remove_gaps_prof(k,im,psq); if (d1_ng == 0) d1_ng = d1;
	d2_ng = segment_remove_gaps_prof(jp,l,psq); if (d2_ng == 0) d2_ng = d2;
	sc = (allow_loop(k,im,covct) && allow_loop(jp,l,covct))?
	    omx->P->dp[l][d+d1+d2] + p->tP[3] + p->l3[d1_ng-1][d2_ng-1] + score_loop_intloop_prof(k, im, p, psq) + score_loop_intloop_prof(jp, l, p, psq) : -eslINFINITY;
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }

    // rule19: M1 -> F0
    //
    sc    = omx->M1->dp[j][d] + p->tM1[1];
    sumsc = e2_FLogsum(sumsc, sc);
 
    break;
    
  case RBG_F5:
    // rule3: F0 -> a F5 a'
    //             
    //      a      F5    a'
    //      i-1 i______j j+1
    //      i-1__________j+1
    //             F0
    //
    if (i > 1 && j < L && d > 2) {
      if (force_bp || allow_bp) {
	sc = omx->F0->dp[jp][d+2] + p->tF0[0] + emitsc_pair1; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
    
    // rule6: F5 -> a F5 a'
    //             
    //      a      F5    a'
    //      i-1 i______j j+1
    //      i-1__________j+1
    //             F5
    //
    if (i > 1 && j < L && d > 2) {
      if (force_bp || allow_bp) {
	sc = omx->F5->dp[jp][d+2] + p->tF5[0] + emitsc_stck1; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
 
    break;
    
  case RBG_P:
    // rule4: F0 -> a P a'
    //
    //      a       P    a'
    //      i-1 i______j j+1
    //      i-1__________j+1
    //             F0
    //
    if (i > 1 && j < L && d > 2) {
      if (force_bp || allow_bp) {
	sc = omx->F0->dp[jp][d+2] + p->tF0[1] + emitsc_pair2; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
    
    // rule7: F5 -> a P a'
    //
    //      a      P    a'
    //      i-1 i______j j+1
    //      i-1__________j+1
    //             F5
    //
    if (i > 1 && j < L && d > 2) {
      if (force_bp || allow_bp) {
	sc = omx->F5->dp[jp][d+2] + p->tF5[1] + emitsc_stck2; 
	sumsc = e2_FLogsum(sumsc, sc);
      }
    }
    
    break;
    
  case RBG_M:
    // rule13: P -> M1 M
    //
    //     M1        M
    //  k_____i-1 i____j
    //  k______________j
    //          P
    //
    for (d1 = 0; d1 <= i; d1++) {
      k = i - d1;
      
      sc    = omx->P->dp[j][d1+d] + imx->M1->dp[im][d1] + p->tP[4];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    
    // rule14: M -> M1 M
    //
    //     M1        M
    //  k_____i-1 i____j
    //  k______________j
    //          M
    //
    for (d1 = 0; d1 <= i; d1++) {
      k = i - d1;
      
      sc    = omx->M->dp[j][d1+d] + imx->M1->dp[im][d1] + p->tM[0];
      sumsc = e2_FLogsum(sumsc, sc);
    }
  
    break;
    
  case RBG_R:
    // rule15: M -> R
    //
    sc    = omx->M->dp[j][d] + p->tM[1];
    sumsc = e2_FLogsum(sumsc, sc);

    // rule16: R -> R a
    //
    //      R     a
    //  i______j j+1
    //  i________j+1
    //       R
    if (j < L && allow_sj) {
      sc    = omx->R->dp[jp][d+1] + emitsc_singj + p->tR[0];
      sumsc = e2_FLogsum(sumsc, sc);
    }
    break;
    
  case RBG_M1:
    // rule13: P -> M1 M
    //
    //     M1        M
    //  i_____j j+1____l
    //  i______________l
    //          P
    //
    for (d2 = 1; d2 <= L-j; d2++) {
      l = j + d2;
      sc    = omx->P->dp[l][d2+d] + imx->M->dp[l][d2] + p->tP[4];
      sumsc = e2_FLogsum(sumsc, sc);
     }


    // rule14: M -> M1 M
    //
    //     M1        M
    //  i_____j j+1____l
    //  i______________l
    //          M
    //
     for (d2 = 1; d2 <= L-j; d2++) {
      l = j + d2;
      sc    = omx->M->dp[l][d2+d] + imx->M->dp[l][d2] + p->tM[0];
      sumsc = e2_FLogsum(sumsc, sc);
     }

    // rule17: R -> M1
    //
    sc    = omx->R->dp[j][d] + p->tR[1];
    sumsc = e2_FLogsum(sumsc, sc);

    // rule18: M1 -> a M1
    //
    //   a    M1
    //  i-1 i____j
    //  i-1______j
    //       M1
    if (d < L && allow_si) {
      sc    = omx->M1->dp[j][d+1] + emitsc_singi + p->tM1[0];
      sumsc = e2_FLogsum(sumsc, sc);
    }
  
    break;
    
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize RBG nt %d\n", w);
  }

  *ret_sc = sumsc;
  return eslOK;

 ERROR:
  return status;
}


static int 
dp_recursion_rbg_posterior_pair(FOLDPARAM *foldparam, RBGparam *p, PSQ *psq, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx, RBG_MX *omx, int w, int j, int d,
				 POST *post, int *ret_nneg, char *errbuf, int verbose)
{
  SCVAL    thispp;        // add to the pp posterior
  SCVAL    pp;
  double   emitsc_pair1, emitsc_pair2;
  double   emitsc_stck1, emitsc_stck2;
  int      allow_bp, force_bp;
  int      L = psq->n;
  int      d1, d2;
  int      i, k, l;
  int      d_ng, d1_ng, d2_ng;
  int      status;

  i = j - d + 1;

  pp = post->pp[i][j];

  // decide on constraints
  force_bp = force_bpair(i, j, covct);
  allow_bp = allow_bpair(foldparam->power_thresh, foldparam->neg_eval_thresh, foldparam->hloop_min, i, j, L, covct, exclude, spair, ret_nneg);
  
  // emission scores
  emitsc_pair1  = emitsc_pair_prof(i, j,    psq, p->e_pair1);
  emitsc_pair2  = emitsc_pair_prof(i, j,    psq, p->e_pair2);
  emitsc_stck1  = emitsc_stck_prof(i, j, L, psq, p->e_pair1, p->e_stck1);
  emitsc_stck2  = emitsc_stck_prof(i, j, L, psq, p->e_pair2, p->e_stck2);
  
   // Follow the grammar
  switch(w) {
  case RBG_S:
  case RBG_P:
  case RBG_M:
  case RBG_M1:
  case RBG_R:
    break;     

  case RBG_F0:
    if (d > 2) {
      if (force_bp || allow_bp) {
	/* rule3: F0 -> a F5 a' */
	thispp = omx->F0->dp[j][d] + imx->F5->dp[j-1][d-2] + p->tF0[0] + emitsc_pair1;
	pp     = e2_FLogsum(pp, thispp);

	/* rule4: F0 -> a P a' */
	thispp = omx->F0->dp[j][d] + imx->P->dp[j-1][d-2] + p->tF0[1] + emitsc_pair2;
	pp     = e2_FLogsum(pp, thispp);
      }
    }
        
    /* rule5: F0 -> a a' */
    if (d == 2) {
      if (force_bp || allow_bp) {
	thispp = omx->F0->dp[j][d] + p->tF0[2] + emitsc_pair2;
	pp = e2_FLogsum(pp, thispp);
      }
    }
    break;
    
  case RBG_F5:
    if (d > 2) {
      if (force_bp || allow_bp) {
	/* rule6: F5 -> a F5 a' */
	thispp = omx->F5->dp[j][d] + imx->F5->dp[j-1][d-2] + p->tF5[0] + emitsc_stck1;
	pp = e2_FLogsum(pp, thispp);

	/* rule7: F5 -> a P a' */
	thispp = omx->F5->dp[j][d] + imx->P->dp[j-1][d-2] + p->tF5[1] + emitsc_stck2;
	pp = e2_FLogsum(pp, thispp);
      }
    }
        
    /* rule8: F5 -> a a' */
    if (d == 2) {
      if (force_bp || allow_bp) {
	thispp =
	  omx->F5->dp[j][d] + p->tF5[2] + emitsc_stck2;
        pp = e2_FLogsum(pp, thispp);
      }
    }
    
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "cannot recognize RBG nt %d\n", w);
  }

  post->pp[i][j] = pp;
  post->pp[j][i] = pp;										
  
  return eslOK;

 ERROR:
  return status;
}

// do not allow if pair is in the exclude list
static int
allow_bpair(double power_thresh, double neg_eval_thresh, int hloop_min, int i, int j, int L, int *covct, COVLIST *exclude, SPAIR *spair, int *ret_nneg) 
{
  double power;
  double eval;
  int    allow = FALSE;
  int    idx;
  int    n;

  if (i == 0 || j == 0 || i == j) return allow;
  
  // check if is has the minimum loop requeriment
  // check if pair is in the excluded list
  for (n = 0; n < exclude->n; n ++) {
    if ((exclude->cov[n].i == i && exclude->cov[n].j == j) ||
	(exclude->cov[n].i == j && exclude->cov[n].j == i)   ) return FALSE;
  }

  // check that the pair can form and does not have too much power or an Evalue in the grey zone
  if (covct[i] == 0 && covct[j] == 0) {
    idx   = INDEX(i-1, j-1, L);
    power = spair[idx].power;
    eval  = spair[idx].Eval;

    if      (power < power_thresh)    allow = TRUE;   // no power,  allow to pair
    else if (eval  < neg_eval_thresh) allow = TRUE;   // has power, but it is in the gray zone for significantly covarying, allow to pair
    else {                                            // power and unlikely to be covarying, a negative pair
      //printf("^^ negative: I %d J %d\n", i, j);
      if (ret_nneg) (*ret_nneg) ++;
    }
   }

  return allow;
}

// a basepair is forced only if a covarying pair, eg, covct[i] = j
static int
force_bpair(int i, int j, int *covct) 
{
  int force = FALSE;

  if (covct[i] == j && covct[j] == i) force = TRUE; // already paired to each other
  
  return force;
}


// a set of residues i..j is allowed to form a hairpin loop if
//
//    (1) there are no cov residues inside
//    (2) if the closing pair is a covaring pair, allow any length
//    (3) if the closing pair is not covaring enforce the hloop_min cutoff
//
static int
allow_hairpin(int hloop_min, int i, int j, int L, int *covct)
{
  int allow = TRUE;
  int iscov;         // TRUE if closing basepair is a covariation
  int hlen;          // hairpin len
  int ibp, jbp;
  int k;

  hlen = j - i + 1;
  ibp  = (i > 0)? i - 1 : 0;
  jbp  = (j < L)? j + 1 : 0;
  
  // first check that there is no covariation inside
  for (k = i; k <= j; k ++) if (covct[k] > 0) return FALSE;
  
  // if closing pair is not a covarying pair, force the hloop_min limit
  iscov = force_bpair(ibp, jbp, covct);
  if (!iscov && hlen < hloop_min) return FALSE;

  return allow;
}

// a set of residues i..j is allowed to form a loop, unless any of the residues
// is involved in a covarying basepair
static int
allow_loop(int i, int j, int *covct)
{
  int allow = TRUE;
  int k;

  // first check that there is no covariation inside
  for (k = i; k <= j; k ++) if (covct[k] > 0) return FALSE;

  return allow;
}

// A residue is allowed to be single stranded unless it is involved in  a covariation
// covct[i] > 0 if in a covarying pair
//
static int
allow_single(int i, int *covct) 
{
  int allow = TRUE;

  if (covct[i] > 0) allow = FALSE; // in a covarying pair

  return allow;
}

static SCVAL
emitsc_stck(int i, int j, int L, ESL_DSQ *dsq, SCVAL e_pair[NP], SCVAL e_stck[NP][NP])
{
  SCVAL sc;
  int   idx;
  int   cdx;
  int   ip = i-1;
  int   jp = j+1;

  /* no stacking on gaps of any kind or a basepair involving the first or/and last position */
  if (dsq[ip] >= NB || dsq[jp] >= NB || i == 1 || j == L) { 
    return emitsc_pair(i, j, dsq, e_pair); 
  }
 
  cdx = dsq[ip]*NB + dsq[jp];

  if (dsq[i] >= NB || dsq[j] >= NB) { // ignore gaps
    sc = -eslINFINITY;
  }
  else {
    idx = dsq[i]*NB + dsq[j];
    sc = e_stck[cdx][idx];
  }
  return sc;
}
static SCVAL
emitsc_stck_prof(int i, int j, int L, PSQ *psq, SCVAL e_pair[NP], SCVAL e_stck[NP][NP])
{
  float *profi  = psq->prof[i];
  float *profj  = psq->prof[j];
  float *profip;
  float *profjp;
  float  factor = 3.0;
  float  pi, pj;
  float  pip, pjp;
  float num     = -eslINFINITY;
  float den     = -eslINFINITY;
  float sc      = -eslINFINITY;
  float logval;
  int   idx;
  int   cdx;
  int   ip = i-1;
  int   jp = j+1;
  int   k1, k2, k3, k4;

  /* no stacking on gaps of any kind or a basepair involving the first or/and last position */
  if (i <= 1 || j >= L) 
    return emitsc_pair_prof(i, j, psq, e_pair); 

  profip = psq->prof[ip];
  profjp = psq->prof[jp];

  e2_FLogsumInit();
    
  for (k1 = 0; k1 < NB; k1 ++) 
    for (k2 = 0; k2 < NB; k2 ++) {
      cdx = k1*NB + k2;
      
      pip = profip[k1];
      pjp = profjp[k2];
      
      if (fabs(pip-pjp) < log(factor)) {
	logval = pip + pjp + (float)e_pair[cdx];
	den    = e2_FLogsum(den, logval);
      
	for (k3 = 0; k3 < NB; k3 ++) 
	  for (k4 = 0; k4 < NB; k4 ++) {
	    
	    idx = k3*NB + k4;
	    
	    pi = profi[k3];
	    pj = profj[k4];
	    
	    if (fabs(pi-pj) < log(factor)) {
	      logval = pi + pj + pip + pjp + (float)e_stck[cdx][idx] + (float)e_pair[cdx];
	      num    = e2_FLogsum(num, logval);
	    }
	  }
      }
    }
  if (den > -eslINFINITY) sc = num - den;
  
  return (SCVAL)sc;
}

static SCVAL
emitsc_pair(int i, int j, ESL_DSQ *dsq, SCVAL e_pair[NP])
{
  SCVAL sc;
  int   idx;

  if (dsq[i] >= NB || dsq[j] >= NB) { // ignore gaps
    sc = -eslINFINITY;
  }
  else {
    idx = dsq[i]*NB + dsq[j];
    sc = e_pair[idx];
  }

  return sc;
}
static SCVAL
emitsc_pair_prof(int i, int j, PSQ *psq, SCVAL e_pair[NP])
{
  float *profi;
  float *profj;
  float  factor = 3.0;
  float  pi, pj;
  float  sc    = -eslINFINITY;
  float  logval;
  int    idx;
  int    k1, k2;

  if (i >= j) return (SCVAL)sc;
  
  e2_FLogsumInit();
  
  profi = psq->prof[i];
  profj = psq->prof[j];
  for (k1 = 0; k1 < NB; k1 ++) 
    for (k2 = 0; k2 < NB; k2 ++) {

      idx = k1*NB + k2;
      pi = profi[k1];
      pj = profj[k2];

      if (fabs(pi-pj) < log(factor)) {
	logval = pi + pj + (float)e_pair[idx];
	sc     = e2_FLogsum(sc, logval);
      }
  }

  return (SCVAL)sc;
}

static SCVAL
emitsc_sing(int i, ESL_DSQ *dsq, SCVAL e_sing[NB])
{
  SCVAL sc;
  
  if (dsq[i] < NB) sc = e_sing[dsq[i]];
  else             sc = log(0.25);

  return sc;
}
static SCVAL
emitsc_sing_prof(int i, PSQ *psq, SCVAL e_sing[NB])
{
  float *prof;
  float  sc = -eslINFINITY;
  float  logval;
  int    k;

  if (i > psq->n) return (SCVAL)sc;

  e2_FLogsumInit();
  
  prof = psq->prof[i];
  for (k = 0; k < NB; k ++) {
    logval = prof[k] + (float)e_sing[k];
    sc = e2_FLogsum(sc, logval);
  }
  // columns with all gaps
  if (sc == -eslINFINITY) sc = log(0.25);
  
  return (SCVAL)sc;
}

static SCVAL
score_loop_hairpin(int i, int j, RBGparam *p, ESL_DSQ *dsq)
{
  SCVAL sc = 0.;
  int   k;

  for (k = i; k <= j; k ++) 
    sc += emitsc_sing(k, dsq, p->e_sing_l1);

  return sc;
}
static SCVAL
score_loop_hairpin_prof(int i, int j, RBGparam *p, PSQ *psq)
{
  SCVAL sc = 0.;
  int   k;

  for (k = i; k <= j; k ++) 
    sc += emitsc_sing_prof(k, psq, p->e_sing_l1);

  return sc;
}

static SCVAL
score_loop_bulge(int i, int j, RBGparam *p, ESL_DSQ *dsq)
{
  SCVAL sc = 0.;
  int   k;

  for (k = i; k <= j; k ++) 
    sc += emitsc_sing(k, dsq, p->e_sing_l2);

  return sc;
}
static SCVAL
score_loop_bulge_prof(int i, int j, RBGparam *p, PSQ *psq)
{
  SCVAL sc = 0.;
  int   k;

  for (k = i; k <= j; k ++) 
    sc += emitsc_sing_prof(k, psq, p->e_sing_l2);
  
  return sc;
}

static SCVAL
score_loop_intloop(int i, int j, RBGparam *p, ESL_DSQ *dsq)
{
  SCVAL sc = 0.;
  int   k;

  for (k = i; k <= j; k ++) 
    sc += emitsc_sing(k, dsq, p->e_sing_l3);

  return sc;
}
static SCVAL
score_loop_intloop_prof(int i, int j, RBGparam *p, PSQ *psq)
{
  SCVAL sc = 0.;
  int   k;

  for (k = i; k <= j; k ++) 
    sc += emitsc_sing_prof(k, psq, p->e_sing_l3);

  return sc;
}

static int
segment_remove_gaps(int i, int j, ESL_DSQ *dsq)
{
  int newlen = 0;
  int x;

  for (x = i; x <= j; x ++) 
    if (dsq[x] < NB) newlen ++;

  return newlen;
}

static int
segment_remove_gaps_prof(int i, int j, PSQ *psq)
{
  int newlen = 0;
  int x;

  for (x = i; x <= j; x ++) 
    if (psq->prof[x][NB] < 1.) newlen ++;

  return newlen;
}

/* Function:  esl_vec_DLogNorm()
 */
static int
vec_SCVAL_LogNorm(SCVAL *scvec, int n)
{
  double *vec = NULL;
  double  denom;
  int     i;
  int     status;
  
  ESL_ALLOC(vec, sizeof(double *) * n);
  for (i = 0; i < n; i ++) vec[i] = (double)scvec[i];
  
  denom = esl_vec_DLogSum(vec, n);
  esl_vec_DIncrement(vec, n, -1.*denom);

  for (i = 0; i < n; i ++) scvec[i] = (SCVAL)vec[i];

  free(vec);
  return eslOK;
  
 ERROR:
  if (vec) free(vec);
  return status;
}

static int
dvec_SCVAL_LogNorm(int n1, int n2, SCVAL dvec[n1][n2])
{
  double *vec = NULL;
  double  denom;
  int     n = n1 * n2;
  int     i1, i2;
  int     status;

  ESL_ALLOC(vec, sizeof(double) * n);
  for (i1 = 0; i1 < n1; i1 ++)
    for (i2 = 0; i2 < n2; i2 ++)
      vec[i1*n2 + i2] = (double)dvec[i1][i2];
  
  denom = esl_vec_DLogSum(vec, n);
  esl_vec_DIncrement(vec, n, -1.*denom);
  
  for (i1 = 0; i1 < n1; i1 ++)
    for (i2 = 0; i2 < n2; i2 ++)
      dvec[i1][i2] = (SCVAL)vec[i1*n2+i2];

  free(vec);
  return eslOK;

 ERROR:
  if (vec) free(vec);
  return status;
			  
}
