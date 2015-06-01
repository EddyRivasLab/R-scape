/* cococyk.c */

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

#include "covgrammars.h"
#include "cococyk.h"


static int   dp_recursion_g6 (G6param  *p, ESL_SQ *sq, int *ct, G6_MX  *cyk, int j, int d, SCVAL *ret_sc, ESL_STACK *alts, char *errbuf, int verbose);
static int   dp_recursion_g6s(G6Sparam *p, ESL_SQ *sq, int *ct, G6_MX  *cyk, int j, int d, SCVAL *ret_sc, ESL_STACK *alts, char *errbuf, int verbose);
static int   dp_recursion_bgr(BGRparam *p, ESL_SQ *sq, int *ct, BGR_MX *cyk, int j, int d, SCVAL *ret_sc, ESL_STACK *alts, char *errbuf, int verbose);
static int   allow_single(int i, int *ct) ;
static int   allow_bpair(int i, int j, int *ct);
static SCVAL emitsc_stck(int i, int j, ESL_DSQ *dsq, SCVAL *e_pair, SCVAL **e_stck);
static SCVAL emitsc_pair(int i, int j, ESL_DSQ *dsq, SCVAL *e_pair);
static SCVAL emitsc_sing(int i,        ESL_DSQ *dsq, SCVAL *e_sing);

int
COCOCYK(ESL_RANDOMNESS *r, enum grammar_e G, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6param  *g6p  = NULL;
  G6Sparam *g6sp = NULL;
  BGRparam *bgrp = NULL;
  int status;

  /* get the grammar parameters and run the corresponding CYK */
  switch(G) {
  case G6:
    /* Transfer scores from static built-in storage */
    status = COCOVYK_G6_GetParam(&g6p, errbuf, verbose);
    if (status != eslOK) goto ERROR;    
    status = COCOCYK_G6(r, g6p, sq, ct, ret_cct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case G6S:
    status = COCOVYK_G6S_GetParam(&g6sp, errbuf, verbose);
    if (status != eslOK) goto ERROR; 
    status = COCOCYK_G6S(r, g6sp, sq, ct, ret_cct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case BGR:
    status = COCOVYK_BGR_GetParam(&bgrp, errbuf, verbose);
    if (status != eslOK) goto ERROR; 
    status = COCOCYK_BGR(r, bgrp, sq, ct, ret_cct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  }

  if (g6p)  free(g6p);
  if (g6sp) free(g6sp);
  if (bgrp) free(bgrp);
  return eslOK;

 ERROR:
  return status;
}

int
COCOVYK_G6_GetParam(G6param **ret_p, char *errbuf, int verbose)
{
  G6param *p = NULL;
  int      x;
  int      status;

 ESL_ALLOC(p, sizeof(G6param));

  p->t1[0] = G6_PRELOADS_TrATrBTrB.t1[0];
  p->t1[1] = G6_PRELOADS_TrATrBTrB.t1[1];
  p->t2[0] = G6_PRELOADS_TrATrBTrB.t2[0];
  p->t2[1] = G6_PRELOADS_TrATrBTrB.t2[1];
  p->t3[0] = G6_PRELOADS_TrATrBTrB.t3[0];
  p->t3[1] = G6_PRELOADS_TrATrBTrB.t3[1];

  for (x = 0; x < 4;  x ++) p->e_sing[x] = G6_PRELOADS_TrATrBTrB.e_sing[x];
  for (x = 0; x < 16; x ++) p->e_pair[x] = G6_PRELOADS_TrATrBTrB.e_pair[x];

  *ret_p = p;
  return eslOK;

 ERROR:
  if (p) free(p);
  return status;
}

int
COCOVYK_G6S_GetParam(G6Sparam **ret_p, char *errbuf, int verbose)
{
 G6Sparam *p = NULL;
 int       x, y;
 int       status;

 ESL_ALLOC(p, sizeof(G6Sparam));

  p->t1[0] = G6S_PRELOADS_TrATrBTrB.t1[0];
  p->t1[1] = G6S_PRELOADS_TrATrBTrB.t1[1];
  p->t2[0] = G6S_PRELOADS_TrATrBTrB.t2[0];
  p->t2[1] = G6S_PRELOADS_TrATrBTrB.t2[1];
  p->t3[0] = G6S_PRELOADS_TrATrBTrB.t3[0];
  p->t3[1] = G6S_PRELOADS_TrATrBTrB.t3[1];

  for (x = 0; x < 4;  x ++)   p->e_sing[x]    = G6S_PRELOADS_TrATrBTrB.e_sing[x];
  for (x = 0; x < 16; x ++)   p->e_pair[x]    = G6S_PRELOADS_TrATrBTrB.e_pair[x];
  for (x = 0; x < 16; x ++) 
    for (y = 0; y < 16; y ++) p->e_stck[x][y] = G6S_PRELOADS_TrATrBTrB.e_stck[x][y];

  *ret_p = p;
  return eslOK;

 ERROR:
  if (p) free(p);
  return status;

  return eslOK;
}

int
COCOVYK_BGR_GetParam(BGRparam **ret_p, char *errbuf, int verbose)
{
 BGRparam *p = NULL;
 int       x, y;
 int       l;
 int       status;

 ESL_ALLOC(p, sizeof(BGRparam));

  p->tS[0] = BGR_PRELOADS_TrATrBTrB.tS[0];
  p->tS[1] = BGR_PRELOADS_TrATrBTrB.tS[1];
  p->tS[2] = BGR_PRELOADS_TrATrBTrB.tS[2];

  p->tF0[0] = BGR_PRELOADS_TrATrBTrB.tF0[0];
  p->tF0[1] = BGR_PRELOADS_TrATrBTrB.tF0[1];

  p->tF5[0] = BGR_PRELOADS_TrATrBTrB.tF5[0];
  p->tF5[1] = BGR_PRELOADS_TrATrBTrB.tF5[1];

  p->tP[0] = BGR_PRELOADS_TrATrBTrB.tP[0];
  p->tP[1] = BGR_PRELOADS_TrATrBTrB.tP[1];
  p->tP[2] = BGR_PRELOADS_TrATrBTrB.tP[2];
  p->tP[3] = BGR_PRELOADS_TrATrBTrB.tP[3];
  p->tP[4] = BGR_PRELOADS_TrATrBTrB.tP[4];

  p->tM[0]  = BGR_PRELOADS_TrATrBTrB.tM[0];
  p->tM[1]  = BGR_PRELOADS_TrATrBTrB.tM[1];

  p->tR[0]  = BGR_PRELOADS_TrATrBTrB.tR[0];
  p->tR[1]  = BGR_PRELOADS_TrATrBTrB.tR[1];

  p->tM1[0] = BGR_PRELOADS_TrATrBTrB.tM1[0];
  p->tM1[1] = BGR_PRELOADS_TrATrBTrB.tM1[1];
  
  for (x = 0; x < 4;  x ++) {
    p->e_sing[x]    = BGR_PRELOADS_TrATrBTrB.e_sing[x];
    p->e_sing_l1[x] = BGR_PRELOADS_TrATrBTrB.e_sing_l1[x];
    p->e_sing_l2[x] = BGR_PRELOADS_TrATrBTrB.e_sing_l2[x];
    p->e_sing_l3[x] = BGR_PRELOADS_TrATrBTrB.e_sing_l3[x];
  }
  
  for (x = 0; x < 16; x ++) {
    p->e_pair1[x] = BGR_PRELOADS_TrATrBTrB.e_pair1[x];
    p->e_pair2[x] = BGR_PRELOADS_TrATrBTrB.e_pair2[x];
    
    for (y = 0; y < 16; y ++) {
      p->e_stck1[x][y] = BGR_PRELOADS_TrATrBTrB.e_stck1[x][y];
      p->e_stck2[x][y] = BGR_PRELOADS_TrATrBTrB.e_stck2[x][y];
    }
  }

  for (l = 0; l < MAXLOOP_H; l ++) p->l1[l] = BGR_PRELOADS_TrATrBTrB.l1[l];
  for (l = 0; l < MAXLOOP_B; l ++) p->l2[l] = BGR_PRELOADS_TrATrBTrB.l2[l];
  for (l = 0; l < MAXLOOP_I; l ++) p->l3[l] = BGR_PRELOADS_TrATrBTrB.l3[l];

  *ret_p = p;
  return eslOK;

 ERROR:
  if (p) free(p);
  return status;

  return eslOK;
}

int
COCOCYK_G6(ESL_RANDOMNESS *r, G6param  *p, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6_MX *gmx = NULL;
  int    status;

  gmx = G6MX_Create(sq->n);

  /* Fill the cyk matrix */
  if ((status = COCOCYK_G6_Fill        (p, sq, ct, gmx, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = COCOCYK_G6_Traceback(r, p, sq, ct, gmx, ret_cct, errbuf, verbose))  != eslOK) goto ERROR;
  
  G6MX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) G6MX_Destroy(gmx);
  return status;
}

int
COCOCYK_G6S(ESL_RANDOMNESS *r, G6Sparam  *p, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6_MX *gmx = NULL;
  int    status;

  gmx = G6MX_Create(sq->n);

  /* Fill the cyk matrix */
  if ((status = COCOCYK_G6S_Fill        (p, sq, ct, gmx, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = COCOCYK_G6S_Traceback(r, p, sq, ct, gmx, ret_cct, errbuf, verbose))  != eslOK) goto ERROR;
  
  G6MX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) G6MX_Destroy(gmx);
  return status;
}

int
COCOCYK_BGR(ESL_RANDOMNESS *r, BGRparam  *p, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  BGR_MX *gmx = NULL;
  int     status;

  gmx = BGRMX_Create(sq->n);

  /* Fill the cyk matrix */
  if ((status = COCOCYK_BGR_Fill        (p, sq, ct, gmx, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = COCOCYK_BGR_Traceback(r, p, sq, ct, gmx, ret_cct, errbuf, verbose))  != eslOK) goto ERROR;
  
  BGRMX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) BGRMX_Destroy(gmx);
  return status;

}

int
COCOCYK_G6_Fill(G6param *p, ESL_SQ *sq, int *ct, G6_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  SCVAL sc = -eslINFINITY;
  int   L = sq->n;
  int   j, d;
  int   status;

 /* G6 grammar
  */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {
	status = dp_recursion_g6(p, sq, ct, cyk, j, d, &(cyk->S->dp[j][d]), NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6 cocoCYK failed");
	if (verbose) printf("\nG6 cocoCYK %f j=%d d=%d L=%d\n", cyk->S->dp[j][d], j, d, L); 
     } 
  sc = cyk->S->dp[L][L];
  if (verbose) printf("G6 cocoCYKscore = %f\n", sc);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
COCOCYK_G6S_Fill(G6Sparam  *p, ESL_SQ *sq, int *ct, G6_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  SCVAL sc = -eslINFINITY;
  int   L = sq->n;
  int   j, d;
  int   status;

 /* G6S grammar
  */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {
	status = dp_recursion_g6s(p, sq, ct, cyk, j, d, &(cyk->S->dp[j][d]), NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "G6 cocoCYK failed");
	if (verbose) printf("\nG6S cocoCYK %f j=%d d=%d L=%d\n", cyk->S->dp[j][d], j, d, L); 
     } 
  sc = cyk->S->dp[L][L];
  if (verbose) printf("G6S cocoCYKscore = %f\n", sc);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
COCOCYK_BGR_Fill(BGRparam  *p, ESL_SQ *sq, int *ct, BGR_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  SCVAL sc = -eslINFINITY;
  int   L = sq->n;
  int   j, d;
  int   status;

 /* BGR grammar
  */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {
	status = dp_recursion_bgr(p, sq, ct, cyk, j, d, &(cyk->S->dp[j][d]), NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "BGR cocoCYK failed");
	if (verbose) printf("\nBGR cocoCYK %f j=%d d=%d L=%d\n", cyk->S->dp[j][d], j, d, L); 
     } 
  sc = cyk->S->dp[L][L];
  if (verbose) printf("BGR cocoCYKscore = %f\n", sc);

  *ret_sc = sc;
  return eslOK;

 ERROR:
  return status;
}

int
COCOCYK_G6_Traceback(ESL_RANDOMNESS *r, G6param  *p, ESL_SQ *sq, int *ct, G6_MX *cyk, int **ret_cct, char *errbuf, int verbose) 
{
  int *cct = NULL;
  int  L = sq->n;
  int  status;

  ESL_ALLOC(cct, sizeof(int) * (L+1));
  esl_vec_ISet(cct, L+1, 0);
       
  *ret_cct = cct;
  return eslOK;

 ERROR:
  if (cct) free(cct);
  return status;
}

int
COCOCYK_G6S_Traceback(ESL_RANDOMNESS *r, G6Sparam  *p, ESL_SQ *sq, int *ct, G6_MX *cyk, int **ret_cct, char *errbuf, int verbose) 
{
   int *cct = NULL;
  int  L = sq->n;
  int  status;

  ESL_ALLOC(cct, sizeof(int) * (L+1));
  esl_vec_ISet(cct, L+1, 0);
       
  *ret_cct = cct;
  return eslOK;

 ERROR:
  if (cct) free(cct);
  return status;
}

int
COCOCYK_BGR_Traceback(ESL_RANDOMNESS *r, BGRparam  *p, ESL_SQ *sq, int *ct, BGR_MX *cyk, int **ret_cct, char *errbuf, int verbose) 
{
  int *cct = NULL;
  int  L = sq->n;
  int  status;

  ESL_ALLOC(cct, sizeof(int) * (L+1));
  esl_vec_ISet(cct, L+1, 0);
       
  *ret_cct = cct;
  return eslOK;

 ERROR:
  if (cct) free(cct);
  return status;
}


static int 
dp_recursion_g6 (G6param  *p, ESL_SQ *sq, int *ct, G6_MX  *cyk, int j, int d, SCVAL *ret_sc, ESL_STACK *alts, char *errbuf, int verbose)
{
  ESL_DSQ *dsq = sq->dsq;
  SCVAL    bestsc = -eslINFINITY;
  SCVAL    sc;
  int      d1;
  int      r;
  int      i, k;

  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;

  if (d == 0) { *ret_sc = bestsc; return eslOK; }

  /* rule0: S -> LS */
  for (d1 = 0; d1 <= d; d1++) {
    k = i + d1 - 1;
    
    sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1+1] + p->t1[0];

    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, r);
	esl_stack_IPush(alts, d1);
      }
    }
  }
  r++;

  /* rule1: S -> L */
  d1 = 0;
  sc = cyk->L->dp[j][d] + p->t1[1];
  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;
  
  /* rule2: L -> a F a' */
  d1 = 0;
  sc = (allow_bpair(i+1, j+1, ct))?  cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pair(i+1, j+1, dsq, &p->e_pair[0]) : -eslINFINITY;
  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;

  /* rule3: L -> a */
  d1 = 0;
  sc = (allow_single(i+1, ct))? p->t2[1] + emitsc_sing(i+1, dsq, p->e_sing) : -eslINFINITY;
  
  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;

  /* rule4: F -> a F a' */
  d1 = 0;
  sc = (allow_bpair(i+1, j+1, ct))?  cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_pair(i+1, j+1, dsq, &p->e_pair[0]) : -eslINFINITY;

  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;

  /* rule5: F -> LS */
  for (d1 = 0; d1 <= d; d1++) {
    k = i + d1 - 1;
    
    sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1+1] + p->t3[1];

    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, r);
	esl_stack_IPush(alts, d1);
      }
    }
  }
  r++;

  *ret_sc = bestsc;
 
  return eslOK;
}

static int 
dp_recursion_g6s(G6Sparam *p, ESL_SQ *sq, int *ct, G6_MX  *cyk, int j, int d, SCVAL *ret_sc, ESL_STACK *alts, char *errbuf, int verbose)
{
  ESL_DSQ *dsq = sq->dsq;
  SCVAL    bestsc = -eslINFINITY;
  SCVAL    sc;
  int      d1;
  int      r;
  int      i, k;
  
  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;

  if (d == 0) { *ret_sc = bestsc; return eslOK; }

  /* rule0: S -> LS */
  for (d1 = 0; d1 <= d; d1++) {
    k = i + d1 - 1;
    
    sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1+1] + p->t1[0];

    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, r);
	esl_stack_IPush(alts, d1);
      }
    }
  }
  r++;

  /* rule1: S -> L */
  d1 = 0;
  sc = cyk->L->dp[j][d] + p->t1[1];
  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;
  
  /* rule2: L -> a F a' */
  d1 = 0;
  sc = (allow_bpair(i+1, j+1, ct))?  cyk->F->dp[j-1][d-2] + p->t2[0] + emitsc_pair(i+1, j+1, dsq, &p->e_pair[0]) : -eslINFINITY;
  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;

  /* rule3: L -> a */
  d1 = 0;
  sc = (allow_single(i+1, ct))? p->t2[1] + emitsc_sing(i+1, dsq, p->e_sing) : -eslINFINITY;
  
  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;

  /* rule4: F -> a F a' */
  d1 = 0;
  sc = (allow_bpair(i+1, j+1, ct))?  cyk->F->dp[j-1][d-2] + p->t3[0] + emitsc_stck(i+1, j+1, dsq, &p->e_pair[0], (SCVAL **)p->e_stck) : -eslINFINITY;

  if (sc >= bestsc) {
    if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
      if (alts) esl_stack_Reuse(alts);
      bestsc = sc;
    }     
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }
  r++;

  /* rule5: F -> LS */
  for (d1 = 0; d1 <= d; d1++) {
    k = i + d1 - 1;
    
    sc = cyk->L->dp[k][d1] + cyk->S->dp[j][d-d1+1] + p->t3[1];

    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }     
      if (alts) {
	esl_stack_IPush(alts, r);
	esl_stack_IPush(alts, d1);
      }
    }
  }
  r++;

  *ret_sc = bestsc;

  return eslOK;
}

static int 
dp_recursion_bgr(BGRparam *p, ESL_SQ *sq, int *ct, BGR_MX *cyk, int j, int d, SCVAL *ret_sc, ESL_STACK *alts, char *errbuf, int verbose)
{
  ESL_DSQ *dsq = sq->dsq;
  SCVAL    bestsc = -eslINFINITY;
  SCVAL    sc;
  int      d1;
  int      r;
  int      i, k;
  
  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;

  /* rule0: S -> a S */
  r++;

  /* rule1: S -> L */
  r++;

  /* rule2: S -> epsilon */
  r++;

  /* rule3: F0 -> a F5 a' */
  r++;

  /* rule4: F0 -> a P a' */
  r++;

  /* rule5: F5 -> a F5^{bb'} a' */
  r++;

  /* rule6: F5 -> a P^{bb'} a' */
  r++;

  /* rule7: P -> m..m */
  r++;

  /* rule8: P -> m..m F0 */
  r++;

  /* rule9: P -> F0 m..m */
  r++;

  /* rule10: P -> m..m F0 m..m */
  r++;

  /* rule11: P -> M1 M */
  r++;

  /* rule12: M -> M M1 */
  r++;

  /* rule13: M -> R */
  r++;

  /* rule14: R -> R a */
  r++;

  /* rule15: R -> M1 */
  r++;

  /* rule16: M1 -> a M1 */
  r++;

  /* rule17: M1 -> F0 */
  r++;


 
 
  *ret_sc = bestsc;

  return eslOK;
}


static int
allow_single(int i, int *ct) 
{
  int allow = FALSE;

  if (ct[i] == 0 ) allow = TRUE; // not paired

  return allow;
}

static int
allow_bpair(int i, int j, int *ct) 
{
  int allow = FALSE;

  if (ct[i] == j && ct[j] == i) allow = TRUE; // already paired to each other
  if (ct[i] == 0 && ct[j] == 0) allow = TRUE; // both unpaired

  return allow;
}



static SCVAL
emitsc_stck(int i, int j, ESL_DSQ *dsq, SCVAL *e_pair, SCVAL **e_stck)
{
  SCVAL sc;
  int   idx;
  int   cdx;
  int   x;
  int   ip = i-1;
  int   jp = j+1;

  /* no stcking on gaps of any kind */
  if (dsq[ip] >= 4 || dsq[jp] >= 4) { 
    return emitsc_pair(i, j, dsq, e_pair); 
  }

  cdx = dsq[ip]*4 + dsq[jp];

  if (dsq[i] >= 4 && dsq[j] >= 4) { // ignore double gaps
    sc = -eslINFINITY;
  }
  else if (dsq[i] >= 4) { // single gaps treat as missing data
    sc = 0.0;
    for (x = 0; x < 4; x++) {
      idx = x*4 + dsq[j];
      sc += e_stck[cdx][idx];
    }
  }
  else if (dsq[j] >= 4) {  // single gaps treat as missing data
    sc = 0.0;
    for (x = 0; x < 4; x++) {
      idx = dsq[i]*4 + x;
      sc += e_stck[cdx][idx];
    }
  }
  else {
    idx = dsq[i]*4  + dsq[j];
    printf("cdx %d idx %d\n", cdx, idx);
    sc = e_stck[cdx][idx];
  }

  return sc;
}

static SCVAL
emitsc_pair(int i, int j, ESL_DSQ *dsq, SCVAL *e_pair)
{
  SCVAL sc;
  int   idx;
  int   x;

  if (dsq[i] >= 4 && dsq[j] >= 4) { // ignore double gaps
    sc = -eslINFINITY;
  }
  else if (dsq[i] >= 4) { // single gaps treat as missing data
    sc = 0.0;
    for (x = 0; x < 4; x++) {
      idx = x*4 + dsq[j];
      sc += e_pair[idx];
    }
  }
  else if (dsq[j] >= 4) {  // single gaps treat as missing data
    sc = 0.0;
    for (x = 0; x < 4; x++) {
      idx = dsq[i]*4 + x;
      sc += e_pair[idx];
    }
  }
  else {
    idx = dsq[i] *4 + dsq[j];
    sc = e_pair[idx];
  }

  return sc;
}

static SCVAL
emitsc_sing(int i, ESL_DSQ *dsq, SCVAL *e_sing)
{
  SCVAL sc;
  
  if (dsq[i] < 4) sc = e_sing[dsq[i]];
  else            sc = 0.0;

  return sc;
}
