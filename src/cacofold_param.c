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
#include "structure.h"

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
CACO_RBG_R3D_GetParam(R3D *r3d, RBGparam **ret_rbgp, R3Dparam **ret_r3dp, char *errbuf, int verbose)
{
  RBGparam *rbgp = NULL;
  R3Dparam *r3dp = NULL;
  SCVAL     tP0;
  SCVAL     tP1;
  SCVAL     tP2;
  SCVAL     tP3;
  int       status;
  
  status = CACO_RBG_GetParam(&rbgp, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  
  status = R3D_GetParam(&r3dp, errbuf, verbose);
  if (status != eslOK) goto ERROR;

  // modify the  P -> t[0] m...m | t[1] m..., F0 | t[2] F0 m...m | t[3] m...m F0 m...m
  //  to
  //             P -> t[0]*(1-pHL)  HL_0  | t[0]*pHL/nHL  HL_1  | ... | t[0]*pHL/nHL  HL_{nHL}  
  //             P -> t[1]*(1-pBL5) BL_0  | t[1]*pBL5/nBL BL_1  | ... | t[1]*pBL5/nBL BL_{nBL}  
  //             P -> t[2]*(1-pBL3) BL_0  | t[2]*pBL3/nBL BL_1  | ... | t[2]*pBL3/nBL BL_{nBL}  
  //             P -> t[3]*(1-pIL)  IL_0  | t[3]*pHL/nIL  IL_1  | ... | t[3]*pHL/nIL  IL_{nIL}  

  if (r3d->nHL > 0) {
    tP0 = rbgp->tP[0];
    rbgp->tP[0]     = tP0 + log(1.0 - exp(r3dp->HLp->pHL));
    r3dp->HLp->pHL += tP0 - log(r3d->nHL);
  }
  if (r3d->nBL > 0) {
    tP1 = rbgp->tP[1];
    tP2 = rbgp->tP[2];
    rbgp->tP[1]      = tP1 + log(1.0 - exp(r3dp->BLp->pBL5));
    rbgp->tP[2]      = tP2 + log(1.0 - exp(r3dp->BLp->pBL3));
    r3dp->BLp->pBL5 += tP1 - log(r3d->nBL);
    r3dp->BLp->pBL3 += tP2 - log(r3d->nBL);
  }
  if (r3d->nIL_total > 0) {
    tP3 = rbgp->tP[3];
    rbgp->tP[3]     = tP3 + log(1.0 - exp(r3dp->ILp->pIL));
    r3dp->ILp->pIL += tP3 - log(r3d->nIL_total);
  }

  if (verbose) {
    printf("RBG_R3D Param\n");
    printf("P0 %f to %f rest %f\n", tP0, rbgp->tP[0], r3dp->HLp->pHL);
    printf("P1 %f to %f rest %f\n", tP1, rbgp->tP[1], r3dp->BLp->pBL5);
    printf("P2 %f to %f rest %f\n", tP2, rbgp->tP[2], r3dp->BLp->pBL3);
    printf("P3 %f to %f rest %f\n", tP3, rbgp->tP[3], r3dp->ILp->pIL);
  }
  
  *ret_rbgp = rbgp;
  *ret_r3dp = r3dp;
  return status;

 ERROR:
  if (rbgp) free(rbgp);
  if (r3dp) R3D_Param_Destroy(r3dp);

  return status;
}
