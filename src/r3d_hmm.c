/* r3d_hmm.c */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_buffer.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "e2_profilesq.h"
#include "logsum.h"
#include "r3d.h"
#include "r3d_hmm.h"

static SCVAL emitsc_prof(int iabs, int K, double **pm, SCVAL *e);

// For a sequence of length M
//
//   S^{0} --> S^{1} | I^{0}
//
//   S^{k} --> x_k S^{k+1} | x_k I^{k} | S^{k+1} | I^{k} ,, k in [1..M-1]
//   S^{M} --> x_M         | x_M I^{M} | e       | I^{M} 
//
//   I^{k} --> x I^{k} | S^{k+1} ,, k in [0..M-1]
//   I^{n} --> x I^{n} | e      
//
// if len(RNAmotif) = 0
//
//   S^{0} -->  x I^{0} | e 
//
//
extern R3D_HMM *
R3D_hmm_Create(const ESL_ALPHABET *abc, char *RMmotif, char *errbuf, int verbose)
{
  R3D_HMM *hmm = NULL;
  float   *p;
  char     c;
  int      M;
  int      k;
  int      status;

  ESL_ALLOC(hmm, sizeof(R3D_HMM));

  // Transition probabilities (all use the same)
  //
  //   S^{0} --> S^{1} | I^{0}
  //             tB[0]   tB[1]
  //
  hmm->tB[0] = log(HMM_tB);
  hmm->tB[1] = log(1. - HMM_tB);
  
  //   S^{k} --> x_k S^{k+1} | x_k I^{k} | S^{k+1} | I^{k}  ,, k in [1..M-1]
  //   S^{M} --> x_M         | x_M I^{M} | e       | I^{M}
  //                tS[0]       tS[1]      tS[2]     tS[3]
  //
  hmm->tS[0] = log(HMM_tM);
  hmm->tS[1] = log(HMM_tMI);
  hmm->tS[2] = log(HMM_tD);
  hmm->tS[3] = log(HMM_tDI);
  esl_vec_FLogNorm(hmm->tS, abc->K);

  //   I^{k} --> x I^{k} | S^{k+1} ,, k in [0..M-1]
  //   I^{n} --> x I^{n} | e      
  //              tI[0]     tI[1]
  //
  hmm->tI[0] = log(HMM_tI);
  hmm->tI[1] = log(1. - HMM_tI);
    
  // Emission probabilities (profiled)
  //
  M = strlen(RMmotif);
  hmm->e  = NULL;
  ESL_ALLOC(hmm->e,    sizeof(SCVAL *) * (M+1)); hmm->e[0] = NULL;
  ESL_ALLOC(hmm->e[0], sizeof(SCVAL)   * (M+1) * abc->K);
  for (k = 1; k <= M; k++)
    hmm->e[k] = hmm->e[0] + k*abc->K;
  
  // assign the emission probabilities for the I state
  // from G6X
  p = hmm->e[0];
  p[0] = -1.013848;
  p[1] = -1.753817;
  p[2] = -1.518912;
  p[3] = -1.406016;

  // assign the emission probabilities for each S state
  for (k = 0; k < M; k++) {
    c = RMmotif[k];
    p = hmm->e[k+1];
    
    if      (c == 'A')   {
      p[0] = 1 - 3*EPSILON;
      p[1] = p[2] = p[3] = EPSILON;
    }
    else if (c == 'C') {
      p[1] = 1 - 3*EPSILON;
      p[0] = p[2] = p[3] = EPSILON;
    }
    else if (c == 'G') {
      p[2] = 1 - 3*EPSILON;
      p[0] = p[1] = p[3] = EPSILON;
    }
    else if (c == 'U') {
      p[3] = 1 - 3*EPSILON;
      p[0] = p[1] = p[2] = EPSILON;
    }
    else if (c == 'T') {
      p[3] = 1 - EPSILON;
      p[0] = p[1] = p[2] = EPSILON;
    }
    else if (c == 'M') { // M = A/C
      p[0] = p[1] = (1 - 2*EPSILON) / 2;
      p[2] = p[3] = EPSILON;
    }
    else if (c == 'R') { // R = A/G
      p[0] = p[2] = (1 - 2*EPSILON) / 2;
      p[1] = p[3] = EPSILON;
    }
    else if (c == 'W') { // W = A/U
      p[0] = p[3] = (1 - 2*EPSILON) / 2;
      p[1] = p[2] = EPSILON;
    }
    else if (c == 'S') { // S = C/G
      p[1] = p[2] = (1 - 2*EPSILON) / 2;
      p[0] = p[3] = EPSILON;
    }
    else if (c == 'Y') { // Y = C/U
      p[1] = p[3] = (1 - 2*EPSILON) / 2;
      p[0] = p[2] = EPSILON;
    }
    else if (c == 'K') { // K = G/U
      p[2] = p[3] = (1 - 2*EPSILON) / 2;
      p[0] = p[1] = EPSILON;
    }
    else if (c == 'B') { // B = C/G/U (not A)
      p[1] = p[2] = p[3] = (1 - EPSILON) / 3;
      p[0] = EPSILON;
    }
    else if (c == 'D') { // D = A/G/U (not C)
      p[0] = p[2] = p[3] = (1 - EPSILON) / 3;
      p[1] = EPSILON;
    }
    else if (c == 'H') { // H = A/C/U (not G)
      p[0] = p[1] = p[3] = (1 - EPSILON) / 3;
      p[2] = EPSILON;
    }
    else if (c == 'V') { // V = A/C/G (not U)
      p[0] = p[1] = p[2] = (1 - EPSILON) / 3;
      p[3] = EPSILON;
    }
    else if (c == 'N') { // N = A/C/G/U
      p[0] = p[1] = p[2] = p[3] = 1 / abc->K;
    }
    else if (c == '-') { // an empty segmentm nothing to add here
      M = 1;
      break;
    }
    else {
      printf("R3D does not allow character %c in a motif (%s)\n", c, RMmotif);
      goto ERROR;
    }
  
    //normalize
    esl_vec_FLog(p, abc->K);
    esl_vec_FLogNorm(p, abc->K);
  }

  hmm->M   = M;
  hmm->K   = abc->K;
  hmm->abc = abc;
  
  return hmm;
  
 ERROR:
  return NULL;
}

void
R3D_hmm_Destroy(R3D_HMM *hmm)
{
  if (!hmm) return;

  if (hmm->e) {
    if (hmm->e[0]) free(hmm->e[0]);
    free(hmm->e);
  }
  free(hmm);
  return;
}

R3D_HMX *
R3D_hmx_Create(int allocL, int allocM)
{
  R3D_HMX *mx = NULL;
  int      i;
  int      status;
  
  ESL_ALLOC(mx, sizeof(R3D_HMX));
  mx->dp_mem = NULL;
  mx->Sdp    = NULL;
  mx->Idp    = NULL;

  ESL_ALLOC(mx->dp_mem, sizeof(SCVAL) * (allocL+1) * allocM);
  mx->ncells = (allocL+1) * allocM;
  
  ESL_ALLOC(mx->Sdp, sizeof (SCVAL *) * (allocL+1));
  ESL_ALLOC(mx->Idp, sizeof (SCVAL *) * (allocL+1));
  mx->allocR = allocL+1;

  for (i = 0; i <= allocL; i++) {
    mx->Sdp[i] = mx->dp_mem + i*allocM;
    mx->Idp[i] = mx->dp_mem + i*allocM;
  }
  mx->validR = allocL+1;
  mx->allocM = allocM;

  mx->L = 0;
  mx->M = 0;
  return mx;

 ERROR:
  R3D_hmx_Destroy(mx);
  return NULL;
}

int
R3D_hmx_GrowTo(R3D_HMX *mx, int L, int M)
{
  uint64_t ncells;
  void    *p;
  int      do_reset = FALSE;
  int      i;
  int      status;

  if (L < mx->allocR && M <= mx->allocM) return eslOK;

  /* Do we have to reallocate the 2D matrix, or can we get away with
   * rejiggering the row pointers into the existing memory? 
   */
  ncells = (L+1) * (M+1);
  if (ncells > mx->ncells) 
    {
      ESL_RALLOC(mx->dp_mem, p, sizeof(SCVAL) * ncells);
      mx->ncells = ncells;
      do_reset   = TRUE;
    }

  /* must we reallocate row pointers? */
  if (L >= mx->allocR)
    {
      ESL_RALLOC(mx->Sdp, p, sizeof(SCVAL *) * (L+1));
      ESL_RALLOC(mx->Idp, p, sizeof(SCVAL *) * (L+1));
      mx->allocR = L+1;
      mx->allocM = M;
      do_reset   = TRUE;
    }

  /* must we widen the rows? */
  if (M > mx->allocM)
    {
      mx->allocM = M;
      do_reset = TRUE;
    }

  /* must we set some more valid row pointers? */
  if (L >= mx->validR)
    do_reset = TRUE;

  /* did we trigger a relayout of row pointers? */
  if (do_reset)
    {
      mx->validR = ESL_MIN(mx->ncells / mx->allocM, mx->allocR);
      for (i = 0; i < mx->validR; i++)
	mx->Sdp[i] = mx->dp_mem + i*mx->allocM;
	mx->Idp[i] = mx->dp_mem + i*mx->allocM;
    }
  mx->M = 0;
  mx->L = 0;
  return eslOK;

 ERROR:
  return status;
}

void
R3D_hmx_Destroy(R3D_HMX *mx)
{
  if (mx == NULL) return;

  if (mx->dp_mem != NULL) free(mx->dp_mem);
  if (mx->Sdp    != NULL) free(mx->Sdp);
  if (mx->Idp    != NULL) free(mx->Idp);
  free(mx);
  return;
}


SCVAL
R3D_hmm_Forward(int i, int j, double **pm, R3D_HMM *hmm, R3D_HMX *fwd, char *errbuf)
{
  SCVAL fwdsc  = -eslINFINITY;
  SCVAL sumsc  = -eslINFINITY;
  SCVAL sc;	        	        /* score for a rule */
  SCVAL emit_k;
  SCVAL emit_i;
  int   L = j - i + 1;
  int   M = fwd->M;
  int   l, la;
  int   k;

  if (L <= 0) return fwdsc;

  e2_FLogsumInit();

  for (l = 0; l < L; l++) {
    la = l + i;

    emit_i = emitsc_prof(la, hmm->abc->K, pm, hmm->e[0]);

    // First S^k then I^k
    for (k = 0; k <= M; k++) {

      emit_k = (k > 0)? emitsc_prof(la, hmm->abc->K, pm, hmm->e[k]) : -eslINFINITY;
      
      // S^{0}[l=0) = log(1)
      // S^{0}[l>0) = log(0)
      //
      if (k == 0) {
	fwd->Sdp[k][l] = (l == 0)? 0 : -eslINFINITY;
      }
      // S^{1}
      //
      else if (k == 1) {
	sumsc = -eslINFINITY;
	//  tB     S^{0}(l-1)
	sc    = (l>0)? fwd->Sdp[0][l-1] + hmm->tB[0] : -eslINFINITY;
	sumsc = e2_FLogsum(sumsc, sc);
 	//  (1-tI) I^{0}(l)
	sc    = fwd->Idp[0][l] + hmm->tI[1];
	sumsc = e2_FLogsum(sumsc, sc);
	//
	fwd->Sdp[k][l] = sumsc;
      }
      // S^{k>1}[l]
      //
      else {
	sumsc = -eslINFINITY;
	// t_M P_k(x_l) S^{k-1}(l-1)
	sc    = (l>0)? fwd->Sdp[k-1][l-1] + hmm->tS[0] + emit_k : -eslINFINITY;
	sumsc = e2_FLogsum(sumsc, sc);
	// t_D          S^{k-1}(l)
	sc    = fwd->Sdp[k-1][l] + hmm->tS[2];
	sumsc = e2_FLogsum(sumsc, sc);
	// (1-tI)       I^{k-1}(l)
	sc    = fwd->Idp[k-1][l] + hmm->tI[1];
	sumsc = e2_FLogsum(sumsc, sc);
	//
	fwd->Sdp[k][l] = sumsc;
      }
      
      // I^{0}
      //
      if (k == 0) {
	sumsc = -eslINFINITY;
	// (1-tB) S^{0}(l)
	sc    = fwd->Sdp[0][l] + hmm->tB[1];
	sumsc = e2_FLogsum(sumsc, sc);
	// tI P_I(x_l) I^{0}(l-1)
	sc    = (l>0)? fwd->Idp[0][l-1] + hmm->tI[0] + emit_i : -eslINFINITY;
	sumsc = e2_FLogsum(sumsc, sc);
	//
	fwd->Idp[0][l] = sumsc;
      }
      // I^{k>0}[l]
      //
      else {
	sumsc = -eslINFINITY;
	// t_MI P_k(x_l) S^{k}(l-1)
	sc    = (l>0)? fwd->Sdp[k-1][l-1] + hmm->tS[1] + emit_k : -eslINFINITY;
	sumsc = e2_FLogsum(sumsc, sc);
	// t_DI          S^{k}(l)
	sc    = fwd->Sdp[k][l] + hmm->tS[3];
	sumsc = e2_FLogsum(sumsc, sc);
	// tI   P_I(x_l) I^{k}(l-1)
	sc    = (l>0)? fwd->Idp[k][l-1] + hmm->tI[0] + emit_i : -eslINFINITY;
	sumsc = e2_FLogsum(sumsc, sc);
	//
	fwd->Idp[k][l] = sumsc;
      }
    }
  }
  
  // End state E(j)
  //
  // t_M P_M(x_j) S^M(j-1)
  sc    = (L>1)? fwd->Sdp[M][L-2] + hmm->tS[0] + emitsc_prof(j, hmm->abc->K, pm, hmm->e[M]) : -eslINFINITY;
  fwdsc = e2_FLogsum(fwdsc, sc);
  // t_D S^M(j)
  sc    = fwd->Sdp[M][L-1]   + hmm->tS[2];
  fwdsc = e2_FLogsum(fwdsc, sc);                 	        
  // (1-t_I) I^M(j)
  sc    = fwd->Idp[M][L-1]   + hmm->tI[1];
  fwdsc = e2_FLogsum(fwdsc, sc);

  if (fwdsc <= -eslINFINITY) esl_fail(errbuf, "R3D HMM failed: Forward sc = -inf.");

  return fwdsc;
}


static SCVAL
emitsc_prof(int i, int K, double **pm, SCVAL *e)
{
  float  sc = -eslINFINITY;
  float  logval;
  int    k;

  if (i < 1) return (SCVAL)sc;

  e2_FLogsumInit();
  
  for (k = 0; k < K; k ++) {
    logval = log(pm[i-1][k]) + (float)e[k];
    sc = e2_FLogsum(sc, logval);
  }
  // columns with all gaps
  if (sc == -eslINFINITY) sc = log(0.25);
  
  return (SCVAL)sc;
}
