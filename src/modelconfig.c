/* Model configuration: 
 * Converting a core model to a fully configured Plan7 search profile.
 * 
 * Contents:
 *     1. Routines in the exposed API.
 *     2. Unit tests.
 *     3. Test driver.
 *     4. Statistics collection driver.
 *     5. Copyright and license
 * 
*/
#include "p7_config.h"

#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2.h"
#include "e2_msa.h"
#include "e2_profile.h"
#include "e2hmmer_profile.h"
#include "logsum.h"
#include "modelconfig.h"


/*****************************************************************
 * 1. Routines in the exposed API.
 *****************************************************************/
 
/* Function:  e2_ProfileConfig()
 * Synopsis:  Configure a search profile.
 *
 * Purpose:   Given two evolutionary modesl, the ancestral
 *            model <bg>, a desired search <mode> (one of 
 *            <e2_GLOBAL> or <e2_LOCAL>), and an
 *            expected ancestral sequence length <L>; configure the
 *            search model in <gm>.
 *            
 * Returns:   <eslOK> on success; the profile <gm> now contains 
 *            scores and is ready for searching target sequences.
 *            
 * Throws:    <eslEMEM> on allocation error.
 */
int 
e2_ProfileConfig(const E1_MODEL *evoml, const E1_MODEL *evomr, float *ancf, E2_PROFILE *gm, float L, E2_ALI e2ali, int verbose)
{
  float  *tsc = gm->tsc;
  float   sc;
  float   p;
  float   ddfactor;
  float   nj;
  float   lj;
  float   pmove, ploop;
  int     e2fix = FALSE;
  int     b, c;
  int     a;
  int     status;

  if (e2ali == E2HMMER || e2ali == E2FHMMER) { status = eslFAIL; goto ERROR; }
  if (e2ali == E2F) e2fix = TRUE;
  
  e2_FLogdiffInit();
  
  /* Contract checks */
  if (evoml->mode      != evomr->mode)      ESL_XEXCEPTION(eslEINVAL, "evom's don't match");
  if (evoml->abc->type != evomr->abc->type) ESL_XEXCEPTION(eslEINVAL, "evom's don't match");
  
  /* Copy some pointer references and other info across from evom's  */
  gm->mode = evoml->mode; 
  if (gm->name != NULL) free(gm->name);
  if (gm->acc  != NULL) free(gm->acc);
  if (gm->desc != NULL) free(gm->desc);
  if ((status = esl_strdup(evoml->name,   -1, &(gm->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(evoml->acc,    -1, &(gm->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(evoml->desc,   -1, &(gm->desc))) != eslOK) goto ERROR;
   
  /* fixed for rev models, set to msa length otherwise */
  p = (evoml->R->p > 0.0 && evoml->R->p < 1.0)? evoml->R->p : (double)L / ((double)L + 1.);
  
  tsc[e2P_DD_DD] = log(evoml->t[e1H_DD]) + log(evomr->t[e1H_DD]);
  ddfactor = (e2fix)? 0.0 : -log(1.0 - p * evoml->t[e1H_DD] * evomr->t[e1H_DD]);

  /* the extra states */
  /*            Parameter <nj> controls the unihit/multihit behavior of
   *            the profile. <nj> is the expected number of J segments.
   *            <nj=0> is unihit. <nj=1> is the standard multihit model,
   *            with <t_{EJ} = 0.5>.
   */
  nj = 0;
  if      (gm->mode == e2_GLOBAL) nj = 0.0;
  else if (gm->mode == e2_LOCAL)  nj = 0.5;
  gm->xsc[e2P_EE][e2P_MOVE] = logf ( 1.0f / (1.0f + nj));       // E->C1
  gm->xsc[e2P_EE][e2P_LOOP] = logf ( nj / (1.0f + nj));         // E->J1

  /* Configure N,J,C transitions so they bear lj/3 of the total
   * sequence length L for local mode. 
   */
  lj = (float)L/0.2;

  pmove = 1.0;
  if      (gm->mode == e2_GLOBAL) pmove = 1.0;                            /* 1.0     for global */
  else if (gm->mode == e2_LOCAL)  pmove = (2.0f + nj) / (lj + 2.0f + nj); /*         for local  */
  ploop = 1.0f - pmove;

  gm->xsc[e2P_N1][e2P_LOOP] = gm->xsc[e2P_C1][e2P_LOOP] = gm->xsc[e2P_J1][e2P_LOOP] = logf(ploop);
  gm->xsc[e2P_N2][e2P_LOOP] = gm->xsc[e2P_C2][e2P_LOOP] = gm->xsc[e2P_J2][e2P_LOOP] = logf(ploop);
  gm->xsc[e2P_N1][e2P_MOVE] = gm->xsc[e2P_C1][e2P_MOVE] = gm->xsc[e2P_J1][e2P_MOVE] = logf(pmove);
  gm->xsc[e2P_N2][e2P_MOVE] = gm->xsc[e2P_C2][e2P_MOVE] = gm->xsc[e2P_J2][e2P_MOVE] = logf(pmove);
  
  /* Transition scores. */
  if (verbose) {
    //e1_rate_Dump(stdout, evoml->R);
    e1_model_DumpTransitions(stdout, evoml);
    //e1_model_DumpEmissions(stdout, evoml);
    //e1_rate_Dump(stdout, evomr->R);
    e1_model_DumpTransitions(stdout, evomr);
  }

  tsc[e2P_BB_IB] = log(evoml->t[e1H_BI]);
  tsc[e2P_IB_IB] = log(evoml->t[e1H_II]);

  tsc[e2P_BB_SS] = log(evoml->t[e1H_BS]) + log(evomr->t[e1H_BS]) + log(p);
  tsc[e2P_IB_SS] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_BS]) + log(p);
  tsc[e2P_SS_SS] = log(evoml->t[e1H_SS]) + log(evomr->t[e1H_SS]) + log(p);
  tsc[e2P_DS_SS] = log(evoml->t[e1H_DS]) + log(evomr->t[e1H_SS]) + log(p);
  tsc[e2P_IS_SS] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_SS]) + log(p);
  tsc[e2P_SD_SS] = log(evoml->t[e1H_SS]) + log(evomr->t[e1H_DS]) + log(p);
  tsc[e2P_DD_SS] = log(evoml->t[e1H_DS]) + log(evomr->t[e1H_DS]) + log(p);
  tsc[e2P_ID_SS] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_DS]) + log(p);
  tsc[e2P_BI_SS] = log(evoml->t[e1H_BS]) + log(evomr->t[e1H_IS]) + log(p);
  tsc[e2P_SI_SS] = log(evoml->t[e1H_SS]) + log(evomr->t[e1H_IS]) + log(p);
  tsc[e2P_DI_SS] = log(evoml->t[e1H_DS]) + log(evomr->t[e1H_IS]) + log(p);
  tsc[e2P_II_SS] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_IS]) + log(p);

  tsc[e2P_BB_DS] = log(evoml->t[e1H_BD]) + log(evomr->t[e1H_BS]) + log(p);
  tsc[e2P_IB_DS] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_BS]) + log(p);
  tsc[e2P_SS_DS] = log(evoml->t[e1H_SD]) + log(evomr->t[e1H_SS]) + log(p);
  tsc[e2P_DS_DS] = log(evoml->t[e1H_DD]) + log(evomr->t[e1H_SS]) + log(p);
  tsc[e2P_IS_DS] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_SS]) + log(p);
  tsc[e2P_SD_DS] = log(evoml->t[e1H_SD]) + log(evomr->t[e1H_DS]) + log(p);
  tsc[e2P_DD_DS] = log(evoml->t[e1H_DD]) + log(evomr->t[e1H_DS]) + log(p);
  tsc[e2P_ID_DS] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_DS]) + log(p);
  tsc[e2P_BI_DS] = log(evoml->t[e1H_BD]) + log(evomr->t[e1H_IS]) + log(p);
  tsc[e2P_SI_DS] = log(evoml->t[e1H_SD]) + log(evomr->t[e1H_IS]) + log(p);
  tsc[e2P_DI_DS] = log(evoml->t[e1H_DD]) + log(evomr->t[e1H_IS]) + log(p);
  tsc[e2P_II_DS] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_IS]) + log(p);

  tsc[e2P_SS_IS] = log(evoml->t[e1H_SI]);
  tsc[e2P_DS_IS] = log(evoml->t[e1H_DI]);  
  tsc[e2P_IS_IS] = log(evoml->t[e1H_II]);

  tsc[e2P_BB_SD] = log(evoml->t[e1H_BS]) + log(evomr->t[e1H_BD]) + log(p);
  tsc[e2P_IB_SD] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_BD]) + log(p);
  tsc[e2P_SS_SD] = log(evoml->t[e1H_SS]) + log(evomr->t[e1H_SD]) + log(p);
  tsc[e2P_DS_SD] = log(evoml->t[e1H_DS]) + log(evomr->t[e1H_SD]) + log(p);
  tsc[e2P_IS_SD] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_SD]) + log(p);
  tsc[e2P_SD_SD] = log(evoml->t[e1H_SS]) + log(evomr->t[e1H_DD]) + log(p);
  tsc[e2P_DD_SD] = log(evoml->t[e1H_DS]) + log(evomr->t[e1H_DD]) + log(p);
  tsc[e2P_ID_SD] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_DD]) + log(p);
  tsc[e2P_BI_SD] = log(evoml->t[e1H_BS]) + log(evomr->t[e1H_ID]) + log(p);
  tsc[e2P_SI_SD] = log(evoml->t[e1H_SS]) + log(evomr->t[e1H_ID]) + log(p);
  tsc[e2P_DI_SD] = log(evoml->t[e1H_DS]) + log(evomr->t[e1H_ID]) + log(p);
  tsc[e2P_II_SD] = log(evoml->t[e1H_IS]) + log(evomr->t[e1H_ID]) + log(p);

  tsc[e2P_IB_DD] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_BD]) + ddfactor;
  tsc[e2P_SS_DD] = log(evoml->t[e1H_SD]) + log(evomr->t[e1H_SD]) + ddfactor;
  tsc[e2P_DS_DD] = log(evoml->t[e1H_DD]) + log(evomr->t[e1H_SD]) + ddfactor;
  tsc[e2P_IS_DD] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_SD]) + ddfactor;
  tsc[e2P_SD_DD] = log(evoml->t[e1H_SD]) + log(evomr->t[e1H_DD]) + ddfactor;
  tsc[e2P_ID_DD] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_DD]) + ddfactor;
  tsc[e2P_BI_DD] = log(evoml->t[e1H_BD]) + log(evomr->t[e1H_ID]) + ddfactor;
  tsc[e2P_SI_DD] = log(evoml->t[e1H_SD]) + log(evomr->t[e1H_ID]) + ddfactor;
  tsc[e2P_DI_DD] = log(evoml->t[e1H_DD]) + log(evomr->t[e1H_ID]) + ddfactor;
  tsc[e2P_II_DD] = log(evoml->t[e1H_ID]) + log(evomr->t[e1H_ID]) + ddfactor;

  tsc[e2P_SD_ID] = log(evoml->t[e1H_SI]);
  tsc[e2P_DD_ID] = log(evoml->t[e1H_DI]);
  tsc[e2P_ID_ID] = log(evoml->t[e1H_II]);

  tsc[e2P_BB_BI] = log(evomr->t[e1H_BI]);
  tsc[e2P_BI_BI] = log(evomr->t[e1H_II]);

  tsc[e2P_SS_SI] = log(evomr->t[e1H_SI]);
  tsc[e2P_SD_SI] = log(evomr->t[e1H_DI]);
  tsc[e2P_SI_SI] = log(evomr->t[e1H_II]);

  tsc[e2P_DS_DI] = log(evomr->t[e1H_SI]);
  tsc[e2P_DD_DI] = log(evomr->t[e1H_DI]);
  tsc[e2P_DI_DI] = log(evomr->t[e1H_II]);
 
  tsc[e2P_IB_II] = log(evomr->t[e1H_BI]);
  tsc[e2P_IS_II] = log(evomr->t[e1H_SI]);
  tsc[e2P_ID_II] = log(evomr->t[e1H_DI]);
  tsc[e2P_II_II] = log(evomr->t[e1H_II]);
 
  tsc[e2P_IB_EE] = log(evoml->t[e1H_IE]) + log(evomr->t[e1H_BE]) + log(1.0-p); 
  tsc[e2P_SS_EE] = log(evoml->t[e1H_SE]) + log(evomr->t[e1H_SE]) + log(1.0-p); 
  tsc[e2P_DS_EE] = log(evoml->t[e1H_DE]) + log(evomr->t[e1H_SE]) + log(1.0-p); 
  tsc[e2P_IS_EE] = log(evoml->t[e1H_IE]) + log(evomr->t[e1H_SE]) + log(1.0-p); 
  tsc[e2P_SD_EE] = log(evoml->t[e1H_SE]) + log(evomr->t[e1H_DE]) + log(1.0-p); 
  tsc[e2P_DD_EE] = log(evoml->t[e1H_DE]) + log(evomr->t[e1H_DE]) + log(1.0-p); 
  tsc[e2P_ID_EE] = log(evoml->t[e1H_IE]) + log(evomr->t[e1H_DE]) + log(1.0-p); 
  tsc[e2P_BI_EE] = log(evoml->t[e1H_BE]) + log(evomr->t[e1H_IE]) + log(1.0-p); 
  tsc[e2P_SI_EE] = log(evoml->t[e1H_SE]) + log(evomr->t[e1H_IE]) + log(1.0-p); 
  tsc[e2P_DI_EE] = log(evoml->t[e1H_DE]) + log(evomr->t[e1H_IE]) + log(1.0-p); 
  tsc[e2P_II_EE] = log(evoml->t[e1H_IE]) + log(evomr->t[e1H_IE]) + log(1.0-p); 

 /* Double substition emission scores. */
  for (b = 0; b < gm->abc->K; b++) 
    for (c = 0; c < gm->abc->K; c++) 
      {
	sc = -eslINFINITY;
	for (a = 0;  a < gm->abc->K; a++) {
	  sc = p7_FLogsum(sc, log(ancf[a]) + log(evoml->sub->mx[a][b]) + log(evomr->sub->mx[a][c]));
	}
	gm->sssc[b][c] = sc;	  
      }
  
  /* Orphan substition emission scores. */
  for (b = 0; b < gm->abc->K; b++) {
    sc = -eslINFINITY;
    for (a = 0;  a < gm->abc->K; a++)
      sc = p7_FLogsum(sc, log(ancf[a]) + log(evoml->sub->mx[a][b]));
    gm->ssc[b][e2P_SL] = sc;	  

    sc = -eslINFINITY;
    for (a = 0;  a < gm->abc->K; a++)
      sc = p7_FLogsum(sc, log(ancf[a]) + log(evomr->sub->mx[a][b]));
    gm->ssc[b][e2P_SR] = sc;	  
  }
  
  /* Insertion emission scores. */
  for (b = 0; b < gm->abc->K; b++) { 
    gm->isc[b][e2P_SL] = log(evoml->ins[b]);
    gm->isc[b][e2P_SR] = log(evomr->ins[b]);
  }
  
  /* Flanking emission scores. */
  for (b = 0; b < gm->abc->K; b++) { 
    gm->fsc[b][e2P_SL] = log(evoml->ins[b]);
    gm->fsc[b][e2P_SR] = log(evomr->ins[b]);
  }
  
 return eslOK;
  
 ERROR:
  return status;
}

/* Function:  e2hmmer_ProfileConfig()
 * Synopsis:  Configure a search profile.
 *
 * Purpose:   Given two evolutionary modesl, the ancestral
 *            model <bg>, a desired search <mode> (one of 
 *            <e2_GLOBAL> or <e2_LOCAL>), and an
 *            expected ancestral sequence length <L>; configure the
 *            search model in <gm>.
 *            
 * Returns:   <eslOK> on success; the profile <gm> now contains 
 *            scores and is ready for searching target sequences.
 *            
 * Throws:    <eslEMEM> on allocation error.
 */
int 
e2hmmer_ProfileConfig(const P7_RATE *R, float tl, float tr, const P7_HMM *evoml, const P7_HMM *evomr, P7_BG *bg, 
		      E2HMMER_PROFILE *gm, float L, E2_ALI e2ali, int mode, int verbose)
{
  float  *tsc;
  float  *trl;
  float  *trr;
  float   sc;
  int     dtl, dtr;      /* discretized times */
  int     e2fix = FALSE;
  int     k;
  int     b, c;
  int     a;
  int     x;
  int     status;

  if (e2ali == E2 || e2ali == E2F) { status = eslFAIL; goto ERROR; }
  if (e2ali == E2FHMMER) e2fix = TRUE;

  e2_FLogdiffInit();

  /* discretize times for substitutions */
  dtl = select_discrete_time(R, tl); if (dtl < 0) { status = eslFAIL; goto ERROR; }
  dtr = select_discrete_time(R, tr); if (dtr < 0) { status = eslFAIL; goto ERROR; }
 
   /* Contract checks */
  if (evoml->abc->type != evomr->abc->type) ESL_XEXCEPTION(eslEINVAL, "evom's don't match");

  /* Copy some pointer references and other info across from evom's  */
  gm->M = evoml->M;
  gm->mode = mode; 
  if (mode != e2_GLOBAL) { printf("only GLOBAL model allowed\n"); status = eslFAIL; goto ERROR; }
 
  if (gm->name != NULL) free(gm->name);
  if (gm->acc  != NULL) free(gm->acc);
  if (gm->desc != NULL) free(gm->desc);
  if ((status = esl_strdup(evoml->name,   -1, &(gm->name))) != eslOK) { status = eslFAIL; goto ERROR; }
  if ((status = esl_strdup(evoml->acc,    -1, &(gm->acc)))  != eslOK) { status = eslFAIL; goto ERROR; }
  if ((status = esl_strdup(evoml->desc,   -1, &(gm->desc))) != eslOK) { status = eslFAIL; goto ERROR; }
   
  /* change the bg model */
  bg->p1 = L / (L + 1.);

  for (k = 0; k <= gm->M; k++) {
    tsc = gm->tsc[k];
    trl = evoml->t[k];
    trr = evomr->t[k];

    /* Transition scores. Non standard */
    tsc[e2HP_BB_IB] = log(0.495/2.);
    tsc[e2HP_BB_BI] = log(0.495/2.);
    tsc[e2HP_BB_SS] = log(0.495) + log(2.0/(gm->M*(gm->M+1))) + log(bg->p1);
    tsc[e2HP_BB_EE] = log(0.01);

    tsc[e2HP_IB_IB] = log(0.495/2.);
    tsc[e2HP_IB_iI] = log(0.495/2.);
    tsc[e2HP_IB_SS] = log(0.495) + log(bg->p1);
    tsc[e2HP_IB_EE] = log(0.01); 
    
    tsc[e2HP_BI_Ii] = log(0.30);
    tsc[e2HP_BI_BI] = log(0.30);
    tsc[e2HP_BI_SS] = log(0.30) + log(bg->p1);
    tsc[e2HP_BI_EE] = log(0.01); 
    
   /* Transition scores. The standar ones */
    tsc[e2HP_SS_SS] = log(trl[p7H_MM]) + log(trr[p7H_MM]) + log(bg->p1);
    tsc[e2HP_DS_SS] = log(trl[p7H_DM]) + log(trr[p7H_MM]) + log(bg->p1);
    tsc[e2HP_IS_SS] = log(trl[p7H_IM]) + log(trr[p7H_MM]) + log(bg->p1);
    tsc[e2HP_SD_SS] = log(trl[p7H_MM]) + log(trr[p7H_DM]) + log(bg->p1);
    tsc[e2HP_DD_SS] = log(trl[p7H_DM]) + log(trr[p7H_DM]) + log(bg->p1);
    tsc[e2HP_ID_SS] = log(trl[p7H_IM]) + log(trr[p7H_DM]) + log(bg->p1);
    tsc[e2HP_SI_SS] = log(trl[p7H_MM]) + log(trr[p7H_IM]) + log(bg->p1);
    tsc[e2HP_DI_SS] = log(trl[p7H_DM]) + log(trr[p7H_IM]) + log(bg->p1);
    tsc[e2HP_Ii_SS] = log(trl[p7H_IM]) + log(trr[p7H_IM]) + log(bg->p1);
    tsc[e2HP_iI_SS] = log(trl[p7H_IM]) + log(trr[p7H_IM]) + log(bg->p1);
    
    tsc[e2HP_SS_DS] = log(trl[p7H_MD]) + log(trr[p7H_MM]) + log(bg->p1);
    tsc[e2HP_DS_DS] = log(trl[p7H_DD]) + log(trr[p7H_MM]) + log(bg->p1);
    tsc[e2HP_SD_DS] = log(trl[p7H_MD]) + log(trr[p7H_DM]) + log(bg->p1);
    tsc[e2HP_DD_DS] = log(trl[p7H_DD]) + log(trr[p7H_DM]) + log(bg->p1);
    tsc[e2HP_SI_DS] = log(trl[p7H_MD]) + log(trr[p7H_IM]) + log(bg->p1);
    tsc[e2HP_DI_DS] = log(trl[p7H_DD]) + log(trr[p7H_IM]) + log(bg->p1);
    
    tsc[e2HP_SS_IS] = log(trl[p7H_MI]);
    tsc[e2HP_IS_IS] = log(trl[p7H_II]);
    
    tsc[e2HP_SS_SD] = log(trl[p7H_MM]) + log(trr[p7H_MD]) + log(bg->p1);
    tsc[e2HP_DS_SD] = log(trl[p7H_DM]) + log(trr[p7H_MD]) + log(bg->p1);
    tsc[e2HP_IS_SD] = log(trl[p7H_IM]) + log(trr[p7H_MD]) + log(bg->p1);
    tsc[e2HP_SD_SD] = log(trl[p7H_MM]) + log(trr[p7H_DD]) + log(bg->p1);
    tsc[e2HP_DD_SD] = log(trl[p7H_DM]) + log(trr[p7H_DD]) + log(bg->p1);
    tsc[e2HP_ID_SD] = log(trl[p7H_IM]) + log(trr[p7H_DD]) + log(bg->p1);
   
    tsc[e2HP_SS_DD] = log(trl[p7H_MD]) + log(trr[p7H_MD]) + log(bg->p1);
    tsc[e2HP_DS_DD] = log(trl[p7H_DD]) + log(trr[p7H_MD]) + log(bg->p1);
    tsc[e2HP_SD_DD] = log(trl[p7H_MD]) + log(trr[p7H_DD]) + log(bg->p1);
    tsc[e2HP_DD_DD] = log(trl[p7H_DD]) + log(trr[p7H_DD]) + log(bg->p1);
    
    tsc[e2HP_SD_ID] = log(trl[p7H_MI]);
    tsc[e2HP_ID_ID] = log(trl[p7H_II]);
   
    tsc[e2HP_SS_SI] = log(trr[p7H_MI]);
    tsc[e2HP_SI_SI] = log(trr[p7H_II]);
    
    tsc[e2HP_DS_DI] = log(trr[p7H_MI]);
    tsc[e2HP_DI_DI] = log(trr[p7H_II]);
    
 
    tsc[e2HP_SI_Ii] = log(trl[p7H_MI]);
    tsc[e2HP_Ii_Ii] = log(trl[p7H_II]);
    tsc[e2HP_iI_Ii] = log(trl[p7H_II]);


    tsc[e2HP_IS_iI] = log(trr[p7H_MI]);
    tsc[e2HP_Ii_iI] = log(trr[p7H_II]);
    tsc[e2HP_iI_iI] = log(trr[p7H_II]);
    
    tsc[e2HP_SS_EE] = log(1.0-bg->p1);  
    tsc[e2HP_DS_EE] = log(1.0-bg->p1); 
    tsc[e2HP_IS_EE] = log(1.0-bg->p1); 
    tsc[e2HP_SD_EE] = log(1.0-bg->p1); 
    tsc[e2HP_DD_EE] = log(1.0-bg->p1); 
    tsc[e2HP_ID_EE] = log(1.0-bg->p1); 
    tsc[e2HP_SI_EE] = log(1.0-bg->p1); 
    tsc[e2HP_DI_EE] = log(1.0-bg->p1); 
    tsc[e2HP_Ii_EE] = log(1.0-bg->p1); 
    tsc[e2HP_iI_EE] = log(1.0-bg->p1); 
    
#if 1
    for (x = 0; x < e2HP_NTRANS; x ++) {
      if (isnan(tsc[x])) { printf("nan transition for k = %d x=%d\n", k, x); exit(1); }
    }
#endif

#if 0
    for (x = 0; x < e2HP_NTRANS; x ++) {
      if (x==e2HP_SS_SI) printf("S->I k %d/%d x %d tsc %f\n", k, gm->M, x, tsc[x]);
    }
#endif

    /* Double substition emission scores. */ 
    for (b = 0; b < gm->abc->K; b++) 
      for (c = 0; c < gm->abc->K; c++) 
	{
	  sc = -eslINFINITY;
	  for (a = 0;  a < gm->abc->K; a++) {
	    sc = p7_FLogsum(sc, log(R->pzero[k][a]) + log(R->Pdt[dtl][k]->mx[a][b]) + log(R->Pdt[dtr][k]->mx[a][c]));
	  }
	  gm->sssc[k][b][c] = sc;
	}
    
    /* Orphan substition emission scores. */
    for (b = 0; b < gm->abc->K; b++) {

      sc = -eslINFINITY;
      for (a = 0;  a < gm->abc->K; a++) {
	sc = p7_FLogsum(sc, log(R->pzero[k][a]) + log(R->Pdt[dtl][k]->mx[a][b]));
      }
      gm->ssc[k][b][e2HP_SL] = sc;	  
      
      sc = -eslINFINITY;
      for (a = 0;  a < gm->abc->K; a++) {
	sc = p7_FLogsum(sc, log(R->pzero[k][a]) + log(R->Pdt[dtr][k]->mx[a][b]));
      }
      gm->ssc[k][b][e2HP_SR] = sc;	  
    }
  
    /* Insertion emission scores. */
    for (b = 0; b < gm->abc->K; b++) { 
      gm->isc[k][b][e2HP_SL] = log(evoml->ins[k][b]);
      gm->isc[k][b][e2HP_SR] = log(evomr->ins[k][b]);
    }
  }
  
 return eslOK;
  
 ERROR:
  return status;
}

int
select_discrete_time(const P7_RATE *R, float time)
{
  int t;
  int dt = -1;

   if (time < 1.0) {
      for (t = 0; t < R->ndt; t++) {
	if (R->dtval[t] >= time) { dt = t; break; }
      }
      if (t == R->ndt) dt = R->ndt;	
    }
    else {
      for (t = R->ndt-1; t >= 0; t--) {
	if (R->dtval[t] <= time) { dt = t; break; }
      }
      if (t < 0) dt = R->ndt;	
    }

   return dt;
}
