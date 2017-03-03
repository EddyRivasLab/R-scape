/* Optimal accuracy alignment; generic version.
 * 
 * Contents:
 *   1. Optimal alignment accuracy fill.
 *   2. Optimal alignment accuracy traceback.
 *   3. Benchmark driver
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 *   7. Copyright and license information
 * 
 */

#include "p7_config.h"

#include <float.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_stack.h"
 #include "esl_vectorops.h"

#include "hmmer.h"
#include "e2.h"
#include "e2hmmer_profile.h"
#include "e2fhmmer_generic_optacc.h"
#include "e2_profilesq.h"
#include "modelconfig.h"

/*****************************************************************
 * 1. Optimal alignment fill and traceback.
 *****************************************************************/

#define TSCDELTA(k,s) ( (tsc[(k)][(s)] == -eslINFINITY) ? FLT_MIN : 1.0)

/* The TSCDELTA is used to make impossible paths impossible in the
 * optimal accuracy decoding algorithm; see Kall et al (2005). What we
 * want to do is multiply by a Kronecker delta that's 1 when the
 * transition probability is finite, and 0 when it's zero (when the
 * log prob is -eslINFINITY). But we can't do that easily, when we're_
 * in log space, because 0 * -eslINFINITY = NaN. Instead, we use a
 * tiny number (FLT_MIN, ~1e-37).
 * 
 * A side concern is that we don't want to put a bunch of if-else
 * branches in the code; compilers should be able to generate more
 * efficient code from the TSCDELTA() construction.
 */


static int isgap(float *p, int K);

/* Function:  e2fhmmer_GOptimalAccuracy()
 * Synopsis:  Optimal accuracy decoding: fill. 
 * Incept:    ER, Thu Oct 10 09:40:51 EDT 2013
 *            adapted from p7_GOptimalAccuracy()
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm \citep{Kall05}.
 *            
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <L> using the query model <gm>
 *            
 *            Caller also provides a DP matrix <gx>. The routine fills this in
 *            with OA scores.
 *            
 * Args:      gm    - query profile      
 *            pp    - posterior decoding matrix created by <e2_GPosteriorDecoding()>
 *            gx    - RESULT: caller provided DP matrix for <gm->M> by <L> 
 *            ret_e - RETURN: expected number of correctly decoded positions 
 *
 * Returns:   <eslOK> on success, and <*ret_e> contains the final OA
 *            score, which is the expected number of correctly decoded
 *            positions in the target sequence (up to <L>).
 *
 * Throws:    (no abnormal error conditions)
 */
int
e2fhmmer_GOptimalAccuracy(const E2HMMER_PROFILE *gm, const E2_GMX *pp, E2_GMX *gx, const PSQ *psql, const PSQ *psqr, float *ret_e, int *ret_k)
{
  float const  **tsc = (float const **)gm->tsc;
  float        **dp  = gx->dp;						
  float         sc;
  float         EE = -eslINFINITY;
  int           L;
  int           isgapl, isgapr;
  int           apair;
  int           x;  
  int           xv; /* linear memory indices */
  int           kmax;
  int           k, kv, kk;
  
  /* OptAcc:
   *           states DD, EE needs to be evaluated last.
   *
   * Order: BB, SS, DS, SD, IB, IS, ID, BI, SI, II, DD, EE
   */
 
  if (gx->Lcol != gx->Lrow) { printf("sqs are not aligned\n"); return eslFAIL; }
  L = gx->Lcol;

  x = 0;
  k = 0;
  BBMXM(x,k) = 0.0;
  SSMXM(x,k) = DSMXM(x,k) = SDMXM(x,k) = DDMXM(x,k) = -eslINFINITY;
  IBMXM(x,k) = ISMXM(x,k) = IDMXM(x,k)              = -eslINFINITY;
  BIMXM(x,k) = SIMXM(x,k) = DIMXM(x,k)              = -eslINFINITY;
  IiMXM(x,k) = iIMXM(x,k)                           = -eslINFINITY;
  
  for (k = 1; k <= gx->M; k++) {
    BBMXM(x,k) = 0.0;
    SSMXM(x,k) = DSMXM(x,k) = SDMXM(x,k) = DDMXM(x,k) = -eslINFINITY;
    IBMXM(x,k) = ISMXM(x,k) = IDMXM(x,k)              = -eslINFINITY;
    BIMXM(x,k) = SIMXM(x,k) = DIMXM(x,k)              = -eslINFINITY;
    IiMXM(x,k) = iIMXM(x,k)                           = -eslINFINITY;
    if (x == L) {
      sc = ESL_MAX(ESL_MAX(TSCDELTA(k,e2HP_BB_EE) * BBMXM(x,k), 
			   TSCDELTA(k,e2HP_IB_EE) * IBMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_SS_EE) * SSMXM(x,k),
			   TSCDELTA(k,e2HP_DS_EE) * DSMXM(x,k)));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(k,e2HP_IS_EE) * ISMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_SD_EE) * SDMXM(x,k), 
			   TSCDELTA(k,e2HP_DD_EE) * DDMXM(x,k)));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(k,e2HP_ID_EE) * IDMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_BI_EE) * BIMXM(x,k), 
			   TSCDELTA(k,e2HP_SI_EE) * SIMXM(x,k)));
      sc = ESL_MAX(ESL_MAX(sc,
			   TSCDELTA(k,e2HP_DI_EE) * DIMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_Ii_EE) * IiMXM(x,k), 
			   TSCDELTA(k,e2HP_iI_EE) * iIMXM(x,k)));
      if (sc > EE) { EE = sc; kmax = k; }
    }
  
#if 0
    if (x==L) printf("e2fhmmer_OA x %d k %d BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		     x, k, 
		     BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), 
		     SDMXM(x,k), DDMXM(x,k), IDMXM(x,k), 
		     BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IIMXM(x,k), EE);
#endif
  }

  k = 0;
  for (x = 1; x <= L; x++) {
    xv = x - 1;

    isgapl = isgap(psql->prof[x], psql->abc->K);
    isgapr = isgap(psqr->prof[x], psql->abc->K);
    if      (!isgapl && !isgapr) apair = PAIR;
    else if ( isgapl && !isgapr) apair = RRES;
    else if (!isgapl &&  isgapr) apair = LRES;
    else                         apair = NONE;

    /* BB state  0 transitions */
    BBMXM(x,k) = -eslINFINITY;
    if (apair == NONE) BBMXM(x,k) = ESL_MAX(BBMXM(x,k), BBMXM(xv,k));
    /* SS state 0 transitions */
    SSMXM(x,k) = -eslINFINITY;	  
    if (apair == NONE) SSMXM(x,k) = ESL_MAX(SSMXM(x,k), SSMXM(xv,k));
    /* DS state 0 transitions */
    DSMXM(x,k) = -eslINFINITY;	  
     if (apair == NONE) DSMXM(x,k) = ESL_MAX(DSMXM(x,k), DSMXM(xv,k));
   /* SD state 0 transitions */
    SDMXM(x,k) = -eslINFINITY;	  
     if (apair == NONE) SDMXM(x,k) = ESL_MAX(SDMXM(x,k), SDMXM(xv,k));
   /* DD state 0 transitions */
    DDMXM(x,k) = -eslINFINITY;
    if (apair == NONE) DDMXM(x,k) = ESL_MAX(DDMXM(x,k), DDMXM(xv,k));

    /* IB state 2 transitions */
    IBMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(TSCDELTA(k,e2HP_BB_IB) * (BBMXM(xv,k) + pp->dp[x][e2G_IB]), 
		   TSCDELTA(k,e2HP_IB_IB) * (IBMXM(xv,k) + pp->dp[x][e2G_IB]));
      IBMXM(x,k) = sc;
    }
    if (apair == NONE) IBMXM(x,k) = ESL_MAX(IBMXM(x,k), IBMXM(xv,k));
    /* IS state 2 transitions  */
    ISMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(TSCDELTA(k,e2HP_SS_IS) * (SSMXM(xv,k) + pp->dp[x][e2G_IS]),
		   TSCDELTA(k,e2HP_IS_IS) * (ISMXM(xv,k) + pp->dp[x][e2G_IS]));
      ISMXM(x,k) = sc;	  
    }
    if (apair == NONE) ISMXM(x,k) = ESL_MAX(ISMXM(x,k), ISMXM(xv,k));
    /* ID state 2 transitions */
    IDMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(TSCDELTA(k,e2HP_SD_ID) * (SDMXM(xv,k) + pp->dp[x][e2G_ID]),
		   TSCDELTA(k,e2HP_ID_ID) * (IDMXM(xv,k) + pp->dp[x][e2G_ID]));
      IDMXM(x,k) = sc;	  
    }
    if (apair == NONE) IDMXM(x,k) = ESL_MAX(IDMXM(x,k), IDMXM(xv,k));
   /* BI state 2 transitions */
    BIMXM(x,k) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(TSCDELTA(k,e2HP_BB_BI) * (BBMXM(xv,k) + pp->dp[x][e2G_BI]), 
		   TSCDELTA(k,e2HP_BI_BI) * (BIMXM(xv,k) + pp->dp[x][e2G_BI]));
      BIMXM(x,k) = sc;
    }
    if (apair == NONE) BIMXM(x,k) = ESL_MAX(BIMXM(x,k), BIMXM(xv,k));
    /* SI state 3 transitions */
    SIMXM(x,k) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(TSCDELTA(k,e2HP_SS_SI) * (SSMXM(xv,k) + pp->dp[x][e2G_SI]),
		   TSCDELTA(k,e2HP_SI_SI) * (SIMXM(xv,k) + pp->dp[x][e2G_SI]));
      SIMXM(x,k) = sc;
    }
    if (apair == NONE) SIMXM(x,k) = ESL_MAX(SIMXM(x,k), SIMXM(xv,k));
    /* DI state 2 transitions */
    DIMXM(x,k) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(TSCDELTA(k,e2HP_DS_DI) * (DSMXM(xv,k) + pp->dp[x][e2G_DI]),
		   TSCDELTA(k,e2HP_DI_DI) * (DIMXM(xv,k) + pp->dp[x][e2G_DI]));
      DIMXM(x,k) = sc;	  
    }
    if (apair == NONE) DIMXM(x,k) = ESL_MAX(DIMXM(x,k), DIMXM(xv,k));
    /* Ii state 4 transitions */
    IiMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(ESL_MAX(TSCDELTA(k,e2HP_BI_Ii) * (BIMXM(xv,k) + pp->dp[x][e2G_ii]),
			   TSCDELTA(k,e2HP_SI_Ii) * (SIMXM(xv,k) + pp->dp[x][e2G_ii])),
		   ESL_MAX(TSCDELTA(k,e2HP_Ii_Ii) * (IiMXM(xv,k) + pp->dp[x][e2G_ii]),
			   TSCDELTA(k,e2HP_iI_Ii) * (iIMXM(xv,k) + pp->dp[x][e2G_ii])));
      IiMXM(x,k) = sc;
    }
    if (apair == NONE) IiMXM(x,k) = ESL_MAX(IiMXM(x,k), IiMXM(xv,k));
    /* iI state 4 transitions */
    iIMXM(x,k) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(ESL_MAX(TSCDELTA(k,e2HP_IB_iI) * (IBMXM(xv,k) + pp->dp[x][e2G_II]),
			   TSCDELTA(k,e2HP_IS_iI) * (ISMXM(xv,k) + pp->dp[x][e2G_II])),
		   ESL_MAX(TSCDELTA(k,e2HP_Ii_iI) * (IiMXM(xv,k) + pp->dp[x][e2G_II]),
			   TSCDELTA(k,e2HP_iI_iI) * (iIMXM(xv,k) + pp->dp[x][e2G_II])));
		   iIMXM(x,k) = sc;
		   }
    if (apair == NONE) iIMXM(x,k) = ESL_MAX(iIMXM(x,k), iIMXM(xv,k));
   
    /* EE state 13 transitions */
    sc = -eslINFINITY;
    if (x == L) {
      sc = ESL_MAX(ESL_MAX(TSCDELTA(k,e2HP_BB_EE) * BBMXM(x,k), 
			   TSCDELTA(k,e2HP_IB_EE) * IBMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_SS_EE) * SSMXM(x,k),
			   TSCDELTA(k,e2HP_DS_EE) * DSMXM(x,k)));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(k,e2HP_IS_EE) * ISMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_SD_EE) * SDMXM(x,k), 
			   TSCDELTA(k,e2HP_DD_EE) * DDMXM(x,k)));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(k,e2HP_ID_EE) * IDMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_BI_EE) * BIMXM(x,k), 
			   TSCDELTA(k,e2HP_SI_EE) * SIMXM(x,k)));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(k,e2HP_DI_EE) * DIMXM(x,k)),
		   ESL_MAX(TSCDELTA(k,e2HP_Ii_EE) * IiMXM(x,k), 
			   TSCDELTA(k,e2HP_iI_EE) * iIMXM(x,k)));
    }
   
    if (sc > EE) { EE = sc; kmax = k; }
  }
  
  /* Recursion. Done as a pull.
   */
    for (x = 1; x <= L; x++) {
      xv = x - 1;
      
      isgapl = isgap(psql->prof[x], psql->abc->K);
      isgapr = isgap(psqr->prof[x], psqr->abc->K);
      if      (!isgapl && !isgapr) apair = PAIR;
      else if ( isgapl && !isgapr) apair = RRES;
      else if (!isgapl &&  isgapr) apair = LRES;
      else                         apair = NONE;
      
      for (k = 1; k <= gx->M; k++) {
	kv  = k - 1;
	kk  = k*e2G_NSCELLS;

	/* BB state  0 transitions */
	BBMXM(x,k) = -eslINFINITY;
	if (apair == NONE) BBMXM(x,k) = ESL_MAX(BBMXM(x,k), BBMXM(xv,k));
	/* SS state 13 transitions */
	SSMXM(x,k) = -eslINFINITY;
	if (apair == PAIR) {
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(kv,e2HP_BB_SS) * (BBMXM(xv,kv) + pp->dp[x][kk+e2G_SS]), 
			       TSCDELTA(kv,e2HP_IB_SS) * (IBMXM(xv,kv) + pp->dp[x][kk+e2G_SS])),
		       ESL_MAX(TSCDELTA(kv,e2HP_SS_SS) * (SSMXM(xv,kv) + pp->dp[x][kk+e2G_SS]),
			       TSCDELTA(kv,e2HP_DS_SS) * (DSMXM(xv,kv) + pp->dp[x][kk+e2G_SS])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(kv,e2HP_IS_SS) * (ISMXM(xv,kv) + pp->dp[x][kk+e2G_SS])),
		       ESL_MAX(TSCDELTA(kv,e2HP_SD_SS) * (SDMXM(xv,kv) + pp->dp[x][kk+e2G_SS]), 
			       TSCDELTA(kv,e2HP_DD_SS) * (DDMXM(xv,kv) + pp->dp[x][kk+e2G_SS])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(kv,e2HP_ID_SS) * (IDMXM(xv,kv) + pp->dp[x][kk+e2G_SS])),
		       ESL_MAX(TSCDELTA(kv,e2HP_BI_SS) * (BIMXM(xv,kv) + pp->dp[x][kk+e2G_SS]), 
			       TSCDELTA(kv,e2HP_SI_SS) * (SIMXM(xv,kv) + pp->dp[x][kk+e2G_SS])));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(kv,e2HP_DI_SS) * (DIMXM(xv,kv) + pp->dp[x][kk+e2G_SS])),
		       ESL_MAX(TSCDELTA(kv,e2HP_Ii_SS) * (IiMXM(xv,kv) + pp->dp[x][kk+e2G_SS]), 
			       TSCDELTA(kv,e2HP_iI_SS) * (iIMXM(xv,kv) + pp->dp[x][kk+e2G_SS])));
	  SSMXM(x,k) = sc;	  
	}
	if (apair == NONE) SSMXM(x,k) = ESL_MAX(SSMXM(x,k), SSMXM(xv,k));

	/* DS state 6 transitions */
	DSMXM(x,k) = -eslINFINITY;
	if (apair == RRES) {
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(kv,e2HP_SS_DS) * (SSMXM(xv,kv) + pp->dp[x][kk+e2G_DS]),
			       TSCDELTA(kv,e2HP_DS_DS) * (DSMXM(xv,kv) + pp->dp[x][kk+e2G_DS])),
		       ESL_MAX(TSCDELTA(kv,e2HP_SD_DS) * (SDMXM(xv,kv) + pp->dp[x][kk+e2G_DS]),
			       TSCDELTA(kv,e2HP_DD_DS) * (DDMXM(xv,kv) + pp->dp[x][kk+e2G_DS])));
	  sc = ESL_MAX(sc,
		       ESL_MAX(TSCDELTA(kv,e2HP_SI_DS) * (SIMXM(xv,kv) + pp->dp[x][kk+e2G_DS]),
			       TSCDELTA(kv,e2HP_DI_DS) * (DIMXM(xv,kv) + pp->dp[x][kk+e2G_DS])));
	  DSMXM(x,k) = sc;	  
	}
	if (apair == NONE) DSMXM(x,k) = ESL_MAX(DSMXM(x,k), DSMXM(xv,k));

	/* SD state 6 transitions */
	SDMXM(x,k) = -eslINFINITY;
	if (apair == LRES) {
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(kv,e2HP_SS_SD) * (SSMXM(xv,kv) + pp->dp[x][kk+e2G_SD]),
			       TSCDELTA(kv,e2HP_DS_SD) * (DSMXM(xv,kv) + pp->dp[x][kk+e2G_SD])),
		       ESL_MAX(TSCDELTA(kv,e2HP_IS_SD) * (ISMXM(xv,kv) + pp->dp[x][kk+e2G_SD]),
			       TSCDELTA(kv,e2HP_SD_SD) * (SDMXM(xv,kv) + pp->dp[x][kk+e2G_SD])));
	  sc = ESL_MAX(sc,
		       ESL_MAX(TSCDELTA(kv,e2HP_DD_SD) * (DDMXM(xv,kv) + pp->dp[x][kk+e2G_SD]),
			       TSCDELTA(kv,e2HP_ID_SD) * (IDMXM(xv,kv) + pp->dp[x][kk+e2G_SD])));
	  SDMXM(x,k) = sc;	  
	}
	if (apair == NONE) SDMXM(x,k) = ESL_MAX(SDMXM(x,k), SDMXM(xv,k));
	
	/* DD state 4 transitions */
	DDMXM(x,k) = -eslINFINITY;
#if 0
	sc = ESL_MAX(ESL_MAX(TSCDELTA(kv,e2HP_SS_DD) * (SSMXM(x,kv) + pp->dp[x][kk+e2G_DD]),
			     TSCDELTA(kv,e2HP_DS_DD) * (DSMXM(x,kv) + pp->dp[x][kk+e2G_DD])),
		     ESL_MAX(TSCDELTA(kv,e2HP_SD_DD) * (SDMXM(x,kv) + pp->dp[x][kk+e2G_DD]),
			     TSCDELTA(kv,e2HP_DD_DD) * (DDMXM(x,kv) + pp->dp[x][kk+e2G_DD])));
	DDMXM(x,k) = ESL_MAX(sc, DDMXM(xv,k));
#endif
	if (apair == NONE) {
#if 0
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(kv,e2HP_SS_DD) * (SSMXM(xv,kv) + pp->dp[x][kk+e2G_DD]),
			       TSCDELTA(kv,e2HP_DS_DD) * (DSMXM(xv,kv) + pp->dp[x][kk+e2G_DD])),
		       ESL_MAX(TSCDELTA(kv,e2HP_SD_DD) * (SDMXM(xv,kv) + pp->dp[x][kk+e2G_DD]),
			       TSCDELTA(kv,e2HP_DD_DD) * (DDMXM(xv,kv) + pp->dp[x][kk+e2G_DD])));
	  DDMXM(x,k) = ESL_MAX(sc, DDMXM(xv,k));
#endif
	  DDMXM(x,k) = ESL_MAX(DDMXM(x,k), DDMXM(xv,k));
	}
	/* IB state 2 transitions */
	IBMXM(x,k) = -eslINFINITY;
	if (apair == LRES) {
	  sc = ESL_MAX(TSCDELTA(k,e2HP_BB_IB) * (BBMXM(xv,k) + pp->dp[x][kk+e2G_IB]), 
		       TSCDELTA(k,e2HP_IB_IB) * (IBMXM(xv,k) + pp->dp[x][kk+e2G_IB]));
	  IBMXM(x,k) = sc;
	}
	if (apair == NONE) IBMXM(x,k) = ESL_MAX(IBMXM(x,k), IBMXM(xv,k));
	/* IS state 2 transitions  */
	ISMXM(x,k) = -eslINFINITY;
	if (apair == LRES) {
	  sc = ESL_MAX(TSCDELTA(k,e2HP_SS_IS) * (SSMXM(xv,k) + pp->dp[x][kk+e2G_IS]),
		       TSCDELTA(k,e2HP_IS_IS) * (ISMXM(xv,k) + pp->dp[x][kk+e2G_IS]));
	  ISMXM(x,k) = sc;
	}
	if (apair == NONE) ISMXM(x,k) = ESL_MAX(ISMXM(x,k), ISMXM(xv,k));
	/* ID state 2 transitions */
	IDMXM(x,k) = -eslINFINITY;
	if (apair == LRES) {
	  sc = ESL_MAX(TSCDELTA(k,e2HP_SD_ID) * (SDMXM(xv,k) + pp->dp[x][kk+e2G_ID]),
		       TSCDELTA(k,e2HP_ID_ID) * (IDMXM(xv,k) + pp->dp[x][kk+e2G_ID]));
	  IDMXM(x,k) = sc;	  
	}
	if (apair == NONE) IDMXM(x,k) = ESL_MAX(IDMXM(x,k), IDMXM(xv,k));
	
	/* BI state 2 transitions */
	BIMXM(x,k) = -eslINFINITY;
	if (apair == RRES) {
	  sc = ESL_MAX(TSCDELTA(k,e2HP_BB_BI) * (BBMXM(xv,k) + pp->dp[x][kk+e2G_BI]), 
		       TSCDELTA(k,e2HP_BI_BI) * (BIMXM(xv,k) + pp->dp[x][kk+e2G_BI]));
	  sc = ESL_MAX(sc, BIMXM(xv,k));
	  BIMXM(x,k) = sc;
	}
	if (apair == NONE) BIMXM(x,k) = ESL_MAX(BIMXM(x,k), BIMXM(xv,k));
	/* SI state 3 transitions */
	SIMXM(x,k) = -eslINFINITY;
	if (apair == RRES) {
	  sc = ESL_MAX(TSCDELTA(k,e2HP_SS_SI) * (SSMXM(xv,k) + pp->dp[x][kk+e2G_SI]),
		       TSCDELTA(k,e2HP_SI_SI) * (SIMXM(xv,k) + pp->dp[x][kk+e2G_SI]));
	  if (apair == NONE) sc = ESL_MAX(sc, SIMXM(xv,k));
	  SIMXM(x,k) = sc;
	}
	if (apair == NONE) SIMXM(x,k) = ESL_MAX(SIMXM(x,k), SIMXM(xv,k));
	/* DI state 2 transitions */
	DIMXM(x,k) = -eslINFINITY;
	if (apair == RRES) {
	  sc = ESL_MAX(TSCDELTA(k,e2HP_DS_DI) * (DSMXM(xv,k) + pp->dp[x][kk+e2G_DI]),
		       TSCDELTA(k,e2HP_DI_DI) * (DIMXM(xv,k) + pp->dp[x][kk+e2G_DI]));
	  if (apair == NONE) sc = ESL_MAX(sc, DIMXM(xv,k));
	  DIMXM(x,k) = sc;
	}	  
	if (apair == NONE) DIMXM(x,k) = ESL_MAX(DIMXM(x,k), DIMXM(xv,k));
	/* Ii state 4 transitions */
	IiMXM(x,k) = -eslINFINITY;
	if (apair == LRES) {
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(k,e2HP_BI_Ii) * (BIMXM(xv,k) + pp->dp[x][kk+e2G_ii]),
			       TSCDELTA(k,e2HP_SI_Ii) * (SIMXM(xv,k) + pp->dp[x][kk+e2G_ii])),
		       ESL_MAX(TSCDELTA(k,e2HP_Ii_Ii) * (IiMXM(xv,k) + pp->dp[x][kk+e2G_ii]),
			       TSCDELTA(k,e2HP_iI_Ii) * (iIMXM(xv,k) + pp->dp[x][kk+e2G_ii])));
	  IiMXM(x,k) = sc;
	}
	if (apair == NONE) IiMXM(x,k) = ESL_MAX(IiMXM(x,k), IiMXM(xv,k));
	/* iI state 4 transitions */
	iIMXM(x,k) = -eslINFINITY;
	if (apair == RRES) {
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(k,e2HP_IB_iI) * (IBMXM(xv,k) + pp->dp[x][kk+e2G_II]),
			       TSCDELTA(k,e2HP_IS_iI) * (ISMXM(xv,k) + pp->dp[x][kk+e2G_II])),
		       ESL_MAX(TSCDELTA(k,e2HP_Ii_iI) * (IiMXM(xv,k) + pp->dp[x][kk+e2G_II]),
			       TSCDELTA(k,e2HP_iI_iI) * (iIMXM(xv,k) + pp->dp[x][kk+e2G_II])));
	  iIMXM(x,k) = sc;
	}
	if (apair == NONE) iIMXM(x,k) = ESL_MAX(iIMXM(x,k), iIMXM(xv,k));
	
	/* EE state 13 transitions */
	if (x == L) {
	  sc = ESL_MAX(ESL_MAX(TSCDELTA(k,e2HP_BB_EE) * BBMXM(x,k), 
			       TSCDELTA(k,e2HP_IB_EE) * IBMXM(x,k)),
		       ESL_MAX(TSCDELTA(k,e2HP_SS_EE) * SSMXM(x,k),
			       TSCDELTA(k,e2HP_DS_EE) * DSMXM(x,k)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(k,e2HP_IS_EE) * ISMXM(x,k)),
		       ESL_MAX(TSCDELTA(k,e2HP_SD_EE) * SDMXM(x,k), 
			       TSCDELTA(k,e2HP_DD_EE) * DDMXM(x,k)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(k,e2HP_ID_EE) * IDMXM(x,k)),
		       ESL_MAX(TSCDELTA(k,e2HP_BI_EE) * BIMXM(x,k), 
			       TSCDELTA(k,e2HP_SI_EE) * SIMXM(x,k)));
	  sc = ESL_MAX(ESL_MAX(sc, 
			       TSCDELTA(k,e2HP_DI_EE) * DIMXM(x,k)),
		       ESL_MAX(TSCDELTA(k,e2HP_Ii_EE) * IiMXM(x,k), 
			       TSCDELTA(k,e2HP_iI_EE) * iIMXM(x,k)));
	  
	  if (sc > EE) { EE = sc; kmax = k; }
	}
#if 0
	if (x==27||x==26) 
	  printf("e2fhmmer_OA x %d  k %d pair %d BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f iI %f EE %f\n", 
				 x, k, apair,
			 BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), 
			 SDMXM(x,k), DDMXM(x,k), IDMXM(x,k), 
			 BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k), EE);
#endif
      }
    }
    
    if (ret_e != NULL) *ret_e = EE;
    if (ret_k != NULL) *ret_k = kmax;
    
    return eslOK;
}
/*---------------------- end, oa fill ---------------------------*/

static inline int   apair(int x, PSQ *sql, PSQ *sqr);
static inline void  ancestral_res_isgap(int y, PSQ *sqa);
static inline float get_postprob(const E2_GMX *pp, int scur, int sprv, int k, int x);
static inline int   find_argmax(ESL_RANDOMNESS *r, const float *v, int n);

static inline int   select_ss_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
						int *ret_k, int *ret_x, int *ret_i, int *ret_j, 
						P7_RATE *R, int dtl, int dtr, PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa);
static inline int   select_ds_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
						int *ret_k, int *ret_x, int *ret_i, int *ret_j,
						P7_RATE *R, int dtl, int dtr, PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa);
static inline int   select_sd_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
						int *ret_k, int *ret_x, int *ret_i, int *ret_j, 
						P7_RATE *R, int dtl, int dtr, PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa);
static inline int   select_dd_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
						int *ret_k, int *ret_x, int *ret_i, int *ret_j,  
						P7_RATE *R,                   PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa);

static inline int   select_ib(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair);
static inline int   select_is(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair);
static inline int   select_id(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair);
static inline int   select_bi(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j, int L, int apair);
static inline int   select_si(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j, int L, int apair);
static inline int   select_di(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j, int L, int apair);
static inline int   select_il(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair);
static inline int   select_ir(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j, int L, int apair);
static inline int   select_ee(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int      x, int L);
static inline int   select_xx(int sprvprv, int *ret_x);

/*****************************************************************
 * 2. Optimal alignment accuracy, traceback
 *****************************************************************/

/* Function:  e2fhmmer_GOATrace()
 * Synopsis:  Optimal accuracy decoding: traceback.
 * Incept:    ER, Thu Oct 10 09:40:51 EDT 2013
 *
 * Purpose:   The traceback stage of the optimal accuracy decoding algorithm
 *            \citep{Kall05}.
 *            
 *            Caller provides the OA DP matrix <gx> that was just
 *            calculated by <e2_GOptimalAccuracy()>, as well as the
 *            posterior decoding matrix <pp>, which was calculated by
 *            Forward/Backward on a target sequence of length <L>
 *            using the query model <gm>.
 *            
 *            Caller provides an empty traceback structure <tr> to
 *            hold the result, allocated to hold optional posterior
 *            probability annotation on residues (with
 *            <e2_trace_CreateWithPP()>, generally).  This will be
 *            internally reallocated as needed for larger traces.
 *
 *            We also construct the profile for the ancestral sequence,
 *            by adding the counts of the descendant. Profiles are
 *            still in counts, not normalized yet.
 *
 * Args:      gm    - query profile      
 *            pp    - posterior decoding matrix created by <e2_PosteriorDecoding()>
 *            gx    - OA DP matrix calculated by  <e2_OptimalAccuracyDP()>
 *            tr    - RESULT: OA traceback, allocated with posterior probs
 *
 * Returns:   <eslOK> on success, and <tr> contains the OA traceback.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
e2fhmmer_GOATrace(ESL_RANDOMNESS *r, P7_RATE *R7, float timel, float timer, 
		  P7_BG *bg, const E2HMMER_PROFILE *gm, int kmax,
		  const E2_GMX *pp, const E2_GMX *gx, E2_TRACE *tr, PSQ *psql, PSQ *psqr, PSQ **ret_psqA)
{
  PSQ       *psqA   = NULL;
  int        L;
  int        x;                 /* position in alingment   */
  int        k;                 /* position in model */
  int        a      = 1;        /* position in the ancestral sequence    */
  int        i, j;
  float      postprob;
  int        sprv, scur;
  int        sprvprv;
  int        dtl, dtr;          /* discretized times */
  int        apairprv, apaircur;
  int        status;

#ifdef e2_DEBUGGING
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  psqA = psq_Create(gm->abc);

  if (psql->n != psqr->n) { printf("sqs are not aligned\n"); return eslFAIL; }

  /* start */
  L = x = psql->n;
  psql->n = i = psq_NResidues(psql);
  psqr->n = j = psq_NResidues(psqr);

  k        = kmax;
  apairprv = NONE;
  if ((status = e2_trace_AppendWithPP(tr, e2T_EE, k, a, i, j, 0.0)) != eslOK) return status;
  printf("EE(kmax %d x %d i %d j %d) = %f\n", kmax, x, i, j, gx->dp[x][k* e2G_NSCELLS + e2G_EE]);
   
  dtl = select_discrete_time(R7, timel); if (dtl < 0) goto ERROR;
  dtr = select_discrete_time(R7, timer); if (dtr < 0) goto ERROR;
 
  sprvprv = sprv = e2T_EE;
  while (sprv != e2T_BB) 
    {      
#if 0
      printf("%s  k %d x %d L %d apair %d i %d j %d a %d N %d || ",
	     e2_model_DecodeStatetype(sprv), k, x, L, apairprv, i, j, a, tr->N);
#endif

      switch (sprv) {
      case e2T_SS: scur = select_ss_and_ancestral_res(r, gm, gx, &k, &x, &i, &j, R7, dtl, dtr, psql, psqr, apairprv, &a, psqA); break;
      case e2T_DS: scur = select_ds_and_ancestral_res(r, gm, gx, &k, &x, &i, &j, R7, dtl, dtr, psql, psqr, apairprv, &a, psqA); break;
      case e2T_SD: scur = select_sd_and_ancestral_res(r, gm, gx, &k, &x, &i, &j, R7, dtl, dtr, psql, psqr, apairprv, &a, psqA); break;
      case e2T_DD: scur = select_dd_and_ancestral_res(r, gm, gx, &k, &x, &i, &j, R7,           psql, psqr, apairprv, &a, psqA); break;

      case e2T_IB: scur = select_ib(r, gm, gx, k, &x, &i,     L, apairprv); ancestral_res_isgap (                     a, psqA); break;
      case e2T_IS: scur = select_is(r, gm, gx, k, &x, &i,     L, apairprv); ancestral_res_isgap (                     a, psqA); break;
      case e2T_ID: scur = select_id(r, gm, gx, k, &x, &i,     L, apairprv); ancestral_res_isgap (                     a, psqA); break;

      case e2T_BI: scur = select_bi(r, gm, gx, k, &x,     &j, L, apairprv); ancestral_res_isgap (                     a, psqA); break;
      case e2T_SI: scur = select_si(r, gm, gx, k, &x,     &j, L, apairprv); ancestral_res_isgap (                     a, psqA); break;
      case e2T_DI: scur = select_di(r, gm, gx, k, &x,     &j, L, apairprv); ancestral_res_isgap (                     a, psqA); break;
      case e2T_ii: scur = select_il(r, gm, gx, k, &x, &i,     L, apairprv); ancestral_res_isgap (                     a, psqA); break;
      case e2T_II: scur = select_ir(r, gm, gx, k, &x,     &j, L, apairprv); ancestral_res_isgap (                     a, psqA); break;
 
      case e2T_XX: scur = select_xx(sprvprv,      &x);                                                                          break;
      case e2T_EE: scur = select_ee(r, gm, gx, k, x, L);                                                                        break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      apaircur = apair(x, psql, psqr);
      postprob = get_postprob(pp, scur, sprv, k, x);

     if ((scur != e2T_XX && apaircur != NONE) || scur ==e2T_BB) {
	if ((status = e2_trace_AppendWithPP(tr, scur, k, a, i, j, postprob)) != eslOK) return status;      
      }
      
 #if 0
      printf("%s  apair %d i %d j %d a %d N %d\n",
	     e2_model_DecodeStatetype(scur), apaircur, i, j, a, tr->N);
#endif

      psq_Grow(psqA, NULL);
      apairprv = apaircur;
      sprvprv = sprv;
      sprv    = scur;
    }
 
  if (i != 0 || j != 0) { printf("e2fhmmer_GOATrace() failed. i=%d j=%d\n", i, j); exit(1); }

  tr->M     = a; /* length of ancestral sequence+1  */
  tr->Lrow  = psql->n;
  tr->Lcol  = psqr->n;
  tr->rowsq = e2P_SL;
  if ((status = e2_trace_Reverse(tr)) != eslOK) goto ERROR; 
  if ((status = psq_Reverse(psqA))    != eslOK) goto ERROR;

  *ret_psqA = psqA;
  return eslOK;

 ERROR:  
  if (psqA != NULL) psq_Destroy(psqA);
  return status;
}

static inline int
apair(int x, PSQ *sql, PSQ *sqr)
{
  int apair;
  int isgapl;
  int isgapr;
  
  isgapl = isgap(sql->prof[x], sql->abc->K);
  isgapr = isgap(sqr->prof[x], sqr->abc->K);
  if      (!isgapl && !isgapr) apair = PAIR;
  else if ( isgapl && !isgapr) apair = RRES;
  else if (!isgapl &&  isgapr) apair = LRES;
  else                         apair = NONE;
  
  return apair;
}

static inline float
get_postprob(const E2_GMX *pp, int scur, int sprv, int k, int x)
{
  float **dp  = pp->dp;

  switch (scur) {
  case e2T_BB: return BBMXM(x,k);
  case e2T_IB: return IBMXM(x,k);
  case e2T_SS: return SSMXM(x,k);
  case e2T_DS: return DSMXM(x,k);
  case e2T_IS: return ISMXM(x,k);
  case e2T_SD: return SDMXM(x,k);
  case e2T_DD: return DDMXM(x,k);
  case e2T_ID: return IDMXM(x,k);
  case e2T_BI: return BIMXM(x,k);
  case e2T_SI: return SIMXM(x,k);
  case e2T_DI: return DIMXM(x,k);
  case e2T_II: return iIMXM(x,k);
  case e2T_ii: return IiMXM(x,k);
  default:     return 0.0;
  }
}

static inline void
ancestral_res_isgap(int y, PSQ *sqa)
{
  int K = sqa->abc->K;

  esl_vec_FSet(sqa->prof[y], K+1, -eslINFINITY); 
  sqa->prof[y][K] = 0.;
}

/* only difference with esl_vec_FArgmax()
 * is that if the arg max is dergenerate,
 * it picks one randomly as suposed to the 
 * first one in the array */
static inline int
find_argmax(ESL_RANDOMNESS *r, const float *v, int n) 
{
  ESL_STACK *alt = NULL;
  float      max;
  int        argmax;
  int        nequiv;
  int        i;
  int        x;

  alt = esl_stack_ICreate();
 
  max = esl_vec_FMax(v, n);
  for (i = 0; i < n; i ++) {
    if (v[i] >= max)  esl_stack_IPush(alt, i);
  }

  /* Choose one of the alternatives at random */
  nequiv = esl_stack_ObjectCount(alt);  /* how many solutions? */
  x = esl_rnd_Roll(r, nequiv);          /* uniformly, 0.nequiv-1 */
  esl_stack_DiscardTopN(alt, x);        /* dig down to choice */
  esl_stack_IPop(alt, &argmax);         /* pop it off */
    
  esl_stack_Destroy(alt);
  return argmax;
}

static inline int
select_ss_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
			    int *ret_k, int *ret_x, int *ret_i, int *ret_j,
			    P7_RATE *R, int dtl, int dtr, PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA(k,) macro works */
  float        path[13];
  int          state[13] = { e2T_BB, e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_ii, e2T_II  };
  int          k = *ret_k;
  int          x = *ret_x;
  int          y = *ret_y;
  int          i = *ret_i;
  int          j = *ret_j;
  float       *pl = sql->prof[x];
  float       *pr = sqr->prof[x];
  int          kv;
  int          xv;
  int          stidx;
  int          K = sql->abc->K;
  int          a, b, c;
 
  kv = k - 1;
  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_BB_SS) * BBMXM(xv,kv);
  path[1]  = TSCDELTA(k,e2HP_IB_SS) * IBMXM(xv,kv);
  path[2]  = TSCDELTA(k,e2HP_SS_SS) * SSMXM(xv,kv);
  path[3]  = TSCDELTA(k,e2HP_DS_SS) * DSMXM(xv,kv);
  path[4]  = TSCDELTA(k,e2HP_IS_SS) * ISMXM(xv,kv);
  path[5]  = TSCDELTA(k,e2HP_SD_SS) * SDMXM(xv,kv);
  path[6]  = TSCDELTA(k,e2HP_DD_SS) * DDMXM(xv,kv);
  path[7]  = TSCDELTA(k,e2HP_ID_SS) * IDMXM(xv,kv);
  path[8]  = TSCDELTA(k,e2HP_BI_SS) * BIMXM(xv,kv);
  path[9]  = TSCDELTA(k,e2HP_SI_SS) * SIMXM(xv,kv);
  path[10] = TSCDELTA(k,e2HP_DI_SS) * DIMXM(xv,kv);
  path[11] = TSCDELTA(k,e2HP_Ii_SS) * IiMXM(xv,kv);
  path[12] = TSCDELTA(k,e2HP_iI_SS) * iIMXM(xv,kv);
  stidx = find_argmax(r, path, 13);
 
  if (apair == NONE) {
    *ret_k = k;
    *ret_x = x;
    *ret_i = i;
    *ret_j = j;
    *ret_y = y;
  }
  else {
    if (apair != PAIR) { printf("ancestral_res_fromSS() you should not be here. x=%d\n", x); exit(1); }
    esl_vec_FSet(sqa->prof[y], K+1, -eslINFINITY); 
    for (a = 0;  a < K; a++) 
      {
	for (b = 0; b < K; b++) 
	  for (c = 0; c < K; c++) 
	    {	  
	      sqa->prof[y][a] = 
		p7_FLogsum(sqa->prof[y][a], log(R->pzero[k][a]) +
			   log(R->Pdt[dtl][k]->mx[a][b]) + 
			   log(R->Pdt[dtr][k]->mx[a][c]) + pl[b] + pr[c]);
#if 0
	      if (k==67&&b==11&&c==11) printf("^^ k %d a %d b %d c %d zero %f sc %f\n", 
					      k, a, b, c, R->pzero[k][a], sqa->prof[y][a]);
#endif
	    }
     }
#if 0
    if (k==67) printf("^^ k %d dtl %d dtr %d b %d c %d sqa %d root %d \n", 
		      k, dtl, dtr,
		      esl_vec_FArgMax(pl,K), 
		      esl_vec_FArgMax(pr,K), 
		      esl_vec_FArgMax(sqa->prof[y],K), 
		      esl_vec_FArgMax(R->pzero[k], K));
#endif
    esl_vec_FLogNorm(sqa->prof[y], K+1); 
    esl_vec_FLog    (sqa->prof[y], K+1); 
    
    *ret_k = kv;
    *ret_x = xv;
    *ret_i = i - 1;
    *ret_j = j - 1;
    *ret_y = y + 1;
    sqa->n ++;
  }
  
  return (apair == NONE)? e2T_XX : state[stidx];
}


static inline int
select_ds_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
			    int *ret_k, int *ret_x, int *ret_i, int *ret_j,
			    P7_RATE *R, int dtl, int dtr, PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA(k,) macro works */
  float        path[6];
  int          state[6] = { e2T_SS, e2T_DS, e2T_SD, e2T_DD, e2T_SI, e2T_DI };
  int          k = *ret_k;
  int          x = *ret_x;
  int          y = *ret_y;
  int          i = *ret_i;
  int          j = *ret_j;
  int          kv;
  int          xv;
  int          stidx;
  int          K = sqa->abc->K;
  int          a, b;
  
  kv = k - 1;
  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_SS_DS) * SSMXM(xv,kv);
  path[1]  = TSCDELTA(k,e2HP_DS_DS) * DSMXM(xv,kv);
  path[2]  = TSCDELTA(k,e2HP_SD_DS) * SDMXM(xv,kv);
  path[3]  = TSCDELTA(k,e2HP_DD_DS) * DDMXM(xv,kv);
  path[4]  = TSCDELTA(k,e2HP_SI_DS) * SIMXM(xv,kv);
  path[5]  = TSCDELTA(k,e2HP_DI_DS) * DIMXM(xv,kv);
  stidx = find_argmax(r, path, 6);

  if (apair == NONE) {
    *ret_k = k;
    *ret_x = x;
    *ret_i = i;
    *ret_j = j;
    *ret_y = y;
  }
  else {
   if (apair != RRES) { printf("ancestral_res_fromDS() you should not be here. x=%d\n", x); exit(1); }
    esl_vec_FSet(sqa->prof[y], K+1, -eslINFINITY); 
     
    for (a = 0;  a < K; a++) 
      for (b = 0; b < K; b++) 
	{	  
	  sqa->prof[y][a] =
	    p7_FLogsum(sqa->prof[y][a], log(R->pzero[k][a]) + log(R->Pdt[dtr][k]->mx[a][b])+ sqr->prof[x][b]);
	}
    
    esl_vec_FLogNorm(sqa->prof[y], K+1); 
    esl_vec_FLog    (sqa->prof[y], K+1); 
 
    *ret_k = kv;
    *ret_x = xv;
    *ret_i = i;
    *ret_j = j - 1;
    *ret_y = y + 1;
    sqa->n ++;
  }

  return (apair == NONE)? e2T_XX : state[stidx];
}


static inline int
select_sd_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
			    int *ret_k, int *ret_x, int *ret_i, int *ret_j,
			    P7_RATE *R, int dtl, int dtr, PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  =(float const **) gm->tsc;	/* so TSCDELTA(k,) macro works */
  float        path[6];
  int          state[6] = { e2T_SS, e2T_DS, e2T_IS, e2T_SD, e2T_DD, e2T_ID };
  int          k = *ret_k;
  int          x = *ret_x;
  int          y = *ret_y;
  int          i = *ret_i;
  int          j = *ret_j;
  int          kv;
  int          xv;
  int          stidx;
  int          K = sqa->abc->K;
  int          a, b;

  kv = k - 1;
  xv = x - 1;
  
  path[0]  = TSCDELTA(k,e2HP_SS_SD) * SSMXM(xv,kv);
  path[1]  = TSCDELTA(k,e2HP_DS_SD) * DSMXM(xv,kv);
  path[2]  = TSCDELTA(k,e2HP_IS_SD) * ISMXM(xv,kv);
  path[3]  = TSCDELTA(k,e2HP_SD_SD) * SDMXM(xv,kv);
  path[4]  = TSCDELTA(k,e2HP_DD_SD) * DDMXM(xv,kv);
  path[5]  = TSCDELTA(k,e2HP_ID_SD) * IDMXM(xv,kv);
  stidx = find_argmax(r, path, 6);

  if (apair == NONE) {
    *ret_k = k;
    *ret_x = x;
    *ret_i = i;
    *ret_j = j;
    *ret_y = y;
  }
  else {
   if (apair != LRES) { printf("ancestral_res_fromSD() you should not be here. x=%d\n", x); exit(1); }

    esl_vec_FSet(sqa->prof[y], K+1, -eslINFINITY); 
    for (a = 0;  a < K; a++) 
      for (b = 0; b < K; b++) 
	{	  
	  sqa->prof[y][a] = 
	    p7_FLogsum(sqa->prof[y][a], log(R->pzero[k][a]) + log(R->Pdt[dtl][k]->mx[a][b]) + sql->prof[x][b]);
	}
    esl_vec_FLogNorm(sqa->prof[y], K+1); 
    esl_vec_FLog    (sqa->prof[y], K+1); 

    *ret_k = kv;
    *ret_x = xv;
    *ret_i = i - 1;
    *ret_j = j;
    *ret_y = y + 1;
    sqa->n ++;
  }
  
  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_dd_and_ancestral_res(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, 
			    int *ret_k, int *ret_x, int *ret_i, int *ret_j,
			    P7_RATE *R, PSQ *sql, PSQ *sqr, int apair, int *ret_y, PSQ *sqa)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA() macro works */
  float        path[4];
  int          state[4] = { e2T_SS, e2T_DS, e2T_SD, e2T_DD };
  int          k = *ret_k;
  int          x = *ret_x;
  int          y = *ret_y;
  int          i = *ret_i;
  int          j = *ret_j;
  int          kv;
  int          xv;
  int          stidx;
  int          K = sqa->abc->K;
  int          a;

  kv = k - 1;
  xv = x - 1;
 
  if (apair == NONE) {
    path[0]  = TSCDELTA(k,e2HP_SS_DD) * SSMXM(xv,kv);
    path[1]  = TSCDELTA(k,e2HP_DS_DD) * DSMXM(xv,kv);
    path[2]  = TSCDELTA(k,e2HP_SD_DD) * SDMXM(xv,kv);
    path[3]  = TSCDELTA(k,e2HP_DD_DD) * DDMXM(xv,kv);
    stidx = find_argmax(r, path, 4);

   /* add residues to the ancestral sequence */
    esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
    
    for (a = 0; a < K; a++) 
      sqa->prof[y][a] = log(R->pzero[k][a]);
    esl_vec_FLogNorm(sqa->prof[y], K+1); 
    esl_vec_FLog    (sqa->prof[y], K+1); 

    *ret_k = kv;
    *ret_x = xv;
    *ret_i = i;
    *ret_j = j;
    *ret_y = y + 1;
    sqa->n ++;
  }

  return state[stidx];
}


static inline int
select_ib(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work  */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA(k,) macro works */
  float        path[2];
  int          state[2] = { e2T_BB, e2T_IB};
  int          x = *ret_x;
  int          i = *ret_i;
  int          xv;
  int          stidx;
  
  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_BB_IB) * BBMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_IB_IB) * IBMXM(xv,k);
  stidx = find_argmax(r, path, 2);

  *ret_x = (apair == NONE)? x : xv;
  *ret_i = (apair == NONE)? i : i - 1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_is(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          state[2] = { e2T_SS, e2T_IS};
  int          x = *ret_x;
  int          i = *ret_i;
  int          xv;
  int          stidx;

  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_SS_IS) * SSMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_IS_IS) * ISMXM(xv,k);
  stidx = find_argmax(r, path, 2);

  *ret_x = (apair == NONE)? x : xv;
  *ret_i = (apair == NONE)? i : i - 1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_id(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          state[2] = { e2T_SD, e2T_ID };
  int          x = *ret_x;
  int          i = *ret_i;
  int          xv;
  int          stidx;
  
  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_SD_ID) * SDMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_ID_ID) * IDMXM(xv,k);
  stidx = find_argmax(r, path, 2);

  *ret_x = (apair == NONE)? x : xv;
  *ret_i = (apair == NONE)? i : i - 1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_bi(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j,int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          state[2] = { e2T_BB, e2T_BI };
  int          x = *ret_x;
  int          j = *ret_j;
  int          xv;
  int          stidx;

  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_BB_BI) * BBMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_BI_BI) * BIMXM(xv,k);
  stidx = find_argmax(r, path, 2);

  *ret_x = (apair == NONE)? x : xv;
  *ret_j = (apair == NONE)? j : j - 1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_si(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j, int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA(k,) macro works */
  float        path[2];
  int          state[2] = { e2T_SS, e2T_SI };
  int          x = *ret_x;
  int          j = *ret_j;
  int          xv;
  int          stidx;

  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_SS_SI) * SSMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_SI_SI) * SIMXM(xv,k);
  stidx = find_argmax(r, path, 2);

  *ret_x = (apair == NONE)? x : xv;
  *ret_j = (apair == NONE)? j : j - 1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_di(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j, int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          state[2] = { e2T_DS, e2T_DI };
  int          x = *ret_x;
  int          j = *ret_j;
  int          xv;
  int          stidx;

  xv = x - 1;
 
  path[0]  = TSCDELTA(k,e2HP_DS_DI) * DSMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_DI_DI) * DIMXM(xv,k);
  stidx = find_argmax(r, path, 2);

  *ret_x = (apair == NONE)? x : xv;
  *ret_j = (apair == NONE)? j : j - 1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_il(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_i, int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work  */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA(k,) macro works */
  float        path[4];
  int          state[4] = { e2T_BI, e2T_SI, e2T_ii, e2T_II };
  int          x = *ret_x;
  int          i = *ret_i;
  int          xv;
  int          stidx;

  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_BI_Ii) * BIMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_SI_Ii) * SIMXM(xv,k);
  path[2]  = TSCDELTA(k,e2HP_Ii_Ii) * IiMXM(xv,k);
  path[3]  = TSCDELTA(k,e2HP_iI_Ii) * iIMXM(xv,k);
  stidx = find_argmax(r, path, 4);

  *ret_x = (apair == NONE)? x : xv;
  *ret_i = (apair == NONE)? i : i-1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_ir(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int *ret_x, int *ret_j, int L, int apair)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work  */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA(k,) macro works */
  float        path[4];
  int          state[4] = { e2T_IB, e2T_IS, e2T_ii, e2T_II };
  int          x = *ret_x;
  int          j = *ret_j;
  int          xv;
  int          stidx;

  xv = x - 1;

  path[0]  = TSCDELTA(k,e2HP_IB_iI) * IBMXM(xv,k);
  path[1]  = TSCDELTA(k,e2HP_IS_iI) * ISMXM(xv,k);
  path[2]  = TSCDELTA(k,e2HP_Ii_iI) * IiMXM(xv,k);
  path[3]  = TSCDELTA(k,e2HP_iI_iI) * iIMXM(xv,k);
  stidx = find_argmax(r, path, 4);

  *ret_x = (apair == NONE)? x : xv;
  *ret_j = (apair == NONE)? j : j-1;

  return (apair == NONE)? e2T_XX : state[stidx];
}

static inline int
select_xx(int sprvprv, int *ret_x)
{
  int          x = *ret_x;
  int          xv;

  xv = x - 1;
  *ret_x = xv;
  return sprvprv;
}

static inline int
select_ee(ESL_RANDOMNESS *r, const E2HMMER_PROFILE *gm, const E2_GMX *gx, int k, int x, int L)
{
  float       **dp   = gx->dp;	/* so {MDI}MXM() macros work       */
  float const **tsc  = (float const **)gm->tsc;	/* so TSCDELTA() macro works */
  float        path[13];
  int          state[13] = { e2T_BB, e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, e2T_DD, e2T_ID, e2T_BI,  e2T_SI, e2T_DI, e2T_ii, e2T_II};

  path[0]  = TSCDELTA(k,e2HP_BB_EE) * BBMXM(x,k);
  path[1]  = TSCDELTA(k,e2HP_IB_EE) * IBMXM(x,k);
  path[2]  = TSCDELTA(k,e2HP_SS_EE) * SSMXM(x,k);
  path[3]  = TSCDELTA(k,e2HP_DS_EE) * DSMXM(x,k);
  path[4]  = TSCDELTA(k,e2HP_IS_EE) * ISMXM(x,k);
  path[5]  = TSCDELTA(k,e2HP_SD_EE) * SDMXM(x,k);
  path[6]  = TSCDELTA(k,e2HP_DD_EE) * DDMXM(x,k);
  path[7]  = TSCDELTA(k,e2HP_ID_EE) * IDMXM(x,k);
  path[8]  = TSCDELTA(k,e2HP_BI_EE) * BIMXM(x,k);
  path[9]  = TSCDELTA(k,e2HP_SI_EE) * SIMXM(x,k);
  path[10] = TSCDELTA(k,e2HP_DI_EE) * DIMXM(x,k);
  path[11] = TSCDELTA(k,e2HP_Ii_EE) * IiMXM(x,k);
  path[12] = TSCDELTA(k,e2HP_iI_EE) * iIMXM(x,k);

  return state[find_argmax(r, path, 13)];
}

static int 
isgap(float *p, int K) 
{
  if (p[K] == 0.) return 1;
  else return 0;
}


/*---- internal functions ---*/

/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef E2OA_TESTDRIVE
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "msatree.h"

static int run_e2f(const PSQ *sql, const PSQ *sqr, PSQ **ret_sqa, E1_RATE *R, P7_BG *bg, float time1, float time2, double tol, int verbose);

static int
utest_e2fhmmer_OptimalAccuracy(ESL_MSA *msa, ESL_TREE *T, char *errbuf, int verbose)
{
  E1_RATE     *R = NULL;                                    /* assumes a unique set of rates for the whole process */
  P7_BG       *bg = NULL;
  ESL_STACK   *vs = NULL;	                            /* node index stack */
  ESL_SQ      *sq = NULL;
  PSQ        **sqn = NULL;                                  /* a profile sequence for each internal node */
  PSQ         *sql = NULL;                                  /* convenience pointer to sq in left branch */
  PSQ         *sqr = NULL;                                  /* convenience pointer to sq in left branch */
  float        tol = 0.0001;
  int          v;
  int          which;
  int          status;

  /* Associate tree leaves to msa sequences */
  if ((status = Tree_ReorderTaxaAccordingMSA(msa, T, errbuf, verbose)) != eslOK) goto ERROR;

  /* allocate the profile sequence for internal nodes */
  ESL_ALLOC(sqn, sizeof(PSQ *) * (T->N-1));
  for (v = 0; v < T->N-1; v ++) 
    sqn[v] = NULL;
 
  /* the background model (ancestral sequence) */
  bg = e1_bg_Create(msa->abc);
  if (bg == NULL) { status = eslEMEM; goto ERROR; }

  /* the evolutionary model (Rates) */
  R = e1_rate_Create(msa->abc);
  if (bg == NULL) { status = eslEMEM; goto ERROR; }
 
  /* some adhoc values */
  R->muD[e1R_B] = R->muD[e1R_S] = R->muD[e1R_D] = R->muD[e1R_I] = 0.1;
  R->muE[e1R_B] = R->muE[e1R_S] = R->muE[e1R_D] = R->muE[e1R_I] = 0.05;
  R->lbd[e1R_B] = R->lbd[e1R_S] = R->lbd[e1R_D] = R->lbd[e1R_I] = 0.05;
  R->xiz = R->rz = 0.1;
  ratematrix_emrate_LoadRate(R->em, "BLOSUM62", errbuf, verbose);
  
 /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslEMEM; goto ERROR; }
  if (esl_stack_IPush(vs, T->N-2) != eslOK) { status = eslEMEM; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      if (T->left[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->left[v];	
	esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	sql = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	esl_sq_Destroy(sq);
       }
      else sql = sqn[T->left[v]];

     if (T->right[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->right[v];
	esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	sqr = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	esl_sq_Destroy(sq);
      }
      else sqr = sqn[T->right[v]];
      
      if (sql != NULL && sqr != NULL) { /* ready to go: find ancestral profile sq running the e2 algorithm */
	if (verbose) printf("\nNODE %d parent %d | %d (%f,len=%d) %d (%f,len=%d)\n", v, T->parent[v], T->left[v], T->ld[v], (int)sql->n, T->right[v], T->rd[v], (int)sqr->n);
	if ((status = run_e2f(sql, sqr, &(sqn[v]), R, bg, T->ld[v], T->rd[v], tol, verbose)) != eslOK) { status = eslEMEM; goto ERROR; }
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslEMEM; goto ERROR; }; /* push parent into stack unless already at the root */
      }
      else if (sql == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])   != eslOK) { status = eslEMEM; goto ERROR; };
      }
      else if (sqr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v])  != eslOK) { status = eslEMEM; goto ERROR; };
      }
   }
  
  for (v = 0; v < T->N-1; v ++) psq_Destroy(sqn[v]); free(sqn);
  e1_bg_Destroy(bg);
  e1_rate_Destroy(R);
  esl_stack_Destroy(vs);
  return eslOK;

 ERROR:
  if (bg) e1_bg_Destroy(bg);
  if (R)  e1_rate_Destroy(R);
  return status;
}

static int
run_e2f(const PSQ *sql, const PSQ *sqr, PSQ **ret_sqa, E1_RATE *R, P7_BG *bg, float time1, float time2, double tol, int verbose)
{
  P7_HMM           *evol = NULL;
  P7_HMM           *evor = NULL;
  E2_TRACE         *tr   = NULL;
  E2_GMX           *gx1  = NULL;
  E2_GMX           *gx2  = NULL;
  E2HMMER_PROFILE  *gm   = NULL;
  ESL_ALPHABET     *abc = (ESL_ALPHABET *)sql->abc;
  PSQ              *sqa = psq_Create(abc);
  float             fsc, bsc;
  float             accscore;
  float             L;         /* average length of the two sequences */
  int               status;

  if (sql->abc->type != sqr->abc->type) {status = eslFAIL; goto ERROR; }

  /* Allocations */
  gx1 = e2hmmer_gmx_Create(sql->n, sqr->n);
  gx2 = e2hmmer_gmx_Create(sql->n, sqr->n);
  tr  = e2_trace_CreateWithPP();
  p7_FLogsumInit();

   /* Evolve the models */
  evol = e1_model_Create(R, time1, bg->f, abc, tol); if (evol == NULL) {status = eslFAIL; goto ERROR; }
  evor = e1_model_Create(R, time2, bg->f, abc, tol); if (evor == NULL) {status = eslFAIL; goto ERROR; }

  /* Configure a profile */
  L = 0.5 * (sql->n + sqr->n); 
  gm = e2_profile_Create(abc);
  e2_ProfileConfig(evol, evor, bg, gm, L, e2_GLOBAL);

  /* set length for background model */
  e1_bg_SetLength(bg, L);

  /* Run Forward, Backward; do OA fill and trace */
  e2fhmmer_GForward (sql, sqr, gm, gx1, &fsc);
  e2fhmmer_GBackward(sql, sqr, gm, gx2, &bsc);
  if (verbose) {
    printf("e2fhmmer_fwd = %.4f nats\n", fsc);
    printf("e2fhmmer_bck = %.4f nats\n", bsc);
  }
  if (fabs(bsc-fsc) > 0.01*L) { status = eslFAIL; goto ERROR; }

  e2fhmmer_GDecoding(gm, gx1, gx2, gx2);                   /* <gx2> is now the posterior decoding matrix */
  e2fhmmer_GOptimalAccuracy(gm, gx2, gx1, &accscore);      /* <gx1> is now the OA matrix */

  e2fhmmer_GOATrace(evol, evor, bg, gm, gx2, gx1, tr, sql, sqr, &sqa);  
  if (verbose) {
    printf("ancestral sequence length %d\n", (int)sqa->n);
    printf("acc = %.4f (%.2f%%)\n", accscore, (sqa->n > 0)? accscore * 100. / (float) sqa->n : accscore * 100);
    e2_trace_Dump(stdout, tr, gm, sql, sqr, sqa);
  }
 
  e1_model_Destroy(evol);
  e1_model_Destroy(evor);
  e2hmmer_gmx_Destroy(gx1);
  e2hmmer_gmx_Destroy(gx2);
  e2_profile_Destroy(gm);
  e2_trace_Destroy(tr);
 
  *ret_sqa = sqa;
  return eslOK;

 ERROR:
  if (evol) e1_model_Destroy(evol);
  if (evor) e1_model_Destroy(evor);
  if (gx1)  e2hmmer_gmx_Destroy(gx1);
  if (gx2)  e2hmmer_gmx_Destroy(gx2);
  if (gm)   e2_profile_Destroy(gm);
  if (tr)   e2_trace_Destroy(tr);
  return status;
}
#endif /* E2OA_TESTDRIVE */

/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef E2FOA_TESTDRIVE
/* gcc -std=gnu99 -Isrc -Ilib/hmmer4/src -Ilib/hmmer4/src/dp_reference -Ilib/hmmer4/src/dp_vector -Ilib/hmmer4/lib/easel -g -Wall -DE2FOA_TESTDRIVE -Isrc -Ilib/hmmer4/src -Ilib/hmmer4/src/dp_reference -Ilib/hmmer4/src/dp_vector -Ilib/hmmer4/lib/easel  -Lsrc -Llib/hmmer4/src -Llib/hmmer4/src/dp_reference -Llib/hmmer4/src/dp_vector -Llib/hmmer4/lib/easel  -o src/e2fhmmer_OA_utest src/e2msa.o src/e1_bg.o src/e1_model.o src/e1_rate.o src/e2_generic_decoding.o src/e2_generic_fwdback.o src/e2_generic_optacc.o src/e2hmmer_gmx.o src/e2_msa.o src/e2_pipeline.o src/e2_profile.o src/e2_profilesq.o src/e2_trace.o src/e2fhmmer_generic_decoding.o src/e2fhmmer_generic_fwdback.o src/e2fhmmer_generic_optacc.o src/branchinference.o src/evohmmer.o src/evomsa.o src/evopipeline.o src/fetchpfamdb.o src/homology.o src/logsum.o src/msatree.o src/mya.o src/muscle.o src/modelconfig.o src/orthologs.o src/plot.o src/ratematrix.o src/ratebuilder.o src/tracealign.o  src/Dickersonrates.o src/ssifile.o src/tkf_model.o   -lm -lhmmer -leasel  

 * ./e2fhmmer_OA_utest ../data/fn3.sto 
 */
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "msatree.h"


static ESL_OPTIONS options[] = {
  /* name           type       default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL,   NULL, NULL, NULL, "show brief help on version and usage",              0 },
  { "-v",         eslARG_NONE,   FALSE, NULL, NULL,   NULL, NULL, NULL, "be verbose",                                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msa>";
static char banner[] = "test driver for e2fhmmer_optacc.c";

int
main(int argc, char **argv)
{ 
  char           *msg = "OPTIMAL_ACCURACY unit test failed";
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char            errbuf[eslERRBUFSIZE];
  char           *msafile;
  ESLX_MSAFILE   *afp = NULL;
  ESL_MSA        *msa = NULL; 
  ESL_ALPHABET   *abc = NULL;
  ESL_TREE       *T = NULL;
  int             status = eslOK;
  int             hstatus = eslOK;
  int             verbose;

  if (esl_opt_ArgNumber(go) != 1)                  { puts("Incorrect number of command line arguments");        exit(1); }
  if ((msafile  = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <msafile> argument on command line");  exit(1); }

  /* Options */
  verbose = esl_opt_GetBoolean(go, "-v");

 /* Open the MSA file */
  status = eslx_msafile_Open(&abc, msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);

  /* read the MSA */
  hstatus = eslx_msafile_Read(afp, &msa);
  if (hstatus != eslOK) eslx_msafile_ReadFailure(afp, status);
  if (verbose) { if (eslx_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal(msg); }

  /* calculate the Tree using FastTree*/
  if (Tree_CalculateExtFromMSA(msa, &T, errbuf, verbose) != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  if (verbose) esl_tree_WriteNewick(stdout, T);

  /* root the Tree */
   if (Tree_InterLeafMaxDistRooted(T, NULL, errbuf, verbose) != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  if (verbose) esl_tree_WriteNewick(stdout, T);

  status = utest_e2fhmmer_OptimalAccuracy(msa, T, errbuf, verbose);
  if (status != eslOK)  { printf("%s\n", errbuf); esl_fatal(msg); }

  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);
  esl_msa_Destroy(msa);
  esl_tree_Destroy(T);
  eslx_msafile_Close(afp);
  return 0;
}
#endif /* E2OA_TESTDRIVE */


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
