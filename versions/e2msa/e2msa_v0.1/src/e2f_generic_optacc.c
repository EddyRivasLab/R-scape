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
#include "e2_gmx.h"
#include "e2f_generic_optacc.h"
#include "e2_profile.h"
#include "e2_profilesq.h"

/*****************************************************************
 * 1. Optimal alignment fill and traceback.
 *****************************************************************/
#define TSCDELTA(s)    ( (tsc[(s)]      == -eslINFINITY) ? FLT_MIN : 1.0)
#define XSCDELTA(s, x) ( (xsc[(s)][(x)] == -eslINFINITY) ? FLT_MIN : 1.0)

/* The TSCDELTA is used to make impossible paths impossible in the
 * optimal accuracy decoding algorithm; see Kall et al (2005). What we
 * want to do is multiply by a Kronecker delta that's 1 when the
 * transition probability is finite, and 0 when it's zero (when the
 * log prob is -eslINFINITY). But we can't do that easily, when we're
 * in log space, because 0 * -eslINFINITY = NaN. Instead, we use a
 * tiny number (FLT_MIN, ~1e-37).
 * 
 * A side concern is that we don't want to put a bunch of if-else
 * branches in the code; compilers should be able to generate more
 * efficient code from the TSCDELTA() construction.
 */


static int isgap(float *p, int K);

/* Function:  e2f_GOptimalAccuracy()
 * Synopsis:  Optimal accuracy decoding: fill. 
 * Incept:    ER, Thu Oct 10 09:40:51 EDT 2013
 *            adapted from p7_GOptimalAccuracy()
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm \citep{Kall05}.
 *            
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <L> using the query model <gm>.
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
e2f_GOptimalAccuracy(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, const E2_GMX *pp, E2_GMX *gx, float *ret_e)
{
  float const  *tsc  = (float const *)gm->tsc;
  float       **dp   = gx->dp;						
  float         sc;
  enum apair_e  apair;
  int           K = psql->abc->K;
  int           isgapl, isgapr;
  int           L;
  int           x;  
  int           xv; /* linear memory indices */
 
  /* OptAcc:
   *           states DD, EE needs to be evaluated last.
   *
   * Order: BB, SS, DS, SD, IB, IS, ID, BI, SI, II, DD, EE
   */
 
  if (gx->Lcol != gx->Lrow) { printf("sqs are not aligned\n"); return eslFAIL; }
  L = gx->Lcol;

  x = 0;
  BBMX(x) = 0.0;
  SSMX(x) = DSMX(x) = SDMX(x) = -eslINFINITY;
  IBMX(x) = ISMX(x) = IDMX(x) = -eslINFINITY;
  BIMX(x) = SIMX(x) = DIMX(x) = IIMX(x) = -eslINFINITY;
  DDMX(x) = -eslINFINITY;
  E2G_XMX(gx, x, e2G_EE) = -eslINFINITY;
  
#if 0
    printf("e2f_OA x %d  BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
	   x, 
	   BBMX(x), IBMX(x), SSMX(x), DSMX(x), ISMX(x), 
	   SDMX(x), DDMX(x), IDMX(x), 
	   BIMX(x), SIMX(x), DIMX(x), IIMX(x), E2G_XMX(gx, x, e2G_EE));
#endif

  /* Recursion. Done as a pull.
   */
  for (x = 1; x <= L; x++) {

    isgapl = isgap(psql->prof[x], K);
    isgapr = isgap(psqr->prof[x], K);
    if      (!isgapl && !isgapr) apair = PAIR;
    else if ( isgapl && !isgapr) apair = RRES;
    else if (!isgapl &&  isgapr) apair = LRES;
    else                         apair = NONE;
    
    xv = x - 1;
    
    /* BB state  0 transitions */
    BBMX(x) = -eslINFINITY;
    /* SS state 12 transitions */
    SSMX(x) = -eslINFINITY;
    if (apair == PAIR) { 
      sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_SS) * (BBMX(xv) + pp->dp[x][e2G_SS]), 
			   TSCDELTA(e2P_IB_SS) * (IBMX(xv) + pp->dp[x][e2G_SS])),
		   ESL_MAX(TSCDELTA(e2P_SS_SS) * (SSMX(xv) + pp->dp[x][e2G_SS]),
			   TSCDELTA(e2P_DS_SS) * (DSMX(xv) + pp->dp[x][e2G_SS])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_IS_SS) * (ISMX(xv) + pp->dp[x][e2G_SS])),
		   ESL_MAX(TSCDELTA(e2P_SD_SS) * (SDMX(xv) + pp->dp[x][e2G_SS]), 
			   TSCDELTA(e2P_DD_SS) * (DDMX(xv) + pp->dp[x][e2G_SS])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_ID_SS) * (IDMX(xv) + pp->dp[x][e2G_SS])),
		   ESL_MAX(TSCDELTA(e2P_BI_SS) * (BIMX(xv) + pp->dp[x][e2G_SS]), 
			   TSCDELTA(e2P_SI_SS) * (SIMX(xv) + pp->dp[x][e2G_SS])));
      sc = ESL_MAX(sc, 
		   ESL_MAX(TSCDELTA(e2P_DI_SS) * (DIMX(xv) + pp->dp[x][e2G_SS]), 
			   TSCDELTA(e2P_II_SS) * (IIMX(xv) + pp->dp[x][e2G_SS])));
      SSMX(x) = sc;	  
    }
    /* DS state 12 transitions */
    DSMX(x) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_DS) * (BBMX(xv) + pp->dp[x][e2G_DS]), 
			   TSCDELTA(e2P_IB_DS) * (IBMX(xv) + pp->dp[x][e2G_DS])),
		   ESL_MAX(TSCDELTA(e2P_SS_DS) * (SSMX(xv) + pp->dp[x][e2G_DS]),
			   TSCDELTA(e2P_DS_DS) * (DSMX(xv) + pp->dp[x][e2G_DS])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_IS_DS) * (ISMX(xv) + pp->dp[x][e2G_DS])),
		   ESL_MAX(TSCDELTA(e2P_SD_DS) * (SDMX(xv) + pp->dp[x][e2G_DS]), 
			   TSCDELTA(e2P_DD_DS) * (DDMX(xv) + pp->dp[x][e2G_DS])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_ID_DS) * (IDMX(xv) + pp->dp[x][e2G_DS])),
		   ESL_MAX(TSCDELTA(e2P_BI_DS) * (BIMX(xv) + pp->dp[x][e2G_DS]), 
			   TSCDELTA(e2P_SI_DS) * (SIMX(xv) + pp->dp[x][e2G_DS])));
      sc = ESL_MAX(sc, 
		   ESL_MAX(TSCDELTA(e2P_DI_DS) * (DIMX(xv) + pp->dp[x][e2G_DS]), 
			   TSCDELTA(e2P_II_DS) * (IIMX(xv) + pp->dp[x][e2G_DS])));
      DSMX(x) = sc;
    }	  
    /* SD state 12 transitions */
   SDMX(x) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_BB_SD) * (BBMX(xv) + pp->dp[x][e2G_SD]), 
			   TSCDELTA(e2P_IB_SD) * (IBMX(xv) + pp->dp[x][e2G_SD])),
		   ESL_MAX(TSCDELTA(e2P_SS_SD) * (SSMX(xv) + pp->dp[x][e2G_SD]),
			   TSCDELTA(e2P_DS_SD) * (DSMX(xv) + pp->dp[x][e2G_SD])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_IS_SD) * (ISMX(xv) + pp->dp[x][e2G_SD])),
		   ESL_MAX(TSCDELTA(e2P_SD_SD) * (SDMX(xv) + pp->dp[x][e2G_SD]), 
			   TSCDELTA(e2P_DD_SD) * (DDMX(xv) + pp->dp[x][e2G_SD])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_ID_SD) * (IDMX(xv) + pp->dp[x][e2G_SD])),
		   ESL_MAX(TSCDELTA(e2P_BI_SD) * (BIMX(xv) + pp->dp[x][e2G_SD]), 
			   TSCDELTA(e2P_SI_SD) * (SIMX(xv) + pp->dp[x][e2G_SD])));
      sc = ESL_MAX(sc, 
		   ESL_MAX(TSCDELTA(e2P_DI_SD) * (DIMX(xv) + pp->dp[x][e2G_SD]), 
			   TSCDELTA(e2P_II_SD) * (IIMX(xv) + pp->dp[x][e2G_SD])));
      SDMX(x) = sc;	  
    }
    
    /* IB state 2 transitions */
    IBMX(x) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(TSCDELTA(e2P_BB_IB) * (BBMX(xv) + pp->dp[x][e2G_IB]), 
		   TSCDELTA(e2P_IB_IB) * (IBMX(xv) + pp->dp[x][e2G_IB]));
      IBMX(x) = sc;
    }
    /* IS state 3 transitions  */
    ISMX(x) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(         TSCDELTA(e2P_SS_IS) * (SSMX(xv) + pp->dp[x][e2G_IS]),
		    ESL_MAX(TSCDELTA(e2P_DS_IS) * (DSMX(xv) + pp->dp[x][e2G_IS]),
			    TSCDELTA(e2P_IS_IS) * (ISMX(xv) + pp->dp[x][e2G_IS])));
      ISMX(x) = sc;
    }	  
    /* ID state 3 transitions */
    BBMX(x) = -eslINFINITY;
    if (apair == LRES) {
      sc = ESL_MAX(        TSCDELTA(e2P_SD_ID) * (SDMX(xv) + pp->dp[x][e2G_ID]),
		   ESL_MAX(TSCDELTA(e2P_DD_ID) * (DDMX(xv) + pp->dp[x][e2G_ID]),
			   TSCDELTA(e2P_ID_ID) * (IDMX(xv) + pp->dp[x][e2G_ID])));
      IDMX(x) = sc;	  
    }

    /* BI state 2 transitions */
   BIMX(x) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(TSCDELTA(e2P_BB_BI) * (BBMX(xv) + pp->dp[x][e2G_BI]), 
		   TSCDELTA(e2P_BI_BI) * (BIMX(xv) + pp->dp[x][e2G_BI]));
      BIMX(x) = sc;
    }
    /* SI state 3 transitions */
    SIMX(x) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(        TSCDELTA(e2P_SS_SI) * (SSMX(xv) + pp->dp[x][e2G_SI]),
		   ESL_MAX(TSCDELTA(e2P_SD_SI) * (SDMX(xv) + pp->dp[x][e2G_SI]),
			   TSCDELTA(e2P_SI_SI) * (SIMX(xv) + pp->dp[x][e2G_SI])));
      SIMX(x) = sc;
    }
    /* DI state 3 transitions */
   DIMX(x) = -eslINFINITY;
   if (apair == RRES) {
     sc = ESL_MAX(        TSCDELTA(e2P_DS_DI) * (DSMX(xv) + pp->dp[x][e2G_DI]),
		  ESL_MAX(TSCDELTA(e2P_DD_DI) * (DDMX(xv) + pp->dp[x][e2G_DI]),
			  TSCDELTA(e2P_DI_DI) * (DIMX(xv) + pp->dp[x][e2G_DI])));
     DIMX(x) = sc;
   }	  
   /* II state 4 transitions */
   IIMX(x) = -eslINFINITY;
    if (apair == RRES) {
      sc = ESL_MAX(ESL_MAX(TSCDELTA(e2P_IB_II) * (IBMX(xv) + pp->dp[x][e2G_II]),
			   TSCDELTA(e2P_IS_II) * (ISMX(xv) + pp->dp[x][e2G_II])),
		   ESL_MAX(TSCDELTA(e2P_ID_II) * (IDMX(xv) + pp->dp[x][e2G_II]),
			   TSCDELTA(e2P_II_II) * (IIMX(xv) + pp->dp[x][e2G_II])));
      IIMX(x) = sc;
    }

    /* DD state 10 transitions */
    DDMX(x) = -eslINFINITY;
    if (apair == NONE) {
      sc = ESL_MAX(        TSCDELTA(e2P_IB_DD) * (IBMX(xv) + pp->dp[x][e2G_DD]),
		   ESL_MAX(TSCDELTA(e2P_SS_DD) * (SSMX(xv) + pp->dp[x][e2G_DD]),
			   TSCDELTA(e2P_DS_DD) * (DSMX(xv) + pp->dp[x][e2G_DD])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_IS_DD) * (ISMX(xv) + pp->dp[x][e2G_DD])),
		   ESL_MAX(TSCDELTA(e2P_SD_DD) * (SDMX(xv) + pp->dp[x][e2G_DD]), 
			   TSCDELTA(e2P_DD_DD) * (DDMX(xv) + pp->dp[x][e2G_DD])));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_ID_DD) * (IDMX(xv) + pp->dp[x][e2G_DD])),
		   ESL_MAX(TSCDELTA(e2P_BI_DD) * (BIMX(xv) + pp->dp[x][e2G_DD]), 
			   TSCDELTA(e2P_SI_DD) * (SIMX(xv) + pp->dp[x][e2G_DD])));
      sc = ESL_MAX(sc, 
		   ESL_MAX(TSCDELTA(e2P_DI_DD) * (DIMX(xv) + pp->dp[x][e2G_DD]), 
			   TSCDELTA(e2P_II_DD) * (IIMX(xv) + pp->dp[x][e2G_DD])));
      DDMX(x) = sc; 
    }

    /* EE state 12 transitions */
    sc = -eslINFINITY;
    if (x == L) {
      sc = ESL_MAX(        TSCDELTA(e2P_IB_EE) * IBMX(x),
		   ESL_MAX(TSCDELTA(e2P_SS_EE) * SSMX(x),
			   TSCDELTA(e2P_DS_EE) * DSMX(x)));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_IS_EE) * ISMX(x)),
		   ESL_MAX(TSCDELTA(e2P_SD_EE) * SDMX(x), 
			   TSCDELTA(e2P_DD_EE) * DDMX(x)));
      sc = ESL_MAX(ESL_MAX(sc, 
			   TSCDELTA(e2P_ID_EE) * IDMX(x)),
		   ESL_MAX(TSCDELTA(e2P_BI_EE) * BIMX(x), 
			   TSCDELTA(e2P_SI_EE) * SIMX(x)));
      sc = ESL_MAX(sc, 
		   ESL_MAX(TSCDELTA(e2P_DI_EE) * DIMX(x), 
			   TSCDELTA(e2P_II_EE) * IIMX(x)));
    }
    E2G_XMX(gx, x, e2G_EE) = sc;

#if 0
    if (x>=0) printf("e2f_OA x %d  apair %d BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		     x, apair,
		     BBMX(x), IBMX(x), SSMX(x), DSMX(x), ISMX(x), 
		     SDMX(x), DDMX(x), IDMX(x), 
		     BIMX(x), SIMX(x), DIMX(x), IIMX(x), E2G_XMX(gx, x, e2G_EE));
#endif
    
  }
  
  if (ret_e != NULL) *ret_e = E2G_XMX(gx, L, e2G_EE);
  
  return eslOK;
}
/*---------------------- end, oa fill ---------------------------*/


static inline void  ancestral_res_fromSS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa);
static inline void  ancestral_res_fromSD(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa);
static inline void  ancestral_res_fromDS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa);
static inline void  ancestral_res_fromDD(const float *frq, const E2_GMX *pp,                           int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa);
static inline void  ancestral_res_isgap(int k, PSQ *sqa);

static inline float get_postprob(const E2_GMX *pp, int scur, int sprv, int x);
static inline int   find_argmax(ESL_RANDOMNESS *r, const float *v, int n);

static inline int   select_ib(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_ss(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_ds(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_is(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_sd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_dd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_id(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_bi(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_si(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_di(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_ii(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);
static inline int   select_ee(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L);

/*****************************************************************
 * 2. Optimal alignment accuracy, traceback
 *****************************************************************/

/* Function:  e2f_GOATrace()
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
e2f_GOATrace(ESL_RANDOMNESS *r, const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, const E2_PROFILE *gm, const E2_GMX *pp, const E2_GMX *gx, E2_TRACE *tr, const PSQ *psql, const PSQ *psqr, PSQ **ret_psqA)
{
  PSQ       *psqA   = NULL;
  int        L;
  int        x;                 /* position in alingment   */
  int        k      = 1;        /* position in the ancestral sequence    */
  float      postprob;
  int        sprv, scur;
  int        status;

#ifdef e2_DEBUGGING
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  psqA = psq_Create(gm->abc);

  if (psql->n != psqr->n) { printf("sqs are not aligned\n"); return eslFAIL; }
  L = psql->n;

  /* start */
  x = L;
  if ((status = e2_trace_AppendWithPP(tr, e2T_EE, 0, k, x, x, 0.0)) != eslOK) return status;
 
  sprv = e2T_EE;
  while (sprv != e2T_BB) 
    {
  
      switch (sprv) {
      case e2T_SS: scur = select_ss(r, gm, gx, x, L); ancestral_res_fromSS(evol, evor, frq,  x, psql, psqr,  k, psqA); x--; k++; psqA->n++; break;
      case e2T_DS: scur = select_ds(r, gm, gx, x, L); ancestral_res_fromDS(evol, evor, frq,  x, psql, psqr,  k, psqA); x--; k++; psqA->n++; break;
      case e2T_SD: scur = select_sd(r, gm, gx, x, L); ancestral_res_fromSD(evol, evor, frq,  x, psql, psqr,  k, psqA); x--; k++; psqA->n++; break;
      case e2T_DD: scur = select_dd(r, gm, gx, x, L); ancestral_res_fromDD(frq, pp,          x, psql, psqr,  k, psqA); x--; k++; psqA->n++; break;

      case e2T_IB: scur = select_ib(r, gm, gx, x, L); ancestral_res_isgap (                                  k, psqA); x--; k++; psqA->n++; break;
      case e2T_IS: scur = select_is(r, gm, gx, x, L); ancestral_res_isgap (                                  k, psqA); x--; k++; psqA->n++; break;
      case e2T_ID: scur = select_id(r, gm, gx, x, L); ancestral_res_isgap (                                  k, psqA); x--; k++; psqA->n++; break;
 
      case e2T_BI: scur = select_bi(r, gm, gx, x, L); ancestral_res_isgap (                                  k, psqA); x--; k++; psqA->n++; break;
      case e2T_SI: scur = select_si(r, gm, gx, x, L); ancestral_res_isgap (                                  k, psqA); x--; k++; psqA->n++; break;
      case e2T_DI: scur = select_di(r, gm, gx, x, L); ancestral_res_isgap (                                  k, psqA); x--; k++; psqA->n++; break;
      case e2T_II: scur = select_ii(r, gm, gx, x, L); ancestral_res_isgap (                                  k, psqA); x--; k++; psqA->n++; break;
      case e2T_EE: scur = select_ee(r, gm, gx, x, L);                                                                                       break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      postprob = get_postprob(pp, scur, sprv, x);
      if ((status = e2_trace_AppendWithPP(tr, scur, 0, k, x, x, postprob)) != eslOK) return status;
      psq_Grow(psqA, NULL);
      sprv = scur;
    }
 
  tr->M     = k; /* length of ancestral sequence+1  */
  tr->Lrow  = L;
  tr->Lcol  = L;
  tr->rowsq = gx->rowsq;
  if ((status = e2_trace_Reverse(tr)) != eslOK) goto ERROR;
 
  if ((status = psq_Reverse(psqA))    != eslOK) goto ERROR;

  *ret_psqA = psqA;
  return eslOK;

 ERROR:  
  if (psqA != NULL) psq_Destroy(psqA);
  return status;
}

static inline void
ancestral_res_fromSS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa)
{
  float *pl = sql->prof[x];
  float *pr = sqr->prof[x];
  int    K = sql->abc->K;
  int    a, b, c;
  
  if (isgap(pl, K) || isgap(pr, K)) { printf("ancestral_res_fromSS() you should not be here. x=%d\n", x); exit(1); }

  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  for (a = 0;  a < K; a++) 
    for (b = 0; b < K; b++) 
      for (c = 0; c < K; c++) 
	{	  
	  sqa->prof[k][a] = 
	    p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evol->sub->mx[a][b]) + log(evor->sub->mx[a][c]) + pl[b] + pr[c]);
	}
  esl_vec_FLogNorm(sqa->prof[k], K+1); 
  esl_vec_FLog    (sqa->prof[k], K+1); 
}

static inline void
ancestral_res_fromSD(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa)
{
  float *pl = sql->prof[x];
  float *pr = sqr->prof[x];
  int    K = sql->abc->K;
  int    a, b;
 
  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  if (!isgap(pr, K)) { printf("ancestral_res_fromSD() you should not be here. x=%d\n", x); exit(1); }

  for (a = 0;  a < K; a++) 
    for (b = 0; b < K; b++) 
      {	  
	sqa->prof[k][a] =  
	  p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evol->sub->mx[a][b]) + pl[b]);
      }
  
  esl_vec_FLogNorm(sqa->prof[k], K+1); 
  esl_vec_FLog    (sqa->prof[k], K+1); 
}

static inline void
ancestral_res_fromDS(const E1_MODEL *evol, const E1_MODEL *evor, const float *frq, int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa)
{
  float *pl = sql->prof[x];
  float *pr = sqr->prof[x];
  int    K  = sql->abc->K;
  int    a, b;
  
  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  if (!isgap(pl, K)) { printf("ancestral_res_fromDS() you should not be here. x=%d\n", x); exit(1); }

  for (a = 0;  a < K; a++) 
    for (b = 0; b < K; b++) 
      {	  
	sqa->prof[k][a] =
	  p7_FLogsum(sqa->prof[k][a], log(frq[a]) + log(evor->sub->mx[a][b]) + pr[b]);
      }
  
  esl_vec_FLogNorm(sqa->prof[k], K+1); 
  esl_vec_FLog    (sqa->prof[k], K+1); 
}

/* the number of ancestral residues involved
 * in a double deletion is governed by a geometric
 * distribution */
static inline void
ancestral_res_fromDD(const float *frq, const E2_GMX *pp, int x, const PSQ *sql, const PSQ *sqr, int k, PSQ *sqa)
{
  float *pl = sql->prof[x];
  float *pr = sqr->prof[x];
  int    K = sqa->abc->K;
  int    a;
  
  if (!isgap(pl, K) || !isgap(pr, K)) { printf("ancestral_res_fromDD() you should not be here. x=%d\n", x); exit(1); }

  /* add residues from the background distribution for the ancestral sequence */
  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  
  for (a = 0; a < K; a++) 
    sqa->prof[k][a] = log(frq[a]);
  esl_vec_FLogNorm(sqa->prof[k], K+1); 
  esl_vec_FLog    (sqa->prof[k], K+1); 
}

static inline float
get_postprob(const E2_GMX *pp, int scur, int sprv, int x)
{
  float **dp  = pp->dp;

  switch (scur) {
  case e2T_BB: return BBMX(x);
  case e2T_IB: return IBMX(x);
  case e2T_SS: return SSMX(x);
  case e2T_DS: return DSMX(x);
  case e2T_IS: return ISMX(x);
  case e2T_SD: return SDMX(x);
  case e2T_DD: return DDMX(x);
  case e2T_ID: return IDMX(x);
  case e2T_BI: return BIMX(x);
  case e2T_SI: return SIMX(x);
  case e2T_DI: return DIMX(x);
  case e2T_II: return IIMX(x);
  case e2T_EE: return E2G_XMX(pp, x, e2G_EE);
  default:    return 0.0;
  }
}

static inline void
ancestral_res_isgap(int k, PSQ *sqa)
{
  int K = sqa->abc->K;

  esl_vec_FSet(sqa->prof[k], K+1, -eslINFINITY); 
  sqa->prof[k][K] = 0.;
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
select_ib(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work  */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          xv;
  int          state[2] = { e2T_BB, e2T_IB };
  
  xv = x - 1;

  path[0]  = TSCDELTA(e2P_BB_IB) * BBMX(xv);
  path[1]  = TSCDELTA(e2P_IB_IB) * IBMX(xv);

  return state[find_argmax(r, path, 2)];
}

static inline int
select_ss(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[12];
  int          xv;
  int          state[12] = { e2T_BB, e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };
  
  xv = x - 1;

  path[0]  = TSCDELTA(e2P_BB_SS) * BBMX(xv);
  path[1]  = TSCDELTA(e2P_IB_SS) * IBMX(xv);
  path[2]  = TSCDELTA(e2P_SS_SS) * SSMX(xv);
  path[3]  = TSCDELTA(e2P_DS_SS) * DSMX(xv);
  path[4]  = TSCDELTA(e2P_IS_SS) * ISMX(xv);
  path[5]  = TSCDELTA(e2P_SD_SS) * SDMX(xv);
  path[6]  = TSCDELTA(e2P_DD_SS) * DDMX(xv);
  path[7]  = TSCDELTA(e2P_ID_SS) * IDMX(xv);
  path[8]  = TSCDELTA(e2P_BI_SS) * BIMX(xv);
  path[9]  = TSCDELTA(e2P_SI_SS) * SIMX(xv);
  path[10] = TSCDELTA(e2P_DI_SS) * DIMX(xv);
  path[11] = TSCDELTA(e2P_II_SS) * IIMX(xv);
  return state[find_argmax(r, path, 12)];
}

static inline int
select_ds(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[12];
  int          xv;
  int          state[12] = { e2T_BB, e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };
  
  xv = x - 1;

  path[0]  = TSCDELTA(e2P_BB_DS) * BBMX(xv);
  path[1]  = TSCDELTA(e2P_IB_DS) * IBMX(xv);
  path[2]  = TSCDELTA(e2P_SS_DS) * SSMX(xv);
  path[3]  = TSCDELTA(e2P_DS_DS) * DSMX(xv);
  path[4]  = TSCDELTA(e2P_IS_DS) * ISMX(xv);
  path[5]  = TSCDELTA(e2P_SD_DS) * SDMX(xv);
  path[6]  = TSCDELTA(e2P_DD_DS) * DDMX(xv);
  path[7]  = TSCDELTA(e2P_ID_DS) * IDMX(xv);
  path[8]  = TSCDELTA(e2P_BI_DS) * BIMX(xv);
  path[9]  = TSCDELTA(e2P_SI_DS) * SIMX(xv);
  path[10] = TSCDELTA(e2P_DI_DS) * DIMX(xv);
  path[11] = TSCDELTA(e2P_II_DS) * IIMX(xv);

  return state[find_argmax(r, path, 12)];
}

static inline int
select_is(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          xv;
  int          state[3] = { e2T_SS, e2T_DS, e2T_IS };
  
  xv = x - 1;

  path[0]  = TSCDELTA(e2P_SS_IS) * SSMX(xv);
  path[1]  = TSCDELTA(e2P_DS_IS) * DSMX(xv);
  path[2]  = TSCDELTA(e2P_IS_IS) * ISMX(xv);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_sd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[12];
  int          xv;
  int          state[12] = { e2T_BB, e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };
  
  xv = x - 1;
 
  path[0]  = TSCDELTA(e2P_BB_SD) * BBMX(xv);
  path[1]  = TSCDELTA(e2P_IB_SD) * IBMX(xv);
  path[2]  = TSCDELTA(e2P_SS_SD) * SSMX(xv);
  path[3]  = TSCDELTA(e2P_DS_SD) * DSMX(xv);
  path[4]  = TSCDELTA(e2P_IS_SD) * ISMX(xv);
  path[5]  = TSCDELTA(e2P_SD_SD) * SDMX(xv);
  path[6]  = TSCDELTA(e2P_DD_SD) * DDMX(xv);
  path[7]  = TSCDELTA(e2P_ID_SD) * IDMX(xv);
  path[8]  = TSCDELTA(e2P_BI_SD) * BIMX(xv);
  path[9]  = TSCDELTA(e2P_SI_SD) * SIMX(xv);
  path[10] = TSCDELTA(e2P_DI_SD) * DIMX(xv);
  path[11] = TSCDELTA(e2P_II_SD) * IIMX(xv);

  return state[find_argmax(r, path, 12)];
}

static inline int
select_dd(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[11];
  int          xv;
  int          state[11] = { e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, 
			     e2T_DD, e2T_ID, e2T_BI, e2T_SI, e2T_DI, e2T_II  };

  xv = x - 1;
 
  path[0]  = TSCDELTA(e2P_IB_DD) * IBMX(xv);
  path[1]  = TSCDELTA(e2P_SS_DD) * SSMX(xv);
  path[2]  = TSCDELTA(e2P_DS_DD) * DSMX(xv);
  path[3]  = TSCDELTA(e2P_IS_DD) * ISMX(xv);
  path[4]  = TSCDELTA(e2P_SD_DD) * SDMX(xv);
  path[5]  = TSCDELTA(e2P_DD_DD) * DDMX(xv);
  path[6]  = TSCDELTA(e2P_ID_DD) * IDMX(xv);
  path[7]  = TSCDELTA(e2P_BI_DD) * BIMX(xv);
  path[8]  = TSCDELTA(e2P_SI_DD) * SIMX(xv);
  path[9]  = TSCDELTA(e2P_DI_DD) * DIMX(xv);
  path[10] = TSCDELTA(e2P_II_DD) * IIMX(xv);

  return state[find_argmax(r, path, 11)];
}

static inline int
select_id(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          xv;
  int          state[3] = { e2T_SD, e2T_DD, e2T_ID  };
  
  xv = x - 1;

  path[0]  = TSCDELTA(e2P_SD_ID) * SDMX(xv);
  path[1]  = TSCDELTA(e2P_DD_ID) * DDMX(xv);
  path[2]  = TSCDELTA(e2P_ID_ID) * IDMX(xv);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_bi(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[2];
  int          xv;
  int          state[2] = { e2T_BB, e2T_BI };
  
  xv = x - 1;

  path[0]  = TSCDELTA(e2P_BB_BI) * BBMX(xv);
  path[1]  = TSCDELTA(e2P_BI_BI) * BIMX(xv);

  return state[find_argmax(r, path, 2)];
}

static inline int
select_si(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          xv;
  int          state[3] = { e2T_SS, e2T_SD, e2T_SI };
  
  xv = x - 1;

  path[0]  = TSCDELTA(e2P_SS_SI) * SSMX(xv);
  path[1]  = TSCDELTA(e2P_SD_SI) * SDMX(xv);
  path[2]  = TSCDELTA(e2P_SI_SI) * SIMX(xv);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_di(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[3];
  int          xv;
  int          state[3] = { e2T_DS, e2T_DD, e2T_DI };
  
  xv = x - 1;
 
  path[0]  = TSCDELTA(e2P_DS_DI) * DSMX(xv);
  path[1]  = TSCDELTA(e2P_DD_DI) * DDMX(xv);
  path[2]  = TSCDELTA(e2P_DI_DI) * DIMX(xv);

  return state[find_argmax(r, path, 3)];
}

static inline int
select_ii(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work  */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[4];
  int          xv;
  int          state[4] = { e2T_IB, e2T_IS, e2T_ID , e2T_II };
  
  xv = x -1;

  path[0]  = TSCDELTA(e2P_IB_II) * IBMX(xv);
  path[1]  = TSCDELTA(e2P_IS_II) * ISMX(xv);
  path[2]  = TSCDELTA(e2P_ID_II) * IDMX(xv);
  path[3]  = TSCDELTA(e2P_II_II) * IIMX(xv);
  return state[find_argmax(r, path, 4)];
}

static inline int
select_ee(ESL_RANDOMNESS *r, const E2_PROFILE *gm, const E2_GMX *gx, int x, int L)
{
  float      **dp   = gx->dp;	/* so {MDI}MX() macros work       */
  float const *tsc  = gm->tsc;	/* so TSCDELTA() macro works */
  float        path[11];
  int          state[121] = { e2T_IB, e2T_SS, e2T_DS, e2T_IS, e2T_SD, e2T_DD, e2T_ID, e2T_BI,  e2T_SI, e2T_DI, e2T_II};

  path[0]  = TSCDELTA(e2P_IB_EE) * IBMX(x);
  path[1]  = TSCDELTA(e2P_SS_EE) * SSMX(x);
  path[2]  = TSCDELTA(e2P_DS_EE) * DSMX(x);
  path[3]  = TSCDELTA(e2P_IS_EE) * ISMX(x);
  path[4]  = TSCDELTA(e2P_SD_EE) * SDMX(x);
  path[5]  = TSCDELTA(e2P_DD_EE) * DDMX(x);
  path[6]  = TSCDELTA(e2P_ID_EE) * IDMX(x);
  path[7]  = TSCDELTA(e2P_BI_EE) * BIMX(x);
  path[8]  = TSCDELTA(e2P_SI_EE) * SIMX(x);
  path[9]  = TSCDELTA(e2P_DI_EE) * DIMX(x);
  path[10] = TSCDELTA(e2P_II_EE) * IIMX(x);

  return state[find_argmax(r, path, 11)];
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

static int run_e2f(const PSQ *sql, const PSQ *sqr, PSQ **ret_sqa, E1_RATE *R, E1_BG *bg, float time1, float time2, double tol, int verbose);

static int
utest_e2f_OptimalAccuracy(ESL_MSA *msa, ESL_TREE *T, char *errbuf, int verbose)
{
  E1_RATE     *R = NULL;                                    /* assumes a unique set of rates for the whole process */
  E1_BG       *bg = NULL;
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
run_e2f(const PSQ *sql, const PSQ *sqr, PSQ **ret_sqa, E1_RATE *R, E1_BG *bg, float time1, float time2, double tol, int verbose)
{
  E1_MODEL       *evol = NULL;
  E1_MODEL       *evor = NULL;
  E2_TRACE       *tr   = NULL;
  E2_GMX         *gx1  = NULL;
  E2_GMX         *gx2  = NULL;
  E2_PROFILE     *gm   = NULL;
  ESL_ALPHABET   *abc = (ESL_ALPHABET *)sql->abc;
  PSQ            *sqa = psq_Create(abc);
  float           fsc, bsc;
  float           accscore;
  float           L;         /* average length of the two sequences */
  int             status;

  if (sql->abc->type != sqr->abc->type) {status = eslFAIL; goto ERROR; }

  /* Allocations */
  gx1 = e2_gmx_Create(sql->n, sqr->n);
  gx2 = e2_gmx_Create(sql->n, sqr->n);
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
  e2f_GForward (sql, sqr, gm, gx1, &fsc);
  e2f_GBackward(sql, sqr, gm, gx2, &bsc);
  if (verbose) {
    printf("e2f_fwd = %.4f nats\n", fsc);
    printf("e2f_bck = %.4f nats\n", bsc);
  }
  if (fabs(bsc-fsc) > 0.01*L) { status = eslFAIL; goto ERROR; }

  e2f_GDecoding(gm, gx1, gx2, gx2);                   /* <gx2> is now the posterior decoding matrix */
  e2f_GOptimalAccuracy(gm, gx2, gx1, &accscore);      /* <gx1> is now the OA matrix */

  e2f_GOATrace(evol, evor, bg, gm, gx2, gx1, tr, sql, sqr, &sqa);  
  if (verbose) {
    printf("ancestral sequence length %d\n", (int)sqa->n);
    printf("acc = %.4f (%.2f%%)\n", accscore, (sqa->n > 0)? accscore * 100. / (float) sqa->n : accscore * 100);
    e2_trace_Dump(stdout, tr, gm, sql, sqr, sqa);
  }
 
  e1_model_Destroy(evol);
  e1_model_Destroy(evor);
  e2_gmx_Destroy(gx1);
  e2_gmx_Destroy(gx2);
  e2_profile_Destroy(gm);
  e2_trace_Destroy(tr);
 
  *ret_sqa = sqa;
  return eslOK;

 ERROR:
  if (evol) e1_model_Destroy(evol);
  if (evor) e1_model_Destroy(evor);
  if (gx1)  e2_gmx_Destroy(gx1);
  if (gx2)  e2_gmx_Destroy(gx2);
  if (gm)   e2_profile_Destroy(gm);
  if (tr)   e2_trace_Destroy(tr);
  return status;
}
#endif /* E2OA_TESTDRIVE */

/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef E2FOA_TESTDRIVE
/* gcc -std=gnu99 -Isrc -Ilib/hmmer4/src -Ilib/hmmer4/src/dp_reference -Ilib/hmmer4/src/dp_vector -Ilib/hmmer4/lib/easel -g -Wall -DE2FOA_TESTDRIVE -Isrc -Ilib/hmmer4/src -Ilib/hmmer4/src/dp_reference -Ilib/hmmer4/src/dp_vector -Ilib/hmmer4/lib/easel  -Lsrc -Llib/hmmer4/src -Llib/hmmer4/src/dp_reference -Llib/hmmer4/src/dp_vector -Llib/hmmer4/lib/easel  -o src/e2f_OA_utest src/e2msa.o src/e1_bg.o src/e1_model.o src/e1_rate.o src/e2_generic_decoding.o src/e2_generic_fwdback.o src/e2_generic_optacc.o src/e2_gmx.o src/e2_msa.o src/e2_pipeline.o src/e2_profile.o src/e2_profilesq.o src/e2_trace.o src/e2f_generic_decoding.o src/e2f_generic_fwdback.o src/e2f_generic_optacc.o src/branchinference.o src/evohmmer.o src/evomsa.o src/evopipeline.o src/fetchpfamdb.o src/homology.o src/logsum.o src/msatree.o src/mya.o src/muscle.o src/modelconfig.o src/orthologs.o src/plot.o src/ratematrix.o src/ratebuilder.o src/tracealign.o  src/Dickersonrates.o src/ssifile.o src/tkf_model.o   -lm -lhmmer -leasel  

 * ./e2f_OA_utest ../data/fn3.sto 
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
static char banner[] = "test driver for e2f_optacc.c";

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

  status = utest_e2f_OptimalAccuracy(msa, T, errbuf, verbose);
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
