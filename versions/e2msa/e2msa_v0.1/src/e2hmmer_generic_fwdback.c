/* Forward/Backward algorithms; generic (non-SIMD) versions.
 * 
 * Contents:
 *   1. Forward, Backward implementations.  
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 */

#include "p7_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "e2.h"
#include "e2_profilesq.h"
#include "e2fhmmer_generic_fwdback.h"
#include "e2_profile.h"
#include "e2hmmer_profile.h"


/*****************************************************************
 * 1. Forward, Backward implementations.
 *****************************************************************/

/* Function:  e2hmmer_GForward()
 * Synopsis:  The Forward algorithm.
 *
 * Purpose:   The Forward dynamic programming algorithm. 
 *
 *            Given two aligned profile sequences <psql> and <psqr>,  
 *            of length L, a 
 *            model profile <gm>, and DP matrix <gx> allocated
 *            in linear time;
 *            calculate the probability of the sequence
 *            given the model using the Forward algorithm; return the
 *            Forward matrix in <gx>, and the Forward score in <ret_sc>.
 *           
 *            The Forward score is in lod score form.
 *           
 *            21 states
 *            110 transtions total
 *
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      psql   - profile sequence, 1..L
 *            psqr   - profile sequence, 1..L
 *            gm     - e2 model profile. 
 *            gx     - DP matrix with linear memory allocation
 *            opt_sc - optRETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
e2hmmer_GForward(const PSQ *psql, const PSQ *psqr, const E2HMMER_PROFILE *gm, E2_GMX *gx, float *opt_sc)
{
  PSQ          *sqrow;                          /* larger  sequence moves by rows */
  PSQ          *sqcol;                          /* shorter sequence moves by cols */
  float const **tsc = (float const **)gm->tsc;
  float       **dp  = gx->dp;						
  float         sc;
  float         SSsc;
  float         Slsc, Srsc;
  float         Ilsc, Irsc;
  float        *pl = NULL;
  float        *pr = NULL;
  float         EE = -eslINFINITY;
  int           M = gm->M;
  int           K = psql->abc->K;
  int           isgapl, isgapr;
  enum apair_e  apair;
  int           b,c;                  /* alphabet index */
  int           rowsq;
  int           i, j;  
  int           ijx, ijv, ipj, ijp;   /* linear memor indices */
  int           status;
  
  p7_FLogsumInit();    
  
  sqcol = (psql->n <= psqr->n)? (PSQ *)psql : (PSQ *)psqr;
  sqrow = (psql->n <= psqr->n)? (PSQ *)psqr : (PSQ *)psql;
  rowsq = (psql->n <= psqr->n)? e2P_SR      : e2P_SL;
  gx->rowsq = rowsq;

  /* reallocation, if needed */
  if ( (status = e2_gmx_GrowTo(gx, M, sqrow->n, sqcol->n)) != eslOK) return status;
  gx->M    = M;
  gx->L    = L;
  gx->type = E2_FORWARD;


  /* Forward:
   *           states DD, EE need to be evaluated last.
   *
   * Order: BB, SS, DS, SD, IB, IS, ID, BI, SI, II, DD, EE
   */

  /* Initialization of k=0 x=0 */
  x = 0;
  k = 0;
  BBMXM(x,k) = 0.0;
  SSMXM(x,k) = DSMXM(x,k) = SDMXM(x,k) = -eslINFINITY;
  IBMXM(x,k) = ISMXM(x,k) = IDMXM(x,k) = -eslINFINITY;
  BIMXM(x,k) = SIMXM(x,k) = DIMXM(x,k) = IiMXM(x,k) = iIMXM(x,k) = -eslINFINITY;
  DDMXM(x,k) = -eslINFINITY; 
  EEMXM(x,k) = -eslINFINITY;

  /* Initialization of the x=0 k > 0 rows */
  for (k = 1; k <= M; k ++) {
    kv = k - 1;

    /* BB state 0 transitions  */
    BBMXM(x,k) = 0.0; 
    /* SS state 0 transitions */
    SSMXM(x,k) = -eslINFINITY; 
    /* DS state 0 transitions */
    DSMXM(x,k) = -eslINFINITY;
    /* SD state 0 transitions */
    SDMXM(x,k) = -eslINFINITY;
    /* DD state 1 transitions */
    DDMXM(x,k) = DDMXM(x,kv) + e2hmmerTSC(kv,e2HP_DD_DD);
    /* IB state 0 transitions */
    IBMXM(x,k) = -eslINFINITY;
    /* IS state 0 transitions  */
    ISMXM(x,k) = -eslINFINITY;
    /* ID state 0 transitions */
    IDMXM(x,k) = -eslINFINITY;
    /* BI state 0 transitions */
    BIMXM(x,k) = -eslINFINITY;
    /* SI state 0 transitions */
    SIMXM(x,k) = -eslINFINITY;
    /* DI state 0 transitions */
    DIMXM(x,k) = -eslINFINITY;
    /* Ii state 0 transitions */
    IiMXM(x,k) = -eslINFINITY;
    /* iI state 0 transitions */
    iIMXM(x,k) = -eslINFINITY;
    /* EE state 13 transitions */
    if (x == L) {
      sc = p7_FLogsum(p7_FLogsum(BBMXM(x,k) + e2hmmerTSC(k,e2HP_BB_EE), 
				 IBMXM(x,k) + e2hmmerTSC(k,e2HP_IB_EE)),
		      p7_FLogsum(SSMXM(x,k) + e2hmmerTSC(k,e2HP_SS_EE),
				 DSMXM(x,k) + e2hmmerTSC(k,e2HP_DS_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 ISMXM(x,k) + e2hmmerTSC(k,e2HP_IS_EE)),
		      p7_FLogsum(SDMXM(x,k) + e2hmmerTSC(k,e2HP_SD_EE), 
				 DDMXM(x,k) + e2hmmerTSC(k,e2HP_DD_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 IDMXM(x,k) + e2hmmerTSC(k,e2HP_ID_EE)),
		      p7_FLogsum(BIMXM(x,k) + e2hmmerTSC(k,e2HP_BI_EE), 
				 SIMXM(x,k) + e2hmmerTSC(k,e2HP_SI_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 DIMXM(x,k) + e2hmmerTSC(k,e2HP_DI_EE)),
		      p7_FLogsum(IiMXM(x,k) + e2hmmerTSC(k,e2HP_Ii_EE), 
				 iIMXM(x,k) + e2hmmerTSC(k,e2HP_iI_EE)));
      EEMXM(x,k) = sc;
      EE = p7_FLogsum(EE, EEMXM(x,k));
    }
    else EEMXM(x,k) = -eslINFINITY;
  }

  /* Initialization of the k=0 row. */
  k = 0;
  for (x = 1; x <= L; x++) {
    xv = x - 1;
    
    pl = psql->prof[x];
    pr = psqr->prof[x];
    
    Ilsc = -eslINFINITY; /* insertion score */
    Irsc = -eslINFINITY; /* insertion score */
    for (b = 0; b < K; b ++) {
      Ilsc = p7_FLogsum(Ilsc, gm->isc[k][b][e2HP_SL] + pl[b]);
      Irsc = p7_FLogsum(Irsc, gm->isc[k][b][e2HP_SR] + pr[b]);
    }
    isgapl = isgap(pl, K);
    isgapr = isgap(pr, K);
    if      (!isgapl && !isgapr) apair = PAIR;
    else if ( isgapl && !isgapr) apair = RRES;
    else if (!isgapl &&  isgapr) apair = LRES;
    else                         apair = NONE;
    
    /* BB state 0 transitions  */
    BBMXM(x,k) = -eslINFINITY;
    if (apair == NONE) BBMXM(x,k) = BBMXM(xv,k);
    /* SS state 0 transitions */
    SSMXM(x,k) = -eslINFINITY;    
    if (apair == NONE) SSMXM(x,k) = SSMXM(xv,k);
    /* DS state 0 transitions */
    DSMXM(x,k) = -eslINFINITY;
    if (apair == NONE) DSMXM(x,k) = DSMXM(xv,k);
    /* SD state 0 transitions */
    SDMXM(x,k) = -eslINFINITY;
    if (apair == NONE) SDMXM(x,k) = SDMXM(xv,k);
    /* DD state 0 transitions */
    DDMXM(x,k) = -eslINFINITY;
    if (apair == NONE) DDMXM(x,k) = DDMXM(xv,k);
    /* IB state 2 transitions */
    IBMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = p7_FLogsum(BBMXM(xv,k) + e2hmmerTSC(k, e2HP_BB_IB), 
		      IBMXM(xv,k) + e2hmmerTSC(k, e2HP_IB_IB));
      IBMXM(x,k) = sc + Ilsc;
    }
    if (apair == NONE) IBMXM(x,k) = IBMXM(xv,k);
    /* IS state 2 transitions  */
    ISMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = p7_FLogsum(SSMXM(xv,k) + e2hmmerTSC(k, e2HP_SS_IS),
		      ISMXM(xv,k) + e2hmmerTSC(k, e2HP_IS_IS));
      ISMXM(x,k) = sc + Ilsc;
    }
    if (apair == NONE) ISMXM(x,k) = ISMXM(xv,k);
    /* ID state 2 transitions */
    IDMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = p7_FLogsum(SDMXM(xv,k) + e2hmmerTSC(k, e2HP_SD_ID),
		      IDMXM(xv,k) + e2hmmerTSC(k, e2HP_ID_ID));
      IDMXM(x,k) = sc + Ilsc;
    }      
    if (apair == NONE) IDMXM(x,k) = IDMXM(xv,k);

    /* BI state 2 transitions */
    BIMXM(x,k) = -eslINFINITY;
    if (apair ==  RRES) {
      sc = p7_FLogsum(BBMXM(xv,k) + e2hmmerTSC(k, e2HP_BB_BI), 
		      BIMXM(xv,k) + e2hmmerTSC(k, e2HP_BI_BI));
      BIMXM(x,k) = sc + Irsc;
    }
    if (apair == NONE) BIMXM(x,k) = BIMXM(xv,k);

   /* SI state 2 transitions */
    SIMXM(x,k) = -eslINFINITY;
    if (apair == RRES) {
      sc = p7_FLogsum(SSMXM(xv,k) + e2hmmerTSC(k, e2HP_SS_SI),
		      SIMXM(xv,k) + e2hmmerTSC(k, e2HP_SI_SI));
      SIMXM(x,k) = sc + Irsc;
    }
    if (apair == NONE) SIMXM(x,k) = SIMXM(xv,k);

    /* DI state 2 transitions */
    DIMXM(x,k) = -eslINFINITY;
    if (apair == RRES) {
      sc = p7_FLogsum(DSMXM(xv,k) + e2hmmerTSC(k,e2HP_DS_DI),
		      DIMXM(xv,k) + e2hmmerTSC(k,e2HP_DI_DI));
      DIMXM(x,k) = sc + Irsc;
    }
   if (apair == NONE) DIMXM(x,k) = DIMXM(xv,k);

    /* Ii state 4 transitions */
    IiMXM(x,k) = -eslINFINITY;
    if (apair == LRES) {
      sc = p7_FLogsum(p7_FLogsum(BIMXM(xv,k) + e2hmmerTSC(k,e2HP_BI_Ii),
				 SIMXM(xv,k) + e2hmmerTSC(k,e2HP_SI_Ii)),
		      p7_FLogsum(IiMXM(xv,k) + e2hmmerTSC(k,e2HP_Ii_Ii),
		                 iIMXM(xv,k) + e2hmmerTSC(k,e2HP_iI_Ii)));
      IiMXM(x,k) = sc + Ilsc;
    }
    if (apair == NONE) IiMXM(x,k) = IiMXM(xv,k);

    /* iI state 4 transitions */
    iIMXM(x,k) = -eslINFINITY;
    if (apair == RRES) {
      sc = p7_FLogsum(p7_FLogsum(IBMXM(xv,k) + e2hmmerTSC(k,e2HP_IB_iI),
				 ISMXM(xv,k) + e2hmmerTSC(k,e2HP_IS_iI)),
		      p7_FLogsum(IiMXM(xv,k) + e2hmmerTSC(k,e2HP_Ii_iI),
		                 iIMXM(xv,k) + e2hmmerTSC(k,e2HP_iI_iI)));
      iIMXM(x,k) = sc + Irsc;
    }
    if (apair == NONE) iIMXM(x,k) = iIMXM(xv,k);
   
   /* EE state 13 transitions */
    EEMXM(x,k) = -eslINFINITY;
    if (x == L) {
      sc = p7_FLogsum(p7_FLogsum(BBMXM(x,k) + e2hmmerTSC(k,e2HP_BB_EE), 
				 IBMXM(x,k) + e2hmmerTSC(k,e2HP_IB_EE)),
		      p7_FLogsum(SSMXM(x,k) + e2hmmerTSC(k,e2HP_SS_EE),
				 DSMXM(x,k) + e2hmmerTSC(k,e2HP_DS_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 ISMXM(x,k) + e2hmmerTSC(k,e2HP_IS_EE)),
		      p7_FLogsum(SDMXM(x,k) + e2hmmerTSC(k,e2HP_SD_EE), 
				 DDMXM(x,k) + e2hmmerTSC(k,e2HP_DD_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 IDMXM(x,k) + e2hmmerTSC(k,e2HP_ID_EE)),
		      p7_FLogsum(BIMXM(x,k) + e2hmmerTSC(k,e2HP_BI_EE), 
				 SIMXM(x,k) + e2hmmerTSC(k,e2HP_SI_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc,
				 DIMXM(x,k) + e2hmmerTSC(k,e2HP_DI_EE)),
		      p7_FLogsum(IiMXM(x,k) + e2hmmerTSC(k,e2HP_Ii_EE), 
				 iIMXM(x,k) + e2hmmerTSC(k,e2HP_iI_EE)));
      EEMXM(x,k) = sc;
      EE = p7_FLogsum(EE, EEMXM(x,k));
    } 

#if 0
    if (x==L)printf("e2fFWD x %d k %d/%d apair %d BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f iI %f  EE %f\n", 
		    x, k, M, apair, 
		    BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), 
		    SDMXM(x,k), DDMXM(x,k), IDMXM(x,k), 
		    BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k), EEMXM(x,k));
#endif
   }
  
  /* main recursion */
  for (k = 1; k <= M; k ++) {
    
    kv = k - 1;

   for (x = 1; x <= L; x++) {
      xv = x - 1;
      
      pl = psql->prof[x];
      pr = psqr->prof[x];
      
      SSsc = -eslINFINITY;
      Slsc = -eslINFINITY; /* orphan substitution score */
      Srsc = -eslINFINITY; /* orphan substitution score */
      Ilsc = -eslINFINITY; /* insertion score */
      Irsc = -eslINFINITY; /* insertion score */
      for (b = 0; b < K; b ++) {
	Slsc = p7_FLogsum(Slsc, gm->ssc[k][b][e2HP_SL] + pl[b]);
	Ilsc = p7_FLogsum(Ilsc, gm->isc[k][b][e2HP_SL] + pl[b]);
	Srsc = p7_FLogsum(Srsc, gm->ssc[k][b][e2HP_SR] + pr[b]);
	Irsc = p7_FLogsum(Irsc, gm->isc[k][b][e2HP_SR] + pr[b]);
	for (c = 0; c < K; c ++) {
	  SSsc = p7_FLogsum(SSsc, gm->sssc[k][b][c] + pl[b] + pr[c]);
	}
      }
      isgapl = isgap(pl, K);
      isgapr = isgap(pr, K);
      if      (!isgapl && !isgapr) apair = PAIR;
      else if ( isgapl && !isgapr) apair = RRES;
      else if (!isgapl &&  isgapr) apair = LRES;
      else                         apair = NONE;
      
      /* BB state 0 transitions  */
      BBMXM(x,k) = -eslINFINITY; 
      if (apair == NONE) BBMXM(x,k) = BBMXM(xv,k);

      /* SS state 13 transitions */
      SSMXM(x,k) = -eslINFINITY;
      if (apair == PAIR) { 
	sc = p7_FLogsum(p7_FLogsum(BBMXM(xv,kv) + e2hmmerTSC(kv, e2HP_BB_SS), 
				   IBMXM(xv,kv) + e2hmmerTSC(kv, e2HP_IB_SS)),
			p7_FLogsum(SSMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SS_SS),
				   DSMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DS_SS)));
	sc = p7_FLogsum(p7_FLogsum(sc, 
				   ISMXM(xv,kv) + e2hmmerTSC(kv, e2HP_IS_SS)),
			p7_FLogsum(SDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SD_SS), 
				   DDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DD_SS)));
	sc = p7_FLogsum(p7_FLogsum(sc, 
				   IDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_ID_SS)),
			p7_FLogsum(BIMXM(xv,kv) + e2hmmerTSC(kv, e2HP_BI_SS), 
				   SIMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SI_SS)));
	sc = p7_FLogsum(p7_FLogsum(sc, 
				   DIMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DI_SS)),
			p7_FLogsum(IiMXM(xv,kv) + e2hmmerTSC(kv, e2HP_Ii_SS), 
				   iIMXM(xv,kv) + e2hmmerTSC(kv, e2HP_iI_SS)));
	SSMXM(x,k) = sc + SSsc;
      }
      if (apair == NONE) SSMXM(x,k) = SSMXM(xv,k);

      /* DS state 6 transitions */
      DSMXM(x,k) = -eslINFINITY;
      if (apair == RRES) {
	sc = p7_FLogsum(p7_FLogsum(SSMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SS_DS), 
				   DSMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DS_DS)),
			p7_FLogsum(SDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SD_DS),
				   DDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DD_DS)));
	sc = p7_FLogsum(sc, 
			p7_FLogsum(SIMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SI_DS), 
				   DIMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DI_DS)));
	DSMXM(x,k) = sc + Srsc;
      }
      if (apair == NONE) DSMXM(x,k) = DSMXM(xv,k);
     
      /* SD state 6 transitions */
      SDMXM(x,k) = -eslINFINITY;
      if (apair == LRES) {
	sc = p7_FLogsum(p7_FLogsum(SSMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SS_SD), 
				   DSMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DS_SD)),
			p7_FLogsum(ISMXM(xv,kv) + e2hmmerTSC(kv, e2HP_IS_SD),
				   SDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_SD_SD)));
	sc = p7_FLogsum(sc, 
			p7_FLogsum(DDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_DD_SD), 
				   IDMXM(xv,kv) + e2hmmerTSC(kv, e2HP_ID_SD)));
	SDMXM(x,k) = sc + Slsc;
      }
      if (apair == NONE) SDMXM(x,k) = SDMXM(xv,k);
      
      /* DD state 4 transitions */
      sc = p7_FLogsum(p7_FLogsum(SSMXM(x,kv) + e2hmmerTSC(kv,e2HP_SS_DD), 
			         SDMXM(x,kv) + e2hmmerTSC(kv,e2HP_SD_DD)),
		      p7_FLogsum(DSMXM(x,kv) + e2hmmerTSC(kv,e2HP_DS_DD),
		  	         DDMXM(x,kv) + e2hmmerTSC(kv,e2HP_DD_DD)));
      DDMXM(x,k) = sc;
      
      /* DD state 4 extra transitions */
      if (apair == NONE) {
	sc = p7_FLogsum(p7_FLogsum(SSMXM(xv,kv) + e2hmmerTSC(kv,e2HP_SS_DD), 
				   SDMXM(xv,kv) + e2hmmerTSC(kv,e2HP_SD_DD)),
			p7_FLogsum(DSMXM(xv,kv) + e2hmmerTSC(kv,e2HP_DS_DD),
				   DDMXM(xv,kv) + e2hmmerTSC(kv,e2HP_DD_DD)));
	sc = p7_FLogsum(DDMXM(x,k),sc);
	DDMXM(x,k) = p7_FLogsum(sc, DDMXM(xv,k));
      }

      /* IB state 2 transitions */
      IBMXM(x,k) = -eslINFINITY;
      if (apair == LRES) {
	sc = p7_FLogsum(BBMXM(xv,k) + e2hmmerTSC(k,e2HP_BB_IB), 
			IBMXM(xv,k) + e2hmmerTSC(k,e2HP_IB_IB));
	IBMXM(x,k) = sc + Ilsc;
      }
      if (apair == NONE) IBMXM(x,k) = IBMXM(xv,k);
      
      /* IS state 2 transitions  */
      ISMXM(x,k) = -eslINFINITY;
      if (apair == LRES) {
	sc = p7_FLogsum(SSMXM(xv,k) + e2hmmerTSC(k,e2HP_SS_IS),
			ISMXM(xv,k) + e2hmmerTSC(k,e2HP_IS_IS));
	ISMXM(x,k) = sc + Ilsc;
      }
      if (apair == NONE) ISMXM(x,k) = ISMXM(xv,k);

      /* ID state 2 transitions */
      IDMXM(x,k) = -eslINFINITY;
      if (apair == LRES) {
	sc = p7_FLogsum(SDMXM(xv,k) + e2hmmerTSC(k,e2HP_SD_ID),
			IDMXM(xv,k) + e2hmmerTSC(k,e2HP_ID_ID));
	IDMXM(x,k) = sc + Ilsc;
      }
      if (apair == NONE) IDMXM(x,k) = IDMXM(xv,k);
      
      /* BI state 2 transitions */
      BIMXM(x,k) = -eslINFINITY;
      if (apair ==  RRES) {
	sc = p7_FLogsum(BBMXM(xv,k) + e2hmmerTSC(k,e2HP_BB_BI), 
			BIMXM(xv,k) + e2hmmerTSC(k,e2HP_BI_BI));
	BIMXM(x,k) = sc + Irsc;
     }
      if (apair == NONE) BIMXM(x,k) =  BIMXM(xv,k);

      /* SI state 2 transitions */
      SIMXM(x,k) = -eslINFINITY;
      if (apair == RRES) {
	sc = p7_FLogsum(SSMXM(xv,k) + e2hmmerTSC(k,e2HP_SS_SI),
		        SIMXM(xv,k) + e2hmmerTSC(k,e2HP_SI_SI));
	SIMXM(x,k) = sc + Irsc;
      }
      if (apair == NONE) SIMXM(x,k) = SIMXM(xv,k);

     /* DI state 2 transitions */
      DIMXM(x,k) = -eslINFINITY;
      if (apair == RRES) {
	sc = p7_FLogsum(DSMXM(xv,k) + e2hmmerTSC(k,e2HP_DS_DI),
			DIMXM(xv,k) + e2hmmerTSC(k,e2HP_DI_DI));
	DIMXM(x,k) = sc + Irsc;
      }
      if (apair == NONE) DIMXM(x,k) = DIMXM(xv,k);

      /* Ii state 4 transitions */
      IiMXM(x,k) = -eslINFINITY;
      if (apair == LRES) {
	sc = p7_FLogsum(p7_FLogsum(BIMXM(xv,k) + e2hmmerTSC(k, e2HP_BI_Ii),
				   SIMXM(xv,k) + e2hmmerTSC(k, e2HP_SI_Ii)),
			p7_FLogsum(IiMXM(xv,k) + e2hmmerTSC(k, e2HP_Ii_Ii),
				   iIMXM(xv,k) + e2hmmerTSC(k, e2HP_iI_Ii)));
	IiMXM(x,k) = sc + Ilsc;
      }
      if (apair == NONE) IiMXM(x,k) = IiMXM(xv,k);

      /* iI state 4 transitions */
      iIMXM(x,k) = -eslINFINITY;
      if (apair == RRES) {
	sc = p7_FLogsum(p7_FLogsum(IBMXM(xv,k) + e2hmmerTSC(k,e2HP_IB_iI),
				   ISMXM(xv,k) + e2hmmerTSC(k,e2HP_IS_iI)),
			p7_FLogsum(IiMXM(xv,k) + e2hmmerTSC(k,e2HP_Ii_iI),
				   iIMXM(xv,k) + e2hmmerTSC(k,e2HP_iI_iI)));
	iIMXM(x,k) = sc + Irsc;
      }
      if (apair == NONE) iIMXM(x,k) = iIMXM(xv,k);
      
      /* EE state 13 transitions */
      EEMXM(x,k) = -eslINFINITY;
      if (x == L) {
	sc = p7_FLogsum(p7_FLogsum(BBMXM(x,k) + e2hmmerTSC(k,e2HP_BB_EE), 
				   IBMXM(x,k) + e2hmmerTSC(k,e2HP_IB_EE)),
			p7_FLogsum(SSMXM(x,k) + e2hmmerTSC(k,e2HP_SS_EE),
				   DSMXM(x,k) + e2hmmerTSC(k,e2HP_DS_EE)));
	sc = p7_FLogsum(p7_FLogsum(sc, 
				   ISMXM(x,k) + e2hmmerTSC(k,e2HP_IS_EE)),
			p7_FLogsum(SDMXM(x,k) + e2hmmerTSC(k,e2HP_SD_EE), 
				   DDMXM(x,k) + e2hmmerTSC(k,e2HP_DD_EE)));
	sc = p7_FLogsum(p7_FLogsum(sc, 
				   IDMXM(x,k) + e2hmmerTSC(k,e2HP_ID_EE)),
			p7_FLogsum(BIMXM(x,k) + e2hmmerTSC(k,e2HP_BI_EE), 
				   SIMXM(x,k) + e2hmmerTSC(k,e2HP_SI_EE)));
	sc = p7_FLogsum(p7_FLogsum(sc,
				   DIMXM(x,k) + e2hmmerTSC(k,e2HP_DI_EE)), 
			p7_FLogsum(IiMXM(x,k) + e2hmmerTSC(k,e2HP_Ii_EE), 
				   iIMXM(x,k) + e2hmmerTSC(k,e2HP_iI_EE)));
  	EEMXM(x,k) = sc;
	EE = p7_FLogsum(EE, EEMXM(x,k));
      }
      
#if 0
      if (x==L)printf("e2fFWD x %d k %d/%d apair %d SSsc %f BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f iI %f  EE %f\n", 
			      x, k, M, apair, SSsc,
		      BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), 
		      SDMXM(x,k), DDMXM(x,k), IDMXM(x,k), 
		      BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k), EEMXM(x,k));
#endif
   }
  }
  
 if (opt_sc != NULL) *opt_sc = EE;
  
  gx->Lrow = L;
  gx->Lcol = L;
  
  return eslOK;
}

/* Function:  e2hmmer_GBackward()
 * Synopsis:  The Backward algorithm.
 *
 * Purpose:   The Backward dynamic programming algorithm.
 * 
 *            Given two profile sequences <psql> and <psqr>,  of
 *            equal length L = <psql->n> = <psqr->n>, a 
 *            model profile <gm>, and DP matrix <gx> allocated
 *            in linear time for L + 3 cells;
 *            calculate the probability of the sequence
 *            given the model using the Backward algorithm; return the
 *            Backward matrix in <gx>, and the Backward score in <ret_sc>.
 *           
  *           
 *            The Backward score is in lod score form. 
 *
 * Args:      psql   - profile sequence, 1..L
 *            psqr   - profile sequence, 1..L
 *            gm     - e2 model profile. 
 *            gx     - DP matrix with linear memory allocation
 *            opt_sc - optRETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
e2hmmer_GBackward(const PSQ *psql, const PSQ *psqr, const E2HMMER_PROFILE *gm, E2_GMX *gx, float *opt_sc)
{
  float const **tsc = (float const **)gm->tsc;
  float       **dp  = gx->dp;						
  float         sc;
  float         SSsc;
  float         Slsc, Srsc;
  float         Ilsc, Irsc;
  float        *pl = NULL;
  float        *pr = NULL;
  float         BB = -eslINFINITY;
  int           L;  
  int           K = psql->abc->K;
  int           isgapl, isgapr;
  enum apair_e  apair;
  int           b,c;                  /* alphabet index */
  int           x, xv;
  int           k, kv;
  
  /* Note: backward calculates the probability we can get *out* of
   * a cell(i,j); exclusive of emitting residue x_i or y_j.
   */
  p7_FLogsumInit();    
  
  /* Backward:
   */
  /* Initialize x = L && k = M */
  x = L;
  k = gm->M;

  /* EE state */
  EEMXM(x,k) = 0.0;
  /* BB state */
  BBMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_BB_EE);  
  /* IB state */
  IBMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_IB_EE);  
  /* SS state */
  SSMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_SS_EE);  
  /* DS state  */
  DSMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_DS_EE);	
  /* IS state  */
  ISMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_IS_EE);	
  /* SD state */
  SDMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_SD_EE);       
  /* DD state */
  DDMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_DD_EE);
  /* ID state  */
  IDMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_ID_EE);	
  /* BI state */
  BIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_BI_EE);
  /* SI state */
  SIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_SI_EE);	
  /* DI state  */
  DIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_DI_EE);
  /* Ii state */
  IiMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_Ii_EE);
  /* iI state */
  iIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_iI_EE);
#if 0
  if (x==91) printf("\ne2fBCW x %d  k %d BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f iI %f EE %f\n", x, k, 
	   BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), SDMXM(x,k), 
	   DDMXM(x,k), IDMXM(x,k), BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k), EEMXM(x,k));
#endif    
 
 
  /* Initialize x = L && k < M */
  x = L;
  for (k = gm->M-1; k >= 0; k --) {
    kv = k + 1;
        
    /* EE state  */
    EEMXM(x,k) = 0.0;
    
    /* BB state 4-3=1 transitions */
    BBMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k, e2HP_BB_EE); 
 					 
    /* IB state 4-3=1 transitions */
    IBMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_IB_EE);
	
    /* SS state (7+1)-6=2 transitions */
    SSMXM(x,k) = p7_FLogsum(EEMXM(x,k)  + e2hmmerTSC(k,e2HP_SS_EE), 
			    DDMXM(x,kv) + e2hmmerTSC(k,e2HP_SS_DD));
	
    /* DS state (6+1)-5=2 transitions */
    DSMXM(x,k) = p7_FLogsum(EEMXM(x,k)  + e2hmmerTSC(k, e2HP_DS_EE), 
			    DDMXM(x,kv) + e2hmmerTSC(k, e2HP_DS_DD));
	
    /* IS state 5-4=1 transitions */
    ISMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_IS_EE);	
    
    /* SD state (6+1)-5=2 transitions */
    SDMXM(x,k) = p7_FLogsum(EEMXM(x,k)  + e2hmmerTSC(k, e2HP_SD_EE), 
			    DDMXM(x,kv) + e2hmmerTSC(k, e2HP_SD_DD));
    
    /* DD state (5+1)-4=2 transitions */
    DDMXM(x,k) = p7_FLogsum(EEMXM(x,k)  + e2hmmerTSC(k,e2HP_DD_EE), 
			    DDMXM(x,kv) + e2hmmerTSC(k,e2HP_DD_DD));
	
    /* ID state 4-3=1 transitions */
    IDMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_ID_EE);
    
    /* BI state 4-3=1 transitions */
    BIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_BI_EE);
	
    /* SI state 5-4=1 transitions */
    SIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k, e2HP_SI_EE);
    
    /* DI state 4-3=1 transitions */
    DIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_DI_EE);
    
    /* Ii state 4-3=1 transitions */
    IiMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_Ii_EE); 

    /* iI state 4-3=1 transitions */
    iIMXM(x,k) = EEMXM(x,k) + e2hmmerTSC(k,e2HP_iI_EE); 
    
#if 0
    if (x==0)printf("e2fBCW x %d  k %d BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f iI %f EE %f\n", x, k, 
	   BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), SDMXM(x,k), 
	   DDMXM(x,k), IDMXM(x,k), BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k), EEMXM(x,k));
#endif    
  }

  /* Initialize k = M && x < L */
  k = gm->M;
  for (x = L-1; x >= 0; x--) {

    xv = x + 1;

    pl = psql->prof[xv];
    pr = psqr->prof[xv];
    SSsc = -eslINFINITY;
    Slsc = -eslINFINITY; /* orphan substitution score */
    Srsc = -eslINFINITY; /* orphan substitution score */
    Ilsc = -eslINFINITY; /* insertion score */
    Irsc = -eslINFINITY; /* insertion score */
    for (b = 0; b < K; b ++) {
      Ilsc = p7_FLogsum(Ilsc, gm->isc[k][b][e2HP_SL] + pl[b]);
      Irsc = p7_FLogsum(Irsc, gm->isc[k][b][e2HP_SR] + pr[b]);
    }
    isgapl = isgap(pl, K);
    isgapr = isgap(pr, K);
    
    if      (!isgapl && !isgapr) apair = PAIR;
    else if ( isgapl && !isgapr) apair = RRES;
    else if (!isgapl &&  isgapr) apair = LRES;
    else                         apair = NONE;
    
   /* EE state  */
    EEMXM(x,k) = -eslINFINITY;
    
    /* BB state 4-1=3 transitions */
    sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IBMXM(xv,k) + e2hmmerTSC(k,e2HP_BB_IB) + Ilsc : -eslINFINITY,
			       (apair == RRES)? BIMXM(xv,k) + e2hmmerTSC(k,e2HP_BB_BI) + Irsc : -eslINFINITY),
		                                EEMXM(x, k) + e2hmmerTSC(k,e2HP_BB_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? BBMXM(xv,k) : -eslINFINITY));
    BBMXM(x,k) = sc;
 	
    /* IB state 4-1=3 transitions */
    sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IBMXM(xv,k) + e2hmmerTSC(k,e2HP_IB_IB) + Ilsc : -eslINFINITY, 
			       (apair == RRES)? iIMXM(xv,k) + e2hmmerTSC(k,e2HP_IB_iI) + Irsc : -eslINFINITY),
		                                EEMXM(x, k) + e2hmmerTSC(k,e2HP_IB_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? IBMXM(xv,k) : -eslINFINITY));
    IBMXM(x,k) = sc;
	
   /* SS state (7+1)-5=3 transitions */
    sc = p7_FLogsum(p7_FLogsum((apair == LRES)? ISMXM(xv,k) + e2hmmerTSC(k,e2HP_SS_IS) + Ilsc : -eslINFINITY,
			       (apair == RRES)? SIMXM(xv,k) + e2hmmerTSC(k,e2HP_SS_SI) + Irsc : -eslINFINITY),
		                                EEMXM(x, k) + e2hmmerTSC(k,e2HP_SS_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? SSMXM(xv,k) : -eslINFINITY));
    SSMXM(x,k) = sc;
	
    /* DS state (6+1)-5=2 transitions */
    sc = p7_FLogsum((apair == RRES)? DIMXM(xv,k) + e2hmmerTSC(k, e2HP_DS_DI) + Irsc : -eslINFINITY,
		                     EEMXM(x, k) + e2hmmerTSC(k, e2HP_DS_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? DSMXM(xv,k) : -eslINFINITY));
    DSMXM(x,k) = sc;
	
    /* IS state 5-2=3 transitions */
    sc = p7_FLogsum((apair == LRES)? ISMXM(xv,k) + e2hmmerTSC(k,e2HP_IS_IS) + Ilsc : -eslINFINITY,
		    (apair == RRES)? iIMXM(xv,k) + e2hmmerTSC(k,e2HP_IS_iI) + Irsc : -eslINFINITY);
    sc = p7_FLogsum(sc, ((apair == NONE)? ISMXM(xv,k) : -eslINFINITY));
    ISMXM(x,k) = p7_FLogsum(sc, EEMXM(x,k) + e2hmmerTSC(k,e2HP_IS_EE));	
    
    /* SD state (6+1)-5=2 transitions */
    sc = p7_FLogsum((apair == LRES)? IDMXM(xv,k) + e2hmmerTSC(k, e2HP_SD_ID) + Ilsc : -eslINFINITY,
		                     EEMXM(x, k) + e2hmmerTSC(k, e2HP_SD_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? SDMXM(xv,k) : -eslINFINITY));
    SDMXM(x,k) = sc;
	
    /* DD state (5+1)-5=1 transitions */
    sc = p7_FLogsum(EEMXM(x,k) + e2hmmerTSC(k,e2HP_DD_EE), ((apair == NONE)? DDMXM(xv,k) : -eslINFINITY));
    DDMXM(x,k) = sc;
	
    /* ID state 4-2=2 transitions */
    sc = p7_FLogsum((apair == LRES)? IDMXM(xv,k) + e2hmmerTSC(k,e2HP_ID_ID) + Ilsc : -eslINFINITY,
		                     EEMXM(x, k) + e2hmmerTSC(k,e2HP_ID_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? IDMXM(xv,k) : -eslINFINITY));
    IDMXM(x,k) = sc;
    
    /* BI state 4-1=3 transitions */
    sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_BI_Ii) + Ilsc : -eslINFINITY,
			       (apair == RRES)? BIMXM(xv,k)  + e2hmmerTSC(k,e2HP_BI_BI) + Irsc : -eslINFINITY),
		                                EEMXM(x, k)  + e2hmmerTSC(k,e2HP_BI_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? BIMXM(xv,k) : -eslINFINITY));
    BIMXM(x,k) = sc;
	
    /* SI state 5-2=3 transitions */
    sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_SI_Ii) + Ilsc : -eslINFINITY,
			       (apair == RRES)? SIMXM(xv,k)  + e2hmmerTSC(k,e2HP_SI_SI) + Irsc : -eslINFINITY),
		                                EEMXM(x, k)  + e2hmmerTSC(k,e2HP_SI_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? SIMXM(xv,k) : -eslINFINITY));
    SIMXM(x,k) = sc;
	
    /* DI state 4-2=2 transitions */
    sc = p7_FLogsum((apair == RRES)? DIMXM(xv,k) + e2hmmerTSC(k,e2HP_DI_DI) + Irsc : -eslINFINITY,
		                     EEMXM(x, k) + e2hmmerTSC(k,e2HP_DI_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? DIMXM(xv,k) : -eslINFINITY));
    DIMXM(x,k) = sc;
	
    /* Ii state 4-1=3 transitions */
    sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_Ii_Ii) + Ilsc : -eslINFINITY,
			       (apair == RRES)? iIMXM(xv,k)  + e2hmmerTSC(k,e2HP_Ii_iI) + Irsc : -eslINFINITY),
		                                EEMXM(x, k)  + e2hmmerTSC(k,e2HP_Ii_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? IiMXM(xv,k) : -eslINFINITY));
    IiMXM(x,k) = sc;
 
    /* iI state 4-1=3 transitions */
    sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_iI_Ii) + Ilsc : -eslINFINITY,
			       (apair == RRES)? iIMXM(xv,k)  + e2hmmerTSC(k,e2HP_iI_iI) + Irsc : -eslINFINITY),
		                                EEMXM(x, k)  + e2hmmerTSC(k,e2HP_iI_EE));
    sc = p7_FLogsum(sc, ((apair == NONE)? iIMXM(xv,k) : -eslINFINITY));
    iIMXM(x,k) = sc;
    
    if (x == 0) BB = p7_FLogsum(BB, BBMXM(x,k));

 #if 0
    if (x==0)printf("e2fBCW x %d  k %d pair %d BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f iI %f EE %f\n", x, k, apair,
	   BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), SDMXM(x,k), 
	   DDMXM(x,k), IDMXM(x,k), BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k), EEMXM(x,k));
#endif    
  }

  /* Main recursion x < L && k < M */
  for (k = gm->M-1; k >= 0; k --) {
    kv = k + 1;

    for (x = L-1; x >= 0; x--)
      {
	/* ahead one */
	xv = x + 1;

	pl = psql->prof[xv];
	pr = psqr->prof[xv];
	SSsc = -eslINFINITY;
	Slsc = -eslINFINITY; /* orphan substitution score */
	Srsc = -eslINFINITY; /* orphan substitution score */
	Ilsc = -eslINFINITY; /* insertion score */
	Irsc = -eslINFINITY; /* insertion score */
	for (b = 0; b < K; b ++) {
	  Slsc = p7_FLogsum(Slsc, gm->ssc[kv][b][e2HP_SL] + pl[b]);
	  Ilsc = p7_FLogsum(Ilsc, gm->isc[k] [b][e2HP_SL] + pl[b]);
	  Srsc = p7_FLogsum(Srsc, gm->ssc[kv][b][e2HP_SR] + pr[b]);
	  Irsc = p7_FLogsum(Irsc, gm->isc[k] [b][e2HP_SR] + pr[b]);
	  for (c = 0; c < K; c ++) {
	    SSsc = p7_FLogsum(SSsc, gm->sssc[kv][b][c] + pl[b] + pr[c]);
	  }
	}
	isgapl = isgap(pl, K);
	isgapr = isgap(pr, K);
	
	if      (!isgapl && !isgapr) apair = PAIR;
	else if ( isgapl && !isgapr) apair = RRES;
	else if (!isgapl &&  isgapr) apair = LRES;
	else                         apair = NONE;
	
	/* EE state  */
	EEMXM(x,k) = -eslINFINITY;
    
	/* BB state 4 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IBMXM(xv,k)  + e2hmmerTSC(k,e2HP_BB_IB) + Ilsc : -eslINFINITY,
				   (apair == RRES)? BIMXM(xv,k)  + e2hmmerTSC(k,e2HP_BB_BI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_BB_SS) + SSsc : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_BB_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? BBMXM(xv,k) : -eslINFINITY));
	BBMXM(x,k) = sc;
	
	/* IB state 4 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IBMXM(xv,k)  + e2hmmerTSC(k,e2HP_IB_IB) + Ilsc : -eslINFINITY, 
				   (apair == RRES)? iIMXM(xv,k)  + e2hmmerTSC(k,e2HP_IB_iI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_IB_SS) + SSsc : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_IB_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? IBMXM(xv,k) : -eslINFINITY));
	IBMXM(x,k) = sc;
	
	/* SS state 7+1 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? ISMXM(xv,k)  + e2hmmerTSC(k,e2HP_SS_IS) + Ilsc : -eslINFINITY,
				   (apair == RRES)? SIMXM(xv,k)  + e2hmmerTSC(k,e2HP_SS_SI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_SS_SS) + SSsc : -eslINFINITY,
				   (apair == LRES)? SDMXM(xv,kv) + e2hmmerTSC(k,e2HP_SS_SD) + Slsc : -eslINFINITY));
	sc = p7_FLogsum(p7_FLogsum(sc,
				   (apair == RRES)? DSMXM(xv,kv) + e2hmmerTSC(k,e2HP_SS_DS) + Srsc : -eslINFINITY),
			p7_FLogsum((apair == NONE)? DDMXM(xv,kv) + e2hmmerTSC(k,e2HP_SS_DD)        : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_SS_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? SSMXM(xv,k) : -eslINFINITY));
	SSMXM(x,k) = p7_FLogsum(sc, DDMXM(x,kv) + e2hmmerTSC(k,e2HP_SS_DD));
	
	/* DS state 6+1 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == RRES)? DIMXM(xv,k)  + e2hmmerTSC(k, e2HP_DS_DI) + Irsc : -eslINFINITY,
				   (apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k, e2HP_DS_SS) + SSsc : -eslINFINITY),
			p7_FLogsum((apair == LRES)? SDMXM(xv,kv) + e2hmmerTSC(k, e2HP_DS_SD) + Slsc : -eslINFINITY,
				   (apair == RRES)? DSMXM(xv,kv) + e2hmmerTSC(k, e2HP_DS_DS) + Srsc : -eslINFINITY));
	sc = p7_FLogsum(sc,
			p7_FLogsum((apair == NONE)? DDMXM(xv,kv) + e2hmmerTSC(k, e2HP_DS_DD)        : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k, e2HP_DS_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? DSMXM(xv,k) : -eslINFINITY));
	DSMXM(x,k) = p7_FLogsum(sc, DDMXM(x,kv) + e2hmmerTSC(k, e2HP_DS_DD));
	
	/* IS state 5 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? ISMXM(xv,k)  + e2hmmerTSC(k,e2HP_IS_IS) + Ilsc : -eslINFINITY,
				   (apair == RRES)? iIMXM(xv,k)  + e2hmmerTSC(k,e2HP_IS_iI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_IS_SS) + SSsc : -eslINFINITY,
				   (apair == LRES)? SDMXM(xv,kv) + e2hmmerTSC(k,e2HP_IS_SD) + Slsc : -eslINFINITY));
	sc = p7_FLogsum(sc, ((apair == NONE)? ISMXM(xv,k) : -eslINFINITY));
	ISMXM(x,k) = p7_FLogsum(sc, EEMXM(x,k) + e2hmmerTSC(k,e2HP_IS_EE));	
	
	/* SD state 6+1 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IDMXM(xv,k)  + e2hmmerTSC(k, e2HP_SD_ID) + Ilsc : -eslINFINITY,
				   (apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k, e2HP_SD_SS) + SSsc : -eslINFINITY),
			p7_FLogsum((apair == LRES)? SDMXM(xv,kv) + e2hmmerTSC(k, e2HP_SD_SD) + Slsc : -eslINFINITY,
				   (apair == RRES)? DSMXM(xv,kv) + e2hmmerTSC(k, e2HP_SD_DS) + Srsc : -eslINFINITY));
	sc = p7_FLogsum(sc,
			p7_FLogsum((apair == NONE)? DDMXM(xv,kv) + e2hmmerTSC(k, e2HP_SD_DD)        : -eslINFINITY,
				                    EEMXM(x,k)   + e2hmmerTSC(k, e2HP_SD_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? SDMXM(xv,k) : -eslINFINITY));
	SDMXM(x,k) = p7_FLogsum(sc, DDMXM(x,kv) + e2hmmerTSC(k, e2HP_SD_DD));
	
	/* DD state 5+1 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_DD_SS) + SSsc : -eslINFINITY,
				   (apair == LRES)? SDMXM(xv,kv) + e2hmmerTSC(k,e2HP_DD_SD) + Slsc : -eslINFINITY),
			p7_FLogsum((apair == RRES)? DSMXM(xv,kv) + e2hmmerTSC(k,e2HP_DD_DS) + Srsc : -eslINFINITY,
				   (apair == NONE)? DDMXM(xv,kv) + e2hmmerTSC(k,e2HP_DD_DD)        : -eslINFINITY));
	sc = p7_FLogsum(sc, 
			p7_FLogsum(DDMXM(x,kv) + e2hmmerTSC(k,e2HP_DD_DD),
				   EEMXM(x,k)  + e2hmmerTSC(k,e2HP_DD_EE)));
 
	DDMXM(x,k) = p7_FLogsum(sc, ((apair == NONE)? DDMXM(xv,k) : -eslINFINITY));
	
	/* ID state 4 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IDMXM(xv,k)  + e2hmmerTSC(k,e2HP_ID_ID) + Ilsc : -eslINFINITY,
				   (apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_ID_SS) + SSsc : -eslINFINITY),
			p7_FLogsum((apair == LRES)? SDMXM(xv,kv) + e2hmmerTSC(k,e2HP_ID_SD) + Slsc : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_ID_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? IDMXM(xv,k) : -eslINFINITY));
	IDMXM(x,k) = sc;
	
	/* BI state 4 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_BI_Ii) + Ilsc : -eslINFINITY,
				   (apair == RRES)? BIMXM(xv,k)  + e2hmmerTSC(k,e2HP_BI_BI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_BI_SS) + SSsc : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_BI_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? BIMXM(xv,k) : -eslINFINITY));
	BIMXM(x,k) = sc;
	
	/* SI state 5 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_SI_Ii) + Ilsc : -eslINFINITY,
				   (apair == RRES)? SIMXM(xv,k)  + e2hmmerTSC(k,e2HP_SI_SI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_SI_SS) + SSsc : -eslINFINITY,
				   (apair == RRES)? DSMXM(xv,kv) + e2hmmerTSC(k,e2HP_SI_DS) + Srsc : -eslINFINITY));
	sc = p7_FLogsum(sc, 
			p7_FLogsum((apair == NONE)? SIMXM(xv,k) : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_SI_EE)));
	SIMXM(x,k) = sc;
	
	/* DI state 4 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == RRES)? DIMXM(xv,k)  + e2hmmerTSC(k,e2HP_DI_DI) + Irsc : -eslINFINITY,
				   (apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_DI_SS) + SSsc : -eslINFINITY),
			p7_FLogsum((apair == RRES)? DSMXM(xv,kv) + e2hmmerTSC(k,e2HP_DI_DS) + Srsc : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_DI_EE)));
	sc = p7_FLogsum(sc, ((apair == NONE)? DIMXM(xv,k) : -eslINFINITY));
	DIMXM(x,k) = sc;
	
	/* Ii state 4 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_Ii_Ii) + Ilsc : -eslINFINITY,
				   (apair == RRES)? iIMXM(xv,k)  + e2hmmerTSC(k,e2HP_Ii_iI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_Ii_SS) + SSsc : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_Ii_EE)));
	sc = p7_FLogsum(sc, (apair == NONE)? IiMXM(xv,k) : -eslINFINITY);
	IiMXM(x,k) = sc;
	
	/* iI state 4 transitions */
	sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IiMXM(xv,k)  + e2hmmerTSC(k,e2HP_iI_Ii) + Ilsc : -eslINFINITY,
				   (apair == RRES)? iIMXM(xv,k)  + e2hmmerTSC(k,e2HP_iI_iI) + Irsc : -eslINFINITY),
			p7_FLogsum((apair == PAIR)? SSMXM(xv,kv) + e2hmmerTSC(k,e2HP_iI_SS) + SSsc : -eslINFINITY,
				                    EEMXM(x, k)  + e2hmmerTSC(k,e2HP_iI_EE)));
	sc = p7_FLogsum(sc, (apair == NONE)? iIMXM(xv,k) : -eslINFINITY);
	iIMXM(x,k) = sc;
	
	if (x == 0) BB = p7_FLogsum(BB, BBMXM(x,k));

 #if 0
	if (x==0) printf("e2fBCW x %d  k %d pair %d BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f iI %f EE %f\n", 
			 x, k, apair, 
			 BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), SDMXM(x,k), 
			 DDMXM(x,k), IDMXM(x,k), BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k), EEMXM(x,k));
#endif
      }

  }

 
  if (opt_sc != NULL) *opt_sc = BB;
  
  gx->Lrow = L;
  gx->Lcol = L;
  
  return eslOK;
}


int 
isgap(float *p, int K) 
{
  if (p[K] == 0.) return 1;
  else return 0;
}
