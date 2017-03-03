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
#include "e2f_generic_fwdback.h"
#include "e2_profilesq.h"


/*****************************************************************
 * 1. Forward, Backward implementations.
 *****************************************************************/
static int isgap(float *p, int K);

/* Function:  e2f_GForward()
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
e2f_GForward(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, E2_GMX *gx, float *opt_sc)
{
  float const  *tsc  = (float const  *)gm->tsc;
  float       **dp   = gx->dp;						
  float         sc;
  float         SSsc;
  float         Slsc, Srsc;
  float         Ilsc, Irsc;
  float        *pl = NULL;
  float        *pr = NULL;
  int           L;
  int           K = psql->abc->K;
  int           isgapl, isgapr;
  enum apair_e  apair;
  int           b,c;                  /* alphabet index */
  int           x, xv;                /* linear memor indices */
  
  if (psql->n != psqr->n) { printf("e2f_GForward(): sqs are not aligned\n"); return eslFAIL; }
  L = psql->n;

  p7_FLogsumInit();    
  
#if 0
  ddfactor = -log(1.0 - bg->p * exp(tsc[e2P_DD_DD]));
#endif

  /* Forward:
   *           states DD, EE need to be evaluated last.
   *
   * Order: BB, SS, DS, SD, IB, IS, ID, BI, SI, II, DD, EE
   */
  /* Initialization of the zero row. */
  x = 0;
  BBMX(x) = 0.0;
  SSMX(x) = DSMX(x) = SDMX(x) = -eslINFINITY;
  IBMX(x) = ISMX(x) = IDMX(x) = -eslINFINITY;
  BIMX(x) = SIMX(x) = DIMX(x) = IIMX(x) = -eslINFINITY;
  DDMX(x) = -eslINFINITY; 
  E2G_XMX(gx, x, e2G_EE) = -eslINFINITY;
  
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
      Slsc = p7_FLogsum(Slsc, gm->ssc[b][e2P_SL] + pl[b]);
      Ilsc = p7_FLogsum(Ilsc, gm->isc[b][e2P_SL] + pl[b]);
      Srsc = p7_FLogsum(Srsc, gm->ssc[b][e2P_SR] + pr[b]);
      Irsc = p7_FLogsum(Irsc, gm->isc[b][e2P_SR] + pr[b]);
      for (c = 0; c < K; c ++) {
	SSsc = p7_FLogsum(SSsc, gm->sssc[b][c] + pl[b] + pr[c]);
      }
    }
    isgapl = isgap(pl, K);
    isgapr = isgap(pr, K);
    if      (!isgapl && !isgapr) apair = PAIR;
    else if ( isgapl && !isgapr) apair = RRES;
    else if (!isgapl &&  isgapr) apair = LRES;
    else                         apair = NONE;

    /* BB state 0 transitions  */
    BBMX(x) = -eslINFINITY;
    /* SS state 12 transitions */
    SSMX(x) = -eslINFINITY;
    if (apair == PAIR) { 
      sc = p7_FLogsum(p7_FLogsum(BBMX(xv) + e2TSC(e2P_BB_SS), 
				 IBMX(xv) + e2TSC(e2P_IB_SS)),
		      p7_FLogsum(SSMX(xv) + e2TSC(e2P_SS_SS),
				 DSMX(xv) + e2TSC(e2P_DS_SS)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 ISMX(xv) + e2TSC(e2P_IS_SS)),
		      p7_FLogsum(SDMX(xv) + e2TSC(e2P_SD_SS), 
				 DDMX(xv) + e2TSC(e2P_DD_SS)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 IDMX(xv) + e2TSC(e2P_ID_SS)),
		      p7_FLogsum(BIMX(xv) + e2TSC(e2P_BI_SS), 
				 SIMX(xv) + e2TSC(e2P_SI_SS)));
      sc = p7_FLogsum(sc, 
		      p7_FLogsum(DIMX(xv) + e2TSC(e2P_DI_SS), 
				 IIMX(xv) + e2TSC(e2P_II_SS)));
      SSMX(x) = sc + SSsc;
    }
    /* DS state 12 transitions */
    DSMX(x) = -eslINFINITY;
    if (apair == RRES) {
      sc = p7_FLogsum(p7_FLogsum(BBMX(xv) + e2TSC(e2P_BB_DS), 
				 IBMX(xv) + e2TSC(e2P_IB_DS)),
		      p7_FLogsum(SSMX(xv) + e2TSC(e2P_SS_DS),
				 DSMX(xv) + e2TSC(e2P_DS_DS)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 ISMX(xv) + e2TSC(e2P_IS_DS)),
		      p7_FLogsum(SDMX(xv) + e2TSC(e2P_SD_DS), 
				 DDMX(xv) + e2TSC(e2P_DD_DS)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 IDMX(xv) + e2TSC(e2P_ID_DS)),
		      p7_FLogsum(BIMX(xv) + e2TSC(e2P_BI_DS), 
				 SIMX(xv) + e2TSC(e2P_SI_DS)));
      sc = p7_FLogsum(sc, 
		      p7_FLogsum(DIMX(xv) + e2TSC(e2P_DI_DS), 
				 IIMX(xv) + e2TSC(e2P_II_DS)));
      DSMX(x) = sc + Srsc;
    }
    /* SD state 12 transitions */
    SDMX(x) = -eslINFINITY;
    if (apair == LRES) {
      sc = p7_FLogsum(p7_FLogsum(BBMX(xv) + e2TSC(e2P_BB_SD), 
				 IBMX(xv) + e2TSC(e2P_IB_SD)),
		      p7_FLogsum(SSMX(xv) + e2TSC(e2P_SS_SD),
				 DSMX(xv) + e2TSC(e2P_DS_SD)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 ISMX(xv) + e2TSC(e2P_IS_SD)),
		      p7_FLogsum(SDMX(xv) + e2TSC(e2P_SD_SD), 
				 DDMX(xv) + e2TSC(e2P_DD_SD)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 IDMX(xv) + e2TSC(e2P_ID_SD)),
		      p7_FLogsum(BIMX(xv) + e2TSC(e2P_BI_SD), 
				 SIMX(xv) + e2TSC(e2P_SI_SD)));
      sc = p7_FLogsum(sc, 
		      p7_FLogsum(DIMX(xv) + e2TSC(e2P_DI_SD), 
				 IIMX(xv) + e2TSC(e2P_II_SD)));
      SDMX(x) = sc + Slsc;
    }
    
    /* IB state 2 transitions */
    IBMX(x) = -eslINFINITY;
    if (apair == LRES) {
      sc = p7_FLogsum(BBMX(xv) + e2TSC(e2P_BB_IB), 
		      IBMX(xv) + e2TSC(e2P_IB_IB));
      IBMX(x) = sc + Ilsc;
    }
    
    /* IS state 3 transitions  */
   ISMX(x) = -eslINFINITY;
   if (apair == LRES) {
     sc = p7_FLogsum(           SSMX(xv) + e2TSC(e2P_SS_IS),
				p7_FLogsum(DSMX(xv) + e2TSC(e2P_DS_IS),
					   ISMX(xv) + e2TSC(e2P_IS_IS)));
     ISMX(x) = sc + Ilsc;
    }
    /* ID state 3 transitions */
   IDMX(x) = -eslINFINITY;
   if (apair == LRES) {
     sc = p7_FLogsum(           SDMX(xv) + e2TSC(e2P_SD_ID),
				p7_FLogsum(DDMX(xv) + e2TSC(e2P_DD_ID),
					   IDMX(xv) + e2TSC(e2P_ID_ID)));
     IDMX(x) = sc + Ilsc;
   }
    
    /* BI state 2 transitions */
   BIMX(x) = -eslINFINITY;
   if (apair ==  RRES) {
     sc = p7_FLogsum(BBMX(xv) + e2TSC(e2P_BB_BI), 
		     BIMX(xv) + e2TSC(e2P_BI_BI));
     BIMX(x) = sc + Irsc;
   }
   /* SI state 3 transitions */
   SIMX(x) = -eslINFINITY;
   if (apair == RRES) {
     sc = p7_FLogsum(           SSMX(xv) + e2TSC(e2P_SS_SI),
				p7_FLogsum(SDMX(xv) + e2TSC(e2P_SD_SI),
					   SIMX(xv) + e2TSC(e2P_SI_SI)));
     SIMX(x) = sc + Irsc;
   }
   /* DI state 3 transitions */
   DIMX(x) = -eslINFINITY;
   if (apair == RRES) {
     sc = p7_FLogsum(           DSMX(xv) + e2TSC(e2P_DS_DI),
				p7_FLogsum(DDMX(xv) + e2TSC(e2P_DD_DI),
					   DIMX(xv) + e2TSC(e2P_DI_DI)));
     DIMX(x) = sc + Irsc;
   }
   /* II state 4 transitions */
   IIMX(x) = -eslINFINITY;
   if (apair == RRES) {
     sc = p7_FLogsum(p7_FLogsum(IBMX(xv) + e2TSC(e2P_IB_II),
				ISMX(xv) + e2TSC(e2P_IS_II)),
		     p7_FLogsum(IDMX(xv) + e2TSC(e2P_ID_II),
				IIMX(xv) + e2TSC(e2P_II_II)));
     IIMX(x) = sc + Irsc;
   }
   /* DD state 12 transitions */
   DDMX(x) = -eslINFINITY;
#if 0
   sc = p7_FLogsum(p7_FLogsum(BBMX(x) + e2TSC(e2P_BB_DD) + ddfactor, 
			      IBMX(x) + e2TSC(e2P_IB_DD) + ddfactor),
		   p7_FLogsum(SSMX(x) + e2TSC(e2P_SS_DD) + ddfactor,
			      DSMX(x) + e2TSC(e2P_DS_DD) + ddfactor));
   sc = p7_FLogsum(sc,
		   p7_FLogsum(ISMX(x) + e2TSC(e2P_IS_DD) + ddfactor,
			      SDMX(x) + e2TSC(e2P_SD_DD) + ddfactor));
   sc = p7_FLogsum(p7_FLogsum(sc, 
			      IDMX(x) + e2TSC(e2P_ID_DD) + ddfactor),
		   p7_FLogsum(BIMX(x) + e2TSC(e2P_BI_DD) + ddfactor, 
			      SIMX(x) + e2TSC(e2P_SI_DD) + ddfactor));
   sc = p7_FLogsum(sc, 
		   p7_FLogsum(DIMX(x) + e2TSC(e2P_DI_DD) + ddfactor, 
			      IIMX(x) + e2TSC(e2P_II_DD) + ddfactor));
   DDMX(x) = sc;
#endif

   if (apair == NONE) {
     sc = p7_FLogsum(           IBMX(xv) + e2TSC(e2P_IB_DD),
		     p7_FLogsum(SSMX(xv) + e2TSC(e2P_SS_DD),
				DSMX(xv) + e2TSC(e2P_DS_DD)));
     sc = p7_FLogsum(sc,
		     p7_FLogsum(ISMX(xv) + e2TSC(e2P_IS_DD),
				SDMX(xv) + e2TSC(e2P_SD_DD)));
     sc = p7_FLogsum(p7_FLogsum(sc, 
				IDMX(xv) + e2TSC(e2P_ID_DD)),
		     p7_FLogsum(BIMX(xv) + e2TSC(e2P_BI_DD), 
				SIMX(xv) + e2TSC(e2P_SI_DD)));
     sc = p7_FLogsum(p7_FLogsum(sc,
				DDMX(xv) + e2TSC(e2P_DD_DD)), 
		     p7_FLogsum(DIMX(xv) + e2TSC(e2P_DI_DD), 
				IIMX(xv) + e2TSC(e2P_II_DD)));
     DDMX(x) = p7_FLogsum(DDMX(x),sc);
   }
  
    /* EE state 12 transitions */
    if (x == L) {
      sc = p7_FLogsum(           IBMX(x) + e2TSC(e2P_IB_EE),
		      p7_FLogsum(SSMX(x) + e2TSC(e2P_SS_EE),
				 DSMX(x) + e2TSC(e2P_DS_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 ISMX(x) + e2TSC(e2P_IS_EE)),
		      p7_FLogsum(SDMX(x) + e2TSC(e2P_SD_EE), 
				 DDMX(x) + e2TSC(e2P_DD_EE)));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 IDMX(x) + e2TSC(e2P_ID_EE)),
		      p7_FLogsum(BIMX(x) + e2TSC(e2P_BI_EE), 
				 SIMX(x) + e2TSC(e2P_SI_EE)));
      sc = p7_FLogsum(sc, 
		      p7_FLogsum(DIMX(x) + e2TSC(e2P_DI_EE), 
				 IIMX(x) + e2TSC(e2P_II_EE)));
      E2G_XMX(gx, x, e2G_EE) = sc;
    }
    else E2G_XMX(gx, x, e2G_EE) = -eslINFINITY;
    
#if 0
    if (x==L)printf("e2fFWD x %d apair %d iiss %d %f BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f  EE %f\n", 
		    x, apair, e2P_II_SS, e2TSC(e2P_II_SS),
		    BBMX(x), IBMX(x), SSMX(x), DSMX(x), ISMX(x), 
		    SDMX(x), DDMX(x), IDMX(x), 
		    BIMX(x), SIMX(x), DIMX(x), IIMX(x), E2G_XMX(gx, x, e2G_EE));
#endif
  }
  
  if (opt_sc != NULL) *opt_sc = E2G_XMX(gx, L, e2G_EE);
  if (*opt_sc <= -eslINFINITY) { printf("-infinity score\n"); return eslFAIL; }

  gx->Lrow = L;
  gx->Lcol = L;
  
  return eslOK;
}

/* Function:  e2_GBackward()
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
e2f_GBackward(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, E2_GMX *gx, float *opt_sc)
{
  float const  *tsc = (float const *)gm->tsc;
  float       **dp  = gx->dp;						
  float         sc;
  float         SSsc;
  float         Slsc, Srsc;
  float         Ilsc, Irsc;
  float        *pl = NULL;
  float        *pr = NULL;
  float         ddfactor;
  int           L;  
  int           K = psql->abc->K;
  int           isgapl, isgapr;
  enum apair_e  apair;
  int           b,c;                  /* alphabet index */
  int           x, xv;
  
  if (psql->n != psqr->n) { printf("sqs are not aligned\n"); return eslFAIL; }
  L = psql->n;

  /* Note: backward calculates the probability we can get *out* of
   * a cell(i,j); exclusive of emitting residue x_i or y_j.
   */
  p7_FLogsumInit();    
  
#if 0
  ddfactor = -log(1.0 - bg->p * exp(tsc[e2P_DD_DD]));
#endif
  
  /* Backward:
   *           states EE, DD need to be evaluated first.
   *
   * Order: EE, DD, BB, SS, DS, SD, IB, IS, ID, BI, SI, II
   */
  /* Initialize the row x = L  */
  x = L;
  
  /* EE state 0 transitions */
  E2G_XMX(gx, x, e2G_EE) = 0.0;
  /* DD state 6 - 5 = 1 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_DD_EE);
  DDMX(x) = sc;
  
 /* BB state 7 - 5 = 2 transitions */
  BBMX(x) = -eslINFINITY;
  /* SS state 7 - 5 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SS_EE);
  SSMX(x) = sc;
  /* DS state 7 - 5 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_DS_EE);
  DSMX(x) = sc;
  /* SD state 7 - 5 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SD_EE);
  SDMX(x) = sc;
  
  /* IB state 7 - 5 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_IB_EE);
  IBMX(x) = sc;      
  /* IS state 7 - 5 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_IS_EE);
  ISMX(x) = sc;
  /* ID state 7 - 5 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_ID_EE);
  IDMX(x) = sc;
  
  /* BI state 6 - 4 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_BI_EE);
  BIMX(x) = sc;
  /* SI state 6 - 4 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SI_EE);
  SIMX(x) = sc;     
  /* DI state 6 - 4 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_DI_EE);
  DIMX(x) = sc;      
  /* II state 6 - 4 = 2 transitions */
  sc = E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_II_EE);
  IIMX(x) = sc;

 #if 0
 /* BB state 7 - 5 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_BB_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_BB_EE));
  BBMX(x) = sc;
  /* SS state 7 - 5 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_SS_DD) + ddfactor, 
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SS_EE));
  SSMX(x) = sc;
  /* DS state 7 - 5 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_DS_DD) + ddfactor, 
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_DS_EE));
  DSMX(x) = sc;
  /* SD state 7 - 5 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_SD_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SD_EE));
  SDMX(x) = sc;
  
  /* IB state 7 - 5 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_IB_DD) + ddfactor, 
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_IB_EE));
  IBMX(x) = sc;      
  /* IS state 7 - 5 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_IS_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_IS_EE));
  ISMX(x) = sc;
  /* ID state 7 - 5 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_ID_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_ID_EE));
  IDMX(x) = sc;
  
  /* BI state 6 - 4 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_BI_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_BI_EE));
  BIMX(x) = sc;
  /* SI state 6 - 4 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_SI_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SI_EE));
  SIMX(x) = sc;     
  /* DI state 6 - 4 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_DI_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_DI_EE));
  DIMX(x) = sc;      
  /* II state 6 - 4 = 2 transitions */
  sc = p7_FLogsum(DDMX(x) + e2TSC(e2P_II_DD) + ddfactor,
		  E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_II_EE));
  IIMX(x) = sc;
#endif 

#if 0
  printf("e2fBCW x %d  BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f EE %f\n", x, 
	 BBMX(x), IBMX(x), SSMX(x), DSMX(x), ISMX(x), SDMX(x), 
	 DDMX(x), IDMX(x), BIMX(x), SIMX(x), DIMX(x), IIMX(x), E2G_XMX(gx, x, e2G_EE));
#endif    
  
  /* Main recursion */
  for (x = L-1; x >= 0; x--)
    {
      pl = psql->prof[x+1];
      pr = psqr->prof[x+1];
      SSsc = -eslINFINITY;
      Slsc = -eslINFINITY; /* orphan substitution score */
      Srsc = -eslINFINITY; /* orphan substitution score */
      Ilsc = -eslINFINITY; /* insertion score */
      Irsc = -eslINFINITY; /* insertion score */
      for (b = 0; b < K; b ++) {
	Slsc = p7_FLogsum(Slsc, gm->ssc[b][e2P_SL] + pl[b]);
	Ilsc = p7_FLogsum(Ilsc, gm->isc[b][e2P_SL] + pl[b]);
	Srsc = p7_FLogsum(Srsc, gm->ssc[b][e2P_SR] + pr[b]);
	Irsc = p7_FLogsum(Irsc, gm->isc[b][e2P_SR] + pr[b]);
	for (c = 0; c < K; c ++) {
	  SSsc = p7_FLogsum(SSsc, gm->sssc[b][c] + pl[b] + pr[c]);
	}
      }
      isgapl = isgap(pl, K);
      isgapr = isgap(pr, K);

      if      (!isgapl && !isgapr) apair = PAIR;
      else if ( isgapl && !isgapr) apair = RRES;
      else if (!isgapl &&  isgapr) apair = LRES;
      else                         apair = NONE;

      /* ahead one */
      xv = x + 1;
      
      /* order: EE, DD then all the rest in any order */
      
      /* EE state no transitions */
      E2G_XMX(gx, x, e2G_EE) = -eslINFINITY;
      /* DD state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_DD_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_DD_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? SDMX(xv) + e2TSC(e2P_DD_SD) + Slsc : -eslINFINITY,
				 (apair == LRES)? IDMX(xv) + e2TSC(e2P_DD_ID) + Ilsc : -eslINFINITY));
      sc = p7_FLogsum(p7_FLogsum(sc,
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_DD_DD)        : -eslINFINITY), 
		      p7_FLogsum((apair == RRES)? DIMX(xv) + e2TSC(e2P_DD_DI) + Irsc : -eslINFINITY,
				                  E2G_XMX(gx, x, e2G_EE)  + e2TSC(e2P_DD_EE)));
      DDMX(x) = sc;
      
      /* BB state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IBMX(xv) + e2TSC(e2P_BB_IB) + Ilsc : -eslINFINITY,
				 (apair == PAIR)? SSMX(xv) + e2TSC(e2P_BB_SS) + SSsc : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? DSMX(xv) + e2TSC(e2P_BB_DS) + Srsc : -eslINFINITY,
				 (apair == LRES)? SDMX(xv) + e2TSC(e2P_BB_SD) + Slsc : -eslINFINITY));
      sc = p7_FLogsum(sc,
		      (apair == RRES)? BIMX(xv) + e2TSC(e2P_BB_BI) + Irsc : -eslINFINITY);
      BBMX(x) = sc;
      /* SS state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_SS_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_SS_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? ISMX(xv) + e2TSC(e2P_SS_IS) + Ilsc : -eslINFINITY,
				 (apair == LRES)? SDMX(xv) + e2TSC(e2P_SS_SD) + Slsc : -eslINFINITY));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_SS_DD)        : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? SIMX(xv) + e2TSC(e2P_SS_SI) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SS_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_SS_DD) + ddfactor);
#endif
      SSMX(x) = sc;
      /* DS state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_DS_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_DS_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? ISMX(xv) + e2TSC(e2P_DS_IS) + Ilsc : -eslINFINITY,
				 (apair == LRES)? SDMX(xv) + e2TSC(e2P_DS_SD) + Slsc : -eslINFINITY));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_DS_DD)        : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? DIMX(xv) + e2TSC(e2P_DS_DI) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_DS_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_DS_DD) + ddfactor);
#endif
      DSMX(x) = sc;
      /* SD state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_SD_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_SD_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? SDMX(xv) + e2TSC(e2P_SD_SD) + Slsc : -eslINFINITY,
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_SD_DD)        : -eslINFINITY));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 (apair == LRES)? IDMX(xv) + e2TSC(e2P_SD_ID) + Ilsc : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? SIMX(xv) + e2TSC(e2P_SD_SI) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SD_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_SD_DD) + ddfactor);
#endif
      SDMX(x) = sc;
      
      /* IB state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == LRES)? IBMX(xv) + e2TSC(e2P_IB_IB) + Ilsc : -eslINFINITY, 
				 (apair == PAIR)? SSMX(xv) + e2TSC(e2P_IB_SS) + SSsc : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? DSMX(xv) + e2TSC(e2P_IB_DS) + Srsc : -eslINFINITY,
				 (apair == LRES)? SDMX(xv) + e2TSC(e2P_IB_SD) + Slsc : -eslINFINITY));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_IB_DD)        : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? IIMX(xv) + e2TSC(e2P_IB_II) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_IB_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_IB_DD) + ddfactor);
#endif
      IBMX(x) = sc;
      /* IS state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_IS_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_IS_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? ISMX(xv) + e2TSC(e2P_IS_IS) + Ilsc : -eslINFINITY,
				 (apair == LRES)? SDMX(xv) + e2TSC(e2P_IS_SD) + Slsc : -eslINFINITY));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_IS_DD)        : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? IIMX(xv) + e2TSC(e2P_IS_II) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_IS_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_IS_DD) + ddfactor);
#endif
      ISMX(x) = sc;
      /* ID state 7 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_ID_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_ID_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? SDMX(xv) + e2TSC(e2P_ID_SD) + Slsc : -eslINFINITY,
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_ID_DD)        : -eslINFINITY));
      sc = p7_FLogsum(p7_FLogsum(sc, 
				 (apair == LRES)? IDMX(xv) + e2TSC(e2P_ID_ID) + Ilsc : -eslINFINITY),
		      p7_FLogsum((apair == RRES)? IIMX(xv) + e2TSC(e2P_ID_II) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_ID_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_ID_DD) + ddfactor);
#endif
      IDMX(x) = sc;
      
      /* BI state 6 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_BI_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_BI_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? SDMX(xv) + e2TSC(e2P_BI_SD) + Slsc : -eslINFINITY,
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_BI_DD)        : -eslINFINITY));
      sc = p7_FLogsum(sc,
		      p7_FLogsum((apair == RRES)? BIMX(xv) + e2TSC(e2P_BI_BI) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_BI_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_BI_DD) + ddfactor);
#endif
      BIMX(x) = sc;
      /* SI state 6 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_SI_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_SI_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? SDMX(xv) + e2TSC(e2P_SI_SD) + Slsc : -eslINFINITY,
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_SI_DD)        : -eslINFINITY));
      sc = p7_FLogsum(sc, 
		      p7_FLogsum((apair == RRES)? SIMX(xv) + e2TSC(e2P_SI_SI) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_SI_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_SI_DD) + ddfactor);
#endif
      SIMX(x) = sc;
      /* DI state 6 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_DI_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_DI_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? SDMX(xv) + e2TSC(e2P_DI_SD) + Slsc : -eslINFINITY,
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_DI_DD)        : -eslINFINITY));
      sc = p7_FLogsum(sc, 
		      p7_FLogsum((apair == RRES)? DIMX(xv) + e2TSC(e2P_DI_DI) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_DI_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_DI_DD) + ddfactor);
#endif
      DIMX(x) = sc;
      /* II state 6 transitions */
      sc = p7_FLogsum(p7_FLogsum((apair == PAIR)? SSMX(xv) + e2TSC(e2P_II_SS) + SSsc : -eslINFINITY, 
				 (apair == RRES)? DSMX(xv) + e2TSC(e2P_II_DS) + Srsc : -eslINFINITY),
		      p7_FLogsum((apair == LRES)? SDMX(xv) + e2TSC(e2P_II_SD) + Slsc : -eslINFINITY,
				 (apair == NONE)? DDMX(xv) + e2TSC(e2P_II_DD)        : -eslINFINITY));
      sc = p7_FLogsum(sc, 
		      p7_FLogsum((apair == RRES)? IIMX(xv) + e2TSC(e2P_II_II) + Irsc : -eslINFINITY,
				 E2G_XMX(gx, x, e2G_EE) + e2TSC(e2P_II_EE)));
#if 0
      sc = p7_FLogsum(sc, DDMX(x) + e2TSC(e2P_II_DD) + ddfactor);
#endif
      IIMX(x) = sc;
#if 0
      if (x==0) printf("e2fBCW x %d  BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f EE %f\n", x, 
		       BBMX(x), IBMX(x), SSMX(x), DSMX(x), ISMX(x), SDMX(x), 
		       DDMX(x), IDMX(x), BIMX(x), SIMX(x), DIMX(x), IIMX(x), E2G_XMX(gx, x, e2G_EE));
#endif
    }
  
  if (opt_sc != NULL) *opt_sc = BBMX(0);
  
  gx->Lrow = L;
  gx->Lcol = L;
  
  return eslOK;
}


static int 
isgap(float *p, int K) 
{
  if (p[K] == 0.) return 1;
  else return 0;
}
