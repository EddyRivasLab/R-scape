/* Posterior decoding algorithms; generic versions.
 * 
 * Contents:
 *   1. Posterior decoding algorithms.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 */
#include "p7_config.h"

#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"
#include "e2.h"
#include "e2f_generic_decoding.h"

/*****************************************************************
 * 1. Posterior decoding algorithms.
 *****************************************************************/
/* Function:  e2f_GDecoding()
 * Synopsis:  Posterior decoding of residue assignments.
 *
 * Purpose:   Calculates a posterior decoding of the residues in a
 *            target sequence, given profile <gm> and filled Forward
 *            and Backward matrices <fwd>, <bck> for the profile
 *            aligned to that target sequence. The resulting posterior
 *            decoding is stored in a DP matrix <pp>, provided by the
 *            caller. 
 *            
 *            There are 11 emitting states:
 *
 *            SS: emits residue <i> and residue <j>
 *            DS: emits                 residue <j>
 *            SD: emits residue <i>
 *
 *            IB: emits residue <i>
 *            IS: emits residue <i>
 *            ID: emits residue <i>
 *            IE: emits residue <i>
 *
 *            BI: emits                 residue <j>
 *            SI: emits                 residue <j>
 *            DI: emits                 residue <j>
 *            II: emits                 residue <j>
 *
 *            where:
 *            <i> index for longer sequence, 
 *            <j> index for shorter sequence.
 *
 *            The sum over all these possibilities for a given 
 *            residue <i> (5 terms) is 1.0. The sum over all 
 *            these possibilities for a given residue <j> 
 *            (7 terms) is 1.0.
 *
 *            The caller may pass the Backward matrix <bck> as <pp>,
 *            in which case <bck> will be overwritten with
 *            <pp>. However, the caller may \emph{not} overwrite <fwd>
 *            this way; an <(i-1)> dependency in the calculation of
 *            NN, CC, JJ transitions prevents this.
 *
 * Args:      gm   - profile (must be the same that was used to fill <fwd>, <bck>).
 *            fwd  - filled Forward matrix 
 *            bck  - filled Backward matrix
 *            pp   - RESULT: posterior decoding matrix.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Note:      Burns time renormalizing each row. If you don't do this,
 *            probabilities will have an error of +/- 0.001 or so, creeping
 *            in from error in FLogsum()'s table approximation and even
 *            in log() and exp() themselves; including "probabilities"
 *            up to  ~1.001. Though this isn't going to break anything
 *            in normal use, it does drive the unit tests wild; the SSE
 *            implementation is more accurate, and unit tests that try
 *            to compare SSE and generic results will see differences,
 *            some sufficient to alter the choice of OA traceback.
 *    
 */
int
e2f_GDecoding(const E2_PROFILE *gm, const E2_GMX *fwd, E2_GMX *bck, E2_GMX *pp, char *errbuf)
{
  float      **dp = pp->dp;
  int          L  = fwd->Lrow;
  float        overall_sc = E2G_XMX(fwd, L, e2G_EE);
  float        sum;
  float        denom;
  int          x;                    /* linear memory index */
  int          status;

  if (overall_sc <= -eslINFINITY) ESL_XFAIL(eslFAIL, errbuf, "-INFINITY score");

  for (x = 0; x <= L; x++) {
    
    BBMX(x) = 0.0;
    SSMX(x) = expf(fwd->dp[x][e2G_SS] + bck->dp[x][e2G_SS] - overall_sc); 
    DSMX(x) = expf(fwd->dp[x][e2G_DS] + bck->dp[x][e2G_DS] - overall_sc); 
    SDMX(x) = expf(fwd->dp[x][e2G_SD] + bck->dp[x][e2G_SD] - overall_sc); 
    DDMX(x) = expf(fwd->dp[x][e2G_DD] + bck->dp[x][e2G_DD] - overall_sc); 
    
    IBMX(x) = expf(fwd->dp[x][e2G_IB] + bck->dp[x][e2G_IB] - overall_sc);
    ISMX(x) = expf(fwd->dp[x][e2G_IS] + bck->dp[x][e2G_IS] - overall_sc); 
    IDMX(x) = expf(fwd->dp[x][e2G_ID] + bck->dp[x][e2G_ID] - overall_sc); 
    
    BIMX(x) = expf(fwd->dp[x][e2G_BI] + bck->dp[x][e2G_BI] - overall_sc); 
    SIMX(x) = expf(fwd->dp[x][e2G_SI] + bck->dp[x][e2G_SI] - overall_sc); 
    DIMX(x) = expf(fwd->dp[x][e2G_DI] + bck->dp[x][e2G_DI] - overall_sc); 
    IIMX(x) = expf(fwd->dp[x][e2G_II] + bck->dp[x][e2G_II] - overall_sc); 
    
  }
  
  /* renormalize */
  for (x = 0; x <= L; x++) {
    sum = 0.0;
    sum += SSMX(x); 
    sum += DSMX(x);
    sum += SDMX(x);
    sum += DDMX(x);
    sum += IBMX(x);
    sum += ISMX(x); 
    sum += IDMX(x);
    sum += BIMX(x);
    sum += SIMX(x); 
    sum += DIMX(x);
    sum += IIMX(x);
  
    denom = (sum > 0.)? 1.0 / sum : 1.0;
    SSMX(x) *= denom; 
    DSMX(x) *= denom;
    SDMX(x) *= denom;
    DDMX(x) *= denom;
    IBMX(x) *= denom;
    ISMX(x) *= denom; 
    IDMX(x) *= denom;
    BIMX(x) *= denom;
    SIMX(x) *= denom; 
    DIMX(x) *= denom;
    IIMX(x) *= denom;

#if 0
      if (x==0) printf("decoding x %d  BB %f IB %f SS %f DS %F IS %f SD %f DD %f ID %f BI %f SI %f DI %f II %f | %f %f %f\n", x, 
		       BBMX(x), IBMX(x), SSMX(x), DSMX(x), ISMX(x), SDMX(x), 
		       DDMX(x), IDMX(x), BIMX(x), SIMX(x), DIMX(x), IIMX(x), fwd->dp[x][e2G_SS], bck->dp[x][e2G_SS], overall_sc);
#endif

  }
  
  return eslOK;

 ERROR:
  return status;
}
