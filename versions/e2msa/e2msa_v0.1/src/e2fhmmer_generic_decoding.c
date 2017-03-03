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
#include "e2fhmmer_generic_decoding.h"

/*****************************************************************
 * 1. Posterior decoding algorithms.
 *****************************************************************/
/* Function:  e2fhmmer_GDecoding()
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
e2fhmmer_GDecoding(const E2HMMER_PROFILE *gm, const E2_GMX *fwd, E2_GMX *bck, float overall_sc, E2_GMX *pp)
{
  float      **dp = pp->dp;
  int          L  = fwd->Lrow;
  int          M  = fwd->M;
  float        sum;
  float        denom;
  int          x;                    /* linear memory index */
  int          k, kk;

  for (x = 0; x <= L; x++) {
    for (k = 0; k <= M; k ++) {
      kk = k*e2G_NSCELLS;
      BBMXM(x,k) = expf(fwd->dp[x][kk+e2G_BB] + bck->dp[x][kk+e2G_BB] - overall_sc); 
      SSMXM(x,k) = expf(fwd->dp[x][kk+e2G_SS] + bck->dp[x][kk+e2G_SS] - overall_sc); 
      DSMXM(x,k) = expf(fwd->dp[x][kk+e2G_DS] + bck->dp[x][kk+e2G_DS] - overall_sc); 
      SDMXM(x,k) = expf(fwd->dp[x][kk+e2G_SD] + bck->dp[x][kk+e2G_SD] - overall_sc); 
      DDMXM(x,k) = expf(fwd->dp[x][kk+e2G_DD] + bck->dp[x][kk+e2G_DD] - overall_sc); 
      
      IBMXM(x,k) = expf(fwd->dp[x][kk+e2G_IB] + bck->dp[x][kk+e2G_IB] - overall_sc);
      ISMXM(x,k) = expf(fwd->dp[x][kk+e2G_IS] + bck->dp[x][kk+e2G_IS] - overall_sc); 
      IDMXM(x,k) = expf(fwd->dp[x][kk+e2G_ID] + bck->dp[x][kk+e2G_ID] - overall_sc); 
      
      BIMXM(x,k) = expf(fwd->dp[x][kk+e2G_BI] + bck->dp[x][kk+e2G_BI] - overall_sc); 
      SIMXM(x,k) = expf(fwd->dp[x][kk+e2G_SI] + bck->dp[x][kk+e2G_SI] - overall_sc); 
      DIMXM(x,k) = expf(fwd->dp[x][kk+e2G_DI] + bck->dp[x][kk+e2G_DI] - overall_sc); 
      IiMXM(x,k) = expf(fwd->dp[x][kk+e2G_II] + bck->dp[x][kk+e2G_II] - overall_sc); 
      iIMXM(x,k) = expf(fwd->dp[x][kk+e2G_ii] + bck->dp[x][kk+e2G_ii] - overall_sc); 
      
    }
     
    /* renormalize */
    sum = 0.0;
    for (k = 0; k <= M; k++) {
      sum += SSMXM(x,k); 
      sum += DSMXM(x,k);
      sum += SDMXM(x,k);
      sum += DDMXM(x,k);
      sum += IBMXM(x,k);
      sum += ISMXM(x,k); 
      sum += IDMXM(x,k);
      sum += BIMXM(x,k);
      sum += SIMXM(x,k); 
      sum += DIMXM(x,k);
      sum += IiMXM(x,k);
      sum += iIMXM(x,k);
    }
    
    denom = (sum > 0.)? 1.0 / sum : 1.0;
    
    for (k = 0; k <= M; k++) {
      SSMXM(x,k) *= denom; 
      DSMXM(x,k) *= denom;
      SDMXM(x,k) *= denom;
      DDMXM(x,k) *= denom;
      IBMXM(x,k) *= denom;
      ISMXM(x,k) *= denom; 
      IDMXM(x,k) *= denom;
      BIMXM(x,k) *= denom;
      SIMXM(x,k) *= denom; 
      DIMXM(x,k) *= denom;
      IiMXM(x,k) *= denom;
      iIMXM(x,k) *= denom;
    }

#if 0
    for (k = 0; k <= M; k ++) {
      printf("e2fhmmer_Decoding x %d  k %d BB %f IB %f SS %f DS %f IS %f SD %f DD %f ID %f BI %f SI %f DI %f Ii %f  iI %f  \n", 
	     x, k, 
	     BBMXM(x,k), IBMXM(x,k), SSMXM(x,k), DSMXM(x,k), ISMXM(x,k), 
	     SDMXM(x,k), DDMXM(x,k), IDMXM(x,k), 
	     BIMXM(x,k), SIMXM(x,k), DIMXM(x,k), IiMXM(x,k), iIMXM(x,k));
    }
#endif
  }
  
  
  return eslOK;
}
