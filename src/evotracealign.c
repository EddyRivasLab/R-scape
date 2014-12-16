/* Construction of multiple alignments from traces.
 * 
 * Contents:
 *   1. API for aligning sequence or MSA traces
 *   2. Internal functions used by the API
 *   3. Example
 *   4. Copyright and license.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "base/p7_alidisplay.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "search/modelconfig.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_viterbi.h"
#include "dp_sparse/sparse_fwdback.h"
#include "dp_sparse/sparse_decoding.h"

#include "p7_evopipeline.h"
#include "evotracealign.h"

/*****************************************************************
 * 1. API for aligning sequence or MSA traces
 *****************************************************************/



/* Function: p7_tracealign_ComputeTraces()
 * Synopsis: Compute traces for an array of sequences.
 *
 * Purpose:  Given an <hmm> and a set of sequences <sq> (along with an
 *           <offset> into the first sequence for which a trace is
 *           desired), calculate the state path (trace) for each of
 *           <N> sequences. The calling function provides a allocated
 *           array of P7_TRACE's (<tr>) into which the results are
 *           placed.
 *
 * Returns:  <eslOK> on success.
 */
int
p7_evotracealign_ComputeTraces(P7_RATE *R, P7_HMM *hmm, ESL_SQ  **sq, int offset, int N, P7_TRACE  **tr, P7_BG *bg, int noevo, float fixtime, float tol, char *errbuf, int verbose)
{
  P7_PROFILE      *gm  = NULL;
  P7_OPROFILE     *om  = NULL;
  P7_CHECKPTMX    *cx  = NULL;
  P7_SPARSEMASK   *sm  = NULL;
  P7_SPARSEMX     *sx1 = NULL;
  P7_SPARSEMX     *sx2 = NULL;
  float            time;
  float            firststep = 1e+0;
  int              idx;
  float            fwdsc;  /* Forward raw score, nats         */
  float            vsc;	   /* Viterbi raw score, nats         */
  int              unihit = TRUE;

  gm = p7_profile_Create (hmm->M, hmm->abc);
  om = p7_oprofile_Create(hmm->M, hmm->abc);

  /* May need to play with config choice 
   * We want the <om> in local multihit mode for sparsification step
   * But we want the <gm> in unihit dual local/glocal mode for alignment
   */
  p7_profile_Config(gm, hmm, bg);                    /* *multihit* glocal/local...   */
  p7_oprofile_Convert(gm, om);	                     /*    ... *multihit* local      */
  p7_profile_ConfigCustom(gm, hmm, bg, 0, 0.0, 0.5); /*    ... *unihit* glocal/local */

  cx  = p7_checkptmx_Create (hmm->M, sq[offset]->n, ESL_MBYTES(p7_RAMLIMIT));
  sm  = p7_sparsemask_Create(hmm->M, sq[offset]->n);
  sx1 = p7_sparsemx_Create  (sm);
  sx2 = p7_sparsemx_Create  (sm);

  for (idx = offset; idx < offset+ N; idx++)
    {
      /* special case: a sequence of length 0. HMMER model can't generate 0 length seq. Set tr->N == 0 as a flag. (bug #h100 fix) */
      if (sq[idx]->n == 0) { tr[idx]->N = 0; continue; }

      time = 1.0;

      p7_oprofile_ReconfigLength(om, sq[idx]->n);
      p7_profile_SetLength      (gm, sq[idx]->n);

      /* Sparse mask is constructed by local multihit alignment to <om> in the filters */
      //p7_ForwardFilter (sq[idx]->dsq, sq[idx]->n, om, cx, &fwdsc);
      p7_evopli_OptimizeForwardFilter(sq[idx]->dsq, sq[idx]->n, &time, R, hmm, gm, om, bg, cx, &fwdsc, firststep, noevo, fixtime, unihit, tol, errbuf, verbose);
      //p7_BackwardFilter(sq[idx]->dsq, sq[idx]->n, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT);

      /* Alignment itself is constructed by unihit dual-mode local/glocal to <gm> */
      //p7_SparseViterbi (sq[idx]->dsq, sq[idx]->n, gm, sm, sx1, tr[idx], &vsc);
      p7_evopli_OptimizeSparseViterbi  (sq[idx]->dsq, sq[idx]->n, &time, R, hmm, gm, om, bg, sm, cx, sx1, tr[idx],   &vsc, firststep, TRUE, fixtime, unihit, tol, errbuf, verbose);
      if (tr[idx]->pp) {
	p7_evopli_OptimizeSparseForward(sq[idx]->dsq, sq[idx]->n, &time, R, hmm, gm, om, bg, sm, cx, sx1,          &fwdsc, firststep, TRUE, fixtime, unihit, tol, errbuf, verbose);
	// p7_SparseForward (sq[idx]->dsq, sq[idx]->n, gm, sm, sx1, &fwdsc);
	p7_SparseBackward(sq[idx]->dsq, sq[idx]->n, gm, sm, sx2, &fwdsc);
	p7_SparseDecoding(sq[idx]->dsq, sq[idx]->n, gm, sx1, sx2, sx2); /* decode sx2 in place */
	p7_sparsemx_TracePostprobs(sx2, tr[idx]);
      }
      p7_trace_Index(tr[idx]);
      
      p7_checkptmx_Reuse(cx);
      p7_sparsemask_Reuse(sm);
      p7_sparsemx_Reuse(sx1);
      p7_sparsemx_Reuse(sx2);
    }

  p7_sparsemx_Destroy(sx1);
  p7_sparsemx_Destroy(sx2);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  return eslOK;
}



/*--------------- end, exposed API ------------------------------*/




/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id: tracealign.c 4514 2013-07-13 14:00:23Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk/src/misc/tracealign.c $
 *****************************************************************/



