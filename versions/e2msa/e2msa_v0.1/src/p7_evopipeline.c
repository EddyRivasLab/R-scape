/* H3's accelerated seq/profile comparison pipeline
 * Contents:
 *   1. P7_PIPELINE: allocation, initiazation, destruction
 *   2. Setting pipeline for next target or model
 *   3. Testing bitscore/E-value against reporting/inclusion thresholds
 *   4. Statistics output from a completed pipeline.
 *   5. The main pipeline call: p7_Pipeline().
 *   6. The pipeline specialized for nhmmer/longtarget.
 *   7. The acceleration pipeline, vastly simplified, as one call. 
 *   8. Example 1: search mode (in a sequence db)
 *   9. Example 2: scan mode (in an HMM db)
 *   10. Copyright and license information
 * 
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "easel.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_domain.h"
#include "base/p7_profile.h"   /* used by the hmmscan workaround */
#include "base/p7_scoredata.h"
#include "base/p7_tophits.h"
#include "base/p7_hmmwindow.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/io.h"
#include "dp_vector/msvfilter.h"
#include "dp_vector/vitfilter.h"
#include "dp_vector/fwdfilter.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_viterbi.h"
#include "dp_sparse/sparse_fwdback.h"
#include "dp_sparse/sparse_decoding.h"
#include "dp_sparse/sparse_masstrace.h"
#include "dp_sparse/sparse_envscore.h"
#include "dp_sparse/sparse_null2.h"


#include "misc/logsum.h"

#include "search/modelconfig.h" /* used by the hmmscan workaround */
#include "search/p7_pipeline.h"

#include "e2_config.h"
#include "e1_rate.h"
#include "evohmmer.h"
#include "minimize.h"
#include "p7_evopipeline.h"
#include "ratematrix.h"

#define RECALIBRATE 0

static int    optimize_msvfilter(const ESL_DSQ *dsq, int n, float *ret_time,
				 P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx, 
				 float *ret_usc, double firststep, float fixtime, int noevo, float tol, char *errbuf, int be_verbose);
static int    optimize_pack_paramvector        (double *p, long np, struct optimize_data *data);
static int    optimize_unpack_paramvector      (double *p, long np, struct optimize_data *data);
static void   optimize_bracket_define_direction(double *p, long np, struct optimize_data *data);
static double optimize_msvfilter_func          (double *p, long np, void *dptr);
static double optimize_viterbifilter_func      (double *p, long np, void *dptr);
static double optimize_forwardfilter_func      (double *p, long np, void *dptr);
static double optimize_sparseviterbi_func      (double *p, long np, void *dptr);
static double optimize_sparseforward_func      (double *p, long np, void *dptr);
static double func_msvfilter(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx,  
			     float time, int noevo, float tol, char *errbuf, int verbose);
static double func_viterbifilter(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx,  
				 float time, int unihit, int noevo, float tol, char *errbuf, int verbose);
static double func_forwardfilter(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_CHECKPTMX *cx, 
				 float time, int unihit, int noevo, float tol, char *errbuf, int verbose);
static double func_sparseviterbi(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxx, P7_TRACE *tr, 
				 float time, int unihit, int noevo, float tol, char *errbuf, int verbose);
static double func_sparseforward(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxf,
				 float time, int unihit, int noevo, float tol, char *errbuf, int verbose);

/* SRE: FIXME 3.1 in progress */
static int workaround_get_profile(P7_PIPELINE *pli, const P7_OPROFILE *om, const P7_BG *bg, P7_PROFILE **ret_gm);

static int workaround_get_starprofile(P7_PIPELINE *pli, const P7_BG *bg, const P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_HMM **ret_ehmm, int noevo);
static int workaround_evolve_profile(double time, int n, const P7_RATE *R, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, int unihit, int verbose);
static int workaround_calibrate(ESL_RANDOMNESS *r, int n, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om);
static int workaroud_scanmodels(P7_PIPELINE *pli, int n, EMRATE *emR, P7_RATE **ret_R, P7_HMM **ret_ehmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, int noevo, float fixtime, float tol, char *errbuf);

/*****************************************************************
 * 5. The main pipeline call: p7_Pipeline()
 *****************************************************************/


/* Function:  p7_Pipeline()
 * Synopsis:  HMMER3's accelerated seq/profile comparison pipeline.
 *
 * Purpose:   Run H3's accelerated pipeline to compare profile <om>
 *            against sequence <sq>. If a significant hit is found,
 *            information about it is added to the <hitlist>. The pipeline 
 *            accumulates beancounting information about how many comparisons
 *            flow through the pipeline while it's active.
 *            
 *            Each time <p7_Pipeline()> is called on a new sequence,
 *            caller must first call <p7_pipeline_NewSeq()> as an
 *            initialization step. Each time it is called on a new
 *            model, caller must first call <p7_pipeline_NewModel()>.
 *
 *            Caller must first configure the length parameterization
 *            of <om> and <bg> as it wants; that will not be done
 *            here. Probably caller wants to do
 *            <p7_oprofile_ReconfigLength(om, sq->n)> and
 *            <p7_bg_SetLength(bg, sq->n)>.
 *            
 *            Caller does not need to worry about the allocation sizes
 *            of the DP matrices inside the <pli> object. DP routines
 *            automatically resize these matrices as needed for the
 *            <gm->M> by <sq->n> comparison. As a result, note that
 *            the internals of these matrices in <pli> may have
 *            changed upon return, because of reallocation; caller
 *            cannot expect any pointers to remain valid.
 *            
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>. 
 *            Statistics in <pli> are updated; DP matrices in <pli>
 *            may have been reallocated. We currently do not guarantee
 *            anything about the state of the data in those matrices;
 *            it is considered to be internal data for <pli>.
 *            
 *            <eslEINVAL> if (in a scan pipeline) we're supposed to
 *            set GA/TC/NC bit score thresholds but the model doesn't
 *            have any.
 *            
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 * 
 * Future:    Right now it needs both <gm> and <om> because the filters
 *            use <om> and sparse DP uses <gm>. Not a problem in hmmsearch,
 *            which typically has both; more of a problem for hmmscan,
 *            which reads <om>. In the future, we need to optimize a bit
 *            more; perhaps we'll have P7_MODEL, holding annotation,
 *            hmm, profile params, and striped vector params for MSV,
 *            Vit, and Fwd/Back, with the ability to create vector params
 *            from profile and vice versa, on demand.
 */
int
p7_EvoPipeline(P7_PIPELINE *pli, EMRATE *emR, P7_RATE *oR, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *hitlist, int noevo, float fixtime, int unihit)
{
  ESL_RANDOMNESS  *r       = NULL;	  /* RNG for E-value calibration simulations */
  P7_HMM          *ehmm    = NULL;
  P7_RATE         *R       = (oR)? oR : NULL;
  P7_DOMAIN       *dcl     = NULL;        /* array of domain data structures <0..ndom-1> */
  P7_HIT          *hit     = NULL;        /* ptr to the current hit output data          */
  float            usc, vitsc, fwdsc;     /* DP scores                                   */
  float            spvitsc, spfwdsc;      /* sparse DP scores                            */
  float            filtersc;              /* HMM null filter score                       */
  float            nullsc;                /* null model score                            */
  float            seqbias;  
  float            seq_score;             /* the corrected per-seq bit score             */
  float            sum_score;             /* the corrected reconstruction score          */
  float            pre_score, pre2_score; /* uncorrected bit scores for seq              */
  double           P;                     /* P-value of a hit                            */
  double           lnP;                   /* log P-value of a hit                        */
  int              Ld;                    /* # of residues in envelopes                  */
  int              d,z,i;
  float            null2[p7_MAXCODE];
  int              noverlaps;
  int              last_ibe;
  int              best_d;
  float            time;                 /* fwdfilter      time */
  float            spvtime;              /* sparse viterbi time */
  float            spftime;              /* sparse forward time */
  float            tol = 0.1;
  float            firststep;
  int              be_verbose = FALSE;
  int              status;
  
 /* init */
  time = 1.0; 
  firststep = 1e+0;
  
#if RECALIBRATE
  r = esl_randomness_CreateFast(seed);
#endif
  
  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */

  /* Base null model score (we could calculate this in NewSeq(), for a scan pipeline) */
  p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

  /* Make a copy of hmm into ehmm which will be evolved */
  workaround_get_starprofile(pli, bg, hmm, gm, om, &ehmm, noevo); 

  /* First level filter: the SSV and MSV filters */
  // noevo = TRUE
  if ((status = optimize_msvfilter(sq->dsq, sq->n, &time, R, ehmm, gm, om, bg, pli->fx, 
				   &usc, firststep, fixtime, TRUE, tol, pli->errbuf, be_verbose)) != eslOK)      
    printf("\nsequence %s msvfilter did not optimize\n", sq->name);  
  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  if (P > pli->F1) goto NOHIT;
  pli->stats.n_past_msv++;

  /* biased composition HMM filtering */
  if (pli->do_biasfilter)
    {
      p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);
      seq_score = (usc - filtersc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > pli->F1) goto NOHIT;
    }
  else filtersc = nullsc;
  pli->stats.n_past_bias++;

  /* In scan mode, if it passes the MSV filter, read the rest of the profile */
  if (pli->mode == p7_SCAN_MODELS)
    {
      workaroud_scanmodels(pli, sq->n, emR, &R, &ehmm, gm, om, bg, noevo, fixtime, 0.0001, pli->errbuf);
    }

  /* Second level filter: ViterbiFilter(), multihit with <om> */
  if (P > pli->F2)
    {
      if ((status = p7_evopli_OptimizeViterbiFilter(sq->dsq, sq->n, &time, R, ehmm, gm, om, bg, pli->fx, 
						    &vitsc, firststep, noevo, fixtime, unihit, tol, pli->errbuf, be_verbose)) != eslOK) 
	printf("\nsequence %s vitfilter did not optimize\n", sq->name);
      
      seq_score = (vitsc-filtersc) / eslCONST_LOG2;
      P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
      if (P > pli->F2) goto NOHIT;
    }
  pli->stats.n_past_vit++;

  /* Checkpointed Forward. Check score as a filter step before proceeding to backwards/decoding. */
  time = 1.0;
  if ((status = p7_evopli_OptimizeForwardFilter(sq->dsq, sq->n, &time, R, ehmm, gm, om, bg, pli->cx, 
						&fwdsc, firststep, noevo, fixtime, unihit, tol, pli->errbuf, be_verbose)) != eslOK)      
    printf("\nsequence %s forwardfilter did not optimize\n", sq->name);

#if RECALIBRATE
  workaround_calibrate(r, sq->n, bg, ehmm, gm, om);
#endif

  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > pli->F3) goto NOHIT;
  pli->stats.n_past_fwd++;

#if 1
  printf("\nrealhit: fwdsc %f seq_score %f time %f P %f F3 %f| SEQ: %s\n", fwdsc, seq_score, time, P, pli->F3, sq->name);
#endif

  /* ok, it's for real; passes vectorized local scoring filters.
   * Finish with Backwards and Decoding, defining a sparse mask in <sm>.
   */
  //p7_BackwardFilter(sq->dsq, sq->n, om, pli->cx, pli->sm, p7_SPARSIFY_THRESH);

#if 0
  /* FIXME 3.1
   * hmmscan needs <gm>; read from the hmm file (.h3m) on disk
   * expect that this is slow and inefficient; come back and refactor later.
   * (in hmmpgmd, hfp is NULL; all our profiles are already in a cache)
   */
  if (pli->mode == p7_SCAN_MODELS && pli->hfp)
    workaround_get_profile(pli, om, bg, &gm);
#endif

  /* Now we can hand it over to sparse DP, with the full glocal/local model */
  //p7_SparseViterbi (sq->dsq, sq->n, gm, pli->sm,  pli->sxx, pli->tr, &spvitsc);
  spvtime = spftime = time;
  p7_evopli_OptimizeSparseViterbi(sq->dsq, sq->n, &spvtime, R, ehmm, gm, om, bg, pli->sm, pli->cx, pli->sxx, pli->tr, &spvitsc, firststep, TRUE, fixtime, unihit, tol, pli->errbuf, be_verbose);
  
  //p7_SparseForward (sq->dsq, sq->n, gm, pli->sm,  pli->sxf,          &spfwdsc);
  p7_evopli_OptimizeSparseForward(sq->dsq, sq->n, &spftime, R, ehmm, gm, om, bg, pli->sm, pli->cx, pli->sxf,          &spfwdsc, firststep, TRUE, fixtime, unihit, tol, pli->errbuf, be_verbose);

  p7_SparseBackward(sq->dsq, sq->n, gm, pli->sm,  pli->sxb,          NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, pli->sxf, pli->sxb, pli->sxd);
  p7_sparsemx_TracePostprobs(pli->sxd, pli->tr); /* annotate the trace. */
  p7_trace_Index(pli->tr);                       /* index domains in the trace */
  p7_sparsemx_Reuse(pli->sxx);                   /* reuse the Viterbi matrix for mass tracing, coming up next */

#if 1
  printf("        sparse_vitsc %f sparse_fwdsc %f ndom %d time_v %f time_f %f| SEQ: %s\n", spvitsc, spfwdsc, pli->tr->ndom, spvtime, spftime, sq->name);
#endif

  dcl = p7_domain_Create(pli->tr->ndom);
  //ESL_ALLOC(dcl, sizeof(P7_DOMAIN) * pli->tr->ndom);
  ESL_REALLOC(pli->wrk,  sizeof(float) * (gm->M+1));
  ESL_REALLOC(pli->n2sc, sizeof(float) * (sq->n+1));
  esl_vec_FSet(pli->n2sc, sq->n+1, 0.0f);

  noverlaps  = 0;
  last_ibe   = -1;
  best_d     = 0;
  for (d = 0; d < pli->tr->ndom; d++)
    {
      /* Determine envelope coords by mass trace. */
      p7_SparseMasstrace(sq->dsq, sq->n, gm, pli->sxf, pli->sxb, pli->tr, pli->tr->anch[d], p7_SPARSIFY_THRESH, pli->sxx, pli->mt, 
                         &(dcl[d].iae), &(dcl[d].ibe), &(dcl[d].kae), &(dcl[d].kbe));
      p7_sparsemx_Reuse (pli->sxx);
      p7_masstrace_Reuse(pli->mt);

      /* Keep track of overlaps */
      if (dcl[d].iae <= last_ibe) noverlaps++;
      last_ibe = dcl[d].ibe;

      /* Transfer alignment coords from trace. [Do we need to do this?] */
      dcl[d].ia = pli->tr->sqfrom[d];  
      dcl[d].ib = pli->tr->sqto[d];  
      dcl[d].ka = pli->tr->hmmfrom[d];  
      dcl[d].kb = pli->tr->hmmto[d];  

      /* Determine envelope score. [We have a fast approximation available, but empirical experiments determined that it rarely holds well.] */
      p7_SparseEnvscore(sq->dsq, sq->n, gm, dcl[d].iae, dcl[d].ibe, dcl[d].kae, dcl[d].kbe, pli->sm, pli->sxx, &(dcl[d].envsc));
      p7_sparsemx_Reuse(pli->sxx);

      /* Determine null2 correction 
       * Upon return, null2[x] = \log f'(x)/f(x).
       * Then we record that in n2sc[iae..ibe]: seqbias correction mustn't overcount overlapped envelopes
       */
      if (pli->do_null2) 
      {
        p7_sparse_Null2ByExpectation(gm, pli->sxd, dcl[d].iae, dcl[d].ibe, dcl[d].kae, dcl[d].kbe, pli->wrk, null2);
        dcl[d].domcorrection = 0.;
        for (i = dcl[d].iae; i <= dcl[d].ibe; i++)
          {
            dcl[d].domcorrection += null2[ sq->dsq[i] ];
            pli->n2sc[i]          = null2[ sq->dsq[i] ];
          }
        dcl[d].dombias = p7_FLogsum(0., log(bg->omega) + dcl[d].domcorrection);
      }
      else 
      {
        dcl[d].domcorrection = 0.;
        dcl[d].dombias       = 0.;
      }

      /* alignment "accuracy score": expected number of correctly aligned positions */
      dcl[d].oasc = 0.;
      for (z = pli->tr->tfrom[d]; z <= pli->tr->tto[d]; z++)
        if (pli->tr->i[z]) dcl[d].oasc += pli->tr->pp[z];

      /* domain bit score = null2-corrected bit score of the sequence if this were the only envelope in it. */
      dcl[d].bitscore = (dcl[d].envsc - (nullsc + dcl[d].dombias)) / eslCONST_LOG2;
      dcl[d].lnP      = esl_exp_logsurv( dcl[d].bitscore, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

      if (dcl[d].bitscore > dcl[best_d].bitscore) best_d = d;

      /* Viterbi alignment of the domain */
      dcl[d].ad = p7_alidisplay_Create(pli->tr, d, om, sq);

      /* We're initializing a P7_DOMAIN structure in dcl[d] by hand, without a Create().
       * We're responsible for initiazing all elements of this structure.
       */
      dcl[d].is_reported = FALSE; /* will get set later by p7_tophits_Threshold() */
      dcl[d].is_included = FALSE; /* ditto */
    }
  
  /* Calculate the null2-corrected per-seq score */
  seqbias = (pli->do_null2 ? p7_FLogsum(0.0, log(bg->omega) + esl_vec_FSum(pli->n2sc, sq->n+1)) : 0.0);
  pre_score =  (spfwdsc - nullsc) / eslCONST_LOG2;              /* BITS */
  seq_score =  (spfwdsc - (nullsc + seqbias)) / eslCONST_LOG2;  /* BITS */

  
  /* Calculate the "reconstruction score": estimated
   * per-sequence score as sum of individual domains,
   * discounting domains that don't contribute positive
   * score after being null2-corrected.
   */
  if (pli->do_null2) 
    {
      float baseline = sq->n * gm->xsc[p7P_C][p7P_LOOP] + gm->xsc[p7P_C][p7P_MOVE];
      sum_score = 0.0f;
      seqbias   = 0.0f;
      Ld        = 0;
      for (d = 0; d < pli->tr->ndom; d++)
      {
        if (dcl[d].envsc - dcl[d].dombias > baseline)
          {
            sum_score += dcl[d].envsc - (sq->n - (dcl[d].ibe - dcl[d].iae + 1)) * gm->xsc[p7P_C][p7P_LOOP] - gm->xsc[p7P_C][p7P_MOVE];
            Ld        += dcl[d].ibe - dcl[d].iae + 1;
            seqbias   += dcl[d].domcorrection; /* NATS */
          }
      }
      seqbias    = p7_FLogsum(0.0, log(bg->omega) + seqbias);      /* NATS */
      sum_score += (sq->n - Ld) * gm->xsc[p7P_C][p7P_LOOP] + gm->xsc[p7P_C][p7P_MOVE];
      pre2_score = (sum_score - nullsc) / eslCONST_LOG2;           /* BITS */
      sum_score  = (sum_score - (nullsc+seqbias)) / eslCONST_LOG2; /* BITS */
    }
  else 
    {
      pre2_score = pre_score;
      sum_score  = seq_score;
    }

  /* SRE TESTING: always use the reconstruction score. The fwd score is not adequately corrected */
  seq_score = sum_score;
  pre_score = pre2_score;


#if 0
  /* Let sum_score override the seq_score when it's better, and it includes at least 1 domain */
  if (Ld > 0 && sum_score > seq_score)
    {
      seq_score = sum_score;
      pre_score = pre2_score;
    }
#endif

  /* Apply thresholding and determine whether to put this
   * target into the hit list. E-value thresholding may
   * only be a lower bound for now, so this list may be longer
   * than eventually reported.
   */
  lnP =  esl_exp_logsurv (seq_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  if (! p7_pipeline_TargetReportable(pli, seq_score, lnP))
    {
      p7_domain_Destroy(dcl, pli->tr->ndom);
    }
  else
    {
      p7_tophits_CreateNextHit(hitlist, &hit);
      if (pli->mode == p7_SEARCH_SEQS) {
        if (                       (status  = esl_strdup(sq->name, -1, &(hit->name)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
        if (sq->acc[0]  != '\0' && (status  = esl_strdup(sq->acc,  -1, &(hit->acc)))   != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
        if (sq->desc[0] != '\0' && (status  = esl_strdup(sq->desc, -1, &(hit->desc)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
      } else {
        if ((status  = esl_strdup(om->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(om->acc,  -1, &(hit->acc)))   != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(om->desc, -1, &(hit->desc)))  != eslOK) esl_fatal("allocation failure");
      } 
      hit->window_length = 0;
      hit->sortkey       = pli->inc_by_E ? -lnP : seq_score; /* per-seq output sorts on bit score if inclusion is by score  */
      
      hit->score      = seq_score;
      hit->pre_score  = pre_score;
      hit->sum_score  = sum_score;

      hit->lnP        = lnP;
      hit->pre_lnP    = esl_exp_logsurv (hit->pre_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
      hit->sum_lnP    = esl_exp_logsurv (hit->sum_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

      hit->fwd_score    = fwdsc;        // from the local filter 
      hit->spfwd_score  = spfwdsc;      // from the sparse global/local model
      hit->time         = (double)time;

      p7_sparsemx_ExpectedDomains(pli->sxd, 1, sq->n, &(hit->nexpected));
      hit->noverlaps  = noverlaps;
      hit->ndom       = pli->tr->ndom;

      hit->flags       = 0;
      hit->nreported   = 0;
      hit->nincluded   = 0;
      hit->best_domain = best_d;

      hit->seqidx       = -1; /* nhmmer/longtarget only */
      hit->subseq_start = -1; /* nhmmer/longtarget only */
      
      hit->dcl        = dcl;
      dcl             = NULL;

      /* If we're using model-specific bit score thresholds (GA | TC |
       * NC) and we're in an hmmscan pipeline (mode = p7_SCAN_MODELS),
       * then we *must* apply those reporting or inclusion thresholds
       * now, because this model is about to go away; we won't have
       * its thresholds after all targets have been processed.
       * 
       * If we're using E-value thresholds and we don't know the
       * search space size (Z_setby or domZ_setby =
       * p7_ZSETBY_NTARGETS), we *cannot* apply those thresholds now,
       * and we *must* wait until all targets have been processed
       * (see p7_tophits_Threshold()).
       * 
       * For any other thresholding, it doesn't matter whether we do
       * it here (model-specifically) or at the end (in
       * p7_tophits_Threshold()). 
       * 
       * What we actually do, then, is to set the flags if we're using
       * model-specific score thresholds (regardless of whether we're
       * in a scan or a search pipeline); otherwise we leave it to 
       * p7_tophits_Threshold(). p7_tophits_Threshold() is always
       * responsible for *counting* the reported, included sequences.
       * 
       * [xref J5/92]
       */
      if (pli->use_bit_cutoffs)
      {
        if (p7_pipeline_TargetReportable(pli, hit->score, hit->lnP))
        {
          hit->flags |= p7_IS_REPORTED;
          if (p7_pipeline_TargetIncludable(pli, hit->score, hit->lnP))
            hit->flags |= p7_IS_INCLUDED;
        }

        for (d = 0; d < hit->ndom; d++)
        {
          if (p7_pipeline_DomainReportable(pli, hit->dcl[d].bitscore, hit->dcl[d].lnP))
          {
            hit->dcl[d].is_reported = TRUE;
            if (p7_pipeline_DomainIncludable(pli, hit->dcl[d].bitscore, hit->dcl[d].lnP))
              hit->dcl[d].is_included = TRUE;
          }
        }
      }
    }
  if (pli->mode == p7_SCAN_MODELS && pli->hfp) p7_profile_Destroy(gm); /* DON'T free it if we're in hmmpgmd, with all models cached */

  p7_hmm_Destroy(ehmm);   
  if (oR == NULL) p7_RateDestroy(R);
  if (r) esl_randomness_Destroy(r);
  return eslOK;
  
 NOHIT:
  p7_hmm_Destroy(ehmm);   
  if (oR == NULL) p7_RateDestroy(R);
  if (r) esl_randomness_Destroy(r);
  return eslOK;
  
 ERROR:
  if (dcl)  free(dcl);
  if (ehmm) p7_hmm_Destroy(ehmm);
  if (oR == NULL) p7_RateDestroy(R);
  if (r)    esl_randomness_Destroy(r);
  return status;
}


/* Temporary workaround for a problem in the hmmscan version of the pipeline.
 * hmmscan reads vectorized profile in two pieces from .h3f, .h3p files.
 * In H3.0 we only needed vectorized profile in the pipeline.
 * In H3.1 we need both <om> and <gm>.
 * hmmscan doesn't have <gm>.
 * Options include:
 *   1. Convert <om> to <gm>.
 *      - We'd have to recalculate the BGMk and DGkE entry/exit wing retractions
 *        from the profile; normally modelconfig does this from the HMM.
 *      - This might be slow.
 *   2. Store the profile on disk too.
 *      - I like this option, but I'd rather do it together with some
 *        other reengineering. We have a lot of redundancy in HMM, PROFILE,
 *        and OPROFILE, particularly in annotation. Should consider
 *        having a P7_MODEL container around P7_MODELINFO (annotation),
 *        P7_HMM, P7_PROFILE, vector MSV, vector VF, vector FB subparts, with ways to
 *        convert amongst HMM, PROFILE, MSV, VF, and FB sections, or read
 *        them from disk. Support delayed read of annotation including
 *        name/desc, allow HMMs (and target seqs?) to be numbered for
 *        efficiency. Measure memory footprints and timing of read, MPI 
 *        transmit, conversion.
 *   3. As it happens, we do have the HMM on disk, in the .h3m file.
 *      We can read it, and convert it to a profile.
 *      This might be slow - especially since we need to alloc/dealloc
 *      the HMM and profile in the pipeline, rather than reusing them.
 *      
 * (3) is the fastest option to implement,
 * and right now the pressure is to get 3.1 compiling and running asap;
 * we can optimize/polish later from a baseline implementation.
 */
static int
workaround_get_profile(P7_PIPELINE *pli, const P7_OPROFILE *om, const P7_BG *bg, P7_PROFILE **ret_gm)
{
  
  P7_HMM     *hmm = NULL;
  P7_PROFILE *gm  = NULL;
  int         status;

  if ( (status = p7_hmmfile_Position(pli->hfp, om->offs[p7_MOFFSET])) != eslOK) goto ERROR; /* {eslESYS | eslEINVAL} */
  if ( (status = p7_hmmfile_Read(pli->hfp, (ESL_ALPHABET **) &(om->abc), &hmm)) != eslOK) goto ERROR; /* eslEOF | eslEINCOMPAT; {eslEMEM | eslESYS} */
  /* the ESL_ALPHABET ** cast was to get rid of a const; safe, but ugly */

  if ( (    gm = p7_profile_Create(hmm->M, om->abc))                  == NULL)  { status = eslEMEM; goto ERROR; }
  if ( (status = p7_profile_Config(gm, hmm, bg))                      != eslOK) goto ERROR;
  if ( (status = p7_profile_SetLength(gm, om->L))                     != eslOK) goto ERROR;
  *ret_gm = gm;

  p7_hmm_Destroy(hmm);
  return eslOK;
  
 ERROR:
  if (hmm) p7_hmm_Destroy(hmm);
  if (gm)  p7_profile_Destroy(gm);
  return status;
}

static int
workaround_get_starprofile(P7_PIPELINE *pli, const P7_BG *bg, const P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_HMM **ret_ehmm, int noevo)
{  
  P7_HMM      *ehmm = NULL;
  int          status;

  if (pli->mode == p7_SCAN_MODELS) return eslOK;
  
  if (!noevo) {
    /* Construct evolved HMM */
    ehmm = p7_hmm_Clone(hmm);
    
    /* Convert to an optimized model */
    if ( (status = p7_profile_Config(gm, ehmm, bg))       != eslOK) goto ERROR;
    if ( (status = p7_profile_SetLength(gm, om->L))       != eslOK) goto ERROR;
    if ( (status = p7_oprofile_Convert(gm, om))           != eslOK) goto ERROR;      
    if ( (status = p7_oprofile_ReconfigLength(om, om->L)) != eslOK) goto ERROR;
  }
  
  *ret_ehmm = ehmm;
  return eslOK;
  
 ERROR:
  if (ret_ehmm) p7_hmm_Destroy(ehmm);
  return status;
}

static int
workaround_evolve_profile(double time, int n, const P7_RATE *R, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, int unihit, int verbose)
{  
  int status;
  int len = 1 * n;

  if (R == NULL) return eslOK;
  
  /* evolved HMM */
  if ( (status = p7_EvolveFromRate(NULL, hmm, R, bg, time, 0.0001, NULL, verbose)) != eslOK) goto ERROR; 
  
  /* evolved profiles gm and om */
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, len);
  if (om != NULL) {
    p7_oprofile_Convert(gm, om);    
    p7_profile_ConfigCustom(gm, hmm, bg, 0, 1.0, 0.5);  /*   ER test */ 
    p7_profile_SetLength(gm, len);                        /* er: need to restet the length here */
    p7_oprofile_ReconfigLength(om, len);
  }
  if (unihit) {
    p7_profile_ConfigCustom(gm, hmm, bg, 0, 0.0, 0.5); /*    gm *unihit* glocal/local */
    p7_profile_SetLength(gm, len);                       /* er: need to restet the length here */
  }
  
  return eslOK;
  
 ERROR:
  return status;
}

static int
workaround_calibrate(ESL_RANDOMNESS *r, int n, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om)
{  
  p7_Calibrate(hmm, NULL, &r, &bg, NULL, NULL);

 /* evolved profiles gm and om */
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, n);
  p7_oprofile_Convert(gm, om);      
  p7_oprofile_ReconfigLength(om, n);
  
  return eslOK;
}

static int
workaroud_scanmodels(P7_PIPELINE *pli, int n, EMRATE *emR, P7_RATE **ret_R, P7_HMM **ret_ehmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, int noevo, float fixtime, float tol, char *errbuf)
{
  P7_RATE *R    = NULL;
  P7_HMM  *ehmm = NULL;
  int      status;

  if (noevo) {
    if (pli->hfp) p7_oprofile_ReadRest(pli->hfp, om);
    p7_oprofile_ReconfigRestLength(om, n);
    if ((status = p7_pipeline_NewModelThresholds(pli, om)) != eslOK) return status; /* pli->errbuf has err msg set */
    return eslOK;
  }

  /* the hmm */
  if ((status = p7_hmmfile_Read(pli->hfp, (ESL_ALPHABET **)(&om->abc), &ehmm)) != eslOK) goto ERROR;

  /* Create the hmm rate */
  if ((status = p7_RateCalculate(NULL, ehmm, bg, emR, NULL, &R, R->evomodel, R->betainf, fixtime, tol, errbuf, FALSE)) != eslOK) goto ERROR;

  if ((status = workaround_get_starprofile(pli, bg, ehmm, gm, om, &ehmm, noevo)) != eslOK) goto ERROR;

  *ret_R    = R;
  *ret_ehmm = ehmm;

  return eslOK;

 ERROR:
  if (ehmm) p7_hmm_Destroy(ehmm);
  if (R)    p7_RateDestroy(R);
  return status;
}


/* --------------------------end --------------------- */

static int
optimize_msvfilter(const ESL_DSQ *dsq, int n, float *ret_time,
		   P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx, 
		   float *ret_usc, double firststep, float fixtime, int noevo,
		   float tol, char *errbuf, int be_verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 sc;
  double                 usc_init;
  float                  time;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;
  
  time = (isfixtime)? fixtime : *ret_time;
  usc_init = func_msvfilter((ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, fx, time, noevo, tol, errbuf, be_verbose);
  if (noevo || isfixtime || usc_init == eslINFINITY) {
    *ret_usc  = usc_init;
    *ret_time = time;
    return eslOK;
  }

  np = 1;     /* variable: time */
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);
  
 /* Copy shared info into the "data" structure
   */
  data.time       = time;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.om         = (P7_OPROFILE *)om;
  data.bg         = bg;
  data.firststep  = firststep;
  data.fx         = fx;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.be_verbose = be_verbose;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, (long)np, &data);
 
  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (long)np, &data);
  status = min_Bracket(p, u, np, data.firststep,
		       &optimize_msvfilter_func,
		       (void *) (&data), 
		       tol, wrk, &sc);
  if (status != eslOK) 
    esl_fatal("optimize_msvfilter(): bad bracket minimization");	
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  data.usc = -sc;
  if (be_verbose) printf("END MSV OPTIMIZATION: time %f usc %f --> %f\n", data.time, usc_init, data.usc);
  
  *ret_usc  = data.usc;
  *ret_time = data.time;
  
  /* clean up */
  if (u   != NULL) free(u);
  if (p   != NULL) free(p);
  if (wrk != NULL) free(wrk);
  return eslOK;

 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;
}

int
p7_evopli_OptimizeViterbiFilter(const ESL_DSQ *dsq, int n, float *ret_time,
				P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx, 
				float *ret_vitsc, double firststep, int noevo, float fixtime, int unihit,
				float tol, char *errbuf, int be_verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 sc;
  double                 vitsc_init;
  float                  time;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;

  time = (isfixtime)? fixtime : *ret_time;
  vitsc_init = func_viterbifilter((ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, fx, time, unihit, noevo, tol, errbuf, be_verbose);
  if (noevo || isfixtime || vitsc_init == eslINFINITY) {
    *ret_vitsc  = vitsc_init;
    *ret_time   = time;
    return eslOK;
  }

  np = 1;     /* variables: time */

  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);
  
 /* Copy shared info into the "data" structure
   */
  data.time       = time;
  data.unihit     = unihit;
  data.noevo      = noevo;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.om         = (P7_OPROFILE *)om;
  data.bg         = bg;
  data.firststep  = firststep;
  data.fx         = fx;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.be_verbose = be_verbose;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, (long)np, &data);
 
  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (long)np, &data);
  status = min_Bracket(p, u, np, data.firststep,
		       &optimize_viterbifilter_func,
		       (void *) (&data), 
		       tol, wrk, &sc);
  if (status != eslOK) 
    esl_fatal("optimize_viterbifilter(): bad bracket minimization");	
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  data.vitsc = -sc;
  if (be_verbose) printf("END VIT OPTIMIZATION: time %f vitsc %f --> %f\n", data.time, vitsc_init, data.vitsc);
  
  *ret_vitsc = data.vitsc;
  *ret_time  = data.time;
  
  /* clean up */
  if (u   != NULL) free(u);
  if (p   != NULL) free(p);
  if (wrk != NULL) free(wrk);
  return eslOK;

 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;
}

int
p7_evopli_OptimizeForwardFilter(const ESL_DSQ *dsq, int n, float *ret_time,
				P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_CHECKPTMX *cx, 
				float *ret_fwdsc, double firststep, int noevo, float fixtime, int unihit,
				float tol, char *errbuf, int be_verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                      */
  double                *u = NULL;             /* max initial step size vector          */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches     */
  double                 sc;
  double                 fwdsc_init;
  float                  time;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;

  time = (isfixtime)? fixtime : *ret_time;
  fwdsc_init = func_forwardfilter((ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, cx, time, unihit, noevo, tol, errbuf, be_verbose);
  if (noevo || isfixtime || fwdsc_init == eslINFINITY) {
    *ret_fwdsc  = fwdsc_init;
    *ret_time   = time;
    return eslOK;
  }
  
  np = 1;     /* variables: time */
  
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);

  /* Copy shared info into the "data" structure
   */
  data.time       = time;
  data.unihit     = unihit;
  data.noevo      = noevo;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.om         = (P7_OPROFILE *)om;
  data.bg         = bg;
  data.firststep  = firststep;
  data.cx         = cx;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.be_verbose = be_verbose;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, (long)np, &data);

  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (long)np, &data);
  status = min_Bracket(p, u, np, data.firststep,
		       &optimize_forwardfilter_func,
		       (void *) (&data), 
		       tol, wrk, &sc);
  if (status != eslOK) 
    esl_fatal("optimize_forwardfilter(): bad bracket minimization");		
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  data.fwdsc = -sc;
  if (be_verbose) printf("END FWD OPTIMIZATION: time %f fwdsc %f --> %f\n", data.time, fwdsc_init, data.fwdsc);
  
  *ret_fwdsc = data.fwdsc;
  *ret_time  = data.time;
  
  /* clean up */
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return eslOK;

 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;
}


int
p7_evopli_OptimizeSparseViterbi(const ESL_DSQ *dsq, int n, float *ret_time,
				P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, 
				P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxx,  P7_TRACE *tr,
				float *ret_vitsc, double firststep, int noevo, float fixtime, int unihit,
				float tol, char *errbuf, int be_verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                      */
  double                *u = NULL;             /* max initial step size vector          */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches     */
  double                 sc;
  double                 vitsc_init;
  float                  time;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;

  time = (isfixtime)? fixtime : *ret_time;
  vitsc_init = func_sparseviterbi((ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, sm, cx, sxx, tr, time, unihit, noevo, tol, errbuf, be_verbose);
  if (noevo || isfixtime || vitsc_init == eslINFINITY) {
    *ret_vitsc  = vitsc_init;
    *ret_time   = time;
    return eslOK;
  }
  if (vitsc_init == -eslINFINITY) {
    *ret_vitsc  = vitsc_init;
    *ret_time   = time;
    return eslOK;
  }

  np = 1;     /* variables: time */
  
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);

  /* Copy shared info into the "data" structure
   */
  data.time       = time;
  data.unihit     = unihit;
  data.noevo      = noevo;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.bg         = bg;
  data.firststep  = firststep;
  data.sm         = sm;
  data.cx         = cx;
  data.sxx        = sxx;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.be_verbose = be_verbose;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, (long)np, &data);

  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (long)np, &data);
  status = min_Bracket(p, u, np, data.firststep,
		       &optimize_sparseviterbi_func,
		       (void *) (&data), 
		       tol, wrk, &sc);
  if (status != eslOK) 
    esl_fatal("optimize_sparseviterbi(): bad bracket minimization");		
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  data.vitsc = -sc;
  if (be_verbose) printf("END SPARSEVIT OPTIMIZATION: time %f vitsc %f --> %f\n", data.time, vitsc_init, data.vitsc);
  
  *ret_vitsc = data.vitsc;
  *ret_time  = data.time;
  
  /* clean up */
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return eslOK;

 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;
}

int
p7_evopli_OptimizeSparseForward(const ESL_DSQ *dsq, int n, float *ret_time,
				P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, 
				P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxf, 
				float *ret_fwdsc, double firststep, int noevo, float fixtime, int unihit,
				float tol, char *errbuf, int be_verbose)
{
  struct optimize_data   data;
  double                *p = NULL;	       /* parameter vector                      */
  double                *u = NULL;             /* max initial step size vector          */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches     */
  double                 sc;
  double                 fwdsc_init;
  float                  time;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;

  time = (isfixtime)? fixtime : *ret_time;
  fwdsc_init = func_sparseforward((ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, sm, cx, sxf, time, unihit, noevo, tol, errbuf, be_verbose);
  if (noevo || isfixtime || fwdsc_init == eslINFINITY) {
    *ret_fwdsc  = fwdsc_init;
    *ret_time   = time;
    return eslOK;
  }
   if (fwdsc_init == -eslINFINITY) {
    *ret_fwdsc  = fwdsc_init;
    *ret_time   = time;
    return eslOK;
  }
 
  np = 1;     /* variables: time */
  
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);

  /* Copy shared info into the "data" structure
   */
  data.time       = time;
  data.unihit     = unihit;
  data.noevo      = noevo;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.bg         = bg;
  data.firststep  = firststep;
  data.sm         = sm;
  data.cx         = cx;
  data.sxf        = sxf;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.be_verbose = be_verbose;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, (long)np, &data);

  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (long)np, &data);
  status = min_Bracket(p, u, np, data.firststep,
		       &optimize_sparseforward_func,
		       (void *) (&data), 
		       tol, wrk, &sc);
  if (status != eslOK) 
    esl_fatal("optimize_sparsefwderbi(): bad bracket minimization");		
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  data.fwdsc = -sc;
  if (be_verbose) printf("END SPARSEFWD OPTIMIZATION: time %f fwdsc %f --> %f\n", data.time, fwdsc_init, data.fwdsc);
  
  *ret_fwdsc = data.fwdsc;
  *ret_time  = data.time;
  
  /* clean up */
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return eslOK;

 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;
}


static int
optimize_pack_paramvector(double *p, long np, struct optimize_data *data)
{
  int   x = 0;
  
  p[x] = (data->time < 1.0)? log(data->time) : data->time - 1.0;

  return eslOK;  
}


static int
optimize_unpack_paramvector(double *p, long np, struct optimize_data *data)
{
  float time;
  float tmax = 10.0;
  int   x = 0;
  
  time = (p[x] < 0.0)? exp(p[x]) : p[x] + 1.0; 
  if (time > tmax) time = tmax;
  
  data->time = time;
  return eslOK;
}

static void
optimize_bracket_define_direction(double *u, long np, struct optimize_data *data)
{
  int x;
  for (x = 0; x < np; x++) u[x] = 0.25;
  u[np] = 0.25;
}

static double
optimize_msvfilter_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, np, data);
  
  data->usc = func_msvfilter(dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->fx, 
			     data->time, data->noevo, data->tol, data->errbuf, data->be_verbose);
  
  if (data->usc == eslINFINITY) data->usc = 1000.;
  return -(double)data->usc;
}

static double
optimize_viterbifilter_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, np, data);
  
  data->vitsc = func_viterbifilter(dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->fx, 
				   data->time, data->unihit, data->noevo, data->tol, data->errbuf, data->be_verbose);
  
  if (data->vitsc == eslINFINITY) data->vitsc = 1000.;
  return -(double)data->vitsc;
}

static double
optimize_forwardfilter_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, np, data);

  data->fwdsc = func_forwardfilter(dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->cx, 
				   data->time, data->unihit, data->noevo, data->tol, data->errbuf, data->be_verbose);
  
  return -(double)data->fwdsc;
}

static double
optimize_sparseviterbi_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, np, data);
  
  data->vitsc = func_sparseviterbi(dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->sm, data->cx, data->sxx, NULL, 
				   data->time, data->unihit, data->noevo, data->tol, data->errbuf, data->be_verbose);
  
  if (data->vitsc == eslINFINITY) data->vitsc = 1000.;
  return -(double)data->vitsc;
}

static double
optimize_sparseforward_func(double *p, long np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, np, data);
  
  data->fwdsc = func_sparseforward(dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->sm, data->cx, data->sxf,
				   data->time, data->unihit, data->noevo, data->tol, data->errbuf, data->be_verbose);
  
  if (data->fwdsc == eslINFINITY) data->fwdsc = 1000.;
  return -(double)data->fwdsc;
}


static double
func_msvfilter(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx, 
	       float time, int noevo, float tol, char *errbuf, int verbose)
{
  float  usc;
  
  if (!noevo) {
    /* Construct evolved profile */
    if (workaround_evolve_profile((double)time, n, R, bg, hmm, gm, om, FALSE, FALSE) != eslOK) exit(1);
    
    /* Convert to an optimized profile */
    p7_filtermx_Reuse(fx);
  }
  
  p7_MSVFilter(dsq, n, om, fx, &(usc));
  
#if 0
  printf("time %f usc %f\n", time, usc);
#endif

  return (double)usc;
 }

static double
func_viterbifilter(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx, 
		   float time, int unihit, int noevo, float tol, char *errbuf, int verbose)
{
  float  vitsc;
  
  if (!noevo) {
    /* Construct evolved profile */
    if (workaround_evolve_profile((double)time, n, R, bg, hmm, gm, om, unihit, FALSE) != eslOK) exit(1);
    
    /* Convert to an optimized profile */
    p7_filtermx_Reuse(fx);
  }

  p7_ViterbiFilter(dsq, n, om, fx, &(vitsc));
  
#if 0
  printf("time %f vitsc %f\n", time, vitsc);
#endif

  return (double)vitsc;
 }

static double
func_forwardfilter(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_CHECKPTMX *cx, 
		   float time, int unihit, int noevo, float tol, char *errbuf, int verbose)
{
  float   fwdsc;

  if (!noevo) {
    /* Construct evolved profile */
    if (workaround_evolve_profile((double)time, n, R, bg, hmm, gm, om, unihit, FALSE) != eslOK) exit(1);
    
    p7_checkptmx_Reuse(cx);
  }
  
  p7_ForwardFilter(dsq, n, om, cx, &(fwdsc)); 

#if 0
  printf("time %f fwdsc %f\n", time, fwdsc);
#endif

  return (double)fwdsc;
 }

static double 
func_sparseviterbi(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, 
		   P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxx, P7_TRACE *tr, 
		   float time, int unihit, int noevo, float tol, char *errbuf, int verbose)
{
  float  vitsc;
    
  if (!noevo) {    
    /* Construct evolved profile */
    if (workaround_evolve_profile((double)time, n, R, bg, hmm, gm, om, unihit, FALSE) != eslOK) exit(1);

    if (tr) p7_trace_Reuse(tr);
    p7_checkptmx_Reuse(cx);
    p7_ForwardFilter(dsq, n, om, cx, NULL); 
  }
  
  p7_BackwardFilter(dsq, n, om, cx, sm, p7_SPARSIFY_THRESH);
  p7_SparseViterbi(dsq, n, gm, sm, sxx, tr, &vitsc);
  
#if 0
  printf("time %f sparsevitsc %f noevo? %d unihit? %d gm->nj %f\n", time, vitsc, noevo, unihit, gm->nj);
#endif

  return (double)vitsc;
}

static double 
func_sparseforward(ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, 
		   P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxf,
		   float time, int unihit, int noevo, float tol, char *errbuf, int verbose)
{
  float  fwdsc;
  
  if (!noevo) {
    /* Construct evolved profile */
    if (workaround_evolve_profile((double)time, n, R, bg, hmm, gm, om, unihit, FALSE) != eslOK) exit(1);
  
    p7_checkptmx_Reuse(cx);
    p7_ForwardFilter (dsq, n, om, cx, NULL); 
    p7_BackwardFilter(dsq, n, om, cx, sm, p7_SPARSIFY_THRESH);
  }
  
  p7_SparseForward(dsq, n, gm, sm, sxf, &fwdsc);
  
#if 0
  printf("time %f sparseforward %f unihit? %d gm->nj %f\n", time, fwdsc, unihit, gm->nj);
#endif

  return (double)fwdsc;
}

/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk/src/search/p7_pipeline.c $
 * SVN $Id: p7_evopipeline.c 4617 2014-02-21 18:06:34Z eddys $
 *****************************************************************/
