/* P7_EVOPIPELINE is the standardized pipeline for one profile/sequence
 * comparison, from the fast filters down through domain postprocessing,
 * alignment, and scoring.
 */

#ifndef P7_EVOPIPELINE_INCLUDED
#define P7_EVOPIPELINE_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "base/p7_bg.h"
#include "base/p7_hmmfile.h"
#include "base/p7_hmmwindow.h"
#include "base/p7_masstrace.h"
#include "base/p7_scoredata.h"
#include "base/p7_tophits.h"
#include "base/p7_trace.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_reference/p7_refmx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/p7_filtermx.h"

#include "evohmmer.h"

struct optimize_data {
  float           time;
  int             unihit;
  int             noevo;
  ESL_DSQ        *dsq;
  int             n;
  P7_RATE        *R;
  P7_HMM         *hmm;
  P7_PROFILE     *gm;
  P7_OPROFILE    *om;
  P7_BG          *bg;
  P7_FILTERMX    *fx;
  P7_CHECKPTMX   *cx;
  P7_SPARSEMASK  *sm;
  P7_SPARSEMX    *sxx;
  P7_SPARSEMX    *sxf;
  P7_SPARSEMX    *sxb;
  float          usc;
  float          vitsc;
  float          fwdsc;
  double         firststep;
  double         tol;
  char          *errbuf;
  int            be_verbose;
};

extern int p7_evopli_OptimizeViterbiFilter(const ESL_DSQ *dsq, int n, float *ret_time, 
					   P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_FILTERMX *fx, 
					   float *ret_vitsc, double firststep, int noevo, float fixtime, int unihit, float tol, char *errbuf, int be_verbose);
extern int p7_evopli_OptimizeForwardFilter(const ESL_DSQ *dsq, int n, float *ret_time, 
					   P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_CHECKPTMX *cx, 
					   float *ret_vitsc, double firststep, int noevo, float fixtime, int unihit, float tol, char *errbuf, int be_verbose);
extern int p7_evopli_OptimizeSparseViterbi(const ESL_DSQ *dsq, int n, float *ret_time,
					   P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, 
					   P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxx,  P7_TRACE *tr,
					   float *ret_vitsc, double firststep, int noevo, float fixtime, int unihit,
					   float tol, char *errbuf, int be_verbose);
extern int p7_evopli_OptimizeSparseForward(const ESL_DSQ *dsq, int n, float *ret_time,
					   P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, 
					   P7_SPARSEMASK *sm, P7_CHECKPTMX *cx, P7_SPARSEMX *sxf,  
					   float *ret_fwdsc, double firststep, int noevo, float fixtime, int unihit,
					   float tol, char *errbuf, int be_verbose);
extern int p7_EvoPipeline(P7_PIPELINE *pli, EMRATE *emR, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *th, 
			  int noevo, float fixtime, int unihit);

#endif /*P7_EVOPIPELINE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
