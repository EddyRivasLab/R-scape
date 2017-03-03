#ifndef P7_EVOTRACEALIGN_INCLUDED
#define P7_EVOTRACEALIGN_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_msa.h"
#include "esl_sq.h"

#include "base/p7_hmm.h"
#include "base/p7_trace.h"

#include "evohmmer.h"

extern int p7_evotracealign_ComputeTraces(P7_RATE *R, P7_HMM *hmm, ESL_SQ  **sq, int offset, int N, P7_TRACE  **tr, P7_BG *bg, int noevo, float fixtime, float tol, char *errbuf, int verbose);

#endif /*P7_EVOTRACEALIGN_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
