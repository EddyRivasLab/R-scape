/* e2_tree.h
 *
 */
#ifndef E2_TREE_INCLUDED
#define E2_TREE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_tree.h"

#include "hmmer.h"

#include "e2.h"
#include "e1_bg.h"
#include "e1_rate.h"
#include "e2_pipeline.h"
#include "evohmmer.h"

extern int e2_tree_UPGMA(ESL_TREE **ret_T, int n, ESL_SQ **seq, ESL_MSA *msa, float *frq, ESL_RANDOMNESS *r, E2_PIPELINE *pli,
			 E1_RATE *R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7, 
			 E2_ALI e2ali, int mode, int do_viterbi, float fixtime, float tinit, double tol, char *errbuf, int verbose);
#endif /*E2_TREE_INCLUDED*/


/************************************************************
 * @LICENSE@
 ************************************************************/
