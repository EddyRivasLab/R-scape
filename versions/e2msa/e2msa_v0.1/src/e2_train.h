/* e2_train
 *
 */
#ifndef E2_TRAIN_INCLUDED
#define E2_TRAIN_INCLUDED


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

#define PARAM2RATE(p) (fabs((p)))
#define RATE2PARAM(r) ((r))

#define PARAM2BERNLOG(p) ( ((p) < 100)? exp((p))/(1.0+exp((p))) : 1.0 )
#define BERN2PARAMLOG(q) ( ((q) < 1.0)? (log(q) - log(1.0-(q))) : -eslINFINITY )

#define PARAM2BERNOULLI(p) ( fabs((p)) / (1.0 + fabs((p))) )
#define BERNOULLI2PARAM(q) ( ((q) < 1.0)? ((q) / (1.0 - (q))) : 0.0 )

struct e2_train_data {
  ESL_RANDOMNESS *r;
  E1_RATE        *R;
  E1_BG          *bg;
  float           sc;
  E2_PIPELINE    *pli;
  E2_ALI          e2ali;
  int             mode;
  int             do_viterbi;
  int             nmsa;
  ESL_TREE     **Tlist;
  ESL_MSA      **msalist;
  float        **msafrq;
  int            it;
  char           *errbuf;
  float           tol;
  int             verbose;
};

extern int e2_train(ESL_RANDOMNESS *r, int nmsa, ESL_TREE **Tlist, ESL_MSA **msalist, float **msafrq, E1_RATE *R, E2_PIPELINE *pli, E1_BG *bg, int mode, int do_viterbi,
		    int amoeba, double tol, char *errbuf, int verbose);
extern int e2_transitions_pack_paramvector  (int *ret_x, double *p, long np, E1_RATE *R, char *errbuf, int verbose);
extern int e2_transitions_unpack_paramvector(int *ret_x, double *p, long np, E1_RATE *R, char *errbuf, int verbose);

#endif /*E2_TRAIN_INCLUDED*/


/************************************************************
 * @LICENSE@
 ************************************************************/
