/* msatree - funtions to build a tree from an alignment
 *
 */
#ifndef BRANCHINFERENCE_INCLUDED
#define BRANCHINFERENCE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_tree.h"

#include "e2.h"
#include "e1_rate.h"
#include "e2_pipeline.h"
#include "fetchpfamdb.h"

struct bf_data {
  ESL_RANDOMNESS *r;
  ESL_MSA        *msa;
  ESL_TREE       *T;
  ESL_TREE       *Ti;
  E2_PIPELINE    *pli;
  E2_ALI          e2ali;
  E1_RATE        *R;
  P7_RATE        *R7;
  float           tsat;
  E1_BG          *bg;
  P7_BG          *bg7;
  float           tol;
  char           *errbuf;
  int             verbose;
};

extern int ehmmBranchInference(FILE *outfp, char *gnuplot, char *distfile, char *histofile, char *ratefile, char *binratefile,
			       ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, SPECIES *SPE, P7_PIPELINE *pli, P7_BG *bg,
			       float evalcutoff, P7_HMM *hmm, int noopt, int do_viterbi, float incbintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose);
extern int e2BranchInference(FILE *outfp, char *gnuplot, char *distfile, char *histofile, char *ratefile, char *binratefile,
			     ESL_RANDOMNESS *r, ESL_MSA *msa, float *msafrq, ESL_TREE *T, SPECIES *SPE, E2_PIPELINE *pli, E1_BG *bg, P7_BG *bg7, E2_ALI e2ali,
			     float evalcutoff, E1_RATE *R, P7_RATE *R7, int mode, int do_viterbi, float incbintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose);
extern int binf_InferedDistances(FILE *outfp,  char *histofile, char *ratefile, char *binratefile, char *gnuplot, char *ssidistfile, 
				 const ESL_MSA *msa, ESL_TREE *T, ESL_TREE *Ti, SPECIES *SPE, 
				 char **taxalist, int ntaxa, float maxbintime, float incbintime, float tlinear, int full, 
				 float evalcutoff, int view, char *errbuf, int verbose); 

#endif /*BINF_INCLUDED*/


/************************************************************
 * @LICENSE@
 ************************************************************/
