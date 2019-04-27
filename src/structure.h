/* structure.h
 *
 *   
*/
#ifndef STRUCTURE_INCLUDED
#define STRUCTURE_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"

#include "covgrammars.h"
#include "correlators.h"

// folding parameters
//
// power threshold
#define POWER_THRESH 0.80

// paramters to include extra helices
#define  INCOMPFRAC     0.51          // max fraction of residues in a helix that overlap with another existing helix in order to be removed
#define  MINHELIX       15            // min length of a helix without any covarying basepairs in order to be reported
#define  LASTFOLD       1
#define  HELIX_UNPAIRED 2             // max number of unpaired residues in the definition of a helix

typedef struct cov_s {
  int64_t i;
  int64_t j;
  
  int64_t nsubs;
  double  power;
  double  score;

  int     isbp;
  
} COV;

typedef struct covlist_s {
  int64_t  n;
  COV     *cov;
  
} COVLIST;

extern int struct_COCOMCYK(struct data_s *data, ESL_MSA *msa, int *ret_nct, int ***ret_cykctlist, int minloop,
			   RANKLIST *ranklist, HITLIST *hitlist, enum grammar_e G, THRESH *thresh);
extern int struct_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int nct, int **ctlist, struct mutual_s *mi, int *msamap, int firstpos,
			  SAMPLESIZE samplesize,  HITLIST *hitlist, int dosvg, int verbose, char *errbuf);
extern int struct_SplitCT(int *ct, int L, int *ret_nct, int ***ret_ctlist, int verbose);
extern int struct_CTMAP(int L, int nct, int **ctlist, int OL, int *msamap, int ***ret_octlist, char ***ret_sslist, FILE *fp, int verbose);
#endif
