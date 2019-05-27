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
// TRUE if we do one last fold without covariations
#define  LASTFOLD    0

// power threshold
#define POWER_THRESH 0.95

// parameters for the main nested structure
#define HLOOP_MIN 1                  // minimum length of a hairpin loop. If i-j is the closing pair: i-x-x-j is minhloop = 3
                                     // unless there are covariations forcing a smaller hairpin loop.

// parameters to break non-nested structures in helices
#define  HELIX_UNPAIRED 2            // max number of unpaired residues in a non-nested helix

// parameters for selecting non-nested helices without covariations
#define  OVERLAPFRAC       0.51      // max fraction of paired residues that overlap with another existing helix in order to be removed
#define  MINHELIX          15        // min length to be reported

// special parameter for selecting helices with covariations
//
// use these setting for maximal display of basepairing even if overlaping or contiguous
// REMOVE_CONTCOV     0
// HELIX_OVERLAP_TRIM 0
//
#define REMOVE_CONTCOV     1         // if TRUE it does not display covariations between contiguous residues
#define HELIX_OVERLAP_TRIM 0         // TRUE for trimming non-nested helices with covariations to remove ovelap with the main non-nested structure

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

typedef struct pair_s {
  int64_t i;
  int64_t j;
} PAIR;
typedef struct pairlist_s {
  int64_t   n;
  PAIR     *pair;
  
} PAIRLIST;

enum fold_e {
  CYK,
  DECODING,
};

extern int struct_CACOFOLD(struct data_s *data, ESL_MSA *msa, int *ret_nct, int ***ret_cykctlist, 
			   RANKLIST *ranklist, HITLIST *hitlist, enum grammar_e G, enum fold_e F, THRESH *thresh);
extern int struct_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int nct, int **ctlist, struct mutual_s *mi, int *msamap, int firstpos,
			  SAMPLESIZE samplesize,  HITLIST *hitlist, int dosvg, int verbose, char *errbuf);
extern int struct_SplitCT(int *ct, int L, int *ret_nct, int ***ret_ctlist, int verbose);
extern int struct_CTMAP(int L, int nct, int **ctlist, int OL, int *msamap, int ***ret_octlist, char ***ret_sslist, FILE *fp, int verbose);
#endif
