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

extern int struct_CYKCOV(struct data_s *data, ESL_MSA *msa, int *ret_nct, int ***ret_cykctlist, int minloop, RANKLIST *ranklist, HITLIST *hitlist, enum grammar_e G, THRESH *thresh);
extern int struct_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int nct, int **ctlist, struct mutual_s *mi, int *msamap, int firstpos,
				    SAMPLESIZE samplesize,  HITLIST *hitlist, int dosvg, int verbose, char *errbuf);
extern int struct_ExpandCT(char *r2rfile, int r2rall,  ESL_RANDOMNESS *r, ESL_MSA *msa, int *ret_nct, int ***ret_ctlist,
				     int minloop, enum grammar_e G, int verbose, char *errbuf);
extern int struct_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, enum grammar_e G, int minloop, int verbose, char *errbuf);
extern int struct_SplitCT(int *ct, int L, int *ret_nct, int ***ret_ctlist, int verbose);

#endif
