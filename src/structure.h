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

// some folding parameters
//
// power threshold
#define POWER_THRESH 0.80

// paramters to include extra helices
#define  INCOMPFRAC  0.20          // max fraction of residues in a helix that overlap with another existing helix in order to be removed
#define  MINHELIX    4             // min length of a helix without any covarying basepairs
#define  LASTFOLD    1

/* LENIENT = 1  The covarying pairs that have been explained already are not allowed to pair with each other again but 
 *              could be pair with some other residue. Would allow to see triplets
 *
 * LENINEN = 0  The already explained covarying residues are forced to remaoin unpared in cacading foldings
 */
#define  LENIENT 0


extern int struct_COCOMCYK(struct data_s *data, ESL_MSA *msa, int *ret_nct, int ***ret_cykctlist, int minloop,
			   RANKLIST *ranklist, HITLIST *hitlist, enum grammar_e G, THRESH *thresh);
extern int struct_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int nct, int **ctlist, struct mutual_s *mi, int *msamap, int firstpos,
			  SAMPLESIZE samplesize,  HITLIST *hitlist, int dosvg, int verbose, char *errbuf);
extern int struct_SplitCT(int *ct, int L, int *ret_nct, int ***ret_ctlist, int verbose);
extern int struct_CTMAP(int L, int nct, int **ctlist, int OL, int *msamap, int ***ret_octlist, char ***ret_sslist, FILE *fp, int verbose);
#endif
