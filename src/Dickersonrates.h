/* msatree - funtions to build a tree from an alignment
 *
 */
#ifndef DICKERSONRATES_INCLUDED
#define DICKERSONRATES_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_tree.h"

#include "e2.h"
#include "fetchpfamdb.h"


extern int e2DickersonChangeRates(FILE *outfp, char *gnuplot, char *distfile, char *histofile, char *ratefile, char *binratefile,
				  ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, SPECIES *SPE, E2_PIPELINE *pli, E2_ALI e2ali, 
				  float evalcutoff, E1_RATE *R,  P7_RATE *R7, E1_BG *bg, P7_BG *bg7, 
				  float incbintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose);
extern int e2DickersonRatesCluster(FILE *outfp,  char *histofile, char *ratefile, char *biratefile, char *gnuplot, char *ssidistfile,
				   const ESL_MSA *msa, const ESL_TREE *T, const E2_TRACE **Tr, const PSQ **sqn, 
				   SPECIES *SPE, char **taxalist, int ntaxa, float maxbintime, float incbintime, float tlinear, int full, 
				   float evalcutoff, int view, char *errbuf, int verbose);
extern int naiveDickersonRates(const ESL_ALPHABET *abc, const ESL_DSQ *asq1, const ESL_DSQ *asq2, int64_t L, float *ret_chrate, float *ret_cchrate);
extern int e2DickersonRates(int n, int m, int lca, const ESL_TREE *T, const ESL_MSA *msa, const E2_TRACE **tr, const PSQ **sqn, 
			    float *ret_e2srate, float *ret_e2drate, float *ret_e2brate, float *ret_e2irate);

#endif /*DICKERSONRATES_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
