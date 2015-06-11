/* covariation.h
 *
 *   
*/
#ifndef COVARIATION_INCLUDED
#define COVARIARITION_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_tree.h"
#include "ribosum_matrix.h"

#include "covgrammars.h"
#include "covmeasures.h"

#define W     0.1     // COV with
#define BMAX  2.0     // max COV score per position
#define BMIN -2.0     // min COV score per position



typedef enum {
  NullNONE = 0,
  Null1    = 1,
  Null2    = 2,
  Null3    = 3,
} NULLTYPE;


typedef struct ranklist_s {
  int       nb;            /* number of bins                                  */
  double    w;	    	   /* fixed width of each bin                         */
  double    bmin, bmax;	   /* sc bounds: all sc satisfy bmin < sc <= bmax     */
  double    scmin, scmax;  /* smallest, largest sample value sc observed      */
 
  double   *covBP;
  double   *covNBP;

} RANKLIST;

typedef struct hit_s {
  int64_t i;
  int64_t j;
  
  double sc;
  double expcovNBP;

  int is_bpair;
  int is_compatible;
} HIT;

typedef struct hitlist_s{
  int  nhit;
  HIT *hit;

}  HITLIST;

extern int              COV_SignificantPairs_Ranking(RANKLIST **ret_ranklist, HITLIST **ret_hitlist, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP,
							int nbpairs, int verbose, char *errbuf);
extern RANKLIST        *COV_CreateRankList(int L, double bmax, double bmin, double w);
extern int              COV_CreateHitList(FILE *fp, HITLIST **ret_hitlist, double threshsc, struct mutual_s *mi, int *msamap, int *ct, RANKLIST *ranklist, 
					     int verbose, char *errbuf);
extern void             COV_FreeRankList(RANKLIST *ranklist);
extern void             COV_FreeHitList(HITLIST *hitlist);
extern int              COV_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int *ct, int verbose, char *errbuf);
extern int              COV_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen);
extern int              COV_CYKCOVCT(FILE *outfp, char *gnuplot, char *dplotfile, char *R2Rcykfile, char *R2Rversion, int R2Rall, ESL_RANDOMNESS *r, 
					ESL_MSA **msa, struct mutual_s *mi, int *msamap, int minloop, enum grammar_e G, 
					int maxFP, double expectFP, int nbpairs, char *errbuf, int verbose);
extern int              COV_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int *ct, struct mutual_s *mi, int *msamap, HITLIST *hitlist, int verbose, char *errbuf);
extern int              COV_R2R(char *r2rfile, char *r2rversion, int r2rall, ESL_MSA **msa, int *ct, int *msamap, HITLIST *hitlist, int makepdf, 
				   int verbose, char *errbuf);
extern int              COV_R2Rpdf(char *r2rfile, char *r2rversion, int verbose, char *errbuf);
extern int              COV_ExpandCT(char *r2rfile, int r2rall,  ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, int minloop, enum grammar_e G, int verbose, char *errbuf);
extern int              COV_ExpandCT_Naive(ESL_MSA *msa, int *ct, int minloop, int verbose, char *errbuf);
extern int              COV_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ct, enum grammar_e G, int minloop, int verbose, char *errbuf);
  

#endif
