/* covariation.h
 *
 *   
*/
#ifndef MUTUALINFO_INCLUDED
#define MUTUALINFO_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_tree.h"
#include "ribosum_matrix.h"

#include "covgrammars.h"


#define NCOVTYPE = 5;
typedef enum {
  CHI   = 0,
  CHIp  = 1,
  CHIa  = 2,

  GT    = 3,
  GTp   = 4,
  GTa   = 5,

  MI    = 6,
  MIp   = 7,
  MIa   = 8,

  MIr   = 9,
  MIrp  = 10,
  MIra  = 11,

  MIg   = 12,
  MIgp  = 13,
  MIga  = 14,

  OMES  = 15,
  OMESp = 16,
  OMESa = 17,

  COVALL  = 18,
  COVNONE = 19,
} COVTYPE;

typedef enum {
  NullNONE = 0,
  Null1    = 1,
  Null2    = 2,
  Null3    = 3,
} NULLTYPE;

typedef enum {
  C16      = 0,
  C2       = 1,
  CSELECT  = 2,
} COVCLASS;

typedef enum {
  APC = 0,
  ASC = 1,
} CORRTYPE;

typedef enum{
  NAIVE  = 0,
  PHYLO  = 1,
  DCA    = 2,
  AKMAEV = 3,
} METHOD;

struct mutual_s {
  int64_t         alen;
  int64_t         nseq;
  double       ***pp;    // joint probability of two position [0,alen-1][0.alen-1][0..15]
  double        **pm;    // marginal probabilities [0,alen-1][0..3]
  int           **nseff; // effective number of sequences  [0,alen-1][0,alen-1]
  int           **ngap;  // number of gaps  [0,alen-1][0,alen-1]

  COVTYPE         type;
  COVCLASS        class;
  ESL_DMATRIX    *COV;    // covariation matrix (MI, MIp, MIr, MIa, CHI,...)  mutual information
 
  double          besthreshCOV;
  double          minCOV;
  double          maxCOV;

  int             ishuffled;
  int             nseqthresh; // if nseq <= nseqthresh use C2 method otherwise use C16

  ESL_ALPHABET   *abc;
};

typedef struct ranklist_s {
  int     N;
  double *sc;
  int    *covBP;
  int    *covNBP;

} RANKLIST;

typedef struct hit_s {
  int64_t i;
  int64_t j;
  
  double sc;
  double expcovNBP;

  double covBP_avg_null1;
  double covBP_avg_null2;
  double covBP_avg_null3;
  double covBP_std_null1;
  double covBP_std_null2;
  double covBP_std_null3;

  int is_bpair;
  int is_compatible;
} HIT;

typedef struct hitlist_s{
  int  nhit;
  HIT *hit;

}  HITLIST;

extern int              Mutual_Calculate(ESL_MSA **omsa, int *msamap, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, COVTYPE covtype, COVCLASS covclass, 
					 int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, char *gnuplot, char *dplotfile, char *R2Rfile, char *r2rversion, int r2rall, int maxFP, double expectFP, 
					 int nbpairs, double tol, int verbose, char *errbuf);
extern int              Mutual_Probs(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf);
extern int              Mutual_ValidateProbs(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateCHI     (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
						 int analyze, HITLIST **ret_hitlist, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateCHI_C16 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateCHI_C2  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateOMES    (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
						 int analyze, HITLIST **ret_hitlist, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateOMES_C16(struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateOMES_C2 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateGT      (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
						 int analyze, HITLIST **ret_hitlist, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateGT_C16  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateGT_C2   (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateMI      (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
						 int analyze, HITLIST **ret_hitlist, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMI_C16  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateMI_C2   (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateMIr     (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
						 int analyze, HITLIST **ret_hitlist, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMIr_C16 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateMIr_C2  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateMIg     (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
						 int analyze, HITLIST **ret_hitlist, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMIg_C16 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateMIg_C2  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_CalculateCOVCorrected(struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP, int nbpairs, 
						     CORRTYPE corrtype, int analyze, HITLIST **ret_hitlist, double tol, int verbose, char *errbuf);
extern int              Mutual_COVTYPEString(char **ret_covtype, COVTYPE type, char *errbuf);
extern int              Mutual_String2COVTYPE(char *covtype, COVTYPE *ret_type, char *errbuf);
extern struct mutual_s *Mutual_Create(int64_t alen, int64_t nseq, int isshuffled, int nseqthresh, ESL_ALPHABET *abc);
extern int              Mutual_ReuseCOV(struct mutual_s *mi, COVTYPE mitype, COVCLASS covclass);
extern void             Mutual_Destroy(struct mutual_s *mi);
extern int              Mutual_NaivePP(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);
extern int              Mutual_SignificantPairs_Ranking(RANKLIST **ret_ranklist, HITLIST **ret_hitlist, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, int maxFP, double expectFP,
							int nbpairs, int verbose, char *errbuf);
extern RANKLIST        *Mutual_CreateRankList();
extern int              Mutual_CreateHitList(FILE *fp, HITLIST **ret_hitlist, double threshsc, struct mutual_s *mi, int *msamap, int *ct, int N, RANKLIST *ranklist, 
					     int verbose, char *errbuf);
extern void             Mutual_FreeRankList(RANKLIST *ranklist);
extern void             Mutual_FreeHitList(HITLIST *hitlist);
extern int              Mutual_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int *ct, int verbose, char *errbuf);
extern int              Mutual_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen);
extern int              Mutual_CYKCOVCT(FILE *outfp, char *gnuplot, char *dplotfile, char *R2Rcykfile, char *R2Rversion, int R2Rall, ESL_RANDOMNESS *r, 
					ESL_MSA **msa, struct mutual_s *mi, int *msamap, int minloop, enum grammar_e G, 
					int maxFP, double expectFP, int nbpairs, char *errbuf, int verbose);
extern int              Mutual_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int *ct, struct mutual_s *mi, int *msamap, HITLIST *hitlist, int verbose, char *errbuf);
extern int              Mutual_R2R(char *r2rfile, char *r2rversion, int r2rall, ESL_MSA **msa, int *ct, int *msamap, HITLIST *hitlist, int makepdf, 
				   int verbose, char *errbuf);
extern int              Mutual_R2Rpdf(char *r2rfile, char *r2rversion, int verbose, char *errbuf);
extern int              Mutual_ExpandCT(char *r2rfile, int r2rall,  ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, int minloop, enum grammar_e G, int verbose, char *errbuf);
extern int              Mutual_ExpandCT_Naive(ESL_MSA *msa, int *ct, int minloop, int verbose, char *errbuf);
extern int              Mutual_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ct, enum grammar_e G, int minloop, int verbose, char *errbuf);
  

#endif
