/* covariation.h
 *
 *   
*/
#ifndef COVARIATION_INCLUDED
#define COVARIATION_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_histogram.h"
#include "esl_msa.h"
#include "esl_tree.h"
#include "ribosum_matrix.h"

#include "covgrammars.h"

#define W     0.1     // COV with

#define NCOVTYPE = 5;

typedef enum {
  CHI   = 0,
  CHIp  = 1,
  CHIa  = 2,
  CHIs  = 3,

  GT    = 4,
  GTp   = 5,
  GTa   = 6,
  GTs   = 7,

  MI    = 8,
  MIp   = 9,
  MIa   = 10,
  MIs   = 11,

  MIr   = 12,
  MIrp  = 13,
  MIra  = 14,
  MIrs  = 15,

  MIg   = 16,
  MIgp  = 17,
  MIga  = 18,
  MIgs  = 19,

  OMES  = 20,
  OMESp = 21,
  OMESa = 22,
  OMESs = 23,

  COVALL  = 24,
  COVNONE = 25,
} COVTYPE;


typedef enum {
  C16      = 0,
  C2       = 1,
  CSELECT  = 2,
} COVCLASS;

typedef enum {
  APC = 0,
  ASC = 1,
  SCA = 2, //spearman correction of attenuation
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

typedef enum {
  NullNONE = 0,
  Null1    = 1,
  Null1b   = 2,
  Null2    = 3,
  Null2b   = 4,
  Null3    = 5,
  Null4    = 6,
} NULLTYPE;


typedef struct ranklist_s {
  ESL_HISTOGRAM *ha;             /* histogram of scores (all pairs) */
  ESL_HISTOGRAM *ht;             /* histogram of scores (truncated pairs == no ss pairs) */
  double        *covBP;
  double        *covNBP;
  double         scthresh;
} RANKLIST;

typedef struct hit_s {
  int64_t i;
  int64_t j;
  
  double sc;
  double Eval;

  double covNBP;
  double covNBPu;
  double covNBPf;

  double covRBP;
  double covRBPu;
  double covRBPf;

  int is_bpair;
  int is_compatible;
} HIT;

typedef struct hitlist_s{
  int     nhit;
  HIT    *hit;

  double  covthresh;
}  HITLIST;

typedef enum {
  covNBP  = 0, // cov NonBPs:    max total        (covNPB)      allowed
  covNBPu = 1, // cov NonBPs;    max per_position (covNBP/alen) allowed
  covNBPf = 2, // cov NonBPs;    max fraction     (covNBP/NBP)  allowed
  covRBP  = 3, // cov RandomBPs: max total        (covRPB)      allowed
  covRBPu = 4, // cov RandomBPs; max per_position (covRBP/alen) allowed
  covRBPf = 5, // cov RandomBPs; max fraction     (covRBP/RBP)  allowed
  Eval    = 6, // Eval;          max expected number of CovNBPs allowed
} THRESHTYPE;

typedef struct thresh_s {
  THRESHTYPE type;
  double     val;  // the actual thershold value
  double     sc;   // cov for the given alignment at threshold
} THRESH;

typedef enum {
  GIVSS = 0,
  CYKSS = 1,
  RANSS = 2,
} MODE;

extern int              cov_Calculate(ESL_RANDOMNESS *r, ESL_MSA **omsa, int *msamap, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				      RANKLIST *ranklist_null, RANKLIST *ranklist_aux, RANKLIST  **ret_ranklist, HITLIST **ret_hitlist, METHOD method, 
				      COVTYPE covtype, COVCLASS covclass, 
				      int *ct, double bmin, double w, double pmass, double *ret_mu, double *ret_lambda, 
				      FILE *outfp, FILE *rocfp, FILE *sumfp, char *gnuplot, char *dplotfile, 
				      char *R2Rfile, char *r2rversion, int r2rall, THRESH *thresh, MODE mode, int nbpairs, int donull2b, 
				      double tol, int verbose, char *errbuf);
extern int              cov_Probs(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, int donull2b, 
				  double tol, int verbose, char *errbuf);
extern int              cov_ValidateProbs(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              cov_CalculateCHI     (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, double bmin, 
					      double w, double pmass, double *ret_mu, double *ret_lambda, 
					      FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
					      int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, 
					      double tol, int verbose, char *errbuf);
extern int              cov_CalculateCHI_C16 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateCHI_C2  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateOMES    (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, double bmin, double w, 
					      double pmass, double *ret_mu, double *ret_lambda, 
					      FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
					      int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, RANKLIST *ranklist_aux,
					      double tol, int verbose, char *errbuf);
extern int              cov_CalculateOMES_C16(struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateOMES_C2 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateGT      (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, double bmin, double w, 
					      double pmass, double *ret_mu, double *ret_lambda, 
					      FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
					      int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, RANKLIST *ranklist_aux,
					      double tol, int verbose, char *errbuf);
extern int              cov_CalculateGT_C16  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateGT_C2   (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateMI      (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, double bmin, double w, 
					      double pmass, double *ret_mu, double *ret_lambda, 
					      FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
					      int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, RANKLIST *ranklist_aux,
					      double tol, int verbose, char *errbuf);
extern int              cov_CalculateMI_C16  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateMI_C2   (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateMIr     (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, double bmin, double w, 
					      double pmass, double *ret_mu, double *ret_lambda, 
					      FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
					      int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, RANKLIST *ranklist_aux,
					      double tol, int verbose, char *errbuf);
extern int              cov_CalculateMIr_C16 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateMIr_C2  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateMIg     (COVCLASS covclass, struct mutual_s *mi, int *msamap, int *ct, double bmin, double w,
					      double pmass, double *ret_mu, double *ret_lambda, 
					      FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
					      int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, 
					      double tol, int verbose, char *errbuf);
extern int              cov_CalculateMIg_C16 (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateMIg_C2  (struct mutual_s *mi, int verbose, char *errbuf);
extern int              cov_CalculateCOVCorrected(struct mutual_s *mi, int *msamap, int *ct, double bmin, double w, 
						  double pmass, double *ret_mu, double *ret_lambda, 
						  FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode, int nbpairs, 
						  CORRTYPE corrtype, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, 
						  RANKLIST *ranklist_null, RANKLIST *ranklist_aux,
						  double tol, int verbose, char *errbuf);
extern int              cov_THRESHTYPEString(char **ret_threshtype, THRESHTYPE type, char *errbuf);
extern int              cov_COVTYPEString(char **ret_covtype, COVTYPE type, char *errbuf);
extern int              cov_String2COVTYPE(char *covtype, COVTYPE *ret_type, char *errbuf);
extern struct mutual_s *cov_Create(int64_t alen, int64_t nseq, int isshuffled, int nseqthresh, ESL_ALPHABET *abc);
extern int              cov_ReuseCOV(struct mutual_s *mi, COVTYPE mitype, COVCLASS covclass);
extern void             cov_Destroy(struct mutual_s *mi);
extern int              cov_NaivePP(ESL_RANDOMNESS *r, ESL_MSA *msa, struct mutual_s *mi, int donull2b, double tol, int verbose, char *errbuf);
extern int              cov_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					double tol, int verbose, char *errbuf);
extern int              cov_SignificantPairs_Ranking(RANKLIST *ranklist_null, RANKLIST *ranklist_aux, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, 
						     struct mutual_s *mi, int *msamap, int *ct, double bmin, double w, 
						     double pmass, double *ret_mu, double *ret_lambda,
						     FILE *outfp, FILE *rocfp, FILE *sumfp, THRESH *thresh, MODE mode,
						     int nbpairs, int verbose, char *errbuf);
extern RANKLIST        *cov_CreateRankList(double bmax, double bmin, double w);
extern int              cov_GrowRankList(RANKLIST **oranklist, double bmax, double  bmin);
extern int              cov_DumpRankList(FILE *fp, RANKLIST *ranklist);
extern int              cov_DumpHistogram(FILE *fp, ESL_HISTOGRAM *h);
extern int              cov_CreateHitList(FILE *fp, HITLIST **ret_hitlist, THRESH *thresh, struct mutual_s *mi, int *msamap, int *ct, 
					  RANKLIST *ranklist, RANKLIST *ranklist_null, double pmass, double mu, double lambda, 
					  int usenull, char *covtype, char *threshtype, MODE mode, int verbose, char *errbuf);
extern void             cov_FreeRankList(RANKLIST *ranklist);
extern void             cov_FreeHitList(HITLIST *hitlist);
extern int              cov_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int *ct, int verbose, char *errbuf);
extern int              cov_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen);
extern int              cov_CYKCOVCT(FILE *outfp, char *gnuplot, char *dplotfile, char *R2Rcykfile, char *R2Rversion, int R2Rall, ESL_RANDOMNESS *r, 
				     ESL_MSA **msa, struct mutual_s *mi, int *msamap, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, 
				     RANKLIST **ret_ranklist, double bmin, double w, double pmass, double *ret_mu, double *ret_lambda, 
				     int minloop, enum grammar_e G, THRESH *thresh,  double covthresh, int nbpairs, char *errbuf, int verbose);
extern int              cov_ExpFitHistogram(ESL_HISTOGRAM *h, double pmass, double *ret_newmass, double *ret_mu, double *ret_lambda, int verbose, char *errbuf);
extern int              cov_WriteHistogram(char *gnuplot, char *covhisfile, char *nullcovhisfile, RANKLIST *ranklist, RANKLIST *ranklist_null, 
					   RANKLIST *ranklist_aux, char *title, double pmass, double mu, double lambda, int verbose, char *errbuf);
extern int              cov_PlotHistogramSurvival(char *gnuplot, char *covhisfile, RANKLIST *ranklist, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, 
						  char *title, double pmass, double mu, double lambda, int dosvg, char *errbuf);
extern int              cov_CreateNullCov(char *gnuplot, char *nullcovfile, int L, int *ct, RANKLIST *ranklist, RANKLIST *ranklist_null, int dosvg, char *errbuf);
extern int              cov_PlotNullCov(char *gnuplot, char *nullcovfile, double maxBP, double maxcovBP, double maxcovRBPf, int dosvg);
extern int              cov_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int *ct, struct mutual_s *mi, int *msamap, HITLIST *hitlist, 
				    int dosvg, int verbose, char *errbuf);
extern int              cov_R2R(char *r2rfile, char *r2rversion, int r2rall, ESL_MSA **msa, int *ct, int *msamap, HITLIST *hitlist, int makepdf, int makesvg,
				int verbose, char *errbuf);
extern int              cov_R2Rpdf(char *r2rfile, char *r2rversion, int verbose, char *errbuf);
extern int              cov_R2Rsvg(char *r2rfile, char *r2rversion, int verbose, char *errbuf);
extern int              cov_ExpandCT(char *r2rfile, int r2rall,  ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, int minloop, enum grammar_e G, int verbose, char *errbuf);
extern int              cov_ExpandCT_Naive(ESL_MSA *msa, int *ct, int minloop, int verbose, char *errbuf);
extern int              cov_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ct, enum grammar_e G, int minloop, int verbose, char *errbuf);
extern int              cov_ranklist_Bin2Bin(int b, ESL_HISTOGRAM *h, ESL_HISTOGRAM *new, int *ret_newb);
#endif
