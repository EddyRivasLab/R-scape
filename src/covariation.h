/* covariation.h
 *
 *   
*/
#ifndef COVARIATION_INCLUDED
#define COVARIATION_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_histogram.h"
#include "esl_msa.h"
#include "esl_tree.h"
#include "ribosum_matrix.h"

#include "covgrammars.h"

#define BMIN -100.0     // minimum COV value (-10.0 is sensible)
#define HPTS  200       // number of point in histogram
#define W     0.5       // default histogram with

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

  RAF   = 18,
  RAFp  = 19,
  RAFa  = 20,

  RAFS  = 21,
  RAFSp = 22,
  RAFSa = 23,

  CCF   = 24,
  CCFp  = 25,
  CCFa  = 26,

  COVNONE = 27,
} COVTYPE;


typedef enum {
  C16      = 0,
  C2       = 1,
  CSELECT  = 2,
} COVCLASS;

typedef enum {
  APC = 0,
  ASC = 1
} CORRTYPE;

typedef enum{
  NAIVE      = 0,
  NULLPHYLO  = 1,
  DCA        = 2,
  AKMAEV     = 3,
} METHOD;

struct mutual_s {
  int64_t         alen;
  int64_t         nseq;
  double       ***pp;          // joint probability of two position [0,alen-1][0.alen-1][0..15]
  double        **pm;          // marginal probabilities [0,alen-1][0..3]
  double        **nseff;       // effective number of sequences  [0,alen-1][0,alen-1]
  double        **ngap;        // number of gaps  [0,alen-1][0,alen-1]

  COVTYPE         type;
  COVCLASS        class;
  ESL_DMATRIX    *COV;    // covariation matrix (MI, MIp, MIr, MIa, CHI,...)  mutual information
 
  double          besthreshCOV;
  double          minCOV;
  double          maxCOV;

  int             ishuffled;
  int             nseqthresh; // if nseq <= nseqthresh use C2 method otherwise use C16
  int             alenthresh; // if alen <= alenthresh use C2 method otherwise use C16

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
  double        *eval;
  double         scthresh;
} RANKLIST;

typedef struct hit_s {
  int64_t i;
  int64_t j;
  
  double sc;
  double Eval;

  int is_bpair;
  int is_compatible;
} HIT;

typedef struct hitlist_s{
  int     nhit;
  HIT  **srthit;
  HIT    *hit;

  double  covthresh;
}  HITLIST;

typedef enum {
  Eval    = 0, // Eval;          max expected number of CovNBPs allowed
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

struct data_s {
  FILE                *outfp;
  FILE                *outsrtfp;
  FILE                *rocfp;
  FILE                *sumfp; 
  char                *gnuplot;
  char                *dplotfile;
  char                *cykdplotfile;
  char                *R2Rfile;
  char                *R2Rcykfile;
  int                  R2Rall;
  ESL_RANDOMNESS      *r;
  RANKLIST            *ranklist_null;
  RANKLIST            *ranklist_aux;
  struct mutual_s     *mi; 
  THRESH              *thresh;
  METHOD               method;
  MODE                 mode;
  COVTYPE              covtype;
  int                  onbpairs;
  int                  nbpairs;
  ESL_TREE            *T;
  struct ribomatrix_s *ribosum;
  int                 *ct;
  int                 *msamap;
  int                  firstpos;
  double               bmin;
  double               w;
  int                  Nfit;       // minimum number of point to fit
  double               pmass;      // actual_pmass = MIN(Nfit/Nc,pmass)
  int                  doexpfit;   // TRUE for an exponential fit, default is chi-square
  double               tau;        // for chi-square fit  of null distribution, effective number of degress of freedom
  double               mu;         // for exponential fit of null distribution
  double               lambda;     // for exponential fit of null distribution
  ESL_DMATRIX         *allowpair;  // allows to propose a non-WC set pairing rules
  double               tol;
  int                  verbose;
  char                *errbuf;
  int                  donull2b;
};


extern int              cov_Calculate(struct data_s *data, ESL_MSA *msa, RANKLIST  **ret_ranklist, HITLIST **ret_hitlist, int analize);
extern int              cov_Probs(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
				   METHOD method, int donull2b, double tol, int verbose, char *errbuf);
extern int              cov_ValidateProbs(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              cov_CalculateCHI     (COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateCHI_C16 (struct mutual_s *mi,                         int verbose, char *errbuf);
extern int              cov_CalculateCHI_C2  (struct mutual_s *mi, ESL_DMATRIX *allowpair, int verbose, char *errbuf);
extern int              cov_CalculateOMES    (COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateOMES_C16(struct mutual_s *mi,                         int verbose, char *errbuf);
extern int              cov_CalculateOMES_C2 (struct mutual_s *mi, ESL_DMATRIX *allowpair, int verbose, char *errbuf);
extern int              cov_CalculateGT      (COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateGT_C16  (struct mutual_s *mi,                         int verbose, char *errbuf);
extern int              cov_CalculateGT_C2   (struct mutual_s *mi, ESL_DMATRIX *allowpair, int verbose, char *errbuf);
extern int              cov_CalculateMI      (COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateMI_C16  (struct mutual_s *mi,                         int verbose, char *errbuf);
extern int              cov_CalculateMI_C2   (struct mutual_s *mi, ESL_DMATRIX *allowpair, int verbose, char *errbuf);
extern int              cov_CalculateMIr     (COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateMIr_C16 (struct mutual_s *mi,                         int verbose, char *errbuf);
extern int              cov_CalculateMIr_C2  (struct mutual_s *mi, ESL_DMATRIX *allowpair, int verbose, char *errbuf);
extern int              cov_CalculateMIg     (COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateMIg_C16 (struct mutual_s *mi,                         int verbose, char *errbuf);
extern int              cov_CalculateMIg_C2  (struct mutual_s *mi, ESL_DMATRIX *allowpair, int verbose, char *errbuf);
extern int              cov_CalculateRAF     (COVCLASS covclass, struct data_s *data, ESL_MSA *msa, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateRAFS    (COVCLASS covclass, struct data_s *data, ESL_MSA *msa, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateCCF     (COVCLASS covclass, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_CalculateCCF_C16 (struct mutual_s *mi,                         int verbose, char *errbuf);
extern int              cov_CalculateCOVCorrected(CORRTYPE corrtype, struct data_s *data, int analyze, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern int              cov_THRESHTYPEString(char **ret_threshtype, THRESHTYPE type, char *errbuf);
extern int              cov_COVTYPEString(char **ret_covtype, COVTYPE type, char *errbuf);
extern int              cov_String2COVTYPE(char *covtype, COVTYPE *ret_type, char *errbuf);
extern struct mutual_s *cov_Create(int64_t alen, int64_t nseq, int isshuffled, int nseqthresh, int thresh, ESL_ALPHABET *abc, COVCLASS covclass);
extern int              cov_ReuseCOV(struct mutual_s *mi, COVTYPE mitype, COVCLASS covclass);
extern void             cov_Destroy(struct mutual_s *mi);
extern int              cov_NaivePP(ESL_RANDOMNESS *r, ESL_MSA *msa, struct mutual_s *mi, int donull2b, double tol, int verbose, char *errbuf);
extern int              cov_Marginals(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              cov_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					double tol, int verbose, char *errbuf);
extern int              cov_SignificantPairs_Ranking(struct data_s *data, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
extern RANKLIST        *cov_CreateRankList(double bmax, double bmin, double w);
extern int              cov_GrowRankList(RANKLIST **oranklist, double bmax, double  bmin);
extern int              cov_DumpRankList(FILE *fp, RANKLIST *ranklist);
extern int              cov_DumpHistogram(FILE *fp, ESL_HISTOGRAM *h);
extern int              cov_CreateHitList(struct data_s *data, struct mutual_s *mi, RANKLIST *ranklist, HITLIST **ret_hitlist,
					  char *covtype, char *threshtype, int usenull);
extern int              cov_WriteHitList(FILE *fp, int nhit, HITLIST *hitlist, int *msamap, int firstpos);
extern int              cov_WriteRankedHitList(FILE *fp, int nhit, HITLIST *hitlist, int *msamap, int firstpos);
extern void             cov_FreeRankList(RANKLIST *ranklist);
extern void             cov_FreeHitList(HITLIST *hitlist);
extern int              cov_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int firstpos, int *ct, int verbose, char *errbuf);
extern int              cov_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen);
extern int              cov_CYKCOVCT(struct data_s *data, ESL_MSA *msa, RANKLIST **ret_ranklist, int minloop, enum grammar_e G, double covthresh);
extern int              cov_NullFitGamma(ESL_HISTOGRAM *h, double pmass, double *ret_newmass, double *ret_mu, double *ret_lambda, double *ret_k, int verbose, char *errbuf);
extern int              cov_NullFitExponential(ESL_HISTOGRAM *h, double pmass, double *ret_newmass, double *ret_mu, double *ret_lambda, int verbose, char *errbuf);
extern int              cov_WriteHistogram(struct data_s *data, char *gnuplot, char *covhisfile, char *covqqfile, RANKLIST *ranklist, char *title);
extern int              cov_PlotHistogramSurvival(struct data_s *data, char *gnuplot, char *covhisfile, RANKLIST *ranklist, char *title, int dosvg);
extern int              cov_PlotHistogramQQ(struct data_s *data, char *gnuplot, char *covqqfile, RANKLIST *ranklist, char *title, int dosvg);
extern int              cov_PlotNullCov(char *gnuplot, char *nullcovfile, double maxBP, double maxcovBP, double maxcovRBPf, int dosvg);
extern int              cov_DotPlot(char *gnuplot, char *dplotfile,  ESL_MSA *msa, int *ct, struct mutual_s *mi, int *msamap, int firstpos, HITLIST *hitlist, 
				    int dosvg, int verbose, char *errbuf);
extern int              cov_R2R(char *r2rfile, int r2rall, ESL_MSA *msa, int *ct, HITLIST *hitlist,
				int makepdf, int makesvg, int verbose, char *errbuf);
extern int              cov_R2Rpdf(char *r2rfile, int verbose, char *errbuf);
extern int              cov_R2Rsvg(char *r2rfile, int verbose, char *errbuf);
extern int              cov_ExpandCT(char *r2rfile, int r2rall,  ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, int minloop, enum grammar_e G, 
				     int verbose, char *errbuf);
extern int              cov_ExpandCT_Naive(ESL_MSA *msa, int *ct, int minloop, int verbose, char *errbuf);
extern int              cov_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ct, enum grammar_e G, int minloop, int verbose, char *errbuf);
extern int              cov_ranklist_Bin2Bin(int b, ESL_HISTOGRAM *h, ESL_HISTOGRAM *new, int *ret_newb);
#endif
