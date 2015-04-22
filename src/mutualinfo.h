/* mutualinfo.h
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

#define NCOVTYPE = 5;
typedef enum {
  CHI   = 0,
  CHIp  = 1,
  CHIa  = 2,

  GT    = 3,
  GTp   = 4,
  GTa   = 5,

  OMES  = 6,
  OMESp = 7,
  OMESa = 8,

  MI    = 9,
  MIp   = 10,
  MIa   = 11,

  MIr   = 12,
  MIrp  = 13,
  MIra  = 14,

  COVNONE  = 15,
} COVTYPE;

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

  COVTYPE         type;
  ESL_DMATRIX    *COV;    // covariation matrix (MI, MIp, MIr, MIa, CHI,...)  mutual information
  double         *H;      // entropy per position
 
  double          besthreshCOV;
  double          minCOV;
  double          maxCOV;

  ESL_ALPHABET   *abc;
};


extern int              Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					 METHOD method, int *ct, FILE *rocfp, int maxFP, int ishuffled, double tol, int verbose, char *errbuf);
extern int              Mutual_Probs(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf);
extern int              Mutual_ValidateProbs(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateH(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateCHI (struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, int analyze, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateOMES(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, int analyze, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateGT  (struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, int analyze, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMI  (struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, int analyze, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMIr (struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, int analyze, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateCOVCorrected(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, CORRTYPE corrtype, double tol, int verbose, char *errbuf);
extern int              Mutual_COVTYPEString(char **ret_covtype, COVTYPE type, char *errbuf);
extern int              Mutual_String2COVTYPE(char *covtype, COVTYPE *ret_type, char *errbuf);
extern struct mutual_s *Mutual_Create(int64_t alen, int64_t nseq, ESL_ALPHABET *abc);
extern int              Mutual_ReuseCOV(struct mutual_s *mi, COVTYPE mitype);
extern void             Mutual_Destroy(struct mutual_s *mi);
extern int              Mutual_NaivePP(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);
extern int              Mutual_SignificantPairs_Ranking(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, int verbose, char *errbuf);
extern int              Mutual_SignificantPairs_ZScore(struct mutual_s *mi, int *ct, int verbose, char *errbuf);


#endif
