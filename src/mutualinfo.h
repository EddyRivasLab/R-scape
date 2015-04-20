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
  GT    = 1,
  OMES  = 2,
  MI    = 3,
  MIr   = 4,

  CHIp  = 5,
  GTp   = 6,
  OMESp = 7,
  MIp   = 8,
  MIrp  = 9,

  CHIa  = 10,
  GTa   = 11,
  OMESa = 12,
  MIa   = 13,
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
  double       ***cp;    // joint  counts for two position [0,alen-1][0.alen-1][0..15]
  double       ***pp;    // joint  probabilities for two position [0,alen-1][0.alen-1][0..15]
  double        **ps;    // single probabilities for  a  position [0,alen-1][0..3]

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
extern int              Mutual_NaivePS(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPS(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);
extern int              Mutual_SignificantPairs_Ranking(struct mutual_s *mi, int *ct, FILE *rocfp, int maxFP, int ishuffled, int verbose, char *errbuf);
extern int              Mutual_SignificantPairs_ZScore(struct mutual_s *mi, int *ct, int verbose, char *errbuf);


#endif
