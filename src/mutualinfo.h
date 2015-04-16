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

typedef enum {
  CHI  = 0,
  MI   = 1,
  MIa  = 2,
  MIp  = 3,
  MIr  = 4,
} MITYPE;

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

  ESL_DMATRIX    *COV;    // covariation matrix (MI, MIp, MIr, MIa, CHI,...)  mutual information
  double         *H;      // entropy per position
 
  double          besthreshCOV;
  double          minCOV;
  double          maxCOV;

  ESL_ALPHABET   *abc;
};


extern int              Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					 METHOD method, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf);
extern int              Mutual_Probs(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, METHOD method, double tol, int verbose, char *errbuf);
extern int              Mutual_ValidateProbs(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateH(struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateCHI(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMI (struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMIa(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMIp(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf);
extern int              Mutual_CalculateMIr(struct mutual_s *mi, int *ct, int plotroc, int maxFP, double tol, int verbose, char *errbuf);
extern struct mutual_s *Mutual_Create(int64_t alen, int64_t nseq, ESL_ALPHABET *abc);
extern int              Mutual_ReuseCOV(struct mutual_s *mi);
extern void             Mutual_Destroy(struct mutual_s *mi);
extern int              Mutual_NaivePP(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_NaivePS(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPS(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);
extern int              Mutual_SignificantPairs_Ranking(struct mutual_s *mi, int *ct, MITYPE whichmi, int plotroc, int maxFP, int verbose, char *errbuf);
extern int              Mutual_SignificantPairs_ZScore(struct mutual_s *mi, int *ct, int verbose, char *errbuf);


#endif
