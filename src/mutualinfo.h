/* mutualinfo.h
 *
 *   
*/
#ifndef MUTUALINFO_INCLUDED
#define MUTUALINFO_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_tree.h"
#include "ribosum_matrix.h"

typedef enum {
  MI  = 0,
  MIa = 1,
  MIp = 2,
  MIr = 3,
} MItype;

struct mutual_s {
  int64_t      alen;
  double       ***pp;    // joint probability dist. of two position at the root [0,alen-1][0.alen-1][0..15]
  double        **ps;    // probability dist. of a position at the root [0,alen-1][0..3]

  ESL_DMATRIX   *MI;     // MI  mutual information
  ESL_DMATRIX   *MIa;    // MIa mutual information
  ESL_DMATRIX   *MIp;    // MIp mutual information
  ESL_DMATRIX   *MIr;    // MIr mutual information
  double        *H;      // entropy per position

  double         threshMI;
  double         threshMIa;
  double         threshMIp;
  double         threshMIr;

  double         minMI;
  double         maxMI;
  double         minMIa;
  double         maxMIa;
  double         minMIp;
  double         maxMIp;
  double         minMIr;
  double         maxMIr;
};

extern int              Mutual_Analyze(int *ct, struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_AnalyzeSignificantPairs(int *ct, struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_AnalyzeRanking(int *ct, struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					 int naive, double tol, int verbose, char *errbuf);
extern struct mutual_s *Mutual_Create(int64_t alen, int K);
extern void             Mutual_Destroy(struct mutual_s *mi);
extern int              Mutual_NaivePP(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_NaivePS(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);
extern int              Mutual_PostOrderPS(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, 
					   double tol, int verbose, char *errbuf);

#endif
