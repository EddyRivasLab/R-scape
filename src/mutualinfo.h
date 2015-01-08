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


struct mutual_s {
  int64_t      alen;
  ESL_DMATRIX *pp;     // joint probability dist. of two position at the root
  double      *ps;     // probability dist. of a position at the root

  ESL_DMATRIX *MI;     // MI  mutual information
  ESL_DMATRIX *MIa;    // MIa mutual information
  ESL_DMATRIX *MIp;    // MIp mutual information
  ESL_DMATRIX *MIr;    // MIr mutual information
  double      *H;      // entropy per position
 
};

extern int              Mutual_Calculate(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf);
extern struct mutual_s *Mutual_Create(int64_t alen);
extern void             Mutual_Destroy(struct mutual_s *mi);
extern int              Mutual_PostOrderPP(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf);
extern int              Mutual_PostOrderPS(ESL_MSA *msa, ESL_TREE *T, struct ribomatrix_s *ribosum, struct mutual_s *mi, int verbose, char *errbuf);

#endif
