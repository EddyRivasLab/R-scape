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


struct mutual_s {
  int64_t      alen;
  ESL_DMATRIX *MI;     // MI  mutual information
  ESL_DMATRIX *MIa;    // MIa mutual information
  ESL_DMATRIX *MIp;    // MIp mutual information
  ESL_DMATRIX *MIr;    // MIr mutual information
  double      *H;      // entropy per position
 
};

extern int              Mutual_Calculate(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf);
extern struct mutual_s *Mutual_Create(int64_t alen);
extern void             Mutual_Destroy(struct mutual_s *mi);

#endif
