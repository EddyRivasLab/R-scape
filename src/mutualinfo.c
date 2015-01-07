/* mutualinfo.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "mutualinfo.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"

int                 
Mutual_Calculate(ESL_MSA *msa, struct mutual_s *mi, double tol, int verbose, char *errbuf)
{
  return eslOK;
}

struct mutual_s *
Mutual_Create(int64_t alen)
{
  struct mutual_s *mi = NULL;
  int              status;

  ESL_ALLOC(mi, sizeof(struct mutual_s));
  mi->alen = alen;
  
  mi->MI  = esl_dmatrix_Create(alen, alen);
  mi->MIa = esl_dmatrix_Create(alen, alen);
  mi->MIp = esl_dmatrix_Create(alen, alen);
  mi->MIr = esl_dmatrix_Create(alen, alen);
  ESL_ALLOC(mi->H, sizeof(double) * alen);
 
  return mi;

 ERROR:
  return NULL;
}

void                
Mutual_Destroy(struct mutual_s *mi)
{
  if (mi) {
    esl_dmatrix_Destroy(mi->MI);
    esl_dmatrix_Destroy(mi->MIa);
    esl_dmatrix_Destroy(mi->MIp);
    esl_dmatrix_Destroy(mi->MIr);
    free(mi->H);
    free(mi);
  }
}
