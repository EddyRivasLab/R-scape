/* cykcov.h
 *
 *   
*/
#ifndef CYKCOV_INCLUDED
#define CYKCOV_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"

#include "covariation.h"

typedef float SCVAL;

typedef struct {
  int       L;    /* msa length */
  SCVAL  **dp;   /* L * L triangular DP matrix */
} GMX;

extern int   CYKCOV(ESL_RANDOMNESS *r, struct mutual_s *mi, int **ret_ct, SCVAL *ret_sc, int minloop, char *errbuf, int verbose);
extern int   CYKCOV_Fill(struct mutual_s *mi, GMX **ret_cyk, SCVAL *ret_sc, int minloop, char *errbuf, int verbose);
extern int   CYKCOV_Traceback(ESL_RANDOMNESS *r, struct mutual_s *mi, GMX *cyk, int **ret_ct, int minloop, char *errbuf, int verbose);
extern GMX  *GMX_Create(int L);
extern void  GMX_Destroy(GMX *gmx);
extern void  GMX_Dump(FILE *fp, GMX *gmx);


#endif
