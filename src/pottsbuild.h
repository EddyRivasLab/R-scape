/* pottsbuild.h
 *
 *   
*/
#ifndef POTTSBUILD_INCLUDED
#define POTTSBUILD_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

typedef enum {
  NONE = 0,
  FULL = 1,
  APLM = 2,
  GINV = 3,
  ACE  = 4,
  BML  = 5,
} POTTSTRAIN;

typedef struct potts_s {
  int64_t   L;       /* length of alignment */
  double  **h;       /* parameter h[0..L-1][0..K-1]         */
  double ***e;       /* couplings e[0..L-1][0..L-1][a*K+b]  */

  POTTSTRAIN    train;
  double        mu;      /* regularization constant */
  
  ESL_ALPHABET *abc;
} PT;

struct optimize_data {
  PT            *pt;
  ESL_MSA       *msa;
  float          logp;
  double         firststep;
  double         tol;
  char          *errbuf;
  int            verbose;
};

extern PT   *potts_Create(int64_t L, ESL_ALPHABET *abc, double mu, POTTSTRAIN pottstrain);
extern void  potts_Destroy(PT *pt);
extern int   potts_Assign(PT *pt, double val);
extern int   potts_Build(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
extern int   potts_OptimizeGD(PT *pt, ESL_MSA *msa, float firststep, float tol, char *errbuf, int verbose);
extern int   potts_GaugeZeroSum(PT *pt, double tol, char *errbuf);
#endif
