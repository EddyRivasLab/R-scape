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
#include "esl_random.h"

typedef enum {
  NONE = 0,
  FULL = 1,
  APLM = 2,
  GINV = 3,
  ACE  = 4,
  BML  = 5,
} PTTRAIN;

typedef enum {
  FROEB = 0,
  AVG   = 1,
  DI    = 2,
} PTSCTYPE;

typedef struct potts_s {
  int64_t   L;       /* length of alignment */
  double  **h;       /* parameter h[0..L-1][0..K-1]         */
  double ***e;       /* couplings e[0..L-1][0..L-1][a*K+b]  */

  PTTRAIN       train;
  PTSCTYPE      sctype;
  double        mu;      /* regularization constant */
  
  ESL_ALPHABET *abc;
} PT;


struct optimize_data {
  ESL_RANDOMNESS *r;
  PT             *pt;
  ESL_MSA        *msa;
  float           logp;
  double          firststep;
  double          tol;
  char           *errbuf;
  int             verbose;
};

extern PT   *potts_Build(ESL_RANDOMNESS *r, ESL_MSA *msa, double ptmu, PTTRAIN pttrain, PTSCTYPE ptsctype, float tol, char *errbuf, int verbose);
extern int   potts_OptimizeGD(PT *pt, ESL_MSA *msa, float firststep, float tol, char *errbuf, int verbose);
extern PT   *potts_Create(int64_t L, ESL_ALPHABET *abc, double mu, PTTRAIN pttrain, PTSCTYPE ptsctype);
extern void  potts_Destroy(PT *pt);
extern void  potts_Dump(FILE *fp, PT *pt);
extern int   potts_InitGaussian(ESL_RANDOMNESS *r, PT *pt, double mu, double sigma);
extern int   potts_InitGT(ESL_RANDOMNESS *r, ESL_MSA *msa, PT *pt, float tol, char *errbuf, int verbose);
extern int   potts_GaugeZeroSum(PT *pt, char *errbuf, int verbose);
#endif
