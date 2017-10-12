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
  ML    = 1,
  PLM   = 2,
  APLM  = 3,
  GINV  = 4,
  ACE   = 5,
  BML   = 6,
} PTTRAIN;

typedef enum {
  SCNONE = 0,
  FROEB  = 1,
  AVG    = 2,
  DI     = 3,
} PTSCTYPE;

typedef struct potts_s {
  int64_t   L;       /* length of alignment */
  double  **h;       /* parameter h[0..L-1][0..K-1]         */
  double ***e;       /* couplings e[0..L-1][0..L-1][a*K+b]  */

  PTTRAIN       train;
  PTSCTYPE      sctype;
  double        muh;      /* regularization constant for hi */
  double        mue;      /* regularization constant for eij */
  
  ESL_ALPHABET *abc;
} PT;


struct optimize_data {
  ESL_RANDOMNESS *r;
  PT             *pt;
  ESL_MSA        *msa;
  int             pos;
  double          tol;
  char           *errbuf;
  int             verbose;
};

extern int   potts_AssignZero(ESL_RANDOMNESS *r, PT *pt, char *errbuf, int verbose);
extern int   potts_AssignGaussian(ESL_RANDOMNESS *r, PT *pt, double mu, double sigma, char *errbuf, int verbose);
extern int   potts_AssignGT(ESL_RANDOMNESS *r, ESL_MSA *msa, PT *pt, float tol, char *errbuf, int verbose);
extern PT   *potts_Build(ESL_RANDOMNESS *r, ESL_MSA *msa, double ptmuh, double ptmue, PTTRAIN pttrain, PTSCTYPE ptsctype, FILE *pottsfp,
			 float tol, char *errbuf, int verbose);
extern PT   *potts_Create(int64_t L, int K, ESL_ALPHABET *abc, double muh, double mue, PTTRAIN pttrain, PTSCTYPE ptsctype);
extern void  potts_Destroy(PT *pt);
extern int   potts_GaugeZeroSum(PT *pt, char *errbuf, int verbose);
extern int   potts_OptimizeCGD_PLM (PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
extern int   potts_OptimizeCGD_APLM(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
extern int   potts_OptimizeLBFGS_APLM(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
extern PT   *potts_Read(char *paramfile, ESL_ALPHABET *abc, char *errbuf);
extern void  potts_Write(FILE *fp, PT *pt);

#endif
