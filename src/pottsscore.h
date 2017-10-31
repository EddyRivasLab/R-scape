/* pottsscore.h
 *
 *   
 */
#ifndef POTTSSCORE_INCLUDED
#define POTTSSCORE_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

#include "correlators.h"
#include "pottsbuild.h"

#define PLMDIM(i,L,K,K2)    ( (i)*(K) + (i)*((L)-1)*(K2) - 0.5*(i)*((i)-1)*(K2) ) // \sum_{j=0}^{i-1} [ K + K*K*(L-1-j) ]
#define PLMIDX(i,j,L,K,K2)  ( PLMDIM(i,L,K,K2) + (K) + ((j)-(i)-1)*(K2) )
#define PLMIDXR(i,j,L,K,K2) ( (K) + ((j)-(i)-1)*(K2) )

#define APLMDIM(L,K,K2)     ( (K) + ((L)-1)*(K2) )   // number of parameters to optimize for a given i
#define APLMIDXL(j,K,K2)    ( (K) + (j)*(K2) )       // j < i
#define APLMIDXG(j,K,K2)    ( (K) + ((j)-1)*(K2) )   // j > i

extern int    potts_NLogp_ML                (PT *pt, ESL_MSA *msa, double *ret_logp,         char *errbuf, int verbose);
extern int    potts_NLogp_PLM               (PT *pt, ESL_MSA *msa, double *ret_logp, PT *gr, char *errbuf, int verbose);
extern int    potts_NLogp_APLM       (int i, PT *pt, ESL_MSA *msa, double *ret_logp, PT *gr, char *errbuf, int verbose);

extern double potts_Hi   (int i, int a, PT *pt, ESL_DSQ *sq);
extern double potts_Logzi(int i,        PT *pt, ESL_DSQ *sq, double *Hi);
extern double potts_Zi   (int i,        PT *pt, ESL_DSQ *sq, double *Hi);

extern int    potts_CalculateCOV         (struct data_s *data);
extern int    potts_CalculateCOVFrobenius(struct data_s *data);
extern int    potts_CalculateCOVAverage  (struct data_s *data);
#endif
