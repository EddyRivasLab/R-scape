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

extern int    potts_ML_NLogp           (PT *pt, ESL_MSA *msa, double *ret_logp,                char *errbuf, int verbose);
extern int    potts_PLM_NLogp          (PT *pt, ESL_MSA *msa, double *ret_logp, double *dlogp, char *errbuf, int verbose);
extern int    potts_APLM_NLogp(int pos, PT *pt, ESL_MSA *msa, double *ret_logp, double *dlogp, char *errbuf, int verbose);
extern int    potts_PLM_NLogp_Packed        (int np, double *p, PT *pt, ESL_MSA *msa, double *ret_nlogp, double *dnlogp, char *errbuf, int verbose);
extern int    potts_APLM_NLogp_Packed(int i, int np, double *p, PT *pt, ESL_MSA *msa, double *ret_nlogp, double *dnlogp, char *errbuf, int verbose);
extern int    potts_CalculateCOV         (struct data_s *data);
extern int    potts_CalculateCOVFrobenius(struct data_s *data);
extern int    potts_CalculateCOVAverage  (struct data_s *data);

extern double potts_APLM_H   (int i, int a, PT *pt, ESL_DSQ *sq);
extern double potts_APLM_logz(int i,        PT *pt, ESL_DSQ *sq);

extern double potts_APLM_H_Packed(int i, int a, double *p, int L, int Kg, ESL_DSQ *sq);
extern double potts_PLM_H_Packed (int i, int a, double *p, int L, int Kg, ESL_DSQ *sq);
#endif
