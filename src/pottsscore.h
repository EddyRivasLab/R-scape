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

extern int potts_SumH             (PT *pt, ESL_MSA *msa, double *ret_sc,   char *errbuf, int verbose);
extern int potts_MLLogp           (PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose);
extern int potts_PLMLogp          (PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose);
extern int potts_APLMLogp(int pos, PT *pt, ESL_MSA *msa, double *ret_logp, char *errbuf, int verbose);
extern int potts_CalculateCOV(struct data_s *data);
extern int potts_CalculateCOVFrobenius(struct data_s *data);
extern int potts_CalculateCOVAverage(struct data_s *data);
#endif
