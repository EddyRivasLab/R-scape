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

extern int potts_UnNormLogp(PT *pt, ESL_MSA *msa, double *ret_sc, char *errbuf, int verbose);
extern int potts_CalculateCOV(PT *pt, struct data_s *data, RANKLIST **ret_ranklist, HITLIST **ret_hitlist);
#endif
