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

extern int   potts_Score(PT *pt, ESL_MSA *msa, float tol, char *errbuf, int verbose);
#endif
