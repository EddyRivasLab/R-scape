/* fastasp 
 *
 */
#ifndef FASTSP_INCLUDED
#define FASTSP_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"

extern int FastSP_Run(const ESL_MSA *msar, const ESL_MSA *msae, int *ret_tph, int *ret_trueh, int *ret_foundh, int *ret_cac, int *ret_ac, 
		      double *ret_sen, double *ret_ppv, double *ret_F, double *ret_TC, int lcu_r, int lcu_e, char *errbuf, int verbose);

#endif /*FASTSP_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/