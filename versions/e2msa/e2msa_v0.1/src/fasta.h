/* fasta
 *
 */
#ifndef FASTA_INCLUDED
#define FASTA_INCLUDED


#include <stdio.h>		/* FILE */

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

extern int SSEARCH_Align (char *version, const ESL_MSA *msa, ESL_MSA **ret_ssearchmsa, float *ret_sc, char *mx, int fasta_f, int fasta_g, char *errbuf, int verbose);
#endif /*FASTA_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
