/* ncbiblast
 *
 */
#ifndef NCBIBLAST_INCLUDED
#define NCBIBLAST_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_tree.h"

extern int NCBIBLAST_Align(const ESL_MSA *msa, int wordsize, ESL_MSA **ret_blastmsa, char *errbuf, int verbose);

#endif /*NCBIBLAST_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
