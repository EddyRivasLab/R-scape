/* MSAProbs
 *
 */
#ifndef MSAPROBS_INCLUDED
#define MSAPROBS_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_tree.h"

extern int MSAProbs_Align(const ESL_MSA *msa, ESL_MSA **ret_musclemsa, char *errbuf, int verbose);

#endif /*MSAPROBS_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
