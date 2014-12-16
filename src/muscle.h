/* muscle 
 *
 */
#ifndef MUSCLE_INCLUDED
#define MUSCLE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_tree.h"

extern int MUSCLE_Align(const ESL_MSA *msa, ESL_MSA **ret_musclemsa, char *errbuf, int verbose);

#endif /*MUSCLE_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
