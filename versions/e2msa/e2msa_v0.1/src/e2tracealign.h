/* tracealign
 *
 */
#ifndef TRACEALIGN_INCLUDED
#define TRACEALIGN_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_tree.h"

#include "e2.h"
#include "e2_profilesq.h"
#include "e2_msa.h"

extern int e2_Tracealign(char *msaname, ESL_TREE *T, PSQ **sq, E2_TRACE **tr, ESL_MSA **ret_msa, E2_ALI e2ali, char *errbuf, int verbose);
extern int e2_TracealignSeqs(E2_TRACE *tr, PSQ *sql, PSQ *sqr, int a_idx, int dl_idx, int dr_idx, ESL_MSA *msa, char *errbuf, int verbose);

#endif /*TREACEALIGN_INCLUDED*/


/************************************************************
 * @LICENSE@
 ************************************************************/
