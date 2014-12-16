/* e2_evolve
 *
 */
#ifndef E2EVOLVE_INCLUDED
#define E2EVOLVE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_tree.h"
#include "esl_msa.h"

#include "e2.h"
#include "e1_bg.h"
#include "e1_rate.h"

#define LENNAME 20

extern int e2_evolve_Alignment(ESL_RANDOMNESS *r, E1_RATE *R, ESL_TREE *T, E1_BG *bg, ESL_MSA **ret_msa, double tol, char *errbuf, int verbose);

#endif /*E2EVOLVE_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
