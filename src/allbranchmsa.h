/* allbranchmsa - funtions to build a tree from an alignment
 *
 */
#ifndef ALLBRANCHMSA_INCLUDED
#define ALLBRANCHMSA_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_tree.h"

typedef struct chit_s {
  int64_t i;
  int64_t j;
  
  int64_t posi;
  int64_t posj;

  int     cont;
  
  double sc;
} CHIT;

typedef struct clist_s{
  int      nhit;
  CHIT   **srthit;
  CHIT    *hit;
} CLIST;


extern int       AllBranchMSA_Plot(char *plotfile, char *gnuplot, ESL_TREE *T, int *msamap, ESL_MSA *allmsa, int *ct, char *errbuf, int verbose);

#endif /*ALLBRANCHMSA_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
