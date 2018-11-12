/* cykcov.h
 *
 *   
*/
#ifndef CYKCOV_INCLUDED
#define CYKCOV_INCLUDED

#include <stdio.h>		/* FILE */


#include "easel.h"

#include "covariation.h"
#include "covgrammars.h"


extern int   CYKCOV(ESL_RANDOMNESS *r, struct mutual_s *mi, CLIST *clist, int **ret_ct, char **ret_ss,
		    SCVAL *ret_sc, int minloop, THRESH *thresh,  char *errbuf, int verbose);
extern int   CYKCOV_Fill(struct mutual_s *mi, CLIST *clist, GMX *cyk, SCVAL *ret_sc, int minloop, THRESH *thresh, char *errbuf, int verbose);
extern int   CYKCOV_Traceback(ESL_RANDOMNESS *r, struct mutual_s *mi, CLIST *clist, GMX *cyk, int **ret_ct, int minloop, THRESH *thresh, char *errbuf, int verbose);

#endif
