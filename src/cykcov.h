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
#include "structure.h"


extern int   CYKCOV(ESL_RANDOMNESS *r, struct mutual_s *mi, CLIST *clist, int *ret_nct, int ***ret_ctlist, COVLIST ***ret_covexclude,
		    int ncvpairs, int minloop, THRESH *thresh,  char *errbuf, int verbose);
extern int   CYKCOV_Structures(ESL_RANDOMNESS *rng, struct mutual_s *mi, CLIST *clist, int *ret_nct, int ***ret_ctlist, COVLIST ***ret_covexclude,
			       int ncvpairs, int minloop, THRESH *thresh, char *errbuf, int verbose) ;
extern int   CYKCOV_Both(ESL_RANDOMNESS *rng, struct mutual_s *mi, CLIST *clist, COVLIST *explained, SCVAL *ret_sc, int **ret_ct, int *ret_ncv,
			 int minloop, THRESH *thresh, char *errbuf, int verbose);
extern int   CYKCOV_Fill(struct mutual_s *mi, CLIST *clist, COVLIST *explained, GMX *cyk, SCVAL *ret_sc, int minloop,
			 THRESH *thresh, char *errbuf, int verbose);
extern int   CYKCOV_Traceback(ESL_RANDOMNESS *r, struct mutual_s *mi, CLIST *clist, COVLIST *explained, GMX *cyk, int **ret_ct, int minloop,
			      THRESH *thresh, char *errbuf, int verbose);

#endif
