/* cococyk_bgr.h
 *
 *   
*/
#ifndef COCOCYK_BGR_INCLUDED
#define COCOCYK_BGR_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"

#include "covgrammars.h"


extern int   COCOCYK_BGR(ESL_RANDOMNESS *r, struct mutual_s *mi, int **ret_ct, SCVAL *ret_sc, int minloop, int maxFP, double target_eval, char *errbuf, int verbose);
extern int   COCOCYK_BGR_Fill(struct mutual_s *mi, GMX *cyk, SCVAL *ret_sc, int minloop, double covthresh, char *errbuf, int verbose);
extern int   COCOCYK_BGR_Traceback(ESL_RANDOMNESS *r, struct mutual_s *mi, GMX *cyk, int **ret_ct, int minloop, double covthresh, char *errbuf, int verbose);

#endif
