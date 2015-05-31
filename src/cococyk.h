/* cococyk.h
 *
 *   
*/
#ifndef COCOCYK_INCLUDED
#define COCOCYK_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "covgrammars.h"


extern int COCOCYK    (ESL_RANDOMNESS *r, enum grammar_e  G,    ESL_SQ *sq, int *ct, int **ret_ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int COCOCYK_G6 (ESL_RANDOMNESS *r, G6param        *g6p,  ESL_SQ *sq, int *ct, int **ret_ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int COCOCYK_G6S(ESL_RANDOMNESS *r, G6Sparam       *g6sp, ESL_SQ *sq, int *ct, int **ret_ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int COCOCYK_BGR(ESL_RANDOMNESS *r, BGRparam       *bgrp, ESL_SQ *sq, int *ct, int **ret_ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int COCOCYK_G6_Fill (G6param  *p, ESL_SQ *sq, int *ct, G6_MX  *gmx, SCVAL *ret_sc, char *errbuf, int verbose);
extern int COCOCYK_G6S_Fill(G6Sparam *p, ESL_SQ *sq, int *ct, G6_MX  *gmx, SCVAL *ret_sc, char *errbuf, int verbose); 
extern int COCOCYK_BGR_Fill(BGRparam *p, ESL_SQ *sq, int *ct, BGR_MX *gmx, SCVAL *ret_sc, char *errbuf, int verbose); 
extern int COCOCYK_G6_Traceback( ESL_RANDOMNESS *r, G6param  *p, ESL_SQ *sq, int *ct, G6_MX  *gmx, int **ret_cct, char *errbuf, int verbose);
extern int COCOCYK_G6S_Traceback(ESL_RANDOMNESS *r, G6Sparam *p, ESL_SQ *sq, int *ct, G6_MX  *gmx, int **ret_cct, char *errbuf, int verbose); 
extern int COCOCYK_BGR_Traceback(ESL_RANDOMNESS *r, BGRparam *p, ESL_SQ *sq, int *ct, BGR_MX *gmx, int **ret_cct, char *errbuf, int verbose); 

#endif
