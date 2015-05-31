/* cococyk.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covgrammars.h"
#include "cococyk.h"


int
COCOCYK(ESL_RANDOMNESS *r, enum grammar_e G, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6param  *g6p  = NULL;
  G6Sparam *g6sp = NULL;
  BGRparam *bgrp = NULL;

  int status;

  /* get the grammar parameters and run the corresponding CYK */
  switch(G) {
  case G6:
    status = COCOCYK_G6(r, g6p, sq, ct, ret_cct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case G6S:
    status = COCOCYK_G6S(r, g6sp, sq, ct, ret_cct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  case BGR:
   status = COCOCYK_BGR(r, bgrp, sq, ct, ret_cct, ret_sc, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    break;
  }

  if (g6p)  free(g6p);
  if (g6sp) free(g6sp);
  if (bgrp) free(bgrp);
  return eslOK;

 ERROR:
  return status;
}


int
COCOCYK_G6(ESL_RANDOMNESS *r, G6param  *p, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6_MX *gmx = NULL;
  int    status;

  gmx = G6MX_Create(sq->n);

  /* Fill the cyk matrix */
  if ((status = COCOCYK_G6_Fill        (p, sq, ct, gmx, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = COCOCYK_G6_Traceback(r, p, sq, ct, gmx, ret_cct, errbuf, verbose))  != eslOK) goto ERROR;
  
  G6MX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) G6MX_Destroy(gmx);
  return status;
}

int
COCOCYK_G6S(ESL_RANDOMNESS *r, G6Sparam  *p, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  G6_MX *gmx = NULL;
  int    status;

  gmx = G6MX_Create(sq->n);

  /* Fill the cyk matrix */
  if ((status = COCOCYK_G6S_Fill        (p, sq, ct, gmx, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = COCOCYK_G6S_Traceback(r, p, sq, ct, gmx, ret_cct, errbuf, verbose))  != eslOK) goto ERROR;
  
  G6MX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) G6MX_Destroy(gmx);
  return status;
}

int
COCOCYK_BGR(ESL_RANDOMNESS *r, BGRparam  *p, ESL_SQ *sq, int *ct, int **ret_cct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  BGR_MX *gmx = NULL;
  int     status;

  gmx = BGRMX_Create(sq->n);

  /* Fill the cyk matrix */
  if ((status = COCOCYK_BGR_Fill        (p, sq, ct, gmx, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;    
  /* Report a traceback */
  if ((status = COCOCYK_BGR_Traceback(r, p, sq, ct, gmx, ret_cct, errbuf, verbose))  != eslOK) goto ERROR;
  
  BGRMX_Destroy(gmx);
  return eslOK;

 ERROR:
  if (gmx) BGRMX_Destroy(gmx);
  return status;

}


int
COCOCYK_G6_Fill(G6param  *p, ESL_SQ *sq, int *ct, G6_MX *gmx, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  return eslOK;
}

int
COCOCYK_G6S_Fill(G6Sparam  *p, ESL_SQ *sq, int *ct, G6_MX *gmx, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  return eslOK;
}

int
COCOCYK_BGR_Fill(BGRparam  *p, ESL_SQ *sq, int *ct, BGR_MX *gmx, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  return eslOK;
}

int
COCOCYK_G6_Traceback(ESL_RANDOMNESS *r, G6param  *p, ESL_SQ *sq, int *ct, G6_MX *gmx, int **ret_cct, char *errbuf, int verbose) 
{
  return eslOK;

}

int
COCOCYK_G6S_Traceback(ESL_RANDOMNESS *r, G6Sparam  *p, ESL_SQ *sq, int *ct, G6_MX *gmx, int **ret_cct, char *errbuf, int verbose) 
{
  return eslOK;

}

int
COCOCYK_BGR_Traceback(ESL_RANDOMNESS *r, BGRparam  *p, ESL_SQ *sq, int *ct, BGR_MX *gmx, int **ret_cct, char *errbuf, int verbose) 
{
  return eslOK;

}
