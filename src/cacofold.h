/* cacofold.h
 *
 *   
*/
#ifndef CACOFOLD_INCLUDED
#define CACOFOLD_INCLUDED

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "covgrammars.h"
#include "correlators.h"
#include "structure.h"

extern int CACO_CYK         (ESL_RANDOMNESS *r, enum grammar_e  G,    ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_DECODING    (ESL_RANDOMNESS *r, enum grammar_e  G,    ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6_CYK      (ESL_RANDOMNESS *r, G6param        *g6p,  ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6_DECODING (ESL_RANDOMNESS *r, G6param        *g6p,  ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6S_CYK     (ESL_RANDOMNESS *r, G6Sparam       *g6sp, ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6S_DECODING(ESL_RANDOMNESS *r, G6Sparam       *g6sp, ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_RBG_CYK(ESL_RANDOMNESS *r, RBGparam       *rbgp, ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_RBG_DECODING(ESL_RANDOMNESS *r, RBGparam       *rbgp, ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6_GetParam (G6param  **ret_p, char *errbuf, int verbose);
extern int CACO_G6S_GetParam(G6Sparam **ret_p, char *errbuf, int verbose);
extern int CACO_RBG_GetParam(RBGparam **ret_p, char *errbuf, int verbose);
extern int CACO_G6_Fill_CYK (G6param  *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6_MX  *cyk, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6S_Fill_CYK(G6Sparam *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6_MX  *cyk, SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_RBG_Fill_CYK(RBGparam *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, RBG_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_G6_Traceback_CYK( ESL_RANDOMNESS *r, G6param  *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6_MX  *cyk, int **ret_cct, char *errbuf, int verbose);
extern int CACO_G6S_Traceback_CYK(ESL_RANDOMNESS *r, G6Sparam *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6_MX  *cyk, int **ret_cct, char *errbuf, int verbose); 
extern int CACO_RBG_Traceback_CYK(ESL_RANDOMNESS *r, RBGparam *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, RBG_MX *cyk, int **ret_cct, char *errbuf, int verbose); 

#endif
