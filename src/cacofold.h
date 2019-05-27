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

extern int CACO_CYK          (ESL_RANDOMNESS *r, enum grammar_e  G,      ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_DECODING     (ESL_RANDOMNESS *r, enum grammar_e  G,      ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6X_CYK      (ESL_RANDOMNESS *r, G6Xparam        *g6xp,  ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6X_DECODING (ESL_RANDOMNESS *r, G6Xparam        *g6xp,  ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6XS_CYK     (ESL_RANDOMNESS *r, G6XSparam       *g6xsp, ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6XS_DECODING(ESL_RANDOMNESS *r, G6XSparam       *g6xsp, ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_RBG_CYK      (ESL_RANDOMNESS *r, RBGparam        *rbgp,  ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_RBG_DECODING (ESL_RANDOMNESS *r, RBGparam        *rbgp,  ESL_SQ *sq, SPAIR *spair, int *ct, int **ret_ct, SCVAL *ret_sc, COVLIST *exclude, char *errbuf, int verbose);
extern int CACO_G6X_GetParam (G6Xparam  **ret_p, char *errbuf, int verbose);
extern int CACO_G6XS_GetParam(G6XSparam **ret_p, char *errbuf, int verbose);
extern int CACO_RBG_GetParam( RBGparam  **ret_p, char *errbuf, int verbose);
extern int CACO_G6X_Fill_CYK (G6Xparam  *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6X_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6XS_Fill_CYK(G6XSparam *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6X_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_RBG_Fill_CYK( RBGparam  *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, RBG_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_G6X_Traceback_CYK (ESL_RANDOMNESS *r, G6Xparam  *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6X_MX *cyk, int **ret_cct, char *errbuf, int verbose);
extern int CACO_G6XS_Traceback_CYK(ESL_RANDOMNESS *r, G6XSparam *p, ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, G6X_MX *cyk, int **ret_cct, char *errbuf, int verbose); 
extern int CACO_RBG_Traceback_CYK (ESL_RANDOMNESS *r, RBGparam *p,  ESL_SQ *sq, SPAIR *spair, int *ct, COVLIST *exclude, RBG_MX *cyk, int **ret_cct, char *errbuf, int verbose); 

#endif
