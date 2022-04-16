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

#include "e2_profilesq.h"
#include "covgrammars.h"
#include "correlators.h"
#include "structure.h"

#define TOLVAL  0.7

extern int CACO_CYK          (ESL_RANDOMNESS *r, enum grammar_e G, FOLDPARAM *foldparam, PSQ  *psq,  struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc,  char *errbuf, int verbose);
extern int CACO_DECODING     (ESL_RANDOMNESS *r, enum grammar_e G, FOLDPARAM *foldparam, PSQ  *psq,  struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_MEA          (ESL_RANDOMNESS *r,                   FOLDPARAM *foldparam, POST *post,                      SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose);

extern int CACO_G6X_CYK      (ESL_RANDOMNESS *r,                   FOLDPARAM *foldparam, G6Xparam *g6xp,   PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6X_DECODING (ESL_RANDOMNESS *r,                   FOLDPARAM *foldparam, G6Xparam *g6xp,   PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, POST *post,             char *errbuf, int verbose);
extern int CACO_G6XS_CYK     (ESL_RANDOMNESS *r,                   FOLDPARAM *foldparam, G6XSparam *g6xsp, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6XS_DECODING(ESL_RANDOMNESS *r,                   FOLDPARAM *foldparam, G6XSparam *g6xsp, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, POST *post,             char *errbuf, int verbose);
extern int CACO_RBG_CYK      (ESL_RANDOMNESS *r,                   FOLDPARAM *foldparam, RBGparam *rbgp,   PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, int *ct, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_RBG_DECODING (ESL_RANDOMNESS *r,                   FOLDPARAM *foldparam, RBGparam *rbgp,   PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, POST *post,             char *errbuf, int verbose);

extern int CACO_G6X_GetParam (G6Xparam  **ret_p, char *errbuf, int verbose);
extern int CACO_G6XS_GetParam(G6XSparam **ret_p, char *errbuf, int verbose);
extern int CACO_RBG_GetParam (RBGparam  **ret_p, char *errbuf, int verbose);

extern int CACO_G6X_Fill_CYK (FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk,              SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6X_Inside   (FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx,              SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_G6X_Outside  (FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6X_Posterior(FOLDPARAM *foldparam, G6Xparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, POST *post,    char *errbuf, int verbose);

extern int CACO_G6XS_Fill_CYK (FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6XS_Inside   (FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx,              SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_G6XS_Outside  (FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *omx, G6X_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_G6XS_Posterior(FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *imx, G6X_MX *omx, POST *post,    char *errbuf, int verbose);

extern int CACO_RBG_Fill_CYK  (FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *cyk, SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_RBG_Inside    (FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx,              SCVAL *ret_sc, char *errbuf, int verbose); 
extern int CACO_RBG_Outside   (FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *omx, RBG_MX *imx, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_RBG_Posterior (FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *imx, RBG_MX *omx, POST *post,    char *errbuf, int verbose);

extern int CACO_G6X_Traceback_CYK (ESL_RANDOMNESS *r, FOLDPARAM *foldparam, G6Xparam  *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk, int *ct, char *errbuf, int verbose);
extern int CACO_G6XS_Traceback_CYK(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, G6XSparam *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *cyk, int *ct, char *errbuf, int verbose); 
extern int CACO_RBG_Traceback_CYK (ESL_RANDOMNESS *r, FOLDPARAM *foldparam, RBGparam  *p, PSQ *psq, struct mutual_s *mi, SPAIR *spair, int *covct, COVLIST *exclude, RBG_MX *cyk, int *ct, char *errbuf, int verbose); 

extern int CACO_G6X_MEA_GetParam(G6Xparam **ret_p, double gamma, char *errbuf, int verbose);
extern int CACO_MEA_Fill_CYK                        (FOLDPARAM *foldparam, G6Xparam *meap, POST *post, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *gmx, SCVAL *ret_sc, char *errbuf, int verbose);
extern int CACO_MEA_Traceback_CYK(ESL_RANDOMNESS *r, FOLDPARAM *foldparam, G6Xparam *meap, POST *post, SPAIR *spair, int *covct, COVLIST *exclude, G6X_MX *gmx, int *ct,       char *errbuf, int verbose);

#endif
