/* homology
 *
 */
#ifndef HOMOLOGY_INCLUDED
#define HOMOLOGY_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"


extern int homol_RunPHMMER(int n, int m, const ESL_MSA *msa, float *ret_sc, float *ret_eval, char *errbuf, int verbose);
extern int homol_ParsePHMMERtblout(char *tblout, float *ret_sc, float *ret_eval, char *errbuf, int verbose);

#endif /* HOMOLOGY_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
