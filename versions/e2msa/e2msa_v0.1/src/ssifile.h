/* create a ssi file
 * 
 * ER, Sat Nov 10 16:07:39 EST 2012 [Janelia]
 * SVN $Id: $
 * SVN $URL: $
 */
#ifndef SSIFILE_INCLUDED
#define SSIFILE_INCLUDED

#include "easel.h"
#include "esl_ssi.h"

#include "fetchpfamdb.h"

extern int ssi_spfile_Create(char *spfile, char **ret_ssispfile);
extern int ssi_spfile_GetSpeciesInfo(char *ssispfile, char *sqname, int n, SPECIES *SPE);
extern int ssi_distfile_Create(char *distfile, char **ret_ssidistfile, char ***ret_taxalist, int *ret_ntaxa);
extern int ssi_distfile_GetDistance(char *ssidistfile, char **taxalist, int ntaxa, char *tx1, char *tx2, float *ret_mya);
extern int ssi_riotbl_Create(char *riotbl, char **ret_ssiriotbl, char ***ret_taxalist, int *ret_ntaxa);
extern int ssi_riotbl_GetOrthsc(char *ssiriotbl, char **taxalist, int ntaxa, char *tx1, char *tx2, float *ret_orthsc, int verbose);

#endif /*SSIFILE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
