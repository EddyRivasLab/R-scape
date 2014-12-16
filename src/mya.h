/* mya
 *
 */
#ifndef MYA_INCLUDED
#define MYA_INCLUDED


#include <stdio.h>	

#include "easel.h"

#include "fetchpfamdb.h"

extern int mya_BetweenSpeciesFromFile(char *distfile, char *tx1, char *tx2, float *ret_mya, int verbose);
extern int mya_BetweenSpeciesFromSSIFile(char *ssidistfile, char **taxalist,  int ntaxa, char *tx1, char *tx2, float *ret_mya, int verbose);
extern int mya_BetweenSpecies(char *sp1, char *sp2, float *ret_mya, char *errbuf, int verbose);
extern int mya_BetweenLevels(int n, int m, SPECIES *SPE, float *ret_mya, char *errbuf, int verbose);
#endif /* MYA_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
