/* rfview.h
 *
 *   
*/
#ifndef RFVIEW_INCLUDED
#define RFVIEW_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"

#include "correlators.h"

extern int rfview_Depict(char *rfviewfile, char *omsafile, char *covfile, int nagg, enum agg_e *agg_method, char *helixcovfile, int makepdf, int makesvg, char *errbuf, int verbose);

#endif
