/* rview_cmap - The R-veiw IPA
 *
 */
#ifndef RVIEW_CMAP_INCLUDED
#define RVIEW_CMAP_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_random.h"

#include "rview_contacts.h"
#include "rview_pdbfile.h"

typedef enum{
  MIN = 0,
  CB  = 1,
} DISTMETHOD;

int rview_ContactMap(char *pdbxfile, int *ret_nct, CLIST **ret_clist, double maxD, int minL, char *errbuf, int verbose);


#endif /*RVIEW_CMAP_INCLUDED*/
