/* power
 *
 */
#ifndef POWER_INCLUDED
#define POWER_INCLUDED


#include <stdio.h>   

#include "easel.h"
#include "correlators.h"
#include "esl_histogram.h"

extern int   power_SPAIR_Create(SPAIR **ret_spair, int alen, int *msamap, POWER *power, CLIST *clist, int *nsubs, char *errbuf, int verbose);
extern void  power_SPAIR_Write(FILE *fp, int64_t dim, SPAIR *spair);
extern void  power_Destroy(POWER *power);
extern int   power_Read(char *powerfile, POWER **ret_power, char *errbuf, int verbose);
extern void  power_Write(FILE *fp, POWER *power, int verbose);
extern void  power_WriteFromHistograms(FILE *fp, ESL_HISTOGRAM *hsubs_bp, ESL_HISTOGRAM *hsubs_cv, int verbose);

#endif /* POWER_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
