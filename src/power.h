/* power
 *
 */
#ifndef POWER_INCLUDED
#define POWER_INCLUDED


#include <stdio.h>   

#include "easel.h"
#include "esl_histogram.h"


typedef struct power_s {
  int64_t  ns;
  double  *subs;
  double  *prob;
} POWER;

extern void power_Destroy(POWER *power);
extern int  power_Read(char *powerfile, POWER **ret_power, char *errbuf, int verbose);
extern void power_Write(FILE *fp, POWER *power, int verbose);
extern void power_WriteFromHistograms(FILE *fp, ESL_HISTOGRAM *hsubs_bp, ESL_HISTOGRAM *hsubs_cv, int verbose);

#endif /* POWER_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
