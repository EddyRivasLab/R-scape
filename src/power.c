/*  power - calculate prediction power based on number of substitutiosn
 * 
 *
 * ER, Sun Jul 29 16:19:56 CEST 2018 [Atlantic ocean]
 */

#include "p7_config.h"

#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>
	
#include "easel.h"
#include "esl_buffer.h"
#include "esl_mem.h"

#include "power.h"

void
power_Destroy(POWER *power)
{
  if (power == NULL) return;
  if (power->subs) free(power->subs);
  if (power->prob) free(power->prob);
}

int
power_Read(char *powerfile, POWER **ret_power, char *errbuf, int verbose)
{
  ESL_BUFFER      *bf    = NULL;
  char            *subs  = NULL;
  char            *prob  = NULL;
  POWER           *power = NULL;
  char            *p;
  esl_pos_t        salloc;
  esl_pos_t        len;
  esl_pos_t        n;
  esl_pos_t        i;
  int              idx;
  int              status;
  
  status = esl_buffer_Open(powerfile, NULL, &bf);
  if      (status == eslENOTFOUND) esl_fatal("open failed: %s",     bf ? bf->errmsg : "(no other diagnostics available)");
  else if (status == eslFAIL)      esl_fatal("gzip -dc failed: %s", bf ? bf->errmsg : "(no other diagnostics available)");
  else if (status != eslOK)        esl_fatal("open failed with error code %d", status);

  ESL_ALLOC(power, sizeof(POWER));
  power->ns   = 0;
  power->subs = NULL;
  power->prob = NULL;
  
  while (( status = esl_buffer_GetLine(bf, &p, &n)) == eslOK) 
    {
      if (power->ns == 0) {
	ESL_ALLOC(power->subs, sizeof(double)*(power->ns+1));
	ESL_ALLOC(power->prob, sizeof(double)*(power->ns+1));
      }
      else {
	ESL_REALLOC(power->subs, sizeof(double)*(power->ns+1));
	ESL_REALLOC(power->prob, sizeof(double)*(power->ns+1));
      }

      idx = 0;
      len = 0;
      if (power->ns == 0) {
	salloc = n + 1;
	ESL_ALLOC(subs, sizeof(char) * salloc);
	ESL_ALLOC(prob, sizeof(char) * salloc);
      }

      for (i = 0; i < n; i++) {
	if (! isspace(p[i])) {
	  if      (idx == 0) subs[len++] = p[i];
	  else if (idx == 1) prob[len++] = p[i];
	}
	else { len = 0; idx ++; }
      }
      
      power->subs[power->ns] = atof(subs);
      power->prob[power->ns] = atof(prob);
      if (power->prob[power->ns] < 0.) power->prob[power->ns] = 0.;
      if (power->prob[power->ns] > 1.) power->prob[power->ns] = 1.;
      
      power->ns ++;
    }

  if (status != eslEOF) esl_fatal("file %s: expected EOF, got code %d", bf->filename, status);
  esl_buffer_Close(bf);

  // we got it. Print if asked
  if (verbose) power_Write(stdout, power, verbose);

  *ret_power = power;
  return eslOK;

 ERROR:
  if (power) power_Destroy(power);
  return status;
}

void
power_Write(FILE *fp, POWER *power, int verbose)
{
  int n;

  for (n = 0; n < power->ns; n ++)
    fprintf(fp, "%f\t%f\n", power->subs[n], power->prob[n]);
}

void
power_WriteFromHistograms(FILE *fp, ESL_HISTOGRAM *hsubs_bp, ESL_HISTOGRAM *hsubs_cv, int verbose)
{
  int            i;
  double         x;
  double         fsubs;

  for (i = hsubs_bp->imin; i <= hsubs_bp->imax; i++)
    {
      x = esl_histogram_Bin2LBound(hsubs_bp,i);
      
      if (hsubs_bp->obs[i] > 0) {
	fsubs = (double) hsubs_cv->obs[i] / (double) hsubs_bp->obs[i];	
	if (fsubs > 1.) esl_fatal("error writing power function");
	if (verbose) printf("%f\t%f\t%lld\t%lld\n", x, fsubs, hsubs_bp->obs[i], hsubs_cv->obs[i]);
	fprintf(fp, "%f\t%f\t%lld\t%lld\n", x, fsubs, hsubs_bp->obs[i], hsubs_cv->obs[i]);
      } 
    }
}
