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

#include "plot.h"
#include "power.h"

POWERHIS *
power_Histogram_Create(int bmin, int bmax, double w)
{
  POWERHIS *powerhis = NULL;
  int       status;

  ESL_ALLOC(powerhis, sizeof(POWERHIS));
  
  powerhis->hsubs_pr = esl_histogram_Create(bmin, bmax, w);
  powerhis->hsubs_ur = esl_histogram_Create(bmin, bmax, w);
  powerhis->hsubs_bp = esl_histogram_Create(bmin, bmax, w);
  powerhis->hsubs_cv = esl_histogram_Create(bmin, bmax, w);
  
  return powerhis;

 ERROR:
  return NULL;
}

void
power_Histogram_Destroy(POWERHIS *powerhis)
{
  if (!powerhis) return;
  
  if (powerhis->hsubs_pr) free(powerhis->hsubs_pr);
  if (powerhis->hsubs_ur) free(powerhis->hsubs_ur);
  if (powerhis->hsubs_bp) free(powerhis->hsubs_bp);
  if (powerhis->hsubs_cv) free(powerhis->hsubs_cv);
  
  free(powerhis);
}

int 
power_SPAIR_Create(int *ret_np, SPAIR **ret_spair, int alen, int *msamap, POWER *power, CLIST *clist, int *nsubs, int *ndouble, char *errbuf, int verbose)
{
  SPAIR   *spair = NULL;
  double   prob;
  int64_t  dim = alen * (alen - 1.) / 2;
  int64_t  subs;
  int64_t  n = 0;
  int64_t  s;
  int      i, j;
  int      ipos;
  int      c;
  int      status;
  
  if (ret_spair == NULL) return eslOK;

  if (!nsubs && !ndouble) esl_fatal("need some kind of substituions!");
  
  ESL_ALLOC(spair, sizeof(SPAIR) * dim);
  for (i = 0; i < alen-1; i ++) 
    for (j = i+1; j < alen; j ++) {
      
      subs            = (nsubs)? nsubs[i] + nsubs[j] : ndouble[i*alen+j];                    
      spair[n].i      = msamap[i]+1;
      spair[n].j      = msamap[j]+1;
      spair[n].nsubs  = subs;
      spair[n].power  = 0.;
      spair[n].bptype = BPNONE;
      
      if (power) {
	prob = 0.;
	for (s = 0; s < power->ns; s ++) {
	  if (subs > power->subs[s]) prob = power->prob[s];
	  else break;
	}
	spair[n].power = prob;
      }
       
      if (clist) {
	for (c = 0; c < clist->ncnt; c++) {
	  if (spair[n].i == clist->cnt[c].posi && spair[n].j == clist->cnt[c].posj) {
	    spair[n].bptype = clist->cnt[c].bptype;
	    break;
	  }
	}
      }
      
      if (spair[n].bptype == WWc) 
	if (verbose) printf("WWc: %lld-%lld nsubs %lld prob %f\n", spair[n].i, spair[n].j, spair[n].nsubs, spair[n].power);
      
      n ++;
    }

  if (ret_np) *ret_np = n;
  *ret_spair = spair;
  return eslOK;
  
 ERROR:
  if (spair) free(spair);
  return status;
}
 
void
power_SPAIR_Write(FILE *fp, int64_t dim, SPAIR *spair)
{
  double     expect    = 0.;
  double     avgsub    = 0;
  int64_t    nbp       = 0;
  int64_t    n;

  fprintf(fp, "# left_pos      right_pos    substitutions      power\n");
  fprintf(fp, "#--------------------------------------------------------\n");
  for (n = 0; n < dim; n ++)
    if (spair[n].bptype == WWc) {
      nbp ++;
      expect    += spair[n].power;
      avgsub    += spair[n].nsubs;
      fprintf(fp, "# %lld\t\t%lld\t\t%lld\t\t%.2f\n", spair[n].i, spair[n].j, spair[n].nsubs, spair[n].power);
    }
  avgsub /= (nbp > 0)? nbp : 1;
  
  fprintf(fp, "#\n# BPAIRS %lld\n", nbp);
  fprintf(fp, "# avg substitutions per BP  %.1f\n", avgsub);
  fprintf(fp, "# BPAIRS expected to covary %.1f\n", expect);
}

void
power_Destroy(POWER *power)
{
  if (power == NULL) return;
  if (power->subs) free(power->subs);
  if (power->prob) free(power->prob);
}

int
power_Read(char *powerfile, int doublesubs, POWER **ret_power, char *errbuf, int verbose)
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
  int              reached = FALSE;
  int              idx;
  int              status;
  
  status = esl_buffer_Open(powerfile, NULL, &bf);
  if      (status == eslENOTFOUND) esl_fatal("open failed: %s",     bf ? bf->errmsg : "(no other diagnostics available)");
  else if (status == eslFAIL)      esl_fatal("gzip -dc failed: %s", bf ? bf->errmsg : "(no other diagnostics available)");
  else if (status != eslOK)        esl_fatal("open failed with error code %d", status);

  ESL_ALLOC(power, sizeof(POWER));
  power->ns   = 0;
  power->type = (doublesubs)? DOUB : SUBS;
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
	else { subs[len] = '\0'; len = 0; idx ++; }
      }
      prob[len] = '\0';
      
      power->subs[power->ns] = atof(subs);
      power->prob[power->ns] = atof(prob);
      if (power->prob[power->ns] < 0.)   power->prob[power->ns] = 0.;
      if (power->prob[power->ns] > 1.) { power->prob[power->ns] = 1.; reached = TRUE; }
      if (reached)                       power->prob[power->ns] = 1.;
      
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
power_WriteFromHistograms(FILE *fp, POWERHIS *powerhis, int verbose)
{
  int            b, bcv;
  int64_t        bp, cv;
  double         x, y;
  double         fsubs;

  for (b = powerhis->hsubs_bp->imin; b <= powerhis->hsubs_bp->imax; b++)
    {
      bp = powerhis->hsubs_bp->obs[b];
      x  = esl_histogram_Bin2LBound(powerhis->hsubs_bp, b);

      if (esl_histogram_Bin2LBound(powerhis->hsubs_cv, powerhis->hsubs_cv->imin) > x) cv = 0;
      else {     
	esl_histogram_Score2Bin(powerhis->hsubs_cv, x, &bcv); bcv ++;
	
	y = esl_histogram_Bin2LBound(powerhis->hsubs_cv, bcv);      
	if (x != y) esl_fatal("power_WriteFromHistogram() error, b %d x %f bcv %d y %f", b, x, bcv, y);
	
	cv = powerhis->hsubs_cv->obs[bcv];
	if (bp < cv) esl_fatal("power_WriteFromHistogram() error, x = %f bp %lld cv %lld", x, bp, cv);
      }
      fsubs = (bp > 0)? (double)cv/(double)bp : 0.;
      
      if (verbose) printf( "%f\t%f\t%lld\t%lld\n", x, fsubs, bp, cv);
      fprintf(fp,          "%f\t%f\t%lld\t%lld\n", x, fsubs, bp, cv);
    }
}

void
power_PlotHistograms(char *gnuplot, char *powerhisfile, FILE *powerhisfp, POWERHIS *powerhis, char *powerfile, int powerdouble, char *errbuf, int verbose)
{
  char *powersubspdf   = NULL;
  char *subspdf        = NULL;
  char *subshisfile[2];
  int   status;
  
  esl_histogram_Plot(powerhisfp, powerhis->hsubs_pr);
  esl_histogram_Plot(powerhisfp, powerhis->hsubs_ur);
  fclose(powerhisfp);
  
  esl_sprintf(&subshisfile[0], "%s.pair",   powerhisfile);
  esl_sprintf(&subshisfile[1], "%s.unpair", powerhisfile);
  plot_write_Histogram(subshisfile[0], powerhis->hsubs_pr);
  plot_write_Histogram(subshisfile[1], powerhis->hsubs_ur);
    
  esl_sprintf(&subspdf, "%s.pdf", powerhisfile);
  if (powerdouble) plot_gplot_Histogram(gnuplot, subspdf, 2, subshisfile, "number of double substitutions", FALSE, errbuf, verbose);
  else             plot_gplot_Histogram(gnuplot, subspdf, 2, subshisfile, "number of substitutions",        FALSE, errbuf, verbose);
  free(subspdf);
  
  esl_sprintf(&powersubspdf,  "%s.pdf", powerfile);
  if (powerhis->hsubs_bp->n > 0) {
    if (powerdouble) 
      status = plot_gplot_XYfile(gnuplot, powersubspdf,   powerfile, 1, 2, "number of double substitutions in basepairs", "fraction of covarying basepairs", errbuf);
    else
      status = plot_gplot_XYfile(gnuplot, powersubspdf,   powerfile, 1, 2, "number of substitutions in basepairs",        "fraction of covarying basepairs", errbuf);
    if (status != eslOK) printf("%s\n", errbuf);
  }
  free(powersubspdf);
}
