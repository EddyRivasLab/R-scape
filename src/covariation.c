/* covariation.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "cococyk.h"
#include "covariation.h"
#include "covgrammars.h"
#include "covmeasures.h"
#include "cykcov.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"

static int cykcov_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop);


int
COV_SignificantPairs_Ranking(RANKLIST **ret_ranklist, HITLIST **ret_hitlist, struct mutual_s *mi, int *msamap, int *ct, FILE *outfp, FILE *rocfp, FILE *sumfp, 
				int maxFP, double expectFP, int nbpairs, int verbose, char *errbuf)
{
  ESL_DMATRIX *mtx = mi->COV;
  char        *covtype = NULL;
  RANKLIST    *ranklist = NULL;
  HITLIST     *hitlist = NULL;
  double       pval;
  double       sen;
  double       ppv;
  double       F;
  double       expect;
  double       bestF = 0.0;
  double       bestsen;
  double       bestppv;
  double       thresh;
  double       besthresh = mi->besthreshCOV;
  double       maxFPF;
  double       maxFPsen;
  double       maxFPppv;
  double       maxFPthresh;
  double       expectFPF;
  double       expectFPsen;
  double       expectFPppv;
  double       expectFPthresh;
  double       oneFPF;
  double       oneFPsen;
  double       oneFPppv;
  double       oneFPthresh;
  double       expectTF_frac_total;   // fraction of covarying basepairs relative to the total number of basepairs
  double       expectTF_frac_surv;    // fraction of covarying basepairs relative to the basepairs that survive the gapthresh
  int          fp, tf, t, f, neg;
  int          oneFP_tf, oneFP_t, oneFP_f, oneFP_fp;
  int          best_tf, best_t, best_f, best_fp;
  int          maxFP_tf, maxFP_t, maxFP_f, maxFP_fp;
  int          expectFP_tf, expectFP_t, expectFP_f, expectFP_fp;
  int          i, j;
  int          x;
  int          status;

  COV_COVTYPEString(&covtype, mi->type, errbuf);

  if (rocfp) {
    fprintf(rocfp, "\n# %s ", covtype);  
    if (mi->ishuffled) fprintf(rocfp, "shuffled thresh fp tf found true negatives sen ppv F\n"); 
    else               fprintf(rocfp, "thresh fp tf found true negatives sen ppv F\n"); 
  }

  ranklist = COV_CreateRankList(mi->alen, BMAX, BMIN, W);
  ranklist->scmin = ESL_MAX(ranklist->bmin,mi->minCOV);
  ranklist->scmax = mi->maxCOV;
  if (ranklist->bmax < ranklist->scmax) ESL_XFAIL(eslFAIL, errbuf, "bmax < scmax");

  for (thresh = ranklist->scmax; thresh > ranklist->scmin-ranklist->w; thresh -= ranklist->w) {

    f = t = tf = 0;
    for (i = 0; i < mi->alen-1; i ++) 
      for (j = i+1; j < mi->alen; j ++) {
	if (mtx->mx[i][j] > thresh)   f  ++;
	if (ct[i+1] == j+1) {         t  ++;
	  if (mtx->mx[i][j] > thresh) tf ++;
	}
      }
    
    fp     = f - tf;
    sen    = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv    = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F      = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;
    expect = (mi->alen > 0)? (double)fp/(double)mi->alen : 0.0;
    
    neg = mi->alen * (mi->alen-1) / 2 - t;
    if (rocfp) fprintf(rocfp, "%.5f %d %d %d %d %d %.2f %.2f %.2f\n", thresh, fp, tf, f, t, neg, sen, ppv, F);
     
    x = (int)((thresh-ranklist->bmin)/ranklist->w);
    ranklist->covBP[x]  = (double)tf;
    ranklist->covNBP[x] = (double)(f - tf);
    
    if (expect <= expectFP) {
      expectFPF      = F;
      expectFPsen    = sen;
      expectFPppv    = ppv;
      expectFP_tf    = tf;
      expectFP_f     = f;
      expectFP_t     = t;
      expectFP_fp    = fp;
      expectFPthresh = thresh;
     }
    
    if (maxFP >= 0) {
      if (fp < 1) {
	oneFPF      = F;
	oneFPsen    = sen;
	oneFPppv    = ppv;
	oneFP_tf    = tf;
	oneFP_f     = f;
	oneFP_t     = t;
	oneFP_fp    = fp;
	oneFPthresh = thresh;
      }
      if (fp <= maxFP) {
	maxFPF      = F;
	maxFPsen    = sen;
	maxFPppv    = ppv;
	maxFP_tf    = tf;
	maxFP_f     = f;
	maxFP_t     = t;
	maxFP_fp    = fp;
	maxFPthresh = thresh;
      }
      if (F > bestF) { 
	bestF     = F; 
	bestsen   = sen;
	bestppv   = ppv;
	best_tf   = tf;
	best_t    = t;
	best_f    = f;
	best_fp   = fp;
	besthresh = thresh; 
      }  
    }
  }
      
  if (maxFP >= 0) {
    if (best_fp < maxFP_fp) {
      if (best_fp == 0 && oneFP_tf > best_tf) {
	if (outfp) fprintf(outfp, "# %s before 1FP %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, oneFPthresh, ranklist->scmin, ranklist->scmax,
			   oneFP_fp, oneFP_tf, oneFP_t, oneFP_f, oneFPsen, oneFPppv, oneFPF);
	status = COV_CreateHitList(outfp, &hitlist, oneFPthresh, mi, msamap, ct, ranklist, verbose, errbuf);
      }
      else {
	if (outfp) fprintf(outfp, "# %s optimalF %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, besthresh, ranklist->scmin, ranklist->scmax,
			   best_fp, best_tf, best_t, best_f, bestsen, bestppv, bestF);
	status = COV_CreateHitList(outfp, &hitlist, besthresh, mi, msamap, ct, ranklist, verbose, errbuf);
      } 
    }
    
    if (outfp) fprintf(outfp, "# %s maxFP=%d %f [%f,%f] [%d | %d %d %d | %f %f %f] \n", covtype, maxFP, maxFPthresh, ranklist->scmin, ranklist->scmax,
		       maxFP_fp, maxFP_tf, maxFP_t, maxFP_f, maxFPsen, maxFPppv, maxFPF);
    status = COV_CreateHitList(outfp, &hitlist, maxFPthresh, mi, msamap, ct, ranklist, verbose, errbuf);
  }
  else {
    expectTF_frac_total = (nbpairs    > 0)? 100.*(double)expectFP_tf/(double)nbpairs    : 0.0;
    expectTF_frac_surv  = (expectFP_t > 0)? 100.*(double)expectFP_tf/(double)expectFP_t : 0.0;
    if (sumfp) fprintf(sumfp, "%s\t%d\t%d\t%d\t%.2f\t%.2f\t", covtype, expectFP_tf, expectFP_t, nbpairs, expectTF_frac_surv, expectTF_frac_total);
    
    status = COV_FisherExactTest(&pval, expectFP_tf, expectFP_fp, expectFP_t, mi->alen);
    
    if (outfp) fprintf(outfp, "# %s expectFP=%f (%d cov nonBPs) pval = %f | cov_BP %d/%d covariations %d | sen %f ppv %f F %f] \n", covtype, expectFP, expectFP_fp, pval,
		       expectFP_tf, expectFP_t, expectFP_f, expectFPsen, expectFPppv, expectFPF);
    status = COV_CreateHitList(outfp, &hitlist, expectFPthresh, mi, msamap, ct, ranklist, verbose, errbuf); 
  }
  if (status != eslOK) goto ERROR;

  if (ret_ranklist) *ret_ranklist = ranklist; else COV_FreeRankList(ranklist);
  if (ret_hitlist)  *ret_hitlist  = hitlist;  else COV_FreeHitList(hitlist);
 
  if (covtype) free(covtype); 
  return eslOK;

 ERROR:
  if (ranklist)  COV_FreeRankList(ranklist);
  if (hitlist)   COV_FreeHitList(hitlist);
  if (covtype)   free(covtype); 
  return status;
}

RANKLIST *
COV_CreateRankList(int L, double bmax, double bmin, double w)
{
  RANKLIST *ranklist = NULL;
  int       status;
  
  ESL_ALLOC(ranklist, sizeof(RANKLIST));

  ranklist->bmax = bmax * (double)L;
  ranklist->bmin = bmin * (double)L;
  ranklist->w    = w;
  ranklist->nb   = (int)((ranklist->bmax-ranklist->bmin) / w);

  ranklist->scmin =  DBL_MAX;
  ranklist->scmax = -DBL_MIN;
  
  ESL_ALLOC(ranklist->covBP,  sizeof(double) * ranklist->nb);
  ESL_ALLOC(ranklist->covNBP, sizeof(double) * ranklist->nb);
 
  esl_vec_DSet(ranklist->covBP,  ranklist->nb, 0.);
  esl_vec_DSet(ranklist->covNBP, ranklist->nb, 0.);

  return ranklist;

 ERROR:
  return NULL;
}

int 
COV_CreateHitList(FILE *outfp, HITLIST **ret_hitlist, double threshsc, struct mutual_s *mi, int *msamap, int *ct, RANKLIST *ranklist, int verbose, char *errbuf)
{
  HITLIST *hitlist = NULL;
  int      alloc_nhit = 5;
  int      nhit;
  int      h = 0;
  int      i, j;
  int      ih, jh;
  int      x;
  int      status;
  
  ESL_ALLOC(hitlist, sizeof(HITLIST));

  nhit = alloc_nhit;
  ESL_ALLOC(hitlist->hit, sizeof(HIT) * nhit);

 for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      if (mi->COV->mx[i][j] > threshsc) {

	if (h == nhit - 1) {
	  nhit += alloc_nhit;
	  ESL_REALLOC(hitlist->hit, sizeof(HIT) * nhit);
	}
	for (x = 0; x < ranklist->nb; x ++) {
	  if (mi->COV->mx[i][j] <= ranklist->bmin+(double)x*ranklist->w) hitlist->hit[h].expcovNBP = (mi->alen > 0)? ranklist->covNBP[x]/(double)mi->alen:0.;
	  else break;
	}
	/* initialize */
	 hitlist->hit[h].is_bpair      = FALSE;
	 hitlist->hit[h].is_compatible = FALSE;

	hitlist->hit[h].i = i;
	hitlist->hit[h].j = j;
	hitlist->hit[h].sc = mi->COV->mx[i][j];
	if (ct[i+1] == j+1) { hitlist->hit[h].is_bpair = TRUE;  }
	else                { 
	  hitlist->hit[h].is_bpair = FALSE; 
	  if (ct[i+1] == 0 && ct[j+1] == 0) hitlist->hit[h].is_compatible = TRUE;
	} 
 
	h ++;
      }
    }
 nhit = h;
 hitlist->nhit = nhit;

 //sort by sc? should go here
  
  for (h = 0; h < nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    if (outfp) {
      if      (hitlist->hit[h].is_bpair)      { fprintf(outfp, "* %d %d\t%.4f\t%.2f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc, hitlist->hit[h].expcovNBP); }
      else if (hitlist->hit[h].is_compatible) { fprintf(outfp, "~ %d %d\t%.4f\t%.2f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc, hitlist->hit[h].expcovNBP); }
      else                                    { fprintf(outfp, "  %d %d\t%.4f\t%.2f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc, hitlist->hit[h].expcovNBP); } 
    }
  }

  if (ret_hitlist) *ret_hitlist = hitlist; else COV_FreeHitList(hitlist);
  return eslOK;

 ERROR:
  if (hitlist) COV_FreeHitList(hitlist);
  return status;
}

void
COV_FreeRankList(RANKLIST *ranklist)
{
  if (ranklist == NULL) return;

  if (ranklist->covBP)  free(ranklist->covBP);
  if (ranklist->covNBP) free(ranklist->covNBP);

  free(ranklist);
}

void
COV_FreeHitList(HITLIST *hitlist)
{
  free(hitlist->hit);
  free(hitlist);
}

int
COV_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int *ct, int verbose, char *errbuf)
{
  double       avgi, avgj;
  double       stdi, stdj;
  double       zscorei, zscorej;
  double       zscore;
  int          i, j, k;
  int          ipair;

  for (i = 0; i < mi->alen; i ++) {
    if (ct[i+1] > 0 && ct[i+1] > i+1) {
      ipair = ct[i+1]-1;

      avgi = 0.0;
      stdi = 0.0;
      
      for (j = 0; j < mi->alen; j ++) {
	if (j != ipair) {   
	  if (j != i) {
	    avgi += mi->COV->mx[i][j];
	    stdi += mi->COV->mx[i][j] * mi->COV->mx[i][j];
	  }
	}
	else {
	  avgj = 0.0;
	  stdj = 0.0;
	  for (k = 0; k < mi->alen; k ++) {
	    if (k != j && k != i) {       
	      avgj += mi->COV->mx[j][k];
	      stdj += mi->COV->mx[j][k] * mi->COV->mx[j][k];
	    }
	  }
	  avgj /= (mi->alen-2);
	  stdj /= (mi->alen-2);
	  stdj -= avgj*avgj;
	  stdj = sqrt(stdj);
	}
      }
      avgi /= (mi->alen-2);
      stdi /= (mi->alen-2);
      stdi -= avgi*avgi;
      stdi = sqrt(stdi);

      zscorei = (mi->COV->mx[i][ipair] - avgi) / stdi;
      zscorej = (mi->COV->mx[i][ipair] - avgj) / stdj;
      zscore  = ESL_MIN(zscorej, zscorej);
      printf("[%d][%d] %f | %f | %f %f | %f %f\n", msamap[i], msamap[ipair], zscore, mi->COV->mx[i][ipair], avgi, stdi, avgj, stdj);
    }
  }  
  return eslOK;
}

int
COV_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen)
{
  double pval = 0.0;
  double add, add_fixed;
  double factorial_total, factorial_BP, factorial_NBP;
  double factorial_cnbp,  factorial_ncnbp;
  double factorial_cbp,   factorial_ncbp;
  double factorial_cov,   factorial_non;
  double tol = 1e-3;
  int    NBP;
  int    total;
  int    cov;
  int    non;
  int    cbp, ncbp;
  int    cnbp, ncnbp;
 
  NBP   = alen * (alen - 1) / 2;
  total = BP    + NBP;
  cov   = cBP   + cNBP;
  non   = total - cov;

  esl_stats_LogGamma(total+1, &factorial_total);
  esl_stats_LogGamma(BP+1,    &factorial_BP);
  esl_stats_LogGamma(NBP+1,   &factorial_NBP);
  esl_stats_LogGamma(cov+1,   &factorial_cov);
  esl_stats_LogGamma(non+1,   &factorial_non);

  add_fixed  = factorial_BP + factorial_NBP + factorial_cov + factorial_non - factorial_total;

  for (cnbp = cNBP; cnbp <= NBP; cnbp ++) {
    cbp   = cov - cnbp; if (cbp < 0) break;
    ncbp  = BP  - cbp;
    ncnbp = non - ncbp;
    
    esl_stats_LogGamma(cbp+1,   &factorial_cbp);
    esl_stats_LogGamma(ncbp+1,  &factorial_ncbp);
    esl_stats_LogGamma(cnbp+1,  &factorial_cnbp);
    esl_stats_LogGamma(ncnbp+1, &factorial_ncnbp);

    add = add_fixed - factorial_cbp - factorial_ncbp - factorial_cnbp - factorial_ncnbp;

    pval += exp(add);
  }
  if (pval > 1.0 && pval < 1.0 + tol) pval = 1.0;

  *ret_pval = pval;

  return eslOK;
}

int
COV_CYKCOVCT(FILE *outfp, char *gnuplot, char *dplotfile, char *R2Rcykfile, char *R2Rversion, int R2Rall,  ESL_RANDOMNESS *r, ESL_MSA **omsa, struct mutual_s *mi, 
		int *msamap, int minloop, enum grammar_e G, int maxFP, double expectFP, int nbpairs, char *errbuf, int verbose)
{
  ESL_MSA *msa = *omsa;
  HITLIST *hitlist = NULL;
  int     *cykct = NULL;
  char    *ss = NULL;	
  SCVAL    sc;
  int      status;
            
  /* calculate the cykcov ct vector */
  status = CYKCOV(r, mi, &cykct, &sc, minloop, maxFP, expectFP, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  if (verbose) printf("cykcov score = %f\n", sc);

  /* impose the ct on the msa GC line 'cons_ss' */
  ESL_ALLOC(ss, sizeof(char) * (msa->alen+1));
  esl_ct2simplewuss(cykct, msa->alen, ss);
  /* replace the 'SS_cons' GC line with the new ss */
  esl_sprintf(&(msa->ss_cons), "%s", ss);  
  
  /* redo the hitlist since the ct has now changed */
  status = COV_SignificantPairs_Ranking(NULL, &hitlist, mi, msamap, cykct, outfp, NULL, NULL, maxFP, expectFP, nbpairs, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  /* R2R */
  status = COV_R2R(NULL, R2Rversion, R2Rall, &msa, cykct, msamap, hitlist, FALSE, verbose, errbuf);
  if (status != eslOK) goto ERROR;
  esl_msa_Digitize(mi->abc, msa, errbuf);

  esl_ct2simplewuss(cykct, msa->alen, ss);
  printf("ss:%s\n", ss);
  
  /* expand the CT with compatible/stacked A:U C:G G:U pairs */
  status = COV_ExpandCT(R2Rcykfile, R2Rall, r, msa, &cykct, minloop, G, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  if (verbose) eslx_msafile_Write(stdout, msa, eslMSAFILE_PFAM);
  
  /* R2Rpdf */
  status = COV_R2Rpdf(R2Rcykfile, R2Rversion, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  /* DotPlot */
  status = COV_DotPlot(gnuplot, dplotfile, msa, cykct, mi, msamap, hitlist, verbose, errbuf);
  if (status != eslOK) goto ERROR;

  *omsa = msa;

  COV_FreeHitList(hitlist);
  free(cykct);
  free(ss);
  return eslOK;
  
 ERROR:
  if (hitlist) COV_FreeHitList(hitlist);
  if (cykct)   free(cykct);
  if (ss)      free(ss);
  return status;
}

int              
COV_DotPlot(char *gnuplot, char *dplotfile, ESL_MSA *msa, int *ct, struct mutual_s *mi, int *msamap, HITLIST *hitlist, int verbose, char *errbuf)
{
  FILE   *pipe;
  char   *filename = NULL;
  double  pointsize;
  double  ps_max = 0.40;
  double  ps_min = 0.0003;
  int     L = msamap[msa->alen-1]+1;
  int     h;           /* index for hitlist */
  int     i, ipair;
  int     ih, jh;
  
  pointsize = (mi->maxCOV > 0.)? ps_max/mi->maxCOV : ps_min;

  esl_FileTail(dplotfile, TRUE, &filename);

  pipe = popen(gnuplot, "w");
  
  //fprintf(pipe, "set terminal postscript color 14 \n");
  fprintf(pipe, "set terminal svg size 350,262 fname 'Verdana' fsize 10 \n");
  fprintf(pipe, "set output '%s'\n", dplotfile);
  /* matlab's 'jet' colormap scale */
  fprintf(pipe, "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0)\n");

  fprintf(pipe, "unset key\n");
  fprintf(pipe, "set size ratio -1\n");
  fprintf(pipe, "set pointsize %f\n", pointsize);
  fprintf(pipe, "set title '%s' \n", filename);
  fprintf(pipe, "set ylabel 'alignment position'\n");
  fprintf(pipe, "set xlabel 'alignment position'\n");

  fprintf(pipe, "set style line 1   lt 1 lc rgb 'grey' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 2   lt 1 lc rgb 'brown' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 3   lt 1 lc rgb 'cyan' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 4   lt 1 lc rgb 'red' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 5   lt 1 lc rgb 'orange' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 6   lt 1 lc rgb 'turquoise' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 7   lt 1 lc rgb 'black' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 8   lt 1 lc rgb 'green' pt 7 lw 2 ps variable\n");
  fprintf(pipe, "set style line 9   lt 1 lc rgb 'blue' pt 7 lw 2 ps variable\n");

  fprintf(pipe, "set yrange [1:%d]\n", L);
  fprintf(pipe, "set xrange [1:%d]\n", L);

  fprintf(pipe, "set multiplot\n");

  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "fun(x)=x\n");
  fprintf(pipe, "plot fun(x) with lines ls 7\n");

  // the actual ct
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 9\n");
  for (i = 1; i <= msa->alen; i ++) {
    ipair = ct[i];
    if (ipair > 0) {
      fprintf(pipe, "%d %d %f\n", msamap[i-1]+1,     msamap[ipair-1]+1, (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
      fprintf(pipe, "%d %d %f\n", msamap[ipair-1]+1, msamap[i-1]+1,     (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
    }	
  } 
  fprintf(pipe, "e\n");

  // the covarying basepairs
  //fprintf(pipe, "plot '-' u 1:2:3 with image \n");
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 8 \n");
  for (h = 0; h < hitlist->nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    if (hitlist->hit[h].is_bpair) {
      fprintf(pipe, "%d %d %f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc);
      fprintf(pipe, "%d %d %f\n", msamap[jh]+1, msamap[ih]+1, hitlist->hit[h].sc);
    }	
  } 
  fprintf(pipe, "e\n");

  // covarying pairs compatible with the given structure
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 5\n");
  for (h = 0; h < hitlist->nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    if (hitlist->hit[h].is_compatible) {
      fprintf(pipe, "%d %d %f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc);	
      fprintf(pipe, "%d %d %f\n", msamap[jh]+1, msamap[ih]+1, hitlist->hit[h].sc);	
    }
  } 
  fprintf(pipe, "e\n");
  
  // covarying pairs incompatible with the given structure
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "plot '-' u 1:2:3 with points ls 7\n");
  for (h = 0; h < hitlist->nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    if (!hitlist->hit[h].is_bpair && !hitlist->hit[h].is_compatible) {
      fprintf(pipe, "%d %d %f\n", msamap[ih]+1, msamap[jh]+1, hitlist->hit[h].sc);	
      fprintf(pipe, "%d %d %f\n", msamap[jh]+1, msamap[ih]+1, hitlist->hit[h].sc);	
    }
  } 
  fprintf(pipe, "e\n");
  
  pclose(pipe);
  
  free(filename);
  return eslOK;
}

int
COV_R2R(char *r2rfile, char *r2rversion, int r2rall, ESL_MSA **ret_msa, int *ct, int *msamap, HITLIST *hitlist, int makepdf, int verbose, char *errbuf)
 {
  ESLX_MSAFILE *afp = NULL;
  FILE         *fp = NULL;
  ESL_MSA      *msa = *ret_msa;
  ESL_MSA      *r2rmsa = NULL;
  char         *r2rpdf = NULL;
  char          tmpinfile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char          covtag[12] = "cov_SS_cons";
  char         *args = NULL;
  char         *s = NULL;
  char         *ssstr = NULL;
  char         *covstr = NULL;
  char         *tok;
  int           found;
  int           i;
  int           h;
  int           ih, jh;
  int           tagidx;
  int           status;
 
  /* first modify the ss to a simple <> format. R2R cannot deal with fullwuss 
   */
  ESL_ALLOC(ssstr, sizeof(char) * (msa->alen+1));
  esl_ct2simplewuss(ct, msa->alen, ssstr);
  
  /* replace the 'SS_cons' GC line with the new ss */
  esl_sprintf(&(msa->ss_cons), "%s", ssstr);  
  
  /* R2R input and output in PFAM format (STOCKHOLM in one single block) */
  if ((status = esl_tmpfile_named(tmpinfile,  &fp))                      != eslOK) ESL_XFAIL(status, errbuf, "failed to create input file");
  if ((status = eslx_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_PFAM)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write PFAM file\n");
  fclose(fp);
  
  /* run R2R */
  if ("R2RDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("R2RDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpoutfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  esl_sprintf(&args, "%s/%s/src/r2r --GSC-weighted-consensus %s %s 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1", s, r2rversion, tmpinfile, tmpoutfile);
  system(args);
  fclose(fp);
 
  /* convert output to r2rmsa */
  if (eslx_msafile_Open(NULL, tmpoutfile, NULL, eslMSAFILE_PFAM, NULL, &afp) != eslOK) eslx_msafile_OpenFailure(afp, status);
  afp->format = eslMSAFILE_PFAM;
  if (eslx_msafile_Read(afp, &r2rmsa) != eslOK) eslx_msafile_ReadFailure(afp, status);
  eslx_msafile_Close(afp);

  /* modify the cov_cons_ss line acording to our hitlist */
  if (msa->alen != r2rmsa->alen) ESL_XFAIL(eslFAIL, errbuf, "r2r has modified the alignment\n");
  for (i = 1; i <= msa->alen; i ++) {
    found = FALSE;
    for (h = 0; h < hitlist->nhit; h ++) {
      ih = hitlist->hit[h].i+1;
      jh = hitlist->hit[h].j+1;

      if ((i == ih || i == jh) && hitlist->hit[h].is_bpair) { 
	esl_sprintf(&tok, "2"); 
	found = TRUE; 
      }
      
      if (found) break;
    }
    if (!found) esl_sprintf(&tok, ".");  

    if (i == 1) esl_sprintf(&covstr, "%s", tok);
    else        esl_sprintf(&covstr, "%s%s", covstr, tok);
  }
 
  /* add line #=GF R2R keep allpairs 
   * so that it does not truncate ss.
   * cannot use the standard esl_msa_addGF:
   *             esl_msa_AddGF(msa, "R2R", -1, " keep allpairs", -1);
   * since it does not parse with r2r
   *
   * turns out the above solution can only deal with the  <> annotation
   */
  if (r2rall) 
    esl_msa_AddGF(r2rmsa, "R2R keep all", -1, "", -1);
  else
    esl_msa_AddGF(r2rmsa, "R2R keep allpairs", -1, "", -1);
  
  /* replace the r2r 'cov_SS_cons' GC line with our own */
  for (tagidx = 0; tagidx < r2rmsa->ngc; tagidx++)
    if (strcmp(r2rmsa->gc_tag[tagidx], covtag) == 0) break;
  if (tagidx == r2rmsa->ngc) {
    ESL_REALLOC(r2rmsa->gc_tag, (r2rmsa->ngc+1) * sizeof(char **));
    ESL_REALLOC(r2rmsa->gc,     (r2rmsa->ngc+1) * sizeof(char **));
    r2rmsa->gc[r2rmsa->ngc] = NULL;
    r2rmsa->ngc++;
  }
  if ((status = esl_strdup(covtag, -1, &(r2rmsa->gc_tag[tagidx]))) != eslOK) goto ERROR;
  esl_sprintf(&(r2rmsa->gc[tagidx]), "%s", covstr);

  if (verbose) eslx_msafile_Write(stdout, r2rmsa, eslMSAFILE_PFAM);
  
  /* write the R2R annotated to PFAM format */
  if (r2rfile) {
    if ((fp = fopen(r2rfile, "w")) == NULL) esl_fatal("Failed to open output file %s", r2rfile);
    eslx_msafile_Write(fp, r2rmsa, eslMSAFILE_PFAM);
    fclose(fp);
    
    /* produce the R2R pdf */
    if (makepdf) {
      status = COV_R2Rpdf(r2rfile, r2rversion, verbose, errbuf);
      if (status != eslOK) goto ERROR;
    }
  }
  
  remove(tmpinfile);
  remove(tmpoutfile);
  
  *ret_msa = r2rmsa;
  esl_msa_Destroy(msa);
  if (r2rpdf) free(r2rpdf);
  if (args) free(args);
  if (ssstr)  free(ssstr);
  if (covstr) free(covstr);
  return eslOK;

 ERROR:
  remove(tmpinfile);
  remove(tmpoutfile);
  
  if (msa)    esl_msa_Destroy(msa);
  if (r2rmsa) esl_msa_Destroy(r2rmsa);
  if (r2rpdf) free(r2rpdf);
  if (args)   free(args);
  if (ssstr)  free(ssstr);
  if (covstr) free(covstr);
  return status;
}

int
COV_R2Rpdf(char *r2rfile, char *r2rversion, int verbose, char *errbuf)
{
  char *r2rpdf = NULL;
  char *args = NULL;
  char *s = NULL;

  /* produce the R2R pdf */
  if ("R2RDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("R2RDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&r2rpdf, "%s.pdf", r2rfile);
  esl_sprintf(&args, "%s/%s/src/r2r %s %s >/dev/null", s, r2rversion, r2rfile, r2rpdf);
  //esl_sprintf(&args, "%s/%s/src/r2r %s %s ", s, r2rversion, r2rfile, r2rpdf);
  system(args);
  
  free(args);
  free(r2rpdf);
  
  return eslOK;
 }

int
COV_ExpandCT(char *r2rfile, int r2rall, ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, int minloop, enum grammar_e G, int verbose, char *errbuf)
{
  FILE *fp = NULL;
  char *ss = NULL;
  int   tagidx;
  int   L = msa->alen;
  int   status;
  
  /* replace any R2R line with a new one not tab delimited 
   * R2R chockes on tab delimited lines
   */
  for (tagidx = 0; tagidx < msa->ngf; tagidx++) 
    if (strcmp(msa->gf_tag[tagidx], "R2R") == 0) {
      esl_sprintf(&(msa->gf_tag[tagidx]), "R2R %s", msa->gf[tagidx]);
      esl_sprintf(&(msa->gf[tagidx]), "");
  }

#if 0 // naive method
  status = COV_ExpandCT_Naive(msa, *ret_ct, minloop, verbose, errbuf);
#else // covariance-constrain CYK using a probabilistic grammar
  status = COV_ExpandCT_CCCYK(r, msa, ret_ct, G, minloop, verbose, errbuf);
#endif
  if (status != eslOK) goto ERROR;

  /* replace the 'SS_cons' GC line with the new ss */
  ESL_ALLOC(ss, sizeof(char) * (L+1));
  esl_ct2simplewuss(*ret_ct, L, ss);
  esl_sprintf(&(msa->ss_cons), "%s", ss);  
  //printf("ss:%s\n", msa->ss_cons);

  if ((fp = fopen(r2rfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to open output file %s", r2rfile);
  eslx_msafile_Write(fp, msa, eslMSAFILE_PFAM);
  fclose(fp);
  
  free(ss);
  return eslOK;

 ERROR:
  if (ss) free(ss);
  return status;
}

int
COV_ExpandCT_Naive(ESL_MSA *msa, int *ct, int minloop, int verbose, char *errbuf)
{
  char  tag[10] = "cons";
  char *cons;
  int   tagidx;
  int   L = msa->alen;
  int   i, j, d;
  
  /* get the line #=GC cons */
  for (tagidx = 0; tagidx < msa->ngc; tagidx++)
    if (strcmp(msa->gc_tag[tagidx], tag) == 0) break;
  if (tagidx == msa->ngc) return eslOK; // no cons line to expand the CT
  cons =  msa->gc[tagidx];

  for (j = 0; j < L; j++)
    for (d = 0; d <= j; d++)
      {
	i = j-d+1;

	if ( d >= minloop && is_stacked_pair(i+1, j+1, L, ct) && is_cannonical_pair(cons[i], cons[j]) )
	  { // add this pair
	    //printf("%c %c | %d %d %d\n", cons[i], cons[j], i, j, d);
	    ct[i+1] = j+1;
	    ct[j+1] = i+1;
	  }	
      }

return eslOK;
}

int
COV_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct,  enum grammar_e G, int minloop, int verbose, char *errbuf)
{
  char  tag[10] = "RF_cons";
  char    *rfline = NULL;
  ESL_SQ  *sq = NULL;
  int     *ct = *ret_ct;
  int     *cct = NULL;
  SCVAL    sc;
  float    idthresh = 0.3;
  int      tagidx;
  int      status;
  
  /* create an RF sequence */
  ESL_ALLOC(rfline, sizeof(char) * (msa->alen+1));
  esl_msa_ReasonableRF(msa, idthresh, TRUE, rfline);
  sq = esl_sq_CreateFrom(msa->name,  rfline, msa->desc, msa->acc, msa->ss_cons); 
  
  printf("sq:%s\n", sq->seq);
  esl_sq_Digitize((const ESL_ALPHABET *)msa->abc, sq);
 
  cykcov_remove_inconsistencies(sq, ct, minloop);
 
 /* calculate the convariance-constrain CYK structure using a probabilistic grammar */
  status = COCOCYK(r, G, sq, ct, &cct, &sc, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  if (verbose) printf("coco-cyk score = %f\n", sc);

  if (cct) {
    free(ct); ct = NULL;
    *ret_ct = cct;
  }
  else *ret_ct = ct;

  if (rfline) free(rfline);
  esl_sq_Destroy(sq);
  return eslOK;

 ERROR:
  if (sq) esl_sq_Destroy(sq);
  if (rfline) free(rfline);
  if (cct) free(cct);
  return status;
}


/*---------------- internal functions --------------------- */


static int
cykcov_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop)
{
  int L = sq->n;
  int n;
  int ipair;
  int i;
  int x;

  /* remove covariation that correspond to gap-gap in
   * the RF sequence */
  for (i = 1; i <= L; i++) {
    ipair = ct[i];
    if (ipair > 0 && sq->dsq[i] >= NB && sq->dsq[ipair] >= NB) { // remove this covariation
      ct[i]     = 0;
      ct[ipair] = 0;
    }
  }

  /* remove covariation that are closer than minloop in RF 
   */
  for (i = 1; i <= L; i++) {
    ipair = ct[i];
    if (ipair == 0 || ipair < i) continue;

    n = 0;
    for (x = i+1; x < ipair; x ++) 
      if (sq->dsq[x] < NB) n ++;

    if (n < minloop-2) { // remove this covariation
      ct[i]     = 0;
      ct[ipair] = 0;
    }
  }

  return eslOK;
}
