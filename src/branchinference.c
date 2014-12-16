/* branch inference.
 *
 * Given a msa and a T and an evolutionary model, infer branch lengths.
 *  
 * 
 */

#include "p7_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h> 
#include <string.h> 
#include <sys/types.h> 

#include <math.h> 
#include <float.h> 

#include "easel.h" 
#include "esl_distance.h"
#include "esl_getopts.h" 
#include "esl_histogram.h" 
#include "esl_msa.h" 
#include "esl_msafile.h" 
#include "esl_sq.h" 
#include "esl_sqio.h" 
#include "esl_stack.h" 
#include "esl_tree.h" 
#include "esl_vectorops.h" 

#include "hmmer.h"
#include "evohmmer.h"
#include "p7_evopipeline.h"

#include "e2.h" 
#include "e2_msa.h" 
#include "e2_tree.h" 
#include "branchinference.h" 
#include "Dickersonrates.h" 
#include "fetchpfamdb.h" 
#include "homology.h"
#include "msatree.h"
#include "mya.h"
#include "plot.h"
#include "ssifile.h"


static int ehmm_tree(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE **ret_Ti, P7_HMM *hm, float evalcutoff, float tol, char *errbuf, int verbose);


int
ehmmBranchInference(FILE *outfp, char *gnuplot, char *distfile, char *histofile, char *ratefile, char *binratefile,
		    ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, SPECIES *SPE, P7_PIPELINE *pli, P7_BG *bg,
		    float evalcutoff, P7_HMM *hmm, int noopt, int do_viterbi, float incbintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose)
{
  char        *ssidistfile = NULL;
  char       **taxalist = NULL;
  FILE        *distfp = NULL;
  ESL_TREE    *Ti = NULL;                                   /* the infered tree */ 
  double       maxdt = 0.;
  float        maxmya = 10000.;
  float        maxbintime;
  int          ntaxa = 0;                                   /* those that actually have a species assigned by ncbi */
  int          x;                                           /* index over ntaxa */
  int          status;

  /* maxdt in between two species in T */
  Tree_GetNodeTime(0, T, NULL, NULL, &maxdt, errbuf, verbose);
  maxdt += 2.;
  maxdt *= 2.;
  maxbintime =  (SPE->n > 0)? maxmya : (float)maxdt;

  /* Create an ssi index of distfile */
  if (distfile != NULL) {
    status = ssi_distfile_Create(distfile, &ssidistfile, &taxalist, &ntaxa);
    if (status != eslOK) { fprintf(outfp, "ssi_dist_Create failed for file %s\n", distfile); status = eslFAIL; goto ERROR; }
    if (verbose) fprintf(outfp, "ssifile %s created\n", ssidistfile); 
    if ((distfp = fopen(distfile, "r")) == NULL) { printf("failed to open %s\n", distfile);  status = eslFAIL; goto ERROR;  }
  }
  
  /* infer the branch lengths from the HMM, fill Ti with values */
  status = ehmm_tree(r, msa, &Ti, hmm, evalcutoff, tol, errbuf, verbose);

  if (status != eslOK) goto ERROR; 
  if (1||verbose) Tree_Dump(stdout, Ti, "Infered Tree");

  /* calculate the distances per pair of sequences */
  status = binf_InferedDistances(outfp, histofile, ratefile, binratefile, gnuplot, ssidistfile,  msa, T, Ti, 
				 SPE, taxalist, ntaxa, maxbintime, incbintime, tlinear, full, evalcutoff, view, errbuf, verbose);  
  if (status != eslOK) goto ERROR; 
  
  if (taxalist) { 
    for (x = 0; x < ntaxa; x ++) if (taxalist[x]) free(taxalist[x]); 
    free(taxalist); 
  } 
  if (distfp) fclose(distfp); 
  if (ssidistfile) { 
    remove(ssidistfile);  
    free (ssidistfile); 
  } 
  esl_tree_Destroy(Ti);
  
  return eslOK;
  
 ERROR:
  if (taxalist) {
    for (x = 0; x < ntaxa; x ++) if (taxalist[x]) free(taxalist[x]);
    free(taxalist);
  }
  if (ssidistfile) {
   remove(ssidistfile); 
   free (ssidistfile);
  }
  if (Ti) esl_tree_Destroy(Ti);
  return status;
}


int
e2BranchInference(FILE *outfp, char *gnuplot, char *distfile, char *histofile, char *ratefile, char *binratefile,
		  ESL_RANDOMNESS *r, ESL_MSA *msa, float *msafrq, ESL_TREE *T, SPECIES *SPE, E2_PIPELINE *pli, E1_BG *bg, P7_BG *bg7, E2_ALI e2ali,
		  float evalcutoff, E1_RATE *R, P7_RATE *R7, int mode, int do_viterbi, float incbintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose)
{
  char        *ssidistfile = NULL;
  char       **taxalist = NULL;
  FILE        *distfp = NULL;
  ESL_TREE    *Ti = NULL;                                   /* the infered tree */ 
  double       maxdt = 0.;
  float        maxmya = 10000.;
  float        maxbintime;
  int          ntaxa = 0;                                   /* those that actually have a species assigned by ncbi */
  int          x;                                           /* index over ntaxa */
  int          status;

  /* maxdt in between two species in T */
  Tree_GetNodeTime(0, T, NULL, NULL, &maxdt, errbuf, verbose);
  maxdt += 2.;
  maxdt *= 2.;
  maxbintime =  (SPE->n > 0)? maxmya : (float)maxdt;

  /* Create an ssi index of distfile */
  if (distfile != NULL) {
    status = ssi_distfile_Create(distfile, &ssidistfile, &taxalist, &ntaxa);
    if (status != eslOK) { fprintf(outfp, "ssi_dist_Create failed for file %s\n", distfile); status = eslFAIL; goto ERROR; }
    if (verbose) fprintf(outfp, "ssifile %s created\n", ssidistfile); 
    if ((distfp = fopen(distfile, "r")) == NULL) { printf("failed to open %s\n", distfile);  status = eslFAIL; goto ERROR;  }
  }
  
  /* infer the branch lengths for the tree */
  status = e2_tree_UPGMA(&Ti, msa, msafrq, r, pli, R, R7, bg, bg7, e2ali, mode, do_viterbi, -1.0, -1.0, tol, errbuf, verbose);
  if (status != eslOK) goto ERROR; 
  if (1||verbose) Tree_Dump(stdout, Ti, "Infered Tree");

  /* calculate the distances per pair of sequences */
  status = binf_InferedDistances(outfp, histofile, ratefile, binratefile, gnuplot, ssidistfile,  msa, T, Ti, 
				 SPE, taxalist, ntaxa, maxbintime, incbintime, tlinear, full, evalcutoff, view, errbuf, verbose);  
  if (status != eslOK) goto ERROR; 
  
  if (taxalist) { 
    for (x = 0; x < ntaxa; x ++) if (taxalist[x]) free(taxalist[x]); 
    free(taxalist); 
  } 
  if (distfp) fclose(distfp); 
  if (ssidistfile) { 
    remove(ssidistfile);  
    free (ssidistfile); 
  } 
  esl_tree_Destroy(Ti);
  
  return eslOK;
  
 ERROR:
  if (taxalist) {
    for (x = 0; x < ntaxa; x ++) if (taxalist[x]) free(taxalist[x]);
    free(taxalist);
  }
  if (ssidistfile) {
   remove(ssidistfile); 
   free (ssidistfile);
  }
  if (Ti) esl_tree_Destroy(Ti);
  return status;
}


int 
binf_InferedDistances(FILE *outfp,  char *histofile, char *ratefile, char *binratefile, char *gnuplot, char *ssidistfile, 
		      const ESL_MSA *msa, ESL_TREE *T, ESL_TREE *Ti, SPECIES *SPE, char **taxalist, int ntaxa, 
		      float maxbintime, float incbintime, float tlinear, int full, float evalcutoff, int view, char *errbuf, int verbose) 
{
  ESL_HISTOGRAM *h = NULL;
  double        *evalues = NULL;
  double        *x_mya = NULL;
  double        *y_rate = NULL;
  double        *y_e2rate = NULL;
  double        *yb_rate_ave = NULL;
  double        *yb_rate_std = NULL;
  double        *yb_e2rate_ave = NULL;
  double        *yb_e2rate_std = NULL;
  float          chrate;
  float          cchrate;
  float          mya;
  float          dti;
  float          dt;
  float          sc;      /* homology score */
  float          eval;    /* homology evalue */
  int            lca;
  int            n, m;
  int            N = 0;
  int            bx;
  int            x;
  int            status;

  h = (SPE->n > 0)? esl_histogram_Create(-1.0, maxbintime, incbintime) : esl_histogram_Create(-1.0, maxbintime, incbintime);

  ESL_ALLOC(yb_rate_ave,   sizeof(double) * h->nb); esl_vec_DSet(yb_rate_ave,   h->nb, 0.0);          
  ESL_ALLOC(yb_rate_std,   sizeof(double) * h->nb); esl_vec_DSet(yb_rate_std,   h->nb, 0.0);        
  ESL_ALLOC(yb_e2rate_ave, sizeof(double) * h->nb); esl_vec_DSet(yb_e2rate_ave, h->nb, 0.0);            
  ESL_ALLOC(yb_e2rate_std, sizeof(double) * h->nb); esl_vec_DSet(yb_e2rate_std, h->nb, 0.0);        
  
  ESL_ALLOC  (evalues,  sizeof(double) * (N+1));
  ESL_ALLOC  (x_mya,    sizeof(double) * (N+1));
  ESL_ALLOC  (y_rate,   sizeof(double) * (N+1));
  ESL_ALLOC  (y_e2rate, sizeof(double) * (N+1));
  
  /* force it to go by (0,0) */
  x_mya[N]    = 0.0;
  y_rate[N]   = 0.0;
  y_e2rate[N] = 0.0;
  evalues[N]  = 0.0;
  esl_histogram_Add(h, 0.0); 
  
  esl_histogram_Score2Bin(h, 0.0, &bx);
  yb_rate_ave[bx]   += y_rate[N];
  yb_rate_std[bx]   += y_rate[N] * y_rate[N];
  yb_e2rate_ave[bx] += y_e2rate[N];
  yb_e2rate_std[bx] += y_e2rate[N] * y_e2rate[N];
  N++;
  
  for (n = 0; n < msa->nseq-1; n ++) {
    for (m = n+1; m < msa->nseq; m ++) { 
      status = Tree_FindLowerCommonAncestor(-n, -m, T,  &lca, &dt);                                  if (status != eslOK) goto ERROR;
      status = Tree_FindLowerCommonAncestor(-n, -m, Ti, NULL, &dti);                                 if (status != eslOK) goto ERROR;
      status = naiveDickersonRates(msa->abc, msa->ax[n], msa->ax[m], msa->alen, &chrate, &cchrate);       if (status != eslOK) goto ERROR;
    
      status = homol_RunPHMMER(n, m, msa, &sc, &eval, errbuf, verbose); if (status != eslOK) goto ERROR;
            
      if (SPE->n > 0) {
	if (ssidistfile) { /* get mya from distfile */
	  status = mya_BetweenSpeciesFromSSIFile(ssidistfile, taxalist, ntaxa, SPE->parent[n][SPE->parentype], SPE->parent[m][SPE->parentype], &mya, FALSE); 
	  if (status != eslOK) goto ERROR;
 	}
	else {          /* calculate mya on the spot */
	  status = mya_BetweenLevels(n, m, SPE, &mya, errbuf, verbose); if (status != eslOK) goto ERROR;
 	}
 	
	if (mya > maxbintime) ESL_XFAIL(eslFAIL, errbuf, "mya %f larger than maxtime in histogram %f\n", mya, maxbintime);

	fprintf(outfp,  "%f %f %f %f %f %f %f %s %s %s %s \n", 
		mya, chrate/100., cchrate/100., dt, dti, sc, eval, msa->sqname[n], msa->sqname[m], SPE->spname[n], SPE->spname[m]);
	if (mya >= 0) {
	  ESL_REALLOC(evalues,  sizeof(double) * (N+1));
	  ESL_REALLOC(x_mya,    sizeof(double) * (N+1));
	  ESL_REALLOC(y_rate,   sizeof(double) * (N+1));
	  ESL_REALLOC(y_e2rate, sizeof(double) * (N+1));
	  x_mya[N]    = (double)mya;
	  y_rate[N]   = (double)cchrate/100.;
	  y_e2rate[N] = (double)dti;
	  evalues[N]  = (double)eval;
	  esl_histogram_Add(h, mya); 
	  
	  esl_histogram_Score2Bin(h, (double)mya, &bx);
	  yb_rate_ave[bx]   += y_rate[N];
	  yb_rate_std[bx]   += y_rate[N] * y_rate[N];
	  yb_e2rate_ave[bx] += y_e2rate[N];
	  yb_e2rate_std[bx] += y_e2rate[N] * y_e2rate[N];
	  N++;
	}
      }
      else {
	if (dt > maxbintime) ESL_XFAIL(eslFAIL, errbuf, "time %f larger than maxtime in histogram %f\n", dt, maxbintime);

	fprintf(outfp,  "%f %f %f %f %f %f %f %s %s \n", 
		dt, chrate/100., cchrate/100., dt, dti, sc, eval, msa->sqname[n], msa->sqname[m]);
	ESL_REALLOC(evalues,  sizeof(double) * (N+1));
	ESL_REALLOC(x_mya,    sizeof(double) * (N+1));
	ESL_REALLOC(y_rate,   sizeof(double) * (N+1));
	ESL_REALLOC(y_e2rate, sizeof(double) * (N+1));
	x_mya[N]    = (double)dt;
	y_rate[N]   = (double)cchrate/100.;
	y_e2rate[N] = (double)dti;
	evalues[N]  = (double)eval;
	esl_histogram_Add(h, (double)dt); 
	
	esl_histogram_Score2Bin(h, (double)dt, &bx);
	yb_rate_ave[bx]   += y_rate[N];
	yb_rate_std[bx]   += y_rate[N] * y_rate[N];
	yb_e2rate_ave[bx] += y_e2rate[N];
	yb_e2rate_std[bx] += y_e2rate[N] * y_e2rate[N];
	N++;
      }      
    }
  }
  
  /* normalize the binned data */
  for (x = 0; x < h->nb; x ++) {
    if (h->obs[x] > 0) {
      yb_rate_ave[x]   /= (double) h->obs[x];
      yb_rate_std[x]    = sqrt((yb_rate_std[x]   - yb_rate_ave[x]  *yb_rate_ave[x]  *(double)h->obs[x]) / ((double)h->obs[x]));
      yb_e2rate_ave[x] /= (double) h->obs[x];
      yb_e2rate_std[x]  = sqrt((yb_e2rate_std[x] - yb_e2rate_ave[x]*yb_e2rate_ave[x]*(double)h->obs[x]) / ((double)h->obs[x]));
    }
  }
  
  if ((status = plot_write_Histogram(histofile, h))                                                                                                   != eslOK) goto ERROR;
  if ((status = plot_write_RatesWithLinearRegression(ratefile, evalues, x_mya, y_rate, y_e2rate, tlinear, N, NULL))                                   != eslOK) goto ERROR;
  if ((status = plot_write_BinnedRatesWithLinearRegression(binratefile, h, yb_rate_ave, yb_rate_std, yb_e2rate_ave, yb_e2rate_std, tlinear, N, NULL)) != eslOK) goto ERROR;

  esl_histogram_Destroy(h);
  free(evalues);
  free(x_mya);
  free(y_rate);
  free(y_e2rate);
  free(yb_rate_ave);
  free(yb_rate_std);
  free(yb_e2rate_ave);
  free(yb_e2rate_std);
  return eslOK;
  
 ERROR:
  if (h)             esl_histogram_Destroy(h);
  if (evalues)       free(evalues);
  if (x_mya)         free(x_mya);
  if (y_rate)        free(y_rate);
  if (y_e2rate)      free(y_e2rate); 
  if (yb_rate_ave)   free(yb_rate_ave);
  if (yb_rate_std)   free(yb_rate_std);
  if (yb_e2rate_ave) free(yb_e2rate_ave);
  if (yb_e2rate_std) free(yb_e2rate_std);  
  return status;
}


static int 
ehmm_tree(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE **ret_T, P7_HMM *hmm, float evalcutoff, float tol, char *errbuf, int verbose)
{
  char          tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char          tmptblfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char          tmphmmfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpseqfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  FILE         *fp  = NULL;
  ESL_SQ       *sq    = NULL;
  ESL_DMATRIX  *D    = NULL;
  ESL_TREE     *T    = NULL;
  float        *time = NULL;
  char         *args = NULL;
  char         *s    = NULL;
  double        Eval    = 200;
  double        incdomE = 200;
  double        incE    = 200;
  int           N = msa->nseq;
  int           n, m;
  int           status;

  if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create out file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmptblfile,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create tbl file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmphmmfile,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create hmm file");
  if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm))  != eslOK) ESL_XFAIL(status, errbuf, "failed to write hmm file");
  fclose(fp);

  /* distances */
  D = esl_dmatrix_Create(N, N);
  esl_dmatrix_Set(D,    0.0);
  
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;

  ESL_ALLOC(time, sizeof(float) * N);  
  for (n = 0; n < N; n ++) {
    esl_sq_FetchFromMSA(msa, n, &sq);
    
    if ((status = esl_tmpfile_named(tmpseqfile,  &fp))        != eslOK) ESL_XFAIL(status, errbuf, "failed to create seq file");
    if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, 0)) != eslOK) ESL_XFAIL(status, errbuf, "failed to write seq file");
    fclose(fp);
    
    esl_sprintf(&args, "%s/src/programs/ehmmsearch              --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s %s %s > %s", 
		s,  Eval, incdomE, incE, tmptblfile, tmphmmfile, tmpseqfile, tmpoutfile);  
    system(args);
    
    remove(tmpseqfile);
    esl_sq_Destroy(sq); sq = NULL;
  }
  
  for (n = 0; n < N; n ++) 
    for (m = n+1; m < N; m ++) 
      D->mx[n][m] = D->mx[m][n] = time[n] + time[m];

  if (verbose) esl_dmatrix_Dump(stdout, D, NULL, NULL);
  if (esl_tree_UPGMA(D, &T) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to create tree\n");
  if (verbose) Tree_Dump(stdout, T, "the guide Tree");

  esl_dmatrix_Destroy(D); 
  free(time);

  remove(tmpoutfile);
  remove(tmptblfile);
  remove(tmphmmfile);

  *ret_T = T;
  return eslOK;

 ERROR:
  if (time) free(time);
  if (D) esl_dmatrix_Destroy(D);
  if (T) esl_tree_Destroy(T);
  ret_T = NULL;
  return status;
}
