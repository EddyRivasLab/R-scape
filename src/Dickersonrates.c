/* Dickerson rates
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
#include "esl_getopts.h" 
#include "esl_histogram.h" 
#include "esl_msa.h" 
#include "esl_msafile.h" 
#include "esl_stack.h" 
#include "esl_tree.h" 
#include "esl_vectorops.h" 

#include "hmmer.h" 

#include "e2.h" 
#include "e2_pipeline.h"
#include "e2_trace.h"

#include "Dickersonrates.h" 
#include "fetchpfamdb.h" 
#include "homology.h"
#include "msatree.h"
#include "mya.h"
#include "plot.h"
#include "ssifile.h"

static int e2_changerates(int c, int a, const ESL_TREE *T, PSQ *sqc, PSQ *sqa, E2_TRACE *tr,  
 			  float *ret_srate, float *ret_drate, float *ret_brate, float *ret_irate); 

int
e2DickersonChangeRates(FILE *outfp, char *gnuplot, char *distfile, char *histofile, char *ratefile, char *binratefile,
		       ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, SPECIES *SPE, E2_PIPELINE *pli, E2_ALI e2ali, 
		       float evalcutoff, E1_RATE *R,  P7_RATE *R7, E1_BG *bg, P7_BG *bg7, 
		       float incbintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose)
{
  char        *ssidistfile = NULL;
  char       **taxalist = NULL;
  FILE        *distfp = NULL;
  ESL_STACK   *vs = NULL;	                            /* node index stack */
  ESL_SQ      *sq = NULL;
  E2_TRACE   **tr = NULL;
  PSQ        **sqn = NULL;                                  /* a profile sequence for each internal node */
  PSQ         *sql = NULL;                                  /* convenience pointer to sq in left branch */
  PSQ         *sqr = NULL;                                  /* convenience pointer to sq in left branch */
  double       maxdt = 0.;
  float        maxmya = 10000.;
  float        maxbintime;
  int          do_viterbi = FALSE;
  int          nnodes;
  int          ntaxa = 0;                                   /* those that actually have a species assigned by ncbi */
  int          x;                                           /* index over ntaxa */
  int          v;
  int          which;
  int          status;

  nnodes = (T->N > 1)? T->N-1 : T->N;
  /* allocate the profile sequence for internal nodes */
  ESL_ALLOC(tr,  sizeof(E2_TRACE *) * (nnodes));
  ESL_ALLOC(sqn, sizeof(PSQ      *) * (nnodes));
  for (v = 0; v < nnodes; v ++) {
    sqn[v] = NULL;
    tr[v]  = NULL;
  }

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
  
  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      if (T->left[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->left[v];	
	status = esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	if (status != eslOK) { printf("esl_sq_FetchFromMSA() failed\n");  goto ERROR; }	
  
	switch(e2ali) {
	case E2:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);	
	  break;
	case E2HMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);	
	  break;
	case E2F:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, msa->ax[which], msa->alen);
	  break;
	case E2FHMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	  break;
	default:
	  status = eslFAIL; printf("unknown alicase\n"); goto ERROR;
	}	
	sql = sqn[nnodes+which];
	esl_sq_Destroy(sq); sq = NULL;
     }
      else sql = sqn[T->left[v]];
        
      if (T->right[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->right[v];
	status = esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	if (status != eslOK) { printf("esl_sq_FetchFromMSA() failed\n");  goto ERROR; }
	
	switch(e2ali) {
	case E2:
	if (e2ali == E2 || e2ali == E2HMMER) 
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	  break;
	case E2HMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	  break;
	case E2F:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, msa->ax[which], msa->alen);
	  break;
	case E2FHMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
	  break;
	default:
	  status = eslFAIL; printf("unknown alicase\n"); goto ERROR;
	}
	
	sqr = sqn[nnodes+which];
	esl_sq_Destroy(sq); sq = NULL;
      }
      else sqr = sqn[T->right[v]];
  

      if (sql != NULL && sqr != NULL) { /* ready to go: find ancestral profile sq running the e2 algorithm */
	if (verbose) 
	  fprintf(outfp, "\nNODE %d parent %d | l:%s %d (%f,len=%d) r:%s %d (%f,len=%d)\n", 
		  v, T->parent[v], sql->name, T->left[v], T->ld[v], (int)sql->n, sqr->name, T->right[v], T->rd[v], (int)sqr->n);
	
	/* use the bg frequencies for insertions */
	status = e2_Pipeline(r, pli, sql, sqr, NULL, R, R7, bg, bg7, T->ld[v], T->rd[v], &(sqn[v]), &(tr[v]), NULL, NULL, e2ali, e2_GLOBAL, TRUE, do_viterbi, tol, errbuf, verbose);
	if (status != eslOK) goto ERROR; 
	
	/* push parent into stack unless already at the root */
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslFAIL; goto ERROR; }; 
      }
      else if (sql == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])   != eslOK) { status = eslFAIL; goto ERROR; };
      }
      else if (sqr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v])  != eslOK) { status = eslFAIL; goto ERROR; };
      }
      
      if (T->left[v]  <= 0 && sql != NULL) { psq_Destroy(sql); sql = NULL; }
      if (T->right[v] <= 0 && sqr != NULL) { psq_Destroy(sqr); sqr = NULL; }
    }
  
   status = e2DickersonRatesCluster(outfp, histofile, ratefile, binratefile, gnuplot, ssidistfile,  msa, T, (const E2_TRACE **)tr, (const PSQ **)sqn, 
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
   e1_bg_Destroy(bg); 
   esl_stack_Destroy(vs); 
   for (v = 0; v < T->N-1; v ++) { psq_Destroy(sqn[v]); e2_trace_Destroy(tr[v]); } 
   free(sqn); 
   free(tr); 

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
  if (vs) esl_stack_Destroy(vs);
  for (v = 0; v < T->N-1; v ++) { if (sqn[v]) psq_Destroy(sqn[v]); if (tr[v]) e2_trace_Destroy(tr[v]); }
  if (tr)  free(tr);
  if (sqn) free(sqn);
  return status;
}

int 
e2DickersonRatesCluster(FILE *outfp,  char *histofile, char *ratefile, char *binratefile, char *gnuplot, char *ssidistfile, 
			const ESL_MSA *msa, const ESL_TREE *T, const E2_TRACE **tr, const PSQ **sqn, SPECIES *SPE, 
			char **taxalist, int ntaxa, float maxbintime, float incbintime, float tlinear, int full, float evalcutoff, int view, char *errbuf, int verbose) 
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
  float          e2srate, e2drate;
  float          e2brate, e2irate;
  float          mya;
  float          dt;
  float          sc;      /* homology score */
  float          eval;    /* homology evalue */
  int            lca;
  int            n, m;
  int            N = 0;
  int            nn = 0;
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
      status = Tree_FindLowerCommonAncestor(-n, -m, T, &lca, &dt);                                   if (status != eslOK) goto ERROR;
      status = naiveDickersonRates(msa->abc, msa->ax[n], msa->ax[m], msa->alen, &chrate, &cchrate);  if (status != eslOK) goto ERROR;
      status = e2DickersonRates(n, m, lca, T, msa, tr, sqn, &e2srate, &e2drate, &e2brate, &e2irate); if (status != eslOK) goto ERROR;
     
      status = homol_RunPHMMER(n, m, msa, &sc, &eval, errbuf, verbose); if (status != eslOK) goto ERROR;
            
      if (SPE->n > 0) {
	if (ssidistfile) { /* get mya from distfile */
	  status = mya_BetweenSpeciesFromSSIFile(ssidistfile, taxalist, ntaxa, SPE->parent[n][SPE->parentype], SPE->parent[m][SPE->parentype], &mya, FALSE); if (status != eslOK) goto ERROR;
 	}
	else {          /* calculate mya on the spot */
	  status = mya_BetweenLevels(n, m, SPE, &mya, errbuf, verbose); if (status != eslOK) goto ERROR;
 	}
 	
	if (mya > maxbintime) ESL_XFAIL(eslFAIL, errbuf, "mya %f larger than maxtime in histogram %f\n", mya, maxbintime);

	fprintf(outfp,  "[%d]%f %f %f %f %f %f %f %f %f %s %s %s %s \n", 
		nn++, mya, chrate, cchrate, e2srate, e2drate, e2brate, e2irate, sc, eval, msa->sqname[n], msa->sqname[m], SPE->spname[n], SPE->spname[m]);
	if (mya >= 0) {
	  ESL_REALLOC(evalues,  sizeof(double) * (N+1));
	  ESL_REALLOC(x_mya,    sizeof(double) * (N+1));
	  ESL_REALLOC(y_rate,   sizeof(double) * (N+1));
	  ESL_REALLOC(y_e2rate, sizeof(double) * (N+1));
	  x_mya[N]    = (double)mya;
	  y_rate[N]   = (double)cchrate;
	  y_e2rate[N] = (double)(e2srate+e2drate+e2brate);
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

	fprintf(outfp,  "[%d]%f %f %f %f %f %f %f %f %f %s %s \n", 
		nn++, dt, chrate, cchrate, e2srate, e2drate, e2brate, e2irate, sc, eval, msa->sqname[n], msa->sqname[m]);
	ESL_REALLOC(evalues,  sizeof(double) * (N+1));
	ESL_REALLOC(x_mya,    sizeof(double) * (N+1));
	ESL_REALLOC(y_rate,   sizeof(double) * (N+1));
	ESL_REALLOC(y_e2rate, sizeof(double) * (N+1));
	x_mya[N]    = (double)dt;
	y_rate[N]   = (double)cchrate;
	y_e2rate[N] = (double)(e2srate+e2drate+e2brate);
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

int  
naiveDickersonRates(const ESL_ALPHABET *abc, const ESL_DSQ *asq1, const ESL_DSQ *asq2, int64_t L, float *ret_chrate, float *ret_cchrate) 
{ 
  float chrate  = 0.0; 
  float cchrate = 0.0; 
  int   i; 
  int   n = 0; 
  
  for (i = 1; i <= L; i++) { 
     /* canonical (x < K)  or gap (x = K) */ 
    if (asq1[i] <= abc->K && asq2[i] <= abc->K) { 
      n ++;
      if (asq1[i] == abc->K && asq2[i] == abc->K) n --; /* discard double gaps */ 
      if (asq1[i] != asq2[i]) chrate += 1.0;
    }
  }
  if (n > 0) chrate *= 100.0 / (float)n;

  cchrate = -100.0 * log (1.0 - chrate/100.);
  
  if (chrate >= 99.999999999) {
    /* an anomaly that makes this method crap, just remove it */
    cchrate = 0.0; 
  }

  *ret_chrate  = chrate;
  *ret_cchrate = cchrate;
  
  return eslOK;
}

int 
e2DickersonRates(int n, int m, int lca, const ESL_TREE *T, const ESL_MSA *msa, const E2_TRACE **tr, const PSQ **sqn, 
		 float *ret_e2srate, float *ret_e2drate, float *ret_e2brate, float *ret_e2irate)
{ 
  ESL_STACK      *vs  = NULL;
  ESL_SQ         *sq  = NULL;
  PSQ            *sq1 = NULL;
  PSQ            *sq2 = NULL;
  PSQ            *sqc;                /* convenience pointer to child sequence */
  float           totsrate = 0.0;     /* substitutions rate */
  float           totdrate = 0.0;     /* deletions rate */
  float           totbrate = 0.0;     /* insert blocks rate */
  float           totirate = 0.0;     /* insertions rate */
  float           srate;
  float           drate;
  float           brate;
  float           irate;
  int             parentn, parentm;
  int             child;
  int             v;
  int             status;

  if (( vs = esl_stack_ICreate()) == NULL) { status = eslEMEM; goto ERROR; };

  esl_sq_FetchFromMSA(msa, n, &sq); /* extract the seqs from the msa */
  sq1 = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
  esl_sq_Destroy(sq); sq = NULL;
  
  esl_sq_FetchFromMSA(msa, m, &sq); /* extract the seqs from the msa */
  sq2 = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);
  esl_sq_Destroy(sq); sq = NULL;
  
  parentn = T->taxaparent[n]; sqc = (PSQ *)sq1;           
  child = -n;
  if (esl_stack_IPush(vs, parentn) != eslOK)     { status = eslEMEM; goto ERROR; };
  while (esl_stack_IPop(vs, &v) == eslOK) 
    {       
      status = e2_changerates(child, v, T, sqc, (PSQ *)sqn[v], (E2_TRACE *)tr[v], &srate, &drate, &brate, &irate); if (status != eslOK) goto ERROR;
      totsrate += srate;
      totdrate += drate;
      totbrate += brate;
      totirate += irate;
      if (v > lca) esl_stack_IPush(vs, T->parent[v]); 
      child = v;
      sqc   = (PSQ *)sqn[v];
    }                         

  parentm = T->taxaparent[m]; sqc = (PSQ *)sq2;           
  child = -m;
  if (esl_stack_IPush(vs, parentm) != eslOK)     { status = eslEMEM; goto ERROR; };
  while (esl_stack_IPop(vs, &v) == eslOK) 
    {       
      status = e2_changerates(child, v, T, sqc, (PSQ *)sqn[v], (E2_TRACE *)tr[v], &srate, &drate, &brate, &irate); if (status != eslOK) goto ERROR;
      totsrate += srate;
      totdrate += drate;
      totbrate += brate;
      totirate += irate;
      if (v > lca) esl_stack_IPush(vs, T->parent[v]); 
      child = v;
      sqc   = (PSQ *)sqn[v];
    }                         
   
  *ret_e2srate = totsrate; 
  *ret_e2drate = totdrate; 
  *ret_e2brate = totbrate; 
  *ret_e2irate = totirate; 

  esl_stack_Destroy(vs);
  if (sq)  esl_sq_Destroy(sq);
  psq_Destroy(sq1); 
  psq_Destroy(sq2); 
  return eslOK;

 ERROR:
  if (vs)  esl_stack_Destroy(vs);
  if (sq)  esl_sq_Destroy(sq);
  if (sq1) psq_Destroy(sq1); 
  if (sq2) psq_Destroy(sq2); 
  return status;
}


int 
e2_changerates(int c, int a, const ESL_TREE *T, PSQ *sqc, PSQ *sqa, E2_TRACE *tr, float *ret_srate, float *ret_drate, float *ret_brate, float *ret_irate)
{
  float *pc = NULL;
  float *pa = NULL;
  float  srate = 0.0;           /* changes in res substitutions */
  float  drate = 0.0;           /* changes in res deletions */
  float  brate = 0.0;           /* rate of inserted blocks */
  float  irate = 0.0;           /* rate of inserted residues */
  int    isleftchild = FALSE;
  int    usei;
  int    z;
  int    prvst;
  int    n = 0;
  int    Kg = sqc->abc->K + 1;
  int    status;

  if       (T->left[a]  == c) isleftchild = TRUE;
  else if  (T->right[a] == c) isleftchild = FALSE;
  else     { status = eslFAIL; goto ERROR; }

  if (isleftchild) usei = (tr->rowsq == e2P_SL)? TRUE : FALSE;
  else             usei = (tr->rowsq == e2P_SR)? TRUE : FALSE;

  ESL_ALLOC(pc, sizeof(float) * Kg);
  ESL_ALLOC(pa, sizeof(float) * Kg);

  for (z = 0; z < tr->N; z++) 
    {
      prvst = (z>0)? tr->st[z-1] : tr->st[0];
      switch(tr->st[z]) {
      case e2T_BB:
	break;
      case e2T_IB:
	if (usei) {
	  irate += 1.0;
	if (prvst == e2T_BB) brate += 1.0;
	}
	break;
      case e2T_SS:
	psq_ProfProbs(tr->k[z],                  sqa, pa);
	psq_ProfProbs((usei)? tr->i[z]:tr->j[z], sqc, pc);
	srate += 1.0 - esl_vec_FDot(pc, pa, Kg);
	n ++;
	break;
      case e2T_DS:
	if (usei) drate += 1.0;
	else {
	  psq_ProfProbs(tr->k[z], sqa, pa);
	  psq_ProfProbs(tr->j[z], sqc, pc);
	  srate += 1.0 - esl_vec_FDot(pc, pa, Kg);
	}
	n ++;
	break;
      case e2T_IS:
	if (usei) {
	  irate += 1.0;
	  if (prvst == e2T_SS || prvst == e2T_DS) brate += 1.0;
	}
	break;
      case e2T_SD:
	if (!usei) drate += 1.0;
	else {
	  psq_ProfProbs(tr->k[z], sqa, pa);
	  psq_ProfProbs(tr->i[z], sqc, pc);
	  srate += 1.0 - esl_vec_FDot(pc, pa, Kg);
	}
	n ++;
	break;
      case e2T_DD:
	drate += 1.0;
	n ++;
	break;
      case e2T_ID:
	if (usei) {
	  irate += 1.0;
	  if (prvst == e2T_SD || prvst == e2T_DD) brate += 1.0;
	}
	break;
      case e2T_BI:
	if (!usei) {
	  irate += 1.0;
	  if (prvst == e2T_BB) brate += 1.0;
	}
	break;
      case e2T_SI:
	if (!usei) {
	  irate += 1.0;
	  if (prvst == e2T_SS || prvst == e2T_SD) brate += 1.0;
	}
	break;
      case e2T_DI:
	if (!usei) {
	  irate += 1.0;
	  if (prvst == e2T_DS || prvst == e2T_DD) brate += 1.0;
	}
	break;
      case e2T_II:
	if (!usei) {
	  irate += 1.0;
	  if (prvst == e2T_IB || prvst == e2T_IS || prvst == e2T_ID) brate += 1.0;
	}
	break;
      case e2T_EE:
	break;
      default: status = eslFAIL; goto ERROR;
      }
   }

  if (n > 0.) { 
    srate *= 100/(float)n;
    drate *= 100/(float)n;
    brate *= 100/(float)n;
    irate *= 100/(float)n;
  }
  
  *ret_srate = srate;
  *ret_drate = drate;
  *ret_brate = brate;
  *ret_irate = irate;
  
  if (pc) free(pc);
  if (pa) free(pa);  
  return eslOK;
  
 ERROR:
  if (pc) free(pc);
  if (pa) free(pa);  
  return status;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
