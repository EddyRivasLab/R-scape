/* structure.c */

#include "rscape_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_gamma.h"
#include "esl_histogram.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "contactmap.h"
#include "covariation.h"
#include "correlators.h"
#include "cococyk.h"
#include "covgrammars.h"
#include "cykcov.h"
#include "ribosum_matrix.h"
#include "r2rdepict.h"
#include "structure.h"

#define LASTFOLD 1

// paramters to include extra helices
#define  INCOMPFRAC  0.51          // max fraction of residues in a helix that overlap with another existing helix
#define  MINHELIX    4             // min length of a helix without any covarying basepairs

static int struct_cocomcyk(char *r2rfile, int r2rall,  ESL_RANDOMNESS *r, ESL_MSA *msa, SPAIR *spair, int *ret_nct, int ***ret_ctlist,
			   int minloop, enum grammar_e G, int verbose, char *errbuf);
static int struct_cocomcyk_expandct( ESL_RANDOMNESS *r, ESL_MSA *msa, SPAIR *spair, int **ret_ct, enum grammar_e G, int minloop, int verbose, char *errbuf);
static int struct_write_ss(FILE *fp, int blqsize, int nss, char **sslist);
static int ct_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop, int verbose);
static int ct_split_helices(int *ct, int L, int *ret_idx, int ***ret_ctlist, int verbose);
static int ctlist_break_in_helices(int *ret_nct, int ***ret_ctlist, int L, int verbose);
static int ctlist_helices_select(int nctcov, int **ctlistcov, int *ret_nct, int ***ret_ctlist, int L, int verbose);
static int ctlist_join(int nct, int **ctlist, int L, int **ret_ct, int verbose);

// covariation-constrained multi-CYK (COCOMCYK)
int
struct_COCOMCYK(struct data_s *data, ESL_MSA *msa, int *ret_cyknct, int ***ret_cykctlist, int minloop,
		RANKLIST *ranklist, HITLIST *hitlist, enum grammar_e G, THRESH *thresh)
{
  HITLIST       *cykhitlist = NULL;
  char          *covtype    = NULL;
  char          *threshtype = NULL;
  int          **cykctlist  = NULL;
  CLIST        **cykexclude = NULL;
  int           *cykct      = NULL;
  SCVAL          sc;
  int            cyknct;
  int            i, j;
  int            s;
  int            h;
  int            found;
  int            status;
            
  /* calculate the cykcov ct vector.
   *
   * I run a nussinov-type algorithm that incorporates as many of the significant pairs as possible.
   * These pairs become constraints for the second part of the folding in struct_ExpandCT()
   */
  status = CYKCOV(data->r, data->mi, data->clist, &cyknct, &cykctlist, &cykexclude, (hitlist)?hitlist->nhit:0, minloop, thresh, data->errbuf, data->verbose);
  if (status != eslOK) {
    for (h = 0; h < hitlist->nhit; h ++) {
      i = hitlist->hit[h].i + 1;
      j = hitlist->hit[h].j + 1;
      found = FALSE;
      for (s = 0; s < cyknct; s ++) 
	if (cykctlist[s][i] == j) { found = TRUE; break; }
      if (!found) printf("hit %d %d not found\n", i, j);
    }
    goto ERROR;
  }

  // Use the CT from covariation to do a cov-constrained multi-cyk folding
  status = struct_cocomcyk(data->R2Rcykfile, data->R2Rall, data->r, msa, data->spair, &cyknct, &cykctlist, minloop, G, data->verbose, data->errbuf);
  if (status != eslOK) goto ERROR;

  // create a new contact list from the cykct
  CMAP_ReuseCList(data->clist);
  for (s = 0; s < cyknct; s ++) {
    status = ContacMap_FromCT(data->clist, msa->alen, cykctlist[s], data->clist->mind, data->msamap, NULL);
    if (status != eslOK) goto ERROR;
  }

  /* redo the hitlist since the ct has now changed */
  corr_COVTYPEString(&covtype, data->mi->type, data->errbuf);
  cov_THRESHTYPEString(&threshtype, data->thresh->type, NULL);

  status = cov_CreateCYKHitList(data, cyknct, cykctlist, ranklist, (data->statsmethod != NAIVE)? hitlist : NULL, &cykhitlist, covtype, threshtype);
  if (status != eslOK) goto ERROR;
  
  for (s = 0; s < cyknct; s ++) {
    cykct = cykctlist[s];
    
    for (i = 1; i <= msa->alen; i ++) 
      if (cykct[i] > 0 && i < cykct[i]) data->nbpairs_cyk ++;
  }

  if (data->nbpairs_cyk > 0) {
    /* R2R */
    status = r2r_Depict(data->R2Rcykfile, data->R2Rall, msa, cyknct, cykctlist, cykhitlist, TRUE, TRUE, data->verbose, data->errbuf);
    if (status != eslOK) goto ERROR;
    
    /* DotPlots (pdf,svg) */
    status = struct_DotPlot(data->gnuplot, data->cykdplotfile, msa, cyknct, cykctlist, data->mi, data->msamap, data->firstpos, data->samplesize, cykhitlist,
			 TRUE,  data->verbose, data->errbuf);
    if (status != eslOK) goto ERROR;
    status = struct_DotPlot(data->gnuplot, data->cykdplotfile, msa, cyknct, cykctlist, data->mi, data->msamap, data->firstpos, data->samplesize, cykhitlist,
			 FALSE, data->verbose, data->errbuf);
    if (status != eslOK) goto ERROR;
  }

  *ret_cyknct    = cyknct;
  *ret_cykctlist = cykctlist;

  free(cykexclude);
  cov_FreeHitList(cykhitlist);
  free(covtype);
  free(threshtype);
  return eslOK;
  
 ERROR:
  for (s = 0; s < cyknct; s ++) { if (cykctlist[s]) free(cykctlist[s]); }
  if (cykctlist)  free(cykctlist); 
  if (cykhitlist) cov_FreeHitList(cykhitlist);
  if (covtype)    free(covtype);
  if (threshtype) free(threshtype);
  if (cykct)      free(cykct);
  return status;
}

int              
struct_DotPlot(char *gnuplot, char *dplotfile, ESL_MSA *msa, int nct, int **ctlist, struct mutual_s *mi, int *msamap, int firstpos, SAMPLESIZE samplesize,
	       HITLIST *hitlist, int dosvg, int verbose, char *errbuf)
{
  FILE    *pipe;
  char    *filename = NULL;
  char    *outplot = NULL;
  int     *ct;
  double   pointsize;
  double   ps_max = 0.40;
  double   ps_min = 0.0003;
  int      isempty;
  int      select;
  int      ileft, iright;
  int      h;           /* index for hitlist */
  int      i, ipair;
  int      ih, jh;
  int      s;
  int      status;

  if (gnuplot   == NULL) return eslOK;
  if (dplotfile == NULL) return eslOK;
  
  esl_FileTail(dplotfile, FALSE, &filename);

  pipe = popen(gnuplot, "w");
  
  if (dosvg) {
    ps_max = 1.00;
    ps_min = 0.3;
    esl_sprintf(&outplot, "%s.svg", dplotfile);
    fprintf(pipe, "set terminal svg font 'Arial,12'\n");
  }
  else {
    ps_max = 1.40;
    ps_min = 0.30;
    esl_sprintf(&outplot, "%s.ps", dplotfile);
    fprintf(pipe, "set terminal postscript color 14\n");
  }
  fprintf(pipe, "set output '%s'\n", outplot);
  pointsize = (mi->maxCOV > 0.)? ps_max/mi->maxCOV : ps_min;

  /* matlab's 'jet' colormap scale */
  fprintf(pipe, "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0)\n");

  fprintf(pipe, "unset key\n");
  fprintf(pipe, "set size ratio -1\n");
  fprintf(pipe, "set pointsize %f\n", pointsize);
  fprintf(pipe, "set title '%s' noenhanced\n", filename);
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

  ileft  = msa->alen;
  iright = 1;
  for (s = 0; s < nct; s ++) {
    ct = ctlist[s];
    for (i = 1; i <= msa->alen; i ++) {
      if (ct[i] > 0) {
	if (ct[i] > i && i < ileft)  ileft  = i;
	if (ct[i] < i && i > iright) iright = i;
      }
    }
  }
  if (ileft == msa->alen && iright == 1) {//no structure
    ileft = 1;
    iright = msa->alen;
  }
  if (ileft > iright)     ESL_XFAIL(eslFAIL, errbuf, "error in struct_DotPlot()");
  if (iright > msa->alen) ESL_XFAIL(eslFAIL, errbuf, "error in struct_DotPlot()");

  if (hitlist) {
    ileft  = (hitlist->nhit > 0)? ESL_MIN(ileft,  hitlist->hit[0].i+1) : ileft;
    iright = (hitlist->nhit > 0)? ESL_MAX(iright, hitlist->hit[hitlist->nhit-1].j+1) : iright;
  }

  fprintf(pipe, "set yrange [%d:%d]\n", msamap[ileft-1]+firstpos, msamap[iright-1]+firstpos);
  fprintf(pipe, "set xrange [%d:%d]\n", msamap[ileft-1]+firstpos, msamap[iright-1]+firstpos);

  fprintf(pipe, "set multiplot\n");

  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "fun(x)=x\n");
  fprintf(pipe, "plot fun(x) with lines ls 7\n");

  // the actual ct
  isempty = TRUE;
  for (s = 0; s < nct; s ++) {
    ct = ctlist[s];
    for (i = 1; i <= msa->alen; i ++) { if (ct[i] > 0) { isempty = FALSE; break; } }
    if (isempty == FALSE) break;
  }
  if (!isempty) {
    fprintf(pipe, "set size 1,1\n");
    fprintf(pipe, "set origin 0,0\n");  
    fprintf(pipe, "plot '-' u 1:2:3 with points ls 9\n");
    for (s = 0; s < nct; s ++) {
      ct = ctlist[s];
      for (i = 1; i <= msa->alen; i ++) {
	ipair = ct[i];
	
	if (ipair > 0) {
	  fprintf(pipe, "%d %d %f\n", msamap[i-1]+firstpos,     msamap[ipair-1]+firstpos, 
		  (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
	  fprintf(pipe, "%d %d %f\n", msamap[ipair-1]+firstpos, msamap[i-1]+firstpos,     
		  (mi->COV->mx[i-1][ipair-1]*pointsize > ps_min)? mi->COV->mx[i-1][ipair-1]:ps_min/pointsize);
	}	
      }
    }
    fprintf(pipe, "e\n");
  }
  
  // the covarying basepairs
  if (hitlist) {
    isempty = TRUE;
    for (h = 0; h < hitlist->nhit; h ++) {
      
      switch(samplesize) {
      case SAMPLE_CONTACTS:
	select = (hitlist->hit[h].bptype < BPNONE)?  TRUE : FALSE;
	break;
      case SAMPLE_BP:
	select = (hitlist->hit[h].bptype < STACKED)? TRUE : FALSE;
	break;
      case SAMPLE_WC:
	select = (hitlist->hit[h].bptype == WWc)?    TRUE : FALSE;
	break;	  
      case SAMPLE_ALL:
	select = FALSE;
	break;	  
      }
      
      if (select) { isempty = FALSE; break; }
    }
  }

  if (!isempty && hitlist) {
    fprintf(pipe, "set size 1,1\n");
    fprintf(pipe, "set origin 0,0\n");  
    fprintf(pipe, "plot '-' u 1:2:3 with points ls 8 \n");
    for (h = 0; h < hitlist->nhit; h ++) {
      ih = hitlist->hit[h].i;
      jh = hitlist->hit[h].j;

      switch(samplesize) {
      case SAMPLE_CONTACTS:
	select = (hitlist->hit[h].bptype < BPNONE)?  TRUE : FALSE;
	break;
      case SAMPLE_BP:
	select = (hitlist->hit[h].bptype < STACKED)? TRUE : FALSE;
	break;
      case SAMPLE_WC:
	select = (hitlist->hit[h].bptype == WWc)?    TRUE : FALSE;
	break;	  
      case SAMPLE_ALL:
	select = FALSE;
	break;	  
      }
      
      if (select) {
	fprintf(pipe, "%d %d %f\n", msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc);
	fprintf(pipe, "%d %d %f\n", msamap[jh]+firstpos, msamap[ih]+firstpos, hitlist->hit[h].sc);
      }	
    } 
    fprintf(pipe, "e\n");
  }

  // covarying pairs compatible with the given structure
  if (hitlist) {
    isempty = TRUE;
    for (h = 0; h < hitlist->nhit; h ++) {
      if (hitlist->hit[h].is_compatible) { isempty = FALSE; break; }
    }
    
    if (!isempty) {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");  
      fprintf(pipe, "plot '-' u 1:2:3 with points ls 5\n");
      for (h = 0; h < hitlist->nhit; h ++) {
	ih = hitlist->hit[h].i;
	jh = hitlist->hit[h].j;
	if (hitlist->hit[h].is_compatible) {
	  fprintf(pipe, "%d %d %f\n", msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc);	
	  fprintf(pipe, "%d %d %f\n", msamap[jh]+firstpos, msamap[ih]+firstpos, hitlist->hit[h].sc);	
	}
      } 
      fprintf(pipe, "e\n");
    }
  }
  
  // covarying pairs incompatible with the given structure
  if (hitlist) {
    isempty = TRUE;
    for (h = 0; h < hitlist->nhit; h ++) {
      if (!hitlist->hit[h].is_compatible) { isempty = FALSE; break; }
    }
    
    if (!isempty) {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");  
      fprintf(pipe, "plot '-' u 1:2:3 with points ls 7\n");
      for (h = 0; h < hitlist->nhit; h ++) {
	ih = hitlist->hit[h].i;
	jh = hitlist->hit[h].j;
	if (!hitlist->hit[h].is_compatible) {
	  fprintf(pipe, "%d %d %f\n", msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc);	
	  fprintf(pipe, "%d %d %f\n", msamap[jh]+firstpos, msamap[ih]+firstpos, hitlist->hit[h].sc);	
	}
      } 
      fprintf(pipe, "e\n");
    }
  }
  
  pclose(pipe);
  
  free(outplot);
  free(filename);
  return eslOK;

 ERROR:
  if (outplot) free(outplot);
  if (filename) free(filename);
  if (pipe) pclose(pipe);
  return status;
}


// take the ctlist and convert to octlist in the coordinates of the original alignment
// using the map function msamap
int
struct_CTMAP(int L, int nct, int **ctlist, int OL, int *msamap, int ***ret_octlist, char ***ret_sslist, FILE *fp, int verbose)
{
  int  **octlist = NULL;
  char **sslist  = NULL;
  int   *ct;
  char  *oss;
  int   *oct;
  int    blqsize = 60;
  int    s;
  int    i;
  int    status;

  // initialize
  ESL_ALLOC(octlist, sizeof(int  *) * nct);
  ESL_ALLOC(sslist,  sizeof(char *) * nct);
  for (s = 0; s < nct; s ++) sslist[s]  = NULL;
  for (s = 0; s < nct; s ++) octlist[s] = NULL;
  
  // the main nested structure (s=0) is annotated as SS_cons
  // The rest of the pseudoknots are annotated as SS_cons_1, SS_cons_2
  //
  // SS_cons_xx is not orthodox stockholm format.
  // I may end up changing this later.
  // It is used by R2R, and I keep it here because some of my
  // secondary/pseudoknot structure may overlap with the main nested structure.
  // That is not wrong, it is just showing our uncertainty about the structure due to lack
  // of more covariations.
  for (s = 0; s < nct; s ++) {
    ESL_ALLOC(octlist[s], sizeof(int)  * (OL+1));
    ESL_ALLOC(sslist[s],  sizeof(char) * (OL+1));
    oct = octlist[s];
    oss = sslist[s];
  
    ct = ctlist[s];
    esl_vec_ISet(oct, OL+1, 0); // initialize the oct
    
    // use the mapping to create the oct from ct
    for (i = 0; i < L; i++) 
      if (ct[i+1] > 0) oct[msamap[i]+1] = msamap[ct[i+1]-1]+1;
     
    // the structure in the coordenates of the original alignment
    esl_ct2wuss(oct, OL, oss);
  }
  if (fp)      struct_write_ss(fp,     blqsize, nct, sslist);
  if (verbose) struct_write_ss(stdout, blqsize, nct, sslist);
  
  if (ret_octlist) *ret_octlist = octlist;
  else
    {
      for (s = 0; s < nct; s ++) free(octlist[s]);
      free(octlist);
    }

  if (ret_sslist)  *ret_sslist =  sslist;
  else
    {
      for (s = 0; s < nct; s ++) free(sslist[s]);
      free(sslist);
    }
 
  return eslOK;

 ERROR:
  for (s = 0; s < nct; s ++) if (sslist[s]) free(sslist[s]);
  if (sslist) free(sslist);
  for (s = 0; s < nct; s ++) if (octlist[s]) free(octlist[s]);
  if (octlist) free(octlist);
  return status;
}

// take a ct vector possibly with pseudoknots, and separate
// into one ct without pseudoknots and additional ct's one with
// each of the pseudoknots
int
struct_SplitCT(int *ct, int L, int *ret_nct, int ***ret_ctlist, int verbose)
{
  int  **ctlist = NULL;
  char  *ss1    = NULL;
  char  *ss2    = NULL;
  int    nct = 1;
  int    use;
  int    n;
  int    c;
  int    s;
  int    status;

  ESL_ALLOC(ss1, sizeof(char)  * (L+1));
  ESL_ALLOC(ss2, sizeof(char)  * (L+1));

  esl_ct2wuss(ct, L, ss1);
  if (verbose) printf("given ss\n%s\n", ss1);
  
  // the nested structure
  ESL_ALLOC(ctlist,    sizeof(int *) * nct);
  ESL_ALLOC(ctlist[0], sizeof(int)   * (L+1));
  for (n = 0; n < L; n ++)
    {
      if (isalpha(ss1[n])) ss2[n] = '.';
      else                 ss2[n] = ss1[n];
    }
  ss2[L] = '\0';
  if (verbose) printf("given main structure\n%s\n", ss2);
  esl_wuss2ct(ss2, L, ctlist[0]);

  // the pseudoknots
  for (c = 'a'; c <= 'z'; c ++) {
    use = FALSE;
    for (n = 0; n < L; n ++)
      {
	if      (ss1[n] == c)          { ss2[n] = '>'; use = TRUE; }
	else if (ss1[n] == toupper(c))   ss2[n] = '<';
	else                             ss2[n] = '.';
      }
    ss2[L] = '\0';
    
    if (use) {
      nct ++;
      ESL_REALLOC(ctlist,      sizeof(int *) * nct);
      ESL_ALLOC(ctlist[nct-1], sizeof(int)   * (L+1));

      esl_wuss2ct(ss2, L, ctlist[nct-1]);
      if (verbose) printf("given pseudoknot %d\n%s\n", nct-1, ss2);
    }    
  }

  status = ctlist_break_in_helices(&nct, &ctlist, L, verbose);
  if (status != eslOK) goto ERROR;

  if (verbose) {
    for (s = 0; s < nct; s ++) {
      esl_ct2wuss(ctlist[s], L, ss1);
      printf("given/broken ss %d\n%s\n", s, ss1);

    }
  }
  
  *ret_nct    = nct;
  *ret_ctlist = ctlist;

  free(ss1);
  free(ss2);
  return eslOK;

 ERROR:
  if (ss1) free(ss1);
  if (ss2) free(ss2);
  for (s = 0; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  return status;
}


/*------------------------------ internal functions -----------------------------*/

static int
struct_cocomcyk(char *r2rfile, int r2rall, ESL_RANDOMNESS *r, ESL_MSA *msa, SPAIR *spair, int *ret_nct, int ***ret_ctlist, int minloop, enum grammar_e G, int verbose, char *errbuf)
{
  FILE  *fp        = NULL;
  char  *ss        = NULL;
  int    nct       = (ret_nct)?*ret_nct:0;
  int  **ctlist    = *ret_ctlist;
  int  **newctlist = NULL;
  int   *newct;
  int    new;
  int    L = msa->alen;
  int    s;
  int    i;
  int    status;
  
  if (r2rfile == NULL) return eslOK;

  ESL_ALLOC(ss,        sizeof(char)  * (L+1));
  ESL_ALLOC(newctlist, sizeof(int *) * (nct+1));
   
  for (s = 0; s < nct; s ++) {
    // covariance-constraint CYK using a probabilistic grammar
    ESL_ALLOC(newctlist[s], sizeof(int) * (L+1));
    esl_vec_ICopy(ctlist[s], L+1, newctlist[s]);
    status = struct_cocomcyk_expandct(r, msa, spair, &newctlist[s], G, minloop, verbose, errbuf);
    if (status != eslOK) goto ERROR;     
  }

#if LASTFOLD
   new = nct + 1;
  // Do one more folding in which we force unpaired all that has been paired before
  ESL_ALLOC(newctlist[nct], sizeof(int) * (L+1));
  esl_vec_ISet(newctlist[nct], L+1, 0);
  for (s = 0; s < nct; s ++) {
    newct = newctlist[s];
    
    for (i = 1; i <= L; i ++) {
      if (newct[i] > 0) newctlist[nct][i] = -newct[i];
      if (newct[i] < 0) newctlist[nct][i] =  newct[i];
    }
  }
  status = struct_cocomcyk_expandct(r, msa, spair, &newctlist[nct], G, minloop, verbose, errbuf);
  if (status != eslOK) goto ERROR;
#else
  new = nct;
#endif
  
  // for the extra structures, break in individual helices
  status = ctlist_break_in_helices(&new, &newctlist, L, verbose);
  if (status != eslOK) goto ERROR;     
  
  // all substructures are tested for including at least 1 cov,
  // if they don't they have to be compatible with the major nested structure
  // otherwise they get removed.
  status = ctlist_helices_select(nct, ctlist, &new, &newctlist, L, verbose);
  if (status != eslOK) goto ERROR;     
  
  // if the primary structure (s==0), replace the 'SS_cons' GC line with the new ss
  // all the other substructures add as GC tag as
  //
  // #GC SS_cons_1 
  // #GC SS_cons_2
  status = r2r_Overwrite_SS_cons(msa, new, newctlist, verbose);
  if (status != eslOK) goto ERROR;     
  
  if ((fp = fopen(r2rfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to open r2rfile %s", r2rfile);
  esl_msafile_Write(fp, msa, eslMSAFILE_PFAM);
  fclose(fp);

  *ret_nct    = new;
  *ret_ctlist = newctlist;
  
  free(ss);
  for (s = 0; s < nct; s ++) free(ctlist[s]);
  free(ctlist);
  return eslOK;

 ERROR:
  if (ss) free(ss);
  for (s = 0; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  free(ctlist);
  return status;
}

static int
struct_cocomcyk_expandct( ESL_RANDOMNESS *r, ESL_MSA *msa, SPAIR *spair, int **ret_ct, enum grammar_e G, int minloop, int verbose, char *errbuf)
{
  char    *rfline = NULL;
  char    *newss  = NULL;
  ESL_SQ  *sq     = NULL;
  int     *cct    = NULL;
  int     *ct     = *ret_ct;
  SCVAL    sc;
  float    idthresh = 0.0;
  int      i;
  int      status;

  /* create an RF sequence */
  ESL_ALLOC(rfline, sizeof(char) * (msa->alen+1));
  esl_msa_ReasonableRF(msa, idthresh, TRUE, rfline);
  if (verbose) printf("\nrfline:\n%s\nss_cons\n%s\n", rfline, msa->ss_cons);
  
  sq = esl_sq_CreateFrom(msa->name, rfline, msa->desc, msa->acc, msa->ss_cons); 
  if (sq == NULL) ESL_XFAIL(eslFAIL, errbuf, "failed to create RF sequence");
  status = esl_sq_Digitize((const ESL_ALPHABET *)msa->abc, sq);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to digitize RF sequence");

  // removes covariations that correspond to gap-gap in the RF sequence
  // it also removes covariation that are closer than minloop in RF 
  ct_remove_inconsistencies(sq, ct, minloop, verbose);

  /* calculate the convariance-constraint CYK structure using a probabilistic grammar */
  status = COCOCYK(r, G, sq, spair, ct, &cct, &sc, errbuf, verbose);
  if (status != eslOK) goto ERROR;

  if (1||verbose) {
    ESL_ALLOC(newss, sizeof(char) * (msa->alen+1));
    esl_ct2wuss(cct, msa->alen, newss);
    printf("coco-cyk score = %f\n%s\n", sc, newss);
  }

  if (cct) {
    free(ct); ct = NULL;
    *ret_ct = cct;
  }
  else 
    *ret_ct = ct;
    
  if (rfline) free(rfline);
  if (newss) free(newss);
  esl_sq_Destroy(sq);
  return eslOK;

 ERROR:
  if (sq) esl_sq_Destroy(sq);
  if (rfline) free(rfline);
  if (newss) free(newss);
  if (cct) free(cct);
  return status;
}

static int
struct_write_ss(FILE *fp, int blqsize, int nss, char **sslist)
{
  char  *buf;
  char  *ss;
  char **tag = NULL;
  int   sslen;
  int   tagsize = 0;
  int   n;
  int   s;
  int   status;

  // alocate
  ESL_ALLOC(tag, sizeof(char *) * nss);
  ESL_ALLOC(buf, sizeof(char) * (blqsize+1));
  buf[blqsize] = '\0';

  for (s = 0; s < nss; s ++) {
    if (s == 0) esl_sprintf(&tag[s], "SS_cons");
    else        esl_sprintf(&tag[s], "SS_cons_%d", s);
    tagsize = ESL_MAX(tagsize, strlen(tag[s]));
  }

  // initialize
  n     = 0;
  sslen = strlen(sslist[0]);
  while (n < sslen) {
      
    for (s = 0; s < nss; s ++) {
      ss = sslist[s];
      
      if (blqsize < 0) 
	fprintf(fp, "# %-*s %s\n", tagsize, tag[s], ss);
      else {
	strncpy(buf, ss+n, blqsize);
	fprintf(fp, "# %-*s %s\n", tagsize, tag[s], buf);
      }      
    }
    
    fprintf(fp, "#\n");
    n += (blqsize > 0)? blqsize : sslen;
  }

  for (s = 0; s < nss; s ++) free(tag[s]);
  free(tag);
  if (buf) free(buf);
  return eslOK;

 ERROR:
  for (s = 0; s < nss; s ++) if (tag[s]) free(tag[s]);
  if (tag) free(tag);
  if (buf) free(buf);
  return status;
}

// removes covariations that correspond to gap-gap in the RF sequence
// it also removes covariation that are closer than minloop in RF 
static int
ct_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop, int verbose)
{
  int L = sq->n;
  int n;
  int ipair;
  int i;
  int x;

  // remove covariation that correspond to gap-gap in the RF sequence
  for (i = 1; i <= L; i++) {
    ipair = ct[i];
    if (ipair > 0 && sq->dsq[i] >= NB && sq->dsq[ipair] >= NB) { // remove this covariation
      ct[i]     = 0;
      ct[ipair] = 0;
    }
  }

  // remove covariation that are closer than minloop in RF 
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

// split in helices
static int
ct_split_helices(int *ct, int L, int *ret_nct, int ***ret_ctlist, int verbose)
{
  ESL_STACK  *pda    = NULL;         // stack for secondary structure 
  int       **ctlist = *ret_ctlist;
  int         nct    = *ret_nct;
  int         nsingle_max = 3;       // break stems when there is more than nsigle_max unpaired residues
  int         nsingle;
  int         nfaces;                // number of faces in a structure 
  int         minface;               // max depth of faces in a structure 
  int         npairs = 0;            // total number of basepairs 
  int         npairs_reached = 0;    // number of basepairs found so far 
  int         found_partner;         // true if we've found left partner of a given base in stack pda */
  int         i, j;
  int         s;
  int         status;

  /* total number of basepairs */
  for (j = 1; j <= L; j ++) { if (ct[j] > 0 && j < ct[j]) npairs ++; }
  
  /* Initialization */
  if ((pda  = esl_stack_ICreate()) == NULL) goto FINISH;

  for (j = 1; j <= L; j++)
    {
      if (ct[j] == 0) {   // unpaired 
	if (esl_stack_IPush(pda, j) != eslOK) goto FINISH;
      }
      else if (ct[j] > j) // left side of a bp: push j
	{
	  if (esl_stack_IPush(pda, j) != eslOK) goto FINISH;
	}
      else // right side of a bp; main routine: fing the left partner 
	{
	  found_partner = FALSE;
	  /* Pop back until we find the left partner of j;
	   * In case this is not a nested structure, finding
	   * the left partner of j will require to put bases 
	   * aside into stack auxpk.
	   *
	   * After we find the left partner of j,
	   * store single stranded residues in auxss;
	   * keep track of #faces and the maximum face depth.
	   */
	  nsingle = 0;
	  nfaces  = 0;
	  minface = -1;
	  
	  while (esl_stack_ObjectCount(pda)) 
	    {
	      if (esl_stack_IPop(pda, &i) != eslOK) goto FINISH;
	      
	      if (i < 0) 		/* a face counter */
		{
		  nfaces++;
		  if (i < minface) minface = i;
		}
	      
	      else if (ct[i] == j)  /* we found the i,j pair. */
		{
		  found_partner = TRUE;
		  npairs_reached ++;	
		  
		  /* Now we know i,j pair; and we know how many faces are
		   * above them; and we know the max depth of those faces.
		   * That's enough to label the pair in WUSS notation.
		   * if nfaces == 0, minface is -1; <> a closing bp of a hairpin.
		   * if nfaces == 1, inherit minface, we're continuing a stem.
		   * if nfaces >  1, bump minface in depth; we're closing a bifurc.
		   */
		  if (nfaces > 1) minface--;
		  if (esl_stack_IPush(pda, minface) != eslOK) goto FINISH;
		  //if (verbose) printf("found pair %d %d %d %d %d\n", i, j, nfaces, minface, nsingle);
		  
		  // a new substructure
		  if ( (nfaces == 0)               ||            // a hairpin loop 
		       (nfaces >  1)               ||            // a multiloop 
		       (nfaces == 1 && nsingle > nsingle_max)  ) // break long stems if they have more than > nsingle_max single stranded residudes
		    {	       
		      nct ++;
		      ESL_REALLOC(ctlist,      sizeof(int *) * nct);
		      ESL_ALLOC(ctlist[nct-1], sizeof(int)   * (L+1));
		      esl_vec_ISet(ctlist[nct-1], L+1, 0);
		  }
		    
		  ctlist[nct-1][i] = j;
		  ctlist[nct-1][j] = i;		  
		  break;
		}
	      else if (ct[i] == 0) 
	      {
		nsingle ++;
	      }
	      else /* ct[i]>0, != j: i is paired, but not to j: pseudoknot! */
		{
		esl_stack_Destroy(pda); 
		ESL_EXCEPTION(eslEINVAL, "should not find pseudonots here");
	      }

	    }
	  if (!found_partner) {
	    esl_stack_Destroy(pda); 
	    ESL_EXCEPTION(eslEINVAL, "Cannot find left partner (%d) of base %d. Likely a triplet", ct[j], j);
	  }
	} /* finished finding the left partner of j */
      
    }

  *ret_nct    = nct;
  *ret_ctlist = ctlist;
      
  esl_stack_Destroy(pda);
  return eslOK;
  
 ERROR:
 FINISH:
  if (npairs != npairs_reached) 		  
    ESL_EXCEPTION(eslFAIL, "found %d out of %d pairs.", npairs_reached, npairs);
  if (pda) esl_stack_Destroy(pda);
  return status;
}

static int
ctlist_break_in_helices(int *ret_nct, int ***ret_ctlist, int L, int verbose)
{
  int   **ctlist = *ret_ctlist;
  int   **ctnew  = NULL;
  char   *ss     = NULL;
  int     nct    = *ret_nct;
  int     new    = 1;
  int     prv;
  int     add;
  int     s;
  int     s1;
  int     status;

  // modify only the additional structures
  if (verbose) ESL_ALLOC(ss, sizeof(char)  * (L+1));
  
  // main structure (s=0) remains intact
  ESL_ALLOC(ctnew,        sizeof(int *) * new);
  ESL_ALLOC(ctnew[new-1], sizeof(int)   * (L+1));
  esl_vec_ICopy(ctlist[new-1], L+1, ctnew[new-1]);
  
  // modify only the additional structures (s>0)
  for (s = 1; s < nct; s ++) {
    prv = new;
    if (verbose) {
      esl_ct2wuss(ctlist[s], L, ss);
      printf("break_in_helices: pknot %d\n%s\n", s, ss);
    }
    status = ct_split_helices(ctlist[s], L, &new, &ctnew, verbose);
    if (status != eslOK) goto ERROR;
    
    if (verbose) {
     add = new - prv;
     printf("pseudoknot breaks in %d structures. total structures so far %d prv %d\n", add, new, prv);
      for (s1 = prv; s1 < new; s1 ++) {
	esl_ct2wuss(ctnew[s1], L, ss);
	printf("broken pknot %d\n%s\n", s1, ss);
      }
    }
    free(ctlist[s]); ctlist[s] = NULL;
  }

  *ret_nct    = new;
  *ret_ctlist = ctnew;
  
  free(ctlist);
  if (ss) free(ss);
  
  return eslOK;
  
 ERROR:
  for (s = 0; s < new; s ++) if (ctnew[s])  free(ctnew[s]);
  for (s = 0; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  if (ss)     free(ss);
  return status;
}

// ctlist1 is the list of original covariations split into a main nestes structure and others
//         that cannot explained by it.
// ctlist  is the current list of the structures.
//
static int
ctlist_helices_select(int nctcov, int **ctlistcov, int *ret_nct, int ***ret_ctlist, int L, int verbose)
{
  char   *ssmain   = NULL;
  char   *ss       = NULL;
  int   **ctlist   = *ret_ctlist;
  int   **ctnew    = NULL;
  int    *ctcum    = NULL;
  int    *ctmain   = ctlist[0];
  int    *useme    = NULL;
  int     nct      = *ret_nct;
  double  incompfrac = INCOMPFRAC;
  int     minhelix   = MINHELIX;
  int     new = 1;
  int     hascov;
  int     iscompatible;
  int     isunique;
  int     idx;
  int     ntot, nincomp, ndup;
  int     s, s1;
  int     i;
  int     status;

  // allocate
  ESL_ALLOC(useme, sizeof(int) * nct);
  ESL_ALLOC(ctcum, sizeof(int) * (L+1));
  esl_vec_ISet(useme,   nct, FALSE);
  esl_vec_ICopy(ctmain, L+1, ctcum);  // start the cumulative ct with the main nested structure
  
  useme[0] = TRUE;      // we maintain untouced the main nested structure
  
  if (verbose) {
    ESL_ALLOC(ssmain, sizeof(char)  * (L+1));
    ESL_ALLOC(ss,     sizeof(char)  * (L+1));
    esl_ct2wuss(ctlist[0], L, ssmain);
  }
  
  for (s = 1; s < nct; s ++) {
    hascov       = FALSE;
    iscompatible = TRUE;
    isunique     = TRUE;
    
    // check for covariations
    for (i = 1; i <= L; i ++)
      if (ctlist[s][i] > 0) {
	for (s1 = 0; s1 < nctcov; s1 ++)
	  if (ctlistcov[s1][i] > 0) { hascov = TRUE; break; }
      }
    
    // check for compatiblity with main structure
    // or for a helix with fewer than minhelix bpairs
    nincomp = 0;
    ntot    = 0;
    for (i = 1; i <= L; i ++) {
      if (ctlist[s][i] > 0)                  ntot ++;
      if (ctlist[s][i] > 0 && ctmain[i] > 0) nincomp ++;
    }
    if (nincomp > incompfrac * ntot || (hascov == FALSE && ntot < 2*minhelix)) iscompatible = FALSE;

    // a final check for duplications
    ndup = 0;
    for (i = 1; i <= L; i ++)
      if (ctlist[s][i] > 0  && ctcum[i] > 0) { ndup ++; }
    if (ndup == ntot) isunique = FALSE;

    // after the check add this ct to the cumulative
    for (i = 1; i <= L; i ++)
      if (ctlist[s][i] > 0) ctcum[i] = ctlist[s][i];
    
    if (verbose) {
      esl_ct2wuss(ctlist[s], L, ss);
      printf("%s\n%s\ntot %d nincompatible %d iscomp %d isunique %d\n", ssmain, ss, ntot, nincomp, iscompatible, isunique);
    }

    // add if it satisfies all conditions
    if ((hascov || iscompatible) && isunique) {
      new ++;
      useme[s] = TRUE;
    }
  }
  
  // Write the final set of structures to ctnew
  ESL_ALLOC(ctnew, sizeof(int *) * new);
  idx = 0;
  for (s = 0; s < nct; s ++) {
    if (useme[s]) {
      ESL_ALLOC(ctnew[idx], sizeof(int) * (L+1));
      esl_vec_ICopy(ctlist[s], L+1, ctnew[idx]);
      idx ++;
    }
  }
  
  *ret_nct    = new;
  *ret_ctlist = ctnew;
  
  for (s = 1; s < nct; s ++) free(ctlist[s]);
  free(ctlist);
  free(ctcum);
  free(useme);
  if (ss)     free(ss);
  if (ssmain) free(ssmain);
  return eslOK;
  
 ERROR:
  for (s = 1; s < new; s ++) if (ctnew[s])  free(ctnew[s]);
  for (s = 1; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  if (ctcum)  free(ctcum);
  if (useme)  free(useme);
  if (ss)     free(ss);
  if (ssmain) free(ssmain);
  return status;
}



static int
ctlist_join(int nct, int **ctlist, int L, int **ret_ct, int verbose)
{
  int *ct = NULL;
  int  status;

  ESL_ALLOC(ct, sizeof(int) * (L+1));
  
  *ret_ct = ct;
  return eslOK;

 ERROR:
  if (ct) free(ct);
  return status;
}


