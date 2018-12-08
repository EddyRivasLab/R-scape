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

static int    cykcov_remove_inconsistencies(ESL_SQ *sq, int *ct, int minloop);
static int    is_stacked_pair(int i, int j, int L, int *ct);
static int    is_cannonical_pair(char nti, char ntj);
static int    break_in_helices(int *ret_nct, int ***ret_ctlist, int L);
static int    ct2helices(int *ct, int L, int *ret_idx, int ***ret_ctlist);
static int    helices_keep_if_cov(int nct1, int **ctlist1, int *ret_nct, int ***ret_ctlist, int L);
static int    esl_wuss_join(int nct, int **ctlist, int L, int **ret_ct, int verbose);

int
struct_CYKCOV(struct data_s *data, ESL_MSA *msa, int *ret_nct, int ***ret_cykctlist, int minloop, RANKLIST *ranklist, HITLIST *hitlist, enum grammar_e G, THRESH *thresh)
{
  HITLIST       *cykhitlist = NULL;
  int          **cykctlist  = NULL;
  int           *cykct;
  char          *covtype    = NULL;
  char          *threshtype = NULL;
  SCVAL          sc;
  int            nct;
  int            i;
  int            s;
  int            status;
            
  /* calculate the cykcov ct vector.
   *
   * I run a nussinov-type algorithm that incorporates as many of the significant pairs as possible.
   * These pairs become constrains for the second part of the folding in struct_ExpandCT()
   */
  status = CYKCOV(data->r, data->mi, data->clist, &nct, &cykctlist, hitlist->nhit, minloop, thresh, data->errbuf, data->verbose);
  if (status != eslOK) goto ERROR;

  /* Use the CT from covariation to do a contrain folding */
  status = struct_ExpandCT(data->R2Rcykfile, data->R2Rall, data->r, msa, &nct, &cykctlist, minloop, G, data->verbose, data->errbuf);
  if (status != eslOK) goto ERROR;

  // create a new contact list from the cykct
  CMAP_ReuseCList(data->clist);
  for (s = 0; s < nct; s ++) {
    status = ContacMap_FromCT(data->clist, msa->alen, cykctlist[s], data->clist->mind, data->msamap, NULL);
    if (status != eslOK) goto ERROR;
  }
  
  /* redo the hitlist since the ct has now changed */
  corr_COVTYPEString(&covtype, data->mi->type, data->errbuf);
  cov_THRESHTYPEString(&threshtype, data->thresh->type, NULL);
  status = cov_CreateCYKHitList(data, ranklist, hitlist, &cykhitlist, covtype, threshtype);
  if (status != eslOK) goto ERROR;
  
  for (s = 0; s < nct; s ++) {
    cykct = cykctlist[s];
    
    for (i = 1; i <= msa->alen; i ++) 
      if (cykct[i] > 0 && i < cykct[i]) data->nbpairs_cyk ++;
  }

  if (data->nbpairs_cyk > 0) {
    /* R2R */
    status = r2r_Depict(data->R2Rcykfile, data->R2Rall, msa, nct, cykctlist, cykhitlist, TRUE, TRUE, data->verbose, data->errbuf);
    if (status != eslOK) goto ERROR;
    
    /* DotPlots (pdf,svg) */
    status = struct_DotPlot(data->gnuplot, data->cykdplotfile, msa, nct, cykctlist, data->mi, data->msamap, data->firstpos, data->samplesize, cykhitlist,
			 TRUE,  data->verbose, data->errbuf);
    if (status != eslOK) goto ERROR;
    status = struct_DotPlot(data->gnuplot, data->cykdplotfile, msa, nct, cykctlist, data->mi, data->msamap, data->firstpos, data->samplesize, cykhitlist,
			 FALSE, data->verbose, data->errbuf);
    if (status != eslOK) goto ERROR;
  }

  *ret_nct       = nct;
  *ret_cykctlist = cykctlist;

  cov_FreeHitList(cykhitlist);
  free(covtype);
  free(threshtype);
  return eslOK;
  
 ERROR:
  if (cykhitlist) cov_FreeHitList(cykhitlist);
  if (covtype)    free(covtype);
  if (threshtype) free(threshtype);
  if (cykct)      free(cykct);
  return status;
}

int              
struct_DotPlot(char *gnuplot, char *dplotfile, ESL_MSA *msa, int nct, int **ctlist, struct mutual_s *mi, int *msamap, int firstpos, SAMPLESIZE samplesize, HITLIST *hitlist, 
	    int dosvg, int verbose, char *errbuf)
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


int
struct_ExpandCT(char *r2rfile, int r2rall, ESL_RANDOMNESS *r, ESL_MSA *msa, int *ret_nct, int ***ret_ctlist, int minloop, enum grammar_e G, int verbose, char *errbuf)
{
  FILE  *fp = NULL;
  char  *ss = NULL;
  char  *tag;
  int    nct       = *ret_nct;
  int  **ctlist    = *ret_ctlist;
  int    new;
  int  **newctlist = NULL;
  int    tagidx;
  int    L = msa->alen;
  int    keep;
  int    s, s1;
  int    i;
  int    status;
  
  if (r2rfile == NULL) return eslOK;

  ESL_ALLOC(ss,        sizeof(char)  * (L+1));
  ESL_ALLOC(newctlist, sizeof(int *) * nct);
  new = nct;
  
  for (s = 0; s < nct; s ++) {
    // covariance-constraint CYK using a probabilistic grammar
    ESL_ALLOC(newctlist[s], sizeof(int) * (L+1));
    esl_vec_ICopy(ctlist[s], L+1, newctlist[s]);
    status = struct_ExpandCT_CCCYK(r, msa, &newctlist[s], G, minloop, verbose, errbuf);
    if (status != eslOK) goto ERROR;     
  }
  
  // for the extra structures, break in individual helices
  // maintain only the helices with covariation support
  status = break_in_helices(&new, &newctlist, L);
  if (status != eslOK) goto ERROR;     
  
  // all substructures are tested for including at least 1 cov,
  // otherwise they get removed.
  status = helices_keep_if_cov(nct, ctlist, &new, &newctlist, L);
  if (status != eslOK) goto ERROR;     
  
  // if the primary structure (s==0), replace the 'SS_cons' GC line with the new ss
  //
  // all the other substructures add as GC tag as
  //
  // #GC SS_cons_1 
  // #GC SS_cons_2
  for (s = 0; s < new; s ++) {
    esl_ct2simplewuss(newctlist[s], L, ss);
    if (1||verbose) printf("^^%d cyk structure\n%s\n", s, ss);
    
    if (s == 0) strcpy(msa->ss_cons, ss);
    else {
      esl_sprintf(&tag, "SS_cons_%d", s);
      esl_msa_AppendGC(msa, tag, ss);
    }
  }

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

int
struct_ExpandCT_CCCYK( ESL_RANDOMNESS *r, ESL_MSA *msa, int **ret_ct, enum grammar_e G, int minloop, int verbose, char *errbuf)
{
  char    *rfline = NULL;
  char    *newss = NULL;
  ESL_SQ  *sq = NULL;
  int     *ct = *ret_ct;
  int     *cct = NULL;
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
  
  cykcov_remove_inconsistencies(sq, ct, minloop);

  /* calculate the convariance-constraint CYK structure using a probabilistic grammar */
  status = COCOCYK(r, G, sq, ct, &cct, &sc, errbuf, verbose);
  if (status != eslOK) goto ERROR;

  if (verbose) {
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

// take a ct vector possibly with pseudoknots, and separate
// into one ct without pseudoknots and additional ct's one with
// each of the pseudoknots
int
struct_SplitCT(int *ct, int L, int *ret_nct, int ***ret_ctlist, int verbose)
{
  int  **ctlist = NULL;
  char  *ss1 = NULL;
  char  *ss2 = NULL;
  int    nct = 1;
  int    use;
  int    n;
  int    c;
  int    s;
  int    status;

  ESL_ALLOC(ss1, sizeof(char)  * (L+1));
  ESL_ALLOC(ss2, sizeof(char)  * (L+1));

  esl_ct2wuss(ct, L, ss1);
  if (verbose) printf("%s\n", ss1);
  
  // the nested structure
  ESL_ALLOC(ctlist,    sizeof(int *) * nct);
  ESL_ALLOC(ctlist[0], sizeof(int)   * (L+1));
  for (n = 0; n < L; n ++)
    {
      if (isalpha(ss1[n])) ss2[n] = '.';
      else                 ss2[n] = ss1[n];
    }
  ss2[L] = '\0';
  if (1||verbose) printf("main structure\n%s\n", ss2);
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
      if (1||verbose) printf("pseudoknot %d\n%s\n", nct-1, ss2);
    }    
  }

  break_in_helices(&nct, &ctlist, L);
  
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

static int
is_stacked_pair(int i, int j, int L, int *ct) 
{
  int is_stacked = FALSE;

  if (ct[i] > 0 || ct[j] > 0) return FALSE; // both have to be unpaired

  if (ct[i+1] == j-1 && ct[j-1] == i+1                  ) is_stacked = TRUE;
  if (ct[i-1] == j+1 && ct[j+1] == i-1 && i > 1 && j < L) is_stacked = TRUE;
	
  return is_stacked;
}

static int
is_cannonical_pair(char nti, char ntj) 
{
  int is_cannonical = FALSE;

  if (nti == 'A' && ntj == 'U') is_cannonical = TRUE;
  if (nti == 'C' && ntj == 'G') is_cannonical = TRUE;
  if (nti == 'G' && ntj == 'C') is_cannonical = TRUE;
  if (nti == 'G' && ntj == 'U') is_cannonical = TRUE;
  if (nti == 'G' && ntj == 'Y') is_cannonical = TRUE;
  if (nti == 'U' && ntj == 'A') is_cannonical = TRUE;
  if (nti == 'U' && ntj == 'G') is_cannonical = TRUE;
  if (nti == 'U' && ntj == 'R') is_cannonical = TRUE;
  if (nti == 'R' && ntj == 'Y') is_cannonical = TRUE;
  if (nti == 'Y' && ntj == 'R') is_cannonical = TRUE;

  return is_cannonical;
}


static int
break_in_helices(int *ret_nct, int ***ret_ctlist, int L)
{
  int **ctlist = *ret_ctlist;
  int **ctnew  = NULL;
  int   nct    = *ret_nct;
  int   new    = 1;
  int   s;
  int   status;
  
  // modify only the additional structures
  ESL_ALLOC(ctnew,    sizeof(int *) * new);
  ESL_ALLOC(ctnew[0], sizeof(int)    * (L+1));
  esl_vec_ICopy(ctlist[0], L+1, ctnew[0]);
   
  for (s = 1; s < nct; s ++) 
    ct2helices(ctlist[s], L, &new, &ctnew);

  *ret_nct    = new;
  *ret_ctlist = ctnew;
  
  for (s = 1; s < nct; s ++) free(ctlist[s]);
  free(ctlist);
  
  return eslOK;
  
 ERROR:
  for (s = 1; s < new; s ++) if (ctnew[s])  free(ctnew[s]);
  for (s = 1; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  return status;
}

static int
ct2helices(int *ct, int L, int *ret_nct, int ***ret_ctlist)
{
  ESL_STACK  *pda    = NULL;         /* stack for secondary structure */
  int       **ctlist = *ret_ctlist;
  int         nct    = *ret_nct;
  int         nfaces;                /* number of faces in a cWW structure */
  int         minface;               /* max depth of faces in a cWW structure */
  int         npairs = 0;            /* total number of basepairs */
  int         npairs_reached = 0;    /* number of basepairs found so far */
  int         found_partner;         /* true if we've found left partner of a given base in stack pda */
  int         i, j;
  int         s;
  int         status;

  /* total number of basepairs */
  for (j = 1; j <= L; j ++) { if (ct[j] > 0 && j < ct[j]) npairs ++; }
  
  /* Initialization */
  if ((pda  = esl_stack_ICreate()) == NULL) goto FINISH;

  for (j = 1; j <= L; j++)
    {
      if (ct[j] == 0) {   /* unparied, do nothing */
      }
      else if (ct[j] > j) /* left side of a bp: push j. */
	{
	  if (esl_stack_IPush(pda, j) != eslOK) goto FINISH;
	}
      else                /* right side of a bp; main routine: fingh the left partner */
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
		   * if nfaces > 1, bump minface in depth; we're closing a bifurc.
		   */
		  if (nfaces > 1) minface--;
		  if (esl_stack_IPush(pda, minface) != eslOK) goto FINISH;

		  if (nfaces == 0) {
		    nct ++;
		    ESL_REALLOC(ctlist, sizeof(int *) * nct);
		    ESL_ALLOC(ctlist[nct-1], sizeof(int) * (L+1));
		    esl_vec_ISet(ctlist[nct-1], L+1, 0);
		  }
		    
		  ctlist[nct-1][i] = j;
		  ctlist[nct-1][j] = i;		  
		  break;
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
      
  return eslOK;
  
 ERROR:
 FINISH:
  if (npairs != npairs_reached) 		  
    ESL_EXCEPTION(eslFAIL, "found %d out of %d pairs.", npairs_reached, npairs);
  if (pda   != NULL) esl_stack_Destroy(pda);
  return status;
}

static int
helices_keep_if_cov(int nct1, int **ctlist1, int *ret_nct, int ***ret_ctlist, int L)
{
  int **ctlist = *ret_ctlist;
  int **ctnew  = NULL;
  int   nct    = *ret_nct;
  int   new    = 1;
  int   keep;
  int   s, s1;
  int   i;
  int   status;
  
  // modify only the additional structures
  ESL_ALLOC(ctnew,    sizeof(int *) * new);
  ESL_ALLOC(ctnew[0], sizeof(int)   * (L+1));
  esl_vec_ICopy(ctlist[0], L+1, ctnew[0]);
  
  for (s = 1; s < nct; s ++) {
    keep = FALSE;
    for (i = 1; i <= L; i ++)
      if (ctlist[s][i] > 0) {
	for (s1 = 0; s1 < nct1; s1 ++)
	  if (ctlist1[s1][i] > 0) { keep = TRUE; break; }
      }
    
    if (keep) {
      new ++;
      ESL_REALLOC(ctnew,      sizeof(int *) * new);
      ESL_ALLOC(ctnew[new-1], sizeof(int)   * (L+1));
      esl_vec_ICopy(ctlist[s], L+1, ctnew[new-1]);
    }
  }

  *ret_nct    = new;
  *ret_ctlist = ctnew;
  
  for (s = 1; s < nct; s ++) free(ctlist[s]);
  free(ctlist);
  
  return eslOK;
  
 ERROR:
  for (s = 1; s < new; s ++) if (ctnew[s])  free(ctnew[s]);
  for (s = 1; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  return status;
}

static int
esl_wuss_join(int nct, int **ctlist, int L, int **ret_ct, int verbose)
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

