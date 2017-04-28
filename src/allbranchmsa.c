/*  allbranchmsa - build a tree from an alignment
 * Contents:
 *
 * ER, Fri Apr 28 09:42:16 EDT 2017 [Harvard] 
 * SVN $Id:$
 */

#include "rscape_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>
	
#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "allbranchmsa.h"
#include "plot.h"

static int     hit_sorted_by_sc(const void *vh1, const void *vh2);
static double  score(int c1, int c2, int dim, int **val);

/* Plots where the changes happen in the tree branches */
int
AllBranchMSA_Plot(char *plotfile, char *gnuplot, ESL_TREE *T, int *msamap, ESL_MSA *allmsa, int *ct, char *errbuf, int verbose)
{
  FILE           *plotfp = NULL;
  char           *psfile = NULL;
  FILE           *pipe = NULL;
  int           **mutb = NULL;
  ESL_DSQ        *ax;
  ESL_DSQ        *axl;
  ESL_DSQ        *axr;
  ESL_DSQ         cc;
  ESL_DSQ         cl, cr;
  CLIST          *clist = NULL;
  int             alloc_nhit = 5;
  int             nhit;
  int             nhit_cutoff;
  int             h;
  int             L = allmsa->alen;
  int             nnodes;
  int             mutbl, mutbr;
  int             v, vl, vr;
  int             c;
  int             ci, cj;
  int             pos;
  int             posi, posj;
  int             dim;
  int             i;
  double          sc;
  int           **pmut = NULL;
  int             K = allmsa->abc->K;
  int             k, k1;
  int             bidx;
  int             tt = 0;
  int             tc = 0;
  int             status;

  if (plotfile == NULL) return eslOK;
  if ((plotfp = fopen(plotfile, "w")) == NULL) esl_fatal("Failed to open plotfile %s", plotfile);

  /* allmsa include the ancestral sequences
   * allmsa->ax[0,N-1] = msa->ax[0,N-1]
   * allmsa->ax[N+n] = ancestral sequence at node 0 <= n < N-1 (n=0 is the root)
   */
  if (T == NULL) nnodes = 1;
  else           nnodes = (T->N > 1)? T->N-1 : T->N;
  dim = nnodes + T->N;
  if (allmsa->nseq != dim) esl_fatal("Tree_AllBranchMSAPlot() not a all-branch msa");

  if (ct) {
  }
  
  pipe = popen(gnuplot, "w");
  fprintf(pipe, "set terminal postscript color 14\n");
  esl_sprintf(&psfile, "%s.ps", plotfile);
  fprintf(pipe, "set output '%s'\n", psfile);
  /* matlab's 'jet' colormap scale */
  //fprintf(pipe, "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0)\n");
  fprintf(pipe, "set palette defined (0 'white', 1 'black')\n");
 
  /* solid squares */
  fprintf(pipe, "set style line 101  lt 1 lc rgb 'gray'      pt 5 lw 2 ps 0.3\n");
  fprintf(pipe, "set style line 102  lt 1 lc rgb 'brown'     pt 5 lw 2 ps 0.7\n");
  fprintf(pipe, "set style line 103  lt 1 lc rgb 'cyan'      pt 5 lw 2 ps 0.7\n");
  fprintf(pipe, "set style line 104  lt 1 lc rgb 'red'       pt 5 lw 2 ps 0.7\n");
  fprintf(pipe, "set style line 105  lt 1 lc rgb 'orange'    pt 5 lw 2 ps 0.7\n");
  fprintf(pipe, "set style line 106  lt 1 lc rgb 'turquoise' pt 5 lw 2 ps 0.7\n");
  fprintf(pipe, "set style line 107  lt 1 lc rgb 'black'     pt 5 lw 2 ps 0.7\n");
  fprintf(pipe, "set style line 108  lt 1 lc rgb 'black'     pt 5 lw 2 ps 0.5\n");
  
  fprintf(pipe, "set ylabel 'Tree Branch'\n");
  fprintf(pipe, "set xlabel 'Alignment position'\n");
  fprintf(pipe, "set xrange [%d:%d]\n", 1, L);
  fprintf(pipe, "set yrange [%d:%d]\n", 1, dim);
  fprintf(pipe, "unset title\n");
  fprintf(pipe, "plot '-' u 1:2:3:3 title '' with point ls 107 palette\n");

  ESL_ALLOC(mutb,    sizeof(int *)*(L+1));
  ESL_ALLOC(mutb[0], sizeof(int)  *(L+1)*dim);
  for (c = 1; c <= L; c ++) {
    mutb[c] = mutb[0]  + c*dim;
    for (i = 0; i < dim; i ++) mutb[c][i] = 0;
  }

  bidx = 1;
  ESL_ALLOC(pmut,    sizeof(int *)*K);
  ESL_ALLOC(pmut[0], sizeof(int)  *K*K);
  for (k = 0; k < K; k ++) {
    pmut[k] = pmut[0]  + k*K;
    for (k1 = 0; k1 < K; k1 ++) pmut[k][k1] = 0;
  }
  for (v = 0; v < nnodes; v ++) {
    vl = T->left[v];
    vr = T->right[v];
    ax  = allmsa->ax[T->N+v];
    axl = (vl >= 0)? allmsa->ax[T->N+vl] : allmsa->ax[-vl-1];
    axr = (vr >= 0)? allmsa->ax[T->N+vr] : allmsa->ax[-vr-1];

    // now test
    for (c = 1; c <= L; c ++) {
      cc = ax[c];
      cl = axl[c];
      cr = axr[c];
     
      pos = (msamap)? msamap[c-1]+1 : c;
      
      mutbl = FALSE;
      if (esl_abc_XIsCanonical(allmsa->abc, cc) && esl_abc_XIsCanonical(allmsa->abc, cl)) {
	if (cc != cl) {
	  if (ct && ct[c] > 0) pmut[cc][cl] ++;
	  mutbl = TRUE;
	}
      }
      mutb[c][bidx] = mutbl;
      
      mutbr = FALSE;
      if (esl_abc_XIsCanonical(allmsa->abc, cc) && esl_abc_XIsCanonical(allmsa->abc, cr)) {
	if (cc != cr) {
	  if (ct && ct[c] > 0) pmut[cc][cr] ++;
	  mutbr = TRUE;
	}
      }
      mutb[c][bidx+1] = mutbr;

      if (ct && ct[c] > 0 ) {
	fprintf(plotfp, "%d %d %d\n", c, bidx,   mutbl); 	
	fprintf(pipe,   "%d %d %d\n", c, bidx,   mutbl);
	fprintf(plotfp, "%d %d %d\n", c, bidx+1, mutbr);
	fprintf(pipe,   "%d %d %d\n", c, bidx+1, mutbr);
      }
    }
    bidx += 2;
  }
  fprintf(pipe, "e\n");


  ESL_ALLOC(clist, sizeof(CLIST));
  clist->hit    = NULL;
  nhit = alloc_nhit;
  h    = 0;
  ESL_ALLOC(clist->hit,    sizeof(CHIT)   * nhit);
  ESL_ALLOC(clist->srthit, sizeof(CHIT *) * nhit);
  clist->srthit[0] = clist->hit;

  for (ci = 1; ci <= L; ci ++) {
    posi = (msamap)? msamap[ci-1]+1 : ci;    
    for (cj = ci+1; cj <= L; cj ++) {
      posj = (msamap)? msamap[cj-1]+1 : cj;

      if (h == nhit - 1) {
	nhit += alloc_nhit;
	ESL_REALLOC(clist->hit,    sizeof(CHIT)   * nhit);
	ESL_REALLOC(clist->srthit, sizeof(CHIT *) * nhit);
     }
      /* assign */
      clist->hit[h].i    = ci;
      clist->hit[h].j    = cj;
      clist->hit[h].posi = posi;
      clist->hit[h].posj = posj;
      clist->hit[h].cont = FALSE;
      if (ct && ct[ci] == cj) clist->hit[h].cont = TRUE;
      clist->hit[h].sc   = score(ci, cj, dim, mutb);
      h ++;
    }
  }
  nhit = h;
  clist->nhit = nhit;
  for (h = 0; h < nhit; h++) clist->srthit[h] = clist->hit + h;
  if (nhit > 1) qsort(clist->srthit, nhit, sizeof(CHIT *), hit_sorted_by_sc);

  
  nhit_cutoff = ESL_MIN(nhit, 0.30*L);
  
  fprintf(pipe, "set title 'N = %d'\n", nhit_cutoff);
  fprintf(pipe, "set xlabel 'Alignment position'\n");
  fprintf(pipe, "set ylabel 'Alignment position'\n");     
  fprintf(pipe, "set size square\n");
  fprintf(pipe, "set xrange [%d:%d]\n", 1, msamap[L-1]+1);
  fprintf(pipe, "set yrange [%d:%d]\n", 1, msamap[L-1]+1);
  fprintf(pipe, "plot '-' u 1:2 title '' with point ls 108\n");
  for (h = 0; h < nhit_cutoff; h ++) {
    posi = clist->srthit[h]->posi;
    posj = clist->srthit[h]->posj;

    if (clist->srthit[h]->cont) { fprintf(stdout, "*%d %d %f\n", posi, posj, clist->srthit[h]->sc); tt ++; tc ++; }
      else                      { fprintf(stdout, " %d %d %f\n", posi, posj, clist->srthit[h]->sc); tt ++; }
    fprintf(pipe,   "%d %d %d\n", posi, posj, 1);
    fprintf(pipe,   "%d %d %d\n", posj, posi, 1);
  }
  fprintf(pipe, "e\n");
  printf("%d/%d\n", tc, tt);
  
  fprintf(pipe, "set xlabel 'residue'\n");
  fprintf(pipe, "set ylabel 'residue'\n");     
  fprintf(pipe, "set size square\n");
  fprintf(pipe, "set xrange [-0.1:%f]\n", K-0.9);
  fprintf(pipe, "set yrange [-0.1:%f]\n", K-0.9);
  //fprintf(pipe, "unset yrange\n");
  fprintf(pipe, "plot '-' u 1:2:3:3 title '' with point ls 108 palette\n");
  for (k = 0; k < K; k ++) {
    for (k1 = 0; k1 < K; k1 ++) {
 	fprintf(pipe,   "%d %d %d\n", k, k1, pmut[k][k1]);
      }
  }
  fprintf(pipe, "e\n");
  
  
  fclose(plotfp);
  pclose(pipe);
  
  //plot_file_ps2pdf(psfile);

  free(mutb[0]);
  free(mutb);
  free(pmut[0]);
  free(pmut);
  return eslOK;

 ERROR:
  if (mutb[0]) free(mutb[0]);
  if (mutb)    free(mutb);
  if (pmut[0]) free(pmut[0]);
  if (pmut)    free(pmut);

  return status;
}



/*--- internal functions ---*/

static int
hit_sorted_by_sc(const void *vh1, const void *vh2)
{
  CHIT *h1 = *((CHIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  CHIT *h2 = *((CHIT **) vh2);

  if      (h1->sc < h2->sc) return  1;
  else if (h1->sc > h2->sc) return -1;
  else {
 
    /* report first pair first */
    int dir1 = (h1->i < h2->i ? 1 : -1);
    int dir2 = (h1->j < h2->j ? 1 : -1);
    if (dir1 != dir2) return dir2; // so if dir1 is pos (1), and dir2 is neg (-1), this will return -1, placing h1 before h2;  otherwise, vice versa
    else              return dir1;

  }
}

static double
score(int c1, int c2, int dim, int **val)
{
  double sc = 0.;
  int    val1;
  int    val2;  
  int    i;

  for (i = 0; i < dim; i ++) {
    val1 = val[c1][i];
    val2 = val[c2][i];
    if      (val1 == 1 && val2 == 1) sc += 2.;
    else if (val1 == 0 && val2 == 0) sc += 1.;
    else                             sc -= 0.;
  }
  
  return sc;
}

