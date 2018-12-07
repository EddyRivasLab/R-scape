/* cykcov.c */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covariation.h"
#include "covgrammars.h"
#include "cykcov.h"

static int  accept_pair(int i, int j, struct mutual_s *mi, CLIST *clist, int *ct_input, THRESH *thresh);
static int  covariations_not_in_structure(int *ct, struct mutual_s *mi, CLIST *clist, THRESH *thresh, int *ret_ncv_in, int verbose);
static int  dp_recursion(struct mutual_s *mi, CLIST *clist, int *ct_input, GMX *cyk, int minloop, THRESH *thresh, int j, int d, SCVAL *ret_sc,  ESL_STACK *alts,
			   char *errbuf, int verbose);
static int  remove_ct_from_ctinput(int *ct_input, int L, const int *ct);

// Include as many of the significantly covarying pairs in one structure.
// We do this by running a nussinov-type algorithm
//
// We return a ct vector.
// We are going to distinguish three situations:
//
//  ct[i] =  j   i-j are a covaring pair that can be put in the structure with the max number of covs
//  ct[i] = -1   i is covarying but cannot be added to the structure with max number of covs
//  ct[i] =  0;  i is not covarying
//
// If there are covaring pairs that cannot be added to one single nested st, then we
// remove the already accounted for pairs, and build more structures untill all
// convarying pairs are accounted for. This allows to follow up pseudoknots.
//
int
CYKCOV(ESL_RANDOMNESS *r, struct mutual_s *mi, CLIST *clist, int *ret_nct, int ***ret_ctlist,
       int ncvpairs, int minloop, THRESH *thresh, char *errbuf, int verbose) 
{
  int n;
  int status;

  if (ESL_MIN(thresh->sc_bp, thresh->sc_nbp) > mi->maxCOV) {
    ESL_ALLOC(*ret_ctlist, sizeof(int *));
    ESL_ALLOC(*ret_ctlist[0], sizeof(int) * (mi->alen+1));
    esl_vec_ISet(*ret_ctlist[0], mi->alen+1, 0);
    *ret_nct = 1;
    return eslOK;
  }
  
  // run cyk/traceback to get the maximal number of cov pairs that form a nested structure
  if ((status = CYKCOV_Structures(r, mi, clist, ret_nct, ret_ctlist, ncvpairs, minloop, thresh, errbuf, verbose)) != eslOK) goto ERROR;    
  
  return eslOK;
  
 ERROR:
  for (n = 0; n < *ret_nct; n ++) if (*ret_ctlist[n]) free(*ret_ctlist[n]);
  if (*ret_ctlist) free(*ret_ctlist);
  return status;
}

int
CYKCOV_Structures(ESL_RANDOMNESS *rng, struct mutual_s *mi, CLIST *clist, int *ret_nct, int ***ret_ctlist, int ncvpairs, int minloop, THRESH *thresh, char *errbuf, int verbose) 
{
  int    **ctlist   = NULL;
  int     *ct_input = NULL; // ct_input=0 allowed ct_input=-1 disallowed
  int     *ct;
  SCVAL    sc;
  int      L = mi->alen;
  int      ncv_in;               // number of cov pair part of a given nested structure
  int      ncv_left = ncvpairs;  // number of cov left to account for
  int      nct = 0;
  int      i;
  int      n;
  int      s; 
  int      status;

  // all allowed at first
  ESL_ALLOC(ct_input, sizeof(int) * (L+1));
  esl_vec_ISet(ct_input, L+1, 0);
  
  while(ncv_left > 0) {
    if (nct == 0) ESL_ALLOC  (ctlist, sizeof(int *) * (nct+1));
    else          ESL_REALLOC(ctlist, sizeof(int *) * (nct+1));
    ctlist[nct] = NULL;

    if ((status = CYKCOV_Both(rng, mi, clist, ct_input, &sc, &ctlist[nct], minloop, thresh, errbuf, verbose)) != eslOK) goto ERROR;
    if (sc == 0) break; // the covariations left don't form a structure
    
    ct = ctlist[nct];
    
    // now look for residues that covary but not in the structure and assign then ct = -1
    // we don't want them to basepair with anything else either
    if (covariations_not_in_structure(ct, mi, clist, thresh, &ncv_in, verbose) != eslOK) goto ERROR;
    ncv_left -= ncv_in;
    
    // this function modifies ct_input to exclude the pairs in ct by assigning them ct_input = -1
    if (remove_ct_from_ctinput(ct_input, mi->alen, ct) != eslOK) goto ERROR;

    if (1||verbose) 
      printf("structure %d [%d cv pairs] CYKscore = %f at covthres %f %f\n", nct+1, ncv_in, sc, thresh->sc_bp, thresh->sc_nbp);
       
    nct ++;
  }

  // go back go all structures.
  // those residues that are marked as -1 in the final ct_input, and are not
  // paired in that particular structure, mark as -1.
  // The next folding step in ExpandCT will keep those -1 residues unpaired
  for (s = 0; s < nct; s ++) {
    ct = ctlist[s];
    for (i = 1; i <= L; i ++) if (ct[i] == 0 && ct_input[i] == -1) ct[i] = -1;
  }
  
  *ret_nct    = nct;
  *ret_ctlist = ctlist;

  free(ct_input);
  return eslOK;

 ERROR:
  if (ct_input) free(ct_input);
  for (s = 0; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  return status;
}

int
CYKCOV_Both(ESL_RANDOMNESS *rng, struct mutual_s *mi, CLIST *clist, int *ct_input, SCVAL *ret_sc, int **ret_ct, int minloop, THRESH *thresh, char *errbuf, int verbose) 
{
  GMX    *cyk = NULL;  // CYK DP matrix: M x (L x L triangular)   
  int     status;

  cyk = GMX_Create(mi->alen);
  
  /* Fill the cyk matrix */
  if ((status = CYKCOV_Fill(mi, clist, ct_input, cyk, ret_sc, minloop, thresh, errbuf, verbose)) != eslOK) goto ERROR;
  
  /* Report a traceback */
  if ((status = CYKCOV_Traceback(rng, mi, clist, ct_input, cyk, ret_ct, minloop, thresh, errbuf, verbose))  != eslOK) goto ERROR;

  GMX_Destroy(cyk);
  return eslOK;

 ERROR:
  if (cyk) GMX_Destroy(cyk);
  return status;
}

int
CYKCOV_Fill(struct mutual_s *mi, CLIST *clist, int *ct_input, GMX *cyk, SCVAL *ret_sc, int minloop, THRESH *thresh, char *errbuf, int verbose) 
{
  SCVAL  sc;
  int    L = mi->alen;
  int    j, d;
  int    status;
  
  /* nussinov grammar: S --> a S | a S a' S | e */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {
	status = dp_recursion(mi, clist, ct_input, cyk, minloop, thresh, j, d, &(cyk->dp[j][d]), NULL, errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "CYK failed");
	if (verbose) printf("\nCYK %f j=%d d=%d L=%d\n", cyk->dp[j][d], j, d, L); 
     } 
  sc = cyk->dp[L][L];
  if (verbose) printf("CYKscore = %f at thresh %f %f\n", sc, thresh->sc_bp, thresh->sc_nbp);

  if (ret_sc)  *ret_sc  = sc;
  
  return eslOK;

 ERROR:
  return status;
}


int
CYKCOV_Traceback(ESL_RANDOMNESS *rng, struct mutual_s *mi, CLIST *clist, int *ct_input, GMX *cyk, int **ret_ct, int minloop, THRESH *thresh, char *errbuf, int verbose) 
{
  ESL_STACK      *ns = NULL;             /* integer pushdown stack for traceback */
  ESL_STACK      *alts = NULL;           /* stack of alternate equal-scoring tracebacks */
  int            *ct = NULL;             /* the ct vector with who is paired to whom */
  SCVAL           bestsc;                /* max score over possible rules */
  int             L = mi->alen;
  int             nequiv;                /* number of equivalent alternatives for a traceback */
  int             x;                     /* a random choice from nequiv */
  int             r;                     /* index of a rule */
  int             d, d1;                 /* optimum values of d1 iterator */
  int             i,j,k;                 /* seq coords */
  float           tol = 0.001;
  int             status;
  
  /* We're going to do a simple traceback that only
   * remembers who was a base pair, and keeps a ct[]
   * array. 
   */
  ESL_ALLOC(ct, sizeof(int) * (L+1));
  esl_vec_ISet(ct, L+1, 0);
  
  /* is sq score is -infty, nothing to traceback */
  if (cyk->dp[L][L] == -eslINFINITY) {
    if (verbose) printf("no traceback.\n");
    return eslOK;
  }

  /* We implement a "stochastic" traceback, which chooses randomly
   * amongst equal-scoring alternative parse trees. This is particularly
   * essential for working with ambiguous grammars, for which 
   * choosing an arbitrary optimal parse tree by order of evaluation
   * can easily result in infinite loops. To do this, we will keep
   * a stack of the alternate solutions.
   */
  alts = esl_stack_ICreate();

  /* Start an integer stack for traversing the traceback.
   * push w,i,j = G->ntS_idx,1,L to init. 
   */
  ns = esl_stack_ICreate();
  esl_stack_IPush(ns, 1);
  esl_stack_IPush(ns, L);

  while (esl_stack_ObjectCount(ns) != 0)
  {
      esl_stack_IPop(ns, &j);
      esl_stack_IPop(ns, &i);
      d = j-i+1;
      
      status = dp_recursion(mi, clist, ct_input, cyk, minloop, thresh, j, d, &bestsc, alts, errbuf, verbose);
       if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "CYK failed");
      
      /* Some assertions.
       */
      if (fabs(bestsc-cyk->dp[j][d]) > tol) 
	ESL_XFAIL(eslFAIL, errbuf, "CYKCOV_Traceback(): that can't happen either. i=%d j=%d d=%d bestsc %f cyk %f", 
		  j-d+1, j, d, bestsc, cyk->dp[j][d]); 
      
      /* Now we know one or more equiv solutions, and they're in
       * the stack <alts>, which keeps 2 numbers (r, d1) for each
       * solution. Choose one of them at random.
       */
      nequiv = esl_stack_ObjectCount(alts) / 2; /* how many solutions? */
      x = esl_rnd_Roll(rng, nequiv);            /* uniformly, 0.nequiv-1 */
      esl_stack_DiscardTopN(alts, x*2);         /* dig down to choice */
      esl_stack_IPop(alts, &d1);
      esl_stack_IPop(alts, &r);

      /* Now we know a best rule; figure out where we came from,
       * and push that info onto the <ns> stack.
       */
      if (verbose) {
        printf("-----------------------------------\n"); 
        printf("i=%d j=%d d=%d d1=%d\n", j-d+1, j, d, d1);
	printf("tracing %f\n", bestsc);
        printf("   rule(%d)\n", r);
      }

      i = j - d  + 1;
      k = i + d1 - 1;

      if (r == 0) {
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, j);
      }
      else if (r == 1) {
	esl_stack_IPush(ns, i+1);
	esl_stack_IPush(ns, k-1);
	esl_stack_IPush(ns, k+1);
	esl_stack_IPush(ns, j);
	ct[i] = k;
	ct[k] = i;
      }
      else if (r == 2) {
      }
      else 	
	ESL_XFAIL(eslFAIL, errbuf, "rule %d disallowed. Max number is 2", r);
  }

  *ret_ct = ct;

  esl_stack_Destroy(ns);
  esl_stack_Destroy(alts);
  return eslOK;
  
 ERROR:
  if (ns)   esl_stack_Destroy(ns);
  if (alts) esl_stack_Destroy(alts);
  if (ct)   free(ct);
  return status;
}



static int
dp_recursion(struct mutual_s *mi, CLIST *clist, int *ct_input, GMX *cyk, int minloop, THRESH *thresh, int j, int d, SCVAL *ret_bestsc,  ESL_STACK *alts, char *errbuf, int verbose)
{
  SCVAL bestsc = -eslINFINITY;
  SCVAL sc;
  int   d1;
  int   r;
  int   i, k;
  
  if (alts) esl_stack_Reuse(alts);
   
  i = j - d + 1;
  /* rule1: a S */
  r  = 0;
  d1 = 0;
  if (d > 0) {
    bestsc = (d > 1)? cyk->dp[j][d-1] : 0.;
    if (alts) {
      esl_stack_IPush(alts, r);
      esl_stack_IPush(alts, d1);
    }
  }

  /* rule2: a S a' S  */
  r = 1;
  for (d1 = minloop; d1 <= d; d1++) {
    k = i + d1 - 1;
    sc = (accept_pair(i, k, mi, clist, ct_input, thresh))? cyk->dp[k-1][d1-2] + cyk->dp[j][d-d1] + mi->COV->mx[i-1][k-1] : -eslINFINITY;
    
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }
      
      if (alts) {
	esl_stack_IPush(alts, r);
	esl_stack_IPush(alts, d1);
      }
    }
  }

  /* rule3: epsilon  */
  r = 2;
  if (d == 0) {
    sc = 0;
    if (sc >= bestsc) {
      if (sc > bestsc) { /* if an outright winner, clear/reinit the stack */
	if (alts) esl_stack_Reuse(alts);
	bestsc = sc;
      }
      
      if (alts) {
	esl_stack_IPush(alts, r);
	esl_stack_IPush(alts, 0);
      }
    }
  }

  *ret_bestsc = bestsc;
  return eslOK;
}

static int
accept_pair(int i, int j, struct mutual_s *mi, CLIST *clist, int *ct_input, THRESH *thresh)
{
  double cov  = mi->COV->mx[i-1][j-1];
  int    isbp = FALSE;
  int    c;

  // is it allowed by ct_input?
  if (ct_input[i] < 0 || ct_input[j] < 0) return FALSE;
  
  for (c = 0; c < clist->ncnt; c++) 
    if (i == clist->cnt[c].i && j == clist->cnt[c].j) isbp = clist->cnt[c].isbp;

  if ( isbp && cov >= thresh->sc_bp)  return TRUE;
  if (!isbp && cov >= thresh->sc_nbp) return TRUE;
  
  return FALSE;
}

// now look for residues that covary but not in the structure and assign then ct = -1
// we don't want them to basepair with anything else either
static int 
covariations_not_in_structure(int *ct, struct mutual_s *mi, CLIST *clist, THRESH *thresh, int *ret_ncv_in, int verbose)
{
  double   cov;
  int      ncv_in = 0;
  int      L = mi->alen;
  int      i, j;
  int      isbp;
  int      c;
  
  for (i = 1; i < L; i ++) {

    if (ct[i] > 0) { if (ct[i] > i) ncv_in ++; continue; } // it is in the ss
 
    for (j = i+1; j <= L; j ++) {
      if (ct[j] > 0) continue; // it is in the ss
      
      isbp = FALSE;
      for (c = 0; c < clist->ncnt; c++) 
	if (i == clist->cnt[c].i && j == clist->cnt[c].j) isbp = clist->cnt[c].isbp; 
      cov = mi->COV->mx[i-1][j-1];
      
      if ( isbp && cov > thresh->sc_bp)  { ct[i] = -1; ct[j] = -1; }
      if (!isbp && cov > thresh->sc_nbp) { ct[i] = -1; ct[j] = -1; }
    }
  }

  *ret_ncv_in = ncv_in;
  
  return eslOK;
}

// this function modifies ct_input so that all pairs in ct get -1
static int 
remove_ct_from_ctinput(int *ct_input, int L, const int *ct)
{
  int i;

  for (i = 1; i <= L; i ++)
    if (ct[i] > 0) ct_input[i] = -1; // already accounted in ct remove from ct_input
     
  return eslOK;
}

