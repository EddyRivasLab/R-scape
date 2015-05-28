/* cykcov.c */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hmmer.h"

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
#include "cykcov.h"

static int dp_recursion(struct mutual_s *mi, GMX *cyk, int minloop, int j, int d, SCVAL *ret_sc, char *errbuf, int verbose);

int
CYKCOV(struct mutual_s *mi, GMX **ret_cyk, int **ret_ct, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  int status;
  /* Fill the cyk matrix */
  if ((status = CYKCOV_Fill(mi, ret_cyk, ret_sc, errbuf, verbose)) != eslOK) goto ERROR;
  
  /* Report a traceback */
  if ((status = CYK_Traceback(*ret_cyk, ret_ct, errbuf, verbose))  != eslOK) goto ERROR;
  
  return eslOK;
  
 ERROR:
  return status;
}

int
CYKCOV_Fill(struct mutual_s *mi, GMX **ret_cyk, SCVAL *ret_sc, char *errbuf, int verbose) 
{
  GMX   *cyk;           /* CYK DP matrix: M x (L x L triangular)     */
  SCVAL  sc;
  int    minloop = 5;
  int    L = mi->alen;
  int    j, d, d1;
  int    status;
  
  /* nussinov grammar: S --> S a | a S a' S | e */
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; d++)
      {
	status = dp_recursion(mi, cyk, minloop, j, d, &(cyk->dp[j][d]), errbuf, verbose);
	if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "CYK failed");
     } 
  sc = cyk->dp[L][L];
  
  if (ret_cyk) *ret_cyk = cyk;
  if (ret_sc)  *ret_sc  = sc;
  
  return eslOK;

 ERROR:
  return status;
}


int
CYKCOV_Traceback(struct mutual_s *mi, GMX *cyk, int **ret_ct, char *errbuf, int verbose) 
{
  ESL_STACK      *ns = NULL;             /* integer pushdown stack for traceback */
  ESL_STACK      *alts = NULL;           /* stack of alternate equal-scoring tracebacks */
  ESL_RANDOMNESS *rng = NULL;            /* random numbers for stochastic traceback */
  int            *ct = NULL;             /* the ct vector with who is paired to whom */
  SCVAL           bestsc;                /* max score over possible rules for nonterminal w  */
  int             nequiv;                /* number of equivalent alternatives for a traceback */
  int             x;                     /* a random choice from nequiv */
  int             r;                     /* index of a rule */
  int             d1,d2;                 /* optimum values of d1, d2 iterators */
  int             i,j,d;                 /* seq coords */
  float           tol = 0.001;
  int             status;

  L = cyk->L;

  /* We're going to do a simple traceback that only
   * remembers who was a base pair, and keeps a ct[]
   * array. 
   */
  ESL_ALLOC(ct, sizeof(int *) * (L+1));
  
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
  rng  = esl_randomness_CreateTimeseeded();

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

      /* Some assertions.
       */
      if (fabs(bestsc-cyk->dp[j][d]) > tol) 
	ESL_XFAIL(eslFAIL, errbuf, "CYKCOV_Traceback(): that can't happen either. i=%d j=%d d=%d bestsc %f cyk %f", 
		  j-d+1, j, d, bestsc, cyk->dp[j][d]); 
      
      /* Now we know one or more equiv solutions, and they're in
       * the stack <alts>, which keeps 3 numbers (r, d1, d2) for each
       * solution. Choose one of them at random.
       */
      nequiv = esl_stack_ObjectCount(alts) / 3; /* how many solutions? */
      x = esl_rnd_Roll(rng, nequiv);            /* uniformly, 0.nequiv-1 */
      esl_stack_DiscardTopN(alts, x*3);         /* dig down to choice */
      esl_stack_IPop(alts, &d2);                /* pop it off, in rev order */
      esl_stack_IPop(alts, &d1);
      esl_stack_IPop(alts, &r);

      /* Now we know a best rule; figure out where we came from,
       * and push that info onto the <ns> stack.
       */
      if (verbose) {
        printf("-----------------------------------\n"); 
        printf("j=%d d=%d d1=%d d2=%d\n", j, d, d1, d2);
	printf("tracing %f\n", bestsc);
        printf("   rule(%d)\n", r+1);
      }

  
  }

  *ret_ct = ct;

  esl_stack_Destroy(ns);
  esl_stack_Destroy(alts);
  esl_randomness_Destroy(rng);
  return eslOK;
  
 ERROR:
  if (ns)   esl_stack_Destroy(ns);
  if (alts) esl_stack_Destroy(alts);
  if (rng)  esl_randomness_Destroy(rng);
  if (ct)   free(ct);
  return status;
}
}

GMX  
*GMX_Create(int L)
{
  GMX *gmx = NULL;
  int  pos = 1;
  int  j, d;
  int  status;

  gmx->L = L;
  
  ESL_ALLOC(gmx->dp,    sizeof(SCVAL *) * (L+1));
  ESL_ALLOC(gmx->dp[0], sizeof(SCVAL  ) * ((L+2)*(L+1))/2);

  for (j = 1; j <= gmx->L; j++)
    {
      gmx->dp[j] = gmx->dp[0] + pos;
      pos += j+1;
    }
  
  /* initialize */
  for (j = 0; j <= gmx->L; j++)
    for (d = 0; d <= j; d++)
      gmx->dp[j][d] = -eslINFINITY;
 
  return gmx;
 ERROR:
  return NULL;
}

void 
GMX_Destroy(GMX *gmx)
{
  if (gmx == NULL) return;
  
  if (gmx->dp) free(gmx->dp);
  free(gmx);
}

void
GMX_Dump(FILE *fp, GMX *gmx)
{
  int j, d;

  for (j = 0; j <= gmx->L; j++)
    {
      for (d = 0; d <= j; d++)
        fprintf(fp, "%f ", gmx->dp[j][d]);
      fputc('\n', fp);
    }
  fputc('\n', fp);
}


static int
dp_recursion(struct mutual_s *mi, GMX *cyk, int minloop, int j, int d, SCVAL *ret_sc, char *errbuf, int verbose)
{
  SCVAL sc = -eslINFINITY;
  int   d1;

  if (d == 0) { cyk->dp[j][d] = 0.; continue; }
  if (d == 1) { cyk->dp[j][d] = 0.; continue; }
  
  ESL_MAX(sc, cyk->dp[j-1][d-1]);
  for (d1 = minloop; d1 <= d; d1++)
    ESL_MAX(sc, cyk->dp[j-d+d1-1][d1-2] + cyk->dp[j][d-d1] + mi->COV->mx[j-d][j-d+d1]);

  *ret_sc = sc;
  return eslOK;
}
