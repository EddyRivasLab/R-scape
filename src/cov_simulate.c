/* cov_simulate - funtions to generate synthetic msa with covariation
 * 
 * Contents:
 *
 * ER, Mon Nov  2 13:24:54 EST 2015 [NorthWest] 
 */

#include "p7_config.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_stats.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "e2.h"
#include "cov_simulate.h"
#include "ribosum_matrix.h"

static int cov_evolve_root(ESL_RANDOMNESS *r, ESL_TREE *T, E1_RATE *e1rate, struct ribomatrix_s *ribosum,
			   ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose);

int 
cov_GenerateAlignment(ESL_RANDOMNESS *r, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
		      ESL_MSA **ret_msafull, int noss, double tol, char *errbuf, int verbose)
{
  ESL_MSA        *msa = NULL;     /* alignment of leaf and node sequences */
  char           *name = NULL;
  int            *nidx = NULL;    /* node index array */
  int            *ct = NULL;
  int             L = root->alen;
  int             idx = 0;        /* node index */
  int             nnode;
  int             n;
  int             status;

  /* indexing of internal nodes */
  ESL_ALLOC(nidx, sizeof(int) * T->N-1);
  for(n = 0; n < T->N-1; n ++)
    nidx[n] = -1;

  /* create the alignment */
  nnode = (T->N > 1)? T->N-1 : T->N;
  msa = esl_msa_Create(nnode+T->N, 0);

  esl_sprintf(&name, "syntethic_L%d_N%d", L, T->N);
  status = esl_strdup(name, -1, &(msa->name)); if (status != eslOK) goto ERROR;
  status = esl_strdup(name, -1, &(msa->acc));  if (status != eslOK) goto ERROR;

  /* we have the root sequence and perhaps a secondary structure */
  ESL_ALLOC(ct, sizeof(int)*(L+1));
  esl_vec_ISet(ct, L+1, 0);
  if (noss || root->ss_cons == NULL) esl_vec_ISet(ct, L+1, 0);
  else                               esl_wuss2ct(root->ss_cons, L, ct);

  /* evolve root */
  if (cov_evolve_root(r, T, e1rate, ribosum, msa, &idx, nidx, tol, errbuf, verbose)!= eslOK) { status = eslFAIL; goto ERROR; }

  *ret_msafull = msa;

  free(ct);
  return eslOK;

 ERROR:
  if (ct) free(ct);
  return status;
}


static int 
cov_evolve_root(ESL_RANDOMNESS *r, ESL_TREE *T, E1_RATE *e1rate, struct ribomatrix_s *ribosum,
		ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose)
{
  double ld, rd;
  int    v;       /* index for internal nodes */
  int    dl, dr;
  int    status;
  
  for (v = 0; v < T->N-1; v++) {
    dl = T->left[v];
    dr = T->right[v];
    
    ld = T->ld[v];
    rd = T->rd[v];
    
     
    if (verbose) {
      printf("%s\n%s\n", msa->sqname[nidx[v]],  msa->aseq[nidx[v]]);
      if (dl > 0) printf("%s\n%s\n", msa->sqname[nidx[dl]],  msa->aseq[nidx[dl]]);      
      else        printf("%s\n%s\n", msa->sqname[T->N-1-dl], msa->aseq[T->N-1-dl]);

      if (dr > 0) printf("%s\n%s\n", msa->sqname[nidx[dr]],  msa->aseq[nidx[dr]]);
      else        printf("%s\n%s\n", msa->sqname[T->N-1-dr], msa->aseq[T->N-1-dr]);
     }
  }
  
  return eslOK;

 ERROR:
  return status;
}
