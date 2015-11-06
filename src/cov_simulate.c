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
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "e2.h"
#include "e1_model.h"
#include "e1_rate.h"
#include "cov_simulate.h"
#include "ribosum_matrix.h"

static int  cov_add_root(ESL_MSA *root, ESL_MSA *msa, int *ret_idx, int *nidx, int noss, char *errbuf, int verbose);
static int  cov_evolve_root_ungapped_tree(ESL_RANDOMNESS *r, ESL_TREE *T, E1_RATE *e1rate, struct ribomatrix_s *ribusum,
					  int *ct, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose);
static int  cov_evolve_root_ungapped_star(ESL_RANDOMNESS *r, int N, double abl, E1_RATE *e1rate, struct ribomatrix_s *ribusum,
					  int *ct, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose);
static int  cov_evolve_ascendant_to_descendant_ungapped(ESL_RANDOMNESS *r, int *ret_idx, int *nidx, int p, int d, double time, int inodes, E1_RATE *e1rate,
							struct ribomatrix_s *ribosum, int *ct, ESL_MSA *msa, double tol, char *errbuf, int verbose);
static int  cov_emit_ungapped(ESL_RANDOMNESS *r, E1_MODEL *e1model, ESL_DMATRIX *riboP, double *riboM, int aidx, int didx, int *ct, ESL_MSA *msa, 
			      char *errbuf, int verbose);
static int  cov_substitute(ESL_RANDOMNESS *r, int pos, int aidx, int didx, E1_MODEL *e1model, ESL_DMATRIX *riboP, int *ct, ESL_MSA *msa, int verbose);
static int  cov_addres(ESL_RANDOMNESS *r, int idx, int pos, double *P, ESL_MSA *msa, int verbose);
static int  cov_addpair(ESL_RANDOMNESS *r, int idx, int pos, int pair, double *riboP, int *ct, ESL_MSA *msa, int verbose);
static void ax_dump(ESL_DSQ *dsq);

int 
cov_GenerateAlignment(ESL_RANDOMNESS *r, int N, double abl, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
		      ESL_MSA **ret_msafull, int noss, double tol, char *errbuf, int verbose)
{
  int status;

  status = cov_GenerateAlignmentUngapped(r, N, abl, T, root, e1rate, ribosum, ret_msafull, noss, tol, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  status = cov_GenerateAlignmentAddGaps(r, T, root, e1rate, ribosum, *ret_msafull, noss, tol, errbuf, verbose);
  if (status != eslOK) goto ERROR;
		      
  return eslOK;

 ERROR:
  return status;
}

int 
cov_GenerateAlignmentUngapped(ESL_RANDOMNESS *r, int N, double abl, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
			      ESL_MSA **ret_msafull, int noss, double tol, char *errbuf, int verbose)
{
  ESL_MSA        *msa = NULL;     /* alignment of leaf and node sequences */
  char           *name = NULL;
  int            *nidx = NULL;    /* node index array */
  int            *ct = NULL;
  int             L = root->alen;
  int             idx = 0;        /* node index */
  int             inode;          // number of internal nodes
  int             tnode;          // number of total nodes
  int             n;
  int             status;

  if (T && N != T->N) {  status = eslFAIL; goto ERROR; }

  /* indexing of internal nodes */
  inode = (T)? T->N - 1 : 1;
  ESL_ALLOC(nidx, sizeof(int) * inode);
  for(n = 0; n < inode-1; n ++)
    nidx[n] = -1;

  /* create the alignment */
  tnode = inode + N;
  msa = esl_msa_CreateDigital(root->abc, tnode, 0);

  esl_sprintf(&name, "%s_syntethic_L%d_N%d", root->name, L, N);
  status = esl_strdup(name, -1, &(msa->name)); if (status != eslOK) goto ERROR;
  status = esl_strdup(name, -1, &(msa->acc));  if (status != eslOK) goto ERROR;
  
  /* The ct vector with the secondary structure */
  ESL_ALLOC(ct, sizeof(int)*(L+1));
  esl_vec_ISet(ct, L+1, 0);
  if (noss || root->ss_cons == NULL) esl_vec_ISet(ct, L+1, 0);
  else                               esl_wuss2ct(root->ss_cons, L, ct);

  /* we have the root and a secondary structure, add both to the growing alignment */
  if (cov_add_root(root, msa, &idx, nidx, noss, errbuf, verbose)!= eslOK) { status = eslFAIL; goto ERROR; }

 /* evolve root */
  if (T) status = cov_evolve_root_ungapped_tree(r, T,      e1rate, ribosum, ct, msa, &idx, nidx, tol, errbuf, verbose);
  else   status = cov_evolve_root_ungapped_star(r, N, abl, e1rate, ribosum, ct, msa, &idx, nidx, tol, errbuf, verbose);
  
    if (status != eslOK) { status = eslFAIL; goto ERROR; }

  *ret_msafull = msa;

  free(ct);
  return eslOK;

 ERROR:
  if (ct) free(ct);
  return status;
}

int 
cov_GenerateAlignmentAddGaps(ESL_RANDOMNESS *r, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
			     ESL_MSA *msafull, int noss, double tol, char *errbuf, int verbose)
{

  return eslOK;
}

/* ---- static functions ---- */

static int 
cov_add_root(ESL_MSA *root, ESL_MSA *msa, int *ret_idx, int *nidx, int noss, char *errbuf, int verbose) 
{
  int L = root->alen;
  int i;
  int idx;
  int status;

  idx     = *ret_idx;
  nidx[0] = idx++;

  ESL_ALLOC(msa->sqname[nidx[0]], sizeof(char) * eslERRBUFSIZE);
  sprintf(msa->sqname[nidx[0]], "v%d", 0);
 
  /* set the root length */
  msa->sqlen[nidx[0]] = msa->alen = L;

  /* copy the root sequence to msa */
  ESL_ALLOC(msa->ax[nidx[0]], sizeof(ESL_DSQ) * (L+2));
  for (i = 0; i <= L; i ++)
    msa->ax[nidx[0]][i] = root->ax[0][i];
  msa->ax[nidx[0]][L+1] = eslDSQ_SENTINEL;
  
  /* copy the secondary structure */
  if (!noss) esl_strdup(root->ss_cons, -1, &(msa->ss_cons));

  if (verbose) {
    for (i = 1; i <= L; i ++)
      printf("%d", msa->ax[nidx[0]][i]);
    printf("\n");
    printf("%s\n", msa->ss_cons);
  }

  *ret_idx = idx;
  return eslOK;

 ERROR:
  return status;
}


static int 
cov_evolve_root_ungapped_star(ESL_RANDOMNESS *r, int N, double abl, E1_RATE *e1rate, struct ribomatrix_s *ribosum,
			      int *ct, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose)
{
  int inodes = 1;  // only the root node
  int v = 0;       // only the root node
  int i;           // sequence index
  int status;

  for (i = 0; i < N; i ++) {
    if (cov_evolve_ascendant_to_descendant_ungapped(r, ret_idx, nidx, v, -i, abl, inodes, e1rate, ribosum, ct, msa, tol, errbuf, verbose) != eslOK)
      ESL_XFAIL(eslFAIL, errbuf, "%s\nfailed to evolve from parent %d to daughther %d after time %f", errbuf, v, i, abl);

    if (verbose) {
      printf("\n%s | len %" PRId64 " \n", msa->sqname[nidx[v]], esl_abc_dsqlen(msa->ax[nidx[v]])); ax_dump(msa->ax[nidx[v]]); 
      printf("%s[%d] | len %" PRId64 " \n", msa->sqname[i], i, esl_abc_dsqlen(msa->ax[i])); ax_dump(msa->ax[i]); 
    }

  }
  return eslOK;

 ERROR:
  return status;
}

static int 
cov_evolve_root_ungapped_tree(ESL_RANDOMNESS *r, ESL_TREE *T, E1_RATE *e1rate, struct ribomatrix_s *ribosum,
			      int *ct, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose)
{
  double ld, rd;
  int    v;       /* index for internal nodes */
  int    dl, dr;
  int    inodes = (T->N>1)? T->N-1:1;
  int    status;
  
  for (v = 0; v < T->N-1; v++) {
    dl = T->left[v];
    dr = T->right[v];
    
    ld = T->ld[v];
    rd = T->rd[v];
    
    if (cov_evolve_ascendant_to_descendant_ungapped(r, ret_idx, nidx, v, dl, ld, inodes, e1rate, ribosum, ct, msa, tol, errbuf, verbose) != eslOK)
      ESL_XFAIL(eslFAIL, errbuf, "%s\nfailed to evolve from parent %d to daughther %d after time %f", errbuf, v, dl, ld);
    if (cov_evolve_ascendant_to_descendant_ungapped(r, ret_idx, nidx, v, dr, rd, inodes, e1rate, ribosum, ct, msa, tol, errbuf, verbose) != eslOK)
      ESL_XFAIL(eslFAIL, errbuf, "%s\nfailed to evolve from parent %d to daughther %d after time %f", errbuf, v, dr, rd);

    if (verbose) {
      printf("\n%s | len %" PRId64 " \n", msa->sqname[nidx[v]], esl_abc_dsqlen(msa->ax[nidx[v]])); ax_dump(msa->ax[nidx[v]]); 
      if (dl > 0) { printf("%s[%d] | len %" PRId64 "\n", msa->sqname[nidx[dl]], nidx[dl],  esl_abc_dsqlen(msa->ax[nidx[dl]]));  ax_dump(msa->ax[nidx[dl]]);  }  
      else        { printf("%s[%d] | len %" PRId64 "\n", msa->sqname[T->N-1-dl], T->N-1-dl, esl_abc_dsqlen(msa->ax[T->N-1-dl])); ax_dump(msa->ax[T->N-1-dl]); } 
      
      if (dr > 0) { printf("%s[%d] | len %" PRId64 "\n", msa->sqname[nidx[dr]],  nidx[dr], esl_abc_dsqlen(msa->ax[nidx[dr]]));  ax_dump(msa->ax[nidx[dr]]);  }
      else        { printf("%s[%d] | len %" PRId64 "\n", msa->sqname[T->N-1-dr], T->N-1-dr, esl_abc_dsqlen(msa->ax[T->N-1-dr])); ax_dump(msa->ax[T->N-1-dr]); }
     }
  }
  
  return eslOK;

 ERROR:
  return status;
}

static int
cov_evolve_ascendant_to_descendant_ungapped(ESL_RANDOMNESS *r, int *ret_idx, int *nidx, int p, int d, double time, int inodes, E1_RATE *e1rate, 
					    struct ribomatrix_s *ribosum, int *ct, ESL_MSA *msa, double tol, char *errbuf, int verbose)
{
  E1_MODEL    *e1model = NULL;
  ESL_DMATRIX *riboP = NULL;
  ESL_DMATRIX *riboQ = ribosum->bprsQ;
  double      *riboM = ribosum->bprsM;
  int          idx;
  int          didx;
  int          status;

  printf("\nEMIT node %d --> %d %f\n", p, d, time);
  /* evolve the e1rate and ribosums to time */
  e1model = e1_model_Create(e1rate, time, NULL, NULL, e2_GLOBAL, esl_abc_dsqrlen(msa->abc, msa->ax[nidx[p]]), msa->abc, tol, errbuf, verbose); 
  if (e1model == NULL) ESL_XFAIL(eslFAIL, errbuf, "failed to evolve e1model to time %f", time);
  if (verbose) esl_dmatrix_Dump(stdout, e1model->sub, NULL, NULL);
  
  /* evolve the riboQ(aa'|bb') rate to time */
  riboP = ratematrix_ConditionalsFromRate(time, riboQ, tol, errbuf, verbose);
  if (riboP == NULL) ESL_XFAIL(eslFAIL, errbuf, "failed to evolve ribosum to time %f", time);
  if (verbose) esl_dmatrix_Dump(stdout, riboP, NULL, NULL);

  /* generate the descendant sequence */
  idx = *ret_idx;
  if (d > 0) {
    didx    = idx++;
    nidx[d] = didx;
  }
  else {
    didx = inodes - d;
  }

  ESL_ALLOC(msa->sqname[didx], sizeof(char)    * eslERRBUFSIZE);
  ESL_ALLOC(msa->ax[didx],     sizeof(ESL_DSQ) * (msa->alen+2));

  if (d <= 0) sprintf(msa->sqname[didx], "T%d", -d+1);
  else        sprintf(msa->sqname[didx], "v%d",  d);

  if ((status = cov_emit_ungapped(r, e1model, riboP, riboM, nidx[p], didx, ct, msa, errbuf, verbose)) != eslOK) goto ERROR;
  msa->sqlen[didx] = esl_abc_dsqrlen(msa->abc, msa->ax[didx]);

  /* check */
  if (msa->alen != msa->sqlen[didx]) {
    printf("seq at node %d len is %d but alen %d\n", d, (int)msa->sqlen[didx], (int)msa->alen);
    ESL_XFAIL(eslFAIL, errbuf, "bad msa at node %d", d);
  }

  if (verbose) {
    printf("ancestral[%d] %s  | len = %" PRId64 " \n", nidx[p], msa->sqname[nidx[p]], esl_abc_dsqlen(msa->ax[nidx[p]]));
    ax_dump(msa->ax[nidx[p]]);
    printf("descendant[%d] %s | len = %" PRId64 " \n", didx, msa->sqname[didx], esl_abc_dsqlen(msa->ax[didx]));
    ax_dump(msa->ax[didx]);
  }
  
  *ret_idx = idx;

  e1_model_Destroy(e1model);
  esl_dmatrix_Destroy(riboP);
  return eslOK;
  
 ERROR:
  if (e1model) e1_model_Destroy(e1model);
  if (riboP)   esl_dmatrix_Destroy(riboP);
  return status;
}

int
cov_emit_ungapped(ESL_RANDOMNESS *r, E1_MODEL *e1model, ESL_DMATRIX *riboP, double *riboM, int aidx, int didx, int *ct, ESL_MSA *msa, char *errbuf, int verbose)
{
  int       L = msa->alen;      /* length of ancestral sq, it includes gaps */
  int       pos = 0;		/* position in alignment 1..alen */
  int       st = e1T_B;  	/* state type */
  int       status;

  /* Renormalize transitions so that T(X->E) = 0 */
  e1_model_RenormNoIndels(e1model);
    
  while (st != e1T_E)
    {
      /* Sample next state type, given current state type */
      switch (st) {
      case e1T_B:
	st = e1T_S; break;
	break;
      case e1T_S:
	st = e1T_S; break;
	break;
      case e1T_D:
	ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	break;
      case e1T_I:
	ESL_XEXCEPTION(eslEINCONCEIVABLE, "impossible.");  	    
	break;
      default: ESL_XEXCEPTION(eslECORRUPT, "impossible state reached during emission");
      }

      /* bump pos if needed */
      if (st == e1T_S) pos ++;
      
      /* a transit to alen+1 is a transit to the E state */
      if (pos == L+1) {
	st = e1T_E;
	msa->ax[didx][pos] = eslDSQ_SENTINEL;
      }

      if (st == e1T_S)  { /* a substitution, sample a residue */
	status = cov_substitute(r, pos, aidx, didx, e1model, riboP, ct, msa, verbose);  
	if (status != eslOK) goto ERROR;
      }
   }

  return eslOK;
   
 ERROR:
   return status;
}


static int
cov_substitute(ESL_RANDOMNESS *r, int pos, int aidx, int didx, E1_MODEL *e1model, ESL_DMATRIX *riboP, int *ct, ESL_MSA *msa, int verbose)
{
  int apos = msa->ax[aidx][pos];
  int apair;
  int pair = ct[pos];
  int K    = e1model->abc->K;
  
  if (apos > K) return eslFAIL;
  
  if (pair == 0) cov_addres (r, didx, pos, e1model->sub->mx[apos], msa, verbose);
  else {
    apair = msa->ax[aidx][pair]; if (apair > K) return eslFAIL;
    cov_addpair(r, didx, pos, pair, riboP->mx[IDX(apos,apair,K)], ct, msa, verbose);
  }
  
  return eslOK;
  
}

static int
cov_addres(ESL_RANDOMNESS *r, int idx, int pos, double *P, ESL_MSA *msa, int verbose)
{
  double pdf = 0.0;
  double x = esl_random(r);
  int    dim = msa->abc->K;
  int    i;

  for (i = 0; i < dim; i++) {
    pdf += P[i];
    if (pdf > x) break;
  }
  if (i == dim) i = dim-1;
  
  msa->ax[idx][pos] = i;

  return eslOK;    
}


static int
cov_addpair(ESL_RANDOMNESS *r, int idx, int pos, int pair, double *riboP, int *ct, ESL_MSA *msa, int verbose)
{
  double pdf = 0.0;
  double x = esl_random(r); 
  int    dim = msa->abc->K;
  int    dim2 = dim*dim;
  int    i;
  
  if (pair == 0)   return eslOK;  // unpaired base
  if (pair <  pos) return eslOK;  // a 3' basepaired position already emitted with ribosum
  
  for (i = 0; i < dim2; i++) {
    pdf += riboP[i];
    if (pdf > x) break;
  }
  if (i == dim2) i = dim2-1;
  
  msa->ax[idx][pos]  = i/dim;
  msa->ax[idx][pair] = i%dim;

  return eslOK;    
}


static void
ax_dump(ESL_DSQ *dsq)
{
  int i;
  for (i = 1; i <= esl_abc_dsqlen(dsq); i ++) 
    printf("%d", dsq[i]);
  printf("\n");
}
