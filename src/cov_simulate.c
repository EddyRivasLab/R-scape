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

static int cov_evolve_root_ungapped(ESL_RANDOMNESS *r, ESL_TREE *T, E1_RATE *e1rate, struct ribomatrix_s *ribusum,
				    int *ct, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose);
static int cov_evolve_ascendant_to_descendant_ungapped(ESL_RANDOMNESS *r, int *ret_idx, int *nidx, int p, int d, double time, ESL_TREE *T, E1_RATE *e1rate,
						       struct ribomatrix_s *ribosum, int *ct, ESL_MSA *msa, double tol, char *errbuf, int verbose);
static int cov_emit_ungapped(ESL_RANDOMNESS *r, E1_MODEL *e1model, ESL_DMATRIX *riboP, double *riboM, int aidx, int didx, int *ct, ESL_MSA *msa, 
			     char *errbuf, int verbose);
static int cov_substitute(ESL_RANDOMNESS *r, int pos, int aidx, int didx, E1_MODEL *e1model, ESL_DMATRIX *riboP, int *ct, ESL_MSA *msa, int verbose);
static int cov_addres(ESL_RANDOMNESS *r, int idx, int pos, double *P, ESL_MSA *msa, int verbose);
static int cov_addpair(ESL_RANDOMNESS *r, int idx, int pos, int pair, double *riboP, int *ct, ESL_MSA *msa, int verbose);
static int asq_rlen(ESL_DSQ *dsq, int alen, int K);

int 
cov_GenerateAlignment(ESL_RANDOMNESS *r, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
		      ESL_MSA **ret_msafull, int noss, double tol, char *errbuf, int verbose)
{
  int status;

  status = cov_GenerateAlignmentUngapped(r, T, root, e1rate, ribosum, ret_msafull, noss, tol, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  status = cov_GenerateAlignmentAddGaps(r, T, root, e1rate, ribosum, *ret_msafull, noss, tol, errbuf, verbose);
  if (status != eslOK) goto ERROR;
		      
  return eslOK;

 ERROR:
  return status;
}

int 
cov_GenerateAlignmentUngapped(ESL_RANDOMNESS *r, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
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
  msa = esl_msa_CreateDigital(root->abc, nnode+T->N, 0);

  esl_sprintf(&name, "syntethic_L%d_N%d", L, T->N);
  status = esl_strdup(name, -1, &(msa->name)); if (status != eslOK) goto ERROR;
  status = esl_strdup(name, -1, &(msa->acc));  if (status != eslOK) goto ERROR;

  /* we have the root sequence and perhaps a secondary structure */
  ESL_ALLOC(ct, sizeof(int)*(L+1));
  esl_vec_ISet(ct, L+1, 0);
  if (noss || root->ss_cons == NULL) esl_vec_ISet(ct, L+1, 0);
  else                               esl_wuss2ct(root->ss_cons, L, ct);

  /* evolve root */
  if (cov_evolve_root_ungapped(r, T, e1rate, ribosum, ct, msa, &idx, nidx, tol, errbuf, verbose)!= eslOK) { status = eslFAIL; goto ERROR; }

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


static int 
cov_evolve_root_ungapped(ESL_RANDOMNESS *r, ESL_TREE *T, E1_RATE *e1rate, struct ribomatrix_s *ribosum,
			 int *ct, ESL_MSA *msa, int *ret_idx, int *nidx, double tol, char *errbuf, int verbose)
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
    
    if (cov_evolve_ascendant_to_descendant_ungapped(r, ret_idx, nidx, v, dl, ld, T, e1rate, ribosum, ct, msa, tol, errbuf, verbose) != eslOK)
      ESL_XFAIL(eslFAIL, errbuf, "%s\nfailed to evolve from parent %d to daughther %d after time %f", errbuf, v, dl, ld);
    if (cov_evolve_ascendant_to_descendant_ungapped(r, ret_idx, nidx, v, dr, rd, T, e1rate, ribosum, ct, msa, tol, errbuf, verbose) != eslOK)
      ESL_XFAIL(eslFAIL, errbuf, "%s\nfailed to evolve from parent %d to daughther %d after time %f", errbuf, v, dr, rd);

     
    if (verbose) {
      printf("%s\n%s\n", msa->sqname[nidx[v]],  msa->ax[nidx[v]]);
      if (dl > 0) printf("%s\n%s\n", msa->sqname[nidx[dl]],  msa->ax[nidx[dl]]);      
      else        printf("%s\n%s\n", msa->sqname[T->N-1-dl], msa->ax[T->N-1-dl]);

      if (dr > 0) printf("%s\n%s\n", msa->sqname[nidx[dr]],  msa->ax[nidx[dr]]);
      else        printf("%s\n%s\n", msa->sqname[T->N-1-dr], msa->ax[T->N-1-dr]);
     }
  }
  
  return eslOK;

 ERROR:
  return status;
}

static int
cov_evolve_ascendant_to_descendant_ungapped(ESL_RANDOMNESS *r, int *ret_idx, int *nidx, int p, int d, double time, ESL_TREE *T, E1_RATE *e1rate, 
					    struct ribomatrix_s *ribosum, int *ct, ESL_MSA *msa, double tol, char *errbuf, int verbose)
{
  E1_MODEL    *e1model = NULL;
  ESL_DMATRIX *riboP = NULL;
  ESL_DMATRIX *riboQ = ribosum->bprsQ;
  double      *riboM = ribosum->bprsM;
  int          idx;
  int          d_idx;
  int          status;
  
  //printf("\nEMIT node %d --> %d %f\n", p, d, time);
  /* evolve the e1rate and ribosums to time */
  e1model = e1_model_Create(e1rate, time, NULL, (float *)ribosum->bg, e2_GLOBAL, asq_rlen(msa->ax[nidx[p]], msa->alen, msa->abc->K), msa->abc, tol, errbuf, verbose); 
  if (e1model == NULL) ESL_XFAIL(eslFAIL, errbuf, "failed to evolve e1model to time %f", time);

  /* evolve the riboQ(aa'|bb') rate to time */
  riboP = ratematrix_ConditionalsFromRate(time, riboQ, tol, errbuf, verbose);
  if (riboP == NULL) ESL_XFAIL(eslFAIL, errbuf, "failed to evolve ribosum to time %f", time);

  /* generate the descendant sequence */
  idx = *ret_idx;
  if (d > 0) {
    d_idx   = idx++;
    nidx[d] = d_idx;
  }
  else {
    d_idx = T->N - 1 - d;
  }

  ESL_ALLOC(msa->sqname[d_idx], sizeof(char) * eslERRBUFSIZE);
  ESL_ALLOC(msa->ax[d_idx],     sizeof(int) * msa->alen);

  if (d <= 0) sprintf(msa->sqname[d_idx], "T%d", -d);
  else        sprintf(msa->sqname[d_idx], "v%d",  d);

  if ((status = cov_emit_ungapped(r, e1model, riboP, riboM, nidx[p], d_idx, ct, msa, errbuf, verbose)) != eslOK) goto ERROR;
  msa->sqlen[d_idx] = strlen(msa->aseq[d_idx]);

  /* check */
  if (msa->alen != msa->sqlen[d_idx]) {
    printf("seq at node %d len is %d but alen %d\n", d, (int)msa->sqlen[d_idx], (int)msa->alen);
    ESL_XFAIL(eslFAIL, errbuf, "bad msa at node %d", d);
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
  int       i = -1;		/* position in the ancestral sequence 0, L-1 */
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

      /* bump i, pos if needed */
      if (st == e1T_S) { i ++; pos ++; }
      
      /* a transit to alen is a transit to the E state */
      if (i == L) { st = e1T_E; pos = 0; }

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
  int apos  = msa->ax[aidx][pos];
  int apair = ct[apos];
  int K   = e1model->abc->K;

  if (apos > K || apair > K) return eslFAIL;

  cov_addres (r, didx, apos,        e1model->sub->mx[apos],           msa, verbose);
  cov_addpair(r, didx, apos, apair, riboP->mx[IDX(apos,apair,K)], ct, msa, verbose);

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


static int
asq_rlen(ESL_DSQ *dsq, int alen, int K)
{
  int rlen = 0;
  int i;

  for (i = 1; i <= alen; i ++) {
    if (dsq[i] < K) rlen ++;
  }

  return rlen;
}
