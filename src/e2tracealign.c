/* Construction of multiple alignments from traces.
 * 
 * Contents:
 *   1. API for aligning sequence or MSA traces
 *   2. Internal functions used by the API
 *   3. Example
 *   4. Copyright and license.
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>


#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "e1_model.h"
#include "e2.h"
#include "e2_profilesq.h"
#include "e2tracealign.h"

static int e2_root(PSQ *sq, ESL_MSA *msa, char *errbuf, int verbose);
static int e2_descendents(int v, E2_TRACE **tr, PSQ **sq, int N, int dl, int dr, ESL_MSA *msa, E2_ALI e2ali, char *errbuf, int verbose);
static int insert(PSQ *psq, int sqidx, int i, int pos, int l, ESL_MSA *msa, int lc, int verbose);
static int add_residue(PSQ *psq, int i, int pos, char *asq, int lc, int verbose);
static int add_gap(int pos, char *asq, int lc, int verbose);

int
e2_TracealignSeqs(E2_TRACE *tr, PSQ *sql, PSQ *sqr, int a_idx, int dl_idx, int dr_idx, ESL_MSA *msa, char *errbuf, int verbose)
{
  char     *ancsq = NULL;
  PSQ      *sqrow;
  PSQ      *sqcol;
  int       pos = 0;		/* position in alignment 0..alen-1 */
  int       apos = 0;		/* position in already existing ancestral msa->aseq */
  int       idxrow, idxcol;
  int       z;
  int       firstmatch = -1;
  int       lastmatch = -1;
  int       unaligned = FALSE;
  int       status;

  sqrow = (tr->rowsq == e2P_SL)? (PSQ *)sql : (PSQ *)sqr;
  sqcol = (tr->rowsq == e2P_SL)? (PSQ *)sqr : (PSQ *)sql;

  idxrow = (tr->rowsq == e2P_SL)? dl_idx :  dr_idx;
  idxcol = (tr->rowsq == e2P_SL)? dr_idx :  dl_idx;
 
  /* make a copy of the ancestral aseq as is. Any gap that appear on it,
   * needs to be passed intact to the descendants; 
   * they are not part of this trace.
   */
  esl_strdup(msa->aseq[a_idx], -1, &ancsq);

  /* first pass to the trace to get the positions
   * before the possible final insert
   */

  for (z = 0; z < tr->N; z++) 
    {
      /* any gaps in the ancestral aseq comes from before, just propagate them, they
       * were not present in this trace 
       */
      if (ancsq[apos] == '-' || ancsq[apos] == '.') { 	
	pos ++; 
	apos ++; 
	z--;
	continue;
      }
      switch(tr->st[z]) {
      case e2T_BB:
	firstmatch = pos;
	break;
     case e2T_S:
       break;
     case e2T_T:
       break;
       /* fall through */
      case e2T_N1: if (z > 0 && tr->st[z-1] == tr->st[z]) pos ++; 
	break; 
      case e2T_J1: if (z > 0 && tr->st[z-1] == tr->st[z]) pos ++; 
	break; 
      case e2T_C1: if (z > 0 && tr->st[z-1] == tr->st[z]) pos ++; 
	break; 
      case e2T_N2: if (z > 0 && tr->st[z-1] == tr->st[z]) pos ++; 
	break; 
      case e2T_J2: if (z > 0 && tr->st[z-1] == tr->st[z]) pos ++; 
	break; 
      case e2T_C2: if (z > 0 && tr->st[z-1] == tr->st[z]) pos ++;
	break; 
      case e2T_XX:
	break;
      case e2T_IB: pos ++;
	firstmatch = pos;
	break;	
      case e2T_SS:
	lastmatch = pos;
	pos ++; apos ++;
	break;	
      case e2T_DS:
	lastmatch = pos;
	pos ++; apos ++;
	break;
      case e2T_IS: pos ++;
	break;	
      case e2T_SD:
	lastmatch = pos;
	pos ++; apos ++;
	break;	
      case e2T_DD:
	lastmatch = pos;
	pos ++;	apos ++;
	break;	
      case e2T_ID: pos ++;
	break;	
     case e2T_ii:  pos ++;	
	break;
      case e2T_BI: pos ++;	
	firstmatch = pos;
	break;
      case e2T_SI: pos ++;	
	break;	
      case e2T_DI: pos ++;	
	break;	
      case e2T_II: pos ++;	
	break;  	
      case e2T_EE: break;
      default: status = eslFAIL; goto ERROR;
      }
    }
  if (lastmatch == -1) unaligned = TRUE;

  pos  = 0;
  apos = 0;
  for (z = 0; z < tr->N; z++) 
    {
      /* any gaps in the ancestral aseq comes from before, just propagate them, they
       * were not present in this trace 
       */
      if (ancsq[apos] == '-') { 	
	add_gap(pos, msa->aseq[dl_idx], FALSE, verbose);
	add_gap(pos, msa->aseq[dr_idx], FALSE, verbose);
	pos ++; 
	apos ++; 
	z--;
	continue;
      }
      if (ancsq[apos] == '.') { 	
	add_gap(pos, msa->aseq[dl_idx], TRUE, verbose);
	add_gap(pos, msa->aseq[dr_idx], TRUE, verbose);
	pos ++; 
	apos ++; 
	z--;
	continue;
      }
 
      switch(tr->st[z]) {
      case e2T_S:
	break;
      case e2T_BB:
	break;
      case e2T_XX:
	break;	
      case e2T_N1:
	if (z > 0 && tr->st[z-1] == tr->st[z]) {
	  insert(sqrow, idxrow, tr->i[z], pos, 1, msa, TRUE, verbose);
	  pos++;
	}
	break;
	
      case e2T_J1:
	if (z > 0 && tr->st[z-1] == tr->st[z]) {
	  insert(sqrow, idxrow, tr->i[z], pos, 1, msa, TRUE, verbose);
	  pos ++;
	}
	break;
	
      case e2T_C1:
	if (z > 0 && tr->st[z-1] == tr->st[z]) {
	  insert(sqrow, idxrow, tr->i[z], pos, 1, msa, TRUE, verbose);
	  pos ++;
	}
	break;
	
      case e2T_N2:
	if (z > 0 && tr->st[z-1] == tr->st[z]) {
	  insert(sqcol, idxcol, tr->j[z], pos, 1, msa, TRUE, verbose);
	  pos ++;
	}
	break;
	
      case e2T_J2:
	if (z > 0 && tr->st[z-1] == tr->st[z]) {
	  insert(sqcol, idxcol, tr->j[z], pos, 1, msa, TRUE, verbose);
	  pos ++;
	}
	break;
	
      case e2T_C2:
	if (z > 0 && tr->st[z-1] == tr->st[z]) {
	  insert(sqcol, idxcol, tr->j[z], pos, 1, msa, TRUE, verbose);
	  pos ++;
	}
	break;
	
      case e2T_IB:
	insert(sqrow, idxrow, tr->i[z], pos, 1, msa, TRUE, verbose);
	pos ++;
	break;
	
      case e2T_SS:
	add_residue(sqrow, tr->i[z], pos, (tr->rowsq == e2P_SL)? msa->aseq[dl_idx] : msa->aseq[dr_idx], FALSE, verbose);
	add_residue(sqcol, tr->j[z], pos, (tr->rowsq == e2P_SL)? msa->aseq[dr_idx] : msa->aseq[dl_idx], FALSE, verbose);
	if (pos == firstmatch) unaligned = FALSE;
	if (pos == lastmatch)  unaligned = TRUE;
	pos ++;
	apos ++;
	break;
	
      case e2T_DS:
	add_gap    (                 pos, (tr->rowsq == e2P_SL)? msa->aseq[dl_idx] : msa->aseq[dr_idx], FALSE, verbose);
	add_residue(sqcol, tr->j[z], pos, (tr->rowsq == e2P_SL)? msa->aseq[dr_idx] : msa->aseq[dl_idx], FALSE, verbose);
	if (pos == firstmatch) unaligned = FALSE;
	if (pos == lastmatch)  unaligned = TRUE;
	pos ++;
	apos ++;
	break;
	
      case e2T_IS:
	insert(sqrow, idxrow, tr->i[z], pos, 1, msa, unaligned, verbose);
	pos ++;
	break;
	
      case e2T_SD:
	add_residue(sqrow, tr->i[z], pos, (tr->rowsq == e2P_SL)? msa->aseq[dl_idx] : msa->aseq[dr_idx], FALSE, verbose);
	add_gap    (                 pos, (tr->rowsq == e2P_SL)? msa->aseq[dr_idx] : msa->aseq[dl_idx], FALSE, verbose);
	if (pos == firstmatch) unaligned = FALSE;
	if (pos == lastmatch)  unaligned = TRUE;
	pos ++;
	apos ++;
	break;
	
      case e2T_DD:
	add_gap(pos, (tr->rowsq == e2P_SL)? msa->aseq[dl_idx] : msa->aseq[dr_idx], FALSE, verbose);
	add_gap(pos, (tr->rowsq == e2P_SL)? msa->aseq[dr_idx] : msa->aseq[dl_idx], FALSE, verbose);
	if (pos == firstmatch) unaligned = FALSE;
	if (pos == lastmatch)  unaligned = TRUE;
	pos ++;
	apos ++;
	break;
	
      case e2T_ID:
	insert(sqrow, idxrow, tr->i[z], pos, 1, msa, unaligned, verbose);
	pos ++;
	break;
	
     case e2T_ii:
	insert(sqrow, idxrow, tr->i[z], pos, 1, msa, unaligned, verbose);
	pos ++;	
	break;

      case e2T_BI:
	insert(sqcol, idxcol, tr->j[z], pos, 1, msa, unaligned, verbose);
	pos ++;	
	break;
      case e2T_SI:
	insert(sqcol, idxcol, tr->j[z], pos, 1, msa, unaligned, verbose);
	pos ++;	
	break;
	
      case e2T_DI:
	insert(sqcol, idxcol, tr->j[z], pos, 1, msa, unaligned, verbose);
	pos ++;	
	break;
	
      case e2T_II:
	insert(sqcol, idxcol, tr->j[z], pos, 1, msa, unaligned, verbose);
	pos ++;	
	break;
  	
      case e2T_EE:
	break;

      case e2T_T:
	break;

      default: status = eslFAIL; goto ERROR;
      }
    }
  
  if (strlen(msa->aseq[dl_idx]) != strlen(msa->aseq[dr_idx]) ||
      strlen(msa->aseq[dl_idx]) != strlen(msa->aseq[a_idx])    )
    {      
      printf("not an alignment for seqs %s (alen=%d) %s (alen=%d). ancestral  %s (alen=%d)\n", 
	     msa->sqname[dl_idx], (int)strlen(msa->aseq[dl_idx]), msa->sqname[dr_idx], (int)strlen(msa->aseq[dr_idx]), msa->sqname[a_idx], (int)strlen(msa->aseq[a_idx]));
      status = eslFAIL; goto ERROR;
    }

  free(ancsq);
  return eslOK;
   
 ERROR:
  if (ancsq) free(ancsq);
  return status;
}

int
e2_Tracealign(char *msaname, ESL_TREE *T, PSQ **sq, E2_TRACE **tr, ESL_MSA **ret_msa, E2_ALI e2ali, char *errbuf, int verbose)
{
  ESL_MSA        *msa = NULL;   /* alignment of leaf and node sequences */
  int             dl, dr; 
  int             v;            /* index for internal nodes */
  int             nnode;
  int             status;
  
  if (ret_msa == NULL) return eslOK;

  /* create the alignment */
  nnode = (T->N > 1)? T->N-1 : T->N;
  msa = esl_msa_Create(nnode+T->N, -1);
  msa->nseq = nnode+T->N;

  status = esl_strdup(msaname, -1, &(msa->name)); if (status != eslOK) goto ERROR;
  status = esl_strdup(msaname, -1, &(msa->acc));  if (status != eslOK) goto ERROR;

  e2_root(sq[0], msa, errbuf, verbose);

  for (v = 0; v < T->N-1; v++) {
    dl = T->left[v];
    dr = T->right[v];
    
    if (e2_descendents(v, tr, sq, T->N, dl, dr, msa, e2ali, errbuf, verbose) != eslOK)
	ESL_XFAIL(eslFAIL, errbuf, "%s\nfailed to grow from parent %d to daughthers %d %d", errbuf, v, dl, dr);
    
    if (verbose) {
      printf("\n%s\n%s\n", msa->sqname[v],  msa->aseq[v]);
      if (dl > 0) printf("%s\n%s\n", msa->sqname[dl],       msa->aseq[dl]);      
      else        printf("%s\n%s\n", msa->sqname[nnode-dl], msa->aseq[nnode-dl]);
      
      if (dr > 0) printf("%s\n%s\n", msa->sqname[dr],      msa->aseq[dr]);
      else        printf("%s\n%s\n", msa->sqname[nnode-dr],msa->aseq[nnode-dr]);
    }
  }
  
  *ret_msa = msa;
  return eslOK;

 ERROR:
  return status;
}


int
e2_root(PSQ *sq, ESL_MSA *msa, char *errbuf, int verbose)
{
  int status;

   ESL_ALLOC(msa->sqname[0], sizeof(char) * eslERRBUFSIZE);
  sprintf(msa->sqname[0], "v%d", 0);
  
  /* set the root length */
  msa->sqlen[0] = sq->n;
  msa->alen     = sq->n;
  
 /* generate the root sequence */
  psq_ConvertToAseq(sq, &(msa->aseq[0]));
 
  if (verbose) {
    printf("ancestral sequence:\n%s\n", msa->aseq[0]);
  }

  return eslOK;

 ERROR:
  return status;
}

int
e2_descendents(int v, E2_TRACE **tr, PSQ **sq, int N, int dl, int dr, ESL_MSA *msa, E2_ALI e2ali, char *errbuf, int verbose)
{
  int dl_idx, dr_idx;            /* indices */
  int nnodes = (N>1)? N-1 : N;
  int status;
  
  if (dl > 0) dl_idx = dl;
  else        dl_idx = nnodes - dl;
  if (dr > 0) dr_idx = dr;
  else        dr_idx = nnodes - dr;
  
  ESL_ALLOC(msa->sqname[dl_idx], sizeof(char) * eslERRBUFSIZE);
  ESL_ALLOC(msa->sqname[dr_idx], sizeof(char) * eslERRBUFSIZE);
  if (dl <= 0) esl_msa_SetSeqName(msa, dl_idx, sq[dl_idx]->name, strlen(sq[dl_idx]->name));
  else         sprintf(msa->sqname[dl_idx], "v%d", dl);
  if (dr <= 0) esl_msa_SetSeqName(msa, dr_idx, sq[dr_idx]->name, strlen(sq[dr_idx]->name));
  else         sprintf(msa->sqname[dr_idx], "v%d", dr);

  ESL_ALLOC(msa->aseq[dl_idx], sizeof(char) * (msa->alen+1));
  ESL_ALLOC(msa->aseq[dr_idx], sizeof(char) * (msa->alen+1));
  msa->aseq[dl_idx][msa->alen] = '\0';
  msa->aseq[dr_idx][msa->alen] = '\0';

  if (e2ali == E2F) {
    psq_ConvertToAseq(sq[dl_idx], &(msa->aseq[dl_idx]));
    psq_ConvertToAseq(sq[dr_idx], &(msa->aseq[dr_idx]));
  }
  else if (e2ali == E2FHMMER) {
    if ((status = e2_TracealignSeqs(tr[v], sq[dl_idx], sq[dr_idx], v, dl_idx, dr_idx, msa, errbuf, verbose)) != eslOK) goto ERROR;
  }
  else {
    if ((status = e2_TracealignSeqs(tr[v], sq[dl_idx], sq[dr_idx], v, dl_idx, dr_idx, msa, errbuf, verbose)) != eslOK) goto ERROR;
  }
  msa->sqlen[dl_idx] = sq[dl_idx]->n;
  msa->sqlen[dr_idx] = sq[dr_idx]->n;

  return eslOK;

 ERROR:
  return status;
}


/* Extend the alignment by l columns starting at 'pos'. 
 */
static  int 
insert(PSQ *psq, int sqidx, int i, int pos, int l, ESL_MSA *msa, int lc, int verbose)
{
  char *new = NULL;
  int   newalen;
  int   x;
  int   n;
  int   status;

  newalen = msa->alen + l;

  for (x = 0; x < msa->nseq; x ++) {
    if (msa->aseq[x] == NULL) continue;
    
    ESL_ALLOC(new, sizeof(char) * (newalen+1));
    new[newalen] = '\0';
    
    /* copy residues before 'pos' */
    for (n = 0; n < pos; n ++) new[n] = msa->aseq[x][n];
    /* move over residues past 'pos' */
    for (n = pos; n < msa->alen; n ++)
      new[n+l] = msa->aseq[x][n];
     
    /* the insertion */
    for (n = 0; n < l; n ++) {
      if (x == sqidx) { add_residue(psq, i+n, pos+n, new, lc, verbose); }
      else            { add_gap(              pos+n, new, lc, verbose); }
    }
    free(msa->aseq[x]);
    msa->aseq[x] = new;
  }
  
  msa->alen = newalen; 
  return eslOK;
  
 ERROR:
  return status;
}

static int
add_residue(PSQ *psq, int i, int pos, char *asq, int lc, int verbose) 
{
  float *p = NULL;
  int    x;
  int    status;
  
  ESL_ALLOC(p, sizeof(float) * (psq->abc->K+1));

  psq_ProfProbs(i, psq, p);
  x = esl_vec_FArgMax(p, psq->abc->K);
  asq[pos] = (lc)? tolower(psq->abc->sym[x]) : psq->abc->sym[x];
  free(p);

  return eslOK;

 ERROR:
  if (p) free(p);
  return status;
}

static int
add_gap(int pos, char *asq, int lc, int verbose) 
{
  asq[pos] = (lc)? '.' : '-';
  return eslOK;
}

