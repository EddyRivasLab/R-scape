/*  msamanip
 *
 * ER, Fri Mar 28 12:42:06 EDT 2014 [Janelia] 
 * SVN $Id:$
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
#include "esl_distance.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"


#include "e2_tree.h"
#include "fastsp.h"
#include "msamanip.h"


static int calculate_Cstats(ESL_MSA *msa, int *ret_maxilen, int *ret_totilen, int *ret_totinum, double *ret_avginum, double *ret_stdinum, double *ret_avgilen, double *ret_stdilen, double *ret_avgsqlen, double *ret_stdsqlen, int *ret_anclen);
static int calculate_Xstats(ESL_MSA *msa, int *ret_maxilen, int *ret_totilen, int *ret_totinum, double *ret_avginum, double *ret_stdinum, double *ret_avgilen, double *ret_stdilen, double *ret_avgsqlen, double *ret_stdsqlen, int *ret_anclen);
static int reorder_msa(ESL_MSA *msa, int *order, char *errbuf);


int 
msamanip_CalculateCT( ESL_MSA *msa, int **ret_ct, int *ret_nbpairs, char *errbuf)
{
  int *ct = NULL;
  int  nbpairs = 0;
  int  i, j;
  int  status;
  
  ESL_ALLOC(ct, sizeof(int) * (msa->alen+1));
  if (msa->ss_cons) esl_wuss2ct(msa->ss_cons, msa->alen, ct);
  else ESL_XFAIL(eslFAIL, errbuf, "no ss for msa");
  
  for (i = 0; i < msa->alen-1; i ++)
    for (j = i+1; j < msa->alen; j ++)
    	if (ct[i+1] == j+1) nbpairs ++;

  *ret_ct      = ct;  
  *ret_nbpairs = nbpairs;
  return eslOK;

 ERROR:
  if (ct) free(ct);
  return status;
}

int
msamanip_NonHomologous(ESL_ALPHABET *abc, ESL_MSA *msar, ESL_MSA *msae, int *ret_nhr, int *ret_nhe, int *ret_hr, int *ret_he, int *ret_hre, char *errbuf)
{
  ESL_SQ *rsq  = NULL;
  ESL_SQ *esq  = NULL;
  int     nsq = msar->nseq;
  int     nhr = 0;      // number of nonhomologous  positions in reference
  int     nhe = 0;      // number of nonhomologous  positions in inferred
  int     hr  = 0;      // number of homologous     positions in reference
  int     he  = 0;      // number of homologous     positions in inferred
  int     hre = 0;      // number of homologous     positions in reference and inferred
  int     i;
  int     x;
  int     status;
  
  if (msar->aseq == NULL) return eslOK;
  
  if (msae == NULL)
    {
      for (i = 0; i < nsq; i++) {
	esl_sq_FetchFromMSA(msar, i, &rsq);
	for (x = 0; x < rsq->n; x++) {	  
	  if (islower(rsq->seq[x])) nhr ++;
	  if (isupper(rsq->seq[x])) hr ++;
	}
	esl_sq_Destroy(rsq); rsq = NULL;
      }
      *ret_nhr = nhr;
      *ret_nhe = 0;
      *ret_hr  = hr;
      *ret_he  = 0;
      *ret_hre = 0;
      return eslOK;
    }

  for (i = 0; i < nsq; i++) {
    esl_sq_FetchFromMSA(msar, i, &rsq);
    esl_sq_FetchFromMSA(msae, i, &esq);
    //printf("rsq %d %s\n", i, rsq->seq);  
    //printf("esq %d %s\n", i, esq->seq); 
    if (rsq->n != esq->n) printf("wrong comparison of sqs rsqlen = %" PRId64 " eqslen = %" PRId64 " \n", rsq->n, esq->n); 
    if (rsq->n != esq->n) ESL_XFAIL(eslFAIL, errbuf, "wrong comparison of sqs rsqlen = %" PRId64 " eqslen = %" PRId64 " \n", rsq->n, esq->n);
    
    /* nhr = number of non-homologous positions in reference msa
     * nhe = number of nhr non-homologous position that are also non-homologous in the inferred msa 
     *
     * hr  = number of homologous positions in reference msa
     * he  = number of homologous positions in inferred  msa
     * hre = number of hr homologous position that are also homologous in the inferred msa
     */
    for (x = 0; x < rsq->n; x++) {
      if (islower(rsq->seq[x])) {
	nhr ++;
	if (islower(esq->seq[x])) nhe ++;
      }

      if (isupper(esq->seq[x])) he ++;
      if (isupper(rsq->seq[x])) hr ++;
      if (isupper(rsq->seq[x]) && isupper(esq->seq[x])) hre ++;
    }
    
    esl_sq_Destroy(rsq); rsq = NULL;
    esl_sq_Destroy(esq); esq = NULL;
  }

  *ret_nhr = nhr;
  *ret_nhe = nhe;
  *ret_hr  = hr;
  *ret_he  = he;
  *ret_hre = hre;

  return eslOK;

 ERROR:
  if (rsq) esl_sq_Destroy(rsq);
  if (esq) esl_sq_Destroy(esq);
  return status;
}


int
msamanip_RemoveGapColumns(double gapthresh, ESL_MSA *msa, char *errbuf, int verbose)
{
  int     *useme = NULL;
  double   gapfreq;
  int      ngaps;
  int      apos;
  int      i;
  int      status;
  
  ESL_ALLOC(useme, sizeof(int) * msa->alen);
 
  for (apos = 0; apos < (int)msa->alen; apos++) {
    /* count the gaps in apos */
    ngaps = 0;
    for (i = 0; i < msa->nseq; i++) 
      if (esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) ngaps ++;
    
    /* apply gapthresh */   
    gapfreq = (double)ngaps / (double) msa->nseq;
    useme[apos] = (gapthresh < gapfreq)? FALSE : TRUE; 
  }
  
  if ((status = esl_msa_RemoveBrokenBasepairs(msa, errbuf, useme)) != eslOK) goto ERROR;
  if ((status = esl_msa_ColumnSubset         (msa, errbuf, useme)) != eslOK) goto ERROR;
  
  free(useme);
  return eslOK;
  
 ERROR:
  if (useme != NULL) free(useme); 
  return status;
}

 int
msamanip_RemoveFragments(float fragfrac, ESL_MSA **msa, int *ret_nfrags, int *ret_seq_cons_len)
{
  ESL_MSA *omsa;
  ESL_MSA *new = NULL;
  ESL_DSQ *dsq = NULL;
  int     *useme = NULL;
  double   flen;
  int64_t  clen;           // seq_cons length 
  int64_t  alen;
  int      len;
  int      n;
  int      i;
  int      x;
  int      nfrags;
  int      status;

  omsa = *msa; /* pointer to original msa */
  if (omsa->ngc == 0) { *ret_seq_cons_len = omsa->alen; return eslOK; }

  for (n = 0; n < omsa->ngc; n++) {
    if (strcmp(omsa->gc_tag[n], "seq_cons") == 0) {
      esl_abc_CreateDsq(omsa->abc, omsa->gc[n], &dsq);
      clen = esl_abc_dsqrlen(omsa->abc, dsq);
      alen = esl_abc_dsqlen(dsq);
      //printf("%s: len %d alen %d\n%s\n\n", omsa->gc_tag[n], (int)clen, (int)alen, omsa->gc[n]);
    }
  }
  if (clen == 0.) { printf("couldn't find 'seq_cons'\n"); status = eslFAIL; goto ERROR; } 
  flen = (double)clen * fragfrac;

  /* for each seq in msa, calculate number of residues that match the seq_cons */
  ESL_ALLOC(useme, sizeof(int) * omsa->nseq);
  for (nfrags = 0, i = 0; i < omsa->nseq; i++)  {
    len = 0;
    for (x = 1; dsq[x] != eslDSQ_SENTINEL; x++) { 
      if (esl_abc_XIsResidue(omsa->abc, dsq[x]) && esl_abc_XIsResidue(omsa->abc, omsa->ax[i][x])) len ++;
    }
    useme[i] = (len < flen) ? 0 : 1;
  }
  if ((status = esl_msa_SequenceSubset(omsa, useme, &new)) != eslOK) goto ERROR;

  *ret_seq_cons_len = clen;
  *ret_nfrags = omsa->nseq - esl_vec_ISum(useme, omsa->nseq);

 /* replace msa */
  esl_msa_Destroy(omsa);
  *msa = new;

  free(dsq);
  free(useme);
  return eslOK;

 ERROR:
  if (new) esl_msa_Destroy(new);
  if (dsq   != NULL) free(dsq);
  if (useme != NULL) free(useme); 
  return status;
}

/* Extract subset with sequences no more than idthesh similar to each other
 */
int
msamanip_SelectSubsetByID(ESL_RANDOMNESS *r, ESL_MSA **msa, float idthresh, int *ret_nremoved)
{      
  ESL_MSA   *omsa  = NULL;
  ESL_MSA   *new  = NULL;
  int       *assignment = NULL;
  int       *nin        = NULL;
  int       *useme      = NULL;
  int        nused      = 0;
  int        nc         = 0;
  int        c;
  int        nskip;
  int        i;
  int        status;

  eslx_msafile_Write(stdout, *msa, eslMSAFILE_STOCKHOLM); 
  omsa = *msa;

  ESL_ALLOC(useme, sizeof(int) * omsa->nseq);
   esl_vec_ISet(useme, omsa->nseq, 0);

  if ((status = esl_msacluster_SingleLinkage(omsa, idthresh, &assignment, &nin, &nc)) != eslOK) goto ERROR;
 
  for (c = 0; c < nc; c++)
    {
      nskip = esl_rnd_Roll(r, nin[c]); /* pick a random seq in this cluster to be the test. */
      for (i = 0; i < omsa->nseq; i++)
	if (assignment[i] == c) {
	  if (nskip == 0) {
	    nused ++;
	    useme[i] = 1;
	    break;
	  } else nskip--; 
	}
    }
  *ret_nremoved = omsa->nseq - nused;

  if ((status = esl_msa_SequenceSubset(omsa, useme, &new))  != eslOK) goto ERROR;
  if ((status = esl_msa_MinimGaps(new, NULL, "-.~", FALSE)) != eslOK) goto ERROR;

 /* replace msa */
  esl_msa_Destroy(omsa);
  *msa = new;

  free(useme);
  free(nin);
  free(assignment);

 return eslOK;

 ERROR:
  if (useme      != NULL) free(useme);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  if (new        != NULL) esl_msa_Destroy(new);
 return status;
}

/* Extract subset with nseq sequences taken at random
 */
int
msamanip_SelectSubset(ESL_RANDOMNESS  *r, int nseq, ESL_MSA **omsa, char **msafile, char *errbuf, int verbose)
{
  FILE         *msafp;
  char         *newfile = NULL;
  ESL_MSA      *new = NULL;
  ESL_MSA      *msa = *omsa;
  ESLX_MSAFILE *submsafp = NULL;
  char         *omsafile = NULL;
  char         *st;
  char         *newacc = NULL;
  int          *array = NULL;
  int          *useme = NULL;
  int           n;
  int           s;
  int           status;

  if (nseq == 0 || nseq >= msa->nseq) return eslOK;
 
  /* the newfile file with submsa */
  if (msafile) {
    omsafile = *msafile;
    if (omsafile) esl_sprintf(&newfile, "%s.sto", omsafile);
    
    /* if newfile exist, read the submsa from existing newfile */
    if (esl_FileExists(newfile)) { 
      free(omsafile);
      *msafile = newfile;
      
      status = eslx_msafile_Open(&msa->abc, newfile, NULL, eslMSAFILE_STOCKHOLM, NULL, &submsafp);
      if (status != eslOK) eslx_msafile_OpenFailure(submsafp, status);
      eslx_msafile_Read(submsafp, &new);
      eslx_msafile_Close(submsafp); 
      esl_msa_Destroy(msa);
      *omsa = new;  
      return eslOK; 
    }
  }

  if (msa) esl_sprintf(&newfile, "%s.sto", msa->name);

  /* otherwise, proceed */
  ESL_ALLOC(array, sizeof(int) * (msa->nseq+1));
  for (n = 0; n <= msa->nseq; n ++) array[n] = n;

  /* randomize array */
  if ((status = esl_vec_IShuffle(r, array, msa->nseq)) != eslOK) ESL_XFAIL(status, errbuf, "failed to randomize array");

  /* get the first 'nseq' */
  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  esl_vec_ISet(useme, msa->nseq, 0);
  for (s = 1; s <= nseq; s ++) useme[array[s]] = 1;
  
  if ((status = esl_msa_SequenceSubset(msa, useme, &new)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create msa subset");

  esl_msa_MinimGaps(new, errbuf, "-_.~", TRUE);

  /* change the accession of the msa to reflect that it is a subset of the original */
  if (msa->acc) {
    st = msa->acc; /* remove the version from the accession */
    esl_strtok(&st, ".", &newacc);
    esl_sprintf(&(new->acc), "%s.random%d.%s", newacc, nseq, st);
  }
  else
    esl_sprintf(&(new->acc), "random%d", nseq);

  /* write the submsa to file */
  if (msafile) esl_sprintf(&newfile, "%s.sto", omsafile);
  
  if ((msafp = fopen(newfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "failed to open %s for writting\n", newfile); 
  if (eslx_msafile_Write(msafp, new, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to write msa to %s", newfile);
  fclose(msafp);

  esl_msa_Destroy(msa);
  if (omsafile) free(omsafile);

  if (msafile) *msafile = newfile;
  if (omsa)    *omsa    = new;

  free(array);
  free(useme);
  return eslOK;

 ERROR:
  if (array)  free(array);
  if (useme)  free(useme);
  if (new)    esl_msa_Destroy(new);
  return status;
}


int
msamanip_SelectRandomSet(ESL_RANDOMNESS *r, ESL_MSA **msa, ESL_MSA **restmsa, int nset)
{      
  ESL_MSA   *omsa  = NULL;
  ESL_MSA   *new   = NULL;
  ESL_MSA   *rest  = NULL;
  int       *useme = NULL;
  int        nuse = 0;
  int        try;
  int        x;
  int        status;

  omsa = *msa;
  
  if (nset == 0) { 
    if (restmsa) {
      *restmsa = omsa;
    }
    else esl_msa_Destroy(omsa); 
    *msa = NULL;
    return eslOK;
  }
  if (omsa->nseq < nset) { printf("not enough sequences %d to get %d\n", omsa->nseq, nset); return eslOK; }
 
  ESL_ALLOC(useme, sizeof(int) * omsa->nseq);
  esl_vec_ISet(useme, omsa->nseq, 0);

  try = esl_rnd_Roll(r, omsa->nseq); // pick a random seq 
  useme[try] = 1;
  nuse ++;
  
  while (nuse < nset) {
    try = esl_rnd_Roll(r, omsa->nseq); // pick more random seq 
    if (useme[try] == 0) {
      useme[try] = 1;
      nuse ++;
    }
  }
  if ((status = esl_msa_SequenceSubset(omsa, useme, &new))  != eslOK) goto ERROR;
  if ((status = esl_msa_MinimGaps(new, NULL, "-.~", FALSE)) != eslOK) goto ERROR;

  if (restmsa) {
    for (x = 0; x < omsa->nseq; x ++) {
      useme[x] = abs(1-useme[x]);
     }
    
    if ((status = esl_msa_SequenceSubset(omsa, useme, &rest))  != eslOK) goto ERROR;
    if ((status = esl_msa_MinimGaps(rest, NULL, "-.~", FALSE)) != eslOK) goto ERROR;
     *restmsa = rest;
  }
  
  /* replace msa */
  esl_msa_Destroy(omsa);
  *msa = new;

  free(useme);
 
 return eslOK;

 ERROR:
  if (useme != NULL) free(useme);
  if (new   != NULL) esl_msa_Destroy(new);
  if (rest  != NULL) esl_msa_Destroy(rest);
 return status;
}

/* Extract trio of sequences a,b,c such that a is not more than idthesh similar to either b or c
 */
int
msamanip_SelectTrio(ESL_RANDOMNESS *r, ESL_MSA **msa, float idthresh1, float idthresh2)
{      
  ESL_MSA   *omsa  = NULL;
  ESL_MSA   *new  = NULL;
  int       *useme      = NULL;
  int       *order      = NULL;
  double     idinc      = 0.05;
  double     id;
  int        nused;
  int        n;
  int        x;
  int        aidx, bidx, cidx;
  int        nit = 100;
  int        it = 0;
  int        failed = TRUE;
  int        status;

  omsa = *msa;

  if (omsa->nseq < 3) { status = eslFAIL; goto ERROR; }
  
  ESL_ALLOC(useme, sizeof(int) * omsa->nseq);
  ESL_ALLOC(order, sizeof(int) * omsa->nseq);

  while (failed && it <= nit) {
    it ++;
    nused = 0;
    esl_vec_ISet(useme, omsa->nseq, 0);
    for (x = 0; x < omsa->nseq; x ++) order[x] = x;

    /* get one sequence at random */
    aidx = esl_rnd_Roll(r, omsa->nseq); 
    useme[0] = 1;
    order[aidx] = 0;
    order[0]    = aidx;
    nused = 1;
    
    /* get two sequences at about idthres from a */
    n = 0;
    while (nused < 2 && n < omsa->nseq) {
      bidx = aidx;
      while (bidx == aidx) bidx = esl_rnd_Roll(r, omsa->nseq);
      
      n ++;
      esl_dst_XPairId(omsa->abc, omsa->ax[bidx], omsa->ax[aidx], &id, NULL, NULL);
      
      if (id > idthresh1-idinc && id < idthresh1+idinc) {
	nused ++;
	useme[1] = 1;
	order[bidx] = 1;
	order[1]    = bidx;
     }
    }
    if (nused < 2) continue;
    
    n = 0;
    while(nused < 3 && n < omsa->nseq) {
      cidx = aidx;
      while (cidx == aidx || cidx == bidx) cidx = esl_rnd_Roll(r, omsa->nseq);
      
      n ++;
      esl_dst_XPairId(omsa->abc, omsa->ax[cidx], omsa->ax[aidx], &id, NULL, NULL);
      if (id > idthresh1-idinc && id < idthresh1+idinc) {
	esl_dst_XPairId(omsa->abc, omsa->ax[cidx], omsa->ax[bidx], &id, NULL, NULL);
	
	if (id > idthresh2-idinc && id < idthresh2+idinc) {
	  esl_dst_XPairId(omsa->abc, omsa->ax[cidx], omsa->ax[bidx], &id, NULL, NULL);	    
	  nused ++;
	  useme[2] = 1;
	  order[cidx] = 2;
	  order[2]    = cidx;
	}
      }
    }
    if (nused < 3) continue;
 
    failed = FALSE;
  }
  if (failed) { status = eslFAIL; goto ERROR; }

  /* replace msa */

  if ((status = reorder_msa(omsa, order, NULL))             != eslOK) goto ERROR;
  if ((status = esl_msa_SequenceSubset(omsa, useme, &new))  != eslOK) goto ERROR;
  if ((status = esl_msa_MinimGaps(new, NULL, "-.~", FALSE)) != eslOK) goto ERROR;

  esl_msa_Destroy(omsa);
  *msa = new;
  
  free(useme);
  free(order);
  
  return eslOK;
  
 ERROR:
  if (useme) free(useme);
  if (order) free(order);
  if (new)   esl_msa_Destroy(new);
  return status;
}

int
msamanip_ShuffleColums(ESL_RANDOMNESS  *r, ESL_MSA *msa, ESL_MSA **ret_shmsa, char *errbuf, int verbose)
{
  ESL_MSA *shmsa = NULL;
  int     *perm = NULL;
  int      i;
  int      n;
  int      status = eslOK;

  /* copy the original alignemt */
  shmsa = esl_msa_Clone(msa);
  if (shmsa == NULL) ESL_XFAIL(eslFAIL, errbuf, "bad allocation of shuffled msa");

  /* colums permutation */
  ESL_ALLOC(perm, sizeof(int) * (msa->alen));
  for (n = 0; n < msa->alen; n ++) perm[n] = n;
  if ((status = esl_vec_IShuffle(r, perm, msa->alen)) != eslOK) ESL_XFAIL(status, errbuf, "failed to randomize perm");
  
  /* aseq[0..nseq-1][0..alen-1] strings, or
   * ax[0..nseq-1][(0) 1..alen (alen+1)] digital seqs 
   */
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      for (i = 0; i < msa->nseq; i++) {
	for (n = 0; n < msa->alen; n++) {
	  shmsa->aseq[i][n] = msa->aseq[i][perm[n]];
	}
      }
    }
#ifdef eslAUGMENT_ALPHABET
  else
    {
      for (i = 0; i < msa->nseq; i++) {
	for (n = 1; n <= msa->alen; n++) 
	  shmsa->ax[i][n] = msa->ax[i][perm[n-1]+1];
      }
    }
#endif

  if (verbose) {
    eslx_msafile_Write(stdout, msa,   eslMSAFILE_STOCKHOLM);
    eslx_msafile_Write(stdout, shmsa, eslMSAFILE_STOCKHOLM);
 }

  *ret_shmsa = shmsa;

  free(perm);
  return status;

 ERROR:
  if (perm) free(perm);
  if (shmsa) esl_msa_Destroy(shmsa);
  return status;
}

int
msamanip_OutfileHeader(char *acc, char **ret_outheader)
{      
  char *outheader = NULL;
  char *st = NULL;

  esl_FileTail(acc, TRUE, &st);
  esl_sprintf(&outheader, "%s", st);
 
  *ret_outheader = outheader;
  free(st);

  return eslOK;
}

int
msamanip_MSALeaves(ESL_MSA **msa, int incnode)
{
  ESL_MSA *msafull = NULL;
  ESL_MSA *new     = NULL;
  int     *useme   = NULL;
  int      i;
  
  msafull = *msa;

  useme = malloc(msafull->nseq * sizeof(int));
  for (i = 0; i < msafull->nseq; i++) {
    if (!incnode && strncmp(msafull->sqname[i], "v", 1) == 0) useme[i] = FALSE; 
    else                                                      useme[i] = TRUE;
  }
  if (esl_msa_SequenceSubset(msafull, useme, &new) != eslOK)
    esl_fatal("failed to generate leaf alignment");
  if (esl_msa_MinimGaps(new, NULL, "-", FALSE) != eslOK) 
    esl_fatal("failed to generate leaf alignment");
  
  /* replace msa */
  esl_msa_Destroy(msafull);
  *msa = new;
  
  free(useme);
  return eslOK;
}

int
msamanip_DumpStats(FILE *ofp, ESL_MSA *msa, MSA_STAT mstat)
{
  fprintf(ofp, "name                 %%id     %%match     alen    avg_indel_num    avg_indel_len    max_indel_len  ancestral_len      avd_seq_len\n");
  fprintf(ofp, "%-20s %3.0f%%    %3.0f%%    %6d     %3.1f +/- %3.1f     %3.1f +/- %3.1f  %6d        %6d                %3.1f +/- %3.1f \n", 
	  msa->name, mstat.avgid, mstat.avgmatch, (int)msa->alen, mstat.avginum, mstat.stdinum,  mstat.avgilen, mstat.stdilen, mstat.maxilen, mstat.anclen, mstat.avgsqlen, mstat.stdsqlen);
  
 return eslOK;
}

int
msamanip_CStats(const ESL_ALPHABET *abc, ESL_MSA *msa, MSA_STAT *ret_mstat)
{
  MSA_STAT mstat;

  if (msa == NULL) {
    (*ret_mstat).avgid    = 0.0;
    (*ret_mstat).avgmatch = 0.0;
    (*ret_mstat).maxilen  = 0;
    (*ret_mstat).totilen  = 0;
    (*ret_mstat).totinum  = 0;
    (*ret_mstat).avginum  = 0.0;
    (*ret_mstat).stdinum  = 0.0;
    (*ret_mstat).avgilen  = 0.0;
    (*ret_mstat).stdilen  = 0.0;
    (*ret_mstat).avgsqlen = 0.0;
    (*ret_mstat).stdsqlen = 0.0;
    (*ret_mstat).anclen   = 0;
    return eslOK;
   }

  esl_dst_CAverageId   (msa->aseq, msa->nseq, 10000, &mstat.avgid);    /* 10000 is max_comparisons, before sampling kicks in */
  esl_dst_CAverageMatch(msa->aseq, msa->nseq, 10000, &mstat.avgmatch); /* 10000 is max_comparisons, before sampling kicks in */
  mstat.avgid    *= 100;
  mstat.avgmatch *= 100;
 
  calculate_Cstats(msa, &mstat.maxilen, &mstat.totilen, &mstat.totinum, &mstat.avginum, &mstat.stdinum, &mstat.avgilen, &mstat.stdilen, &mstat.avgsqlen, &mstat.stdsqlen, &mstat.anclen);

  if (ret_mstat) *ret_mstat = mstat;
  return eslOK;
}

int
msamanip_XStats(ESL_MSA *msa, MSA_STAT *ret_mstat)
{
  MSA_STAT mstat;
  
  if (msa == NULL) {
    (*ret_mstat).avgid    = 0.0;
    (*ret_mstat).avgmatch = 0.0;
    (*ret_mstat).maxilen  = 0;
    (*ret_mstat).totilen  = 0;
    (*ret_mstat).totinum  = 0;
    (*ret_mstat).avginum  = 0.0;
    (*ret_mstat).stdinum  = 0.0;
    (*ret_mstat).avgilen  = 0.0;
    (*ret_mstat).stdilen  = 0.0;
    (*ret_mstat).avgsqlen = 0.0;
    (*ret_mstat).stdsqlen = 0.0;
    (*ret_mstat).anclen   = 0;
    return eslOK;
   }

  esl_dst_XAverageId   (msa->abc, msa->ax, msa->nseq, 10000, &mstat.avgid);    /* 10000 is max_comparisons, before sampling kicks in */
  esl_dst_XAverageMatch(msa->abc, msa->ax, msa->nseq, 10000, &mstat.avgmatch); /* 10000 is max_comparisons, before sampling kicks in */
  mstat.avgid    *= 100;
  mstat.avgmatch *= 100;
 
  calculate_Xstats(msa, &mstat.maxilen, &mstat.totilen, &mstat.totinum, &mstat.avginum, &mstat.stdinum, &mstat.avgilen, &mstat.stdilen, &mstat.avgsqlen, &mstat.stdsqlen, &mstat.anclen);
  
  if (ret_mstat) *ret_mstat = mstat;
  return eslOK;
}

int
msamanip_CBaseComp(const ESL_ALPHABET *abc, ESL_MSA *msa, float *prior, float **ret_frq)
{
  char    *aseq;
  float   *frq = NULL;
  float    a = 100.0;            // dirichlet constant
  int      K = abc->K;
  int      n;
  int      i;
  int      status;

  ESL_ALLOC(frq, sizeof(float) * K);
  esl_vec_FSet(frq, K, 0.0);

  for (n = 0; n < msa->nseq; n ++) {
    aseq = msa->aseq[n];

    for (i = 0; aseq[i] != '\0'; i++) { 
      if (esl_abc_CIsCanonical(abc, aseq[i])) {
	frq[esl_abc_DigitizeSymbol(abc, aseq[i])] ++;
      }
    }
  }
  
  /* add prior */
  if (prior) esl_vec_FAddScaled(frq, prior, a, K);
  else       esl_vec_FIncrement(frq, K, 1.0);            // laplace prior

  /* normalize */
  esl_vec_FNorm(frq, K);
  
  *ret_frq = frq;
  return eslOK;

 ERROR:
  if (frq) free(frq);
  return status;
}

int
msamanip_XBaseComp(ESL_MSA *msa, float *prior, float **ret_frq)
{
  ESL_DSQ *ax;
  float   *frq = NULL;
  float    a = 100.0;            // dirichlet constant
  int      K = msa->abc->K;
  int      n;
  int      i;
  int      status;

  ESL_ALLOC(frq, sizeof(float) * K);
  esl_vec_FSet(frq, K, 0.0);

  for (n = 0; n < msa->nseq; n ++) {
    ax = msa->ax[n];

    for (i = 1; ax[i] != eslDSQ_SENTINEL; i++) 
      if (esl_abc_XIsCanonical(msa->abc, ax[i]))
	{
	  frq[ax[i]] ++;
	}
  }

  /* add prior */
  if (prior) esl_vec_FAddScaled(frq, prior, a, K);
  else       esl_vec_FIncrement(frq, K, 1.0);            // laplace prior

  /* normalize */
  esl_vec_FNorm(frq, K);
  
  *ret_frq = frq;
  return eslOK;

 ERROR:
  if (frq) free(frq);
  return status;
}

int
msamanip_Benchmark(FILE *benchfp, char *msaname, char *method, ESL_ALPHABET *abc, ESL_MSA *rmsa, MSA_STAT mrstat, ESL_MSA *emsa, MSA_STAT mestat, float sc, 
		   float treeavgt, float time, int lcu_r, int lcu_e, char *errbuf, int verbose)
{
  double   sen    = 0.0;
  double   ppv    = 0.0;
  double   F      = 0.0;
  double   TC     = 0.0;
  double   SPE    = 0.0;
  int      alen;
  int      trueh  = 0;
  int      foundh = 0;
  int      tph    = 0;
  int      cac    = 0; // correctly aligned columns
  int      ac     = 0; // aligned columns in reference alignment
  int      nhr    = 0; // number of nonhomologous  positions in reference
  int      nhe    = 0; // number of nonhomologous  positions in inferred
  int      hr     = 0; // number of homologous  positions in reference
  int      he     = 0; // number of homologous  positions in inferred
  int      hre    = 0; // number of homologous  positions in reference and inferred
  int      status; 
  
  status = FastSP_Run(rmsa, emsa, &tph, &trueh, &foundh, &cac, &ac, &sen, &ppv, &F, &TC, lcu_r, lcu_e, errbuf, verbose);
  if (status != eslOK) goto ERROR;

  alen = (emsa)? (int)emsa->alen : 0.0;

  /* nonhomologies */
  status = msamanip_NonHomologous(abc, rmsa, emsa, &nhr, &nhe, &hr, &he, &hre, errbuf);
  if (status != eslOK) goto ERROR;

  SPE = (nhr > 0)? (float)nhe/(float)nhr : 0.0;
 
  if (method) fprintf(benchfp, "%s         %10d %10d %10d %10d %10d  %2.2f %2.2f %2.2f %2.2f    %f    %f    %s    %" PRId64 "  %d  %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d\n", 
		      msaname, tph, trueh, foundh, nhe, nhr, mrstat.avgid, mestat.avgid, mrstat.avgmatch, mestat.avgmatch, treeavgt, time, method, 
		      rmsa->alen, alen, mrstat.avgsqlen, mrstat.stdsqlen, mestat.avgsqlen, mestat.stdsqlen, 
		      mrstat.avginum, mrstat.stdinum, mestat.avginum, mestat.stdinum, mrstat.avgilen, mrstat.stdilen, mestat.avgilen, mestat.stdilen, sc, hr, he, hre);
  else        fprintf(benchfp, "%s         %10d %10d %10d %10d %10d   %2.2f %2.2f %2.2f %2.2f    %f    %f    NA    %" PRId64 "  %d  %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d\n",      
		      msaname, tph, trueh, foundh, nhe, nhr, mrstat.avgid, mestat.avgid, mrstat.avgmatch, mestat.avgmatch, treeavgt, time, 
		      rmsa->alen, alen, mrstat.avgsqlen, mrstat.stdsqlen, mestat.avgsqlen, mestat.stdsqlen,  
		      mrstat.avginum, mrstat.stdinum, mestat.avginum, mestat.stdinum, mrstat.avgilen, mrstat.stdilen, mestat.avgilen, mestat.stdilen, sc, hr, he, hre);

  if (1||verbose) fprintf(stdout, "#%s         %10d %10d %10d %10d %10d |  %.4f   %.4f   %.4f   %.4f | %2.2f %2.2f %2.2f %2.2f | %.4f     %.2f | %" PRId64 "  %d  %.2f +/- %.2f %.2f +/- %.2f %.4f | %d %d %d\n",      
			  msaname, tph, trueh, foundh, nhe, nhr, sen, ppv, F, SPE, mrstat.avgid, mestat.avgid, mrstat.avgmatch, mestat.avgmatch, 
			  treeavgt, time, 
			  rmsa->alen, alen, mrstat.avgsqlen, mrstat.stdsqlen, mestat.avgsqlen, mestat.stdsqlen, sc, hr, he, hre);
   return eslOK;

 ERROR:
  return status;
}



static int
calculate_Cstats(ESL_MSA *msa, int *ret_maxilen, int *ret_totilen, int *ret_totinum, double *ret_avginum, double *ret_stdinum, double *ret_avgilen, double *ret_stdilen, double *ret_avgsqlen, 
		   double *ret_stdsqlen, int *ret_anclen) 
{
  char   *aseq;
  double  avginum = 0.;
  double  stdinum = 0.;
  double  avgilen = 0.;
  double  stdilen = 0.;
  double  avgsqlen = 0.;
  double  stdsqlen = 0.;
  int     totinum = 0;
  int     totilen = 0;
  int     maxilen = 0;
  int     inum;
  int     ilen = 0;
  int     anclen = 0;
  int     sqlen = 0;
  int     n;
  int     ni = 0;
  int     i;
  
  for (n = 0; n < msa->nseq; n ++) { 
    aseq =  msa->aseq[n];

    inum = 0;
    ilen = 0;

    for (i = 0; aseq[i] != '\0'; i++) { 
      
      /* a new insert */
      if( (i == 0                       && !isalpha(aseq[i])) ||
	  (i >  0 && isalpha(aseq[i-1]) && !isalpha(aseq[i]))   ) {
	ni ++;
	ilen ++;
      }
      /* continue insert */
       else if (!isalpha(aseq[i])) {
	ilen ++;
      }
      /* end an insert */
      else if (i > 0 && !isalpha(aseq[i-1]) && isalpha(aseq[i]) ) {
	if (ilen > maxilen) maxilen = ilen;
	inum ++;
	totilen += ilen;
	avgilen += (double)ilen;
	stdilen += (double)ilen * (double)ilen;
	//printf("ni %d ilen %d avgilen %f prod %f\n", ni, ilen, avgilen, stdilen);
	ilen = 0;
      }     
    }
   
    /* possible last case */
    if (ilen > 0) {
      if (ilen > maxilen) maxilen = ilen;
      inum ++;
      totilen += ilen;
      avgilen += (double)ilen;
      stdilen += (double)ilen * (double)ilen;
    }
    
    totinum += inum;
    avginum += (double)inum;
    stdinum += (double)inum * (double)inum;

    /* ancestral? */
    if (!strcmp(msa->sqname[n], "v0")) {
      for (i = 0; aseq[i] != '\0'; i++) { if (isalpha(aseq[i])) anclen ++; }
    }

    /* avg length */
    sqlen = 0;
    for (i = 0; aseq[i] != '\0'; i++) { if (isalpha(aseq[i])) sqlen ++; }
    avgsqlen += (double)sqlen;
    stdsqlen += (double)sqlen * (double)sqlen;
  }
  
  avginum = (msa->nseq > 0)? avginum/(double)msa->nseq : 0.;
  stdinum = (msa->nseq > 1)? sqrt(( stdinum - avginum * avginum * (double)msa->nseq) / ((double)msa->nseq-1.) ) : 0.;
  avgilen = (ni        > 0)? avgilen/(double)ni : 0.;
  stdilen = (ni        > 1)? sqrt(( stdilen - avgilen *avgilen * (double)ni) / ((double)ni-1.) ) : 0.;

  avgsqlen = (msa->nseq > 0)? avgsqlen/(double)msa->nseq : 0.;
  stdsqlen = (msa->nseq > 1)? sqrt(( stdsqlen - avgsqlen * avgsqlen * (double)msa->nseq) / ((double)msa->nseq-1.) ) : 0.;
  
  if (ret_maxilen)  *ret_maxilen  = maxilen;
  if (ret_totilen)  *ret_totilen  = totilen;
  if (ret_totinum)  *ret_totinum  = totinum;
  if (ret_avginum)  *ret_avginum  = avginum;
  if (ret_stdinum)  *ret_stdinum  = stdinum;
  if (ret_avgilen)  *ret_avgilen  = avgilen;
  if (ret_stdilen)  *ret_stdilen  = stdilen;
  if (ret_avgsqlen) *ret_avgsqlen = avgsqlen;
  if (ret_stdsqlen) *ret_stdsqlen = stdsqlen;
  if (ret_anclen)   *ret_anclen   = anclen;
  return eslOK;
}

static int
calculate_Xstats(ESL_MSA *msa, int *ret_maxilen, int *ret_totilen, int *ret_totinum, double *ret_avginum, double *ret_stdinum, double *ret_avgilen, double *ret_stdilen, double *ret_avgsqlen, double *ret_stdsqlen, int *ret_anclen) 
{
  ESL_DSQ *ax;
  double   avginum = 0.;
  double   stdinum = 0.;
  double   avgilen = 0.;
  double   stdilen = 0.;
  double   avgsqlen = 0.;
  double   stdsqlen = 0.;
  int      maxilen = 0;
  int      totinum = 0;
  int      totilen = 0;
  int      inum = 0;
  int      ilen = 0;
  int      anclen = 0;
  int      sqlen = 0;
  int      n;
  int      ni = 0;
  int      i;
  
  for (n = 0; n < msa->nseq; n ++) { 
    ax =  msa->ax[n];
 
    inum = 0.0;
    ilen = 0;

    for (i = 1; ax[i] != eslDSQ_SENTINEL; i++) { 
      
      /* a new insert */
      if( (i == 1                                            && !esl_abc_XIsCanonical(msa->abc, ax[i])) ||
	  (i >  1 && esl_abc_XIsCanonical(msa->abc, ax[i-1]) && !esl_abc_XIsCanonical(msa->abc, ax[i]))   ) {
	ni ++;
	ilen ++;
      }
     /* continue insert */
       else if (!esl_abc_XIsCanonical(msa->abc, ax[i])) {
	ilen ++;
       }
      /* end an insert */
       else if (i > 1 && !esl_abc_XIsCanonical(msa->abc, ax[i-1]) && esl_abc_XIsCanonical(msa->abc, ax[i])) {
	 if (ilen > maxilen) maxilen = ilen;
	 inum ++;
	 totilen += ilen;
	 avgilen += (double)ilen;
	 stdilen += (double)ilen * (double)ilen;
	 ilen = 0;
       }   
    }
    
    /* possible last case */
    if (ilen > 0) {
      if (ilen > maxilen) maxilen = ilen;
      inum ++;
      totilen += ilen;
      avgilen += (double)ilen;
      stdilen += (double)ilen * (double)ilen;
    }
    
    totinum += inum;
    avginum += (double)inum;
    stdinum += (double)inum * (double)inum;

    /* ancestral? */
    if (!strcmp(msa->sqname[n], "v0")) {
      for (i = 0; ax[i] != eslDSQ_SENTINEL; i++) { if (esl_abc_XIsCanonical(msa->abc, ax[i])) anclen ++; }
    }

    /* avg length */
    sqlen = 0;
    for (i = 1; ax[i] != eslDSQ_SENTINEL; i++) { if (esl_abc_XIsCanonical(msa->abc, ax[i])) sqlen ++; }
    avgsqlen += (double)sqlen;
    stdsqlen += (double)sqlen * (double)sqlen;

 }
  
  avginum = (msa->nseq > 0)? avginum/(double)msa->nseq : 0.;
  stdinum = (msa->nseq > 1)? sqrt(( stdinum - avginum * avginum * (double)msa->nseq) / ((double)msa->nseq-1.) ) : 0.;
  avgilen = (ni        > 0)? avgilen/(double)ni : 0.;
  stdilen = (ni        > 1)? sqrt(( stdilen - avgilen *avgilen * (double)ni) / ((double)ni-1.) ) : 0.;

  avgsqlen = (msa->nseq > 0)? avgsqlen/(double)msa->nseq : 0.;
  stdsqlen = (msa->nseq > 1)? sqrt(( stdsqlen - avgsqlen * avgsqlen * (double)msa->nseq) / ((double)msa->nseq-1.) ) : 0.;
  
  if (ret_maxilen)  *ret_maxilen  = maxilen;
  if (ret_totilen)  *ret_totilen  = totilen;
  if (ret_totinum)  *ret_totinum  = totinum;
  if (ret_avginum)  *ret_avginum  = avginum;
  if (ret_stdinum)  *ret_stdinum  = stdinum;
  if (ret_avgilen)  *ret_avgilen  = avgilen;
  if (ret_stdilen)  *ret_stdilen  = stdilen;
  if (ret_avgsqlen) *ret_avgsqlen = avgsqlen;
  if (ret_stdsqlen) *ret_stdsqlen = stdsqlen;
  if (ret_anclen)   *ret_anclen   = anclen;
  return eslOK;
}

/* lifted from easel/miniapps/esl-alimanip.c */
/* reorder_msa
 *                   
 * Given an array specifying a new order for the sequences in
 * the MSA, reorder it by swapping pointers.
 */
static int
reorder_msa(ESL_MSA *msa, int *order, char *errbuf)
{
  int status;
  char **tmp; 
  ESL_ALLOC(tmp, sizeof(char *) * msa->nseq);
  int i, a;

  /* contract check */
  /* 'order' must be have nseq elements, elements must be in range [0..nseq-1], no duplicates  */
  int *covered;
  ESL_ALLOC(covered, sizeof(int) * msa->nseq);
  esl_vec_ISet(covered, msa->nseq, 0);
  for(i = 0; i < msa->nseq; i++) { 
    /* printf("order[i:%4d]: %4d\n", i, order[i]);
       printf("covered[order[i:%4d]]: %4d\n", i, covered[order[i]]);
    */
    if(covered[order[i]]) ESL_FAIL(eslEINVAL, errbuf, "reorder_msa() order array has duplicate entries for i: %d\n", i);
    covered[order[i]] = 1;
  }
  free(covered);

  /* swap aseq or ax (one or the other must be non-NULL) */
  if(msa->flags & eslMSA_DIGITAL) { /* digital MSA */
    ESL_DSQ **tmp_dsq; 
    ESL_ALLOC(tmp_dsq, sizeof(ESL_DSQ *) * msa->nseq);
    for(i = 0; i < msa->nseq; i++) tmp_dsq[i] = msa->ax[i];
    for(i = 0; i < msa->nseq; i++) msa->ax[i] = tmp_dsq[order[i]];
    free(tmp_dsq);
  }
  else { /* text MSA */
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->aseq[i];
    for(i = 0; i < msa->nseq; i++) msa->aseq[i] = tmp[order[i]];
  }

  /* swap sqnames (mandatory) */
  for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqname[i];
  for(i = 0; i < msa->nseq; i++) msa->sqname[i] = tmp[order[i]];

  /* swap sqacc, if they exist */
  if(msa->sqacc != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqacc[i];
    for(i = 0; i < msa->nseq; i++) msa->sqacc[i] = tmp[order[i]];
  }

  /* swap sqdesc, if they exist */
  if(msa->sqdesc != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sqdesc[i];
    for(i = 0; i < msa->nseq; i++) msa->sqdesc[i] = tmp[order[i]];
  }

  /* swap ss, if they exist */
  if(msa->ss != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->ss[i];
    for(i = 0; i < msa->nseq; i++) msa->ss[i] = tmp[order[i]];
  }

  /* swap sa, if they exist */
  if(msa->sa != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->sa[i];
    for(i = 0; i < msa->nseq; i++) msa->sa[i] = tmp[order[i]];
  }

  /* swap pp, if they exist */
  if(msa->pp != NULL) { 
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->pp[i];
    for(i = 0; i < msa->nseq; i++) msa->pp[i] = tmp[order[i]];
  }

  /* swap gs annotation, if it exists */
  for(a = 0; a < msa->ngs; a++) {
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->gs[a][i];
    for(i = 0; i < msa->nseq; i++) msa->gs[a][i] = tmp[order[i]];
  }

  /* swap gr annotation, if it exists */
  for(a = 0; a < msa->ngr; a++) {
    for(i = 0; i < msa->nseq; i++) tmp[i] = msa->gr[a][i];
    for(i = 0; i < msa->nseq; i++) msa->gr[a][i] = tmp[order[i]];
  }
  free(tmp);
  return eslOK;

 ERROR: 
  return status;
}
