/* contactmap - funtions to create a contact map (both for RNA and protein)
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
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "contactmap.h"
#include "msamanip.h"

static int read_pdbmap(char *pdbmapfile, int L, int *msa2pdb, int *omsa2msa, int *ret_pdblen, char *errbuf);
static int read_pdbcontacts(char *pdbcfile, int *msa2pdb, int *omsa2msa, int hasss, int *ct, CLIST *clist, char *errbuf);
static int isnewcontact(int posi,int  posj, CLIST *clist);

int
ContactMap(char *pdbfile, char *msafile, char *gnuplot, ESL_MSA *msa, int *msa2omsa, int *omsa2msa, int abcisRNA,
	   int **ret_ct, int *ret_nbpairs, CLIST **ret_clist, int **ret_msa2pdb,
	   double cntmaxD, int cntmind, char *errbuf, int verbose)
{
  char   tmpcfile[16]   = "esltmpXXXXXX"; /* template for contacts*/
  char   tmpmapfile[16] = "esltmpXXXXXX"; /* template for msa2pdb map */
  char  *ss = NULL;
  FILE  *tmpfp  = NULL;
  char  *cmd  = NULL;
  char  *args = NULL;
  CLIST *clist = NULL;
  int   *ct = NULL;
  int   *msa2pdb = NULL;
  int    hasss = FALSE;
  int    L = msa->alen;
  int    ct_nbpairs;
  int    alloc_ncnt = 5;
  int    ncnt = 0;
  int    i, j;
  int    status;
  
  status = msamanip_CalculateCT(msa, &ct, &ct_nbpairs, -1.0, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "%s. Failed to calculate ct vector", errbuf);
  if (msa->ss_cons && abcisRNA) hasss = TRUE;
  
   ESL_ALLOC(msa2pdb, sizeof(int) * L);
   esl_vec_ISet(msa2pdb, L, -1);
   
   clist = CMAP_CreateCList(alloc_ncnt);
   if (clist == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to allocate clist");
   clist->mind = cntmind;

   // Start with the ss_cons structure if there is one
   if (hasss)
     {
       for (i = 0; i < L; i ++) {

	 msa2pdb[i] = i;

	 if (ncnt == clist->alloc_ncnt - 1) {
	   clist->alloc_ncnt += alloc_ncnt;
	   ESL_REALLOC(clist->cnt,    sizeof(CNT)   * clist->alloc_ncnt);
	   ESL_REALLOC(clist->srtcnt, sizeof(CNT *) * clist->alloc_ncnt);
	 }
	 
	 if (ct[i+1] > i+1 && (ct[i+1] - i) >= cntmind) {
	   /* assign */
	   clist->cnt[ncnt].i      = i+1;
	   clist->cnt[ncnt].j      = ct[i+1];
	   clist->cnt[ncnt].posi   = msa2omsa[i]+1;
	   clist->cnt[ncnt].posj   = msa2omsa[ct[i+1]-1]+1;
	   clist->cnt[ncnt].pdbi   = -1;
	   clist->cnt[ncnt].pdbj   = -1;
	   clist->cnt[ncnt].isbp   = TRUE;
	   clist->cnt[ncnt].bptype = WWc;
	   clist->cnt[ncnt].D      = +eslINFINITY;
	   clist->cnt[ncnt].sc     = -eslINFINITY;
	   clist->nbps ++;
	   clist->nwwc ++;
	   ncnt ++;
	 }
       }
       clist->ncnt = ncnt;
    }
   if (verbose) CMAP_Dump(stdout, clist);
   
   // Add the pdb annotation if there is one
   if (pdbfile != NULL)
     {
       if ((status = esl_tmpfile_named(tmpmapfile, &tmpfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create pdbmapfile");
       fclose(tmpfp);
       if ((status = esl_tmpfile_named(tmpcfile,   &tmpfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create pdbcfile");
       fclose(tmpfp);
       
       // read the contact from the pdbfile
       if (RSCAPE_BIN) esl_sprintf(&cmd, "%s/pdb_parse.pl", RSCAPE_BIN);  
       else            ESL_XFAIL(status, errbuf, "Failed to find program pdb_parse.pl\n");
       if (abcisRNA)  // run rnaview as well
	 esl_sprintf(&args, "%s -D %f -L %d -W MIN -C %s -M %s -R %s %s %s %s > /dev/null",
		     cmd, cntmaxD, cntmind, tmpmapfile, tmpcfile, pdbfile, msafile, RSCAPE_BIN, gnuplot);
       else 
	 esl_sprintf(&args, "%s -D %f -L %d -W MIN -C %s -M %s %s %s %s %s > /dev/null",
		     cmd, cntmaxD, cntmind, tmpmapfile, tmpcfile, pdbfile, msafile, RSCAPE_BIN, gnuplot);
       printf("%s\n", args);
       status = system(args);
       if (status == -1) ESL_XFAIL(status, errbuf, "Failed to run pdb_parse.pl\n");

       status = read_pdbmap(tmpmapfile, L, msa2pdb, omsa2msa, &(clist->pdblen), errbuf);
       if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "%s. Failed reading pdbmap", errbuf);
       remove(tmpmapfile);
       
       status = read_pdbcontacts(tmpcfile, msa2pdb, omsa2msa, hasss, ct, clist, errbuf);
       if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "%s. Failed reading contacts", errbuf);
       remove(tmpcfile);
     }
    
#if 0
   for (h = 0; h < ncnt; h++) clist->srtcnt[h] = clist->cnt + h;
   if (ncnt > 1) qsort(clist->srtcnt, clist->ncnt, sizeof(CCNT *), cnt_sorted_by_sc);
#endif

   if (abcisRNA) {
     /* Impose the new ct on the msa GC line 'cons_ss' */
     ESL_ALLOC(ss, sizeof(char) * (msa->alen+1));
     esl_ct2wuss(ct, msa->alen, ss);
     /* Replace the 'SS_cons' GC line with the new ss */
     esl_sprintf(&(msa->ss_cons), "%s", ss);
   }
   
   clist->maxD  = cntmaxD;
   *ret_nbpairs = clist->nbps;
   
   if (1||verbose) CMAP_Dump(stdout, clist);

   if (ret_ct)      *ret_ct      = ct;      else free(ct);
   if (ret_clist)   *ret_clist   = clist;   else free(clist);
   if (ret_msa2pdb) *ret_msa2pdb = msa2pdb; else free(msa2pdb);
   
   free(ss);   
   if (cmd) free(cmd);
   if (args) free(args);
   return eslOK;
   
 ERROR:
   if (ss) free(ss);
   if (msa2pdb) free(msa2pdb);
   if (cmd) free(cmd);
   if (args) free(args);
   if (ct) free(ct);
   if (clist) CMAP_FreeCList(clist);
   return status;
}

CLIST *
CMAP_CreateCList(int alloc_ncnt)
{
  CLIST *clist = NULL;
  int    status;
  
  ESL_ALLOC(clist,         sizeof(CLIST));
  ESL_ALLOC(clist->cnt,    sizeof(CNT)   * alloc_ncnt);
  ESL_ALLOC(clist->srtcnt, sizeof(CNT *) * alloc_ncnt);
  clist->srtcnt[0] = clist->cnt;

  clist->alloc_ncnt = alloc_ncnt;
  clist->ncnt   = 0;
  clist->nbps   = 0;
  clist->nwwc   = 0;
  clist->maxD   = -1;
  clist->mind   = -1;
  clist->pdblen = -1;

  return clist;

 ERROR:
  return NULL;
}

void
CMAP_FreeCList(CLIST *list)
{
  if (list == NULL) return;

  if (list->srtcnt) free(list->srtcnt);
  if (list->cnt)    free(list->cnt);
  free(list);
}

BPTYPE 
CMAP_GetBPTYPE(int i, int j, CLIST *clist)
{
  BPTYPE bptype = BPNONE;
  int    h;
  int    ii = (i<j)? i : j;
  int    jj = (i<j)? j : i;
  
  for (h = 0; h < clist->ncnt; h ++) {
    if (ii == clist->cnt[h].i && jj == clist->cnt[h].j) return clist->cnt[h].bptype;
  }
  return bptype;
}

int
CMAP_IsContactLocal(int i, int j, CLIST *clist)
{
  int h;
  int ii = (i<j)? i : j;
  int jj = (i<j)? j : i;
  
  for (h = 0; h < clist->ncnt; h ++) {
    if (ii == clist->cnt[h].i && jj == clist->cnt[h].j) return TRUE;
  }
  return FALSE;
}

int
CMAP_IsBPLocal(int i, int j, CLIST *clist)
{
  int h;
  int ii = (i<j)? i : j;
  int jj = (i<j)? j : i;
  
  for (h = 0; h < clist->ncnt; h ++) {
    if (ii == clist->cnt[h].i && jj == clist->cnt[h].j && clist->cnt[h].bptype < STACKED) return TRUE;
  }
  return FALSE;
}

int
CMAP_IsWCLocal(int i, int j, CLIST *clist)
{
  int h;
  int ii = (i<j)? i : j;
  int jj = (i<j)? j : i;
  
  for (h = 0; h < clist->ncnt; h ++) {
    if (ii == clist->cnt[h].i && jj == clist->cnt[h].j && clist->cnt[h].bptype == WWc) return TRUE;
  }
  return FALSE;
}

int
CMAP_Dump(FILE *fp, CLIST *clist)
{
  int   h;
  int   nbp = 0;
  int   nwc = 0;
  char *bptype = NULL;
  int   status = eslOK;


  for (h = 0; h < clist->ncnt; h ++) {
    if (clist->cnt[h].isbp) nbp ++;
    if (clist->cnt[h].bptype == WWc) { nwc ++; }
    CMAP_BPTYPEString(&bptype, clist->cnt[h].bptype, NULL);
  
    fprintf(fp, "# %d %d | bptype %s\n", (int)clist->cnt[h].posi, (int)clist->cnt[h].posj, bptype);
    free(bptype); bptype = NULL;
  }
  if (nbp != clist->nbps) status = eslFAIL;
  if (nwc != clist->nwwc) status = eslFAIL;

  CMAP_DumpShort(fp, clist);

  if (bptype) free(bptype);
  return status;
}

int
CMAP_DumpShort(FILE *fp, CLIST *clist)
{
  fprintf(fp, "# contacts  %d (%d bpairs %d wc bpairs)\n", clist->ncnt, clist->nbps, clist->nwwc);
  fprintf(fp, "# maxD      %.2f\n", clist->maxD);
  fprintf(fp, "# mind      %.d\n",  clist->mind);
  if (clist->pdblen > 0) { fprintf(fp, "# pdblen    %.d\n",  clist->pdblen); }
  return eslOK;
}


int 
CMAP_BPTYPEString(char **ret_bptype, BPTYPE type, char *errbuf)
{
  int status;

  switch(type) {
  case WWc:       esl_sprintf(ret_bptype, "WWc");       break;
  case WWt:       esl_sprintf(ret_bptype, "WWt");       break;
  case WHc:       esl_sprintf(ret_bptype, "WHc");       break;
  case WHt:       esl_sprintf(ret_bptype, "WHt");       break;
  case WSc:       esl_sprintf(ret_bptype, "WSc");       break;
  case WSt:       esl_sprintf(ret_bptype, "WSt");       break;

  case HWc:       esl_sprintf(ret_bptype, "HWc");       break;
  case HWt:       esl_sprintf(ret_bptype, "HWt");       break;
  case HHc:       esl_sprintf(ret_bptype, "HHc");       break;
  case HHt:       esl_sprintf(ret_bptype, "HHt");       break;
  case HSc:       esl_sprintf(ret_bptype, "HSc");       break;
  case HSt:       esl_sprintf(ret_bptype, "HSt");       break;

  case SWc:       esl_sprintf(ret_bptype, "SWc");       break;
  case SWt:       esl_sprintf(ret_bptype, "SWt");       break;
  case SHc:       esl_sprintf(ret_bptype, "SHc");       break;
  case SHt:       esl_sprintf(ret_bptype, "SHt");       break;
  case SSc:       esl_sprintf(ret_bptype, "SSc");       break;
  case SSt:       esl_sprintf(ret_bptype, "SSt");       break;
    
  case STACKED:   esl_sprintf(ret_bptype, "STACKED");   break;
  case CONTACT:   esl_sprintf(ret_bptype, "CONTACT");   break;
  case BPNONE:    esl_sprintf(ret_bptype, "BPNONE");    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong BPTYPE");  break;
  }

  return eslOK;
  
 ERROR:
  return status;
}

int 
CMAP_String2BPTYPE(char *bptype, BPTYPE *ret_type, char *errbuf)
{
  BPTYPE type;
  int     status;

  if      (!esl_strcmp(bptype, "WWc"))       type = WWc;
  else if (!esl_strcmp(bptype, "WWt"))       type = WWt;
  else if (!esl_strcmp(bptype, "WHc"))       type = WHc;
  else if (!esl_strcmp(bptype, "WHt"))       type = WHt;
  else if (!esl_strcmp(bptype, "WSc"))       type = WSc;
  else if (!esl_strcmp(bptype, "WSt"))       type = WSt;
  
  else if (!esl_strcmp(bptype, "HWc"))       type = HWc;
  else if (!esl_strcmp(bptype, "HWt"))       type = HWt;
  else if (!esl_strcmp(bptype, "HHc"))       type = HHc;
  else if (!esl_strcmp(bptype, "HHt"))       type = HHt;
  else if (!esl_strcmp(bptype, "HSc"))       type = HSc;
  else if (!esl_strcmp(bptype, "HSt"))       type = HSt;
  
  else if (!esl_strcmp(bptype, "SWc"))       type = SWc;
  else if (!esl_strcmp(bptype, "SWt"))       type = SWt;
  else if (!esl_strcmp(bptype, "SHc"))       type = SHc;
  else if (!esl_strcmp(bptype, "SHt"))       type = SHt;
  else if (!esl_strcmp(bptype, "SSc"))       type = SSc;
  else if (!esl_strcmp(bptype, "SSt"))       type = SSt;
  
  else if (!esl_strcmp(bptype, "STACKED"))   type = STACKED;
  else if (!esl_strcmp(bptype, "CONTACT"))   type = CONTACT;
  else if (!esl_strcmp(bptype, "BPNONE"))    type = BPNONE;
  else ESL_XFAIL(eslFAIL, errbuf, "wrong BYTYPE %s", bptype);

  *ret_type = type;
  return eslOK;
  
 ERROR:
  return status;
}


static int
read_pdbmap(char *pdbmapfile, int L, int *msa2pdb, int *omsa2msa, int *ret_pdblen, char *errbuf)
{
  ESL_FILEPARSER  *efp   = NULL;
  char            *tok;
  int              nchain = 0;
  char           **mapfile = NULL;
  int              posi;
  int              pdbi;
  int              pdb_min = L;
  int              pdb_max = 0;
  int              i;
  int              c;
  int              status;

  // Read the names of the .cor files to use
  if (esl_fileparser_Open(pdbmapfile, NULL, &efp) != eslOK)  ESL_XFAIL(eslFAIL, errbuf, "file open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  ESL_ALLOC(mapfile, sizeof(char *));
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbmapfile);
       esl_sprintf(&mapfile[nchain], tok);
 
      nchain ++;
      ESL_REALLOC(mapfile, sizeof(char *)*(nchain+1));
    }
  esl_fileparser_Close(efp);
  
  // Add from all chains if unique
  for (c = 0; c < nchain; c ++) {
    if (esl_fileparser_Open(mapfile[c], NULL, &efp) != eslOK)  ESL_XFAIL(eslFAIL, errbuf, "file open failed");
    esl_fileparser_SetCommentChar(efp, '#');
    while (esl_fileparser_NextLine(efp) == eslOK)
      {
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", mapfile[c]);
	pdbi = atoi(tok);
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", mapfile[c]);
	posi = atoi(tok);
	
	if (posi > 0) {
	  if (pdbi < pdb_min) pdb_min = pdbi;
	  if (pdbi > pdb_max) pdb_max = pdbi;
	  i = omsa2msa[posi-1]+1;
	  msa2pdb[i-1] = pdbi-1;
	}
      }
  }

  *ret_pdblen = pdb_max - pdb_min + 1;
  
  free(mapfile);
  return eslOK;

 ERROR:
  for (c = 0; c < nchain; c++) { if (mapfile[c]) free(mapfile[c]); }
  if (mapfile) free(mapfile);
  return status;
}

static int
read_pdbcontacts(char *pdbcfile, int *msa2pdb, int *omsa2msa, int hasss, int *ct, CLIST *clist, char *errbuf)
{
  ESL_FILEPARSER  *efp   = NULL;
  char            *tok;
  int              ncnt = clist->ncnt;
  int              alloc_ncnt = 5;
  int              posi, posj;
  int              i, j;
  BPTYPE           bptype;
  double           D;
  int              status;

  if (esl_fileparser_Open(pdbcfile, NULL, &efp) != eslOK)  ESL_XFAIL(eslFAIL, errbuf, "file open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbcfile);
      posi = atoi(tok);
      i    = omsa2msa[posi-1]+1;
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbcfile);
      
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbcfile);
      posj = atoi(tok);
      j    = omsa2msa[posj-1]+1;
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbcfile);
      
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbcfile);
      CMAP_String2BPTYPE(tok, &bptype, errbuf);
 
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbcfile);
      D = atof(tok);

      if (i > 0 && j > 0 && (j-i+1) >= clist->mind && isnewcontact(posi, posj, clist)) {
	
	if (ncnt == clist->alloc_ncnt - 1) {
	  clist->alloc_ncnt += alloc_ncnt;
	  ESL_REALLOC(clist->cnt,    sizeof(CNT)   * clist->alloc_ncnt);
	  ESL_REALLOC(clist->srtcnt, sizeof(CNT *) * clist->alloc_ncnt);
	}
	
	clist->cnt[ncnt].posi   = posi;
	clist->cnt[ncnt].posj   = posj;
	clist->cnt[ncnt].i      = i;
	clist->cnt[ncnt].j      = j;
	clist->cnt[ncnt].pdbi   = msa2pdb[i-1]+1;
	clist->cnt[ncnt].pdbj   = msa2pdb[j-1]+1;
	clist->cnt[ncnt].D      = D;
	clist->cnt[ncnt].bptype = bptype;
	clist->ncnt          = ncnt;
	clist->cnt[ncnt].isbp   = (bptype < STACKED)? TRUE : FALSE;
	if (bptype <  STACKED) clist->nbps ++;
	if (bptype == WWc)     clist->nwwc ++;

	if (!hasss) { // ct = 0 not paired, ct[i]=i is contact ct[i]=j a base pair
	  if (bptype == WWc) {
	    ct[i] = j; 
	    ct[j] = i;
	  }
	  else {
	    ct[i] = i; 
	    ct[j] = j;
	  }
	}
	else {
	  if (ct[i] == j && clist->cnt[ncnt].bptype == CONTACT)
	    clist->cnt[ncnt].bptype = WWc;
	  if (bptype == WWc && ct[i] == 0 && ct[j] == 0) {
	    ct[i] = j; 
	    ct[j] = i;
	  }
	}
	
#if 0
	if (bptype == WWc)
	  printf("ncnt %d posi %d %d %d posj %d %d %d |%f\n", ncnt+1,
		 (int)clist->cnt[ncnt].i, (int)clist->cnt[ncnt].posi, (int)clist->cnt[ncnt].pdbi,
		 (int)clist->cnt[ncnt].j, (int)clist->cnt[ncnt].posj, (int)clist->cnt[ncnt].pdbj, clist->cnt[ncnt].D);
#endif
	
	ncnt ++;
      }
    }
  esl_fileparser_Close(efp);
  
  clist->ncnt = ncnt;
  return eslOK;

 ERROR:
  return status;
}

static int
isnewcontact(int posi, int  posj, CLIST *clist)
{
  int h;
  for (h = 0; h < clist->ncnt; h ++) {
    if (posi == clist->cnt[h].posi && posj == clist->cnt[h].posj) return FALSE;
  }
  return TRUE;
}
