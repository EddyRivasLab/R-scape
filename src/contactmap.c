/*  contactmap - funtions to create a contact map (both for RNA and protein)
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

#include "contactmap.h"
#include "msamanip.h"

static int read_pdbmap(char *pdbmapfile, int L, int *msa2pdb, int *omsa2msa, char *errbuf);
static int read_pdbcontacts(char *pdbcfile, int *msa2pdb, int *omsa2msa, int *ct, CLIST *clist, char *errbuf);
static int isnewcontact(int posi,int  posj, CLIST *clist);

int
ContactMap(char *pdbfile, char *msafile, char *gnuplot, ESL_MSA *msa, int *msa2omsa, int *omsa2msa,
	   int **ret_ct, int *ret_nbpairs, CLIST **ret_clist, int **ret_msa2pdb,
	   double cntmaxD, int cntmind, char *errbuf, int verbose)
{
  char   tmpcfile[16]  = "esltmpXXXXXX"; /* template for contacts*/
  char   tmpmapfile[16] = "esltmpXXXXXX"; /* template for msa2pdb map */
  FILE  *tmpfp  = NULL;
  char  *cmd  = NULL;
  char  *args = NULL;
  CLIST *clist = NULL;
  int   *ct = NULL;
  int   *msa2pdb = NULL;
  int    L = msa->alen;
  int    alloc_ncnt = 5;
  int    ncnt = alloc_ncnt;
  int    h = 0;
  int    i, j;
  int    status;

  status = msamanip_CalculateCT(msa, &ct, ret_nbpairs, -1.0, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "%s. Failed to calculate ct vector", errbuf);

   ESL_ALLOC(msa2pdb, sizeof(int) * L);
   esl_vec_ISet(msa2pdb, L, 0);
   
   clist = CMAP_CreateCList(alloc_ncnt);
   if (clist == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to allocate clist");
  
  if (pdbfile == NULL)
    {
      clist->mind = 1; // the cntmind only affects the pbd contact maps, not the RNA structures

      for (i = 0; i < L; i ++) {
	
	msa2pdb[i] = i;
	
	if (h == ncnt - 1) {
	  ncnt += alloc_ncnt;
	  ESL_REALLOC(clist->cnt,    sizeof(CNT)   * ncnt);
	  ESL_REALLOC(clist->srtcnt, sizeof(CNT *) * ncnt);
	}
	
	if (ct[i+1] > i+1) {
	  /* assign */
	  clist->cnt[h].i    = i+1;
	  clist->cnt[h].j    = ct[i+1];
	  clist->cnt[h].posi = msa2omsa[i]+1;
	  clist->cnt[h].posj = msa2omsa[ct[i+1]-1]+1;
	  clist->cnt[h].pdbi = -1;
	  clist->cnt[h].pdbj = -1;
	  clist->cnt[h].isbp = TRUE;
	  clist->cnt[h].D    = +eslINFINITY;
	  clist->cnt[h].sc   = -eslINFINITY;
	  h ++;
	}
      }
      clist->ncnt = h;

    }
  else if (pdbfile != NULL)
    {
      clist->mind = cntmind; // the cntmind affects the pbd contact maps, not the RNA structures
      
      if ((status = esl_tmpfile_named(tmpmapfile, &tmpfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create pdbmapfile");
      fclose(tmpfp);
      if ((status = esl_tmpfile_named(tmpcfile,   &tmpfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create pdbcfile");
      fclose(tmpfp);
      
      if (RSCAPE_BIN)                         // look for the installed executable
	esl_sprintf(&cmd, "%s/pdb_parse.pl", RSCAPE_BIN);  
      else
	ESL_XFAIL(status, errbuf, "Failed to find pdb_parse.pl file\n");
      esl_sprintf(&args, "%s -D %f -W ALL -C %s -M %s %s %s %s %s > /dev/null", cmd, cntmaxD, tmpmapfile, tmpcfile, pdbfile, msafile, RSCAPE_BIN, gnuplot);
      printf("%s\n", args);
      system(args);
      
      read_pdbmap(tmpmapfile, L, msa2pdb, omsa2msa, errbuf);
      remove(tmpmapfile);
     
      read_pdbcontacts(tmpcfile, msa2pdb, omsa2msa, ct, clist, errbuf);      
      remove(tmpcfile);
    }
  else ESL_XFAIL(eslFAIL, errbuf, "could not create contact map");

#if 0
  for (h = 0; h < ncnt; h++) clist->srtcnt[h] = clist->cnt + h;
  if (ncnt > 1) qsort(clist->srtcnt, clist->ncnt, sizeof(CCNT *), cnt_sorted_by_sc);
#endif

  clist->maxD  = cntmaxD;

  if (1||verbose) CMAP_Dump(stdout, clist);
  
  if (ret_ct)      *ret_ct      = ct;      else free(ct);
  if (ret_clist)   *ret_clist   = clist;   else free(clist);
  if (ret_msa2pdb) *ret_msa2pdb = msa2pdb; else free(msa2pdb);


  if (cmd) free(cmd);
  if (args) free(args);
  return eslOK;
  
 ERROR:
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
  clist->ncnt = 0;
  clist->maxD = -1;
  clist->mind = -1;

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

int
CMAP_IsContactLocal(int i, int j, CLIST *clist)
{
  int h;
   
  for (h = 0; h < clist->ncnt; h ++) {
    if (i == clist->cnt[h].i && j == clist->cnt[h].j) return TRUE;
  }
  return FALSE;
}

int
CMAP_IsBPLocal(int i, int j, CLIST *clist)
{
  int h;
   
  for (h = 0; h < clist->ncnt; h ++) {
    if (i == clist->cnt[h].i && j == clist->cnt[h].j && clist->cnt[h].isbp) return TRUE;
  }
  return FALSE;
}

int
CMAP_Dump(FILE *fp, CLIST *clist)
{
  int h;

  fprintf(fp, "#Ncontacts %d\n",   clist->ncnt);
  fprintf(fp, "#maxD      %.2f\n", clist->maxD);
  fprintf(fp, "#mind      %.d\n",  clist->mind);
  for (h = 0; h < clist->ncnt; h ++) {
    fprintf(fp, "%d %d | isbp? %d\n", (int)clist->cnt[h].posi, (int)clist->cnt[h].posj, clist->cnt[h].isbp);
  }
  return FALSE;
}


static int
read_pdbmap(char *pdbmapfile, int L, int *msa2pdb, int *omsa2msa, char *errbuf)
{
  ESL_FILEPARSER  *efp   = NULL;
  char            *tok;
  int              nchain = 0;
  char           **mapfile = NULL;
  int              posi;
  int              pdbi;
  int              i;
  int              c;
  int              status;
  
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
  
  // Add contacts from all chains if unique
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
	  i = omsa2msa[posi-1]+1;
	  msa2pdb[i-1] = pdbi-1;
	}
      }
  }
  
  free(mapfile);
  return eslOK;

 ERROR:
  for (c = 0; c < nchain; c++) { if (mapfile[c]) free(mapfile[c]); }
  if (mapfile) free(mapfile);
  return status;
}

static int
read_pdbcontacts(char *pdbcfile, int *msa2pdb, int *omsa2msa, int *ct, CLIST *clist, char *errbuf)
{
  ESL_FILEPARSER  *efp   = NULL;
  char            *tok;
  int              ncnt = clist->alloc_ncnt;
  int              posi, posj;
  int              i, j;
  double           D;
  int              h = 0;
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
      D = atof(tok);
      
      if (i > 0 && j > 0 && j-i >= clist->mind && isnewcontact(posi, posj, clist)) {
	
	if (h == ncnt - 1) {
	  ncnt += clist->alloc_ncnt;
	  ESL_REALLOC(clist->cnt,    sizeof(CNT)   * ncnt);
	  ESL_REALLOC(clist->srtcnt, sizeof(CNT *) * ncnt);
	}
	
	clist->cnt[h].posi = posi;
	clist->cnt[h].posj = posj;
	clist->cnt[h].i    = i;
	clist->cnt[h].j    = j;
	clist->cnt[h].pdbi = msa2pdb[i-1]+1;
	clist->cnt[h].pdbj = msa2pdb[j-1]+1;
	clist->cnt[h].D    = D;
	clist->cnt[h].isbp = FALSE;
	clist->ncnt        = h;
	ct[i] = i; // ct = 0 not paired, ct[i]=i is contact ct[i]=j a base pair
	ct[j] = j;
	
#if 0
	printf("h %d posi %d %d %d posj %d %d %d |%f\n", h+1,
	       (int)clist->cnt[h].i, (int)clist->cnt[h].posi, (int)clist->cnt[h].pdbi,
	       (int)clist->cnt[h].j, (int)clist->cnt[h].posj, (int)clist->cnt[h].pdbj, clist->cnt[h].D);
#endif
	
	h ++;
      }
    }
  esl_fileparser_Close(efp);

  clist->ncnt = h;
  return eslOK;

 ERROR:
  return status;
}

static int
isnewcontact(int posi,int  posj, CLIST *clist)
{
  int h;
  for (h = 0; h < clist->ncnt; h ++) {
    if (posi == clist->cnt[h].posi && posj == clist->cnt[h].posj) return FALSE;
  }
  return TRUE;
}
