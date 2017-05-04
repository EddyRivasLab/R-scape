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

static int read_pdfcontacts(char *pdbcfile, int *revmap, int *ct, CLIST *clist, char *errbuf);
static int isnewcontact(int posi,int  posj, CLIST *clist);

int
ContactMap(char *msafile, char *gnuplot, ESL_MSA *msa, int *msamap, int *msarevmap,
	   int **ret_ct, int *ret_nbpairs, CLIST **ret_clist,
	   double cntmaxD, char *pdbfile, char *pdbcfile, char *errbuf, int verbose)
{
  char   tmpfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  FILE  *tmpfp = NULL;
  char  *cmd  = NULL;
  char  *args = NULL;
  CLIST *clist = NULL;
  int   *ct = NULL;
  int    L = msa->alen;
  int    alloc_ncnt = 5;
  int    ncnt;
  int    h;
  int    i, j;
  int    status;

  status = msamanip_CalculateCT(msa, &ct, ret_nbpairs, -1.0, errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "%s. Failed to calculate ct vector", errbuf);

  clist = CMAP_CreateCList(alloc_ncnt);
  if (clist == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to allocate clist");
  
  if (pdbfile == NULL)
    {
      h = 0;
      for (i = 0; i < L; i ++) {
	
	if (h == ncnt - 1) {
	  ncnt += alloc_ncnt;
	  ESL_REALLOC(clist->cnt,    sizeof(CNT)   * ncnt);
	  ESL_REALLOC(clist->srtcnt, sizeof(CNT *) * ncnt);
	}
	
	if (ct[i] > i) {
	  /* assign */
	  clist->cnt[h].i    = i+1;
	  clist->cnt[h].j    = ct[i]+1;
	  clist->cnt[h].posi = msamap[i]+1;
	  clist->cnt[h].posj = msamap[ct[i]]+1;
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
      if ((status = esl_tmpfile_named(tmpfile, &tmpfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create contactfile");
      fclose(tmpfp);
      
      if (RSCAPE_BIN)                         // look for the installed executable
	esl_sprintf(&cmd, "%s/pdb_parse.pl", RSCAPE_BIN);  
      else
	ESL_XFAIL(status, errbuf, "Failed to find pdb_parse.pl file\n");
      esl_sprintf(&args, "%s -D %f -W ALL %s %s %s > %s", cmd, cntmaxD, pdbfile, msafile, gnuplot, tmpfile);
      printf("%s\n", args);
      system(args);
      
      read_pdfcontacts(tmpfile, msarevmap, ct, clist, errbuf);
      
      //remove(tmpfile);
    }
  else if (pdbcfile != NULL)
    {
      read_pdfcontacts(pdbcfile, msarevmap, ct, clist, errbuf);
    }
  else ESL_XFAIL(eslFAIL, errbuf, "could not create contact map");

  clist->maxD = cntmaxD;
  *ret_ct     = ct;
  *ret_clist  = clist;
  printf("NC %d NBP %d\n", clist->ncnt, *ret_nbpairs);
  
#if 0
  for (h = 0; h < ncnt; h++) clist->srtcnt[h] = clist->cnt + h;
  if (ncnt > 1) qsort(clist->srtcnt, clist->ncnt, sizeof(CCNT *), cnt_sorted_by_sc);
#endif

  if (cmd) free(cmd);
  if (args) free(args);
  return eslOK;
  
 ERROR:
  if (cmd) free(cmd);
  if (args) free(args);
  if (ct) free(ct);
  if (clist) CMAP_FreeCList(clist);
  return status;
}

void
CMAP_FreeCList(CLIST *list)
{
  if (list == NULL) return;

  if (list->srtcnt) free(list->srtcnt);
  if (list->cnt)    free(list->cnt);
  free(list);
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

  return clist;

 ERROR:
  return NULL;
}


static int
read_pdfcontacts(char *pdbcfile, int *revmap, int *ct, CLIST *clist, char *errbuf)
{
  ESL_FILEPARSER  *efp   = NULL;
  char            *tok;
  int              nchain = 0;
  char           **mapfile = NULL;
  int              ncnt = clist->alloc_ncnt;
  int              posi, posj;
  int              i, j;
  double           D;
  int              c;
  int              h = 0;
  int              status;
  
  if (esl_fileparser_Open(pdbcfile, NULL, &efp) != eslOK)  ESL_XFAIL(eslFAIL, errbuf, "file open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  ESL_ALLOC(mapfile, sizeof(char *));
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbcfile);
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
	posi = atoi(tok)+1;
	i    = revmap[posi-1]+1;
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", mapfile[c]);
	
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", mapfile[c]);
	posj = atoi(tok)+1;
	j    = revmap[posj-1]+1;
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", mapfile[c]);
	
	if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", mapfile[c]);
	D = atof(tok);

	printf("^^ posi %d %d posj %d %d |%f\n", i, posi, j, posj, D);
	if (i > 0 && j > 0 && isnewcontact(posi, posj, clist)) {
	  
	  if (h == ncnt - 1) {
	    ncnt += clist->alloc_ncnt;
	    ESL_REALLOC(clist->cnt,    sizeof(CNT)   * ncnt);
	    ESL_REALLOC(clist->srtcnt, sizeof(CNT *) * ncnt);
	  }
	  
	  clist->cnt[h].posi = posi;
	  clist->cnt[h].posj = posj;
	  clist->cnt[h].i    = i;
	  clist->cnt[h].j    = j;
	  clist->cnt[h].D    = D;
	  clist->cnt[h].isbp = FALSE;
	  clist->ncnt        = h;
	  ct[i] = i; // ct = 0 not paired, ct[i]=i is contact ct[i]=j a base pair
	  ct[j] = j;
	  
	  printf("h %d posi %d %d posj %d %d |%f\n", h+1,
		 (int)clist->cnt[h].i, (int)clist->cnt[h].posi,
		 (int)clist->cnt[h].j, (int)clist->cnt[h].posj, clist->cnt[h].D);
	  h ++;
	}
      }
    esl_fileparser_Close(efp);
  }

  clist->ncnt = h;
  for (c = 0; c < nchain; c++) { free(mapfile[c]); }
  free(mapfile);
  return eslOK;

 ERROR:
  for (c = 0; c < nchain; c++) { if (mapfile[c]) free(mapfile[c]); }
  if (mapfile) free(mapfile);
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
