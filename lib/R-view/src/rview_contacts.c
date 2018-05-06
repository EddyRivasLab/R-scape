/* rview_contacts - functions to infer contacts and basepairs
 *
 */
#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_fileparser.h"

#include "rview_contacts.h"
#include "rview_pdbfile.h"

int
rview_CreateContacts(int nchain, struct chain_s *chain, int *ret_nct, CLIST **ret_clist, char *errbuf, int verbose)
{
  CLIST *clist = NULL;
  int    nct = 0;
  int    status;
  
  clist = CMAP_CreateCList(ALLOC_NCT);
  if (clist == NULL) ESL_XFAIL(eslFAIL, errbuf, "Failed to allocate clist");

  *ret_nct   = nct;
  *ret_clist = clist;

  return eslOK;
  
 ERROR:
  if (clist) CMAP_FreeCList(clist);
  return status;
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
  if (clist->maxD   > 0) fprintf(fp, "# maxD      %.2f\n", clist->maxD);
  if (clist->mind   > 0) fprintf(fp, "# mind      %.d\n",  clist->mind);
  if (clist->L      > 0) fprintf(fp, "# L         %.d\n",  clist->L); 
  if (clist->alen   > 0) fprintf(fp, "# alen      %.d\n",  clist->alen); 
  if (clist->pdblen > 0) fprintf(fp, "# pdblen    %.d\n",  clist->pdblen); 
  return eslOK;
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
  clist->L      = -1;
  clist->alen   = -1;
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
CMAP_IsNewContact(int posi, int  posj, BPTYPE bptype, CLIST *clist)
{
  int h;
  for (h = 0; h < clist->ncnt; h ++) {
    if (posi == clist->cnt[h].posi && posj == clist->cnt[h].posj && clist->cnt[h].bptype == bptype) return FALSE;
  }
  return TRUE;
}


int
CMAP_ReuseCList(CLIST *clist)
{
  if (clist == NULL) return eslOK;
  
  clist->ncnt   = 0;
  clist->nbps   = 0;
  clist->nwwc   = 0;
  clist->maxD   = -1;
  clist->mind   = -1;
  clist->L      = -1;
  clist->alen   = -1;
  clist->pdblen = -1;
  return eslOK;
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
	  if (i > 0) msa2pdb[i-1] = pdbi-1;
	}
      }
      esl_fileparser_Close(efp);
      free(mapfile[c]);
      remove(mapfile[c]);
  }

  *ret_pdblen = pdb_max - pdb_min + 1;
  
  free(mapfile);
  return eslOK;

 ERROR:
  for (c = 0; c < nchain; c++) { if (mapfile[c]) free(mapfile[c]); }
  if (mapfile) free(mapfile);
  return status;
}

