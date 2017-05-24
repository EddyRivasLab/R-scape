/* erfiles.c */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"

#include "rna_header.h"
#include "nrutil.h"
#include "rna.h"
#include "erfiles.h"


int er_PrintSeqs(long **chain_idx,long nchain, long num_residue, char *bseq, long **seidx, char *ChainID, long *ResSeq, char *errbuf)
{
  long c;
  long icum = 1;
  long i;
  long from, to;
  long ib, ie;
  int  status;
  
  for (c = 1; c <= nchain; c ++){ 
    ib = chain_idx[c][1];
    ie = chain_idx[c][2];

    if ((ie - ib) <= 0) continue;

    from = ResSeq[ seidx[ib][1] ];
    to   = ResSeq[ seidx[ie][1] ];
    i    = 0;
    fprintf(stdout, "# RNA/DNA chain_ID\t%c\t%ld\t%ld\n", ChainID[ seidx[ib][1]], from, to);
    fprintf(stdout, "# seq_%c ", ChainID[seidx[ib][1]]);
    while (i < to-from+1) {
      fprintf(stdout, "%c", bseq[icum]);
      i ++;
      icum ++;
    }
    fprintf(stdout, "\n");
  }

  if (icum != num_residue) ESL_XFAIL(eslFAIL, errbuf, "er_PrintSeqs() error.");
  return eslOK;
  
 ERROR:
  return status;
}
  

int er_Contacts(long num_residue, char *bseq, long **seidx, long *RY, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq,char **Miscs, double **xyz)
{
  char **seq = NULL;
  long   i;
  long   ir;
  
  fprintf(stdout, "# seq ");
  for (i = 0; i < num_residue; i ++) {
    ir = seidx[i][1];
    fprintf(stdout, "%c%ld", bseq[i], ir);
  }
  fprintf(stdout, "\n");
  
  return eslOK;
}

  
LIST *
er_CreateList(int alloc_np)
{
  LIST *list = NULL;
  int    status;
  
  ESL_ALLOC(list,         sizeof(LIST));
  ESL_ALLOC(list->pair,    sizeof(PAIR)   * alloc_np);

  list->alloc_np = alloc_np;
  list->np = 0;
  list->maxD = -1;
  list->mind = -1;

  return list;

 ERROR:
  return NULL;
}

void
er_FreeList(LIST *list)
{
  if (list == NULL) return;

   if (list->pair) free(list->pair);
  free(list);
}

int
er_ListDump(FILE *fp, LIST *list)
{
  int p;
  int nbp = 0;
  int nwc = 0;
  char *bptype = NULL;
  
  for (p = 0; p < list->np; p ++) {
    CMAP_BPTYPEString(&bptype, list->pair[p].bptype, NULL);
    
    if (list->pair[p].isbp) nbp ++;
    if (list->pair[p].bptype == WWc) nwc ++;
    
    fprintf(fp, "%ld\t%ld\t%c\t%5ld\t%c\t%c\t%5ld\t%c\t%7s\n",
	    list->pair[p].i,  list->pair[p].j,
	    list->pair[p].chi, list->pair[p].ir, list->pair[p].ic,
	    list->pair[p].jc,  list->pair[p].jr, list->pair[p].chj,
	    bptype);
  }
  fprintf(fp, "#Ncontacts %d (%d bpairs %d wc pairs)\n", list->np, nbp, nwc);
  fprintf(fp, "#maxD      %.2f\n", list->maxD);
  fprintf(fp, "#mind      %.d\n",  list->mind);
  
  free(bptype);
  return eslOK;
}

