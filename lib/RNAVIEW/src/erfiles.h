/* erfiles.h
 *
 *   
*/
#ifndef ERFILES_INCLUDED
#define ERFILES_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"

#include "contactmap.h"

#include "rna_header.h"
#include "nrutil.h"
#include "rna.h"

typedef struct pair_s {
  long i;       // positions in chain
  long j;
  
  long ir;      // positions in pdb sequence
  long jr;

  char ic;    // the character
  char jc;
  char chi;   // chain
  char chj;

  int     isbp;
  BPTYPE  bptype;
  
  double  D;
  
} PAIR;

typedef struct list_s{
  int      alloc_np;
  int      np;
  double   maxD;
  int      mind; // min(j-i+1)
  
  PAIR    *pair;
} LIST;

extern int er_PrintSeqs(long **chain_idx,long nchain, long num_residue, char *bseq, long **seidx, char *ChainID, long *ResSeq, char *errbuf);
extern int er_Contacts(long num_residue, char *bseq, long **seidx, long *RY, char **AtomName,
		       char **ResName, char *ChainID, long *ResSeq,char **Miscs, double **xyz);

extern LIST *er_CreateList(int alloc_np);
extern void  er_FreeList(LIST *list);
extern int   er_ListDump(FILE *fp, LIST *list);

#endif
