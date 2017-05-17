/* contactmap - funtions to create a contact map (both for RNA and protein)
 *
 */
#ifndef CONTACTMAP_INCLUDED
#define CONTACTMAP_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_tree.h"

typedef struct cnt_s {
  int64_t i;        // positions in the analyzed alignment
  int64_t j;
  
  int64_t posi;     // positions in the input alignment
  int64_t posj;

  int64_t pdbi;     // positions in the pdb sequence
  int64_t pdbj;

  int     isbp;
  double  D;
  
  double sc;
} CNT;

typedef struct clist_s{
  int      alloc_ncnt;
  int      ncnt;
  double   maxD;
  int      mind; // min(j-i+1)
  
  CNT    **srtcnt;
  CNT     *cnt;
} CLIST;


extern int    ContactMap(char *pdbfile, char *msafile, char *gnuplot, ESL_MSA *msa, int *msa2omsa, int *omsa2msa,
			 int **ret_ct, int *ret_nbpairs, CLIST **ret_clist, int **ret_msa2pdb,
			 double contD, int cntmind, char *errbuf, int verbose);
extern CLIST *CMAP_CreateCList(int alloc_ncnt);
extern void   CMAP_FreeCList(CLIST *list);
extern int    CMAP_IsContactLocal(int i, int j, CLIST *list);
extern int    CMAP_IsBPLocal(int i, int j, CLIST *list);
extern int    CMAP_Dump(FILE *fp, CLIST *clist);

#endif /*CONTACTMAP_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
