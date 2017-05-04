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
  int64_t i;
  int64_t j;
  
  int64_t posi;
  int64_t posj;

  int     isbp;
  double  D;
  
  double sc;
} CNT;

typedef struct clist_s{
  int      alloc_ncnt;
  int      ncnt;
  double   maxD;
  
  CNT    **srtcnt;
  CNT     *cnt;
} CLIST;


extern int    ContactMap(char *msafile, char *gnuplot, ESL_MSA *msa, int *msamap, int *msarevmap,
			 int **ret_ct, int *ret_nbpairs, CLIST **ret_clist, double contD,
			 char *pdbfile, char *pdbcfile, char *errbuf, int verbose);
extern CLIST *CMAP_CreateCList(int alloc_ncnt);
extern void   CMAP_FreeCList(CLIST *list);
extern int    CMAP_IsContactLocal(int i, int j, CLIST *list);
extern int    CMAP_IsBPLocal(int i, int j, CLIST *list);

#endif /*CONTACTMAP_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
