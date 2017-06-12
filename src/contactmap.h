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

typedef enum {
  WWc       = 0,
  WWt       = 1,
  HHc       = 2,
  HHt       = 3,
  SSc       = 4,
  SSt       = 5,
  WHc       = 6,
  WHt       = 7,
  WSc       = 8,
  WSt       = 9,
  HSc       = 10,
  HSt       = 11,
  STACKED   = 12, // two bases in proximity that are involved in other basepairs
  CONTACT   = 13, // just in proximity in the 3D structure
  BPNONE    = 14,
} BPTYPE;

typedef struct cnt_s {
  int64_t i;        // positions in the analyzed alignment
  int64_t j;
  
  int64_t posi;     // positions in the input alignment
  int64_t posj;

  int64_t pdbi;     // positions in the pdb sequence
  int64_t pdbj;

  int     isbp;
  BPTYPE  bptype;
  
  double  D;
  
  double  sc;
} CNT;

typedef struct clist_s{
  int      alloc_ncnt;
  int      ncnt;    // total number of contacts
  int      nbps;    // total number of basepairs (all 12 types) (RNA only)
  int      nwwc;    // total number of WWc basepairs            (RNA only)

                    // conditions to annotate a pair as "contact"
  double   maxD;    // maximum Eucledian distance 
  int      mind;    // minimum backbone distance min(j-i+1) in pdb sequence
  
  CNT    **srtcnt;
  CNT     *cnt;

  int      pdblen;  // total length of the pdb fragment that is homologous to the alignment
} CLIST;


extern int    ContactMap(char *pdbfile, char *msafile, char *gnuplot, ESL_MSA *msa, int *msa2omsa, int *omsa2msa, int abcisRNA,
			 int **ret_ct, int *ret_nbpairs, CLIST **ret_clist, int **ret_msa2pdb,
			 double contD, int cntmind, char *errbuf, int verbose);
extern CLIST *CMAP_CreateCList(int alloc_ncnt);
extern void   CMAP_FreeCList(CLIST *list);
extern BPTYPE CMAP_GetBPTYPE(int i, int j, CLIST *clist);
extern int    CMAP_IsContactLocal(int i, int j, CLIST *list);
extern int    CMAP_IsBPLocal(int i, int j, CLIST *list);
extern int    CMAP_IsWCLocal(int i, int j, CLIST *clist);
extern int    CMAP_Dump(FILE *fp, CLIST *clist);
extern int    CMAP_DumpShort(FILE *fp, CLIST *clist);
extern int    CMAP_BPTYPEString(char **ret_bptype, BPTYPE type, char *errbuf);
extern int    CMAP_String2BPTYPE(char *bptype, BPTYPE *ret_type, char *errbuf);

#endif /*CONTACTMAP_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
