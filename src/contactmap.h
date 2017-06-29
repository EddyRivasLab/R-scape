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
  WHc       = 2,
  WHt       = 3,
  WSc       = 4,
  WSt       = 5,
  
  HWc       = 6,
  HWt       = 7,
  HHc       = 8,
  HHt       = 9,
  HSc       = 10,
  HSt       = 11,
  
  SWc       = 12,
  SWt       = 13,
  SHc       = 14,
  SHt       = 15,
  SSc       = 16,
  SSt       = 17,
  
  STACKED   = 18, // two bases in proximity that are involved in other basepairs
  CONTACT   = 19, // just in proximity in the 3D structure
  BPNONE    = 20,
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
extern int    ContacMap_FromCT(CLIST *clist, int L, int *ct, int cntmind, int *msa2omsa, int *msa2pdb);
extern int    ContactMap_FromPDB(char *pdbfile, char *msafile, ESL_MSA *msa, int *omsa2msa, int abcisRNA, int *ct, CLIST *clist, int *msa2pdb,
				 double cntmaxD, int cntmind, char *errbuf, int verbose);
extern CLIST *CMAP_CreateCList(int alloc_ncnt);
extern void   CMAP_FreeCList(CLIST *list);
extern int    CMAP_ReuseCList(CLIST *list);
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
