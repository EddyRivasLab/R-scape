/* rview_contacts - funtions to create a contact map (both for RNA and protein)
 *
 */
#ifndef RVIEW_CONTACTS_INCLUDED
#define RVIEW_CONTACTS_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_random.h"

#include "rview_pdbfile.h"

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

typedef struct clist_s {
  int      alloc_ncnt;
  int      ncnt;    // total number of contacts
  int      nbps;    // total number of basepairs (all 12 types) (RNA only)
  int      nwwc;    // total number of WWc basepairs            (RNA only)

                    // conditions to annotate a pair as "contact"
  double   maxD;    // maximum Eucledian distance 
  int      mind;    // minimum backbone distance min(j-i+1) in pdb sequence
  
  CNT    **srtcnt;
  CNT     *cnt;

  int      L;       // total length of the analysed alignment
  int      alen;    // total length of the input alignment
  int      pdblen;  // total length of the pdb fragment that is homologous to the alignment
} CLIST;




#endif /*RVIEW_CONTACTS_INCLUDED*/
