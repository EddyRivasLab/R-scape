/* e2hmmer_profile
 *
 *   
*/
#ifndef E2HMMER_PROFILE_INCLUDED
#define E2HMMER_PROFILE_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"	/* ESL_DSQ, ESL_ALPHABET */
#include "esl_dmatrix.h"	/* ESL_DMATRIX           */
#include "esl_getopts.h"	/* ESL_GETOPTS           */
#include "esl_tree.h"

#include "e2_config.h"
#include "e2.h"
#include "e2_profile.h"

/*******************************************************************************
 * 3. E2_PROFILE: a scoring profile for two descendant from a common ancestral
 *******************************************************************************/
/* Indices for state types
 */
enum e2hmmerp_states_e { 
  e2HP_BB = 0,
  e2HP_IB = 1,
  e2HP_SS = 2,
  e2HP_DS = 3,
  e2HP_IS = 4,
  e2HP_SD = 5,
  e2HP_DD = 6,
  e2HP_ID = 7,
  e2HP_BI = 8,
  e2HP_SI = 9,
  e2HP_DI = 10,
  e2HP_Ii = 11,
  e2HP_iI = 12,
  e2HP_EE = 13
};
#define e2HP_NSTATES 14

/* Indices for transition scores gm->tsc[] */
/* order is optimized for dynamic programming */
enum e2hmmerp_tsc_e {
  e2HP_BB_IB = 0,  // 2 to IB
  e2HP_IB_IB = 1,
  e2HP_BB_SS = 2,  // 13 to SS
  e2HP_IB_SS = 3,
  e2HP_SS_SS = 4,
  e2HP_DS_SS = 5,
  e2HP_IS_SS = 6,
  e2HP_SD_SS = 7,
  e2HP_DD_SS = 8,
  e2HP_ID_SS = 9,
  e2HP_BI_SS = 10,
  e2HP_SI_SS = 11,
  e2HP_DI_SS = 12,
  e2HP_Ii_SS = 13,
  e2HP_iI_SS = 14,
  e2HP_SS_DS = 15,  // 6 to DS
  e2HP_DS_DS = 16,
  e2HP_SD_DS = 17,
  e2HP_DD_DS = 18,
  e2HP_SI_DS = 19,
  e2HP_DI_DS = 20,
  e2HP_SS_IS = 21,  // 2 to IS
  e2HP_IS_IS = 22, 
  e2HP_SS_SD = 23,  // 6 to SD
  e2HP_DS_SD = 24,
  e2HP_IS_SD = 25,
  e2HP_SD_SD = 26,
  e2HP_DD_SD = 27,
  e2HP_ID_SD = 28,
  e2HP_SS_DD = 29, // 4 to DD
  e2HP_DS_DD = 30,
  e2HP_SD_DD = 31,
  e2HP_DD_DD = 32,
  e2HP_SD_ID = 33,  // 2 to ID
  e2HP_ID_ID = 34, 
  e2HP_BB_BI = 35,  // 2 to BI
  e2HP_BI_BI = 36,
  e2HP_SS_SI = 37,  // 2 to SI  
  e2HP_SI_SI = 38, 
  e2HP_DS_DI = 39,  // 2 to DI
  e2HP_DI_DI = 40,
  e2HP_BI_Ii = 41,  // 4 to Ii
  e2HP_SI_Ii = 42,  
  e2HP_Ii_Ii = 43,  
  e2HP_iI_Ii = 44,  
  e2HP_IB_iI = 45,  // 4 to iI
  e2HP_IS_iI = 46,  
  e2HP_Ii_iI = 47,  
  e2HP_iI_iI = 48,  
  e2HP_BB_EE = 49,  // 13 to EE
  e2HP_IB_EE = 50,  
  e2HP_SS_EE = 51,  
  e2HP_DS_EE = 52,  
  e2HP_IS_EE = 53,
  e2HP_SD_EE = 54,
  e2HP_DD_EE = 55,
  e2HP_ID_EE = 56,
  e2HP_BI_EE = 57,
  e2HP_SI_EE = 58,
  e2HP_DI_EE = 59,
  e2HP_Ii_EE = 60,
  e2HP_iI_EE = 61,
};
#define e2HP_NTRANS 62

/* Indices for sequences
 */
enum e2hmmerp_sq_e {
  e2HP_SL = 0, 
  e2HP_SR = 1
};
#define e2HP_NS 2

/* Accessing transition, emission scores */
#define e2HP_TSC(gm, s)     ((gm)->tsc[(s)])
#define e2HP_SSSC(gm, x, y) ((gm)->sssc[(x)][(y)])
#define e2HP_SLSC(gm, x)    ((gm)->ssc[(x)][e2HP_SL])
#define e2HP_SRSC(gm, y)    ((gm)->ssc[(y)][e2HP_SR])
#define e2HP_ILSC(gm, x)    ((gm)->isc[(x)][e2HP_SL])
#define e2HP_IRSC(gm, y)    ((gm)->isc[(y)][e2HP_SR])

typedef struct e2hmmer_profile_s {
  int    M;                                /* number of HMM states */

  float   **tsc;     /* transitions  [0..M][0..e2HP_NTRANS-1], hand-indexed                                */
  float  ***sssc;    /* seq1/seq2 joint subtitution emissions [0..M][0..K-1][0..K-1], hand-indexed        */
  float  ***ssc;     /* orphan subtitution emissions [0..M][0..K-1][0..e2P_NR-1], hand-indexed            */
  float  ***isc;     /* insertion emissions [0..M][0..K-1][0..e2P_NR-1], hand-indexed                     */

  int    mode;                                 /* configured algorithm mode (e.g. e2_LOCAL) */ 

  /* Info, most of which is a copy from parent HMM:                                       */
  char  *name;			/* unique name of model                                   */
  char  *acc;			/* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */

  const ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet                */

} E2HMMER_PROFILE;

/* e2hmmer_profile.c */
extern E2HMMER_PROFILE *e2hmmer_profile_Create(int M, const ESL_ALPHABET *abc);
extern int              e2hmmer_profile_Copy(const E2HMMER_PROFILE *src, E2HMMER_PROFILE *dst);
extern E2HMMER_PROFILE *e2hmmer_profile_Clone(const E2HMMER_PROFILE *gm);
extern int              e2hmmer_profile_Reuse(E2HMMER_PROFILE *gm);
extern void             e2hmmer_profile_Destroy(E2HMMER_PROFILE *gm);

#endif
