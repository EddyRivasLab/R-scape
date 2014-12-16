/* msamanip
 *
 */
#ifndef MSAMANIP_INCLUDED
#define MSAMANIP_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_distance.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msafile.h"

#include "e2_tree.h"
 
typedef struct msa_stat_s {
  double avgid;
  double avgmatch;

  int    totinum;
  double avginum;
  double stdinum;
 
  int    maxilen;
  int    totilen; 
  double avgilen; 
  double stdilen;

  double avgsqlen; 
  double stdsqlen;

  int    anclen;

} MSA_STAT;


extern int msamanip_NonHomologous(ESL_ALPHABET *abc, ESL_MSA *msar, ESL_MSA *msae, int *ret_nhr, int *ret_nhe, int *ret_hr, int *ret_he, int *ret_hre, char *errbuf);
extern int msamanip_RemoveFragments(float fragfrac, ESL_MSA **msa, int *ret_nfrags, int *ret_seq_cons_len);
extern int msamanip_SelectSubset(ESL_RANDOMNESS *r, ESL_MSA **msa, float idthresh, int *ret_nremoved);
extern int msamanip_SelectTrio(ESL_RANDOMNESS *r, ESL_MSA **msa, float idthresh1, float idthresh2);
extern int msamanip_SelectRandomSet(ESL_RANDOMNESS *r, ESL_MSA **msa, ESL_MSA **restmsa, int nset);
extern int msamanip_OutfileHeader(char *acc, char **ret_outheader);
extern int msamanip_MSALeaves(ESL_MSA **msa, int incnode);
extern int msamanip_DumpStats(FILE *ofp, ESL_MSA *msa, MSA_STAT mstat);
extern int msamanip_CStats(const ESL_ALPHABET *abc, ESL_MSA *msa, MSA_STAT *ret_mstat);
extern int msamanip_XStats(ESL_MSA *msa, MSA_STAT *ret_mstat);
extern int msamanip_CBaseComp(const ESL_ALPHABET *abc, ESL_MSA *msa, float *prior, float **ret_msafreq);
extern int msamanip_XBaseComp(ESL_MSA *msa, float *prior, float **ret_msafreq);
extern int msamanip_Benchmark(FILE *benchfp, char *msaname, char *method, ESL_ALPHABET *abc, ESL_MSA *rmsa, MSA_STAT mrstat, ESL_MSA *emsa, MSA_STAT mestat, float sc, 
			      float treeavgt, float time, int lcu_r, int lcu_e, char *errbuf, int verbose);

#endif /*MSAMANIP_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
