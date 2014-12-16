/* phmmer
 *
 */
#ifndef PHMMER_INCLUDED
#define PHMMER_INCLUDED


#include <stdio.h>		/* FILE */

#include "hmmer.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

extern int PHMMER_Align (const ESL_MSA *msa, int isphmmer3,               ESL_MSA **ret_phmmermsa,                     float *ret_sc, P7_BG *bg, int EmL, int EmN, int EvL, int EvN, int EfL, int EfN, float Eft, float popen, float pextend, char *mx, int max, int dostats, char *errbuf, int verbose);
extern int EPHMMER_Align(const ESL_MSA *msa, float time,    float fixpid, ESL_MSA **ret_phmmermsa, float *ret_usetime, float *ret_sc, P7_BG *bg, int EmL, int EmN, int EvL, int EvN, int EfL, int EfN, float Eft, float popen, float pextend, char *mx, int max, int dostats, double tol, char *errbuf, int verbose);
extern int PHMMER_gethmm(ESL_SQ *sq, P7_HMM **ret_hmm, P7_BG *bg, int EmL, int EmN, int EvL, int EvN, int EfL, int EfN, float Eft, float popen, float pextend, char *mx, char *errbuf);
extern int EPHMMER_pid2time(P7_HMM *hmm, float targetpid, float *ret_time, double tol, char *errbuf, int verbose);
extern int HMMER_ParseHMMemit(char *file, char **ret_hmmname, float *ret_hmmME, float *ret_hmmMRE, char *errbuf);
extern int HMMER_EmitStats(char *hmmfile, ESL_ALPHABET *abc, float time, char **ret_hmmname, int *ret_hmmALEN, float *ret_hmmSQLEN, float *ret_hmmPID, float *ret_hmmPMATCH, float *ret_hmmME, float *ret_hmmMRE, char *errbuf);
#endif /*PHMMER_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
