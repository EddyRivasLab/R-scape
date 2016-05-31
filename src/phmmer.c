/*  phmmer -- function to integrate PHMMER msa output
 *
 * ER, Wed Oct 22 09:35:54 EDT 2014 [Janelia] 
 * SVN $Id:$
 */

#include "p7_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "msamanip.h"
#include "minimize.h"
#include "miscellaneous.h"
#include "phmmer.h"

struct hmmpid_data {
  char    *hmmfile;
  int        N;      // number of samples in hmmemit
  double   time;
  double   pid_target;
  double   pid_obs;
  double   firststep;
  double   tol;
  char    *errbuf;
  int      verbose;
};

static int    optimize_pack_paramvector        (double *p, long np, struct hmmpid_data *data);
static int    optimize_unpack_paramvector      (double *p, long np, struct hmmpid_data *data);
static void   optimize_bracket_define_direction(double *u, long np, struct hmmpid_data *data);
static double optimize_pid2time_func           (double *p, long np, void *dptr);
static double func_pid2time                    (char *hmmfile, double time, int N, char *errbuf);

int
PHMMER_Align(const ESL_MSA *msa, int isphmmer3, ESL_MSA **ret_phmmermsa, float *ret_sc, 
	     P7_BG *bg, int EmL, int EmN, int EvL, int EvN, int EfL, int EfN, float Eft, float popen, float pextend, char *mx,
	     int max, int dostats, char *errbuf, int verbose)
{
  ESL_MSAFILE    *afp = NULL;
  FILE            *fp = NULL;
  ESL_MSA         *phmmermsa = NULL;
  P7_HMM          *hmm =  NULL;
  ESL_SQ          *sq = NULL;
  char             tmpdbfile[16]    = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpqueryfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmptblfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpmsafile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmphmmfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char            *args = NULL;
  struct domain_s *dom = NULL;
  int              ndom;
  char            *phmmeropts = NULL;
  char            *s = NULL;
  int              hmmALEN   = -1;
  float            hmmSQLEN  = -1.0;
  float            hmmPID    = -1.0;
  float            hmmPMATCH = -1.0;
  float            hmmME     = -1.0;
  float            hmmMRE    = -1.0;
  int              bestd = 0;
  int              d;
  int              status;
  
  if (msa->nseq != 2)  ESL_XFAIL(status, errbuf, "EPHMMER_align() only works for pair of sequences for now");
  
  /* db in FASTA format */
  if ((status = esl_tmpfile_named(tmpdbfile,  &fp))   != eslOK) ESL_XFAIL(status, errbuf, "failed to create db file");
  esl_sq_FetchFromMSA(msa, 0, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write db  file\n");
  PHMMER_gethmm(sq, &hmm, bg, EmL, EmN, EvL, EvN, EfL, EfN, Eft, popen, pextend, mx, errbuf);

  fclose(fp); 
  esl_sq_Destroy(sq); sq = NULL;

  /* query in FASTA format */
  if ((status = esl_tmpfile_named(tmpqueryfile,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create query file");
  esl_sq_FetchFromMSA(msa, 1, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write query file\n");
  fclose(fp);
  esl_sq_Destroy(sq); sq = NULL;
  
  /* run PHMMER */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmptblfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create tbl file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpoutfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  fclose(fp);
  if (max) esl_sprintf(&phmmeropts, "--cpu 1 -E 10000 --incE 10000 --incdomE 10000 --max --domtblout %s", tmptblfile);
  else     esl_sprintf(&phmmeropts, "--cpu 1 -E 10000 --incE 10000 --incdomE 10000       --domtblout %s", tmptblfile);

  if (isphmmer3)  esl_sprintf(&args, "%s/lib/hmmer-3.1b1/src/phmmer     %s %s %s > %s", s, phmmeropts, tmpqueryfile, tmpdbfile, tmpoutfile);
  else            esl_sprintf(&args, "%s/lib/hmmer4/src/programs/phmmer %s %s %s > %s", s, phmmeropts, tmpqueryfile, tmpdbfile, tmpoutfile);
  system(args);
  printf("%s\n", args);

  /* parse phmmer output to get alignment */ 
  if ((status = esl_tmpfile_named(tmpmsafile, &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create msa file");
  fclose(fp);
  // report best alignment of best hit (-a) and add unaligned flanking regions (-f)
  esl_sprintf(&args, "%s/benchmarks/scripts/hmmer2afa.pl -a -f -q %s -t %s  %s %s ", s, tmpqueryfile, tmpdbfile, tmpoutfile, tmpmsafile);
  system(args);

  /* convert to msa */
  if (esl_msafile_Open(NULL, tmpmsafile, NULL, eslMSAFILE_AFA, NULL, &afp) != eslOK) ESL_XFAIL(status, errbuf, "Failed to open AFA file\n");
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  if (esl_msafile_Read(afp, &phmmermsa) != eslOK) ESL_XFAIL(status, errbuf, "Failed to read AFA file\n");
  if (status != eslOK) esl_msafile_ReadFailure(afp, status);
  esl_msafile_Close(afp);
  
  status = misc_ParseHMMERdomtbl(tmptblfile, &dom, &ndom, errbuf); if (status != eslOK) goto ERROR;
  if (verbose) misc_DomainDump(stdout, dom, ndom);
  
  /* the hmm tmpfile */
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))  != eslOK) ESL_XFAIL(status, errbuf, "failed to create tmphmm file");
  if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm))  != eslOK) ESL_XFAIL(status, errbuf, "failed to write hmm file");
  fclose(fp);   
  
  bestd = 0;
  for (d = 0; d < ndom; d++) {
    if (dostats) {
      /* ehmm stats using ehmmemeit */
      status = HMMER_EmitStats(tmphmmfile, (ESL_ALPHABET *)hmm->abc, dom[d].time, NULL, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, &hmmME, &hmmMRE, errbuf);
      if (status != eslOK) goto ERROR;
      misc_DomainWrite(stdout, &dom[d], hmmALEN, hmmSQLEN, hmmPID, hmmPMATCH, hmmME, hmmMRE);
    }
    if (dom[d].E < dom[bestd].E) bestd = d;
  }

   *ret_sc       = dom[bestd].fwd_score;
  *ret_phmmermsa = phmmermsa;
 
  remove(tmpdbfile);
  remove(tmpqueryfile);
  remove(tmpoutfile);
  remove(tmptblfile);
  remove(tmpmsafile);
  remove(tmphmmfile);
  
  if (dom) free(dom);
  if (hmm) p7_hmm_Destroy(hmm);
  if (sq) esl_sq_Destroy(sq);
  if (phmmeropts) free(phmmeropts);
  if (args) free(args);
  return eslOK;
  
 ERROR:
  remove(tmpdbfile);
  remove(tmpqueryfile);
  remove(tmpoutfile);
  remove(tmptblfile);
  remove(tmpmsafile);
  remove(tmphmmfile);
  if (dom) free(dom);
  if (sq) esl_sq_Destroy(sq);
  if (hmm) p7_hmm_Destroy(hmm);
  if (phmmeropts) free(phmmeropts);
  if (args) free(args);
  return status;  
}

int
EPHMMER_Align(const ESL_MSA *msa, float time, float fixpid, ESL_MSA **ret_phmmermsa, float *ret_usetime, float *ret_sc,
	      P7_BG *bg, int EmL, int EmN, int EvL, int EvN, int EfL, int EfN, float Eft, float popen, float pextend, char *mx, 
	      int max, int dostats, double tol, char *errbuf, int verbose)
{
  ESL_MSAFILE    *afp = NULL;
  FILE            *fp = NULL;
  ESL_MSA         *phmmermsa = NULL;
  P7_HMM          *hmm =  NULL;
  ESL_SQ          *sq = NULL;
  char             tmpdbfile[16]    = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpqueryfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char             tmptblfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpmsafile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmphmmfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  struct domain_s *dom = NULL;
  int              ndom;
  char            *args = NULL;
  char            *phmmeropts = NULL;
  char            *s = NULL;
  float            usetime;
  int              hmmALEN   = -1;
  float            hmmSQLEN  = -1.0;
  float            hmmPID    = -1.0;
  float            hmmPMATCH = -1.0;
  float            hmmME     = -1.0;
  float            hmmMRE    = -1.0;
  int              bestd;
  int              d;
  int              status;
  
  if (msa->nseq != 2)  ESL_XFAIL(status, errbuf, "EPHMMER_align() only works for pair of sequences for now");
  
  /* db in FASTA format */
  if ((status = esl_tmpfile_named(tmpdbfile,  &fp))   != eslOK) ESL_XFAIL(status, errbuf, "failed to create db file");
  esl_sq_FetchFromMSA(msa, 0, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write db  file\n");

  PHMMER_gethmm(sq, &hmm, bg, EmL, EmN, EvL, EvN, EfL, EfN, Eft, popen, pextend, mx, errbuf);
   if (fixpid > 0.) EPHMMER_pid2time(hmm, fixpid, &usetime, tol, errbuf, verbose);
   else             usetime = time;
    
  fclose(fp); 
  esl_sq_Destroy(sq); sq = NULL;
  
  /* query in FASTA format */
  if ((status = esl_tmpfile_named(tmpqueryfile,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create query file");
  esl_sq_FetchFromMSA(msa, 1, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write query file\n");
  fclose(fp);
  esl_sq_Destroy(sq); sq = NULL;  
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmptblfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create tbl file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpoutfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  fclose(fp);
  if (max) esl_sprintf(&phmmeropts, "--cpu 1 -E 10000 --incE 10000 --incdomE 10000 --max --domtblout %s", tmptblfile);
  else     esl_sprintf(&phmmeropts, "--cpu 1 -E 10000 --incE 10000 --incdomE 10000       --domtblout %s", tmptblfile);

  if (usetime < 0.0)  esl_sprintf(&args, "%s/src/programs/ephmmer              %s %s %s > %s", s, phmmeropts, tmpqueryfile, tmpdbfile, tmpoutfile);
  else                esl_sprintf(&args, "%s/src/programs/ephmmer --fixtime %f %s %s %s > %s", s, usetime, phmmeropts, tmpqueryfile, tmpdbfile, tmpoutfile);
  system(args);
  printf("%s\n", args);
  
  /* parse phmmer output to get alignment */ 
  if ((status = esl_tmpfile_named(tmpmsafile, &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create msa file");
  fclose(fp);
  // report best alignment of best hit (-a) and add unaligned flanking regions (-f)
  esl_sprintf(&args, "%s/benchmarks/scripts/hmmer2afa.pl -a -f -q %s -t %s  %s %s ", s, tmpqueryfile, tmpdbfile, tmpoutfile, tmpmsafile);
  system(args);

  status = misc_ParseHMMERedomtbl(tmptblfile, &dom, &ndom, errbuf); if (status != eslOK) goto ERROR;
  if (verbose) misc_DomainDump(stdout, dom, ndom);
  
  /* the hmm tmpfile */
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))  != eslOK) ESL_XFAIL(status, errbuf, "failed to create tmphmm file");
  if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm))  != eslOK) ESL_XFAIL(status, errbuf, "failed to write hmm file");
  fclose(fp);   
  
  bestd = 0;
  for (d = 0; d < ndom; d++) {
    if (dostats) {
      /* ehmm stats using ehmmemeit */
      status = HMMER_EmitStats(tmphmmfile, (ESL_ALPHABET *)hmm->abc, dom[d].time, NULL, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, &hmmME, &hmmMRE, errbuf);
      if (status != eslOK) goto ERROR;
      misc_DomainWrite(stdout, &dom[d], hmmALEN, hmmSQLEN, hmmPID, hmmPMATCH, hmmME, hmmMRE);
    }
    if (dom[d].E < dom[bestd].E) bestd = d;
  }
 
  /* convert to msa */
  if (esl_msafile_Open(NULL, tmpmsafile, NULL, eslMSAFILE_AFA, NULL, &afp) != eslOK) ESL_XFAIL(status, errbuf, "Failed to open AFA file\n");
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  if (esl_msafile_Read(afp, &phmmermsa) != eslOK) ESL_XFAIL(status, errbuf, "Failed to read AFA file\n");
  if (status != eslOK) esl_msafile_ReadFailure(afp, status);
  esl_msafile_Close(afp);
  
  remove(tmpdbfile);
  remove(tmpqueryfile);
  remove(tmptblfile);
  remove(tmpoutfile);
  remove(tmpmsafile);
  remove(tmphmmfile);
  
  *ret_sc        = dom[bestd].fwd_score;
  *ret_usetime   = dom[bestd].time;
  *ret_phmmermsa = phmmermsa;

  if (dom) free(dom);
  if (hmm) p7_hmm_Destroy(hmm);
  if (sq) esl_sq_Destroy(sq);
  if (phmmeropts) free(phmmeropts);
  if (args) free(args);
  return eslOK;
  
 ERROR:
  remove(tmpdbfile);
  remove(tmpqueryfile);
  remove(tmptblfile);
  remove(tmpoutfile);
  remove(tmpmsafile);
  remove(tmphmmfile);
  if (dom) free(dom);
  if (hmm) p7_hmm_Destroy(hmm);
  if (sq) esl_sq_Destroy(sq);
  if (phmmeropts) free(phmmeropts);
  if (args) free(args);
  return status;  
}


int
PHMMER_gethmm(ESL_SQ *sq, P7_HMM **ret_hmm, P7_BG *bg, int EmL, int EmN, int EvL, int EvN, int EfL, int EfN, float Eft, float popen, float pextend, char *mx, char *errbuf) 
{
  P7_BUILDER *bld  = NULL;    
  P7_HMM     *hmm = NULL;
  int         status;
  
  bld = p7_builder_Create(NULL, bg->abc);
  bld->EmL = EmL;
  bld->EmN = EmN;
  bld->EvL = EvL;
  bld->EvN = EvN;
  bld->EfL = EfL;
  bld->EfN = EfN;
  bld->Eft = Eft;
  
  /* Default is stored in the --mx option */
  status = p7_builder_LoadScoreSystem(bld, mx, popen, pextend, bg); 
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "Failed to set single query seq score system:\n%s\n", errbuf);
  
  esl_sq_Digitize(bg->abc, sq);
  p7_SingleBuilder(bld, sq, bg, &hmm, NULL, NULL, NULL); /* bypass everything else - only need HMM model */  
  
  *ret_hmm = hmm;
  p7_builder_Destroy(bld);
  return eslOK;
  
 ERROR:
  if (bld) p7_builder_Destroy(bld);
  if (hmm) p7_hmm_Destroy(hmm);
  return status;
}

int 
EPHMMER_pid2time(P7_HMM *hmm, float targetpid, float *ret_time, double tol, char *errbuf, int verbose) 
{
  struct hmmpid_data  data;
  char                tmphmmfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  FILE               *fp  = NULL;
  double             *p;	          /* parameter vector                  */
  double             *u;                  /* max initial step size vector      */
  double             *wrk; 	          /* 4 tmp vectors of length nbranches */
  double              pid_diff;
  double              firststep = 1e+0;
  int                 np;
  int                 status;
  
  if (targetpid < 0.) return eslOK;
  
  /* the hmm tmpfile */
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))  != eslOK) ESL_XFAIL(status, errbuf, "failed to create tmphmm file");
  if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm))  != eslOK) ESL_XFAIL(status, errbuf, "failed to write hmm file");
  fclose(fp);

  np = 1;     /* variable: time */
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  ESL_ALLOC(wrk, sizeof(double) * (np+1) * 4);
  
  /* Copy shared info into the "data" structure
   */
  data.N          = 500;
  data.time       = 0.4;
  data.hmmfile    = tmphmmfile;
  data.firststep  = firststep;
  data.pid_target = (double)targetpid;
  data.tol        = tol;
  data.errbuf     = errbuf;
  data.verbose    = verbose;

  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, (long)np, &data);
 
  /* pass problem to the optimizer
   */
  optimize_bracket_define_direction(u, (long)np, &data);
  status = min_Bracket(p, u, np, data.firststep,
		       &optimize_pid2time_func,
		       (void *) (&data), 
		       data.tol, wrk, &pid_diff);
  if (status != eslOK) 
    esl_fatal("ehmm_pid2time(): bad bracket minimization");	
  if (verbose) printf("END OPTIMIZATION: time %f pid %f diff %f\n", data.time, data.pid_target, pid_diff);

  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, (long)np, &data);
  if (data.verbose) printf("END PID OPTIMIZATION: pidtarget %f (found %f) diff %f time %f \n", data.pid_target, data.pid_obs, pid_diff, data.time);

  *ret_time = data.time;

  remove(tmphmmfile);
  return eslOK;

 ERROR:
  return status;
}

static int
optimize_pack_paramvector(double *p, long np, struct hmmpid_data *data)
{
  int   x = 0;
  
  p[x] = (data->time < 1.0)? log(data->time) : data->time - 1.0;

  return eslOK;  
}


static int
optimize_unpack_paramvector(double *p, long np, struct hmmpid_data *data)
{
  float time;
  float tmax = 10.0;
  int   x = 0;
  
  time = (p[x] < 0.0)? exp(p[x]) : p[x] + 1.0; 
  if (time > tmax) time = tmax;
  
  data->time = time;
  return eslOK;
}

static void
optimize_bracket_define_direction(double *u, long np, struct hmmpid_data *data)
{
  int x;
  for (x = 0; x < np; x++) u[x] = 0.25;
  u[np] = 0.25;
}

static double
optimize_pid2time_func(double *p, long np, void *dptr)
{
  struct hmmpid_data *data = (struct hmmpid_data *) dptr;
  
  optimize_unpack_paramvector(p, np, data);
  
  data->pid_obs = func_pid2time(data->hmmfile, data->time, data->N, data->errbuf);
  //printf("\nTIME %f pid %f | target pid %f\n", data->time, data->pid_obs, data->pid_target);

  return fabs(data->pid_obs - data->pid_target);
}

static double
func_pid2time(char *hmmfile, double time, int N, char *errbuf)
{
  char    tmpmsafile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char    tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  FILE   *fp  = NULL;
  char   *args = NULL;
  char   *s = NULL;
  float   pid_obs;

  esl_tmpfile_named(tmpmsafile, &fp);
  fclose(fp);
  esl_tmpfile_named(tmpoutfile, &fp);
  fclose(fp);

 /* run ehmmemit and esl-alistat to get stats on the actual hmm */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&args, "%s/src/programs/ehmmemit --time %f -N %d -a -o %s %s >/dev/null", s, time, N, tmpmsafile, hmmfile);  
  system(args);

  /* run esl-alistat */
  esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-alistat %s > %s", s, tmpmsafile, tmpoutfile);
  system(args);
  /* parse esl-alistat output file to get hmmPID */
  misc_ParseESLalistat(tmpoutfile, NULL, NULL, &pid_obs, NULL, errbuf);

  remove(tmpmsafile);
  remove(tmpoutfile);
  if (args) free(args); 
  return (double)pid_obs;
}

int
HMMER_EmitStats(char *hmmfile, ESL_ALPHABET *abc, float time, char **ret_hmmname, int *ret_hmmALEN, float *ret_hmmSQLEN, float *ret_hmmPID, float *ret_hmmPMATCH, 
		float *ret_hmmME, float *ret_hmmMRE, char *errbuf)
{ 
  char          tmpmsafile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpout1file[16] = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpout2file[16] = "esltmpXXXXXX"; /* tmpfile template */
  char         *args = NULL;
  char         *hmmname = NULL;
  char         *s = NULL;
  FILE         *fp = NULL;
  ESL_MSAFILE *afp = NULL;
  ESL_MSA      *msa = NULL;
  MSA_STAT      mstat;    
  int           hmmALEN;
  float         hmmSQLEN;
  float         hmmPID;
  float         hmmPMATCH;
  float         hmmME;
  float         hmmMRE;
  int           N = 1000;
  int           status;

  /* the tmpfiles */
  if ((status = esl_tmpfile_named(tmpmsafile,   &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create msa file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpout1file,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create out file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpout2file,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create out file");
  fclose(fp);

  /* run ehmmemit and esl-alistat to get stats on the actual hmm */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&args, "%s/src/programs/ehmmemit --time %f -N %d -a -o %s --statfile %s %s", s, time, N, tmpmsafile, tmpout1file, hmmfile);  
  system(args);

  /* run esl-alistat */
  esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-alistat %s > %s", s, tmpmsafile, tmpout2file);
  system(args);

  /* parse hmmemit output file to get hmmME and hmmMRE */
  HMMER_ParseHMMemit(tmpout1file, &hmmname, &hmmME, &hmmMRE, errbuf);

  /* msa stats */
  /* Open the MSA file */
  status = esl_msafile_Open(NULL, tmpmsafile, NULL, eslMSAFILE_STOCKHOLM, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  /* read the MSA */
  esl_msafile_Read(afp, &msa);
  msamanip_CStats(abc, msa, &mstat);
  
  misc_ParseESLalistat(tmpout2file, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, errbuf);
  
  if (ret_hmmname)   *ret_hmmname   = hmmname;
  if (ret_hmmME)     *ret_hmmME     = hmmME;
  if (ret_hmmMRE)    *ret_hmmMRE    = hmmMRE;
  if (ret_hmmALEN)   *ret_hmmALEN   = hmmALEN;
  if (ret_hmmSQLEN)  *ret_hmmSQLEN  = hmmSQLEN;
  if (ret_hmmPID)    *ret_hmmPID    = hmmPID;
  if (ret_hmmPMATCH) *ret_hmmPMATCH = mstat.avgmatch;
 
  esl_msafile_Close(afp);
  remove(tmpmsafile);
  remove(tmpout1file);
  remove(tmpout2file);
  if (args) free(args);
  esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (afp) esl_msafile_Close(afp);
  remove(tmpmsafile);
  remove(tmpout1file);
  remove(tmpout2file);
  if (args) free(args);
  if (msa) esl_msa_Destroy(msa);
  return status;
}


int
HMMER_ParseHMMemit(char *file, char **ret_hmmname, float *ret_hmmME, float *ret_hmmMRE, char *errbuf)
{
  ESL_FILEPARSER *fp  = NULL;
  char           *hmmname = NULL;
  float           ME  = -1.0;
  float           MRE = -1.0;
  char           *tok;
  float           time;
  int             status;

  if (esl_fileparser_Open(file, NULL, &fp) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to open file %s", file);
  esl_fileparser_SetCommentChar(fp, '#');
  
  while (esl_fileparser_NextLine(fp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      esl_strdup(tok, -1, &hmmname);
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      time = atof(tok);
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      ME = atof(tok);
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      MRE = atof(tok);
      if (ME < 0. || MRE < 0.) ESL_XFAIL(eslFAIL, errbuf, "error extracting ME or MRE from file %s", file);
    }
  
  *ret_hmmname = hmmname;
  *ret_hmmME   = ME;
  *ret_hmmMRE  = MRE;
  esl_fileparser_Close(fp);
  return eslOK;

 ERROR:
  if (fp) esl_fileparser_Close(fp);
  return status;
}
