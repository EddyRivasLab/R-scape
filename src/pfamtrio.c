/* pfamtrio -- select 3 way alignments such that 
 *             one of the seqs is not more that x% identical to either of the other two.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esl_getopts.h"
#include "esl_distance.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msafile.h"
#include "esl_stopwatch.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2.h"
#include "e2_msa.h"
#include "e2_train.h"
#include "e2_tree.h"
#include "msamanip.h"
#include "msatree.h"


/* Exclusive options for evolutionary model choice */

/* struct cfg_s : "Global" application configuration shared by all threads/processes.
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  /* Shared configuration in masters & workers */
  int                  argc;
  char               **argv;
  
  char                 errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS      *r;	               /* random numbers for stochastic sampling (grm-emit) */
  ESL_ALPHABET        *abc;                    /* the alphabet */
  double               fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double               idthresh1;	       /* fractional identity threshold for (a,b) or (a,c) */
  double               idthresh2;	       /* fractional identity threshold for (b,c) */
  
  char                *bmarkname;

  char                *triofile;
  FILE                *triofp;
  char                *gnuplot;

  int                  voutput;                 /* TRUE for verbose output */
  char                *outheader;              /* header for all output files */
  char                *msaheader;              /* header for all msa-specific output files */
 
  int                  nmsa;
  char                *msafile;
  ESL_MSA             *msa;
  MSA_STAT             mstat;                   /* statistics of the individual alignment */
  int                  ntrio;

  int                  infmt;
  float                tol;
  int                  verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  /* options for msa */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "--I1",           eslARG_REAL,    "0.50",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require trio to have (a,b) && (a,c) < x id",                                                0 },
  { "--I2",           eslARG_REAL,    "0.90",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require trio to have (b,c) < x id",                                                         0 },
  { "--ntrio",        eslARG_INT,      "0",      NULL,     "n>=0",   NULL,    NULL,  NULL,               "max number of trio alignments in set",                                                      0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
  /* msa format */
  { "--informat",  eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                             0 },
  /* other options */
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,       "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "pfamtrio [-options] pfamdb bmarkname";
static char banner[] = "pfam alignments of 3 sequences";

static int extract_trio(struct cfg_s *cfg, int *ret_ntrio);
static int write_msa(FILE *fp, ESL_MSA *msa, int verbose);

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, struct cfg_s *ret_cfg)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  struct cfg_s cfg;
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  cfg.argc = argc;
  cfg.argv = argv;
  
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, cfg.argv[0], banner);
      esl_usage(stdout, cfg.argv[0], usage);
      if (puts("\noptions:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/
      exit(0);
    }
  
  cfg.msafile = NULL;
  if (esl_opt_ArgNumber(go) != 2) { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  if ((cfg.msafile  = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.bmarkname  = esl_opt_GetArg(go, 2)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));
  
  /* outheader for all output files */
  cfg.outheader = NULL;
  esl_sprintf(&cfg.outheader, cfg.bmarkname);

  /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = eslx_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nmsa = 0;
  
  /* alphabet */
  cfg.abc = esl_alphabet_Create(eslAMINO);
  
  /* other options */
  cfg.fragfrac  = esl_opt_IsOn(go, "-F")? esl_opt_GetReal(go, "-F") : -1.0;
  cfg.idthresh1 = esl_opt_GetReal(go, "--I1");
  cfg.idthresh2 = esl_opt_GetReal(go, "--I2");
  cfg.tol       = esl_opt_GetReal   (go, "--tol");
  cfg.verbose   = esl_opt_GetBoolean(go, "-v");
  cfg.voutput   = esl_opt_GetBoolean(go, "--voutput");
  
  cfg.nmsa   = 0;
  cfg.msa    = NULL;
  cfg.ntrio  = esl_opt_GetInteger(go, "--ntrio");
  
  /* open outfile */
  cfg.triofile = NULL;
  cfg.triofp   = NULL;
  esl_sprintf(&cfg.triofile, "%s.trio.sto", cfg.outheader, cfg.ntrio); 
  if ((cfg.triofp = fopen(cfg.triofile, "w")) == NULL) esl_fatal("Failed to open trio file %s", cfg.triofile);
  
  *ret_go  = go;
  *ret_cfg = cfg;
  free(cfg.gnuplot);
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, cfg.argv[0], usage);
  if (puts("\nwhere options are:")                                                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 0, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", cfg.argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

int
main(int argc, char **argv)
{ 
  char           *msg = "e2train failed";
  ESL_GETOPTS    *go;
  struct cfg_s    cfg;
  ESLX_MSAFILE   *afp = NULL;
  float           frac;
  int             seq_cons_len = 0;
  int             nfrags = 0;	  	  /* # of fragments removed */
  int             ntrio = 0;
  int             status = eslOK;
  int             hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* outheader for all msa-output files */
  msamanip_OutfileHeader(cfg.msafile, &cfg.msaheader); 
  
  /* Open the MSA file */
  status = eslx_msafile_Open(&(cfg.abc), cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);
    while ((hstatus = eslx_msafile_Read(afp, &cfg.msa)) != eslEOF) {
    if (hstatus != eslOK) eslx_msafile_ReadFailure(afp, status);
    cfg.nmsa ++;
  }
  eslx_msafile_Close(afp);

  frac = 1.2 * (float)cfg.ntrio/(float)cfg.nmsa;

  /* Open the MSA file */
  status = eslx_msafile_Open(&(cfg.abc), cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);

  /* read the training MSAs */
  while ((hstatus = eslx_msafile_Read(afp, &cfg.msa)) != eslEOF) {
    if (hstatus != eslOK) eslx_msafile_ReadFailure(afp, status);
    
    esl_msa_ConvertDegen2X(cfg.msa); 
    esl_msa_Hash(cfg.msa);
    
    if (esl_opt_IsOn(go, "-F") && msamanip_RemoveFragments(cfg.fragfrac, &cfg.msa, &nfrags, &seq_cons_len) != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
     
    /* print some info */
    if (cfg.voutput) {
      msamanip_XStats(cfg.msa, &cfg.mstat); //  msa aveid and avematch 
      fprintf(stdout, "%s %6d sqs \n", cfg.msa->acc, cfg.msa->nseq);
      msamanip_DumpStats(stdout, cfg.msa, cfg.mstat); 
    }       

    /* extract a trio alignment */
    if (cfg.ntrio == 0 || ( ntrio < cfg.ntrio && esl_random(cfg.r) < frac/(frac+1.) ) ) {
      extract_trio(&cfg, &ntrio); 
    }

    esl_msa_Destroy(cfg.msa); cfg.msa = NULL;

    if (cfg.ntrio > 0 && ntrio == cfg.ntrio) break;
  }
 
  eslx_msafile_Close(afp);
  
  if (1||cfg.verbose) {
    printf("trio %d/%d\n", ntrio, cfg.nmsa);
  }
 
  /* cleanup */          
  fclose(cfg.triofp);
  free(cfg.outheader);
  free(cfg.triofile);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);

  return 0;
}


static int
extract_trio(struct cfg_s *cfg, int *ret_ntrio)
{
  double id;
  int    ntrio = *ret_ntrio;
  int    n, m;
  int    status;

  /* select tree sequences a, b, c, such that a is not more that x% 
   * identical to b and c */  
  if ((status = msamanip_SelectTrio(cfg->r, &cfg->msa, cfg->idthresh1, cfg->idthresh2)) != eslOK) return eslOK;
  ntrio ++;
  
  /* report % id */
  for (n = 0; n < cfg->msa->nseq; n ++) 
    for (m = n+1; m < cfg->msa->nseq; m ++) {
      if ((status = esl_dst_XPairId(cfg->abc, cfg->msa->ax[n], cfg->msa->ax[m], &id, NULL, NULL)) != eslOK) return status;
      printf("%d> id[%d-%d] %f\n", ntrio, n, m, id);
    }
  
  write_msa(cfg->triofp, cfg->msa, TRUE); 
  
  *ret_ntrio = ntrio;

  return eslOK;
}




static int
write_msa(FILE *fp, ESL_MSA *msa, int verbose)
{
  MSA_STAT msastat;

  if (eslx_msafile_Write(fp, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write train msa to file"); 
  if (verbose) {
    msamanip_XStats(msa, &msastat); //  msa aveid and avematch 
    msamanip_DumpStats(stdout, msa, msastat); 
  }
  
  return eslOK;
}
