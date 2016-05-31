/* e2train -- train parameters for e2mas (progressive alignment using an evolutionary model)
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

#define ALPHOPTS "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */

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
  
  ESL_STOPWATCH       *w;
  char                 errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS      *r;	               /* random numbers for stochastic sampling (grm-emit) */
  ESL_ALPHABET        *abc;                    /* the alphabet */
  double               fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double               idthresh;	       /* fractional identity threshold for selecting subset of sequences */

  char                *paramfile;
  char                *outfile;
  FILE                *outfp;

  int                  voutput;                 /* TRUE for verbose output */
  char                *benchfile;
  FILE                *benchfp;

  char                *gnuplot;
  char                *outheader;              /* header for all output files */
  char                *msaheader;              /* header for all msa-specific output files */
 
  int                  nmsa;
  int                  allocmsa;
  char                *msafile;
  ESL_MSA            **msalist;
  ESL_TREE           **Tlist;
  float              **msafrq;
  MSA_STAT            *mstat;                   /* statistics of the individual alignment */
  
  E2_PIPELINE         *pli;
  E1_BG               *bg;	   	       /* null model (copies made of this into threads) */
  int                  mode;                   /* e2_GLOBAL, e2_LOCAL */ 
  int                  do_viterbi;             // default is optimal accuracy

  E1_RATE             *R;
  char                *subsmx;                 /* BLUSUM62, ... */
  EMRATE              *emR;
  EVOM                 evomodel;
  struct rateparam_s   rateparam;

  char                *method;
  int                  probdist;            /* explore the distribution as a function of time */
  int                  amoeba;              /* Needle-Mead optimization */
  int                  infmt;
  float                tol;
  int                  verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  { "--probdist",     eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "probability distribution landscape",                                                        0 },
 /* options for optimization */
   { "--amoeba",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "if TRUE usel Needle-Mead optimization (default is gradient descent)",                       0 },
  /* options for msa */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "-I",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < x id",                                                               0 },
  /* Control of scoring system - substitutions */
  { "--mx",           eslARG_STRING,"BLOSUM62",  NULL,       NULL,   NULL,    NULL,  "--mxfile",         "substitution rate matrix choice (of some built-in matrices)",                               0 },
  { "--mxfile",       eslARG_INFILE,      NULL,  NULL,       NULL,   NULL,    NULL,  "--mx",             "read substitution rate matrix from file <f>",                                               0 },
  /* Control of scoring system - indels */ 
  { "--evomodel",     eslARG_STRING,    "AFG",   NULL,       NULL,   NULL,    NULL,  NULL,               "evolutionary model used",                                                                   0 },
  { "--paramfile",    eslARG_STRING,     NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "file with rate parameters (overrides the individual options below)",                        0 },
  { "--muAM",         eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of ancestral sequences",                                                   0 },
  { "--muAD",         eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of ancestral sequences",                                                   0 },
  { "--muAI",         eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of ancestral sequences",                                                   0 },
  { "--muEM",         eslARG_REAL,     "0.18",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of inserted residues",                                                     0 },
  { "--ldEM",         eslARG_REAL,     "0.176",  NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of adding a res to an insert",                                                         0 },
  { "--muED",         eslARG_REAL,     "0.18",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of inserted residues",                                                     0 },
  { "--ldED",         eslARG_REAL,     "0.176",  NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of adding a res to an insert",                                                         0 },
  { "--muI",          eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of inserted sequences",                                                    0 },
  { "--ldI",          eslARG_REAL,     "0.07",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of adding a whole insert",                                                             0 },
  { "--rI",           eslARG_REAL,     "0.89",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "fragment parameter for I",                                                                  0 },
  { "--rM",           eslARG_REAL,     "0.90",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "fragment parameter for M",                                                                  0 },
  { "--rD",           eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "fragment parameter for D",                                                                  0 },
  { "--sI",           eslARG_REAL,     "0.89",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "mic: adding an insert",                                                                     0 },
  { "--sD",           eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "mic: removing an insert",                                                                   0 },
  { "--vI",           eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "mic: adding to an existing insert",                                                         0 },
  { "--vD",           eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "mic: removing from an existing an insert",                                                  0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
 /* Selecting the alphabet rather than autoguessing it */
  { "--amino",        eslARG_NONE,      TRUE,    NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is protein sequence data",                                                  2 },
  { "--dna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is DNA sequence data",                                                      2 },
  { "--rna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is RNA sequence data",                                                      2 },
  /* msa format */
  { "--informat",  eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                             0 },
  /* other options */
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,       "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "e2train [-options] <msa>";
static char banner[] = "parameter training for e2 evolutionary model";

static int write_output (struct cfg_s *cfg);
static int run_voutput(struct cfg_s *cfg);
static int explore_distribution(struct cfg_s *cfg);
static int move_rate     (int which, int idx, int np, double *param, struct cfg_s *cfg);
static int move_bernoulli(int which, int idx, int np, double *param, struct cfg_s *cfg);

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, struct cfg_s *ret_cfg)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  struct cfg_s cfg;
  int          n;
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
  if (esl_opt_ArgNumber(go) != 1) { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
 
  if ((cfg.msafile  = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
    cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));

  /* outheader for all output files */
  cfg.outheader = NULL;
  msamanip_OutfileHeader(cfg.msafile, &cfg.outheader); 
  
   /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nmsa = 0;

  /* alphabet */
  if      (esl_opt_GetBoolean(go, "--amino"))   cfg.abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg.abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg.abc = esl_alphabet_Create(eslRNA);
  else                                          cfg.abc = NULL;
 
  cfg.w = esl_stopwatch_Create(); 

   /* other options */
  cfg.fragfrac   = esl_opt_IsOn(go, "-F")? esl_opt_GetReal(go, "-F") : -1.0;
  cfg.idthresh   = esl_opt_IsOn(go, "-I")? esl_opt_GetReal(go, "-I") : -1.0;
  cfg.tol        = esl_opt_GetReal   (go, "--tol");
  cfg.verbose    = esl_opt_GetBoolean(go, "-v");
  cfg.voutput    = esl_opt_GetBoolean(go, "--voutput");
  cfg.probdist   = esl_opt_GetBoolean(go, "--probdist"); 
  cfg.amoeba     = esl_opt_GetBoolean(go, "--amoeba"); 
  
  cfg.subsmx     = esl_opt_GetString(go, "--mx");
  cfg.mode       = e2_GLOBAL;
  cfg.do_viterbi = FALSE;

  cfg.evomodel       = e1_rate_Evomodel(esl_opt_GetString(go, "--evomodel"));
  cfg.rateparam.muAM = esl_opt_GetReal  (go, "--muAM");
  cfg.rateparam.muAD = esl_opt_GetReal  (go, "--muAD");
  cfg.rateparam.muAI = esl_opt_GetReal  (go, "--muAI");
  cfg.rateparam.muI  = esl_opt_GetReal  (go, "--muI");
  cfg.rateparam.ldI  = esl_opt_GetReal  (go, "--ldI");
  cfg.rateparam.muEM = esl_opt_GetReal  (go, "--muEM");
  cfg.rateparam.ldEM = esl_opt_GetReal  (go, "--ldEM");
  cfg.rateparam.muED = esl_opt_GetReal  (go, "--muED");
  cfg.rateparam.ldED = esl_opt_GetReal  (go, "--ldED");
  cfg.rateparam.rI   = esl_opt_GetReal  (go, "--rI");
  cfg.rateparam.rM   = esl_opt_GetReal  (go, "--rM");
  cfg.rateparam.rD   = esl_opt_GetReal  (go, "--rD");
  cfg.rateparam.sI   = esl_opt_GetReal  (go, "--sI");
  cfg.rateparam.sD   = esl_opt_GetReal  (go, "--sD");
  cfg.rateparam.vI   = esl_opt_GetReal  (go, "--vI");
  cfg.rateparam.vD   = esl_opt_GetReal  (go, "--vD");

  /* the paramfile takes precedent over individually feed values above */
  cfg.paramfile = NULL;
  if (esl_opt_IsOn(go, "--paramfile")) {
    cfg.paramfile = esl_opt_GetString(go, "--paramfile");
    status = e1_rate_ReadParamfile(cfg.paramfile, &cfg.rateparam, &cfg.evomodel, cfg.errbuf, cfg.verbose);
    if (status != eslOK) esl_fatal("Failed to read paramfile %s\n%s", cfg.paramfile, cfg.errbuf);
  }

  /* open outfile */
  cfg.outfile = NULL;
  cfg.outfp   = NULL;
  if (cfg.amoeba) esl_sprintf(&cfg.outfile, "%sNM.%s.param", cfg.outheader, e1_rate_EvomodelType(cfg.evomodel)); 
  else            esl_sprintf(&cfg.outfile, "%sGD.%s.param", cfg.outheader, e1_rate_EvomodelType(cfg.evomodel)); 
  if ((cfg.outfp = fopen(cfg.outfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outfile);
 
  cfg.R   = NULL;
  cfg.bg  = e1_bg_Create(cfg.abc);
  cfg.emR = NULL;

  cfg.allocmsa = 20;
  cfg.nmsa = 0;
  ESL_ALLOC(cfg.Tlist,   sizeof(ESL_TREE *) * cfg.allocmsa);
  ESL_ALLOC(cfg.msalist, sizeof(ESL_MSA  *) * cfg.allocmsa);
  ESL_ALLOC(cfg.msafrq,  sizeof(float    *) * cfg.allocmsa);
  for (n = 0; n < cfg.allocmsa; n ++) { 
    cfg.Tlist[n]   = NULL;
    cfg.msalist[n] = NULL;
    cfg.msafrq[n]  = NULL;
  }

  cfg.benchfile = NULL;
  cfg.benchfp   = NULL;

  esl_sprintf(&cfg.method, "e2train.%s", e1_rate_EvomodelType(cfg.evomodel)); 

  *ret_go = go;
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
  ESL_MSAFILE   *afp = NULL;
  int             n;
  int             seq_cons_len = 0;
  int             nfrags = 0;	  	  /* # of fragments removed */
  int             nremoved = 0;	          /* # of identical sequences removed */
  int             status = eslOK;
  int             hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* Open the MSA file */
  status = esl_msafile_Open(&(cfg.abc), cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  /* outheader for all msa-output files */
  msamanip_OutfileHeader(cfg.msafile, &cfg.msaheader); 
 
  /* read the training MSAs */
  while ((hstatus = esl_msafile_Read(afp, &cfg.msalist[cfg.nmsa])) != eslEOF) {
    if (hstatus != eslOK) esl_msafile_ReadFailure(afp, status);

    esl_msa_ConvertDegen2X(cfg.msalist[cfg.nmsa]); 
    esl_msa_Hash(cfg.msalist[cfg.nmsa]);
   
    if (esl_opt_IsOn(go, "-F") && msamanip_RemoveFragments(cfg.fragfrac, &cfg.msalist[cfg.nmsa], &nfrags, &seq_cons_len) != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-I")) msamanip_SelectSubsetBymaxID(cfg.r, &cfg.msalist[cfg.nmsa], cfg.idthresh, &nremoved);
    
     msamanip_XBaseComp(cfg.msalist[cfg.nmsa], cfg.bg->f, &cfg.msafrq[cfg.nmsa]);
     
     /* Initializations */
     cfg.pli = e2_pipeline_Create(go, 1, 100, 100);
     
     /* Create the evolutionary rate model */
     cfg.R = e1_rate_CreateWithValues(cfg.abc, cfg.evomodel, cfg.rateparam, cfg.subsmx, NULL, NULL, TRUE, cfg.tol, cfg.errbuf, cfg.verbose);
     if (cfg.R == NULL) { printf("Bad rate model.\n"); esl_fatal(cfg.errbuf); }
 
     /* calculate the tree */
     e2_tree_UPGMA(&cfg.Tlist[cfg.nmsa], cfg.msalist[cfg.nmsa], cfg.msafrq[cfg.nmsa], cfg.r, cfg.pli, cfg.R, NULL, cfg.bg, NULL, E2F, cfg.mode, cfg.do_viterbi, -1.0, -1.0, cfg.tol, cfg.errbuf, cfg.verbose);
     
     /* print some info */
     if (0&&cfg.voutput) {
       msamanip_XStats(cfg.msalist[cfg.nmsa], &cfg.mstat); //  msa aveid and avematch 
       fprintf(stdout, "Alignment\n");
       fprintf(stdout, "%6d          %s\n", cfg.msalist[cfg.nmsa]->nseq, cfg.msafile);
       if (esl_msafile_Write(stdout, cfg.msalist[cfg.nmsa], eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write msa to file"); 
       msamanip_DumpStats(stdout, cfg.msalist[cfg.nmsa], cfg.mstat); 
     }    
     cfg.nmsa ++;

     if (cfg.nmsa >= cfg.allocmsa) {
       cfg.allocmsa = cfg.nmsa + 1;
       ESL_REALLOC(cfg.Tlist,   sizeof(ESL_TREE *) * cfg.allocmsa);
       ESL_REALLOC(cfg.msalist, sizeof(ESL_MSA  *) * cfg.allocmsa);
       ESL_REALLOC(cfg.msafrq,  sizeof(float    *) * cfg.allocmsa);
     }
     
  } 
  if (1||cfg.verbose) printf("nMSA %d\n", cfg.nmsa);
  
  
  /* explore the parameter space */
  if (cfg.probdist) {
    status = explore_distribution(&cfg);
    if (status != eslOK)  { esl_fatal(cfg.errbuf); }
    return eslOK;
  }
  
  /* main function */
  esl_stopwatch_Start(cfg.w);
  status = e2_train(cfg.r, cfg.nmsa, cfg.Tlist, cfg.msalist, cfg.msafrq, cfg.R, cfg.pli, cfg.bg, cfg.mode, cfg.do_viterbi, cfg.amoeba, cfg.tol, cfg.errbuf, cfg.verbose);
  if (status != eslOK)  { esl_fatal(cfg.errbuf); }
  esl_stopwatch_Stop(cfg.w);

  status = write_output(&cfg);
  if (status != eslOK)  { esl_fatal(cfg.errbuf); }
 
  if (cfg.voutput) {
    status = run_voutput(&cfg);
    if (status != eslOK)  { esl_fatal(cfg.errbuf); }
  }  

  /* cleanup */
  esl_stopwatch_Destroy(cfg.w);
  if (cfg.bg)  e1_bg_Destroy(cfg.bg);          
  fclose(cfg.outfp);
  free(cfg.outheader);
  free(cfg.outfile);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  esl_msafile_Close(afp);
  if (cfg.R) e1_rate_Destroy(cfg.R);  
  e2_pipeline_Destroy(cfg.pli);
  if (cfg.emR) ratematrix_emrate_Destroy(cfg.emR, 1);
  for (n = 0; n < cfg.nmsa; n ++) {
    esl_tree_Destroy(cfg.Tlist[n]); 
    esl_msa_Destroy(cfg.msalist[n]); 
    free(cfg.msafrq[n]); 
  }
  free(cfg.Tlist);
  free(cfg.msalist);
  free(cfg.msafrq);
  free(cfg.method);

  return 0;
  
 ERROR:
  return 0;
}

static int
write_output (struct cfg_s *cfg)
{
  double time;
  
  time = cfg->w->user + cfg->w->sys;

  printf("#training set: %s\n", cfg->outheader);
  printf("rI %f rM %f rD %f | sI %f sD %f vI %f vD %f | ldEM %f muEM %f ldED %f muED %f ldI %f muI %f muAM %f muAD %f muAI %f | model %s\n",  
	 cfg->R->rI, cfg->R->rM, cfg->R->rD, 
	 cfg->R->sI, cfg->R->sD, cfg->R->vI, cfg->R->vD, 
	 cfg->R->ldE[e1R_S], cfg->R->muE[e1R_S], cfg->R->ldE[e1R_D], cfg->R->muE[e1R_D], 
	 cfg->R->ldE[e1R_I], cfg->R->muE[e1R_I], 
	 cfg->R->muA[e1R_S], cfg->R->muA[e1R_D], cfg->R->muA[e1R_I], 
	 e1_rate_EvomodelType(cfg->R->evomodel));
  printf("#CPU time (s): %.2f\n", time);

  fprintf(cfg->outfp, "#training set: %s\n", cfg->outheader);
  fprintf(cfg->outfp, "rI %f rM %f rD %f | sI %f sD %f vI %f vD %f | ldEM %f muEM %f ldED %f muED %f ldI %f muI %f muAM %f muAD %f muAI %f | model %s\n",  
	  cfg->R->rI, cfg->R->rM, cfg->R->rD, 
	  cfg->R->sI, cfg->R->sD, cfg->R->vI, cfg->R->vD, 
	  cfg->R->ldE[e1R_S], cfg->R->muE[e1R_S], cfg->R->ldE[e1R_D], cfg->R->muE[e1R_D], 
	  cfg->R->ldE[e1R_I], cfg->R->muE[e1R_I], 
	  cfg->R->muA[e1R_S], cfg->R->muA[e1R_D], cfg->R->muA[e1R_I], 
	  e1_rate_EvomodelType(cfg->R->evomodel));
  fprintf(cfg->outfp, "#CPU time (s): %.2f\n", time);

  return eslOK;
}

static int
run_voutput (struct cfg_s *cfg)
{
  char           *msaname = NULL;
  ESL_MSA        *e2msa = NULL;
  MSA_STAT       *msastat = NULL;
  MSA_STAT       *alle2msastat = NULL;
  MSA_STAT       *e2msastat = NULL;
  double          time;
  double          ttime = 0.0;
  float           sc;
  float           treeavgt;
  int             n;
  int             status;
  
  /* open outfile with benchmarking results */
  
  if (cfg->amoeba) esl_sprintf(&cfg->benchfile, "%sNM.%s.bench", cfg->outheader, e1_rate_EvomodelType(cfg->evomodel));  
  else             esl_sprintf(&cfg->benchfile, "%sGD.%s.bench", cfg->outheader, e1_rate_EvomodelType(cfg->evomodel));   

  if ((cfg->benchfp = fopen(cfg->benchfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg->benchfile);
  fprintf(cfg->benchfp, "#MSA_NAME           TP_hom   TRUE_hom  FOUND_hom      \n");
  
  /* align the training set with the optimal parameters */
  for (n = 0; n < cfg->nmsa; n ++) { 
    esl_sprintf(&msaname, "%s.%d", cfg->msaheader, n); 

    esl_stopwatch_Start(cfg->w);
    status = e2_msa(cfg->r, cfg->R, NULL, cfg->msalist[n], cfg->msafrq[n], cfg->Tlist[n], &e2msa, &sc, cfg->pli, cfg->bg, NULL, E2, OPTNONE, 
		    cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
    esl_stopwatch_Stop(cfg->w);
    
    time   = cfg->w->user + cfg->w->sys;
    ttime += time;
    
    treeavgt = esl_tree_er_AverageBL(cfg->Tlist[n]);

    msamanip_XStats(cfg->msalist[n], &msastat);                  //  reference msa aveid and avematch 
    msamanip_CStats(cfg->msalist[n]->abc, e2msa, &alle2msastat); //  optimal   msa aveid and avematch (with ancestral sequences) 
    
    /* The leaves-only alignment */
    msamanip_MSALeaves(&e2msa, FALSE);   
    msamanip_CStats(cfg->msalist[n]->abc, e2msa, &e2msastat);
    e2msastat->anclen = alle2msastat->anclen; // transfer the len of the ancestral sequence (cannot be calculated from the leaves-only msa)
    
    //if (esl_msafile_Write(stdout, e2msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write msa to file"); 
    
    status = FastSP_Benchmark(cfg->benchfp, msaname, cfg->method, cfg->abc, cfg->msalist[n], msastat, e2msa, e2msastat, sc, treeavgt, time, FALSE, FALSE, cfg->errbuf, cfg->verbose);
    if (status != eslOK) goto ERROR;
    esl_msa_Destroy(e2msa); e2msa = NULL;
  }

  if (e2msa) esl_msa_Destroy(e2msa);
  if (msastat) free(msastat);
  if (e2msastat) free(e2msastat);
  if (alle2msastat) free(alle2msastat);
  fclose(cfg->benchfp);
  free(cfg->benchfile);
   return eslOK;
  
 ERROR:
  if (e2msa) esl_msa_Destroy(e2msa);
  if (msastat) free(msastat);
  if (e2msastat) free(e2msastat);
  if (alle2msastat) free(alle2msastat);
  if (cfg->benchfp) fclose(cfg->benchfp);
  if (cfg->benchfile) free(cfg->benchfile);
   return status;
}

static int
explore_distribution(struct cfg_s *cfg)
{
  double *param = NULL;
  int     which;
  int     x = 0;
  int     np = cfg->R->nrate + cfg->R->nbern;
  int     p;
  int     status; 
  
  /* select a msa */
  which = (int)(esl_random(cfg->r) * cfg->nmsa);
  	
  /* P(msa|param,t) as a function of one of the parameters */ 
  ESL_ALLOC(param, sizeof(double) * (np));
 
  status = e2_transitions_pack_paramvector(&x, param, np, cfg->R, cfg->errbuf, cfg->verbose);
  if (status != eslOK) goto ERROR;
  
  for (p = 0; p < np; p ++) {
    if (p < cfg->R->nrate) status = move_rate     (which, p, np, param, cfg);
    else                   status = move_bernoulli(which, p, np, param, cfg);
   if (status != eslOK) goto ERROR;
  }
  
  free(param);
  return eslOK;

 ERROR:
  if (param) free(param);
  return status;
}

static int
move_rate(int which, int idx, int np, double *param, struct cfg_s *cfg)
{
  float   tinit  = -1.0;
  double  ratez = param[idx];
  double  rate;
  double  ratemax = 0.50;
  double  ratemin = 0.01;
  double  rateinc = 0.01;
  float   sc;
  float   sign = 1.0;
  int     x;
  int     status;
  
  rate = ratez;
  while (rate <= ratemax && rate <= ratemin) {
    x = 0;
    status = e2_transitions_unpack_paramvector(&x, param, np, cfg->R, cfg->errbuf, cfg->verbose);
    if (status != eslOK) goto ERROR;
    
    tinit = esl_tree_er_AverageBL(cfg->Tlist[which]);
    if (cfg->Tlist[which]) { esl_tree_Destroy(cfg->Tlist[which]); cfg->Tlist[which] = NULL; }
    status = e2_tree_UPGMA(&cfg->Tlist[which], cfg->msalist[which], cfg->msafrq[which], cfg->r, cfg->pli, cfg->R, NULL, cfg->bg, NULL, E2F, 
			   cfg->mode, cfg->do_viterbi, -1.0, tinit, cfg->tol, cfg->errbuf, cfg->verbose);     
    if (status != eslOK) goto ERROR;
    
    status = e2_msa(cfg->r, cfg->R, NULL, cfg->msalist[which], cfg->msafrq[which], cfg->Tlist[which], NULL, &sc, cfg->pli, cfg->bg, NULL, E2F, OPTNONE, 
		    cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
    if (status != eslOK) goto ERROR;
    
    printf("%f %f\n", rate, sc);
    param[idx] += sign * RATE2PARAM(rateinc); 
    rate = PARAM2RATE(param[idx]);
    if (rate > ratemax) { param[idx] = RATE2PARAM(ratez-rateinc); sign = -1.0; }
  }
  printf("&\n");
  
  return eslOK;

 ERROR:
  return status;
}

static int
move_bernoulli(int which, int idx, int np, double *param, struct cfg_s *cfg)
{
  float   tinit  = -1.0;
  double  bzero = param[idx];
  double  bernoulli;
  double  binc = 0.01;
  double  sign = 1.0;
  float   sc;
  int     x;
  int     status;

  bernoulli = bzero;
  while (bernoulli < 1.0 && bernoulli > 0.0) {
    x = 0;
    status = e2_transitions_unpack_paramvector(&x, param, np, cfg->R, cfg->errbuf, cfg->verbose);
    if (status != eslOK) goto ERROR;
    
    tinit = esl_tree_er_AverageBL(cfg->Tlist[which]);
    if (cfg->Tlist[which]) { esl_tree_Destroy(cfg->Tlist[which]); cfg->Tlist[which] = NULL; }
    status = e2_tree_UPGMA(&cfg->Tlist[which], cfg->msalist[which], cfg->msafrq[which], cfg->r, cfg->pli, cfg->R, NULL, cfg->bg, NULL, E2F, 
			   cfg->mode, cfg->do_viterbi, -1.0, tinit, cfg->tol, cfg->errbuf, cfg->verbose);     
    if (status != eslOK) goto ERROR;
    
    status = e2_msa(cfg->r, cfg->R, NULL, cfg->msalist[which], cfg->msafrq[which], cfg->Tlist[which], NULL, &sc, cfg->pli, cfg->bg, NULL, E2F, OPTNONE, 
		    cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
    if (status != eslOK) goto ERROR;
    
    printf("%f %f\n", bernoulli, sc);
    param[idx] = BERN2PARAMLOG(bernoulli + sign*binc); 
    bernoulli  = PARAM2BERNLOG(param[idx]);
    if (bernoulli >= 1.0) { param[idx] = RATE2PARAM(bzero-binc); sign = -1.0; }
  }
  printf("&\n");
  
  return eslOK;

 ERROR:
  return status;
}
