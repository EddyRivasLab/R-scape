/* e2sim -- e2msa alignment by simulation
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esl_getopts.h"
#include "esl_distance.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_stopwatch.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2.h"
#include "e2_evolve.h"
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
 
  int                  ntaxa;
  double               fixabl;
  double               minabl;
  double               maxabl;
  double               incabl;
  ESL_TREE            *T;

  int                  nmsa;
  int                  L;                      /* expected length of the ancestral sequence */
  char                *msaname;
  ESL_MSA             *msa;
  float               *msafrq;
  MSA_STAT            *msastat;

  char                *outfile;
  FILE                *outfp;
  int                  voutput;                 /* TRUE for verbose output */

  char                *paramfile;
  EVOM                 evomodel;
  struct rateparam_s   rateparam;
  
  char                *benchfile;
  FILE                *benchfp;
  
  char                *gnuplot;
  char                *outheader;              /* header for all output files */
  
  E2_PIPELINE         *pli;
  E1_BG               *bg;	   	       /* null model (copies made of this into threads) */
  int                  mode;                   /* e2_GLOBAL, e2_LOCAL */ 
  int                  do_viterbi;             // default is optimal accuracy

  E1_RATE             *R;
  char                *subsmx;                 /* BLUSUM62, ... */
  EMRATE              *emR;
  
  char                *method;
  float                tol;
  int                  verbose;
};

static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
 /* parameters to control the simulation */
  { "--nmsa",         eslARG_INT,      "100",    NULL,      "n>0",   NULL,    NULL,  NULL,               "set nmsa to <n>",                                                                           0 },
  { "-L",             eslARG_INT,     "1000",    NULL,      "n>0",   NULL,    NULL,  NULL,               "expected length of ancestral sequences",                                                    0 },
  /* tree options */
  { "--ntaxa",        eslARG_INT,       "2",    NULL,      "n> 0",   NULL,    NULL,  NULL,               "set ntaxa to <n>",                                                                          0 },
  { "--fixabl",       eslARG_REAL,     NULL,    NULL,       "x>0",   NULL,    NULL,  NULL,               "set for a fixed tree average branch lenght",                                                0 },
  { "--minabl",       eslARG_REAL,   "0.01",    NULL,       "x>0",   NULL,    NULL,  NULL,               "minimum tree average branch lenght",                                                        0 },
  { "--maxabl",       eslARG_REAL,   "2.00",    NULL,       "x>0",   NULL,    NULL,  NULL,               "maximum tree average branch lenght",                                                        0 },
  { "--incabl",       eslARG_REAL,   "0.01",    NULL,       "x>0",   NULL,    NULL,  NULL,               "increment tree average branch lenght",                                                      0 },
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
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
 /* Selecting the alphabet rather than autoguessing it */
  { "--amino",        eslARG_NONE,      TRUE,    NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is protein sequence data",                                                  2 },
  { "--dna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is DNA sequence data",                                                      2 },
  { "--rna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is RNA sequence data",                                                      2 },
   /* other options */
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,       "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "e2msa [-options] ";
static char banner[] = "simulate e2msa";

static int generate_random_tree(struct cfg_s  *cfg, double abl, ESL_TREE **ret_T);
static int generate_msa(struct cfg_s *cfg, ESL_MSA **ret_msa, int incnode);
static int run_voutput (struct cfg_s *cfg, int n);

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

  cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));

  /* alphabet */
  if      (esl_opt_GetBoolean(go, "--amino"))   cfg.abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg.abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg.abc = esl_alphabet_Create(eslRNA);
  else                                          cfg.abc = NULL;

   /* tree options */
  cfg.fixabl = (esl_opt_IsOn(go, "--fixabl"))? esl_opt_GetReal(go, "--fixabl") : -1.0;
  cfg.minabl = esl_opt_GetReal   (go, "--minabl");
  cfg.maxabl = esl_opt_GetReal   (go, "--maxabl");
  cfg.incabl = esl_opt_GetReal   (go, "--incabl");
  cfg.ntaxa  = esl_opt_GetInteger(go, "--ntaxa");
  cfg.T      = NULL;
  
  /* msa options */
  cfg.L      = esl_opt_GetInteger(go, "-L");
  cfg.nmsa   = esl_opt_GetInteger(go, "--nmsa");
  cfg.msa    = NULL;
  cfg.msafrq = NULL;

   /* other options */
  cfg.tol     = esl_opt_GetReal   (go, "--tol");
  cfg.verbose = esl_opt_GetBoolean(go, "-v");
  cfg.voutput = esl_opt_GetBoolean(go, "--voutput");
  
  cfg.subsmx    = esl_opt_GetString(go, "--mx");
  cfg.mode      = e2_GLOBAL;
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
  esl_sprintf(&cfg.method, "e2sim.%s", e1_rate_EvomodelType(cfg.evomodel)); 

  cfg.paramfile = NULL;
  if (esl_opt_IsOn(go, "--paramfile")) {
    cfg.paramfile = esl_opt_GetString(go, "--paramfile");
    status = e1_rate_ReadParamfile(cfg.paramfile, &cfg.rateparam, &cfg.evomodel, cfg.errbuf, cfg.verbose);
    if (status != eslOK) esl_fatal("Failed to read paramfile %s\n%s", cfg.paramfile, cfg.errbuf);
  }

  /* open outfile */
  cfg.outfile = NULL;
  cfg.outfp   = NULL;
  if (esl_opt_IsOn(go, "-o"))
    cfg.outfile = esl_opt_GetString(go, "-o");
  else 
    esl_sprintf(&cfg.outfile, "e2sim.%s.sto", e1_rate_EvomodelType(cfg.evomodel)); 
  if ((cfg.outfp = fopen(cfg.outfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outfile);
   
  cfg.R   = NULL;
  cfg.bg  = e1_bg_Create(cfg.abc);
  e1_bg_SetLength(cfg.bg, cfg.L);
  cfg.emR = NULL;

  /* outheader for all output files */
  cfg.outheader = cfg.outfile;
  msamanip_OutfileHeader(cfg.outheader, &cfg.outheader); 

  cfg.benchfile = NULL;
  cfg.pli       = NULL;
  if (cfg.voutput) {
    esl_sprintf(&cfg.benchfile, "%s.%s.bench", cfg.outheader, e1_rate_EvomodelType(cfg.evomodel));
    if ((cfg.benchfp = fopen(cfg.benchfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.benchfile);
    fprintf(cfg.benchfp, "#MSA_NAME           TP_hom   TRUE_hom  FOUND_hom      \n");
    /* Initializations */
    cfg.pli = e2_pipeline_Create(go, 1, 100, 100);
  }
  
  *ret_go = go;
  *ret_cfg = cfg;
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
  char           *msg = "e2sim failed";
  ESL_GETOPTS    *go;
  struct cfg_s    cfg;
  double          abl;
  int             n;
  int             status;
 
  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  if (cfg.fixabl > 0.) abl = cfg.fixabl;
  else                 abl = cfg.minabl;

  /* Create the evolutionary rate model */
  cfg.R = e1_rate_CreateWithValues(cfg.abc, cfg.evomodel, cfg.rateparam, cfg.subsmx, NULL, TRUE, cfg.tol, cfg.errbuf, cfg.verbose);
  if (cfg.R == NULL) { printf("Bad rate model.\n"); esl_fatal(cfg.errbuf); }
  
  for (n = 0; n < cfg.nmsa; n ++) {
    
    status = generate_random_tree(&cfg, abl, &cfg.T);  if (status != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
    status = generate_msa(&cfg, &cfg.msa, FALSE);      if (status != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
    while (cfg.msa->alen == 0) { //redo if alignment of leave ends up having length zero 
      esl_msa_Destroy(cfg.msa); cfg.msa = NULL;
      status = generate_msa(&cfg, &cfg.msa, FALSE);    if (status != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
    }
    
    /* write the msa */
    if (eslx_msafile_Write(cfg.outfp, cfg.msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write to file %s", cfg.outfile); 
    
    if (cfg.voutput) run_voutput (&cfg, n);
     
    esl_tree_Destroy(cfg.T);   cfg.T   = NULL;
    esl_msa_Destroy (cfg.msa); cfg.msa = NULL;
    if (cfg.msafrq) free(cfg.msafrq); cfg.msafrq = NULL;

    if (cfg.fixabl < 0.) 
      abl = (abl+cfg.incabl <= cfg.maxabl)? abl+cfg.incabl : cfg.minabl;
  }

  /* cleanup */
  if (cfg.T)   esl_tree_Destroy(cfg.T);   
  if (cfg.msa) esl_msa_Destroy (cfg.msa); 
  if (cfg.bg)  e1_bg_Destroy(cfg.bg);          
  fclose(cfg.outfp);
  if (!esl_opt_IsOn(go, "-o")) free(cfg.outfile);
  if (cfg.voutput) {
    fclose(cfg.benchfp);
    free(cfg.benchfile);
    e2_pipeline_Destroy(cfg.pli);
  }
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  free(cfg.gnuplot);
  free(cfg.method);

  return 0;
}

static int
generate_random_tree(struct cfg_s  *cfg, double abl, ESL_TREE **ret_T)
{
  ESL_TREE *T = NULL;
  int       status;

  if (esl_tree_Simulate(cfg->r, cfg->ntaxa, &T)                                  != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to simulate the tree");
  if (esl_tree_er_RescaleAverageBL(abl, &T, cfg->tol, cfg->errbuf, cfg->verbose) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to rescale the tree");
  if (esl_tree_Validate(T, NULL)                                                 != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to validate the tree");

  if (1||cfg->verbose) esl_tree_WriteNewick(stdout, T);
  
  *ret_T = T;
  return eslOK;

 ERROR:
  if (T) esl_tree_Destroy(T);
  return status;
}


static int
generate_msa(struct cfg_s  *cfg, ESL_MSA **ret_msa, int incnode)
{
  ESL_MSA *msa     = NULL;     /* alignment of leaves sequences */
  ESL_MSA *msafull = NULL;     /* alignment of leaves and internal node sequences */
  int     *useme   = NULL;
  int      i;
  int      status;

  /* Generate the alignment */
  if (e2_evolve_Alignment(cfg->r, cfg->R, cfg->T, cfg->bg, &msafull, cfg->tol, cfg->errbuf, cfg->verbose) != eslOK)
    ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to generate the alignment");
  if (cfg->verbose) 
    eslx_msafile_Write(stdout, msafull, eslMSAFILE_STOCKHOLM); 
  
  /* The leaves-only alignment */
  useme = malloc(msafull->nseq * sizeof(int));
  for (i = 0; i < msafull->nseq; i++) {
    if (!incnode && !strncmp(msafull->sqname[i], "v", 1)) useme[i] = FALSE; 
    else                                                  useme[i] = TRUE;
  }
  if (esl_msa_SequenceSubset(msafull, useme, &msa) != eslOK)
    ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to generate leaf alignment");
  if (esl_msa_MinimGaps(msa, NULL, "-", FALSE) != eslOK) 
    ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to generate leaf alignment");
  
  *ret_msa = msa;

  free(useme);
  esl_msa_Destroy(msafull);
  return eslOK;

 ERROR:
  if (msa)     esl_msa_Destroy(msa);
  if (msafull) esl_msa_Destroy(msafull);
  if (useme) free(useme);
  return status;
}


 
static int
run_voutput (struct cfg_s *cfg, int n)
{
  ESL_MSA        *e2msa = NULL;
  MSA_STAT       *alle2msastat = NULL;
  MSA_STAT       *e2msastat = NULL;
  float           sc;
  float           treeavgt;
  int             status;
  
  /* stats of the simulated alignment */
  esl_sprintf(&cfg->msaname, "%s.%s.%d", cfg->outheader, e1_rate_EvomodelType(cfg->evomodel), n); 
  fprintf(stdout, "%6d          %s\n", cfg->msa->nseq, cfg->msaname);
  eslx_msafile_Write(stdout, cfg->msa, eslMSAFILE_STOCKHOLM); 
  msamanip_CStats(cfg->abc, cfg->msa, &cfg->msastat);
  msamanip_DumpStats(stdout, cfg->msa, cfg->msastat); 
  msamanip_CBaseComp(cfg->abc, cfg->msa, cfg->bg->f, &cfg->msafrq);
  
  /* rum e2msa */
  esl_msa_Digitize(cfg->abc, cfg->msa, cfg->errbuf);
  status = e2_msa(cfg->r, cfg->R, NULL, cfg->msa, cfg->msafrq, cfg->T, &e2msa, &sc, cfg->pli, cfg->bg, NULL, E2F, OPTNONE, 
		  cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
    
  msamanip_CStats(cfg->msa->abc, e2msa, &alle2msastat); //  optimal   msa aveid and avematch (with ancestral sequences) 
  
 /* The leaves-only alignment */
  msamanip_MSALeaves(&e2msa, FALSE);   
  msamanip_CStats(cfg->msa->abc, e2msa, &e2msastat);
  e2msastat->anclen = alle2msastat->anclen; // transfer the len of the ancestral sequence (cannot be calculated from the leaves-only msa)
  
  //if (eslx_msafile_Write(stdout, e2msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write msa to file"); 
  
  treeavgt = esl_tree_er_AverageBL(cfg->T);
  status = FastSP_Benchmark(cfg->benchfp, cfg->msaname, cfg->method, cfg->abc, cfg->msa, cfg->msastat, e2msa, e2msastat, sc, treeavgt, 0.0, FALSE, FALSE, cfg->errbuf, cfg->verbose);
  if (status != eslOK) goto ERROR;
  esl_msa_Destroy(e2msa); 
  free(e2msastat);
  free(alle2msastat);
  return eslOK;
  
 ERROR:
  if (e2msa) esl_msa_Destroy(e2msa); 
  if (e2msastat) free(e2msastat);
  if (alle2msastat) free(alle2msastat);

  return status;
}

