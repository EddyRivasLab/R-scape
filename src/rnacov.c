/* rnacov -- statistical test for RNA covatiation in an alignment
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



#define ALPHOPTS "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
/* Exclusive options for evolutionary model choice */

/* struct cfg_s : "Global" application configuration shared by all threads/processes.
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  /* Shared configuration in masters & workers */
  int              argc;
  char           **argv;
  
  ESL_STOPWATCH   *w;

  char             errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS  *r;	               /* random numbers for stochastic sampling (grm-emit) */
  ESL_ALPHABET    *abc;                /* the alphabet */
  double           fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double           idthresh;	       /* fractional identity threshold for selecting subset of sequences */

  int              nmsa;
  char            *msafile;
  
  float              tol;
  int                verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "-I",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < x id",                                                               0 },
  { "--submsa",       eslARG_INT,       FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   0 },
  { "--minid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      0 },
  { "--maxid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",                                                      0 },
  /* options for evolutionary method */
  { "--fixpid",       eslARG_REAL,       NULL,   NULL,"0<=x<=100",   NULL,    NULL,  NULL,               "TRUE for using a fix pid for the ehmm model",                                               0 },
  { "--fixtime",      eslARG_REAL,       NULL,   NULL,     "x>=0",      NULL, NULL,  NULL,               "TRUE for using a fix time for the evolutionary models of a pair",                           0 },
  /* Control of scoring system - substitutions */
  { "--mx",           eslARG_STRING,"BLOSUM62",  NULL,       NULL,   NULL,    NULL,  "--mxfile",         "substitution rate matrix choice (of some built-in matrices)",                               0 },
  { "--mxfile",       eslARG_INFILE,      NULL,  NULL,       NULL,   NULL,    NULL,  "--mx",             "read substitution rate matrix from file <f>",                                               0 },
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
static char usage[]  = "[-options] <msa>";
static char banner[] = "rnacov - statistical test for RNA covatiation in an alignment";

static int create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);

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
      (cfg.infmt = eslx_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nmsa = 0;

  /* alphabet */
  if      (esl_opt_GetBoolean(go, "--amino"))   cfg.abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg.abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg.abc = esl_alphabet_Create(eslRNA);
  else                                          cfg.abc = NULL;

  cfg.w = esl_stopwatch_Create(); 

  /*  output file */
  if ( esl_opt_IsOn(go, "-o") ) {
    if ((cfg.outfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else cfg.outfp = stdout;
  
  if (esl_opt_IsOn(go, "--submsa")) cfg.submsa = esl_opt_GetInteger(go, "--submsa");
  else                              cfg.submsa = 0;
  
  /* other options */
  cfg.fragfrac   = esl_opt_IsOn(go, "-F")? esl_opt_GetReal(go, "-F") : -1.0;
  cfg.idthresh   = esl_opt_IsOn(go, "-I")? esl_opt_GetReal(go, "-I") : -1.0;
  cfg.tol        = esl_opt_GetReal   (go, "--tol");
  cfg.verbose    = esl_opt_GetBoolean(go, "-v");
  cfg.voutput    = esl_opt_GetBoolean(go, "--voutput");

  /* branch length options */
  cfg.fixpid     = esl_opt_IsOn(go, "--fixpid")?  esl_opt_GetReal(go, "--fixpid") : -1.0; 
  cfg.fixtime    = esl_opt_IsOn(go, "--fixtime")? esl_opt_GetReal(go, "--fixtime") : -1.0;

  cfg.T = NULL;

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
  char           *msg = "e2msa failed";
  ESL_GETOPTS    *go;
  struct cfg_s    cfg;
  ESLX_MSAFILE   *afp = NULL;
  ESL_MSA        *msa = NULL;          /* the input alignment  */
  int             seq_cons_len = 0;
  int             nfrags = 0;	  	  /* # of fragments removed */
  int             nremoved = 0;	          /* # of identical sequences removed */
  int             status = eslOK;
  int             hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* Open the MSA file */
  status = eslx_msafile_Open(NULL, cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);

  /* read the MSA */
  while ((hstatus = eslx_msafile_Read(afp, &msa)) != eslEOF) {
    if (hstatus != eslOK) eslx_msafile_ReadFailure(afp, status);
    cfg.nmsa ++;

    /* select submsa */
    if (cfg.submsa) {
      if (MSA_Subset(cfg.r, cfg.submsa, &msa, NULL, cfg.errbuf, cfg.verbose) != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
    }
    
    /* outheader for all msa-output files */
    msamanip_OutfileHeader((msa->acc)?msa->acc:cfg.msafile, &cfg.msaheader); 
   
    esl_msa_Hash(msa);
   
    if (esl_opt_IsOn(go, "-F") && msamanip_RemoveFragments(cfg.fragfrac, &msa, &nfrags, &seq_cons_len) != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-I"))   msamanip_SelectSubset(cfg.r, &msa, cfg.idthresh, &nremoved);
    
    /* given msa aveid and avematch */
    msamanip_CStats(cfg.abc, msa, &cfg.mstat);
    msamanip_CBaseComp(cfg.abc, msa, cfg.bg->f, &cfg.msafrq);

    if (esl_opt_IsOn(go, "--minid") && cfg.mstat.avgid < 100.*esl_opt_GetReal(go, "--minid")) continue;
    if (esl_opt_IsOn(go, "--maxid") && cfg.mstat.avgid > 100.*esl_opt_GetReal(go, "--maxid")) continue;

    /* print some info */
    if (cfg.voutput) {
      esl_sprintf(&cfg.msarfname, "%s.%s", cfg.msaheader, e1_rate_EvomodelType(cfg.evomodel)); 
      fprintf(cfg.outfp, "Given alignment\n");
      fprintf(cfg.outfp, "%6d          %s\n", msa->nseq, cfg.msafile);
      if (eslx_msafile_Write(cfg.outfp, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write to file %s", cfg.msarfname); 
      msamanip_DumpStats(cfg.outfp, msa, cfg.mstat); 
    }
    
    status = run_rnacov(go, &cfg, msa);
    if (status != eslOK) esl_fatal("%s Failed to run e2msa", cfg.errbuf);
    
    esl_msa_Destroy(msa); msa = NULL;
    free(cfg.msafrq); cfg.msafrq = NULL;
  }

  /* cleanup */
  esl_stopwatch_Destroy(cfg.w);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  eslx_msafile_Close(afp);
  return 0;
}


static int
create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  int      status;
  
  /* the TREE */
  status = e2_tree_UPGMA(&cfg->T, msa, cfg->msafrq, cfg->r, cfg->pli, cfg->R1[0], cfg->R7, cfg->bg, cfg->bg7, cfg->e2ali, cfg->mode, cfg->do_viterbi, cfg->fixtime, -1.0,
			 cfg->tol, cfg->errbuf, cfg->verbose);
  if (cfg->verbose) Tree_Dump(cfg->outfp, cfg->T, "Tree");
  
  cfg->treeavgt = esl_tree_er_AverageBL(cfg->T);
  
  return eslOK;
}


static int
run_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  float      sc;
  int        nnodes;
  int        r;
  int        status;

  esl_stopwatch_Start(cfg->w);

  dmsa = esl_msa_Clone(msa);
  esl_msa_Digitize(cfg->abc, dmsa, NULL);
  
  /* get the tree
   */
  status = create_tree(go, cfg, dmsa);
  if (status != eslOK)  { esl_fatal(cfg->errbuf); }
  nnodes = (cfg->T->N > 1)? cfg->T->N-1 : cfg->T->N;

 /* main function */
  status = e2_msa(cfg->r, cfg->R1[0], cfg->R7, dmsa, cfg->msafrq, cfg->T, &e2msa, &sc, cfg->pli, cfg->bg, cfg->bg7, cfg->e2ali, cfg->e2optimize, 
		  cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
  if (status != eslOK)  { esl_fatal(cfg->errbuf); }
  esl_stopwatch_Stop(cfg->w);
  msamanip_CStats(cfg->abc, e2msa, &alle2mstat);
  
#if 0
  if (cfg->voutput) { /* write output alignment */
    fprintf(cfg->outfp, "\ne2msa-%s alignment with nodes\n", e1_rate_EvomodelType(cfg->evomodel));
    eslx_msafile_Write(cfg->outfp, e2msa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, e2msa, alle2mstat);
  }
#endif
  
  if (cfg->T) esl_tree_Destroy(cfg->T);
  
  return eslOK;

 ERROR:
 return status;
}

