/* rnacov -- statistical test for covariation in an RNA structural alignment
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
 #include "esl_wuss.h"
 
#include "msamanip.h"
#include "msatree.h"
#include "mutualinfo.h"
#include "ribosum_matrix.h"



#define ALPHOPTS "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
#define METHODOPTS "--naive,--phylo,--dca,--akmaev"              

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
  double           gapthresh;          /* only keep columns with <= gapthresh fraction of gaps in them */

  int              nmsa;
  char            *msafile;
  char            *filename;

  FILE            *outmsafp;
 
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */
  char            *msaheader;          /* header for all msa-specific output files */
  int              infmt;
  
  int              domsa;
  int              doshuffle;

  METHOD           method;
  ESL_TREE        *T;
  double           treeavgt;
 
  char                *ribofile;
  struct ribomatrix_s *ribosum;

  char            *gnuplot;
  
  int              submsa;              /* set to the number of random seqs taken from original msa.
					 * Set to 0 if we are taking all */
  MSA_STAT         omstat;              /* statistics of the original alignment */
  MSA_STAT         mstat;               /* statistics of the analyzed alignment */
  float           *msafrq;
  int             *ct;
  int              onbpairs;
  int              nbpairs;

  int              voutput;
  char            *rocfile;
  FILE            *rocfp; 
  char            *sumfile;
  FILE            *sumfp; 
  char            *shsumfile;
  FILE            *shsumfp; 
  int              maxFP;
  double           ratioFP;

  float            tol;
  int              verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "--maxFP",        eslARG_INT,       FALSE,   NULL,     "n>=0",   NULL,    NULL,  NULL,               "maximum number of FP allowed",                                                              0 },
  { "--ratioFP",      eslARG_REAL,      "0.2",   NULL,   "x>=0.0",   NULL,    NULL,  NULL,               "maximum ration of FP allowed",                                                              0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  /* method */
  { "--naive",        eslARG_NONE,       TRUE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "naive MI calculations",                                                                     0 },
  { "--phylo",        eslARG_NONE,      FALSE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "phylo MI calculations",                                                                     0 },
  { "--dca",          eslARG_NONE,      FALSE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "direct coupling analysis (DCA) MI calculations",                                            0 },
  { "--akmaev",       eslARG_NONE,      FALSE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "akmaev-style MI calculations",                                                              0 },
  /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "-I",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < <x> id",                                                             0 },
  { "--submsa",       eslARG_INT,       FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   0 },  
  { "--gapthresh",    eslARG_REAL,      FALSE,   NULL,  "0<=x<=1",   NULL,    NULL,  NULL,               "keep columns with < <x> fraction of gaps",                                                  0 },
  { "--minid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      0 },
  { "--maxid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",                                                      0 },
  /* Control of scoring system - ribosum */
  { "--ribofile",     eslARG_INFILE,    NULL,    NULL,       NULL,   NULL,    NULL,  "--mx",             "read ribosum structure from file <f>",                                                      0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--outmsa",       eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "write actual msa used to file <f>,",                                                        0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
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
static int run_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, int ishuffled);

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
      esl_banner(stdout, cfg.argv[0], banner);
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
  
  cfg.msaheader = NULL;

   /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = eslx_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nmsa = 0;
  cfg.msafrq = NULL;

  /* alphabet */
  cfg.abc = esl_alphabet_Create(eslRNA);
  esl_alphabet_SetEquiv(cfg.abc, '=', '-');     /* allow = as a gap character too */
  esl_alphabet_SetEquiv(cfg.abc, '.', '-');     /* allow . as a gap character too */

  cfg.w = esl_stopwatch_Create(); 
  
  /*  output file */
  if ( esl_opt_IsOn(go, "-o") ) {
    if ((cfg.outfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else cfg.outfp = stdout;

   /*  outmsa file */
  cfg.outmsafp = NULL;
  if ( esl_opt_IsOn(go, "--outmsa") ) {
    if ((cfg.outmsafp = fopen(esl_opt_GetString(go, "--outmsa"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "--outmsa"));
  } 
  
  esl_FileTail(cfg.msafile, TRUE, &cfg.outheader);  
  if (esl_opt_IsOn(go, "--submsa")) { cfg.submsa = esl_opt_GetInteger(go, "--submsa"); esl_sprintf(&cfg.outheader, "%s_random%d", cfg.outheader, cfg.submsa); }
  else                              { cfg.submsa = 0; }
    
  /* other options */
  cfg.domsa      = TRUE;
  cfg.doshuffle  = TRUE;
  cfg.fragfrac   = esl_opt_IsOn(go, "-F")?          esl_opt_GetReal   (go, "-F")          : -1.0;
  cfg.idthresh   = esl_opt_IsOn(go, "-I")?          esl_opt_GetReal   (go, "-I")          : -1.0;
  cfg.gapthresh  = esl_opt_IsOn(go, "--gapthresh")? esl_opt_GetReal   (go, "--gapthresh") : -1.0;
  cfg.maxFP      = esl_opt_IsOn(go, "--maxFP")?     esl_opt_GetInteger(go, "--maxFP")     : -1;
  cfg.ratioFP    = esl_opt_GetReal   (go, "--ratioFP");
  cfg.tol        = esl_opt_GetReal   (go, "--tol");
  cfg.verbose    = esl_opt_GetBoolean(go, "-v");
  cfg.voutput    = esl_opt_GetBoolean(go, "--voutput");

  if      (esl_opt_GetBoolean(go, "--naive"))  cfg.method = OPTNONE;
  else if (esl_opt_GetBoolean(go, "--phylo"))  cfg.method = PHYLO;
  else if (esl_opt_GetBoolean(go, "--dca"))    cfg.method = DCA;
  else if (esl_opt_GetBoolean(go, "--akmaev")) cfg.method = AKMAEV;
 
 /*  rocplot file */
  esl_sprintf(&cfg.rocfile, "%s.g%.1f.r%.1f.roc", cfg.outheader, cfg.gapthresh, cfg.ratioFP); 
  if ((cfg.rocfp = fopen(cfg.rocfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.rocfile);
  printf("rocfile %s\n", cfg.rocfile);

  /*  summary file */
  esl_sprintf(&cfg.sumfile, "%s.g%.1f.r%.1f.sum", cfg.outheader, cfg.gapthresh, cfg.ratioFP); 
  if ((cfg.sumfp = fopen(cfg.sumfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.sumfile);
  printf("sumfile %s\n", cfg.sumfile);
  
  cfg.shsumfile = NULL;
  cfg.shsumfp = NULL;
  if (cfg.doshuffle) {
    /*  sh-summary file */
    esl_sprintf(&cfg.shsumfile, "%s.g%.1f.r%.1f.shsum", cfg.outheader, cfg.gapthresh, cfg.ratioFP); 
    if ((cfg.shsumfp = fopen(cfg.shsumfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.shsumfile);
    printf("sh-sumfile %s\n", cfg.shsumfile);
  }
  
  cfg.T  = NULL;
  cfg.ct = NULL;
  cfg.nbpairs = 0;

  /* the ribosum matrices */
  cfg.ribofile = NULL;
  cfg.ribosum  = NULL;
  if (cfg.method == PHYLO || cfg.method == AKMAEV) {
    if ( esl_opt_IsOn(go, "--ribofile") ) { cfg.ribofile = esl_opt_GetString(go, "--ribofile"); }
    else esl_sprintf(&cfg.ribofile, "ssu-lsu.ribosum");
    
    cfg.ribosum = Ribosum_matrix_Read(cfg.ribofile, cfg.abc, FALSE, cfg.errbuf);
    if (cfg.ribosum == NULL) esl_fatal("%s\nfailed to create ribosum matrices from file %s\n", cfg.errbuf, cfg.ribofile);
    if (cfg.verbose) Ribosum_matrix_Write(stdout, cfg.ribosum);
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
  char            *msg = "e2msa failed";
  ESL_GETOPTS     *go;
  struct cfg_s     cfg;
  ESLX_MSAFILE    *afp = NULL;
  ESL_MSA         *msa = NULL;            /* the input alignment  */
  ESL_MSA         *shmsa = NULL;          /* the shuffled alignment  */
  int              seq_cons_len = 0;
  int              nfrags = 0;	  	  /* # of fragments removed */
  int              nremoved = 0;	  /* # of identical sequences removed */
  int              status = eslOK;
  int              hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* Open the MSA file */
  status = eslx_msafile_Open(NULL, cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);
  eslx_msafile_SetDigital(afp, cfg.abc);

  /* read the MSA */
  while ((hstatus = eslx_msafile_Read(afp, &msa)) != eslEOF) {
    if (hstatus != eslOK) eslx_msafile_ReadFailure(afp, status);
    cfg.nmsa ++;

    /* stats of the orignial alignment */
    msamanip_XStats(msa, &cfg.omstat);
    msamanip_CalculateCT(msa, NULL, &cfg.onbpairs, cfg.errbuf);

   /* select submsa and then apply msa filters 
    */
    if (cfg.submsa                      && msamanip_SelectSubset(cfg.r, cfg.submsa, &msa, &cfg.outheader, cfg.errbuf, cfg.verbose) != eslOK) { printf("%s\n", cfg.errbuf);          esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-F")          && msamanip_RemoveFragments(cfg.fragfrac, &msa, &nfrags, &seq_cons_len)                    != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-I")          && msamanip_SelectSubsetByID(cfg.r, &msa, cfg.idthresh, &nremoved)                         != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "--gapthresh") && msamanip_RemoveGapColumns(cfg.gapthresh, msa, cfg.errbuf, cfg.verbose)                  != eslOK) { printf("RemoveGapColumns\n");        esl_fatal(msg); }
 
    esl_msa_Hash(msa);
    esl_msa_ConvertDegen2X(msa);
    if (esl_msa_MinimGaps(msa, NULL, "-.~=", FALSE) != eslOK) esl_fatal("Failed to remove minim gaps");
    
    /* given msa aveid and avematch */
    msamanip_XStats(msa, &cfg.mstat);
    
    if (esl_opt_IsOn(go, "--minid") && cfg.mstat.avgid < 100.*esl_opt_GetReal(go, "--minid")) continue;
    if (esl_opt_IsOn(go, "--maxid") && cfg.mstat.avgid > 100.*esl_opt_GetReal(go, "--maxid")) continue;
    
    /* output the actual file used if requested */
    if ( esl_opt_IsOn(go, "--outmsa") ) eslx_msafile_Write(cfg.outmsafp, msa, eslMSAFILE_STOCKHOLM);;
 
    /* print some info */
    if (cfg.voutput) {
      fprintf(cfg.outfp, "Used alignment\n");
      fprintf(cfg.outfp, "%6d          %s\n", msa->nseq, cfg.msafile);
      if (eslx_msafile_Write(cfg.outfp, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write msa"); 
      msamanip_DumpStats(cfg.outfp, msa, cfg.mstat); 
    }
    
    /* the ct vector */
    status = msamanip_CalculateCT(msa, &cfg.ct, &cfg.nbpairs, cfg.errbuf);
    if (status != eslOK) esl_fatal("%s. Failed to calculate ct vector", cfg.errbuf);
    
    /* main function */
    if (cfg.domsa) {
      status = run_rnacov(go, &cfg, msa, FALSE);
      if (status != eslOK) esl_fatal("Failed to run rnacov");
    }

    if (cfg.doshuffle) {
      msamanip_ShuffleColums(cfg.r, msa, &shmsa, cfg.errbuf, cfg.verbose);
      status = run_rnacov(go, &cfg, shmsa, TRUE);
      if (status != eslOK) esl_fatal("%s. Failed to run rnacov shuffled", cfg.errbuf);
    }

    esl_msa_Destroy(msa); msa = NULL;
    free(cfg.ct); cfg.ct = NULL;
    if (shmsa) esl_msa_Destroy(shmsa); shmsa = NULL;
    if (cfg.msafrq) free(cfg.msafrq); cfg.msafrq = NULL;
    if (cfg.T) esl_tree_Destroy(cfg.T); cfg.T = NULL;
    if (cfg.msaheader) free(cfg.msaheader); cfg.msaheader = NULL;
  }

  /* cleanup */
  esl_stopwatch_Destroy(cfg.w);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  if (cfg.ct) free(cfg.ct);
  eslx_msafile_Close(afp);
  fclose(cfg.outfp);
  fclose(cfg.rocfp);
  fclose(cfg.sumfp);
  if (cfg.shsumfp) fclose(cfg.shsumfp);
  free(cfg.outheader);
  free(cfg.gnuplot);
  if (cfg.ribosum) Ribosum_matrix_Destroy(cfg.ribosum);
  return 0;
}


static int
create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  int      status;
  
  /* the TREE */
  status = Tree_CalculateExtFromMSA(msa, &cfg->T, TRUE, cfg->errbuf, cfg->verbose);
  if (status != eslOK) { printf("%s\n", cfg->errbuf); esl_fatal(cfg->errbuf); }

  if (cfg->verbose) Tree_Dump(cfg->outfp, cfg->T, "Tree");
  
  cfg->treeavgt = esl_tree_er_AverageBL(cfg->T);
  
  return eslOK;
}


static int
run_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, int ishuffled)
{
  struct mutual_s *mi = NULL;
  char            *name = NULL;
  int              nnodes;
  int              status;

  esl_stopwatch_Start(cfg->w);

  /* produce a tree
   */
  if (cfg->method != NAIVE) {
    status = create_tree(go, cfg, msa);
    if (status != eslOK)  { esl_fatal(cfg->errbuf); }
    nnodes = (cfg->T->N > 1)? cfg->T->N-1 : cfg->T->N;
  }

  /* create the MI structure */
  mi = Mutual_Create(msa->alen, msa->nseq, cfg->abc);

  /* write MSA info to the sumfile */
  if (msa->acc && msa->name) esl_sprintf(&name, "%s_%s", msa->acc, msa->name);
  else if (msa->acc)         esl_sprintf(&name, "%s", msa->acc);
  else                       esl_sprintf(&name, "%s", cfg->outheader);

  if (!ishuffled) 
    fprintf(cfg->sumfp, "%f\t%s\t%d\t%.2f\t", cfg->ratioFP, name, msa->nseq, cfg->mstat.avgid); 
  else
    fprintf(cfg->shsumfp, "%f\t%s\t%d\t%.2f\t", cfg->ratioFP, name, msa->nseq, cfg->mstat.avgid); 
  
  /* write MSA info to the rocfile */
  fprintf(cfg->rocfp, "# MSA nseq %d alen %" PRId64 " avgid %f\n", msa->nseq, msa->alen, cfg->mstat.avgid);  
 
  /* main function */
  status = Mutual_Calculate(msa, cfg->T, cfg->ribosum, mi, cfg->method, cfg->ct, cfg->rocfp, (!ishuffled)?cfg->sumfp:cfg->shsumfp, cfg->maxFP,
			    cfg->ratioFP, cfg->onbpairs, ishuffled, cfg->tol, cfg->verbose, cfg->errbuf);   
  if (status != eslOK)  { goto ERROR; }

  /* print to stdout */
  fprintf(stdout,     "# MSA %s nseq %d (%d) alen %" PRId64 " (%" PRId64 ") avgid %.2f (%.2f) nbpairs %d (%d)\n", 
	  name, msa->nseq, cfg->omstat.nseq, msa->alen, cfg->omstat.alen, 
	  cfg->mstat.avgid, cfg->omstat.avgid, cfg->nbpairs, cfg->onbpairs);  

  Mutual_Destroy(mi); mi = NULL;
  free(name);
  return eslOK;

 ERROR:
  if (cfg->T) esl_tree_Destroy(cfg->T);
  if (mi) Mutual_Destroy(mi);
  if (name) free(name);
 return status;
}

