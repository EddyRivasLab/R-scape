/* Phcar -- Phylogenetic Covariation Above Random
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "esl_getopts.h"
#include "esl_distance.h"
#include "esl_fileparser.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msafile.h"
#include "esl_msaweight.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "e2.h"

#include "phcar_config.h"

#include "msamanip.h"
#include "msatree.h"
#include "correlators.h"
#include "covariation.h"

#define ALPHOPTS     "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
#define COVTYPEOPTS  "--GTp,--MIp"              
#define THRESHOPTS   "-E"                                          

/* Exclusive options for evolutionary model choice */

/* struct cfg_s : "Global" application configuration shared by all threads/processes.
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
* of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s { /* Shared configuration in masters & workers */
  int              argc;
  char           **argv;
  ESL_STOPWATCH   *watch;
  char             errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS  *r;	               /* random numbers */
  ESL_ALPHABET    *abc;                /* the alphabet */
  int              abcisRNA;
  
  double           fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double           idthresh;	       /* fractional identity threshold for selecting subset of sequences */
  double           gapthresh;          /* only keep columns with <= gapthresh fraction of gaps in them */
  int              singlelink;         /* if TRUE use single linkage clustering to remove highly identical sequences */
  
  double           minidthresh;	       /* fractional minimal identity threshold for selecting subset of sequences */
  
  COVTYPE          covtype;
  COVCLASS         covclass;

  METHOD           method;

  int              onemsa;
  int              nmsa;
  char            *msafile;
  char            *filename;
  char            *msaname;
  int             *msamap;
  int              firstpos;
  int              maxsq_gsc;        /* maximum number of seq to switch from GSC to PG sequence weighting */
  int              tstart;
  int              tend;

  char            *outmsafile;
  FILE            *outmsafp;
 
  char            *smsafile;
  FILE            *smsafp;
 
  char            *outnullfile;
  FILE            *outnullfp;
 
  int              R2Rall;
  char            *R2Rfile;

  char            *outdir;
  char            *outfile;
  char            *outsrtfile;
  FILE            *outfp;
  FILE            *outsrtfp; 
  char            *outheader;          /* header for all output files */
  int              infmt;
   
  char            *covhisfile;
  char            *dplotfile;

  int              nshuffle;

  char            *gnuplot;
  
  int              nseqmin;
  int              consensus;           /* if TRUE, analyze only consensus positions in the alignment */
  int              submsa;              /* set to the number of random seqs taken from original msa.
					 * Set to 0 if we are taking all */
  MSA_STAT        *omstat;              /* statistics of the original alignment */
  MSA_STAT        *mstat;               /* statistics of the analyzed alignment */
  float           *msafrq;
  int             *ct;
  int              onbpairs;
  int              nbpairs;

  int              voutput;
  int              doroc;
  char            *rocfile;
  FILE            *rocfp; 
  char            *sumfile;
  FILE            *sumfp; 

  double           hpts;    /* number of points in the histogram */
  double           bmin;    /* score histograms */
  double           w;
  double           pmass;
  double           fracfit;
  int              doexpfit;  // do an exponential fit, defautl is chi-square
  
  THRESH          *thresh;
  MODE             mode;

  int              window;
  int              slide;

  int              nofigures;

  ESL_MSA         *smsa;
  
  float            tol;
  int              verbose;
};

static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      1 },
  { "--outdir",     eslARG_STRING,       NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "specify a directory for all output files",                                                  1 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                1 },
  { "--r2rall",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "make R2R plot all position in the alignment",                                               1 },
  { "--window",       eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "window size",                                                                               1 },
  { "--slide",        eslARG_INT,      "50",     NULL,      "n>0",   NULL,    NULL,  NULL,               "window slide",                                                                              1 },
  { "--onemsa",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "if file has more than one msa, analyze only the first one",                                 1 },
  { "--nofigures",    eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "write .out and .sum files only",                                                            1 },
  { "--roc",          eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "write .roc file",                                                                           1 },
  { "--expo",         eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "true to do an exponential fit (default is gamma)",                                          0},
  { "--singlelink",   eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "true to use single linkage clustering (default is esl_msaweight_IDFilter)",                 0},
  /* E-values to assess significance */
  { "-E",            eslARG_REAL,      "0.05",   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "Eval: max expected number of covNBPs allowed",                                             1 },
 /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      1 },
  { "-I",             eslARG_REAL,     "1.0",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < <x> id",                                                             1 },
  { "-i",             eslARG_REAL,      NULL,    NULL, "0<=x<1.0",   NULL,    NULL,  NULL,               "require seqs to have >= <x> id",                                                            1 },
  { "--tstart",        eslARG_INT,      FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "min alignment position to analyze [1..alen]",                                               1 },
  { "--tend",          eslARG_INT,      FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "max alignment position to analyze [1..alen]",                                               1 },
  { "--consensus",    eslARG_NONE,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "analyze only consensus (seq_cons) positions",                                               1 },
  { "--submsa",       eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   1 },
  { "--nseqmin",      eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "minimum number of sequences in the alignment",                                              1 },  
  { "--gapthresh",    eslARG_REAL,     "0.5",    NULL,  "0<=x<=1",   NULL,    NULL,  NULL,               "keep columns with < <x> fraction of gaps",                                                  1 },
  { "--minid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      1 },
  { "--maxid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",                                                      1 },
  /* msa format */
  { "--informat",   eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            1 },
  /* null hypothesis */
  { "--nshuffle",      eslARG_INT,       NULL,   NULL,      "n>0",   NULL,    NULL,  NULL,               "number of shuffled sequences",                                                              1 },   
  /* covariation measures */
  { "--GTp",          eslARG_NONE,     "TRUE",   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   APC corrected statistic",                                                              1 },
  { "--MIp",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   APC corrected statistic",                                                              1 },
  /* alphabet type */
  { "--dna",          eslARG_NONE,      FALSE,   NULL,       NULL,      NULL, NULL,  NULL,               "use DNA alphabet",                                                                          0 },
  { "--rna",          eslARG_NONE,      FALSE,   NULL,       NULL,      NULL, NULL,  NULL,               "use RNA alphabet",                                                                          0 },
  { "--amino",        eslARG_NONE,      FALSE,   NULL,       NULL,      NULL, NULL,  NULL,               "use protein alphabet",                                                                      0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       1 },
  { "--outmsa",       eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "write actual msa used to file <f>,",                                                        1 },
  { "--outnull",      eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "write null alignments to file <f>,",                                                        1 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            1 },
 /* other options */  
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,      "42",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>. Use 0 for a random seed.",                                             1 },
  { "--fracfit",      eslARG_REAL,    "1.00",    NULL,   "0<x<=1",   NULL,    NULL,  NULL,               "pmass for censored histogram of cov scores",                                                0 },
  { "--pmass",        eslARG_REAL,    "0.05",    NULL,   "0<x<=1",   NULL,    NULL,  NULL,               "pmass for censored histogram of cov scores",                                                0 },
  { "--scmin",        eslARG_REAL,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "minimum score value considered",                                                            0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "Phylogenetic Covariation Above Random";

static int phcar_banner(FILE *fp, char *progname, char *banner);
static int MSA_banner(FILE *fp, char *msaname, MSA_STAT *mstat, MSA_STAT *omstat, int nbpairs, int onbpairs);
static int get_msaname(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int original_msa_manipulate(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **msa);
static int phcar_for_msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int run_phcar(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, RANKLIST **ret_ranklist, int analyze);
static int calculate_width_histo(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int null_phcar (ESL_GETOPTS *go, struct cfg_s *cfg, int nshuffle, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null2_phcar (ESL_GETOPTS *go, struct cfg_s *cfg, int nshuffle, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null_add2cumranklist(RANKLIST *ranklist, RANKLIST **ocumranklist, int verbose, char *errbuf);
static int create_msa_significant(struct cfg_s *cfg, ESL_MSA *msa, HITLIST *hitlist, ESL_MSA **ret_smsa);
static int msa_pid(int idx, ESL_MSA *msa, double **ret_pid, int **ret_nid, int **ret_n);

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, struct cfg_s *ret_cfg)
{
  ESL_GETOPTS  *go = esl_getopts_Create(options);
  struct cfg_s  cfg;
  char         *path;
  char         *s;
  char         *tok;
  char         *tok1 = NULL;
  struct stat  info;
  char         *outname = NULL;
  int           status;


  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  cfg.argc = argc;
  cfg.argv = argv;
 
  /* Phcar banner */
  phcar_banner(stdout, cfg.argv[0], banner);

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_usage(stdout,  cfg.argv[0], usage);
      if (puts("\noptions:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
     exit(0);
    }

  cfg.msafile = NULL;  
  if (esl_opt_ArgNumber(go) != 1) { 
    if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
 

  if ((cfg.msafile  = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));

  /* find gnuplot */
  cfg.gnuplot = NULL;
  path = getenv("PATH"); // chech for an executable in the path
  s = path; 
  while (esl_strcmp(s, "")) {
    esl_strtok(&s, ":", &tok);
    esl_sprintf(&tok1, "%s/gnuplot", tok);
    if ((stat(tok1, &info) == 0) && (info.st_mode & S_IXOTH)) {
      esl_sprintf(&cfg.gnuplot, "%s -persist", tok1);    
      break;
    }
    free(tok1); tok1 = NULL;
  }
  if (cfg.gnuplot == NULL && "GNUPLOT" && (s = getenv("GNUPLOT"))) { // check for an envvar
    if ((stat(s, &info) == 0) && (info.st_mode & S_IXOTH)) {
      esl_sprintf(&cfg.gnuplot, "%s -persist", s);
    }
  }
  
  /* outheader for all output files */
  cfg.outheader = NULL;
  msamanip_OutfileHeader(cfg.msafile, &cfg.outheader); 
  
  /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nmsa = 0;
  cfg.msafrq = NULL;
  cfg.msaname = NULL;
  cfg.msamap = NULL;

  /* alphabet */
  cfg.abc = NULL;
  cfg.abcisRNA = FALSE;
  if      (esl_opt_GetBoolean(go, "--rna"))   { cfg.abc = esl_alphabet_Create(eslRNA);   cfg.abcisRNA = TRUE;  }
  else if (esl_opt_GetBoolean(go, "--dna"))   { cfg.abc = esl_alphabet_Create(eslDNA);   cfg.abcisRNA = TRUE;  }
  else if (esl_opt_GetBoolean(go, "--amino")) { cfg.abc = esl_alphabet_Create(eslAMINO);                       }
  
  cfg.watch = esl_stopwatch_Create(); 
  
  cfg.outdir = NULL;
  if (esl_opt_IsOn(go, "--outdir")) esl_sprintf( &cfg.outdir, "%s", esl_opt_GetString(go, "--outdir"));
 
  esl_FileTail(cfg.msafile, TRUE, &cfg.filename);
  if ( cfg.outdir ) esl_sprintf( &cfg.outheader, "%s/%s", cfg.outdir, cfg.filename);
  
  if (esl_opt_IsOn(go, "--submsa")) { 
    cfg.submsa = esl_opt_GetInteger(go, "--submsa"); 
    esl_sprintf(&outname, "%s.select%d", cfg.outheader, cfg.submsa);  
    free(cfg.outheader); cfg.outheader = NULL;
    esl_sprintf(&cfg.outheader, "%s", outname); 
  }
  else { cfg.submsa = 0; }
    
  /* other options */
  cfg.consensus   = esl_opt_IsOn(go, "--consensus")?                                   TRUE : FALSE;
  cfg.maxsq_gsc   = 1000;
  cfg.nshuffle    = esl_opt_IsOn(go, "--nshuffle")?   esl_opt_GetInteger(go, "--nshuffle")  : -1.0;
  cfg.fragfrac    = esl_opt_IsOn(go, "-F")?           esl_opt_GetReal   (go, "-F")          :  0.0; // remove sequences with no residues always
  cfg.idthresh    = esl_opt_IsOn(go, "-I")?           esl_opt_GetReal   (go, "-I")          :  0.0;
  cfg.minidthresh = esl_opt_IsOn(go, "-i")?           esl_opt_GetReal   (go, "-i")          : -1.0;
  cfg.nseqmin     = esl_opt_IsOn(go, "--nseqmin")?    esl_opt_GetInteger(go, "--nseqmin")   : -1;
  cfg.gapthresh   = esl_opt_GetReal   (go, "--gapthresh");
  cfg.tstart      = esl_opt_IsOn(go, "--tstart")?     esl_opt_GetInteger(go, "--tstart")    : 0;
  cfg.tend        = esl_opt_IsOn(go, "--tend")?       esl_opt_GetInteger(go, "--tend")      : 0;
  cfg.tol         = esl_opt_GetReal   (go, "--tol");
  cfg.verbose     = esl_opt_GetBoolean(go, "-v");
  cfg.voutput     = esl_opt_GetBoolean(go, "--voutput");
  cfg.window      = esl_opt_IsOn(go, "--window")?     esl_opt_GetInteger(go, "--window")    : -1;
  cfg.slide       = esl_opt_IsOn(go, "--slide")?      esl_opt_GetInteger(go, "--slide")     : -1;
  cfg.onemsa      = esl_opt_IsOn(go, "--onemsa")?     esl_opt_GetBoolean(go, "--onemsa")    : FALSE;
  cfg.nofigures   = esl_opt_IsOn(go, "--nofigures")?  esl_opt_GetBoolean(go, "--nofigures") : FALSE;
  cfg.doroc       = esl_opt_IsOn(go, "--roc")?        esl_opt_GetBoolean(go, "--roc")       : FALSE;
  cfg.doexpfit    = esl_opt_IsOn(go, "--expo")?       esl_opt_GetBoolean(go, "--expo")      : FALSE;
  cfg.R2Rall      = esl_opt_GetBoolean(go, "--r2rall");
  cfg.singlelink  = esl_opt_GetBoolean(go, "--singlelink");
  cfg.method      = NULLPHYLO;
  
  if (cfg.minidthresh > cfg. idthresh) esl_fatal("minidthesh has to be smaller than idthresh");

  ESL_ALLOC(cfg.thresh, sizeof(THRESH));
  if (esl_opt_IsOn(go, "-E")) { cfg.thresh->type = Eval;    cfg.thresh->val = esl_opt_GetReal(go, "-E"); 
  }

  if      (esl_opt_GetBoolean(go, "--GTp"))   cfg.covtype = GTp;
  else if (esl_opt_GetBoolean(go, "--MIp"))   cfg.covtype = MIp;
  
  /* for the cov histograms */
  cfg.hpts    = HPTS;                                                               /* number of points in the histogram */
  cfg.bmin    = esl_opt_IsOn(go, "--scmin")? esl_opt_GetReal(go, "--scmin") : BMIN; /* lowest cov score to bound the histogram */
  cfg.w       = W;                                                                  /* default. histogram step, will be determined for each msa */
  cfg.pmass   = esl_opt_GetReal(go, "--pmass");
  cfg.fracfit = esl_opt_GetReal(go, "--fracfit");

  cfg.mstat  = NULL;
  cfg.omstat = NULL;
  cfg.smsa   = NULL;

  /* output file */
  if ( esl_opt_IsOn(go, "-o") ) {
    esl_sprintf(&cfg.outfile, "%s", esl_opt_GetString(go, "-o"));
    if ((cfg.outfp    = fopen(cfg.outfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outfile);
    esl_sprintf(&cfg.outsrtfile, "%s.sorted", esl_opt_GetString(go, "-o"));
    if ((cfg.outsrtfp = fopen(cfg.outsrtfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outsrtfile);
  } 
  else {
    if (cfg.window <= 0) esl_sprintf(&cfg.outfile, "%s.out", cfg.outheader);
    else                 esl_sprintf(&cfg.outfile, "%s.w%d.s%d.out", cfg.outheader, cfg.window, cfg.slide);
    if ((cfg.outfp    = fopen(cfg.outfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outfile);
    if (cfg.window <= 0) esl_sprintf(&cfg.outsrtfile, "%s.sorted.out", cfg.outheader);
    else                 esl_sprintf(&cfg.outsrtfile, "%s.w%d.s%d.sorted.out", cfg.outheader, cfg.window, cfg.slide);
    if ((cfg.outsrtfp = fopen(cfg.outsrtfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outsrtfile);
  }
  
  /*  rocplot file */
  cfg.rocfile = NULL;
  cfg.rocfp = NULL;
  if (cfg.doroc) {
    if (cfg.window <= 0) esl_sprintf(&cfg.rocfile, "%s.roc", cfg.outheader); 
    else                 esl_sprintf(&cfg.rocfile, "%s.w%d.s%d.roc", cfg.outheader, cfg.window, cfg.slide);
    if ((cfg.rocfp = fopen(cfg.rocfile, "w")) == NULL) esl_fatal("Failed to open rocfile %s", cfg.rocfile);
  }
  
  /*  summary file */
  if (cfg.window <= 0) esl_sprintf(&cfg.sumfile, "%s.sum", cfg.outheader); 
  else                 esl_sprintf(&cfg.sumfile, "%s.w%d.s%d.sum", cfg.outheader, cfg.window, cfg.slide);
  if ((cfg.sumfp = fopen(cfg.sumfile, "w")) == NULL) esl_fatal("Failed to open sumfile %s", cfg.sumfile);
  
  /* file with the msa actually used */
  cfg.outmsafile = NULL;
  cfg.outmsafp   = NULL;
  if (esl_opt_IsOn(go, "--outmsa")) {
    esl_sprintf(&cfg.outmsafile, "%s", esl_opt_GetString(go, "--outmsa"));
    if ((cfg.outmsafp = fopen(cfg.outmsafile, "w")) == NULL) esl_fatal("Failed to open outmsa file %s", cfg.outmsafile);
  } 
  esl_sprintf(&cfg.smsafile, "%s.smsa.sto", cfg.outheader);
  if ((cfg.smsafp = fopen(cfg.smsafile, "w")) == NULL) esl_fatal("Failed to open smsa file %s", cfg.smsafile);
    
  /* file with the null alignments */
  cfg.outnullfile = NULL;
  cfg.outnullfp   = NULL;
  if (esl_opt_IsOn(go, "--outnull")) {
    esl_sprintf(&cfg.outnullfile, "%s", esl_opt_GetString(go, "--outnull"));
    if ((cfg.outnullfp = fopen(cfg.outnullfile, "w")) == NULL) esl_fatal("Failed to open outnull file %s", cfg.outnullfile);
  } 
  
  /* msa-specific files */
  cfg.R2Rfile    = NULL;

  /* covhis file */
  cfg.covhisfile    = NULL;
  
  /* dotplot file */
  cfg.dplotfile    = NULL;

  cfg.onbpairs = 0;
  cfg.nbpairs  = 0;
  cfg.covhisfile = NULL;
    
  *ret_go  = go;
  *ret_cfg = cfg;

  if (outname) free(outname);
  if (tok1) free(tok1);
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, banner, usage);
  if (puts("\nwhere options are:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  esl_getopts_Destroy(go);
  if (outname) free(outname);
  if (tok1) free(tok1);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  if (outname) free(outname);
  if (tok1) free(tok1);
  exit(status);
}


static int
phcar_banner(FILE *fp, char *progname, char *banner)
{
  char *appname = NULL;
  int   status;

  if ((status = esl_FileTail(progname, FALSE, &appname)) != eslOK) return status;

  if (fprintf(fp, "# %s :: %s\n", appname, banner)                                               < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(fp, "# R-scape %s (%s)\n", PHCAR_VERSION, PHCAR_DATE)                              < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(fp, "# %s\n", PHCAR_COPYRIGHT)                                                     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(fp, "# %s\n", PHCAR_LICENSE)                                                       < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");

  if (appname) free(appname);
  return eslOK;

 ERROR:
  if (appname) free(appname);
  return status;
}


static int
MSA_banner (FILE *fp, char *msaname, MSA_STAT *mstat, MSA_STAT *omstat, int nbpairs, int onbpairs)
{
  if (omstat) 
    fprintf(fp, "# MSA %s nseq %d (%d) alen %" PRId64 " (%" PRId64 ") avgid %.2f (%.2f) nbpairs %d (%d)\n", 
	    msaname, mstat->nseq, omstat->nseq, mstat->alen, omstat->alen, 
	    mstat->avgid, omstat->avgid, nbpairs, onbpairs);
  else
    fprintf(fp, "# wMSA %s nseq %d alen %" PRId64 " avgid %.2f nbpairs %d\n", 
	    msaname, mstat->nseq, mstat->alen, mstat->avgid, nbpairs);

  return eslOK;
}

int
main(int argc, char **argv)
{ 
  ESL_GETOPTS     *go;
  struct cfg_s     cfg;
  char            *omsaname = NULL;
  ESL_MSAFILE    *afp = NULL;
  ESL_MSA         *msa = NULL;            /* the input alignment    */
  ESL_MSA         *wmsa = NULL;           /* the window alignment   */
  int             *useme = NULL;
  int              first, last;
  int              i;
  int              status = eslOK;
  int              hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* Open the MSA file */
  status = esl_msafile_Open(&cfg.abc, cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  esl_alphabet_SetEquiv(cfg.abc, '=', '-');     /* allow = as a gap character too */
  esl_alphabet_SetEquiv(cfg.abc, '.', '-');     /* allow . as a gap character too */

  /* read the MSA */
  while ((hstatus = esl_msafile_Read(afp, &msa)) != eslEOF) {
    if (hstatus != eslOK) { printf("%s\n", afp->errmsg) ; esl_msafile_ReadFailure(afp, status); }
    cfg.nmsa ++;
    if (cfg.onemsa && cfg.nmsa > 1) break;

    /* the msaname */
    status = get_msaname(go, &cfg, msa);
    if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to manipulate alignment"); }
    cfg.abcisRNA = FALSE;
    if (msa->abc->type == eslDNA || msa->abc->type == eslRNA) { cfg.abcisRNA = TRUE; }
    cfg.covclass = C16;
    
    if (cfg.window > 0) {
 
      esl_sprintf(&omsaname, "%s", cfg.msaname);

      useme = malloc(sizeof(int)*(msa->alen));
      for (first = 1; first <= ESL_MAX(1, msa->alen); first += cfg.slide) {

	esl_vec_ISet(useme, msa->alen, FALSE);
	wmsa = esl_msa_Clone(msa);

	last = ESL_MIN(first+cfg.window-1, msa->alen);	
	for (i = first-1; i < last; i ++) useme[i] = TRUE;
	status = esl_msa_ColumnSubset(wmsa, cfg.errbuf, useme);
	if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to manipulate alignment"); }

	/* add the name including the from-to information */
	free(cfg.msaname); cfg.msaname = NULL;
	esl_sprintf(&cfg.msaname, "%s_%d-%d", omsaname, first, last);
	esl_msa_SetName(wmsa, cfg.msaname, -1);

	status = original_msa_manipulate(go, &cfg, &wmsa);
	if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to manipulate alignment"); }
	if (wmsa == NULL) continue;
	if (wmsa->alen <= 0) {
	  MSA_banner(cfg.outfp, cfg.msaname, cfg.mstat, cfg.omstat, 0.0, 0.0);
	  continue;
	}

	cfg.firstpos = first;

	status = phcar_for_msa(go, &cfg, wmsa);
	if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to run phcar"); }
 
	esl_msa_Destroy(wmsa); wmsa = NULL;
	if (last >= msa->alen) break;
      }
    }
    else {      
      status = original_msa_manipulate(go, &cfg, &msa);
      if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to manipulate alignment"); }
      if (msa == NULL) continue;
      if (msa->alen == 0) {
	MSA_banner(cfg.outfp, cfg.msaname, cfg.mstat, cfg.omstat, 0.0, 0.0);
	continue;
      }

      cfg.firstpos = 1;
      status = phcar_for_msa(go, &cfg, msa);
      if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to run phcar"); }
    }

    if (omsaname) free(omsaname); omsaname = NULL;
    if (useme) free(useme); useme = NULL;
    if (msa) esl_msa_Destroy(msa); msa = NULL;
    if (cfg.msaname) free(cfg.msaname); cfg.msaname = NULL;
    if (cfg.msamap) free(cfg.msamap); cfg.msamap = NULL;
    if (cfg.omstat) free(cfg.omstat); cfg.omstat = NULL;
    if (cfg.mstat) free(cfg.mstat); cfg.mstat = NULL;
    if (cfg.smsa) esl_msa_Destroy(cfg.smsa); cfg.smsa = NULL;
  }

  /* cleanup */
  fclose(cfg.outfp);
  fclose(cfg.outsrtfp);
  if (cfg.rocfp) fclose(cfg.rocfp);
  fclose(cfg.sumfp);
  if (cfg.outmsafp) fclose(cfg.outmsafp);
  fclose(cfg.smsafp);
  if (cfg.outnullfp) fclose(cfg.outnullfp);
  free(cfg.filename);
  esl_stopwatch_Destroy(cfg.watch);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  esl_msafile_Close(afp);
  if (cfg.msaname) free(cfg.msaname);
  if (cfg.outfile) free(cfg.outfile);
  if (cfg.outsrtfile) free(cfg.outsrtfile);
  if (cfg.outdir) free(cfg.outdir);
  free(cfg.outheader);
  if (cfg.rocfile) free(cfg.rocfile);
  free(cfg.sumfile);
  free(cfg.gnuplot);
  if (cfg.outmsafile) free(cfg.outmsafile);
   free(cfg.smsafile);
  if (cfg.outnullfile) free(cfg.outnullfile);
  if (cfg.thresh) free(cfg.thresh);
  if (useme) free(useme);
  if (cfg.msamap) free(cfg.msamap); 
  if (cfg.omstat) free(cfg.omstat);
  if (cfg.mstat) free(cfg.mstat);
  if (cfg.smsa) esl_msa_Destroy(cfg.smsa);

  return 0;
}

static int
get_msaname(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  char    *msg = "get_msaname failed";
  char *type = NULL;
  char *tp;
  char *tok1;
  char *tok2;
  char *tok = NULL;
  char *submsaname = NULL;
  int   t;
  
  /* the msaname */
  for (t = 0; t < msa->ngf; t++) {
    if (!esl_strcmp(msa->gf_tag[t], "TP")) {
      tp = msa->gf[t];	
      while (*tp != '\0') {
	if (type) free(type); type = NULL;
	if (tok) free(tok); tok = NULL;
	if (esl_strtok(&tp,   ";", &tok1) != eslOK) esl_fatal(msg);
	if (esl_strtok(&tok1, " ", &tok2) != eslOK) esl_fatal(msg);
	esl_sprintf(&tok, "_%s", tok2);
	esl_strcat(&type, -1, tok, -1);
      }
    }
  }
  if      (msa->acc && msa->name) esl_sprintf(&cfg->msaname, "%s_%s", msa->acc, msa->name, type);
  else if (msa->name)             esl_sprintf(&cfg->msaname, "%s",              msa->name, type);
  else if (msa->acc)              esl_sprintf(&cfg->msaname, "%s",    msa->acc);
  else if (cfg->onemsa)           esl_sprintf(&cfg->msaname, "%s",    cfg->filename);
  else                            esl_sprintf(&cfg->msaname, "%s_%d", cfg->filename, cfg->nmsa);
  if (esl_opt_IsOn(go, "--submsa")) {					       
    esl_sprintf(&submsaname, "%s.select%d", cfg->msaname, esl_opt_GetInteger(go, "--submsa"));
    free(cfg->msaname); cfg->msaname = NULL;
    esl_sprintf(&cfg->msaname, "%s", submsaname);
  }

  if (tok) free(tok);
  if (type) free(type);
  if (submsaname) free(submsaname);
  return eslOK;
}

static int
original_msa_manipulate(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa)
{
  ESL_MSA *msa = *omsa;
  int     *useme = NULL;
  char    *msg = "original_msa_manipulate failed";
  char    *type = NULL;
  char    *tok = NULL;
  char    *submsaname = NULL;
  int      alen = msa->alen;
  int64_t  tstart, tend;
  int64_t  startpos, endpos;
  int      seq_cons_len = 0;
  int      nremoved = 0;	  /* # of identical sequences removed */
  int      nfrags = 0;	          /* # of fragments removed */

  /* stats of the original alignment */
  msamanip_XStats(msa, &cfg->omstat);
  msamanip_CalculateCT(msa, NULL, &cfg->onbpairs, -1., cfg->errbuf);
  /* print some info */
  if (cfg->verbose) {
    fprintf(stdout, "Given alignment\n");
    fprintf(stdout, "%6d %d          %s\n", msa->nseq, (int)msa->alen, cfg->msafile);
    if (esl_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write msa"); 
    msamanip_DumpStats(stdout, msa, cfg->omstat); 
  }
  
  /* apply msa filters and then select submsa
   * none of these functions reduce the number of columns in the alignemnt
   */
  if (cfg->fragfrac >= 0.    && msamanip_RemoveFragments(cfg->fragfrac, omsa, &nfrags, &seq_cons_len)                 != eslOK) {
    printf("%s\nremove_fragments failed\n", cfg->errbuf);                 esl_fatal(msg); }
  if (esl_opt_IsOn(go, "-I") && msamanip_SelectSubsetBymaxID(cfg->r, omsa, cfg->idthresh, cfg->singlelink, &nremoved) != eslOK) {
    printf("%s\n", cfg->errbuf); printf("select_subsetBymaxID failed\n"); esl_fatal(msg); }
  if (esl_opt_IsOn(go, "-i") && msamanip_SelectSubsetByminID(cfg->r, omsa, cfg->minidthresh, &nremoved)               != eslOK) {
    printf("%s\n", cfg->errbuf); printf("select_subsetByminID failed\n"); esl_fatal(msg); }
  if (cfg->submsa            && msamanip_SelectSubset(cfg->r, cfg->submsa, omsa, NULL, cfg->errbuf, cfg->verbose)     != eslOK) {
    printf("%s\n", cfg->errbuf);                                          esl_fatal(msg); }
  
  msa = *omsa;
  if (msa->alen != alen) {
    printf("filtering altered the length of the alignemnt!\n");
    esl_fatal(msg);
  }
  
  if (msa == NULL) {
    if (submsaname) free(submsaname);
    if (type) free(type); type = NULL;
    if (tok) free(tok); tok = NULL;
    free(cfg->msaname); cfg->msaname = NULL;
    return eslOK;
  }
  if (msa->nseq < cfg->nseqmin || msa->alen == 0) {
    esl_msa_Destroy(msa); msa = NULL;
    if (submsaname) free(submsaname);
    if (type) free(type); type = NULL;
    if (tok) free(tok); tok = NULL;
    free(cfg->msaname); cfg->msaname = NULL;
    return eslOK;
  }
  
  /* define [tstart,tend]  */
  if      (cfg->tstart == 0 && cfg->tend == 0) { tstart = 1;           tend = msa->alen; }
  else if (cfg->tstart >  0 && cfg->tend == 0) { tstart = cfg->tstart; tend = msa->alen; }
  else if (cfg->tstart == 0 && cfg->tend >  0) { tstart = 1;           tend = cfg->tend; }
  else                                         { tstart = cfg->tstart; tend = cfg->tend; }
  
  /* add tstat..tend information */
  if (tstart > 1 || tend < msa->alen) 
    esl_sprintf(&cfg->msaname, "%s_%d-%d", cfg->msaname, tstart, tend);
  
  /* Now apply [tstart,tend] restriction if given */
  cfg->omstat->alen = tend - tstart + 1;
  if (msamanip_Truncate(msa, tstart, tend, &startpos, &endpos, cfg->errbuf) != eslOK) {
    printf("%s\nTruncate failed\n", cfg->errbuf);        esl_fatal(msg); }
  // check if empty sequences can be removed after the truncation
  if (msamanip_RemoveFragments(cfg->fragfrac, omsa, &nfrags, &seq_cons_len)           != eslOK) {
    printf("%s\nremove_fragments failed\n", cfg->errbuf); esl_fatal(msg); }

  msa = *omsa;
   /* remove columns with gaps.
   * Important: the mapping is done here; cannot remove any other columns beyond this point.
   */
  if (cfg->consensus) {
    if (msamanip_SelectConsensus(msa, &useme, cfg->verbose) != eslOK) {
      printf("%s\nconsensus selection fails\n", cfg->errbuf); esl_fatal(msg);
    }
  }
  if (msamanip_RemoveGapColumns(cfg->gapthresh, msa, startpos, endpos, alen, &cfg->msamap, NULL, &useme, cfg->errbuf, cfg->verbose) != eslOK) {
    printf("%s\n", cfg->errbuf); esl_fatal(msg);
  }
  /* convert degenerates to N, and Missing/Nonresidues to Gap */
  msamanip_ConvertDegen2N(msa);
  msamanip_ConvertMissingNonresidue2Gap(msa);
  
  /* given msa aveid and avematch */
  msamanip_XStats(msa, &cfg->mstat);
  
  if ((esl_opt_IsOn(go, "--minid") && cfg->mstat->avgid < 100.*esl_opt_GetReal(go, "--minid")) ||
      (esl_opt_IsOn(go, "--maxid") && cfg->mstat->avgid > 100.*esl_opt_GetReal(go, "--maxid"))   ) {
    esl_msa_Destroy(msa); msa = NULL;
    if (type) free(type); type = NULL;
    if (tok) free(tok); tok = NULL;
    free(cfg->msaname); cfg->msaname = NULL;
    if (submsaname) free(submsaname);
    return eslOK;  
  }

  /* add the name */
  esl_msa_SetName(msa, cfg->msaname, -1);
  
  /* print some info */
  if (cfg->verbose) {
    fprintf(stdout, "Used alignment\n");
    fprintf(stdout, "%6d          %s\n", msa->nseq, cfg->msafile);
    if (msa->alen > 0 && esl_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write msa"); 
    msamanip_DumpStats(stdout, msa, cfg->mstat); 
  }

  if (tok) free(tok);
  if (type) free(type);
  if (submsaname) free(submsaname);
  if (useme) free(useme);
  return eslOK;
}

static int
phcar_for_msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  RANKLIST *ranklist_null = NULL;
  RANKLIST *ranklist_aux  = NULL;
  int       nshuffle;
  int       analyze;
  int       status;
  
  if (msa == NULL) return eslOK;

  if (msa->nseq <= 1) {
    MSA_banner(cfg->outfp,    cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
    MSA_banner(cfg->outsrtfp, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
    return eslOK; 
  }
 
  /* outmsa file if requested */
  if (cfg->outmsafp) esl_msafile_Write(cfg->outmsafp, msa, eslMSAFILE_STOCKHOLM);

  if (cfg->outdir) {
    /* covhis file */
    esl_sprintf(&cfg->covhisfile,    "%s/%s.surv",     cfg->outdir, cfg->msaname);
  }
  else {
    /* covhis file */
    esl_sprintf(&cfg->covhisfile,    "%s.surv",     cfg->msaname);
  }
  
  /* R2R annotated sto file */
  if (cfg->outdir && !cfg->nofigures) {
    esl_sprintf(&cfg->R2Rfile,    "%s/%s.R2R.sto",     cfg->outdir, cfg->msaname);
    
    /* dotplot file */
    esl_sprintf(&cfg->dplotfile,    "%s/%s.dplot",     cfg->outdir, cfg->msaname);
  }
  else if (!cfg->nofigures) {
    esl_sprintf(&cfg->R2Rfile,    "%s.R2R.sto",     cfg->msaname);
        
    /* dotplot file */
    esl_sprintf(&cfg->dplotfile,    "%s.dplot",     cfg->msaname);
  }

  /* the ct vector */
  status = msamanip_CalculateCT(msa, &cfg->ct, &cfg->nbpairs, -1.0, cfg->errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s. Failed to calculate ct vector", cfg->errbuf);

#if 0
  msamanip_CalculateBC(msa, cfg->ct, &cfg->ft, &cfg->fbp, &cfg->fnbp, cfg->errbuf);
  esl_vec_DDump(stdout, cfg->ft,   cfg->abc->K, "total BC");
  esl_vec_DDump(stdout, cfg->fbp,  cfg->abc->K, "basepairs BC");
  esl_vec_DDump(stdout, cfg->fnbp, cfg->abc->K, "nonbasepairs BC");
#endif
  
  // Print some alignment information
  MSA_banner(stdout, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);

  /* the null model first */
  nshuffle = cfg->nshuffle;
  if (nshuffle < 0) {
    nshuffle = 20;
    if (msa->nseq*msa->alen < 1e3) { nshuffle = 200; }
  }
  if (msa->nseq*msa->alen < 1e3) { cfg->fracfit = 0.3; }
  
  cfg->mode = RANSS;
  status = null_phcar(go, cfg, nshuffle, msa, &ranklist_null);
  if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null phcar", cfg->errbuf);
  
  /* main function */
  cfg->mode = GIVSS;
  analyze = TRUE;
  status = run_phcar(go, cfg, msa, ranklist_null, ranklist_aux, NULL, analyze);
  if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s\n", cfg->errbuf);

  if (cfg->ct) free(cfg->ct); cfg->ct = NULL;
  if (cfg->msafrq) free(cfg->msafrq); cfg->msafrq = NULL;
  if (ranklist_null) cov_FreeRankList(ranklist_null); ranklist_null = NULL;
  if (ranklist_aux) cov_FreeRankList(ranklist_aux); ranklist_aux = NULL;
  
  if (cfg->covhisfile) free(cfg->covhisfile); 
  if (cfg->dplotfile) free(cfg->dplotfile);
  if (cfg->R2Rfile) free(cfg->R2Rfile);
 
  return eslOK;

 ERROR:
  if (cfg->ct) free(cfg->ct); 
  if (cfg->msafrq) free(cfg->msafrq); 
  if (cfg->msaname) free(cfg->msaname);
  if (ranklist_null) cov_FreeRankList(ranklist_null); 
  if (ranklist_aux) cov_FreeRankList(ranklist_aux);
  if (cfg->covhisfile) free(cfg->covhisfile); 
  if (cfg->dplotfile) free(cfg->dplotfile);
  if (cfg->R2Rfile) free(cfg->R2Rfile);
  return status;
}

static int
calculate_width_histo(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  struct mutual_s *mi = NULL;
  struct data_s    data;
  int              status;

  /* weigth the sequences */
  if (msa->nseq <= cfg->maxsq_gsc) esl_msaweight_GSC(msa, NULL);
  else                             esl_msaweight_PB(msa);

  /* create the MI structure */
  mi = corr_Create(msa->alen, msa->nseq, TRUE, msa->nseq, msa->alen, cfg->abc, cfg->covclass);

  /* main function */
  data.outfp         = NULL;
  data.outsrtfp      = NULL;
  data.rocfp         = NULL;
  data.sumfp         = NULL;
  data.dplotfile     = NULL;
  data.R2Rfile       = NULL;
  data.R2Rall        = FALSE;
  data.gnuplot       = NULL;
  data.r             = cfg->r;
  data.ranklist_null = NULL;
  data.ranklist_aux  = NULL;
  data.mi            = mi;
  data.covtype       = cfg->covtype;
  data.thresh        = cfg->thresh;
  data.method        = cfg->method;
  data.mode          = cfg->mode;
  data.ct            = cfg->ct;
  data.onbpairs      = cfg->onbpairs;
  data.nbpairs       = cfg->nbpairs;
  data.msamap        = cfg->msamap;
  data.bmin          = cfg->bmin;
  data.w             = cfg->w;
  data.fracfit       = cfg->fracfit;
  data.pmass         = cfg->pmass;
  data.doexpfit      = cfg->doexpfit;
  data.tol           = cfg->tol;
  data.nofigures     = cfg->nofigures;
  data.verbose       = cfg->verbose;
  data.errbuf        = cfg->errbuf;
  data.donull2b      = FALSE;
  data.ignorebps     = TRUE;

  status = cov_Calculate(&data, msa, NULL, NULL, FALSE);   
  if (status != eslOK) goto ERROR; 
  if (mi->maxCOV <= cfg->bmin) ESL_XFAIL(eslFAIL, cfg->errbuf, "bmin %f should be larger than maxCOV %f\n", cfg->bmin, mi->maxCOV);

  cfg->w = (mi->maxCOV - ESL_MAX(data.bmin,mi->minCOV)) / (double) cfg->hpts;
  
  if (cfg->w < cfg->tol) cfg->w = cfg->tol;
  if (cfg->verbose) printf("w %f minCOV %f bmin %f maxCOV %f\n", cfg->w, mi->minCOV, cfg->bmin, mi->maxCOV);
  if (cfg->w <= 0) return eslFAIL;

  corr_Destroy(mi);
  return eslOK;

 ERROR:
  if (mi) corr_Destroy(mi);
  return status;
}

static int
run_phcar(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, RANKLIST **ret_ranklist, int analyze)
{
  char            *title = NULL;
  struct data_s    data;
  struct mutual_s *mi   = NULL;
  RANKLIST        *ranklist = NULL;
  HITLIST         *hitlist = NULL;
  double          *pid = NULL;
  int             *nid = NULL;
  int             *n = NULL;
  int              nnodes;
  int              status;

  esl_stopwatch_Start(cfg->watch);

  /* weigth the sequences */
  if (msa->nseq <= cfg->maxsq_gsc) esl_msaweight_GSC(msa, NULL);
  else                             esl_msaweight_PB(msa);

  /* print to stdout */
  if (cfg->verbose) {
    MSA_banner(stdout, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
  }
   
  if (cfg->mode != RANSS) {
    MSA_banner(cfg->outfp,    cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
    MSA_banner(cfg->outsrtfp, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
    esl_sprintf(&title, "%s (seqs %d alen %" PRId64 " avgid %d bpairs %d)", 
		cfg->msaname, msa->nseq, msa->alen, (int)ceil(cfg->mstat->avgid), cfg->nbpairs);
  }

  /* create the MI structure */
  mi = corr_Create(msa->alen, msa->nseq, (cfg->mode == RANSS)? TRUE : FALSE, msa->nseq, msa->alen, cfg->abc, cfg->covclass);

  /* write MSA info to the sumfile */
  if (cfg->mode != RANSS) {
    fprintf(cfg->sumfp, "#target_E-val\tMSA\tnseq\talen\tavgid\tTP\tTrue\tFound\tSEN\tPPV\n");
    fprintf(cfg->sumfp, "%g\t\t%s\t%d\t%d\t%.2f\t", cfg->thresh->val, cfg->msaname, msa->nseq, (int)msa->alen, cfg->mstat->avgid); 
  }
  
  /* write MSA info to the rocfile */
  if (cfg->doroc) {
    if (cfg->mode == GIVSS)
      fprintf(cfg->rocfp, "# MSA nseq %d alen %" PRId64 " avgid %f nbpairs %d (%d)\n",
	      msa->nseq, msa->alen, cfg->mstat->avgid, cfg->nbpairs, cfg->onbpairs);
  }
  
  /* main function */
  data.outfp         = (cfg->mode == RANSS)? NULL : cfg->outfp;
  data.outsrtfp      = (cfg->mode == RANSS)? NULL : cfg->outsrtfp;
  data.rocfp         = cfg->rocfp;
  data.sumfp         = (cfg->mode == RANSS)? NULL : cfg->sumfp;
  data.dplotfile     = cfg->dplotfile;
  data.R2Rfile       = cfg->R2Rfile;
  data.R2Rall        = cfg->R2Rall;
  data.gnuplot       = cfg->gnuplot;
  data.r             = cfg->r;
  data.ranklist_null = ranklist_null;
  data.ranklist_aux  = ranklist_aux;
  data.mi            = mi;
  data.covtype       = cfg->covtype;
  data.thresh        = cfg->thresh;
  data.method        = cfg->method;
  data.mode          = cfg->mode;
  data.ct            = cfg->ct;
  data.onbpairs      = cfg->onbpairs;
  data.nbpairs       = cfg->nbpairs;
  data.msamap        = cfg->msamap;
  data.firstpos      = cfg->firstpos;
  data.bmin          = cfg->bmin;
  data.w             = cfg->w;
  data.fracfit       = cfg->fracfit;
  data.pmass         = cfg->pmass;
  data.doexpfit      = cfg->doexpfit;
  data.tol           = cfg->tol;
  data.nofigures     = cfg->nofigures;
  data.verbose       = cfg->verbose;
  data.errbuf        = cfg->errbuf;
  data.donull2b      = FALSE;
  data.ignorebps     = TRUE;

  status = cov_Calculate(&data, msa, &ranklist, &hitlist, analyze);   
  if (status != eslOK) goto ERROR; 
  if (cfg->mode == GIVSS && (cfg->verbose)) cov_DumpRankList(stdout, ranklist);

  if (cfg->mode == GIVSS) {
    if (cfg->verbose) {
      printf("score total distribution\n");
      printf("imin %d imax %d xmax %f xmin %f width %f\n",
	     ranklist->ha->imin, ranklist->ha->imax, ranklist->ha->xmax, ranklist->ha->xmin, ranklist->ha->w);
      printf("score truncated distribution\n");
      printf("imin %d imax %d xmax %f xmin %f width %f\n",
	     ranklist->ht->imin, ranklist->ht->imax, ranklist->ht->xmax, ranklist->ht->xmin, ranklist->ht->w);
    }
    status = cov_WriteHistogram(&data, cfg->gnuplot, cfg->covhisfile, NULL, ranklist, title);
    if (status != eslOK) goto ERROR;

    status = create_msa_significant(cfg, msa, hitlist, &cfg->smsa);
    if (status != eslOK) goto ERROR;
    esl_msafile_Write(stdout, cfg->smsa, eslMSAFILE_STOCKHOLM);
    esl_msafile_Write(cfg->smsafp, cfg->smsa, eslMSAFILE_STOCKHOLM);
    msa_pid(0, cfg->smsa, &pid, &nid, &n);
  }
 
  if (ret_ranklist) *ret_ranklist = ranklist; else if (ranklist) cov_FreeRankList(ranklist);
  if (hitlist) cov_FreeHitList(hitlist); hitlist = NULL;
  if (title) free(title);
  if (pid) free(pid);
  if (nid) free(nid);
  if (n) free(n);
  corr_Destroy(mi); mi = NULL;

  return eslOK;
  
 ERROR:
  if (mi)       corr_Destroy(mi);
  if (ranklist) cov_FreeRankList(ranklist);
  if (hitlist)  cov_FreeHitList(hitlist);
  if (title)    free(title);
  if (pid)      free(pid);
  if (nid)      free(nid);
  if (n)        free(n);
  return status;
}

/* Fitch algorithm on a star topology */
static int
null2_phcar(ESL_GETOPTS *go, struct cfg_s *cfg, int nshuffle, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *allmsa = NULL;
  ESL_MSA   *shmsa = NULL;
  MSA_STAT  *shmstat = NULL;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int       *usecol = NULL;
  int       sc;
  int       s;
  int       status;
  
  // use all columns to do the shuffling
  ESL_ALLOC(usecol, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(usecol, msa->alen+1, TRUE);

  for (s = 0; s < nshuffle; s ++) {

    status = Tree_FitchAlgorithmAncenstral(cfg->r, NULL, msa, &allmsa, &sc, cfg->errbuf, cfg->verbose);
    if (status != eslOK) goto ERROR;
    if (cfg->verbose) {
      esl_msafile_Write(stdout, allmsa, eslMSAFILE_STOCKHOLM); 
      printf("fitch sc %d\n", sc);
    }
    
    status = msamanip_ShuffleTreeSubstitutions(cfg->r, NULL, msa, allmsa, usecol, &shmsa, cfg->errbuf, cfg->verbose);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s.\nFailed to run null phcar", cfg->errbuf);

    /* output null msas to file if requested */
    if (cfg->outnullfp) esl_msafile_Write(cfg->outnullfp, shmsa, eslMSAFILE_STOCKHOLM);
    
    if (cfg->verbose) {
       esl_msafile_Write(stdout, shmsa, eslMSAFILE_STOCKHOLM); 
      //msamanip_XStats(shmsa, &shmstat);
      //msamanip_DumpStats(stdout, shmsa, shmstat);
    }

    if (s == 0) {
      status = calculate_width_histo(go, cfg, shmsa);
      if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s.\nFailed to calculate the width of the histogram", cfg->errbuf);
    }
    status = run_phcar(go, cfg, shmsa, NULL, NULL, &ranklist, TRUE);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s.\nFailed to run null2 phcar", cfg->errbuf);
    if (shmsa == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "error creating shmsa");
   
    status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
    if (status != eslOK) goto ERROR;

    esl_msa_Destroy(allmsa); allmsa = NULL;
    esl_msa_Destroy(shmsa); shmsa = NULL;
    cov_FreeRankList(ranklist); ranklist = NULL;
  }
  
  if (cfg->verbose) {
    printf("null distribution - cumulative\n");
    printf("imin %d imax %d xmax %f xmin %f\n", 
	   cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
    //esl_histogram_Plot(stdout, cumranklist->ha);
    esl_histogram_PlotSurvival(stdout, cumranklist->ha);
  }
  
  if (cfg->verbose) cov_DumpRankList(stdout, cumranklist);
  
  *ret_cumranklist = cumranklist;

  if (allmsa) esl_msa_Destroy(allmsa);
  free(usecol);
  free(shmstat);
  return eslOK;
  
 ERROR:
  if (allmsa) esl_msa_Destroy(allmsa);
  if (shmsa) esl_msa_Destroy(shmsa);
  if (usecol) free(usecol);
  if (shmstat) free(shmstat);
  if (ranklist) cov_FreeRankList(ranklist);
  if (cumranklist) cov_FreeRankList(cumranklist);
  return status;
}

/* suffle the residues in each sequence independently */
static int
null_phcar(ESL_GETOPTS *go, struct cfg_s *cfg, int nshuffle, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *shmsa = NULL;
  MSA_STAT  *shmstat = NULL;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int        n;
  int        sc;
  int        s;
  int        status;
  
  for (s = 0; s < nshuffle; s ++) {

    shmsa = esl_msa_Clone(msa);
    for (n = 0; n < msa->nseq; n ++) {
      esl_rsq_XShuffle(cfg->r, shmsa->ax[n], shmsa->alen, shmsa->ax[n]);
    }

    /* output null msas to file if requested */
    if (cfg->outnullfp) esl_msafile_Write(cfg->outnullfp, shmsa, eslMSAFILE_STOCKHOLM);
    
    if (cfg->verbose) {
      esl_msafile_Write(stdout, shmsa, eslMSAFILE_STOCKHOLM); 
      //msamanip_XStats(shmsa, &shmstat);
      //msamanip_DumpStats(stdout, shmsa, shmstat);
    }

    if (s == 0) {
      status = calculate_width_histo(go, cfg, shmsa);
      if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s.\nFailed to calculate the width of the histogram", cfg->errbuf);
    }
    status = run_phcar(go, cfg, shmsa, NULL, NULL, &ranklist, TRUE);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s.\nFailed to run null phcar", cfg->errbuf);
    if (shmsa == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "error creating shmsa");
   
    status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
    if (status != eslOK) goto ERROR;

    esl_msa_Destroy(shmsa); shmsa = NULL;
    cov_FreeRankList(ranklist); ranklist = NULL;
  }
  
  if (cfg->verbose) {
    printf("null distribution - cumulative\n");
    printf("imin %d imax %d xmax %f xmin %f\n", 
	   cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
    //esl_histogram_Plot(stdout, cumranklist->ha);
    esl_histogram_PlotSurvival(stdout, cumranklist->ha);
  }
  
  if (cfg->verbose) cov_DumpRankList(stdout, cumranklist);
  
  *ret_cumranklist = cumranklist;

  free(shmstat);
  return eslOK;
  
 ERROR:
  if (shmsa) esl_msa_Destroy(shmsa);
  if (shmstat) free(shmstat);
  if (ranklist) cov_FreeRankList(ranklist);
  if (cumranklist) cov_FreeRankList(cumranklist);
  return status;
}

static int
null_add2cumranklist(RANKLIST *ranklist, RANKLIST **ocumranklist, int verbose, char *errbuf)
{
  RANKLIST *cumranklist = *ocumranklist;
  int       b;
  int       cumb;
  
  if (*ocumranklist == NULL) {
    *ocumranklist = cov_CreateRankList(ranklist->ha->bmax, ranklist->ha->bmin, ranklist->ha->w);
    cumranklist = *ocumranklist;
    cumranklist->ha->n     = ranklist->ha->n;
    cumranklist->ha->xmin  = ranklist->ha->xmin;
    cumranklist->ha->xmax  = ranklist->ha->xmax;
    cumranklist->ha->imin  = ranklist->ha->imin;
    cumranklist->ha->imax  = ranklist->ha->imax;
  }
  else {                    
    cov_GrowRankList(ocumranklist, ranklist->ha->bmax, ranklist->ha->bmin);
    cumranklist = *ocumranklist;
    cumranklist->ha->n    += ranklist->ha->n;
    cumranklist->ha->xmin  = ESL_MIN(cumranklist->ha->xmin, ranklist->ha->xmin);
    cumranklist->ha->xmax  = ESL_MAX(cumranklist->ha->xmax, ranklist->ha->xmax);
    cumranklist->ha->imin  = ESL_MIN(cumranklist->ha->imin, ranklist->ha->imin);
    cumranklist->ha->imax  = ESL_MAX(cumranklist->ha->imax, ranklist->ha->imax);
  }
 
  for (b = ranklist->ha->imin; b <= ranklist->ha->imax; b ++) {
    cov_ranklist_Bin2Bin(b, ranklist->ha, cumranklist->ha, &cumb);
    
    if (cumb < cumranklist->ha->nb) {
      cumranklist->ha->obs[cumb] += ranklist->ha->obs[b];
      cumranklist->ha->Nc        += ranklist->ha->obs[b];
      cumranklist->ha->No        += ranklist->ha->obs[b];
    }
  }
  
  if (verbose) {
    printf("null distribution \n");
    printf("imin %d imax %d xmax %f xmin %f | imin %d imax %d xmax %f xmin %f\n", 
	   ranklist->ha->imin, ranklist->ha->imax, ranklist->ha->xmax, ranklist->ha->xmin,
	   cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
    //esl_histogram_Plot(stdout, ranklist->ha);
    //esl_histogram_PlotSurvival(stdout, ranklist->ha);
  }

 return eslOK;
}


static int
create_msa_significant(struct cfg_s *cfg, ESL_MSA *msa, HITLIST *hitlist, ESL_MSA **ret_smsa)
{
  ESL_MSA *smsa = NULL;
  int     *useme = NULL;
  int      h;
  int      status;

  ESL_ALLOC(useme, sizeof(int)*msa->alen);
  esl_vec_ISet(useme, msa->alen, FALSE);

  for (h = 0; h < hitlist->nhit; h++) {
    useme[hitlist->hit[h].i-1] = TRUE;
    useme[hitlist->hit[h].j-1] = TRUE;
  }
  smsa = esl_msa_Clone(msa);
  esl_msa_ColumnSubset(smsa, cfg->errbuf, useme);
 
  free(useme);
  *ret_smsa = smsa;
  return eslOK;

 ERROR:
  return status;
}

static int
msa_pid(int idx, ESL_MSA *msa, double **ret_pid, int **ret_nid, int **ret_n)
{
  ESL_DSQ *ref;
  double  *pid = NULL;
  double  *cid = NULL;
  int     *nid = NULL;
  int     *n = NULL;
  int     *useme = NULL;
  int     *perm = NULL;
  int      s, s1, s2;
  int      status;

  ESL_ALLOC(pid,   sizeof(double)*msa->nseq);
  ESL_ALLOC(cid,   sizeof(double)*msa->nseq);
  ESL_ALLOC(nid,   sizeof(int)   *msa->nseq);
  ESL_ALLOC(n,     sizeof(int)   *msa->nseq);
  ESL_ALLOC(perm,  sizeof(int)   *msa->nseq);
  ESL_ALLOC(useme, sizeof(int)   *msa->nseq);
  ref = msa->ax[idx];
  printf("reference: %s\n", msa->sqname[idx]);
  for (s = 0; s < msa->nseq; s ++) {
    perm[s] = s;
    useme[s] = TRUE;
    esl_dst_XPairId(msa->abc, ref, msa->ax[s], &pid[s], &nid[s], &n[s]);
  }
  esl_vec_DCopy(pid, msa->nseq, cid);
  esl_vec_DSortDecreasing(cid, msa->nseq);
 
  for (s1 = 0; s1 < msa->nseq; s1 ++) {
    for (s2 = 0; s2 < msa->nseq; s2 ++) {
      if (pid[s1] == cid[s2] && useme[s2]) { perm[s2] = s1; useme[s2] = FALSE; break; }
    }
  }
  for (s = 0; s < msa->nseq; s ++) {
    if (s != idx) {
      printf(" pid %.2f nid %d n %d | %s\n", 100*cid[s], nid[perm[s]], n[perm[s]], msa->sqname[perm[s]]);
    }
  }
  *ret_pid = pid;
  *ret_nid = nid;
  *ret_n   = n;
  free(cid);
  free(perm);
  free(useme);
  return eslOK;

 ERROR:
  if (pid)  free(pid);
  if (cid)  free(pid);
  if (nid)  free(nid);
  if (n)    free(n);
  if (perm) free(perm);
  if (useme) free(useme);
  return status;
}
