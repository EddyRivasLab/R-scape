/* rnacov -- significan covarying pairs in an RNA structural alignment
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
#include "covariation.h"
#include "covgrammars.h"
#include "ribosum_matrix.h"

#define ALPHOPTS     "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
#define METHODOPTS   "--naive,--phylo,--dca,--akmaev"              
#define COVTYPEOPTS  "--CHI,--CHIa,--CHIp,--CHIs,--GT,--GTa,--GTp,--GTs,--MI,--MIp,--MIa,--MIs,--MIr,--MIrp,--MIra,--MIrs,--MIg,--MIgp,--MIga,--MIga,--OMES,--OMESp,--OMESa,--OMESa"              
#define COVCLASSOPTS "--C16,--C2"                                          
#define NULLOPTS     "--null1,--null1b,--null2,--null2b,--null3,--null4"                                          
#define THRESHOPTS   "--covNBP,--covNBPu,--covNBPf,--covRBP,--covRBPu,--covRBPf,-E"                                          

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
  double           fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double           idthresh;	       /* fractional identity threshold for selecting subset of sequences */
  double           gapthresh;          /* only keep columns with <= gapthresh fraction of gaps in them */

  double           minidthresh;	       /* fractional minimal identity threshold for selecting subset of sequences */

  COVTYPE          covtype;
  COVCLASS         covclass;

  int              nseqthresh;

  int              onemsa;
  int              nmsa;
  char            *msafile;
  char            *filename;
  char            *msaname;
  int             *msamap;

  char            *outmsafile;
  FILE            *outmsafp;
 
  char            *outdir;
  char            *outfile;
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */
   int              infmt;
 
  char            *R2Rversion;
  int              R2Rall;
  char            *R2Rfile;
  FILE            *R2Rfp;

  int              docyk;
  int              cykLmax;
  char            *R2Rcykfile;
  FILE            *R2Rcykfp;
  int              minloop;
  enum grammar_e   grammar;
  
  int              donullcov;
  char            *covhisfile;
  char            *nullcovhisfile;

  char            *cykcovhisfile;
  char            *cyknullcovhisfile;

  char            *nullcovfile;
  char            *dplotfile;
  char            *cykdplotfile;

  int              nshuffle;

  double          *ft;
  double          *fbp;
  double          *fnbp;

  METHOD           method;
  ESL_TREE        *T;
  double           treeavgt;
 
  NULLTYPE         nulltype;

  char                *ribofile;
  struct ribomatrix_s *ribosum;

  char            *gnuplot;
  
  int              nseqmin;
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

  double           bmin;    /* score histograms */
  double           w;
  double           pmass;
  double           mu;
  double           lambda;

  THRESH          *thresh;
  MODE             mode;

  int              window;
  int              slide;

  float            tol;
  int              verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      1 },
  { "--outdir",     eslARG_STRING,       NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "specify a directory for all output files",                                                  1 },
  { "--cyk",          eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "obtain the structure with maximum covariation",                                             1 },
  { "--r2rall",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "make R2R plot all position in the alignment",                                               1 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                1 },
  { "--window",       eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "window size",                                                                               1 },
  { "--slide",        eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "window slide",                                                                              1 },
  { "--onemsa",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "if file has more than one msa, analyze only the first one",                                 1 },
 /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      1 },
  { "-I",             eslARG_REAL,    "0.97",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < <x> id",                                                             1 },
  { "-i",             eslARG_REAL,      NULL,    NULL, "0<=x<1.0",   NULL,    NULL,  NULL,               "require seqs to have >= <x> id",                                                            1 },
  { "--submsa",       eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   1 },
  { "--nseqmin",      eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "minimum number of sequences in the alignment",                                              1 },  
  { "--gapthresh",    eslARG_REAL,     "0.5",    NULL,  "0<=x<=1",   NULL,    NULL,  NULL,               "keep columns with < <x> fraction of gaps",                                                  1 },
  { "--minid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      1 },
  { "--maxid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",                                                      1 },
  /* msa format */
  { "--informat",   eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            1 },
   /* different ways to assess significance */
  { "-E",            eslARG_REAL,      "0.05",   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "Eval: max expected number of covNBPs allowed",                                              1 },
  { "--covNBP",      eslARG_REAL,       NULL,   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov NonBPs:    max total        (covNPB)      allowed",                                     1 },
  { "--covNBPu",     eslARG_REAL,       NULL,   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov NonBPs;    max per_position (covNBP/alen) allowed",                                     1 },
  { "--covNBPf",     eslARG_REAL,       NULL,   NULL,    "0<x<=1",THRESHOPTS, NULL,  NULL,               "cov NonBPs;    max fraction     (covNBP/NBP)  allowed",                                     1 },
  { "--covRBP",      eslARG_REAL,       NULL,   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov RandomBPs: max total        (covRPB)      allowed",                                     0 },
  { "--covRBPu",     eslARG_REAL,       NULL,   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov RandomBPs; max per_position (covRBP/alen) allowed",                                     0 },
  { "--covRBPf",     eslARG_REAL,       NULL,   NULL,   "0<=x<=1",THRESHOPTS, NULL,  NULL,               "cov RandomBPs; max fraction     (covRBP/RBP)  allowed",                                     0 },
  /* null hypothesis */
  { "--nshuffle",      eslARG_INT,       "20",   NULL,      "n>0",   NULL,    NULL,  NULL,               "number of shuffled sequences",                                                              1 },   
  { "--null1",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null1:  shuffle alignment columns",                                                         0 },
  { "--null1b",       eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null1b: shuffle bpaired_columns and nonbpaired_columns independently",                      0 },
  { "--null2",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null2:  shuffle residues within a column",                                                  0 },
  { "--null2b",       eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null2b: shuffle residues within a column that appear to be canonical (i,j) pairs ",                                                  0 },
  { "--null3",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null3:  null1(b)+null2",                                                                    0 },
  { "--null4",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null4: ",                                                                                   0 },
  /* covariation measures */
  { "--CHIa",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "CHI  ACS corrected calculation",                                                            0 },
  { "--CHIp",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "CHI  APS corrected calculation",                                                            0 },
  { "--CHIs",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "CHI  SCA corrected calculation",                                                            0 },
  { "--CHI",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "CHI  calculation",                                                                          0 },
  { "--GTa",          eslARG_NONE,       TRUE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   ACS corrected calculation",                                                            0 },
  { "--GTp",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   APS corrected calculation",                                                            0 },
  { "--GTs",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   SCA corrected calculation",                                                            0 },
  { "--GT",           eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   calculation",                                                                          0 },
  { "--MIa",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   ACS corrected calculation",                                                            0 },
  { "--MIp",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   APS corrected calculation",                                                            0 },
  { "--MIs",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   SCA corrected calculation",                                                            0 },
  { "--MI",           eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   calculation",                                                                          0 },
  { "--MIra",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIr  ACS corrected calculation",                                                            0 },
  { "--MIrp",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIr  APS corrected calculation",                                                            0 },
  { "--MIrs",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIr  SCA corrected calculation",                                                            0 },
  { "--MIr",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIr  calculation",                                                                          0 },
  { "--MIga",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIg  ACS corrected calculation",                                                            0 },
  { "--MIgp",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIg  APS corrected calculation",                                                            0 },
  { "--MIgs",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIg  SCA corrected calculation",                                                            0 },
  { "--MIg",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIg  calculation",                                                                          0 },
  { "--OMESa",        eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "OMES ACS corrected calculation",                                                            0 },
  { "--OMESp",        eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "OMES APS corrected calculation",                                                            0 },
  { "--OMESs",        eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "OMES SCA corrected calculation",                                                            0 },
  { "--OMES",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "OMES calculation",                                                                          0 },
  /* covariation class */
  { "--C16",         eslARG_NONE,      FALSE,    NULL,       NULL,COVCLASSOPTS,NULL,  NULL,              "use 16 covariation classes",                                                                0 },
  { "--C2",          eslARG_NONE,      FALSE,    NULL,       NULL,COVCLASSOPTS,NULL,  NULL,              "use 2 covariation classes",                                                                 0 },
  { "--nseqthresh",  eslARG_INT,         "8",    NULL,      "n>=0",    NULL,   NULL,"--C2--C16",         "use C2 if nseq <= nseqthresh otherwise use C16",                                            0 },   
   /* phylogenetic method */
  { "--naive",        eslARG_NONE,       TRUE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "naive calculations",                                                                        0 },
  { "--phylo",        eslARG_NONE,      FALSE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "phylo calculations",                                                                        0 },
  { "--dca",          eslARG_NONE,      FALSE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "direct coupling analysis (DCA) MI calculations",                                            0 },
  { "--akmaev",       eslARG_NONE,      FALSE,   NULL,       NULL,METHODOPTS, NULL,  NULL,               "akmaev-style MI calculations",                                                              0 },
   /* Control of scoring system - ribosum */
  { "--ribofile",     eslARG_INFILE,    NULL,    NULL,       NULL,   NULL,    NULL,  "--mx",             "read ribosum structure from file <f>",                                                      0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       1 },
  { "--outmsa",       eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "write actual msa used to file <f>,",                                                        1 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            1 },
 /* other options */  
  { "--cykLmax",       eslARG_INT,     "800",    NULL,      "n>0",   NULL,    NULL, NULL,                "max length to do cykcov calculation",                                                       0 },   
  { "--minloop",       eslARG_INT,       "5",    NULL,      "n>0",   NULL,    NULL, NULL,                "minloop in cykcov calculation",                                                             0 },   
  { "--grammar",    eslARG_STRING,     "BGR",    NULL,       NULL,   NULL,"--cyk",  NULL,                "grammar used for cococyk calculation",                                                      0 },   
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,      "42",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       0 },
  { "--pmass",        eslARG_REAL,    "0.0001",  NULL,       NULL,   NULL,    NULL,  NULL,               "pmass for censored histogram of cov scores",                                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msa>";
static char banner[] = "rnacov - significant covarying pairs in RNA alignments";

static int MSA_banner (FILE *fp, char *msaname, MSA_STAT mstat, MSA_STAT omstat, int nbpairs, int onbpairs);
static int original_msa_manipulate(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **msa);
static int rnacov_for_msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **msa);
static int create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int run_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **msa, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, RANKLIST **ret_ranklist);
static int null1_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null1b_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null2_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null2b_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null3_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null4_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null_add2cumranklist(RANKLIST *ranklist, RANKLIST **ocumranklist, int verbose, char *errbuf);

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, struct cfg_s *ret_cfg)
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
      esl_usage(stdout,  cfg.argv[0], usage);
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
  cfg.msafrq = NULL;
  cfg.msaname = NULL;
  cfg.msamap = NULL;

  /* alphabet */
  cfg.abc = esl_alphabet_Create(eslRNA);
  esl_alphabet_SetEquiv(cfg.abc, '=', '-');     /* allow = as a gap character too */
  esl_alphabet_SetEquiv(cfg.abc, '.', '-');     /* allow . as a gap character too */

  cfg.watch = esl_stopwatch_Create(); 
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));
  
  cfg.outdir = NULL;
  if (esl_opt_IsOn(go, "--outdir")) esl_sprintf( &cfg.outdir, "%s", esl_opt_GetString(go, "--outdir"));
 
 esl_FileTail(cfg.msafile, TRUE, &cfg.filename);
  if ( cfg.outdir ) esl_sprintf( &cfg.outheader, "%s/%s", cfg.outdir, cfg.filename);
  
  if (esl_opt_IsOn(go, "--submsa")) { cfg.submsa = esl_opt_GetInteger(go, "--submsa"); esl_sprintf(&cfg.outheader, "%s.select%d", cfg.outheader, cfg.submsa); }
  else                              { cfg.submsa = 0; }
    
  /* other options */
  cfg.nshuffle    = esl_opt_GetInteger(go, "--nshuffle");
  cfg.nseqthresh  = esl_opt_GetInteger(go, "--nseqthresh");
  cfg.fragfrac    = esl_opt_IsOn(go, "-F")?           esl_opt_GetReal   (go, "-F")           : -1.0;
  cfg.idthresh    = esl_opt_IsOn(go, "-I")?           esl_opt_GetReal   (go, "-I")           : -1.0;
  cfg.minidthresh = esl_opt_IsOn(go, "-i")?           esl_opt_GetReal   (go, "-i")           : -1.0;
  cfg.nseqmin     = esl_opt_IsOn(go, "--nseqmin")?    esl_opt_GetInteger(go, "--nseqmin")    : -1;
  cfg.gapthresh   = esl_opt_GetReal   (go, "--gapthresh");
  cfg.tol         = esl_opt_GetReal   (go, "--tol");
  cfg.verbose     = esl_opt_GetBoolean(go, "-v");
  cfg.voutput     = esl_opt_GetBoolean(go, "--voutput");
  cfg.minloop     = esl_opt_GetInteger(go, "--minloop");
  cfg.docyk       = esl_opt_IsOn(go, "--cyk")? TRUE:FALSE;
  cfg.cykLmax     = esl_opt_GetInteger(go, "--cykLmax");
  cfg.window      = esl_opt_IsOn(go, "--window")?    esl_opt_GetInteger(go, "--window") : -1;
  cfg.slide       = esl_opt_IsOn(go, "--slide")?     esl_opt_GetInteger(go, "--slide")  : -1;
  cfg.onemsa      = esl_opt_IsOn(go, "--onemsa")?    esl_opt_GetBoolean(go, "--onemsa") : FALSE;
  
  if ( esl_opt_IsOn(go, "--grammar") ) {
    if      (esl_strcmp(esl_opt_GetString(go, "--grammar"), "G6")  == 0) cfg.grammar = G6;
    else if (esl_strcmp(esl_opt_GetString(go, "--grammar"), "G6S") == 0) cfg.grammar = G6S;
    else if (esl_strcmp(esl_opt_GetString(go, "--grammar"), "BGR") == 0) cfg.grammar = BGR;
    else esl_fatal("Grammar %s has not been implemented", esl_opt_GetString(go, "--grammar"));
  }
  
  esl_sprintf(&cfg.R2Rversion, "R2R-1.0.4");			
  cfg.R2Rall      = esl_opt_GetBoolean(go, "--r2rall");

  if (cfg.minidthresh > cfg. idthresh) esl_fatal("minidthesh has to be smaller than idthresh");

  cfg.nulltype = NullNONE;
  if      (esl_opt_GetBoolean(go, "--null1"))  cfg.nulltype = Null1;
  else if (esl_opt_GetBoolean(go, "--null1b")) cfg.nulltype = Null1b;
  else if (esl_opt_GetBoolean(go, "--null2"))  cfg.nulltype = Null2;
  else if (esl_opt_GetBoolean(go, "--null2b")) cfg.nulltype = Null2b;
  else if (esl_opt_GetBoolean(go, "--null3"))  cfg.nulltype = Null3;
  else if (esl_opt_GetBoolean(go, "--null4"))  cfg.nulltype = Null4;

  ESL_ALLOC(cfg.thresh, sizeof(THRESH));
  if      (esl_opt_IsOn(go, "--covNBP") )  { cfg.thresh->type = covNBP;  cfg.thresh->val = esl_opt_GetReal(go, "--covNBP");  }
  else if (esl_opt_IsOn(go, "--covNBPu"))  { cfg.thresh->type = covNBPu; cfg.thresh->val = esl_opt_GetReal(go, "--covNBPu"); }
  else if (esl_opt_IsOn(go, "--covNBPf"))  { cfg.thresh->type = covNBPf; cfg.thresh->val = esl_opt_GetReal(go, "--covNBPf"); }
  else if (esl_opt_IsOn(go, "--covRBP") )  { cfg.thresh->type = covRBP;  cfg.thresh->val = esl_opt_GetReal(go, "--covRBP");  }
  else if (esl_opt_IsOn(go, "--covRBPu"))  { cfg.thresh->type = covRBPu; cfg.thresh->val = esl_opt_GetReal(go, "--covRBPu"); }
  else if (esl_opt_IsOn(go, "--covRBPf"))  { cfg.thresh->type = covRBPf; cfg.thresh->val = esl_opt_GetReal(go, "--covRBPf"); }
  else if (esl_opt_IsOn(go, "-E"))         { cfg.thresh->type = Eval;    cfg.thresh->val = esl_opt_GetReal(go, "-E"); if (cfg.nulltype == NullNONE) cfg.nulltype = Null4; }

  if      (esl_opt_GetBoolean(go, "--CHIa"))  cfg.covtype = CHIa;
  else if (esl_opt_GetBoolean(go, "--CHIp"))  cfg.covtype = CHIp;
  else if (esl_opt_GetBoolean(go, "--CHIs"))  cfg.covtype = CHIs;
  else if (esl_opt_GetBoolean(go, "--CHI"))   cfg.covtype = CHI;
  else if (esl_opt_GetBoolean(go, "--GTa"))   cfg.covtype = GTa;
  else if (esl_opt_GetBoolean(go, "--GTp"))   cfg.covtype = GTp;
  else if (esl_opt_GetBoolean(go, "--GTs"))   cfg.covtype = GTs;
  else if (esl_opt_GetBoolean(go, "--GT"))    cfg.covtype = GT;
  else if (esl_opt_GetBoolean(go, "--MIa"))   cfg.covtype = MIa;
  else if (esl_opt_GetBoolean(go, "--MIp"))   cfg.covtype = MIp;
  else if (esl_opt_GetBoolean(go, "--MIs"))   cfg.covtype = MIs;
  else if (esl_opt_GetBoolean(go, "--MI"))    cfg.covtype = MI;
  else if (esl_opt_GetBoolean(go, "--MIra"))  cfg.covtype = MIra;
  else if (esl_opt_GetBoolean(go, "--MIrp"))  cfg.covtype = MIrp;
  else if (esl_opt_GetBoolean(go, "--MIrs"))  cfg.covtype = MIrs;
  else if (esl_opt_GetBoolean(go, "--MIr"))   cfg.covtype = MIr;
  else if (esl_opt_GetBoolean(go, "--MIga"))  cfg.covtype = MIga;
  else if (esl_opt_GetBoolean(go, "--MIgp"))  cfg.covtype = MIgp;
  else if (esl_opt_GetBoolean(go, "--MIgs"))  cfg.covtype = MIgs;
  else if (esl_opt_GetBoolean(go, "--MIg"))   cfg.covtype = MIg;
  else if (esl_opt_GetBoolean(go, "--OMESa")) cfg.covtype = OMESa;
  else if (esl_opt_GetBoolean(go, "--OMESp")) cfg.covtype = OMESp;
  else if (esl_opt_GetBoolean(go, "--OMESs")) cfg.covtype = OMESs;
  else if (esl_opt_GetBoolean(go, "--OMES"))  cfg.covtype = OMES;

 
  if      (esl_opt_GetBoolean(go, "--C16"))   cfg.covclass = C16;
  else if (esl_opt_GetBoolean(go, "--C2"))    cfg.covclass = C2;
  else                                        cfg.covclass = CSELECT;

  if      (esl_opt_GetBoolean(go, "--naive"))  cfg.method = OPTNONE;
  else if (esl_opt_GetBoolean(go, "--phylo"))  cfg.method = PHYLO;
  else if (esl_opt_GetBoolean(go, "--dca"))    cfg.method = DCA;
  else if (esl_opt_GetBoolean(go, "--akmaev")) cfg.method = AKMAEV;
 
  cfg.bmin  = -300.0; /* a guess for lowest cov score */
  cfg.w     = W;      /* histogram step */
  cfg.pmass = esl_opt_GetReal   (go, "--pmass");

  /* output file */
  if ( esl_opt_IsOn(go, "-o") ) {
    esl_sprintf(&cfg.outfile, "%s", esl_opt_GetString(go, "-o"));
    if ((cfg.outfp = fopen(cfg.outfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outfile);
  } 
  else {
    esl_sprintf(&cfg.outfile, "%s.out", cfg.outheader);
    if ((cfg.outfp = fopen(cfg.outfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outfile);
  }
 
  /*  rocplot file */
  esl_sprintf(&cfg.rocfile, "%s.roc", cfg.outheader); 
  if ((cfg.rocfp = fopen(cfg.rocfile, "w")) == NULL) esl_fatal("Failed to open rocfile %s", cfg.rocfile);

  /*  summary file */
  esl_sprintf(&cfg.sumfile, "%s.sum", cfg.outheader); 
  if ((cfg.sumfp = fopen(cfg.sumfile, "w")) == NULL) esl_fatal("Failed to open sumfile %s", cfg.sumfile);
  
  /* msa-specific files */
  cfg.outmsafile = NULL;
  cfg.outmsafp   = NULL;

  cfg.R2Rfile    = NULL;
  cfg.R2Rfp      = NULL;
  cfg.R2Rcykfile = NULL;
  cfg.R2Rcykfp   = NULL;

  /* covhis file */
  cfg.covhisfile    = NULL;
  cfg.cykcovhisfile = NULL;
  
  /* nullcovhis file */
  cfg.donullcov         = FALSE;
  cfg.nullcovhisfile    = NULL;
  cfg.cyknullcovhisfile = NULL;
  
  /* nullcovplot file */
  cfg.nullcovfile = NULL;

  /* dotplot file */
  cfg.dplotfile    = NULL;
  cfg.cykdplotfile = NULL;
  
  cfg.shsumfile = NULL;
  cfg.shsumfp   = NULL;
  
  if (cfg.nulltype != NullNONE) {
    /*  sh-summary file */
    esl_sprintf(&cfg.shsumfile, "%s.shsum", cfg.outheader); 
    if ((cfg.shsumfp = fopen(cfg.shsumfile, "w")) == NULL) esl_fatal("Failed to open shsumfile %s", cfg.shsumfile);
  }
  
  cfg.T  = NULL;
  cfg.ct = NULL;
  cfg.nbpairs = 0;

  cfg.ft   = NULL;
  cfg.fbp  = NULL;
  cfg.fnbp = NULL;

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

  *ret_go  = go;
  *ret_cfg = cfg;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, banner, usage);
  if (puts("\nwhere options are:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
MSA_banner (FILE *fp, char *msaname, MSA_STAT mstat, MSA_STAT omstat, int nbpairs, int onbpairs)
{
  fprintf(fp, "# MSA %s nseq %d (%d) alen %" PRId64 " (%" PRId64 ") avgid %.2f (%.2f) nbpairs %d (%d)\n", 
	  msaname, mstat.nseq, omstat.nseq, mstat.alen, omstat.alen, 
	  mstat.avgid, omstat.avgid, nbpairs, onbpairs);
  return eslOK;
}

int
main(int argc, char **argv)
{ 
  ESL_GETOPTS     *go;
  struct cfg_s     cfg;
  char            *omsaname = NULL;
  ESLX_MSAFILE    *afp = NULL;
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
  status = eslx_msafile_Open(NULL, cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);
  eslx_msafile_SetDigital(afp, cfg.abc);

  /* read the MSA */
  while ((hstatus = eslx_msafile_Read(afp, &msa)) != eslEOF) {
    if (hstatus != eslOK) eslx_msafile_ReadFailure(afp, status);
    cfg.nmsa ++;
    if (cfg.onemsa && cfg.nmsa > 1) break;

    status = original_msa_manipulate(go, &cfg, &msa);
    if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to manipulate alignment"); }
    if (msa == NULL) continue;
    
    if (cfg.window > 0) {
      esl_sprintf(&omsaname, "%s", cfg.msaname);

      useme = malloc(sizeof(int)*(msa->alen+1));
      for (first = 1; first <= msa->alen-cfg.slide; first += cfg.slide) {
	esl_vec_ISet(useme, msa->alen+1, FALSE);
	wmsa = esl_msa_Clone(msa);

	last = ESL_MIN(first+cfg.window-1, msa->alen);
	esl_sprintf(&cfg.msaname, "%s_%d-%d", omsaname, first, last);

	for (i = first; i <= last; i ++) useme[i] = TRUE;
	status = esl_msa_ColumnSubset(wmsa, cfg.errbuf, useme);

	status = rnacov_for_msa(go, &cfg, &wmsa);
	if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to run rnacov"); }
 
	esl_msa_Destroy(wmsa); wmsa = NULL;
      }
    }
    else {
      status = rnacov_for_msa(go, &cfg, &msa);
      if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to run rnacov"); }
    }

    if (omsaname) free(omsaname); omsaname = NULL;
    if (useme) free(useme); useme = NULL;
    if (msa) esl_msa_Destroy(msa); msa = NULL;
    if (cfg.msaname) free(cfg.msaname); cfg.msaname = NULL;
  }

  /* cleanup */
  esl_stopwatch_Destroy(cfg.watch);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  if (cfg.ct) free(cfg.ct);
  eslx_msafile_Close(afp);
  if (cfg.msaname) free(cfg.msaname);
  if (cfg.outfile) free(cfg.outfile);
  if (cfg.outdir) free(cfg.outdir);
  free(cfg.outheader);
  fclose(cfg.outfp);
  fclose(cfg.rocfp);
  fclose(cfg.sumfp);
  free(cfg.gnuplot);
  if (cfg.shsumfp) fclose(cfg.shsumfp);
  if (cfg.ribosum) Ribosum_matrix_Destroy(cfg.ribosum);
  if (cfg.outmsafile) free(cfg.outmsafile);
  if (cfg.outmsafp) fclose(cfg.outmsafp);
  if (cfg.R2Rfile) free(cfg.R2Rfile); 
  if (cfg.R2Rversion) free(cfg.R2Rversion); 
  if (cfg.R2Rfp) fclose(cfg.R2Rfp); 
  if (cfg.covhisfile) free(cfg.covhisfile); 
  if (cfg.nullcovhisfile) free(cfg.nullcovhisfile);
  if (cfg.cykcovhisfile) free(cfg.cykcovhisfile);
  if (cfg.cyknullcovhisfile) free(cfg.cyknullcovhisfile);
  if (cfg.nullcovfile) free(cfg.nullcovfile);
  if (cfg.dplotfile) free(cfg.dplotfile);
  if (cfg.cykdplotfile) free(cfg.cykdplotfile);
  if (cfg.ft) free(cfg.ft);
  if (cfg.fbp) free(cfg.fbp);
  if (cfg.fnbp) free(cfg.fnbp);
  if (cfg.thresh) free(cfg.thresh);
  if (useme) free(useme);
  return 0;
}


static int
original_msa_manipulate(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa)
{
  ESL_MSA *msa = *omsa;
  char    *msg = "original_msa_manipulate failed";
  char    *type = NULL;
  char    *tp;
  char    *tok;
  int      seq_cons_len = 0;
  int      nremoved = 0;	  /* # of identical sequences removed */
  int      nfrags = 0;	  /* # of fragments removed */
  int      t;

  /* stats of the original alignment */
  msamanip_XStats(msa, &cfg->omstat);
  msamanip_CalculateCT(msa, NULL, &cfg->onbpairs, cfg->errbuf);
  
  /* the msaname */
  for (t = 0; t < msa->ngf; t++) {
    if (!esl_strcmp(msa->gf_tag[t], "TP")) {
      tp = msa->gf[t];	
      while (*tp != '\0') {
	if (esl_strtok(&tp,  ";", &tok) != eslOK) esl_fatal(msg);
	if (esl_strtok(&tok, " ", &tok) != eslOK) esl_fatal(msg);
	esl_sprintf(&tok, "_%s", tok);
	esl_strcat(&type, -1, tok, -1);
      }
    }
  }
  if      (msa->acc  && msa->name && type) esl_sprintf(&cfg->msaname, "%s_%s%s", msa->acc, msa->name, type);
  else if (msa->acc  && msa->name)         esl_sprintf(&cfg->msaname, "%s_%s", msa->acc, msa->name);
  else if (msa->name && type)              esl_sprintf(&cfg->msaname, "%s%s", msa->name, type);
  else if (msa->acc)                       esl_sprintf(&cfg->msaname, "%s", msa->acc);
  else if (msa->name)                      esl_sprintf(&cfg->msaname, "%s", msa->name);
  else                                     esl_sprintf(&cfg->msaname, "%s_%d", cfg->filename, cfg->nmsa);
  if (esl_opt_IsOn(go, "--submsa")) esl_sprintf(&cfg->msaname, "%s.select%d", cfg->msaname, esl_opt_GetInteger(go, "--submsa"));
  
  /* apply msa filters and than select submsa
   */
  if (cfg->fragfrac > 0.     && msamanip_RemoveFragments(cfg->fragfrac, &msa, &nfrags, &seq_cons_len)             != eslOK) { printf("remove_fragments failed\n");     esl_fatal(msg); }
  if (esl_opt_IsOn(go, "-I") && msamanip_SelectSubsetBymaxID(cfg->r, &msa, cfg->idthresh, &nremoved)              != eslOK) { printf("select_subsetBymaxID failed\n"); esl_fatal(msg); }
  if (esl_opt_IsOn(go, "-i") && msamanip_SelectSubsetByminID(cfg->r, &msa, cfg->minidthresh, &nremoved)           != eslOK) { printf("select_subsetByminID failed\n"); esl_fatal(msg); }
  if (cfg->submsa            && msamanip_SelectSubset(cfg->r, cfg->submsa, &msa, NULL, cfg->errbuf, cfg->verbose) != eslOK) { printf("%s\n", cfg->errbuf);              esl_fatal(msg); }
  if (msa == NULL) {
    free(type); type = NULL;
    free(cfg->msaname); cfg->msaname = NULL;
    return eslOK;
  }
  if (msa->nseq < cfg->nseqmin) {
    esl_msa_Destroy(msa); msa = NULL;
    free(type); type = NULL;
    free(cfg->msaname); cfg->msaname = NULL;
    return eslOK;
  }
  if (msamanip_RemoveGapColumns(cfg->gapthresh, msa, &cfg->msamap, cfg->errbuf, cfg->verbose) != eslOK) { printf("RemoveGapColumns\n"); esl_fatal(msg); }
  
  esl_msa_Hash(msa);
  msamanip_ConvertDegen2RandomCanonical(cfg->r, msa);
  if (esl_msa_MinimGaps(msa, NULL, "-.~=", FALSE) != eslOK) esl_fatal("Failed to remove minim gaps");
  
  /* given msa aveid and avematch */
  msamanip_XStats(msa, &cfg->mstat);
  
  if ((esl_opt_IsOn(go, "--minid") && cfg->mstat.avgid < 100.*esl_opt_GetReal(go, "--minid")) ||
      (esl_opt_IsOn(go, "--maxid") && cfg->mstat.avgid > 100.*esl_opt_GetReal(go, "--maxid"))   ) {
    esl_msa_Destroy(msa); msa = NULL;
    free(type); type = NULL;
    free(cfg->msaname); cfg->msaname = NULL;
    return eslOK;  
  }
  
  /* print some info */
  if (cfg->verbose) {
    fprintf(cfg->outfp, "Used alignment\n");
    fprintf(cfg->outfp, "%6d          %s\n", msa->nseq, cfg->msafile);
    if (eslx_msafile_Write(cfg->outfp, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write msa"); 
    msamanip_DumpStats(cfg->outfp, msa, cfg->mstat); 
  }
  
  *omsa = msa;
  return eslOK;
}

static int
rnacov_for_msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa)
{
  ESL_MSA  *msa = *omsa;
  RANKLIST *ranklist_null = NULL;
  RANKLIST *ranklist_aux  = NULL;
  int       status;
  
  if (msa->nseq <= 1) {
    MSA_banner(cfg->outfp, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
    return eslOK;
  }

  /* outmsa file if requested */
  if (esl_opt_IsOn(go, "--outmsa")) {
    esl_sprintf(&cfg->outmsafile, "%s/%s_%s.sto",     cfg->outdir, esl_opt_GetString(go, "--outmsa"), cfg->msaname);
    if ((cfg->outmsafp = fopen(cfg->outmsafile, "w")) == NULL) esl_fatal("Failed to open outmsa file %s", cfg->outmsafile);
    eslx_msafile_Write(cfg->outmsafp, msa, eslMSAFILE_STOCKHOLM);
    fclose(cfg->outmsafp);
  } 

  /* R2R annotated sto file */
  if (cfg->outdir) {
    esl_sprintf(&cfg->R2Rfile,    "%s/%s.R2R.sto",     cfg->outdir, cfg->msaname);
    esl_sprintf(&cfg->R2Rcykfile, "%s/%s.cyk.R2R.sto", cfg->outdir, cfg->msaname);
    
    /* covhis file */
    esl_sprintf(&cfg->covhisfile,    "%s/%s.his",     cfg->outdir, cfg->msaname);
    esl_sprintf(&cfg->cykcovhisfile, "%s/%s.cyk.his", cfg->outdir, cfg->msaname);
    
    /* nullcovhis file */
    if (cfg->donullcov) {
      esl_sprintf(&cfg->nullcovhisfile,    "%s/%s.nullhis",     cfg->outdir, cfg->msaname);
      esl_sprintf(&cfg->cyknullcovhisfile, "%s/%s.cyk.nullhis", cfg->outdir, cfg->msaname);
      
      /* nullcovplot file */
      esl_sprintf(&cfg->nullcovfile, "%s/%s.nullcov", cfg->outdir, cfg->msaname);
    }
    
    /* dotplot file */
    esl_sprintf(&cfg->dplotfile,    "%s/%s.dplot",     cfg->outdir, cfg->msaname);
    esl_sprintf(&cfg->cykdplotfile, "%s/%s.cyk.dplot", cfg->outdir, cfg->msaname);
  }
  else {
    esl_sprintf(&cfg->R2Rfile,    "%s.R2R.sto",     cfg->msaname);
    esl_sprintf(&cfg->R2Rcykfile, "%s.cyk.R2R.sto", cfg->msaname);
    
    /* covhis file */
    esl_sprintf(&cfg->covhisfile,    "%s.his",     cfg->msaname);
    esl_sprintf(&cfg->cykcovhisfile, "%s.cyk.his", cfg->msaname);
    
    /* nullcovhis file */
    if (cfg->donullcov) {
      esl_sprintf(&cfg->nullcovhisfile,    "%s.nullhis",     cfg->msaname);
      esl_sprintf(&cfg->cyknullcovhisfile, "%s.cyk.nullhis", cfg->msaname);
    
      /* nullcovplot file */
      esl_sprintf(&cfg->nullcovfile, "%s.nullcov", cfg->msaname);
    }
    
    /* dotplot file */
    esl_sprintf(&cfg->dplotfile,    "%s.dplot",     cfg->msaname);
    esl_sprintf(&cfg->cykdplotfile, "%s.cyk.dplot", cfg->msaname);
  }

  /* the ct vector */
  status = msamanip_CalculateCT(msa, &cfg->ct, &cfg->nbpairs, cfg->errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s. Failed to calculate ct vector", cfg->errbuf);
#if 0
  msamanip_CalculateBC(msa, cfg->ct, &cfg->ft, &cfg->fbp, &cfg->fnbp, cfg->errbuf);
  esl_vec_DDump(stdout, cfg->ft,   cfg->abc->K, "total BC");
  esl_vec_DDump(stdout, cfg->fbp,  cfg->abc->K, "basepairs BC");
  esl_vec_DDump(stdout, cfg->fnbp, cfg->abc->K, "nonbasepairs BC");
#endif
  if (cfg->nbpairs == 0)        cfg->docyk = TRUE;  // calculate the cyk-cov structure if no given one
  if (msa->alen > cfg->cykLmax) cfg->docyk = FALSE; // unless alignment is too long
  
  if (1||cfg->verbose) 
    MSA_banner(stdout, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
  
  /* the null model first */
  cfg->mode = RANSS;
  if (cfg->nulltype == Null1) {
    status = null1_rnacov(go, cfg, msa, &ranklist_null);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null1 rnacov", cfg->errbuf);
  }
  else if (cfg->nulltype == Null1b) {
    status = null1b_rnacov(go, cfg, msa, &ranklist_null);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null1b rnacov", cfg->errbuf);
  }
  else if (cfg->nulltype == Null2) {
    status = null2_rnacov(go, cfg, msa, &ranklist_null);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null2 rnacov", cfg->errbuf);
  }
  else if (cfg->nulltype == Null2b) {
    cfg->nulltype = Null2;
    status = null2_rnacov(go, cfg, msa, &ranklist_aux);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null2 rnacov", cfg->errbuf);
    cfg->nulltype = Null2b;
    status = null2b_rnacov(go, cfg, msa, &ranklist_null);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null2b rnacov", cfg->errbuf);
  }
  else if (cfg->nulltype == Null3) {
    status = null3_rnacov(go, cfg, msa, &ranklist_null);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null3 rnacov", cfg->errbuf);
  }
  else if (cfg->nulltype == Null4) {
#if 0
    cfg->nulltype = Null2;
    status = null2_rnacov(go, cfg, msa, &ranklist_aux);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null2 rnacov", cfg->errbuf);
#endif
    cfg->nulltype = Null4;
    status = null4_rnacov(go, cfg, msa, &ranklist_null);
    if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run null4 rnacov", cfg->errbuf);
  }
  
  /* main function */
  cfg->mode = GIVSS;
  status = run_rnacov(go, cfg, &msa, ranklist_null, ranklist_aux, NULL);
  if (status != eslOK) ESL_XFAIL(status, cfg->errbuf, "%s.\nFailed to run rnacov", cfg->errbuf);

  *omsa = msa;

  free(cfg->ct); cfg->ct = NULL;
  if (cfg->msafrq) free(cfg->msafrq); cfg->msafrq = NULL;
  if (cfg->T) esl_tree_Destroy(cfg->T); cfg->T = NULL;
  cov_FreeRankList(ranklist_null); ranklist_null = NULL;
  if (ranklist_aux) cov_FreeRankList(ranklist_aux); ranklist_aux = NULL;

  return eslOK;

 ERROR:
  if (cfg->ct) free(cfg->ct);
  if (cfg->msafrq) free(cfg->msafrq); 
  if (cfg->T) esl_tree_Destroy(cfg->T); 
  if (cfg->msaname) free(cfg->msaname);
  if (ranklist_null) cov_FreeRankList(ranklist_null); 
  if (ranklist_aux) cov_FreeRankList(ranklist_aux);
  return status;
}

static int
create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  int  status;
  
  /* the TREE */
  status = Tree_CalculateExtFromMSA(msa, &cfg->T, TRUE, cfg->errbuf, cfg->verbose);
  if (status != eslOK) { printf("%s\n", cfg->errbuf); esl_fatal(cfg->errbuf); }
  if (cfg->T) {
    cfg->treeavgt = esl_tree_er_AverageBL(cfg->T); 
    if (cfg->verbose) Tree_Dump(stdout, cfg->T, "Tree");
  }
  
  return eslOK;
}


static int
run_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa, RANKLIST *ranklist_null, RANKLIST *ranklist_aux, RANKLIST **ret_ranklist)
{
  char            *title = NULL;
  struct mutual_s *mi   = NULL;
  ESL_MSA         *msa = *omsa;
  RANKLIST        *ranklist = NULL;
  RANKLIST        *cykranklist = NULL;
  HITLIST         *hitlist = NULL;
  int              donull2b;
  int              nnodes;
  int              status;

  esl_stopwatch_Start(cfg->watch);
  
  /* print to stdout */
  if (cfg->verbose) 
    MSA_banner(stdout, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
  
  if (cfg->mode != RANSS) {
    MSA_banner(cfg->outfp, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
    esl_sprintf(&title, "%s (seqs %d alen %" PRId64 " avgid %d bpairs %d)", 
		cfg->msaname, msa->nseq, msa->alen, (int)ceil(cfg->mstat.avgid), cfg->nbpairs);
  }

  /* produce a tree
   */
  if (cfg->method != NAIVE) {
    status = create_tree(go, cfg, msa);
    if (status != eslOK)  { esl_fatal(cfg->errbuf); }
    nnodes = (cfg->T->N > 1)? cfg->T->N-1 : cfg->T->N;
  }

  /* create the MI structure */
  mi = cov_Create(msa->alen, msa->nseq, (cfg->mode == RANSS)?TRUE:FALSE, cfg->nseqthresh, cfg->abc);
  
  /* write MSA info to the sumfile */
  if (cfg->mode != RANSS) 
    fprintf(cfg->sumfp, "%f\t%s\t%d\t%d\t%.2f\t", cfg->thresh->val, cfg->msaname, msa->nseq, (int)msa->alen, cfg->mstat.avgid); 
  else
    fprintf(cfg->shsumfp, "%f\t%s\t%d\t%d\t%.2f\t", cfg->thresh->val, cfg->msaname, msa->nseq, (int)msa->alen, cfg->mstat.avgid); 
  
  /* write MSA info to the rocfile */
  if (cfg->mode == GIVSS || cfg->mode == CYKSS)
  fprintf(cfg->rocfp, "# MSA nseq %d alen %" PRId64 " avgid %f nbpairs %d (%d)\n", msa->nseq, msa->alen, cfg->mstat.avgid, cfg->nbpairs, cfg->onbpairs);  
  
  /* main function */
  donull2b = (cfg->mode == RANSS && cfg->nulltype == Null2b)? TRUE:FALSE;
  status = cov_Calculate(cfg->r, &msa, cfg->msamap, cfg->T, cfg->ribosum, mi, ranklist_null, ranklist_aux, &ranklist, &hitlist, cfg->method, 
			 cfg->covtype, cfg->covclass, cfg->ct, 
			 cfg->bmin, cfg->w, cfg->pmass, &cfg->mu, &cfg->lambda,
			 (cfg->mode == RANSS)?NULL:cfg->outfp, cfg->rocfp, 
			 (cfg->mode == RANSS)?cfg->shsumfp:cfg->sumfp, cfg->gnuplot, cfg->dplotfile, cfg->R2Rfile, cfg->R2Rversion, cfg->R2Rall, 
			 cfg->thresh, cfg->mode, cfg->onbpairs, donull2b, cfg->tol, cfg->verbose, cfg->errbuf);   
  if (status != eslOK) goto ERROR; 
  if (cfg->mode == GIVSS && (cfg->verbose)) cov_DumpRankList(stdout, ranklist);
    
  if (cfg->mode == GIVSS) {
    if (cfg->verbose) {
      printf("score total distribution\n");
      printf("imin %d imax %d xmax %f xmin %f\n", ranklist->ha->imin, ranklist->ha->imax, ranklist->ha->xmax, ranklist->ha->xmin);
      printf("score truncated distribution\n");
      printf("imin %d imax %d xmax %f xmin %f\n", ranklist->ht->imin, ranklist->ht->imax, ranklist->ht->xmax, ranklist->ht->xmin);
    }
    status = cov_WriteHistogram(cfg->gnuplot, cfg->covhisfile, cfg->nullcovhisfile, ranklist, ranklist_null, ranklist_aux, title, 
				cfg->pmass, cfg->mu, cfg->lambda, cfg->verbose, cfg->errbuf);
    if (status != eslOK) goto ERROR; 
  }
  
  if (cfg->donullcov) {
    status = cov_CreateNullCov(cfg->gnuplot, cfg->nullcovfile, msa->alen, cfg->ct, ranklist, ranklist_null, FALSE, cfg->errbuf);
    if (status != eslOK) goto ERROR; 
  }

  /* find the cykcov structure, and do the cov analysis on it */
  if (cfg->docyk && cfg->mode != RANSS) {
    status = cov_CYKCOVCT(cfg->outfp, cfg->gnuplot, cfg->cykdplotfile, cfg->R2Rcykfile, cfg->R2Rversion, cfg->R2Rall, cfg->r, 
			  &msa, mi, cfg->msamap, ranklist_null, ranklist_aux, &cykranklist, cfg->bmin, cfg->w,
			  cfg->pmass, &cfg->mu, &cfg->lambda, cfg->minloop, cfg->grammar, 
			  cfg->thresh, cfg->thresh->sc, cfg->onbpairs, cfg->errbuf, cfg->verbose);
    if (status != eslOK) goto ERROR;
    
    if (cfg->verbose) {
      printf("score cyk truncated distribution\n");
      printf("imin %d imax %d xmax %f xmin %f\n", cykranklist->ht->imin, cykranklist->ht->imax, cykranklist->ht->xmax, cykranklist->ht->xmin);
      //esl_histogram_Plot(stdout, ranklist->ht);
    }
    status = cov_WriteHistogram(cfg->gnuplot, cfg->cykcovhisfile, cfg->cyknullcovhisfile, cykranklist, ranklist_null, ranklist_aux, title, cfg->pmass, 
				cfg->mu, cfg->lambda, cfg->verbose, cfg->errbuf);
    if (status != eslOK) goto ERROR; 
  }
 
  *omsa = msa;
  if (ret_ranklist) *ret_ranklist = ranklist; else cov_FreeRankList(ranklist);
  if (cykranklist) cov_FreeRankList(cykranklist);
  if (hitlist) cov_FreeHitList(hitlist);
  if (title) free(title);
  cov_Destroy(mi); mi = NULL;
  
  return eslOK;
  
 ERROR:
  *omsa = NULL;
  if (mi)          cov_Destroy(mi);
  if (ranklist)    cov_FreeRankList(ranklist);
  if (cykranklist) cov_FreeRankList(cykranklist);
  if (hitlist)     cov_FreeHitList(hitlist);
  if (title)       free(title);
  return status;
}

/* shuffle all the alignment columns, leave the ss intact */
static int
null1_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *shmsa = NULL;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int       *useme = NULL;
  int        n;                   // index for alignment positions
  int        s;                   // index for shuffles
  int        b;                   // index for ranked list
  int        status;

  /* shuffle all the columns */
  ESL_ALLOC(useme, sizeof(int) * msa->alen);
  for (n = 0; n < msa->alen; n ++) useme[n] = TRUE;

  for (s = 0; s < cfg->nshuffle; s ++) {
    msamanip_ShuffleColumns(cfg->r, msa, &shmsa, useme, cfg->errbuf, cfg->verbose);
    status = run_rnacov(go, cfg, &shmsa, NULL, NULL, &ranklist);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s. Failed to run null1_rnacov", cfg->errbuf);

     status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
     if (status != eslOK) goto ERROR;
     
    esl_msa_Destroy(shmsa); shmsa = NULL;
    cov_FreeRankList(ranklist); ranklist = NULL;
  }
  
  for (b = cumranklist->ha->imin; b <= cumranklist->ha->imax; b ++) {
    cumranklist->covBP[b]  /= (double)cfg->nshuffle;
    cumranklist->covNBP[b] /= (double)cfg->nshuffle;
  }
  
   if (cfg->verbose) {
     printf("null1 distribution - cummulative\n");
     printf("imin %d imax %d xmax %f xmin %f\n", 
	    cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
     //esl_histogram_Plot(stdout, cumranklist->h);
     //esl_histogram_PlotSurvival(stdout, cumranklist->h);
   }

  if (cfg->verbose) cov_DumpRankList(stdout, cumranklist);

  *ret_cumranklist = cumranklist;
  free(useme);
  return eslOK;
  
 ERROR:
  if (useme) free(useme);
  if (shmsa) esl_msa_Destroy(shmsa);
  if (ranklist) cov_FreeRankList(ranklist);
  if (cumranklist) cov_FreeRankList(cumranklist);
  return status;
}

/* shuffle the paired alignment columns and the unpaired columns independently, leave the ss intact */
static int
null1b_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *shmsa = NULL;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int       *useme1 = NULL;
  int       *useme2 = NULL;
  int        n;                   // index for alignment positions
  int        s;                   // index for shuffles
  int        b;                   // index for ranked list
  int        status;

  /* shuffle all the columns */
  ESL_ALLOC(useme1, sizeof(int) * msa->alen);
  ESL_ALLOC(useme2, sizeof(int) * msa->alen);
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] >  0) useme1[n] = TRUE; else useme1[n] = FALSE;
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] == 0) useme2[n] = TRUE; else useme2[n] = FALSE;

  for (s = 0; s < cfg->nshuffle; s ++) {
    msamanip_ShuffleColumns(cfg->r, msa,   &shmsa, useme1, cfg->errbuf, cfg->verbose);
    msamanip_ShuffleColumns(cfg->r, shmsa, &shmsa, useme2, cfg->errbuf, cfg->verbose);
    status = run_rnacov(go, cfg, &shmsa, NULL, NULL, &ranklist);
    if (status != eslOK) ESL_XFAIL(eslFAIL, "%s.\nFailed to run null1b rnacov", cfg->errbuf);
   
     status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
     if (status != eslOK) goto ERROR;
     
     esl_msa_Destroy(shmsa); shmsa = NULL;
     cov_FreeRankList(ranklist); ranklist = NULL;
  }
  
  for (b = cumranklist->ha->imin; b <= cumranklist->ha->imax; b ++) {
    cumranklist->covBP[b]  /= (double)cfg->nshuffle;
    cumranklist->covNBP[b] /= (double)cfg->nshuffle;
  }
  
  if (cfg->verbose) {
    printf("null1b distribution - cummulative\n");
    printf("imin %d imax %d xmax %f xmin %f\n", 
	   cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
    //esl_histogram_Plot(stdout, cumranklist->h);
    //esl_histogram_PlotSurvival(stdout, cumranklist->h);
  }

  *ret_cumranklist = cumranklist;
  free(useme1);
  free(useme2);
  return eslOK;
  
 ERROR:
  if (useme1) free(useme1);
  if (useme2) free(useme2);
  if (shmsa) esl_msa_Destroy(shmsa);
  if (ranklist) cov_FreeRankList(ranklist);
  if (cumranklist) cov_FreeRankList(cumranklist);
  return status;
}


/* shuffle residues within each column */
static int
null2_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *shmsa = NULL;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int       s;
  int       b;
  int       status;

   for (s = 0; s < cfg->nshuffle; s ++) {
     msamanip_ShuffleWithinColumn(cfg->r, msa, &shmsa, cfg->errbuf, cfg->verbose);
     status = run_rnacov(go, cfg, &shmsa, NULL, NULL, &ranklist);
     if (status != eslOK) ESL_XFAIL(eslFAIL, "%s.\nFailed to run null2 rnacov", cfg->errbuf);
     
     status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
     if (status != eslOK) goto ERROR;

     esl_msa_Destroy(shmsa); shmsa = NULL;
     cov_FreeRankList(ranklist); ranklist = NULL;
   }
   
   for (b = cumranklist->ha->imin; b <= cumranklist->ha->imax; b ++) {
     cumranklist->covBP[b]  /= (double)cfg->nshuffle;
     cumranklist->covNBP[b] /= (double)cfg->nshuffle;
   }
   if (cfg->verbose) {
     printf("null2 distribution - cummulative\n");
     printf("imin %d imax %d xmax %f xmin %f\n", 
	    cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
     //esl_histogram_Plot(stdout, cumranklist->h);
     //esl_histogram_PlotSurvival(stdout, cumranklist->h);
   }
   
   if (cfg->verbose) cov_DumpRankList(stdout, cumranklist);
   
   *ret_cumranklist = cumranklist;
   return eslOK;

 ERROR:
   if (shmsa) esl_msa_Destroy(shmsa);
   if (ranklist) cov_FreeRankList(ranklist);
   if (cumranklist) cov_FreeRankList(cumranklist);
   return status;
}

/* shuffle residues within each column in a pair-dependent manner:
 * for two columns i,j shuffle the residues of column i (j) that have
 * a partner in column j (i) forming a canonical pair.
 */
static int
null2b_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *shmsa = NULL;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int       s;
  int       b;
  int       status;

  /* this shuffle is comparisond dependent, cannot create a general shmsa,
   * shuffling needs to be done on the spot.
   * We make a copy of the original msa anyway.
   */
   for (s = 0; s < cfg->nshuffle; s ++) {
     shmsa = esl_msa_Clone(msa);
     status = run_rnacov(go, cfg, &shmsa, NULL, NULL, &ranklist);
     if (status != eslOK) ESL_XFAIL(eslFAIL, "%s.\nFailed to run null2b rnacov", cfg->errbuf);

     status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
     if (status != eslOK) goto ERROR;

     esl_msa_Destroy(shmsa); shmsa = NULL;
     cov_FreeRankList(ranklist); ranklist = NULL;
   }
   
   for (b = cumranklist->ha->imin; b <= cumranklist->ha->imax; b ++) {
     cumranklist->covBP[b]  /= (double)cfg->nshuffle;
     cumranklist->covNBP[b] /= (double)cfg->nshuffle;
   }
   if (cfg->verbose) {
     printf("null2b distribution - cummulative\n");
     printf("imin %d imax %d xmax %f xmin %f\n", 
	    cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
     //esl_histogram_Plot(stdout, cumranklist->h);
     //esl_histogram_PlotSurvival(stdout, cumranklist->h);
   }
   
   if (cfg->verbose) cov_DumpRankList(stdout, cumranklist);
   
   *ret_cumranklist = cumranklist;
   return eslOK;

 ERROR:
   if (shmsa) esl_msa_Destroy(shmsa);
   if (ranklist) cov_FreeRankList(ranklist);
   if (cumranklist) cov_FreeRankList(cumranklist);
   return status;
}

/* mull1(b) + null2 */
static int
null3_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *shmsa = NULL;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int       *useme1 = NULL;
  int       *useme2 = NULL;
  int        n;                   // index for alignment positions
  int        s;                   // index for shuffles
  int        b;                   // index for ranked list
  int        status;

  /* shuffle all the columns */
  ESL_ALLOC(useme1, sizeof(int) * msa->alen);
  ESL_ALLOC(useme2, sizeof(int) * msa->alen);
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] >  0) useme1[n] = TRUE; else useme1[n] = FALSE;
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] == 0) useme2[n] = TRUE; else useme2[n] = FALSE;

  for (s = 0; s < cfg->nshuffle; s ++) {
    msamanip_ShuffleColumns     (cfg->r, msa,   &shmsa, useme1, cfg->errbuf, cfg->verbose);
    msamanip_ShuffleColumns     (cfg->r, shmsa, &shmsa, useme2, cfg->errbuf, cfg->verbose);
    msamanip_ShuffleWithinColumn(cfg->r, shmsa, &shmsa,         cfg->errbuf, cfg->verbose);

    status = run_rnacov(go, cfg, &shmsa, NULL, NULL, &ranklist);
    if (status != eslOK) ESL_XFAIL(eslFAIL, "%s.\nFailed to run rnacov shuffled", cfg->errbuf);

     status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
     if (status != eslOK) goto ERROR;

    esl_msa_Destroy(shmsa); shmsa = NULL;
    cov_FreeRankList(ranklist); ranklist = NULL;
  }
  
  for (b = cumranklist->ha->imin; b <= cumranklist->ha->imax; b ++) {
    cumranklist->covBP[b]  /= (double)cfg->nshuffle;
    cumranklist->covNBP[b] /= (double)cfg->nshuffle;
  }
  
   if (cfg->verbose) {
     printf("null3 distribution - cummulative\n");
     printf("imin %d imax %d xmax %f xmin %f\n", 
	    cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
     //esl_histogram_Plot(stdout, cumranklist->h);
     //esl_histogram_PlotSurvival(stdout, cumranklist->h);
   }

  *ret_cumranklist = cumranklist;
  free(useme1);
  free(useme2);
  return eslOK;
  
 ERROR:
  if (useme1) free(useme1);
  if (useme2) free(useme2);
  if (shmsa) esl_msa_Destroy(shmsa);
  if (ranklist) cov_FreeRankList(ranklist);
  if (cumranklist) cov_FreeRankList(cumranklist);
  return status;
}


/* use a tree to generate residues independently for each alignment column */
static int
null4_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  ESL_MSA   *allmsa = NULL;
  ESL_MSA   *shmsa = NULL;
  MSA_STAT   shmstat;
  RANKLIST  *cumranklist = NULL;
  RANKLIST  *ranklist = NULL;
  int       sc;
  int       s;
  int       b;
  int       status;
  
  status = create_tree(go, cfg, msa);
  if (status != eslOK) goto ERROR;
  if (cfg->T == NULL) {
    if (msa->nseq == 1) 
      return eslOK;
    else 
      return eslFAIL;
  }
 
  status = Tree_FitchAlgorithmAncenstral(cfg->r, cfg->T, msa, &allmsa, &sc, cfg->errbuf, cfg->verbose);
  if (status != eslOK) goto ERROR;
  if (cfg->verbose) {
    eslx_msafile_Write(stdout, allmsa, eslMSAFILE_STOCKHOLM); 
    printf("fitch sc %d\n", sc);
  }

  for (s = 0; s < cfg->nshuffle; s ++) {
    status = msamanip_ShuffleTreeSubstitutions(cfg->r, cfg->T, msa, allmsa, &shmsa, cfg->errbuf, cfg->verbose);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s.\nFailed to run null4 rnacov", cfg->errbuf);
    
    if (cfg->verbose) {
      msamanip_DumpStats(stdout, msa, cfg->mstat);
      //eslx_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM); 
 
      msamanip_XStats(shmsa, &shmstat);
      msamanip_DumpStats(stdout, shmsa, shmstat);
      //eslx_msafile_Write(stdout, shmsa, eslMSAFILE_STOCKHOLM); 
    }

    status = run_rnacov(go, cfg, &shmsa, NULL, NULL, &ranklist);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "%s.\nFailed to run null4 rnacov", cfg->errbuf);
    if (shmsa == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "error creating shmsa");
    
    status = null_add2cumranklist(ranklist, &cumranklist, cfg->verbose, cfg->errbuf);
    if (status != eslOK) goto ERROR;
    
    esl_msa_Destroy(shmsa); shmsa = NULL;
    cov_FreeRankList(ranklist); ranklist = NULL;
  }
  
  for (b = cumranklist->ha->imin; b <= cumranklist->ha->imax; b ++) {
    cumranklist->covBP[b]  /= (double)cfg->nshuffle;
    cumranklist->covNBP[b] /= (double)cfg->nshuffle;
  }
  if (cfg->verbose) {
    printf("null4 distribution - cumulative\n");
    printf("imin %d imax %d xmax %f xmin %f\n", 
	   cumranklist->ha->imin, cumranklist->ha->imax, cumranklist->ha->xmax, cumranklist->ha->xmin);
    //esl_histogram_Plot(stdout, cumranklist->h);
    //esl_histogram_PlotSurvival(stdout, cumranklist->h);
  }
  
  if (cfg->verbose) cov_DumpRankList(stdout, cumranklist);
  
  *ret_cumranklist = cumranklist;

  esl_msa_Destroy(allmsa);
  return eslOK;
  
 ERROR:
  if (allmsa) esl_msa_Destroy(allmsa);
  if (shmsa) esl_msa_Destroy(shmsa);
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
  
  if (cumranklist == NULL) {
    cumranklist = cov_CreateRankList(ranklist->ha->bmax, ranklist->ha->bmin, ranklist->ha->w);
    cumranklist->ha->n     = ranklist->ha->n;

    cumranklist->ha->xmin  = ranklist->ha->xmin;
    cumranklist->ha->xmax  = ranklist->ha->xmax;
       cumranklist->ha->imin  = ranklist->ha->imin;
       cumranklist->ha->imax  = ranklist->ha->imax;
       cumranklist->scthresh  = ranklist->scthresh;
  }
  else {                    
    cov_GrowRankList(&cumranklist, ranklist->ha->bmax, ranklist->ha->bmin);
    cumranklist->scthresh = ESL_MIN(cumranklist->scthresh, ranklist->scthresh);
    cumranklist->ha->n    += ranklist->ha->n;
  }
  for (b = ranklist->ha->imin; b <= ranklist->ha->imax; b ++) {
    cov_ranklist_Bin2Bin(b, ranklist->ha, cumranklist->ha, &cumb);
    
    if (cumb >= cumranklist->ha->imin && cumb <= cumranklist->ha->imax) {
      cumranklist->covBP[cumb]  += ranklist->covBP[b];
      cumranklist->covNBP[cumb] += ranklist->covNBP[b];
      if (b >= ranklist->ha->imin && b <= ranklist->ha->imax) {
	cumranklist->ha->obs[cumb] += ranklist->ha->obs[b];
	cumranklist->ha->Nc        += ranklist->ha->obs[b];
	cumranklist->ha->No        += ranklist->ha->obs[b];
      }
    }
  }
  
  if (verbose) {
    printf("null distribution \n");
    printf("imin %d imax %d xmax %f xmin %f\n", 
	   ranklist->ha->imin, ranklist->ha->imax, ranklist->ha->xmax, ranklist->ha->xmin);
    //esl_histogram_Plot(stdout, ranklist->h);
    //esl_histogram_PlotSurvival(stdout, ranklist->h);
  }

  *ocumranklist = cumranklist;
  return eslOK;
}

