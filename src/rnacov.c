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
#include "covariation.h"
#include "covgrammars.h"
#include "ribosum_matrix.h"



#define ALPHOPTS     "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
#define METHODOPTS   "--naive,--phylo,--dca,--akmaev"              
#define COVTYPEOPTS  "--CHI,--CHIa,--CHIp,--GT,--GTa,--GTp,--MI,--MIp,--MIa,--MIr,--MIrp,--MIra,--MIg,--MIgp,--MIga,--OMES,--OMESp,--OMESa,--ALL"              
#define COVCLASSOPTS "--C16,--C2"                                          
#define NULLOPTS     "--null1,--null1b,--null2,--null3,--null4"                                          
#define THRESHOPTS   "--covNBP,--covNBPu,--covNBPf,--covRBP,--covRBPu,--covRBPf"                                          

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

  double           minidthresh;	       /* fractional minimal identity threshold for selecting subset of sequences */

  COVTYPE          covtype;
  COVCLASS         covclass;

  int              nseqthresh;

  int              nmsa;
  char            *msafile;
  char            *filename;
  char            *msaname;
  int             *msamap;

  FILE            *outmsafp;
 
  char            *outfile;
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */
   int              infmt;
 
  char            *R2Rversion;
  int              R2Rall;
  char            *R2Rfile;
  FILE            *R2Rfp;

  char            *R2Rcykfile;
  FILE            *R2Rcykfp;
  int              minloop;
  enum grammar_e   grammar;
  
  char            *dplotfile;
  char            *cykdplotfile;

  int              domsa;
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

  THRESH          *thresh;

  float            tol;
  int              verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      1 },
  { "--outdir",     eslARG_STRING,       NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "specify a directory for all output files",                                                  1 },
  { "--cykcov",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "obtain the structure with maximum covariation",                                             1 },
  { "--r2rall",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "make R2R plot all position in the alignment",                                               1 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                1 },
 /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      1 },
  { "-I",             eslARG_REAL,    "0.97",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < <x> id",                                                             1 },
  { "-i",             eslARG_REAL,      NULL,    NULL, "0<=x<1.0",   NULL,    NULL,  NULL,               "require seqs to have >= <x> id",                                                            1 },
  { "--submsa",       eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   1 },
  { "--nseqmin",      eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "minimum number of sequences in the alignment",                                              1 },  
  { "--gapthresh",    eslARG_REAL,      NULL,    NULL,  "0<=x<=1",   NULL,    NULL,  NULL,               "keep columns with < <x> fraction of gaps",                                                  1 },
  { "--minid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      1 },
  { "--maxid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",                                                      1 },
  /* msa format */
  { "--informat",   eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            1 },
   /* different ways to assess significance */
  { "--covNBP",      eslARG_REAL,       NULL,   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov NonBPs:    max total        (covNPB)      allowed",                                     1 },
  { "--covNBPu",     eslARG_REAL,     "0.05",   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov NonBPs;    max per_position (covNBP/alen) allowed",                                     1 },
  { "--covNBPf",     eslARG_REAL,       NULL,   NULL,    "0<x<=1",THRESHOPTS, NULL,  NULL,               "cov NonBPs;    max fraction     (covNBP/NBP)  allowed",                                     1 },
  { "--covRBP",      eslARG_REAL,       NULL,   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov RandomBPs: max total        (covRPB)      allowed",                                     1 },
  { "--covRBPu",     eslARG_REAL,       NULL,   NULL,      "x>=0",THRESHOPTS, NULL,  NULL,               "cov RandomBPs; max per_position (covRBP/alen) allowed",                                     1 },
  { "--covRBPf",     eslARG_REAL,       NULL,   NULL,   "0<=x<=1",THRESHOPTS, NULL,  NULL,               "cov RandomBPs; max fraction     (covRBP/RBP)  allowed",                                     1 },
  /* null hypothesis */
  { "--nshuffle",      eslARG_INT,       "20",   NULL,      "n>0",   NULL,    NULL,  NULL,               "number of shuffled sequences",                                                              1 },   
  { "--null1",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null1:  shuffle alignment columns",                                                         0 },
  { "--null1b",       eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null1b: shuffle bpaired_columns and nonbpaired_columns independently",                      0 },
  { "--null2",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null2:  shuffle residues within a column",                                                  0 },
  { "--null3",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null3:  null1(b)+null2",                                                                    0 },
  { "--null4",        eslARG_NONE,      FALSE,   NULL,       NULL,  NULLOPTS, NULL,  NULL,               "null4: ",                                                                                   0 },
  /* covariation measures */
  { "--CHIa",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "CHI  ACS corrected calculation",                                                            0 },
  { "--CHIp",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "CHI  APS corrected calculation",                                                            0 },
  { "--CHI",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "CHI  calculation",                                                                          0 },
  { "--GTa",          eslARG_NONE,       TRUE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   ACS corrected calculation",                                                            0 },
  { "--GTp",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   APS corrected calculation",                                                            0 },
  { "--GT",           eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "GT   calculation",                                                                          0 },
  { "--MIa",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   ACS corrected calculation",                                                            0 },
  { "--MIp",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   APS corrected calculation",                                                            0 },
  { "--MI",           eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MI   calculation",                                                                          0 },
  { "--MIra",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIr  ACS corrected calculation",                                                            0 },
  { "--MIrp",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIr  APS corrected calculation",                                                            0 },
  { "--MIr",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIr  calculation",                                                                          0 },
  { "--MIga",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIg  ACS corrected calculation",                                                            0 },
  { "--MIgp",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIg  APS corrected calculation",                                                            0 },
  { "--MIg",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "MIg  calculation",                                                                          0 },
  { "--OMESa",        eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "OMES ACS corrected calculation",                                                            0 },
  { "--OMESp",        eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "OMES APS corrected calculation",                                                            0 },
  { "--OMES",         eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "OMES calculation",                                                                          0 },
  { "--ALL",          eslARG_NONE,      FALSE,   NULL,       NULL,COVTYPEOPTS, NULL,  NULL,              "calculate all types",                                                                       0 },
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
  { "--minloop",       eslARG_INT,       "5",    NULL,      "n>0",   NULL,    NULL, NULL,                "minloop in cykcov calculation",                                                             0 },   
  { "--grammar",    eslARG_STRING,     "BGR",    NULL,       NULL,   NULL,"--cykcov", NULL,              "grammar used for cococyk calculation",                                                      0 },   
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,      "42",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msa>";
static char banner[] = "rnacov - statistical test for covatiation in RNA alignments";

static int create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int run_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **msa, RANKLIST *ranklist_null, RANKLIST **ret_ranklist);
static int null1_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null1b_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null2_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null3_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);
static int null4_rnacov (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_ranklist_null);

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

  cfg.w = esl_stopwatch_Create(); 
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));

   /* outmsa file */
  cfg.outmsafp = NULL;
  if (esl_opt_IsOn(go, "--outmsa")) {
    if ((cfg.outmsafp = fopen(esl_opt_GetString(go, "--outmsa"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "--outmsa"));
  } 
  
  esl_FileTail(cfg.msafile, TRUE, &cfg.outheader);
  if ( esl_opt_IsOn(go, "--outdir") ) esl_sprintf( &cfg.outheader, "%s/%s", esl_opt_GetString(go, "--outdir"), cfg.outheader);
  
  if (esl_opt_IsOn(go, "--submsa")) { cfg.submsa = esl_opt_GetInteger(go, "--submsa"); esl_sprintf(&cfg.outheader, "%s.select%d", cfg.outheader, cfg.submsa); }
  else                              { cfg.submsa = 0; }
    
  /* other options */
  cfg.domsa       = TRUE;
  cfg.nshuffle    = esl_opt_GetInteger(go, "--nshuffle");
  cfg.nseqthresh  = esl_opt_GetInteger(go, "--nseqthresh");
  cfg.fragfrac    = esl_opt_IsOn(go, "-F")?           esl_opt_GetReal   (go, "-F")           : -1.0;
  cfg.idthresh    = esl_opt_IsOn(go, "-I")?           esl_opt_GetReal   (go, "-I")           : -1.0;
  cfg.minidthresh = esl_opt_IsOn(go, "-i")?           esl_opt_GetReal   (go, "-i")           : -1.0;
  cfg.nseqmin     = esl_opt_IsOn(go, "--nseqmin")?    esl_opt_GetInteger(go, "--nseqmin")    : -1;
  cfg.gapthresh   = esl_opt_IsOn(go, "--gapthresh")?  esl_opt_GetReal   (go, "--gapthresh")  :  1.0;
  cfg.tol         = esl_opt_GetReal   (go, "--tol");
  cfg.verbose     = esl_opt_GetBoolean(go, "-v");
  cfg.voutput     = esl_opt_GetBoolean(go, "--voutput");
  cfg.minloop     = esl_opt_GetInteger(go, "--minloop");
  
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
  else if (esl_opt_GetBoolean(go, "--null3"))  cfg.nulltype = Null3;
  else if (esl_opt_GetBoolean(go, "--null4"))  cfg.nulltype = Null4;

  ESL_ALLOC(cfg.thresh, sizeof(THRESH));
  if      (esl_opt_IsOn(go, "--covNBP") )  { cfg.thresh->type = covNBP;  cfg.thresh->val = esl_opt_GetReal(go, "--covNBP");  }
  else if (esl_opt_IsOn(go, "--covNBPu"))  { cfg.thresh->type = covNBPu; cfg.thresh->val = esl_opt_GetReal(go, "--covNBPu"); }
  else if (esl_opt_IsOn(go, "--covNBPf"))  { cfg.thresh->type = covNBPf; cfg.thresh->val = esl_opt_GetReal(go, "--covNBPf"); }
  else if (esl_opt_IsOn(go, "--covRBP") )  { cfg.thresh->type = covRBP;  cfg.thresh->val = esl_opt_GetReal(go, "--covRBP");  if (cfg.nulltype == NullNONE) cfg.nulltype = Null3; }
  else if (esl_opt_IsOn(go, "--covRBPu"))  { cfg.thresh->type = covRBPu; cfg.thresh->val = esl_opt_GetReal(go, "--covRBPu"); if (cfg.nulltype == NullNONE) cfg.nulltype = Null3; }
  else if (esl_opt_IsOn(go, "--covRBPf"))  { cfg.thresh->type = covRBPf; cfg.thresh->val = esl_opt_GetReal(go, "--covRBPf"); if (cfg.nulltype == NullNONE) cfg.nulltype = Null3; }

  if      (esl_opt_GetBoolean(go, "--CHIa"))  cfg.covtype = CHIa;
  else if (esl_opt_GetBoolean(go, "--CHIp"))  cfg.covtype = CHIp;
  else if (esl_opt_GetBoolean(go, "--CHI"))   cfg.covtype = CHI;
  else if (esl_opt_GetBoolean(go, "--GTa"))   cfg.covtype = GTa;
  else if (esl_opt_GetBoolean(go, "--GTp"))   cfg.covtype = GTp;
  else if (esl_opt_GetBoolean(go, "--GT"))    cfg.covtype = GT;
  else if (esl_opt_GetBoolean(go, "--MIa"))   cfg.covtype = MIa;
  else if (esl_opt_GetBoolean(go, "--MIp"))   cfg.covtype = MIp;
  else if (esl_opt_GetBoolean(go, "--MI"))    cfg.covtype = MI;
  else if (esl_opt_GetBoolean(go, "--MIra"))  cfg.covtype = MIra;
  else if (esl_opt_GetBoolean(go, "--MIrp"))  cfg.covtype = MIrp;
  else if (esl_opt_GetBoolean(go, "--MIr"))   cfg.covtype = MIr;
  else if (esl_opt_GetBoolean(go, "--MIga"))  cfg.covtype = MIga;
  else if (esl_opt_GetBoolean(go, "--MIgp"))  cfg.covtype = MIgp;
  else if (esl_opt_GetBoolean(go, "--MIg"))   cfg.covtype = MIg;
  else if (esl_opt_GetBoolean(go, "--OMESa")) cfg.covtype = OMESa;
  else if (esl_opt_GetBoolean(go, "--OMESp")) cfg.covtype = OMESp;
  else if (esl_opt_GetBoolean(go, "--OMES"))  cfg.covtype = OMES;
  else if (esl_opt_GetBoolean(go, "--ALL"))   cfg.covtype = COVALL;

 
  if      (esl_opt_GetBoolean(go, "--C16"))   cfg.covclass = C16;
  else if (esl_opt_GetBoolean(go, "--C2"))    cfg.covclass = C2;
  else                                        cfg.covclass = CSELECT;

  if      (esl_opt_GetBoolean(go, "--naive"))  cfg.method = OPTNONE;
  else if (esl_opt_GetBoolean(go, "--phylo"))  cfg.method = PHYLO;
  else if (esl_opt_GetBoolean(go, "--dca"))    cfg.method = DCA;
  else if (esl_opt_GetBoolean(go, "--akmaev")) cfg.method = AKMAEV;
 
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
  if ((cfg.rocfp = fopen(cfg.rocfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.rocfile);

  /*  summary file */
  esl_sprintf(&cfg.sumfile, "%s.sum", cfg.outheader); 
  if ((cfg.sumfp = fopen(cfg.sumfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.sumfile);
  
  /* R2R annotated sto file */
  esl_sprintf(&cfg.R2Rfile, "%s.R2R.sto", cfg.outheader);
  cfg.R2Rfp = NULL;

  /* dotplot file */
  esl_sprintf(&cfg.dplotfile, "%s.%s", cfg.outheader, "dplot.svg");
  /* dotplot file */
  cfg.cykdplotfile = NULL;
  if (esl_opt_IsOn(go, "--cykcov")) 
    esl_sprintf(&cfg.cykdplotfile, "%s.%s", cfg.outheader, "cyk.dplot.svg");
 
  cfg.R2Rcykfile = NULL;
  if (esl_opt_IsOn(go, "--cykcov")) {
    esl_sprintf(&cfg.R2Rcykfile, "%s.%s", cfg.outheader, "cyk.R2R.sto");
    cfg.R2Rcykfp = NULL;
  }


  cfg.shsumfile = NULL;
  cfg.shsumfp = NULL;

  if (cfg.nulltype != NullNONE) {
    /*  sh-summary file */
    esl_sprintf(&cfg.shsumfile, "%s.shsum", cfg.outheader); 
    if ((cfg.shsumfp = fopen(cfg.shsumfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.shsumfile);
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

int
main(int argc, char **argv)
{ 
  char            *msg = "rnacov failed";
  ESL_GETOPTS     *go;
  struct cfg_s     cfg;
  ESLX_MSAFILE    *afp = NULL;
  ESL_MSA         *msa = NULL;            /* the input alignment    */
  char            *type = NULL;
  char            *tp;
  char            *tok;
  RANKLIST        *ranklist_null = NULL;
  int              nmsa_noss = 0;
  int              seq_cons_len = 0;
  int              nfrags = 0;	  	  /* # of fragments removed */
  int              nremoved = 0;	  /* # of identical sequences removed */
  int              t;
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

    if (cfg.onbpairs == 0) {
      nmsa_noss ++;
      printf("msa: %s_%s  has no secondary structure.\n", msa->acc, msa->name);
      esl_msa_Destroy(msa); msa = NULL;
      continue;
    }

    /* write MSA info to the sumfile */
    for (t = 0; t < msa->ngf; t++) {
      if (!esl_strcmp(msa->gf_tag[t], "TP")) {
	tp = msa->gf[t];	
	while (*tp != '\0') {
	  if (esl_strtok(&tp, " ", &tok) != eslOK) esl_fatal(msg);
	  esl_strcat(&type, -1, tok, -1);
	}
      }
    }
    if      (msa->acc && msa->name && type) esl_sprintf(&cfg.msaname, "%s_%s_%s", msa->acc, msa->name, type);
    else if (msa->acc && msa->name)         esl_sprintf(&cfg.msaname, "%s_%s", msa->acc, msa->name);
    else if (msa->acc)                      esl_sprintf(&cfg.msaname, "%s", msa->acc);
    else                                    esl_sprintf(&cfg.msaname, "%s", cfg.outheader);

   /* select submsa and then apply msa filters 
    */
    if (cfg.fragfrac > 0.      && msamanip_RemoveFragments(cfg.fragfrac, &msa, &nfrags, &seq_cons_len)          != eslOK) { printf("remove_fragments failed\n");     esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-I") && msamanip_SelectSubsetBymaxID(cfg.r, &msa, cfg.idthresh, &nremoved)            != eslOK) { printf("select_subsetBymaxID failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-i") && msamanip_SelectSubsetByminID(cfg.r, &msa, cfg.minidthresh, &nremoved)         != eslOK) { printf("select_subsetByminID failed\n"); esl_fatal(msg); }
    if (cfg.submsa             && msamanip_SelectSubset(cfg.r, cfg.submsa, &msa, NULL, cfg.errbuf, cfg.verbose) != eslOK) { printf("%s\n", cfg.errbuf);              esl_fatal(msg); }
    if (msa == NULL) {
      free(type); type = NULL;
      free(cfg.msaname); cfg.msaname = NULL;
      continue;
    }
    if (msa->nseq < cfg.nseqmin) {
      esl_msa_Destroy(msa); msa = NULL;
      free(type); type = NULL;
      free(cfg.msaname); cfg.msaname = NULL;
      continue;  
    }
    if (msamanip_RemoveGapColumns(cfg.gapthresh, msa, &cfg.msamap, cfg.errbuf, cfg.verbose) != eslOK) { printf("RemoveGapColumns\n"); esl_fatal(msg); }
 
    esl_msa_Hash(msa);
    esl_msa_ConvertDegen2X(msa);
    if (esl_msa_MinimGaps(msa, NULL, "-.~=", FALSE) != eslOK) esl_fatal("Failed to remove minim gaps");
    
    /* given msa aveid and avematch */
    msamanip_XStats(msa, &cfg.mstat);
    
    if ((esl_opt_IsOn(go, "--minid") && cfg.mstat.avgid < 100.*esl_opt_GetReal(go, "--minid")) ||
	(esl_opt_IsOn(go, "--maxid") && cfg.mstat.avgid > 100.*esl_opt_GetReal(go, "--maxid"))   ) {
      esl_msa_Destroy(msa); msa = NULL;
      free(type); type = NULL;
      free(cfg.msaname); cfg.msaname = NULL;
      continue;  
    }

    /* output the actual file used if requested */
    if (cfg.outmsafp) eslx_msafile_Write(cfg.outmsafp, msa, eslMSAFILE_STOCKHOLM);
 
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
#if 0
    msamanip_CalculateBC(msa, cfg.ct, &cfg.ft, &cfg.fbp, &cfg.fnbp, cfg.errbuf);
    esl_vec_DDump(stdout, cfg.ft,   cfg.abc->K, "total BC");
    esl_vec_DDump(stdout, cfg.fbp,  cfg.abc->K, "basepairs BC");
    esl_vec_DDump(stdout, cfg.fnbp, cfg.abc->K, "nonbasepairs BC");
#endif
    
    /* the null model first */
    if (cfg.nulltype == Null1) {
      status = null1_rnacov(go, &cfg, msa, &ranklist_null);
     if (status != eslOK) esl_fatal("%s.\nFailed to run rnacov", cfg.errbuf);
    }
    if (cfg.nulltype == Null1b) {
      status = null1b_rnacov(go, &cfg, msa, &ranklist_null);
     if (status != eslOK) esl_fatal("%s.\nFailed to run rnacov", cfg.errbuf);
    }
    if (cfg.nulltype == Null2) {
      status = null2_rnacov(go, &cfg, msa, &ranklist_null);
     if (status != eslOK) esl_fatal("%s.\nFailed to run rnacov", cfg.errbuf);
    }
    if (cfg.nulltype == Null3) {
      status = null3_rnacov(go, &cfg, msa, &ranklist_null);
     if (status != eslOK) esl_fatal("%s.\nFailed to run rnacov", cfg.errbuf);
    }
    if (cfg.nulltype == Null4) {
      status = null4_rnacov(go, &cfg, msa, &ranklist_null);
     if (status != eslOK) esl_fatal("%s.\nFailed to run rnacov", cfg.errbuf);
    }
 
    /* main function */
    if (cfg.domsa) {
      status = run_rnacov(go, &cfg, &msa, ranklist_null, NULL);
      if (status != eslOK) esl_fatal("%s.\nFailed to run rnacov", cfg.errbuf);
    }

    if (msa) esl_msa_Destroy(msa); msa = NULL;
    free(cfg.ct); cfg.ct = NULL;
    if (cfg.msafrq) free(cfg.msafrq); cfg.msafrq = NULL;
    if (cfg.T) esl_tree_Destroy(cfg.T); cfg.T = NULL;
    if (cfg.msaname) free(cfg.msaname); cfg.msaname = NULL;
    if (type) free(type); type = NULL;
    COV_FreeRankList(ranklist_null); ranklist_null = NULL;
  }

  if (nmsa_noss > 0) printf("%d msa's without any secondary structure\n", nmsa_noss);

  /* cleanup */
  esl_stopwatch_Destroy(cfg.w);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  if (cfg.ct) free(cfg.ct);
  eslx_msafile_Close(afp);
  if (cfg.msaname) free(cfg.msaname);
  if (type) free(type);
  if (cfg.outfile) free(cfg.outfile);
  fclose(cfg.outfp);
  fclose(cfg.rocfp);
  fclose(cfg.sumfp);
  if (cfg.outmsafp) fclose(cfg.outmsafp);
  if (cfg.shsumfp) fclose(cfg.shsumfp);
  if (cfg.R2Rfp) fclose(cfg.R2Rfp);
  free(cfg.outheader);
  free(cfg.dplotfile);
  if (cfg.cykdplotfile) free(cfg.cykdplotfile);
  free(cfg.R2Rfile);
  free(cfg.R2Rversion);
  free(cfg.gnuplot);
  if (cfg.ribosum) Ribosum_matrix_Destroy(cfg.ribosum);
  if (cfg.ft)   free(cfg.ft);
  if (cfg.fbp)  free(cfg.fbp);
  if (cfg.fnbp) free(cfg.fnbp);
  if (cfg.thresh) free(cfg.thresh);
  if (ranklist_null) COV_FreeRankList(ranklist_null);

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
run_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa, RANKLIST *ranklist_null, RANKLIST **ret_ranklist)
{
  struct mutual_s *mi   = NULL;
  ESL_MSA         *msa = *omsa;
  int              nnodes;
  int              status;
  int              ishuffled = (ret_ranklist)? TRUE:FALSE;

  esl_stopwatch_Start(cfg->w);
  
  /* print to stdout */
  if (cfg->verbose) fprintf(stdout, "# MSA %s nseq %d (%d) alen %" PRId64 " (%" PRId64 ") avgid %.2f (%.2f) nbpairs %d (%d)\n", 
			   cfg->msaname, msa->nseq, cfg->omstat.nseq, msa->alen, cfg->omstat.alen, 
			   cfg->mstat.avgid, cfg->omstat.avgid, cfg->nbpairs, cfg->onbpairs);  
  if (!ishuffled) 
    fprintf(cfg->outfp, "# MSA %s nseq %d (%d) alen %" PRId64 " (%" PRId64 ") avgid %.2f (%.2f) nbpairs %d (%d)\n", 
	    cfg->msaname, msa->nseq, cfg->omstat.nseq, msa->alen, cfg->omstat.alen, 
	    cfg->mstat.avgid, cfg->omstat.avgid, cfg->nbpairs, cfg->onbpairs);  
  
  /* produce a tree
   */
  if (cfg->method != NAIVE) {
    status = create_tree(go, cfg, msa);
    if (status != eslOK)  { esl_fatal(cfg->errbuf); }
    nnodes = (cfg->T->N > 1)? cfg->T->N-1 : cfg->T->N;
  }
  
  /* create the MI structure */
  mi = COV_Create(msa->alen, msa->nseq, ishuffled, cfg->nseqthresh, cfg->abc);
  
  /* write MSA info to the sumfile */
  if (!ishuffled) 
    fprintf(cfg->sumfp, "%f\t%s\t%d\t%d\t%.2f\t", cfg->thresh->val, cfg->msaname, msa->nseq, (int)msa->alen, cfg->mstat.avgid); 
  else
    fprintf(cfg->shsumfp, "%f\t%s\t%d\t%d\t%.2f\t", cfg->thresh->val, cfg->msaname, msa->nseq, (int)msa->alen, cfg->mstat.avgid); 
  
  /* write MSA info to the rocfile */
  fprintf(cfg->rocfp, "# MSA nseq %d alen %" PRId64 " avgid %f nbpairs %d (%d)\n", msa->nseq, msa->alen, cfg->mstat.avgid, cfg->nbpairs, cfg->onbpairs);  
  
  /* main function */
  status = COV_Calculate(&msa, cfg->msamap, cfg->T, cfg->ribosum, mi, ranklist_null, ret_ranklist, cfg->method, cfg->covtype, cfg->covclass, cfg->ct, 
			 (ret_ranklist)?NULL:cfg->outfp, cfg->rocfp, 
			 (ret_ranklist)?cfg->shsumfp:cfg->sumfp, cfg->gnuplot, cfg->dplotfile, cfg->R2Rfile, cfg->R2Rversion, cfg->R2Rall, 
			 cfg->thresh, cfg->onbpairs, cfg->tol, cfg->verbose, cfg->errbuf);   
  if (status != eslOK)  { goto ERROR; }
  
  /* find the cykcov structure, and do the cov analysis on it */
  if (cfg->R2Rcykfile && !ishuffled) {
    status = COV_CYKCOVCT(cfg->outfp, cfg->gnuplot, cfg->cykdplotfile, cfg->R2Rcykfile, cfg->R2Rversion, cfg->R2Rall, cfg->r, 
			  &msa, mi, cfg->msamap, cfg->minloop, cfg->grammar, cfg->thresh, cfg->onbpairs, cfg->errbuf, cfg->verbose);
    if (status != eslOK)  { goto ERROR; }
  }
 
  *omsa = msa;
  COV_Destroy(mi); mi = NULL;
  
  return eslOK;
  
 ERROR:
  if (cfg->T) esl_tree_Destroy(cfg->T);
  if (mi) COV_Destroy(mi);
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
  int        x;                   // index for ranked list
  int        cumx;
  int        status;

  /* shuffle all the columns */
  ESL_ALLOC(useme, sizeof(int) * msa->alen);
  for (n = 0; n < msa->alen; n ++) useme[n] = TRUE;

  for (s = 0; s < cfg->nshuffle; s ++) {
    msamanip_ShuffleColumns(cfg->r, msa, &shmsa, useme, cfg->errbuf, cfg->verbose);
    status = run_rnacov(go, cfg, &shmsa, NULL, &ranklist);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run rnacov shuffled");
    esl_msa_Destroy(shmsa); shmsa = NULL;
    
    if (cumranklist == NULL) cumranklist = COV_CreateRankList(ranklist->bmax, ranklist->bmin, ranklist->w);
    else                     COV_GrowRankList(&cumranklist, ranklist->bmax, ranklist->bmin);
    
    for (x = 0; x < ranklist->nb; x ++) {
      cumx = (int)((x*ranklist->w + ranklist->bmin - cumranklist->bmin)/cumranklist->w);
      cumranklist->covBP[cumx]  += ranklist->covBP[x];
      cumranklist->covNBP[cumx] += ranklist->covNBP[x];
    }
    
    COV_FreeRankList(ranklist); ranklist = NULL;
  }
  
  for (x = 0; x < cumranklist->nb; x ++) {
    cumranklist->covBP[x]  /= (double)cfg->nshuffle;
    cumranklist->covNBP[x] /= (double)cfg->nshuffle;
  }
  
  *ret_cumranklist = cumranklist;
  free(useme);
  return eslOK;
  
 ERROR:
  if (useme) free(useme);
  if (shmsa) esl_msa_Destroy(shmsa);
  if (ranklist) COV_FreeRankList(ranklist);
  if (cumranklist) COV_FreeRankList(cumranklist);
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
  int        x;                   // index for ranked list
  int        cumx;
  int        status;

  /* shuffle all the columns */
  ESL_ALLOC(useme1, sizeof(int) * msa->alen);
  ESL_ALLOC(useme2, sizeof(int) * msa->alen);
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] >  0) useme1[n] = TRUE; else useme1[n] = FALSE;
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] == 0) useme2[n] = TRUE; else useme2[n] = FALSE;

  for (s = 0; s < cfg->nshuffle; s ++) {
    msamanip_ShuffleColumns(cfg->r, msa,   &shmsa, useme1, cfg->errbuf, cfg->verbose);
    msamanip_ShuffleColumns(cfg->r, shmsa, &shmsa, useme2, cfg->errbuf, cfg->verbose);
    status = run_rnacov(go, cfg, &shmsa, NULL, &ranklist);
    if (status != eslOK) goto ERROR;
    esl_msa_Destroy(shmsa); shmsa = NULL;
    
    if (cumranklist == NULL) cumranklist = COV_CreateRankList(ranklist->bmax, ranklist->bmin, ranklist->w);
    else                     COV_GrowRankList(&cumranklist, ranklist->bmax, ranklist->bmin);

    for (x = 0; x < ranklist->nb; x ++) {
      cumx = (int)((x*ranklist->w + ranklist->bmin - cumranklist->bmin)/cumranklist->w);
      cumranklist->covBP[cumx]  += ranklist->covBP[x];
      cumranklist->covNBP[cumx] += ranklist->covNBP[x];
    }
    
    COV_FreeRankList(ranklist); ranklist = NULL;
  }
  
  for (x = 0; x < cumranklist->nb; x ++) {
    cumranklist->covBP[x]  /= (double)cfg->nshuffle;
    cumranklist->covNBP[x] /= (double)cfg->nshuffle;
  }
  
  *ret_cumranklist = cumranklist;
  free(useme1);
  free(useme2);
  return eslOK;
  
 ERROR:
  if (useme1) free(useme1);
  if (useme2) free(useme2);
  if (shmsa) esl_msa_Destroy(shmsa);
  if (ranklist) COV_FreeRankList(ranklist);
  if (cumranklist) COV_FreeRankList(cumranklist);
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
  int       x;
  int       cumx;
  int       status;

   for (s = 0; s < cfg->nshuffle; s ++) {

     msamanip_ShuffleWithinColumn(cfg->r, msa, &shmsa, cfg->errbuf, cfg->verbose);
     status = run_rnacov(go, cfg, &shmsa, NULL, &ranklist);
     if (status != eslOK) ESL_XFAIL(eslFAIL, "%s.\nFailed to run rnacov shuffled", cfg->errbuf);
     esl_msa_Destroy(shmsa); shmsa = NULL;
     
     if (cumranklist == NULL) cumranklist = COV_CreateRankList(ranklist->bmax, ranklist->bmin, ranklist->w);
     else                     COV_GrowRankList(&cumranklist, ranklist->bmax, ranklist->bmin);
     
     for (x = 0; x < ranklist->nb; x ++) {
       cumx = (int)((x*ranklist->w + ranklist->bmin - cumranklist->bmin)/cumranklist->w);
       cumranklist->covBP[cumx]  += ranklist->covBP[x];
       cumranklist->covNBP[cumx] += ranklist->covNBP[x];
     }
     
     COV_FreeRankList(ranklist); ranklist = NULL;
   }
   
   for (x = 0; x < cumranklist->nb; x ++) {
     cumranklist->covBP[x]  /= (double)cfg->nshuffle;
     cumranklist->covNBP[x] /= (double)cfg->nshuffle;
   }

   *ret_cumranklist = cumranklist;
   return eslOK;

 ERROR:
   if (shmsa) esl_msa_Destroy(shmsa);
   if (ranklist) COV_FreeRankList(ranklist);
   if (cumranklist) COV_FreeRankList(cumranklist);
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
  int        x;                   // index for ranked list
  int        cumx;
  int        status;

  /* shuffle all the columns */
  ESL_ALLOC(useme1, sizeof(int) * msa->alen);
  ESL_ALLOC(useme2, sizeof(int) * msa->alen);
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] >  0) useme1[n] = TRUE; else useme1[n] = FALSE;
  for (n = 0; n < msa->alen; n ++) if (cfg->ct[n+1] == 0) useme2[n] = TRUE; else useme2[n] = FALSE;

  for (s = 0; s < cfg->nshuffle; s ++) {
    msamanip_ShuffleWithinColumn(cfg->r, msa,   &shmsa,         cfg->errbuf, cfg->verbose);
    msamanip_ShuffleColumns     (cfg->r, shmsa, &shmsa, useme1, cfg->errbuf, cfg->verbose);
    msamanip_ShuffleColumns     (cfg->r, shmsa, &shmsa, useme2, cfg->errbuf, cfg->verbose);
    status = run_rnacov(go, cfg, &shmsa, NULL, &ranklist);
    if (status != eslOK) ESL_XFAIL(eslFAIL, "%s.\nFailed to run rnacov shuffled", cfg->errbuf);
    esl_msa_Destroy(shmsa); shmsa = NULL;
    
    if (cumranklist == NULL) cumranklist = COV_CreateRankList(ranklist->bmax, ranklist->bmin, ranklist->w);
    else                     COV_GrowRankList(&cumranklist, ranklist->bmax, ranklist->bmin);
    
    for (x = 0; x < ranklist->nb; x ++) {
      cumx = (int)((x*ranklist->w + ranklist->bmin - cumranklist->bmin)/cumranklist->w); 
      cumranklist->covBP[cumx]  += ranklist->covBP[x];
      cumranklist->covNBP[cumx] += ranklist->covNBP[x];
    }
    
    COV_FreeRankList(ranklist); ranklist = NULL;
  }
  
  for (x = 0; x < cumranklist->nb; x ++) {
    cumranklist->covBP[x]  /= (double)cfg->nshuffle;
    cumranklist->covNBP[x] /= (double)cfg->nshuffle;
  }
  
  *ret_cumranklist = cumranklist;
  free(useme1);
  free(useme2);
  return eslOK;
  
 ERROR:
  if (useme1) free(useme1);
  if (useme2) free(useme2);
  if (shmsa) esl_msa_Destroy(shmsa);
  if (ranklist) COV_FreeRankList(ranklist);
  if (cumranklist) COV_FreeRankList(cumranklist);
  return status;
}


/* use a tree to generate residues independently for each alignment column */
static int
null4_rnacov(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, RANKLIST **ret_cumranklist)
{
  
  printf("null4 not implemented yet\n");
  exit(1);
  return eslOK;
  
}
