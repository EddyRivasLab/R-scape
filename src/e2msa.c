/* e2msa -- progressive alignment using an evolutionary model
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
#include "e2_profilesq.h"
#include "e2_msa.h"
#include "e2_tree.h"
#include "evohmmer.h"
#include "miscellaneous.h"
#include "msatree.h"
#include "msaprobs.h"
#include "msamanip.h"
#include "muscle.h"
#include "ncbiblast.h"
#include "fasta.h"
#include "phmmer.h"
#include "fastsp.h"

#define ALPHOPTS "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
#define ALIOPTS  "--e2,--e2f,--e2hmmer,--e2fhmmer"          /* Exclusive options for alignment choice */
#define OPTOPTS  "--optnone,--opttime,--optpara,--optboth"   /* Exclusive options for msa optimization */
#define EVOMOPTS "--LI,--LR,--AFG,--AFGR,--AFR,--AIF,-GG,--G2,--TKF91,--TKF92,--FID,--GTKF92,--GRTKF92,--ITKF92"              
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

  char            *benchfile;
  FILE            *benchfp;
  
  char            *paramfile;

  int              voutput;             /* TRUE for verbose output */

  int              nmsa;
  char            *msafile;
  char            *gnuplot;
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */
  char            *msaheader;          /* header for all msa-specific output files */

  char            *msarfname;
  char            *msaefname;
  char            *msaefname_ephmmer;
  char            *msaefname_phmmer;
  char            *msaefname_phmmer3;
  char            *msaefname_ncbiblast;
  char            *msaefname_ssearch;
  char            *msaefname_muscle;
  char            *msaefname_msaprobs;
  FILE            *msaefp; 
  FILE            *msaefp_ephmmer; 
  FILE            *msaefp_phmmer; 
  FILE            *msaefp_phmmer3; 
  FILE            *msaefp_ncbiblast; 
  FILE            *msaefp_ssearch; 
  FILE            *msaefp_muscle; 
  FILE            *msaefp_msaprobs; 

  char            *hmmfile;            /* for e2hmmer case: open input HMM file */
  
  double           treeavgt;
  ESL_TREE        *T;

  P7_HMMFILE      *hfp;          /* input HMM file */
  P7_HMM          *hmm;          /* HMM            */
  EMRATE          *emR;

  E2_PIPELINE     *pli;
  E1_BG           *bg;	   	       /* null model (copies made of this into threads) */
  P7_BG           *bg7;	   	       /* null model (copies made of this into threads) */
  int              mode;               /* e2_GLOBAL, e2_LOCAL */
  int              do_viterbi;         /* default is fwd/bck and optimal accuracy */
 
  E1_RATE        **R1;
  P7_RATE         *R7;
  int              nr;
  char            *subsmx;             /* BLUSUM62, ... */
  int              userates;           /* if TRUE construct the R from given rate values */
  
  E2_ALI             e2ali;
  EVOM               evomodel;
  struct rateparam_s rateparam;
  double             popen;
  double             pextend;
  double             pcross;
  float              etainf;

  int                submsa;              /* set to the number of random seqs taken from original msa.
					   * Set to 0 if we are taking all */
  MSA_STAT          *mstat;               /* statistics of the input alignment */
  float             *msafrq;

  int                infmt;
  E2_OPT             e2optimize;          /* type of optimization (parameters and/or time) */

  float              fixpid;
  float              fixtime;

  int                EmL;
  int                EmN;
  int                EvL;
  int                EvN;
  int                EfL;
  int                EfN;
  float              Eft;

  char              *mx;
  int                dostats;

  int                doephmmer;           /* TRUE to run ephmmer at the end to compare alignments */
  int                dophmmer;            /* TRUE to run phmmer at the end to compare alignments */
  int                dophmmer3;           /* TRUE to run phmmer at the end to compare alignments */
  int                doncbiblast;         /* TRUE to run ncbiblast at the end to compare alignments */
  int                dossearch;           /* TRUE to run ssearch at the end to compare alignments */
  int                domuscle;            /* TRUE to run muscle at the end to compare alignments */
  int                domsaprobs;          /* TRUE to run MSAProbs at the end to compare alignments */
  int                doe2msa;

  int                dofastsp;            /* compare to original tree */
  int                lcu_r;               /* TRUE for lc in reference unaligned with FastSP */
  int                lcu_e;               /* TRUE for lc in inferrend unaligned with FastSP */

  int                wordsize;            /* for running blast -1 == default */
  int                max;                 /* for running --max with phmmer or hmmsearch */

  char              *fastaversion;        /* fasta-36.3.6d, or ... for running ssearch */
  char              *fasta_s;             /* fasta scoring matrix */
  int                fasta_f;             /* gap-open penalty */
  int                fasta_g;             /* gap-extend penalty */

  char              *method;
  int                probdist;            /* explore the distribution as a function of time */
  int                train;               /* train on given msas */
  float              tol;
  int                verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  { "--ephmmer",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run ephmmer    at the end to compare alignments",                                   0 }, 
  { "--phmmer",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run phmmer     at the end to compare alignments",                                   0 }, 
  { "--phmmer3",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run phmmer3    at the end to compare alignments",                                   0 }, 
  { "--ncbiblast",    eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run ncbiblast  at the end to compare alignments",                                   0 }, 
  { "--ssearch",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run ssearch  at the end to compare alignments",                                     0 }, 
  { "--muscle",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run muscle at the end to compare alignments",                                       0 }, 
  { "--msaprobs",     eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run MSAProbs at the end to compare alignments",                                     0 }, 
  { "--onlyephmmer",  eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run ephmmer alone",                                                                 0 }, 
  { "--onlyphmmer",   eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run phmmer alone",                                                                  0 }, 
  { "--onlyphmmer3",  eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run phmmer3 alone",                                                                 0 }, 
  { "--onlyncbiblast",eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run ncbiblast alone",                                                               0 }, 
  { "--onlyssearch",  eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run sserach alone",                                                                 0 }, 
  { "--onlymuscle",   eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run muscle alone",                                                                  0 }, 
  { "--onlymsaprobs", eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to run MSAProbs alone",                                                                0 }, 
  /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "-I",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < x id",                                                               0 },
  { "--submsa",       eslARG_INT,       FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   0 },
  { "--minid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      0 },
  { "--maxid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",                                                      0 },
  /* options for msa */
  { "--local",        eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  "--e2f,--e2fhmmer", "TRUE for a local alignment, default is global",                                             0 }, 
  { "--viterbi",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  "--e2f,--e2fhmmer", "TRUE for getting the alignment by a viterbi path, default is optimal accuracy",             0 }, 
  { "--e2",           eslARG_NONE,       TRUE,   NULL,       NULL,   ALIOPTS, NULL,  NULL,               "for the e2 algorithm, realign pairs of sequences",                                          0 }, 
  { "--e2f",          eslARG_NONE,      FALSE,   NULL,       NULL,   ALIOPTS, NULL,  NULL,               "for the e2 algorithm, use alignment as is, do not realign",                                 0 }, 
  { "--e2hmmer",      eslARG_STRING,    FALSE,   NULL,       NULL,   ALIOPTS, NULL,  NULL,               "for the e2 algorithm, realign pairs of sequences, use position-specific evomodel",          0 }, 
  { "--e2fhmmer",     eslARG_STRING,    FALSE,   NULL,       NULL,   ALIOPTS, NULL,  NULL,               "for the e2 algorithm, use alignment as is, do not realign, use position-specific evomodel", 0 }, 
  /* options for evolutionary method */
  { "--fixpid",       eslARG_REAL,       NULL,   NULL,"0<=x<=100",   NULL,    NULL,  NULL,               "TRUE for using a fix pid for the ehmm model",                                               0 },
  { "--fixtime",      eslARG_REAL,       NULL,   NULL,     "x>=0",      NULL, NULL,  NULL,               "TRUE for using a fix time for the evolutionary models of a pair",                           0 },
    /* options for optimization */
  { "--optnone",      eslARG_NONE,       TRUE,   NULL,       NULL,   OPTOPTS, NULL,  NULL,               "TRUE for no optimization of parameters or time",                                            0 }, 
  { "--opttime",      eslARG_NONE,      FALSE,   NULL,       NULL,   OPTOPTS, NULL,  NULL,               "TRUE for optimization of time",                                                             0 }, 
  { "--optpara",      eslARG_NONE,      FALSE,   NULL,       NULL,   OPTOPTS, NULL,  NULL,               "TRUE for optimization of parameters",                                                       0 }, 
  { "--optboth",      eslARG_NONE,      FALSE,   NULL,       NULL,   OPTOPTS, NULL,  NULL,               "TRUE for optimization of time and parameters",                                              0 }, 
  /* Control of scoring system - substitutions */
  { "--mx",           eslARG_STRING,"BLOSUM62",  NULL,       NULL,   NULL,    NULL,  "--mxfile",         "substitution rate matrix choice (of some built-in matrices)",                               0 },
  { "--mxfile",       eslARG_INFILE,      NULL,  NULL,       NULL,   NULL,    NULL,  "--mx",             "read substitution rate matrix from file <f>",                                               0 },
  /* Control ncbiblast */
  { "--wordsize",     eslARG_INT,       FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "ncbiblast wordsize, if wordsize < 0, use default",                                          0 },
  /* Control fasta/ssearch */
  { "--fastaversion", eslARG_STRING,    FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "fasta version: fasta-36.3.6d, or ...",                                                      0 },
  { "-f",             eslARG_INT,       FALSE,   NULL,      "n<0",   NULL,    NULL,  NULL,               "gap-open penalty",                                                                          0 },
  { "-g",             eslARG_INT,       FALSE,   NULL,      "n<0",   NULL,    NULL,  NULL,               "gap-extend penalty",                                                                        0 },
  { "-s",             eslARG_STRING,    FALSE,   NULL,       NULL,   NULL, "-f,-g",  NULL,               "Scoring matrix: (protein) BL50, BP62 (sets -f -11 -g -1); P250, OPT5, VT200,VT160, P120, VT120, BL80, VT80, MD40, VT40, MD20, VT20, MD10, VT10; scoring matrix file name; -s ?BL50 adjusts matrix for short queries;",                                          0 },
  /* Control phmmer/hmmsearch */
  { "--max",          eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "use --max option for phmmer or hmmsearch",                                                  0 },
  { "--dostats",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "use --max option for phmmer or hmmsearch",                                                  0 },
  /* Control of scoring system - indels */ 
  { "--evomodel",     eslARG_STRING,    "AFG",   NULL,       NULL,   NULL,    NULL,  NULL,               "evolutionary model used",                                                                   0 },
  { "--paramfile",    eslARG_STRING,     NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "file with rate parameters (overrides the individual options below)",                        0 },
  { "--muAM",         eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of ancestral sequences",                                                   0 },
  { "--muAD",         eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of ancestral sequences",                                                   0 },
  { "--muAI",         eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of ancestral sequences",                                                   0 },
  { "--muEM",         eslARG_REAL,     "0.18",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of inserted residues",                                                     0 },
  { "--ldEM",         eslARG_REAL,     "0.176",   NULL,     "x>=0",   NULL,    NULL,  NULL,              "rate of adding a res to an insert",                                                         0 },
  { "--muED",         eslARG_REAL,     "0.18",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of inserted residues",                                                     0 },
  { "--ldED",         eslARG_REAL,     "0.176",   NULL,     "x>=0",   NULL,    NULL,  NULL,              "rate of adding a res to an insert",                                                         0 },
  { "--muI",          eslARG_REAL,     "0.08",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of deletion of inserted sequences",                                                    0 },
  { "--ldI",          eslARG_REAL,     "0.07",   NULL,     "x>=0",   NULL,    NULL,  NULL,               "rate of adding a whole insert",                                                             0 },
  { "--rI",           eslARG_REAL,     "0.89",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "fragment parameter for I",                                                                  0 },
  { "--rM",           eslARG_REAL,     "0.90",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "fragment parameter for M",                                                                  0 },
  { "--rD",           eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "fragment parameter for D",                                                                  0 },
  { "--userates",     eslARG_NONE,      TRUE,    NULL,       NULL,   NULL,    NULL,  NULL,               "if true start with raw rates, phmmer style otherwise",                                      0 },
  { "--popen",        eslARG_REAL,    "0.020",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "popen in phmmer",                                                                           0 },
  { "--pextend",      eslARG_REAL,    "0.400",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "pextend in phmmer",                                                                         0 },
  { "--pcross",       eslARG_REAL,    "0.001",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,               "pcross, zero in phmm",                                                                      0 },
  /* Control of running options */
  { "--nofastsp",     eslARG_NONE,       TRUE,   NULL,       NULL,   NULL,    NULL,  NULL,               "compare to original alignment using FastSP",                                                0 },
  { "--lcur",         eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE for lower case in reference msa unaligned using FastSP",                               0 },
  { "--lcue",         eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE for lower case in inferred msa unaligned using FastSP",                                0 },
  { "--probdist",     eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "probability distribution landscape",                                                        0 },
  { "--train",        eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "train on msas",                                                                             0 },
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
static char banner[] = "progressive alignment using e2 evolutionary model";

static int create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int run_e2msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int run_phmmer(struct cfg_s *cfg, ESL_MSA *msa, int isphmmer3, int isephmmer, float time);
static int run_ncbiblast(struct cfg_s *cfg, ESL_MSA *msa);
static int run_ssearch(struct cfg_s *cfg, ESL_MSA *msa);
static int run_muscle(struct cfg_s *cfg, ESL_MSA *msa);
static int run_msaprobs(struct cfg_s *cfg, ESL_MSA *msa);
static int run_benchmark(struct cfg_s *cfg, char *method, ESL_MSA *rmsa, MSA_STAT *mrstat, ESL_MSA *emsa, MSA_STAT *mestat, float sc);
static int train(struct cfg_s *cfg, ESL_MSA *msa);
static int explore_distribution(struct cfg_s *cfg, ESL_MSA *msa);

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
  cfg.probdist   = esl_opt_GetBoolean(go, "--probdist"); 
  cfg.train      = esl_opt_GetBoolean(go, "--train"); 

  cfg.mstat = NULL;

  /* branch length options */
  cfg.fixpid     = esl_opt_IsOn(go, "--fixpid")?  esl_opt_GetReal(go, "--fixpid") : -1.0; 
  cfg.fixtime    = esl_opt_IsOn(go, "--fixtime")? esl_opt_GetReal(go, "--fixtime") : -1.0;
 
  cfg.doephmmer   = esl_opt_GetBoolean(go, "--ephmmer");
  cfg.dophmmer    = esl_opt_GetBoolean(go, "--phmmer");
  cfg.dophmmer3   = esl_opt_GetBoolean(go, "--phmmer3");
  cfg.doncbiblast = esl_opt_GetBoolean(go, "--ncbiblast");
  cfg.dossearch   = esl_opt_GetBoolean(go, "--ssearch");
  cfg.domuscle    = esl_opt_GetBoolean(go, "--muscle");
  cfg.domsaprobs  = esl_opt_GetBoolean(go, "--msaprobs");
  cfg.dofastsp    = esl_opt_GetBoolean(go, "--nofastsp");
  cfg.lcu_r       = esl_opt_GetBoolean(go, "--lcur");
  cfg.lcu_e       = esl_opt_GetBoolean(go, "--lcue");
  cfg.doe2msa     = TRUE;
  if (esl_opt_GetBoolean(go, "--onlyephmmer")) {
    cfg.doe2msa     = FALSE;
    cfg.doephmmer   = TRUE;
  }
  if (esl_opt_GetBoolean(go, "--onlyphmmer")) {
    cfg.doe2msa     = FALSE;
    cfg.dophmmer    = TRUE;
  }
  if (esl_opt_GetBoolean(go, "--onlyphmmer3")) {
    cfg.doe2msa     = FALSE;
    cfg.dophmmer3   = TRUE;
  }
  if (esl_opt_GetBoolean(go, "--onlyncbiblast")) {
    cfg.doe2msa     = FALSE;
    cfg.doncbiblast = TRUE;
  }
  if (esl_opt_GetBoolean(go, "--onlyssearch")) {
    cfg.doe2msa     = FALSE;
    cfg.dossearch   = TRUE;
  }
  if (esl_opt_GetBoolean(go, "--onlymuscle")) {
    cfg.doe2msa   = FALSE;
    cfg.domuscle  = TRUE;
  }
   if (esl_opt_GetBoolean(go, "--onlymsaprobs")  ) {
    cfg.doe2msa    = FALSE;
    cfg.domsaprobs = TRUE;
  }
 
   /* ncbiblast options */
  if (esl_opt_IsOn(go, "--wordsize"))      cfg.wordsize = esl_opt_GetInteger(go, "--wordsize");
  else                                     cfg.wordsize = -1;

  /* fasta/ssearch options */
  if (esl_opt_IsOn(go, "-f"))              cfg.fasta_f = esl_opt_GetInteger(go, "-f");
  else                                     cfg.fasta_f = 0;
  if (esl_opt_IsOn(go, "-g"))              cfg.fasta_g = esl_opt_GetInteger(go, "-g");
  else                                     cfg.fasta_g = 0;
  if (esl_opt_IsOn(go, "-s"))              esl_sprintf(&cfg.fasta_s, esl_opt_GetString(go, "-s"));
  else                                     esl_sprintf(&cfg.fasta_s, "default");
  if (esl_opt_IsOn(go, "--fastaversion"))  esl_sprintf(&cfg.fastaversion, esl_opt_GetString(go, "--fastaversion"));
  else                                     esl_sprintf(&cfg.fastaversion, "fasta-36.3.6d");

  /* phmmer */
  if (esl_opt_IsOn(go, "--max"))           cfg.max = TRUE;
  else                                     cfg.max = FALSE;
  if (esl_opt_IsOn(go, "--dostats"))       cfg.dostats = TRUE;
  else                                     cfg.dostats = FALSE;

  cfg.subsmx     = esl_opt_GetString(go, "--mx");
  cfg.msafrq     = NULL;
  cfg.mode       = esl_opt_IsOn(go, "--local")? e2_LOCAL : e2_GLOBAL;
  cfg.do_viterbi = esl_opt_IsOn(go, "--viterbi")? TRUE : FALSE;

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
  esl_sprintf(&cfg.method, "e2msa.%s", e1_rate_EvomodelType(cfg.evomodel)); 

  /* the paramfile takes precedent over individually feed values above */
  cfg.paramfile = NULL;
  if (esl_opt_IsOn(go, "--paramfile")) {
    cfg.paramfile = esl_opt_GetString(go, "--paramfile");
    status = e1_rate_ReadParamfile(cfg.paramfile, &cfg.rateparam, &cfg.evomodel, cfg.errbuf, cfg.verbose);
    if (status != eslOK) esl_fatal("Failed to read paramfile %s\n%s", cfg.paramfile, cfg.errbuf);
  }

  cfg.userates       = esl_opt_GetBoolean(go, "--userates");
  cfg.popen          = esl_opt_GetReal   (go, "--popen");
  cfg.pextend        = esl_opt_GetReal   (go, "--pextend");
  cfg.pcross         = esl_opt_GetReal   (go, "--pcross");
  
  if      (esl_opt_GetBoolean(go, "--e2"))       cfg.e2ali = E2;
  else if (esl_opt_GetBoolean(go, "--e2f"))      cfg.e2ali = E2F;
  else if (esl_opt_GetString (go, "--e2hmmer"))  cfg.e2ali = E2HMMER;
  else if (esl_opt_GetString (go, "--e2fhmmer")) cfg.e2ali = E2FHMMER;
  else                                           cfg.e2ali = E2NONE;
 
  if      (esl_opt_GetBoolean(go, "--optnone")) cfg.e2optimize = OPTNONE;
  else if (esl_opt_GetBoolean(go, "--opttime")) cfg.e2optimize = OPTTIME;
  else if (esl_opt_GetBoolean(go, "--optpara")) cfg.e2optimize = OPTPARA;
  else if (esl_opt_GetBoolean(go, "--optboth")) cfg.e2optimize = OPTBOTH;

  if (cfg.e2ali != E2 && cfg.mode == e2_LOCAL) {
    esl_fatal("local alignment not implemented yet\n");
  }

  /* open outfile with benchmarking results */
  cfg.benchfile = NULL;
  cfg.benchfp   = NULL;
  esl_sprintf(&cfg.benchfile, "%s.%s.bench", cfg.outheader, e1_rate_EvomodelType(cfg.evomodel)); 
  if ((cfg.benchfp = fopen(cfg.benchfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.benchfile);
  fprintf(cfg.benchfp, "#MSA_NAME           TP_hom   TRUE_hom  FOUND_hom      \n");

  cfg.msarfname           = NULL;
  cfg.msaefname           = NULL;
  cfg.msaefname_ephmmer   = NULL;
  cfg.msaefname_phmmer    = NULL;
  cfg.msaefname_phmmer3   = NULL;
  cfg.msaefname_ncbiblast = NULL;
  cfg.msaefname_ssearch   = NULL;
  cfg.msaefname_muscle    = NULL;
  cfg.msaefname_msaprobs  = NULL;
  
  cfg.T = NULL;

  /* if e2hmmer or e2fhmmer, read the input HMM */
  if (cfg.e2ali == E2FHMMER &&
      (cfg.hmmfile = esl_opt_GetString(go, "--e2fhmmer")) == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (cfg.e2ali == E2HMMER &&
      (cfg.hmmfile = esl_opt_GetString(go, "--e2hmmer"))  == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  cfg.R1 = NULL;
  cfg.R7 = NULL;
  cfg.bg  = NULL;
  cfg.bg7 = NULL;
  cfg.bg  = e1_bg_Create(cfg.abc);
  cfg.bg7 = p7_bg_Create(cfg.abc);
 
  cfg.hfp = NULL;
  cfg.hmm = NULL;
  cfg.emR = NULL;

  cfg.EmL     = 200;
  cfg.EmN     = 200;
  cfg.EvL     = 200;
  cfg.EvN     = 200;
  cfg.EfL     = 100;
  cfg.EfN     = 200;
  cfg.Eft     = 0.04;

  cfg.popen   = 0.04;
  cfg.pextend = 0.2;

  esl_sprintf(&cfg.mx, "BLOSUM62");

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
      if (msamanip_SelectSubset(cfg.r, cfg.submsa, &msa, NULL, cfg.errbuf, cfg.verbose) != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
    }
    
    /* outheader for all msa-output files */
    msamanip_OutfileHeader((msa->acc)?msa->acc:cfg.msafile, &cfg.msaheader); 
   
    esl_msa_Hash(msa);
   
    if (esl_opt_IsOn(go, "-F") && msamanip_RemoveFragments(cfg.fragfrac, &msa, &nfrags, &seq_cons_len) != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-I"))   msamanip_SelectSubsetBymaxID(cfg.r, &msa, cfg.idthresh, &nremoved);
    
    /* given msa aveid and avematch */
    msamanip_CStats(cfg.abc, msa, &cfg.mstat);
    msamanip_CBaseComp(cfg.abc, msa, cfg.bg->f, &cfg.msafrq);

    if (esl_opt_IsOn(go, "--minid") && cfg.mstat->avgid < 100.*esl_opt_GetReal(go, "--minid")) continue;
    if (esl_opt_IsOn(go, "--maxid") && cfg.mstat->avgid > 100.*esl_opt_GetReal(go, "--maxid")) continue;

    /* print some info */
    if (cfg.voutput) {
      esl_sprintf(&cfg.msarfname, "%s.%s", cfg.msaheader, e1_rate_EvomodelType(cfg.evomodel)); 
      fprintf(cfg.outfp, "Given alignment\n");
      fprintf(cfg.outfp, "%6d          %s\n", msa->nseq, cfg.msafile);
      if (eslx_msafile_Write(cfg.outfp, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write to file %s", cfg.msarfname); 
      msamanip_DumpStats(cfg.outfp, msa, cfg.mstat); 
    }
    
    if (cfg.train) {
      cfg.doe2msa     = FALSE;
      cfg.doephmmer   = FALSE;
      cfg.dophmmer    = FALSE;
      cfg.dophmmer3   = FALSE;
      cfg.doncbiblast = FALSE;
      cfg.dossearch   = FALSE;
      cfg.domuscle    = FALSE;
      cfg.domsaprobs  = FALSE;
      status = train(&cfg, msa);
      if (status != eslOK)  { esl_fatal(cfg.errbuf); }
      continue;
    }

    /* the main function */
    if (cfg.doe2msa) {
      status = run_e2msa(go, &cfg, msa);
      if (status != eslOK) esl_fatal("%s Failed to run e2msa", cfg.errbuf);
    }
    
    /* run ephmmer */
    if (cfg.doephmmer) {
      status = run_phmmer(&cfg, msa, FALSE, TRUE, cfg.fixtime);
      if (status != eslOK) esl_fatal("%s Failed to run ephmmer", cfg.errbuf);
    }

    /* run phmmer */
    if (cfg.dophmmer) {
      status = run_phmmer(&cfg, msa, FALSE, FALSE, -1.0);
      if (status != eslOK) esl_fatal("%s Failed to run phmmer", cfg.errbuf);
    }

    /* run phmmer3 */
    if (cfg.dophmmer3) {
      status = run_phmmer(&cfg, msa, TRUE, FALSE, -1.0);
      if (status != eslOK) esl_fatal("%s Failed to run phmmer3", cfg.errbuf);
    }
    
    /* run ncbiblast */
    if (cfg.doncbiblast) {
      status = run_ncbiblast(&cfg, msa);
      if (status != eslOK) esl_fatal("%s Failed to run ncbiblast", cfg.errbuf);
    }
    
    /* run ssearch */
    if (cfg.dossearch) {
      status = run_ssearch(&cfg, msa);
      if (status != eslOK) esl_fatal("%s Failed to run ssearch", cfg.errbuf);
    }
    
    /* run muscle */
    if (cfg.domuscle) {
      status = run_muscle(&cfg, msa);
      if (status != eslOK) esl_fatal("%s Failed to run muscle", cfg.errbuf);
    }
    
    /* run msaprobs */
    if (cfg.domsaprobs) {
      status = run_msaprobs(&cfg, msa);
      if (status != eslOK) esl_fatal("%s\nFailed to run msaprobs", cfg.errbuf);
    }
    
    esl_msa_Destroy(msa); msa = NULL;
    free(cfg.msafrq); cfg.msafrq = NULL;
  }

  /* cleanup */
  esl_stopwatch_Destroy(cfg.w);
  if (cfg.bg)  e1_bg_Destroy(cfg.bg);     
  if (cfg.bg7) p7_bg_Destroy(cfg.bg7);     
  fclose(cfg.outfp);
  free(cfg.outheader);
  fclose(cfg.benchfp);
  free(cfg.benchfile);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  eslx_msafile_Close(afp);
  free(cfg.mstat);
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
run_e2msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  ESL_MSA   *dmsa  = NULL;
  ESL_MSA   *e2msa = NULL;
  MSA_STAT  *alle2mstat = NULL;
  MSA_STAT  *e2mstat = NULL;
  float      sc;
  int        nnodes;
  int        r;
  int        status;

  esl_stopwatch_Start(cfg->w);

  dmsa = esl_msa_Clone(msa);
  esl_msa_Digitize(cfg->abc, dmsa, NULL);

  /* Create the evolutionary rate model */
  cfg->nr    = 1;
  cfg->R1    = malloc(sizeof(E1_RATE) * cfg->nr);
  cfg->R1[0] = NULL;
  if (cfg->e2ali == E2 || cfg->e2ali == E2F) {
    for (r = 0; r < cfg->nr; r ++) {
      if (cfg->userates) cfg->R1[r] = e1_rate_CreateWithValues(cfg->abc, cfg->evomodel, cfg->rateparam, cfg->subsmx, NULL, TRUE, cfg->tol, cfg->errbuf, cfg->verbose);
      else               cfg->R1[r] = e1_rate_CreateFromCosts(cfg->abc, cfg->evomodel, cfg->popen, cfg->pextend, cfg->pcross, cfg->subsmx, 
							      NULL, TRUE, cfg->tol, cfg->errbuf, cfg->verbose);
      if (cfg->R1[r] == NULL) { printf("Bad rate model.\n"); esl_fatal(cfg->errbuf); }
      if (1||cfg->verbose)   printf("rI %f rM %f rD %f | ldEM %f muEM %f ldED %f muED %f ldI %f muI %f muA %f %f %f| model %s\n",  
				    cfg->R1[r]->rI, cfg->R1[r]->rM, cfg->R1[r]->rD, 
				    cfg->R1[r]->ldE[e1R_S], cfg->R1[r]->muE[e1R_S], cfg->R1[r]->ldE[e1R_D], cfg->R1[r]->muE[e1R_D], 
				    cfg->R1[r]->ldE[e1R_I], cfg->R1[r]->muE[e1R_I], 
				    cfg->R1[r]->muA[e1R_S], cfg->R1[r]->muA[e1R_D], cfg->R1[r]->muA[e1R_I], 
				    e1_rate_EvomodelType(cfg->R1[r]->evomodel));
    }
  }
  else {  
    /* Open the query profile HMM file */
    status = p7_hmmfile_OpenE(cfg->hmmfile, NULL, &cfg->hfp, cfg->errbuf);
    if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, cfg->errbuf);
    else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, cfg->errbuf);
    else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, cfg->errbuf);  
    
    /* <abc> is not known 'til  HMM is read. */
    status = p7_hmmfile_Read(cfg->hfp, &(cfg->abc), &cfg->hmm);
    if (status != eslOK) p7_Fail("error reading HMM file %s.\n%s\n",  cfg->hmmfile, cfg->errbuf);  
    
    /* Calculate the hmm rate */
    cfg->emR = ratematrix_emrate_Create(cfg->abc, 1);
    ratematrix_emrate_Set(cfg->subsmx, NULL, cfg->bg7->f, cfg->emR, TRUE, cfg->tol, cfg->errbuf, FALSE);

    int betainf;
    int etainf;
    betainf = (cfg->rateparam.muEM > 0.)? ((cfg->rateparam.ldEM < cfg->rateparam.muEM)? cfg->rateparam.ldEM/cfg->rateparam.muEM : 1.0 ) : -1.;
    etainf  = (cfg->rateparam.muI  > 0.)? ((cfg->rateparam.ldI  < cfg->rateparam.muI)?  cfg->rateparam.ldI/cfg->rateparam.muI : 1.0 ) : -1.;
    if (cfg->evomodel != GG) cfg->rateparam.rI = -1.0; /* does not have to be provided */
    if (cfg->evomodel != GG) etainf   = -1.0; /* does not have to be provided */
    if (p7_RateCalculate(stdout, cfg->hmm, cfg->bg7, cfg->emR, NULL, &(cfg->R7), cfg->evomodel, betainf, cfg->fixtime, 0.001, cfg->errbuf, FALSE) != eslOK)  
      esl_fatal("%s", cfg->errbuf);
  }

  /* Initializations */
  cfg->pli = e2_pipeline_Create(go, (cfg->hmm)? cfg->hmm->M:1, 100, 100);
  
  if (cfg->probdist) {
    cfg->doephmmer   = FALSE;
    cfg->dophmmer    = FALSE;
    cfg->dophmmer3   = FALSE;
    cfg->doncbiblast = FALSE;
    cfg->dossearch   = FALSE;
    cfg->domuscle    = FALSE;
    cfg->domsaprobs  = FALSE;
    status = explore_distribution(cfg, msa);
    if (status != eslOK)  { esl_fatal(cfg->errbuf); }
    return eslOK;
  }
  
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

  /* The leaves-only alignment */
  msamanip_MSALeaves(&e2msa, FALSE);   
  msamanip_CStats(cfg->abc, e2msa, &e2mstat);
  e2mstat->anclen = alle2mstat->anclen; // transfer the len of the ancestral sequence (cannot be calculated from the leaves-only msa)

  if (cfg->voutput) { /* write output alignment */
    /* print output info */
    fprintf(cfg->outfp, "\ne2msa-%s alignment\n", e1_rate_EvomodelType(cfg->evomodel));
    eslx_msafile_Write(cfg->outfp, e2msa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, e2msa, e2mstat);
  }

  if (cfg->voutput) { /* the output alignment to an independent file */
    esl_sprintf(&cfg->msaefname, "%s.e2msa.%s", cfg->msaheader, e1_rate_EvomodelType(cfg->evomodel)); 
    if ((cfg->msaefp = fopen(cfg->msaefname, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg->msaefname);
    if (eslx_msafile_Write(cfg->msaefp, e2msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write to file %s", cfg->msaefname); 
    fclose (cfg->msaefp);
  }
  
  if (cfg->dofastsp) {
    if ((status = run_benchmark(cfg, cfg->method, msa, cfg->mstat, e2msa, e2mstat, sc)) != eslOK) goto ERROR;
  }
  
  /* cleanup */
  if (cfg->R1) {
    for (r = 0; r < cfg->nr; r ++) e1_rate_Destroy(cfg->R1[r]);
    free(cfg->R1);
  }
  if (cfg->R7) p7_RateDestroy(cfg->R7);
  if (cfg->msarfname)          free(cfg->msarfname);
  if (cfg->msaefname)          free(cfg->msaefname);
  
  e2_pipeline_Destroy(cfg->pli);
  if (cfg->hmm) p7_hmm_Destroy(cfg->hmm); 
  if (cfg->hfp) p7_hmmfile_Close(cfg->hfp);
  if (cfg->emR) ratematrix_emrate_Destroy(cfg->emR, 1);
 
  if (dmsa)   esl_msa_Destroy(dmsa);
  if (e2msa)  esl_msa_Destroy(e2msa);
  if (e2mstat) free(e2mstat);
  if (alle2mstat) free(alle2mstat);
  if (cfg->T) esl_tree_Destroy(cfg->T);
  
  return eslOK;

 ERROR:
 if (dmsa)  esl_msa_Destroy(dmsa);
 if (e2msa) esl_msa_Destroy(e2msa);
 return status;
}

static int
run_phmmer(struct cfg_s *cfg, ESL_MSA *msa, int isphmmer3, int isephmmer, float time)
{
  char     *name = NULL;
  ESL_MSA  *phmmermsa = NULL;
  MSA_STAT *phmmermstat = NULL;
  float     usetime = -1.0;
  float     sc = -eslINFINITY;
  int       status;

  esl_stopwatch_Start(cfg->w);
  if (isephmmer) EPHMMER_Align(msa, time, cfg->fixpid,       &phmmermsa, &usetime, &sc, cfg->bg7, cfg->EmL, cfg->EmN, cfg->EvL, cfg->EvN, cfg->EfL, cfg->EfN, cfg->Eft, cfg->popen, cfg->pextend, cfg->mx, cfg->max, cfg->dostats, cfg->tol,  cfg->errbuf, cfg->verbose);
  else           PHMMER_Align (msa, (isphmmer3)? TRUE:FALSE, &phmmermsa,           &sc, cfg->bg7, cfg->EmL, cfg->EmN, cfg->EvL, cfg->EvN, cfg->EfL, cfg->EfN, cfg->Eft, cfg->popen, cfg->pextend, cfg->mx, cfg->max, cfg->dostats,            cfg->errbuf, cfg->verbose);
  esl_stopwatch_Stop(cfg->w);
  cfg->treeavgt = usetime;
 
  msamanip_CStats(cfg->abc, phmmermsa, &phmmermstat);
  if (phmmermsa && (cfg->voutput)) {/* the PHMMER output alignment */
    fprintf(cfg->outfp,"\nPHMMER alignment\n");
    eslx_msafile_Write(cfg->outfp, phmmermsa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, phmmermsa, phmmermstat);
  }
    
  if (cfg->voutput && phmmermsa) {/* write the output alignment to an independent file */
    if (isephmmer) {
      esl_sprintf(&cfg->msaefname_ephmmer,  "%s.ephmmer", cfg->msaheader); 
      if ((cfg->msaefp_phmmer = fopen(cfg->msaefname_ephmmer, "w")) == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to open output file %s", cfg->msaefname_ephmmer);
      if (eslx_msafile_Write(cfg->msaefp_ephmmer, phmmermsa, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to write to file %s", cfg->msaefname_ephmmer);
      fclose (cfg->msaefp_ephmmer);
    }   
    else if (isphmmer3) {
      esl_sprintf(&cfg->msaefname_phmmer3, "%s.phmmer3", cfg->msaheader); 
      if ((cfg->msaefp_phmmer3 = fopen(cfg->msaefname_phmmer3, "w")) == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to open output file %s", cfg->msaefname_phmmer3);
      if (eslx_msafile_Write(cfg->msaefp_phmmer3, phmmermsa, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to write to file %s", cfg->msaefname_phmmer3);
      fclose (cfg->msaefp_phmmer3);
    }
    else  {
      esl_sprintf(&cfg->msaefname_phmmer,  "%s.phmmer", cfg->msaheader); 
      if ((cfg->msaefp_phmmer = fopen(cfg->msaefname_phmmer, "w")) == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to open output file %s", cfg->msaefname_phmmer);
      if (eslx_msafile_Write(cfg->msaefp_phmmer, phmmermsa, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to write to file %s", cfg->msaefname_phmmer);
      fclose (cfg->msaefp_phmmer);
    }
  }  
  
  if (cfg->dofastsp) {
    if      (isephmmer) esl_sprintf(&name, "ePHMMER");
    else if (isphmmer3) esl_sprintf(&name, "PHMMER3");
    else                esl_sprintf(&name, "PHMMER");
    if (run_benchmark(cfg, name,  msa, cfg->mstat, phmmermsa, phmmermstat, sc) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run FastSP");
  }

  if (name) free(name);
  if (phmmermsa) esl_msa_Destroy(phmmermsa);
  if (phmmermstat) free(phmmermstat); phmmermstat = NULL;
  if (cfg->msaefname_phmmer3) free(cfg->msaefname_phmmer3); cfg->msaefname_phmmer3 = NULL;
  if (cfg->msaefname_phmmer)  free(cfg->msaefname_phmmer);  cfg->msaefname_phmmer  = NULL;
  return eslOK;
  
 ERROR:
  if (name) free(name);
  if (phmmermsa) esl_msa_Destroy(phmmermsa);
  if (phmmermstat) free(phmmermstat); phmmermstat = NULL;
  if (cfg->msaefname_phmmer3) free(cfg->msaefname_phmmer3); cfg->msaefname_phmmer3 = NULL;
  if (cfg->msaefname_phmmer)  free(cfg->msaefname_phmmer);  cfg->msaefname_phmmer  = NULL;
  return status;
}

static int
run_ncbiblast(struct cfg_s *cfg, ESL_MSA *msa)
{
  ESL_MSA  *blastmsa = NULL;
  MSA_STAT *blastmstat = NULL;
  float     sc = -eslINFINITY;
  int       status;

  esl_stopwatch_Start(cfg->w);
  NCBIBLAST_Align(msa, cfg->wordsize, &blastmsa, cfg->errbuf, cfg->verbose);
  esl_stopwatch_Stop(cfg->w);
  
  msamanip_CStats(cfg->abc, blastmsa, &blastmstat);
  if (blastmsa && (cfg->voutput)) {/* the BLAST output alignment */
    fprintf(cfg->outfp,"\nBLAST alignment\n");
    eslx_msafile_Write(cfg->outfp, blastmsa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, blastmsa, blastmstat);
  }
    
  if (cfg->voutput && blastmsa) {/* write the output alignment to an independent file */
    esl_sprintf(&cfg->msaefname_ncbiblast, "%s.ncbiblast", cfg->msaheader); 
    if ((cfg->msaefp_ncbiblast = fopen(cfg->msaefname_ncbiblast, "w")) == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to open output file %s", cfg->msaefname_ncbiblast);
    if (eslx_msafile_Write(cfg->msaefp_ncbiblast, blastmsa, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to write to file %s", cfg->msaefname_ncbiblast);
    fclose (cfg->msaefp_ncbiblast);
  }  
  
  if (cfg->dofastsp) {
    if (run_benchmark(cfg, "NCBIBLAST", msa, cfg->mstat, blastmsa, blastmstat, sc) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run FastSP");
  }

  if (blastmsa) esl_msa_Destroy(blastmsa);
  if (blastmstat) free(blastmstat);
  if (cfg->msaefname_ncbiblast) free(cfg->msaefname_ncbiblast);
   return eslOK;

 ERROR:
  if (blastmsa) esl_msa_Destroy(blastmsa);
  if (blastmstat) free(blastmstat);
  if (cfg->msaefname_ncbiblast) free(cfg->msaefname_ncbiblast);
  return status;
}

static int
run_ssearch(struct cfg_s *cfg, ESL_MSA *msa)
{
  ESL_MSA  *ssearchmsa = NULL;
  MSA_STAT *ssearchmstat = NULL;			
  float     sc;
  int       status;

  esl_stopwatch_Start(cfg->w);
  SSEARCH_Align(cfg->fastaversion, msa, &ssearchmsa, &sc, cfg->fasta_s, cfg->fasta_f, cfg->fasta_g, cfg->errbuf, cfg->verbose);
  esl_stopwatch_Stop(cfg->w);
  
  msamanip_CStats(cfg->abc, ssearchmsa, &ssearchmstat);
  if (ssearchmsa && (cfg->voutput)) {/* the SSEARCH output alignment */
    fprintf(cfg->outfp,"\nSSEARCH alignment\n");
    eslx_msafile_Write(cfg->outfp, ssearchmsa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, ssearchmsa, ssearchmstat);
  }
    
  if (cfg->voutput && ssearchmsa) {/* write the output alignment to an independent file */
    esl_sprintf(&cfg->msaefname_ssearch, "%s.ssearch", cfg->msaheader); 
    if ((cfg->msaefp_ssearch = fopen(cfg->msaefname_ssearch, "w")) == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to open output file %s", cfg->msaefname_ssearch);
    if (eslx_msafile_Write(cfg->msaefp_ssearch, ssearchmsa, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to write to file %s", cfg->msaefname_ssearch);
    fclose (cfg->msaefp_ssearch);
  }  
  
  if (cfg->dofastsp) {
    if (run_benchmark(cfg, "SSEARCH", msa, cfg->mstat, ssearchmsa, ssearchmstat, sc) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run FastSP");
  }

  if (ssearchmsa) esl_msa_Destroy(ssearchmsa);
  if (ssearchmstat) free(ssearchmstat);
  if (cfg->msaefname_ssearch) free(cfg->msaefname_ssearch);
   return eslOK;

 ERROR:
  if (ssearchmsa) esl_msa_Destroy(ssearchmsa);
  if (cfg->msaefname_ssearch) free(cfg->msaefname_ssearch);
  return status;
}

static int
run_muscle(struct cfg_s *cfg, ESL_MSA *msa)
{
  ESL_MSA  *musclemsa = NULL;
  MSA_STAT *musclemstat = NULL;
  float     sc = -eslINFINITY;
  int       status;

  esl_stopwatch_Start(cfg->w);
  MUSCLE_Align(msa, &musclemsa, cfg->errbuf, cfg->verbose);
  esl_stopwatch_Stop(cfg->w);
  
  msamanip_CStats(cfg->abc, musclemsa, &musclemstat);
  if (cfg->voutput) {/* the MUSCLE output alignment */
    fprintf(cfg->outfp,"\nMUSCLE alignment\n");
    eslx_msafile_Write(cfg->outfp, musclemsa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, musclemsa, musclemstat);
  }
  
  if (cfg->voutput) {/* write the output alignment to an independent file */
    esl_sprintf(&cfg->msaefname_muscle, "%s.muscle", cfg->msaheader); 
    if ((cfg->msaefp_muscle = fopen(cfg->msaefname_muscle, "w")) == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to open output file %s", cfg->msaefname_muscle);
    if (eslx_msafile_Write(cfg->msaefp_muscle, musclemsa, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to write to file %s", cfg->msaefname_muscle);
    fclose (cfg->msaefp_muscle);
  }  

  if (cfg->dofastsp) {
    if (run_benchmark(cfg, "MUSCLE", msa, cfg->mstat, musclemsa, musclemstat, sc) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run FastSP");
  }

  if (musclemsa)             esl_msa_Destroy(musclemsa);
  if (musclemstat)           free(musclemstat);
  if (cfg->msaefname_muscle) free(cfg->msaefname_muscle);
   return eslOK;

 ERROR:
  if (musclemsa)             esl_msa_Destroy(musclemsa);
  if (musclemstat)           free(musclemstat);
  if (cfg->msaefname_muscle) free(cfg->msaefname_muscle);
  return status;
}

static int
run_msaprobs(struct cfg_s *cfg, ESL_MSA *msa)
{
  ESL_MSA  *msaprobsmsa = NULL;  
  MSA_STAT *msaprobsmstat = NULL;
  float     sc = -eslINFINITY;
  int       status;
  
  esl_stopwatch_Start(cfg->w);
  status = MSAProbs_Align(msa, &msaprobsmsa, cfg->errbuf, cfg->verbose);
  esl_stopwatch_Stop(cfg->w);
  if (status != eslOK) goto ERROR;

  msamanip_CStats(cfg->abc, msaprobsmsa, &msaprobsmstat);
  
  if (cfg->voutput) {/* the MSAProbs output alignment */
    fprintf(cfg->outfp, "\nMSAProbs alignment\n");
    eslx_msafile_Write(cfg->outfp, msaprobsmsa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, msaprobsmsa, msaprobsmstat);
  }
  
  if (cfg->voutput) {/* write the output alignment to an independent file */
    esl_sprintf(&cfg->msaefname_msaprobs, "%s.msaprobs", cfg->msaheader); 
    if ((cfg->msaefp_msaprobs = fopen(cfg->msaefname_msaprobs, "w")) == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to open output file %s", cfg->msaefname_msaprobs);
    if (eslx_msafile_Write(cfg->msaefp_msaprobs, msaprobsmsa, eslMSAFILE_STOCKHOLM) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to write to file %s", cfg->msaefname_msaprobs);
    fclose (cfg->msaefp_msaprobs);
  }
  
  if (cfg->dofastsp) {
    if (run_benchmark(cfg, "MSAProbs", msa, cfg->mstat, msaprobsmsa, msaprobsmstat, sc) != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run FastSP");
  }

 if (msaprobsmsa)             esl_msa_Destroy(msaprobsmsa);
 if (msaprobsmstat)           free(msaprobsmstat);
 if (cfg->msaefname_msaprobs) free(cfg->msaefname_msaprobs);
 return eslOK;

 ERROR:
 if (msaprobsmsa)             esl_msa_Destroy(msaprobsmsa);
 if (msaprobsmstat)           free(msaprobsmstat);
 if (cfg->msaefname_msaprobs) free(cfg->msaefname_msaprobs);
 return status;
}


static int 
run_benchmark(struct cfg_s *cfg, char *method, ESL_MSA *rmsa, MSA_STAT *mrstat, ESL_MSA *emsa, MSA_STAT *mestat, float sc)
{
  char   *msaname = NULL;
  double  time;
  int     status; 
  
  time = cfg->w->user + cfg->w->sys;
  
  esl_sprintf(&msaname, "%s.%d", cfg->msaheader, cfg->nmsa);

  status = FastSP_Benchmark(cfg->benchfp, (rmsa->name)? rmsa->name : msaname, method, cfg->abc, rmsa, mrstat, emsa, mestat, sc, cfg->treeavgt, time, cfg->lcu_r, cfg->lcu_e, cfg->errbuf, cfg->verbose);
  if (status != eslOK) goto ERROR;

  if (msaname) free(msaname);
  return eslOK;

 ERROR:
  if (msaname) free(msaname);
  return status;
}

static int
train(struct cfg_s *cfg, ESL_MSA *msa)
{
  ESL_MSA    *pairmsa = NULL;
  MSA_STAT   *pairmstat = NULL;
  int        *useme = NULL;
  PSQ       **psq = NULL;
  int         N = msa->nseq; 
  int         i, j; 
  int         status; 
  
  ESL_ALLOC(useme, sizeof(int) * N); 
  for (i = 0; i < N; i++) useme[i] = FALSE; 
  
  ESL_ALLOC(psq, sizeof(PSQ *) * N); 
  for (i = 0; i < N; i++) { 
    psq[i] = psq_CreateFrom(msa->sqname[i], msa->desc, msa->acc, msa->abc, msa->ax[i], msa->alen); 
  }
  
  for (i = 0; i < N; i++) { 
    useme[i] = TRUE; 
    for (j = i+1; j < N; j++) 
      { 
   	useme[j] = TRUE; 
    	esl_msa_SequenceSubset(msa, useme, &pairmsa); 
   	msamanip_XStats(pairmsa, &pairmstat); 
	printf("# id %f match %f alen %" PRId64 "  inserts %d ilen %d\n", pairmstat->avgid, pairmstat->avgmatch, pairmsa->alen, pairmstat->totinum, pairmstat->totilen);

	useme[j] = FALSE;
	esl_msa_Destroy(pairmsa); pairmsa = NULL;
      }
    useme[i] = FALSE;
  }
  
  if (psq) {
    for (i = 0; i < N; i++) 
      psq_Destroy(psq[i]);
    free(psq);
  }
  
  return eslOK;

 ERROR:
  if (useme)   free(useme);  
  if (pairmsa) esl_msa_Destroy(pairmsa); 
  if (pairmstat) free(pairmstat);
   if (psq) {
    for (i = 0; i < N; i++) psq_Destroy(psq[i]);
    free(psq);
  }
  return status;

}

static int
explore_distribution(struct cfg_s *cfg, ESL_MSA *msa)
{
  ESL_MSA    *pairmsa = NULL;
  MSA_STAT   *pairmstat = NULL;
  int        *useme = NULL;
  PSQ       **psq = NULL;
  ESL_SQ     *sq = NULL;
  float       sc, accsc;
  double      timel, timer;
  double      optime;
  double      time;
  double      tmax = 1e+4; 
  double      tinc = 0.01; 
  double      delta; 
  int         N = msa->nseq; 
  int         i, j; 
  int         status; 
  
   ESL_ALLOC(useme, sizeof(int) * N); 
   for (i = 0; i < N; i++) useme[i] = FALSE; 

   ESL_ALLOC(psq, sizeof(PSQ *) * N); 
   for (i = 0; i < N; i++) { 
     status = esl_sq_FetchFromMSA(msa, i, &sq); /* extract the seqs from the msa */ 

     switch(cfg->e2ali) { 
     case E2: 
       psq[i] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);	 
       break; 
     case E2HMMER: 
       psq[i] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n);	 
       break; 
     case E2F: 
       psq[i] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, msa->ax[i], msa->alen); 
       break; 
     case E2FHMMER: 
       psq[i] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, sq->dsq, sq->n); 
       break; 
     default: 
       ESL_XFAIL(eslFAIL, cfg->errbuf, "unknown alicase\n"); goto ERROR; 
     }	 
    
     esl_sq_Destroy(sq); sq = NULL; 
   }	 
  
   for (i = 0; i < N; i++) { 
     useme[i] =  TRUE; 
     for (j = i+1; j < N; j++) 
       { 
   	useme[j] =  TRUE; 
	
    	esl_msa_SequenceSubset(msa, useme, &pairmsa); 
   	msamanip_XStats(pairmsa, &pairmstat); 
	
   	/* optimal P(s1,s2|t) */ 
   	timel = timer = 0.5; 
   	status = e2_Optimize(cfg->r, cfg->pli, psq[i], psq[j], cfg->msafrq, cfg->R1[0], cfg->R7, cfg->bg, cfg->bg7, &timel, &timer, NULL, NULL, &sc, cfg->e2ali, OPTTIME, 
   			     cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose); 
   	optime = timel+timer; 
   	printf("# id %f match %f L %" PRId64 " optime %f optsc %f\n", pairmstat->avgid, pairmstat->avgmatch, pairmsa->alen, optime, sc); 

   	/* P(s1,s2|t) as a function of time */ 
   	time = 1e-6; 
	while (time < tmax) {
	  e2_Pipeline(cfg->r, cfg->pli, psq[i], psq[j], cfg->msafrq, cfg->R1[0], cfg->R7, cfg->bg, cfg->bg7, 0.5*(float)time, 0.5*(float)time, NULL, NULL, &sc, &accsc, cfg->e2ali, cfg->mode, 
		      TRUE, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
	  printf("%f %f\n", time, sc);
	  delta = (time < 0.5)? tinc : ( (time < 3.0)? tinc * 10 : ( (time < 10.0)? tinc * 100 : tinc * 1000 ));
	  time += delta;
	}
	printf("&\n");
	useme[j] = FALSE;
	
	esl_msa_Destroy(pairmsa); pairmsa = NULL; 
	if (pairmstat) free(pairmstat); pairmstat = NULL;
      }
    useme[i] = FALSE;
  }
  
  if (psq) {
    for (i = 0; i < N; i++) 
      psq_Destroy(psq[i]);
    free(psq);
  }
  
  return eslOK;

 ERROR:
  if (useme)   free(useme);  
  if (pairmsa) esl_msa_Destroy(pairmsa);
  if (pairmstat) free(pairmstat);
  if (sq)      esl_sq_Destroy(sq);
  if (psq) {
    for (i = 0; i < N; i++) psq_Destroy(psq[i]);
    free(psq);
  }
  return status;
}


