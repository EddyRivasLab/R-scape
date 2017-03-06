/* triomsa -- align (b,c) by aligning (a,b) and (a,c)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esl_composition.h"
#include "esl_getopts.h"
#include "esl_distance.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2.h"
#include "e2_profilesq.h"
#include "e2_msa.h"
#include "e2_tree.h"
#include "evohmmer.h"
#include "msatree.h"
#include "msaprobs.h"
#include "msamanip.h"
#include "muscle.h"
#include "fastsp.h"

#define ALPHOPTS "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
#define ALIOPTS  "--e2,--e2f,--e2hmmer,--e2fhmmer"          /* Exclusive options for alignment choice */
#define OPTOPTS  "--optnone,--opttime,--optpara,--optboth"   /* Exclusive options for msa optimization */
#define EVOMOPTS "--LI,--LR,--AFG,--AFGR,--AFR,--AIF,-GG,--AG,--AGA,--TKF91,--TKF92,--FID,--GTKF92,--GRTKF92,--ITKF92" /* Exclusive options for evolutionary model choice */
#define SHUF_OPTS "--mono,--di,--markov0,--markov1,--reverse"   /* toggle group, seq shuffling options          */

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

  char            *benchfile1;
  char            *benchfile2;
  FILE            *benchfp1;
  FILE            *benchfp2;
  
  char            *paramfile;

  int              voutput;             /* TRUE for verbose output */

  int              truncate;
  int              partial;

  char            *dbfile;	       /* name of seq db file             */
  ESL_SQFILE      *dbfp;   	       /* source database for negatives   */
  int              dbfmt;	       /* format code for dbfile          */
  int              db_nseq;  	       /* # of sequences in the db        */
  int              db_maxL;	       /* maximum seq length in db_lens   */
  double           fq[20];  	       /* background frequency distribution, if we're making iid negatives */

  int              nmsa;
  char            *msafile;
  char            *gnuplot;
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */
  char            *msaheader;          /* header for all msa-specific output files */

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
  int              do_viterbi;         // default is optimal accuracy
 
  E1_RATE        **R1;
  P7_RATE         *R7;
  int              nr;
  float            nr_scale;
  char            *subsmx;             /* BLUSUM62, ... */
  float            scaleR;             /* scale substitutions rate matrix. if 0, do not rescale */
  int              userates;           /* if TRUE construct the R from given rate values */
  
  E2_ALI             e2ali;
  EVOM               evomodel;
  struct rateparam_s rateparam;
 
  MSA_STAT          *mstat;               /* statistics of the input alignment */
  float             *msafrq;

  int                infmt;
  E2_OPT             e2optimize;          /* type of optimization (parameters and/or time) */

  float              fixtime;

  char              *method;
  float              tol;
  int                verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  /* options for msa */
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
 /* options for msa */
  { "--local",        eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE for a local alignment, default is global",                                             0 },  
  { "--truncate",     eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to truncate the b,c sequences",                                                        0 }, 
  { "--partial",      eslARG_STRING,     NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to truncate+add other random nts to the b,c sequences",                                0 }, 
   /* options for tree */
  { "--fixtime",      eslARG_REAL,       NULL,   NULL,     "x>=0",      NULL, NULL,  NULL,               "TRUE for using a fix time for the evolutionary models of a pair",                           0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
 /* Selecting the alphabet rather than autoguessing it */
  { "--amino",        eslARG_NONE,    "TRUE",    NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is protein sequence data",                                                  2 },
  { "--dna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is DNA sequence data",                                                      2 },
  { "--rna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is RNA sequence data",                                                      2 },
  /* msa format */
  { "--informat",  eslARG_STRING,       NULL,    NULL,       NULL,     NULL,  NULL,  NULL,               "specify format",                                                                            0 },
  
/* Options controlling negative segment randomization method  */
  { "--mono",        eslARG_NONE,   "default",   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "shuffle preserving monoresidue composition",                                                2 },
  { "--di",          eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "shuffle preserving mono- and di-residue composition",                                       2 },
  { "--markov0",     eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "generate with 0th order Markov properties per input",                                       2 },
  { "--markov1",     eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "generate with 1st order Markov properties per input",                                       2 },
  { "--reverse",     eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "reverse each input",                                                                        2 },
  { "--iid",         eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "generate random iid sequence for negatives",                                                2 },

  /* other options */
  { "--single",      eslARG_NONE,       FALSE,   NULL,      NULL,      NULL,  NULL, NULL,                "embed one, not two domains in each positive",                                               4 },
  { "--minDPL",       eslARG_INT,       "100",   NULL,      NULL,      NULL,  NULL, NULL,                "minimum segment length for DP shuffling",                                                   4 },
  { "--tol",         eslARG_REAL,      "1e-3",   NULL,      NULL,      NULL,  NULL, NULL,                "tolerance",                                                                                 0 },
  { "--seed",         eslARG_INT,         "0",   NULL,    "n>=0",      NULL,  NULL, NULL,                "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msa>";
static char banner[] = "progressive alignmen using e2 evolutionary model";

static int create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);

static int process_dbfile(struct cfg_s *cfg, char *dbfile, int dbfmt);
static int doctor_triomsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa, 
			  ESL_MSA **ret_msa_ab, ESL_MSA **ret_msa_ac, ESL_MSA **ret_msa_bc, ESL_MSA **ret_rmsa_bc);
static int truncate_triomsa(struct cfg_s *cfg, ESL_MSA *msa, int *ret_L1, int *ret_L2);
static int select_domain(struct cfg_s *cfg, ESL_MSA *msa, int *ret_ldom, int *ret_L1, int *ret_L2);
static int expand_partialtriomsa(ESL_GETOPTS *go, struct cfg_s *cfg,  ESL_MSA **omsa, ESL_MSA **omsa_bc, int L1, int L2);
static int domain_notallgaps(ESL_ALPHABET *abc, ESL_DSQ *dsq, int L);
static int set_gaps(ESL_ALPHABET *abc, ESL_DSQ *dsq, int L);
static int set_random_segment(ESL_GETOPTS *go, struct cfg_s *cfg, FILE *logfp, ESL_DSQ *dsq, int L);
static int merge_alignments(struct cfg_s *cfg, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **ret_msa);
static int run_e2msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *rmsa, ESL_MSA *msa);
static int run_triomsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *rmsa, ESL_MSA *msa1, ESL_MSA *msa2);
static int run_benchmark(struct cfg_s *cfg, FILE *benchfp, char *method, ESL_MSA *rmsa, MSA_STAT *mrstat, ESL_MSA *emsa, MSA_STAT *mestat, float sc);

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, struct cfg_s *ret_cfg)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  struct cfg_s cfg;
  int          r;
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

   /*  output file */
  if ( esl_opt_IsOn(go, "-o") ) {
    if ((cfg.outfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else cfg.outfp = stdout;
  
  /* other options */
  cfg.tol        = esl_opt_GetReal   (go, "--tol");
  cfg.verbose    = esl_opt_GetBoolean(go, "-v");
  cfg.voutput    = esl_opt_GetBoolean(go, "--voutput");
  cfg.truncate   = esl_opt_GetBoolean(go, "--truncate");
  cfg.dbfile     = NULL;
  cfg.dbfmt      = eslSQFILE_FASTA;

  if (esl_opt_IsOn(go, "--partial")) {
    cfg.dbfile = esl_opt_GetString(go, "--partial");
    
    if (cfg.abc->type == eslAMINO) esl_composition_SW34(cfg.fq);
    else                           esl_vec_DSet(cfg.fq, cfg.abc->K, 1.0 / (double) cfg.abc->K);
    
    /* Open and process the dbfile; make sure it's in the same alphabet */
    process_dbfile(&cfg, cfg.dbfile, cfg.dbfmt);
  }


    /* tree options */
  cfg.fixtime    = esl_opt_IsOn(go, "--fixtime")? esl_opt_GetReal(go, "--fixtime") : -1.0;
 
  cfg.msafrq     = NULL;
  cfg.subsmx     = esl_opt_GetString(go, "--mx");
  cfg.mode       = esl_opt_IsOn(go, "--local")? e2_LOCAL : e2_GLOBAL;
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
  esl_sprintf(&cfg.method, "e2msa.%s", e1_rate_EvomodelType(cfg.evomodel)); 

  /* the paramfile takes precedent over individually feed values above */
  cfg.paramfile = NULL;
  if (esl_opt_IsOn(go, "--paramfile")) {
    cfg.paramfile = esl_opt_GetString(go, "--paramfile");
    status = e1_rate_ReadParamfile(cfg.paramfile, &cfg.rateparam, &cfg.evomodel, cfg.errbuf, cfg.verbose);
    if (status != eslOK) esl_fatal("Failed to read paramfile %s\n%s", cfg.paramfile, cfg.errbuf);
  }

  cfg.e2ali = E2;
  cfg.e2optimize = OPTNONE;
 
  /* open outfile with benchmarking results */
  cfg.benchfile1 = NULL;
  cfg.benchfile2 = NULL;
  cfg.benchfp1   = NULL;
  cfg.benchfp2   = NULL;
  esl_sprintf(&cfg.benchfile1, "%s.%s.bench1", cfg.outheader, e1_rate_EvomodelType(cfg.evomodel)); 
  esl_sprintf(&cfg.benchfile2, "%s.%s.bench2", cfg.outheader, e1_rate_EvomodelType(cfg.evomodel)); 
  if ((cfg.benchfp1 = fopen(cfg.benchfile1, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.benchfile1);
  if ((cfg.benchfp2 = fopen(cfg.benchfile2, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.benchfile2);
  fprintf(cfg.benchfp1, "#MSA_NAME           TP_hom   TRUE_hom  FOUND_hom      \n");
  fprintf(cfg.benchfp2, "#MSA_NAME           TP_hom   TRUE_hom  FOUND_hom      \n");

  cfg.T = NULL;

  cfg.R1 = NULL;
  cfg.bg = e1_bg_Create(cfg.abc);

  /* Create the evolutionary rate model */
    cfg.nr    = 1;
    cfg.R1    = malloc(sizeof(E1_RATE) * cfg.nr);
    cfg.R1[0] = NULL;
    for (r = 0; r < cfg.nr; r ++) {
      cfg.R1[r] = e1_rate_CreateWithValues(cfg.abc, cfg.evomodel, cfg.rateparam, cfg.subsmx, NULL, NULL, TRUE, cfg.tol, cfg.errbuf, cfg.verbose);
      if (cfg.R1[r] == NULL) { printf("Bad rate model.\n"); esl_fatal(cfg.errbuf); }
      if (1||cfg.verbose)   printf("rI %f rM %f rD %f | ldE %f muE %f ldI %f muI %f muA %f | model %s\n",  
				    cfg.R1[r]->rI, cfg.R1[r]->rM, cfg.R1[r]->rD, 
				    cfg.R1[r]->ldE[e1R_S], cfg.R1[r]->muE[e1R_S], cfg.R1[r]->ldE[e1R_I], cfg.R1[r]->muE[e1R_I], cfg.R1[r]->muA[e1R_S], 
				    e1_rate_EvomodelType(cfg.R1[r]->evomodel));
    }
    
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
  ESL_GETOPTS    *go;
  struct cfg_s    cfg;
  ESL_MSAFILE   *afp = NULL;
  ESL_MSA        *msa = NULL;             /* the input alignment      */
  ESL_MSA        *msa_ab = NULL;          /* pairwise (a,b)           */
  ESL_MSA        *msa_ac = NULL;          /* pairwise (a,c)           */
  ESL_MSA        *msa_bc = NULL;          /* pairwise (b,c)           */
  ESL_MSA        *rmsa_bc = NULL;         /* reference pairwise (b,c) */
  int             r;
  int             status = eslOK;
  int             hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* Open the MSA file */
  status = esl_msafile_Open(&(cfg.abc), cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  /* read the MSA */
  while ((hstatus = esl_msafile_Read(afp, &msa)) != eslEOF) {
    if (hstatus != eslOK) esl_msafile_ReadFailure(afp, status);
    cfg.nmsa ++;
    if (msa->nseq != 3) esl_fatal("this works only for 3seq alignments");

   /* outheader for all msa-output files */
    msamanip_OutfileHeader((msa->acc)?msa->acc:cfg.msafile, &cfg.msaheader); 
   
    esl_msa_ConvertDegen2X(msa); 
    esl_msa_Hash(msa);

    /* Initializations */
    cfg.pli = e2_pipeline_Create(go, 1, 100, 100);
    
    status = doctor_triomsa(go, &cfg, &msa, &msa_ab, &msa_ac, &msa_bc, &rmsa_bc);
    if (1||cfg.verbose) {/* the given trio alignment */
      fprintf(stdout, "\ntrimsa\n");
      esl_msafile_Write(stdout, msa,    eslMSAFILE_STOCKHOLM);
      esl_msafile_Write(stdout, msa_ab, eslMSAFILE_STOCKHOLM);
      esl_msafile_Write(stdout, msa_ac, eslMSAFILE_STOCKHOLM);
      esl_msafile_Write(stdout, msa_bc, eslMSAFILE_STOCKHOLM);
    }

    /* reference msa aveid and avematch */
    msamanip_CStats(cfg.abc, rmsa_bc, &cfg.mstat);
    msamanip_CBaseComp(cfg.abc, rmsa_bc, cfg.bg->f, &cfg.msafrq);
    
    if (cfg.voutput) {/* the MSAProbs output alignment */
      fprintf(stdout, "\nMSAProbs alignment\n");
      esl_msafile_Write(stdout, rmsa_bc, eslMSAFILE_STOCKHOLM);
    }

    /* fixedtime  indirect (b,c) alignment */
    status = run_triomsa(go, &cfg, rmsa_bc, msa_ab, msa_ac);
    if (status != eslOK) esl_fatal("%s Failed to run triomsa", cfg.errbuf);
   
    /* optimal time direct pairwise alignment (b,c) */
    status = run_e2msa(go, &cfg, rmsa_bc, msa_bc);
    if (status != eslOK) esl_fatal("%s Failed to run e2msa", cfg.errbuf);
    
    esl_msa_Destroy(msa_ab);  msa_ab  = NULL;
    esl_msa_Destroy(msa_ac);  msa_ac  = NULL;
    esl_msa_Destroy(msa_bc);  msa_bc  = NULL;
    esl_msa_Destroy(rmsa_bc); rmsa_bc = NULL;
    esl_msa_Destroy(msa);     msa     = NULL;
    free(cfg.msafrq); cfg.msafrq = NULL;
  }
  
  /* cleanup */
  if (cfg.R1) {
    for (r = 0; r < cfg.nr; r ++) e1_rate_Destroy(cfg.R1[r]);
    free(cfg.R1);
  }   
  e2_pipeline_Destroy(cfg.pli);
  if (cfg.T) esl_tree_Destroy(cfg.T);
  esl_stopwatch_Destroy(cfg.w);
  if (cfg.bg)  e1_bg_Destroy(cfg.bg);     
  fclose(cfg.outfp);
  free(cfg.outheader);
  fclose(cfg.benchfp1);
  fclose(cfg.benchfp2);
  free(cfg.benchfile1);
  free(cfg.benchfile2);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  esl_msafile_Close(afp);
  return 0;
}

/* Open the source sequence database for negative subseqs;
 * upon return, cfg->dbfp is open (digital, SSI indexed);
 * cfg->db_maxL and cfg->db_nseq are set.
 */
static int
process_dbfile(struct cfg_s *cfg, char *dbfile, int dbfmt)
{
  ESL_SQ     *sq    = esl_sq_CreateDigital(cfg->abc);
  int         status;
  
  /* Open the sequence file in digital mode */
  status = esl_sqfile_OpenDigital(cfg->abc, dbfile, dbfmt, NULL, &(cfg->dbfp));
  if      (status == eslENOTFOUND) esl_fatal("No such file %s", dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of seqfile %s unrecognized.", dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  /* Read info on each sequence */
  cfg->db_nseq   = 0;
  cfg->db_maxL   = 0;
  while ((status = esl_sqio_ReadInfo(cfg->dbfp, sq)) == eslOK) {
    cfg->db_maxL = ESL_MAX(sq->L, cfg->db_maxL);
    cfg->db_nseq++;
    esl_sq_Reuse(sq);
  }
  if (status != eslEOF) esl_fatal("Something went wrong with reading the seq db");

  /* Open SSI index */
  if (esl_sqfile_OpenSSI(cfg->dbfp, NULL) != eslOK) esl_fatal("Failed to open SSI index file");
  if (cfg->dbfp->data.ascii.ssi->nprimary != cfg->db_nseq)     esl_fatal("oops, nprimary != nseq");

  esl_sq_Destroy(sq);
  return eslOK;
}


static int
doctor_triomsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa, ESL_MSA **ret_msa_ab, ESL_MSA **ret_msa_ac, ESL_MSA **ret_msa_bc, ESL_MSA **ret_rmsa_bc)
{
  ESL_MSA *msa_ab  = NULL;
  ESL_MSA *msa_ac  = NULL;
  ESL_MSA *msa_bc  = NULL;
  ESL_MSA *rmsa_bc = NULL;
  ESL_MSA *msa;
  int     *useme = NULL;
  int      L1, L2;
  int      status;

  msa = *omsa;

  /* first truncate and get the reference alignment and stats */
  if (cfg->truncate || cfg->partial) {
    status = truncate_triomsa(cfg, msa, &L1, &L2);
    if (status != eslOK) goto ERROR;
  }
  
  useme = malloc(msa->nseq * sizeof(int));
  useme[0] = FALSE;
  useme[1] = TRUE;
  useme[2] = TRUE;
  esl_msa_SequenceSubset(msa, useme, &msa_bc);

  /* the msaprobs alignment of bc */
  status = MSAProbs_Align(msa_bc, &rmsa_bc, cfg->errbuf, cfg->verbose);
  if (status != eslOK) esl_fatal("%s Failed to run msaprobs", cfg->errbuf);
  esl_msa_Destroy(msa_bc); msa_bc = NULL;

  /* then (if partial) add the non-homologous regions both to the tmsa and rmsa */
  if (cfg->partial) {
    expand_partialtriomsa(go, cfg, &msa, &rmsa_bc, L1, L2);
  }
  
  /* the three pairwise alignments */
  useme[0] = TRUE;
  useme[1] = TRUE;
  useme[2] = FALSE;
  esl_msa_SequenceSubset(msa, useme, &msa_ab);
  useme[1] = FALSE;
  useme[2] = TRUE;
  esl_msa_SequenceSubset(msa, useme, &msa_ac);
  useme[0] = FALSE;
  useme[1] = TRUE;
  useme[2] = TRUE;
  esl_msa_SequenceSubset(msa, useme, &msa_bc);
  
  esl_msa_MinimGaps(msa_ab, NULL, "-", FALSE);
  esl_msa_MinimGaps(msa_ac, NULL, "-", FALSE);
  esl_msa_MinimGaps(msa_bc, NULL, "-", FALSE);
  
  *omsa = msa;
  *ret_msa_ab  = msa_ab;
  *ret_msa_ac  = msa_ac;
  *ret_msa_bc  = msa_bc;
  *ret_rmsa_bc = rmsa_bc;
  free(useme);
  return eslOK;
  
 ERROR:
  if (msa_ab)  esl_msa_Destroy(msa_ab);
  if (msa_ac)  esl_msa_Destroy(msa_ac);
  if (msa_bc)  esl_msa_Destroy(msa_bc);
  if (rmsa_bc) esl_msa_Destroy(rmsa_bc);
  if (useme) free(useme);
  return status;
}

static int
select_domain(struct cfg_s *cfg, ESL_MSA *msa, int *ret_ldom, int *ret_L1, int *ret_L2)
{
  ESL_DSQ *dsq;
  int      L    = msa->alen;
  int      ldom = -1;
  int      L1   = -1;
  int      L2   = -1;
  int      i;
  int      n;
  int      okdomain = FALSE;

  /* select the length of the domains from a Gaussian (mu=L/3,sigma=L/5)
   */
  do {
    ldom  = (int)fabs(esl_rnd_Gaussian(cfg->r, L/3., L/5.));
  } while (L - ldom < (int)((float)L/5.) || ldom < (int)((float)L/5.));
  
  /* select where domain stats */
  do {
    i = esl_rnd_Roll(cfg->r, L - ldom + 1 ); /* i = 0..L' */

    /* make sure domains are not all gaps */
    okdomain = TRUE;
    for (n = 0; n < msa->nseq; n ++) {
      dsq = msa->ax[n];
      
      okdomain = domain_notallgaps(cfg->abc, dsq+1+L1, ldom);
      if (okdomain == FALSE) break;
    }

  } while (okdomain == FALSE);

  /* now 1           .. i         = gaps (if i==0, there's none); 
   *     i+1         .. i+ldom    = domain
   *     i+ldom+1     .. L        = gaps (if i+ldom==L, there's none);
   */
  L1 = i;			
  L2 = L - ldom - L1;

  *ret_ldom = ldom;
  *ret_L1   = L1;
  *ret_L2   = L2;

  return eslOK;
}

static int
domain_notallgaps(ESL_ALPHABET *abc, ESL_DSQ *dsq, int L)
{
  int l;

  for (l = 0; l < L; l ++)
    if (dsq[l] < abc->K) break;
  if (l == L) return FALSE;

  return TRUE;
}

static int
truncate_triomsa(struct cfg_s *cfg,  ESL_MSA *msa, int *ret_L1, int *ret_L2)
{
  ESL_DSQ *dsq;
  int      ldom;
  int      L1, L2;
  int      n;
  int      status;
  
 select_domain(cfg, msa, &ldom, &L1, &L2);
 printf("ldom %d L1 %d L2 %d alen %"PRId64"\n", ldom, L1, L2, msa->alen);
 if (ldom < 0 || L1 < 0 || L2 < 0) ESL_XFAIL(eslFAIL, cfg->errbuf, "select_domain() failed");

  /* for the last 2 seqs in the triomsa */
  for (n = 1; n < msa->nseq; n ++) {
    dsq = msa->ax[n];
    
    set_gaps(cfg->abc, dsq+1,         L1);
    set_gaps(cfg->abc, dsq+L1+ldom+1, L2);
  }
   
  if (ret_L1) *ret_L1 = L1;
  if (ret_L2) *ret_L2 = L2;
  return eslOK;

 ERROR:
  return status;
}

static int
expand_partialtriomsa(ESL_GETOPTS *go, struct cfg_s *cfg,  ESL_MSA **omsa, ESL_MSA **omsa_bc, int L1, int L2)
{
  ESL_MSA  *new = NULL;
  ESL_MSA  *new_bc = NULL;
  ESL_MSA  *msa;
  ESL_MSA  *msa_bc;
  ESL_SQ  **sql = NULL;
  ESL_SQ  **sqr = NULL;
  int       n;
  int       x;
  int       i;
  int       status;

  msa    = *omsa;
  msa_bc = *omsa_bc;
  if (msa->nseq    != 3) { status = eslFAIL; goto ERROR; }
  if (msa_bc->nseq != 2) { status = eslFAIL; goto ERROR; }

  if (msa->aseq    == NULL) esl_msa_Textize(msa);
  if (msa_bc->aseq == NULL) esl_msa_Textize(msa_bc);

  // create the flanking regions
  ESL_ALLOC(sql, sizeof(ESL_SQ) * 2);
  ESL_ALLOC(sqr, sizeof(ESL_SQ) * 2);
  for (n = 0; n < 2; n ++) {
    sql[n] = esl_sq_CreateDigital(cfg->abc);  
    sqr[n] = esl_sq_CreateDigital(cfg->abc);  
    
    esl_sq_GrowTo(sql[n], L1);
    esl_sq_GrowTo(sqr[n], L2);
    sql[n]->n = L1;
    sqr[n]->n = L2;
    
    set_random_segment(go, cfg, NULL, sql[n]->dsq+1, L1);
    set_random_segment(go, cfg, NULL, sqr[n]->dsq+1, L2);
    
    esl_sq_Textize(sql[n]);
    esl_sq_Textize(sqr[n]);
  }
  
  // then add the flanking regions to the triomsa
  for (x = 0; x < 2; x ++) {
    new = esl_msa_Create(msa->nseq, msa->alen+L1+L2);  
    
    for (n = 0; n < msa->nseq; n ++) {
      esl_msa_SetSeqName(new, n, msa->sqname[n], -1);
     
      for (i = 0; i < new->alen; i ++) {	  
       if (n==x+1) {
	 if      (i <  L1)           new->aseq[n][i] = tolower(sql[x]->seq[i]);
	 else if (i >= L1+msa->alen) new->aseq[n][i] = tolower(sqr[x]->seq[i-(L1+msa->alen)]);
	 else                        new->aseq[n][i] = msa->aseq[n][i-L1];
       }
       else {
	 if (i < L1 || i >= L1+msa->alen) new->aseq[n][i] = '-';
	 else                             new->aseq[n][i] = msa->aseq[n][i-L1];
       }
      }
    }
    
    esl_msa_Destroy(msa); msa = NULL;
    msa = new;
  }
  
  // then add the flanking regions to the rmsa_bc
  for (x = 0; x < 2; x ++) {
   new_bc = esl_msa_Create(2, msa_bc->alen+L1+L2);  

    for (n = 0; n < 2; n ++) {
      esl_msa_SetSeqName(new_bc, n, msa_bc->sqname[n], -1);

      for (i = 0; i < new_bc->alen; i ++) {	  
	if (n==x) {
	  if      (i <  L1)              new_bc->aseq[n][i] = tolower(sql[x]->seq[i]);
	  else if (i >= L1+msa_bc->alen) new_bc->aseq[n][i] = tolower(sqr[x]->seq[i-(L1+msa_bc->alen)]);
	  else                           new_bc->aseq[n][i] = msa_bc->aseq[n][i-L1];
	}
	else {
	  if (i < L1 || i >= L1+msa_bc->alen) new_bc->aseq[n][i] = '-';
	  else                                new_bc->aseq[n][i] = msa_bc->aseq[n][i-L1];
	}
      }
    }

    esl_msa_Destroy(msa_bc); msa_bc = NULL;
    msa_bc = new_bc;
   }
  
  if (sql) {
    for (n = 0; n < 2; n ++) 
      if (sql[n]) esl_sq_Destroy(sql[n]);
    free(sql);
  }
  if (sqr) {
    for (n = 0; n < 2; n ++) 
      if (sqr[n]) esl_sq_Destroy(sqr[n]);
    free(sqr);
  }

  *omsa    = msa;
  *omsa_bc = msa_bc;
  
  //printf("O1 %d\n", x);
  //esl_msafile_Write(stdout, *omsa1, eslMSAFILE_STOCKHOLM);
  //printf("O2 %d\n", x);
  //esl_msafile_Write(stdout, *omsa2, eslMSAFILE_STOCKHOLM);
  
  return eslOK;

 ERROR:
  if (sql) {
    for (n = 0; n < 2; n ++) 
      if (sql[n]) esl_sq_Destroy(sql[n]);
    free(sql);
  }
  if (sqr) {
    for (n = 0; n < 2; n ++) 
      if (sqr[n]) esl_sq_Destroy(sqr[n]);
    free(sqr);
  }
  if (new) esl_msa_Destroy(new);
  if (new_bc) esl_msa_Destroy(new_bc);
  return status;
}


static int
set_gaps(ESL_ALPHABET *abc, ESL_DSQ *dsq, int L)
{
  int x;

  for (x = 0; x < L; x ++) { 
    dsq[x] = abc->K;
  }

  return eslOK;
}

/* Fetch in a random sequence of length <L> from the the pre-digitized
 * concatenated sequence database, select a random subseq, shuffle it
 * by the chosen algorithm; set dsq[1..L] to the resulting randomized
 * segment.
 * 
 * If <logfp> is non-NULL, append one or more "<sqname> <from> <to>"
 * fields to current line, to record where the random segment was
 * selected from. This is useful in cases where we want to track back
 * the origin of a high-scoring segment, in case the randomization
 * wasn't good enough to obscure the identity of a segment.
 * 
 */
static int
set_random_segment(ESL_GETOPTS *go, struct cfg_s *cfg, FILE *logfp, ESL_DSQ *dsq, int L)
{
  ESL_SQ  *sq           = esl_sq_CreateDigital(cfg->abc);
  int      minDPL       = esl_opt_GetInteger(go, "--minDPL");
  int      db_dependent = (esl_opt_GetBoolean(go, "--iid") == TRUE ? FALSE : TRUE);
  char    *pkey         = NULL;
  int      start, end;
  int64_t  Lseq;
  int      status;

  if (L==0) return eslOK;
  if (L > cfg->db_maxL) esl_fatal("can't fetch a segment of length %d; database max is %d\n", L, cfg->db_maxL);

  /* fetch a random subseq from the source database */
  esl_sq_GrowTo(sq, L);
  if (db_dependent) 
    {
      do {                                                     
	if (pkey != NULL) free(pkey);
	if (esl_ssi_FindNumber(cfg->dbfp->data.ascii.ssi, esl_rnd_Roll(cfg->r, cfg->db_nseq), NULL, NULL, NULL, &Lseq, &pkey) != eslOK)
	  esl_fatal("failed to look up a random seq");
      } while (Lseq < L);

      start = 1 + esl_rnd_Roll(cfg->r, Lseq-L);              
      end   = start + L - 1;
      if (esl_sqio_FetchSubseq(cfg->dbfp, pkey, start, end, sq) != eslOK) esl_fatal("failed to fetch subseq");
      esl_sq_ConvertDegen2X(sq);
    }

  /* log sequence source info: <name> <start> <end> */
  if (logfp != NULL && db_dependent) 
    fprintf(logfp, " %-25s %5d %5d", pkey, start, end); 

  /* Now apply the appropriate randomization algorithm */
  if      (esl_opt_GetBoolean(go, "--mono"))    status = esl_rsq_XShuffle  (cfg->r, sq->dsq, L, sq->dsq);
  else if (esl_opt_GetBoolean(go, "--di")) {
    if (L < minDPL)                             status = esl_rsq_XShuffle  (cfg->r, sq->dsq, L, sq->dsq);
    else                                        status = esl_rsq_XShuffleDP(cfg->r, sq->dsq, L, cfg->abc->Kp, sq->dsq);
  } 
  else if (esl_opt_GetBoolean(go, "--markov0")) status = esl_rsq_XMarkov0  (cfg->r, sq->dsq, L, cfg->abc->Kp, sq->dsq);
  else if (esl_opt_GetBoolean(go, "--markov1")) status = esl_rsq_XMarkov1  (cfg->r, sq->dsq, L, cfg->abc->Kp, sq->dsq);
  else if (esl_opt_GetBoolean(go, "--reverse")) status = esl_rsq_XReverse  (sq->dsq, L, sq->dsq);
  else if (esl_opt_GetBoolean(go, "--iid"))     status = esl_rsq_xIID      (cfg->r, cfg->fq, cfg->abc->K, L, sq->dsq);
  if (status != eslOK) esl_fatal("esl's shuffling failed");

  memcpy(dsq, sq->dsq+1, sizeof(ESL_DSQ) * L);
  esl_sq_Destroy(sq);
  free(pkey);
  return eslOK;
}
  

static int
create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  int status;

  /* the TREE */
  status = e2_tree_UPGMA(&cfg->T, 0, NULL, msa, cfg->msafrq, cfg->r, cfg->pli, cfg->R1[0], NULL, cfg->bg, NULL, cfg->e2ali,
			 cfg->mode, cfg->do_viterbi, cfg->fixtime, -1.0, cfg->tol, cfg->errbuf, cfg->verbose);
  if (1||cfg->verbose) Tree_Dump(cfg->outfp, cfg->T, "Tree");
  
  cfg->treeavgt = esl_tree_er_AverageBL(cfg->T);

   return eslOK;
}

static int
run_triomsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *rmsa, ESL_MSA *msa1, ESL_MSA *msa2)
{
  ESL_MSA  *e2msa1 = NULL;
  ESL_MSA  *e2msa2 = NULL;
  ESL_MSA  *e2msa  = NULL;
  MSA_STAT *m1stat = NULL;
  MSA_STAT *m2stat = NULL;
  MSA_STAT *e2mstat = NULL;
  float     sc1, sc2, sc;
  int       status;

  cfg->fixtime = esl_opt_IsOn(go, "--fixtime")? esl_opt_GetReal(go, "--fixtime") : -1.0;

  esl_msa_Digitize(cfg->abc, msa1, NULL);
  esl_msa_Digitize(cfg->abc, msa2, NULL);

  /* Create the tree once is enough */
  status = create_tree(go, cfg, msa1);
  if (status != eslOK)  { esl_fatal(cfg->errbuf); }
  
  status = e2_msa(cfg->r, cfg->R1[0], NULL, 0, NULL, msa1, cfg->msafrq, cfg->T, &e2msa1, &sc1, cfg->pli, cfg->bg, NULL, cfg->e2ali, cfg->e2optimize, 
		  cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
  
  status = e2_msa(cfg->r, cfg->R1[0], NULL, 0, NULL, msa2, cfg->msafrq, cfg->T, &e2msa2, &sc2, cfg->pli, cfg->bg, NULL, cfg->e2ali, cfg->e2optimize, 
		  cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose); 
     
  msamanip_CStats(cfg->abc, e2msa1, &m1stat);
  msamanip_CStats(cfg->abc, e2msa2, &m2stat);
  
  /* total score per position */
  sc = 0.5 * (sc1/m1stat->avgsqlen + sc2/m2stat->avgsqlen);

  /* The leaves-only alignments */
  msamanip_MSALeaves(&e2msa1, FALSE);   
  msamanip_MSALeaves(&e2msa2, FALSE);   
  
  /* merge the two alignments */
  status = merge_alignments(cfg, e2msa1, e2msa2, &e2msa);
  if (status != eslOK) goto ERROR;
  msamanip_CStats(cfg->abc, e2msa, &e2mstat);
  
   if (cfg->voutput) { /* write output alignment */
    /* print output info */
    fprintf(cfg->outfp, "\ntriomsa-%s alignment\n", e1_rate_EvomodelType(cfg->evomodel));
    esl_msafile_Write(cfg->outfp, e2msa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, e2msa, e2mstat);
  }
  
  if ((status = run_benchmark(cfg, cfg->benchfp1, cfg->method, rmsa, cfg->mstat, e2msa, e2mstat, sc)) != eslOK) goto ERROR;
  
  if (e2msa)  esl_msa_Destroy(e2msa);  
  if (e2msa1) esl_msa_Destroy(e2msa1);  
  if (e2msa2) esl_msa_Destroy(e2msa2);  
  return eslOK;

 ERROR:
  if (e2msa)  esl_msa_Destroy(e2msa);  
  if (e2msa1) esl_msa_Destroy(e2msa1);  
  if (e2msa2) esl_msa_Destroy(e2msa2);  
  return status;
}

static int
merge_alignments(struct cfg_s *cfg, ESL_MSA *msa1, ESL_MSA *msa2, ESL_MSA **ret_msa)
{
  ESL_SQ   **sq         = NULL;
  P7_TRACE **tr         = NULL;
  P7_TRACE **tr1        = NULL;
  P7_TRACE **tr2        = NULL;
  ESL_MSA   *msa        = NULL;
  int       *matassign1 = NULL;
  int       *matassign2 = NULL;
  int        optflags   = p7_DEFAULT;
  int        apos;
  int        M;
  int        nseq;
  int        status;

  if (1||cfg->verbose) {
    fprintf(stdout, "\n A-B alignment\n");
    esl_msafile_Write(stdout, msa1, eslMSAFILE_STOCKHOLM);
    fprintf(stdout, "\n A-C alignment\n");
    esl_msafile_Write(stdout, msa2, eslMSAFILE_STOCKHOLM);
  }
 
  ESL_ALLOC(matassign1, sizeof(int) * (msa1->alen+1));
  ESL_ALLOC(matassign2, sizeof(int) * (msa2->alen+1));
  esl_vec_ISet(matassign1, msa1->alen+1, FALSE);
  esl_vec_ISet(matassign2, msa2->alen+1, FALSE);

  /* assign M to positions with residues in the top sequence */
  for (apos = 0; msa1->aseq[0][apos] != '\0'; apos++) {
    if (esl_abc_CIsResidue(cfg->abc, msa1->aseq[0][apos])) matassign1[apos+1] = TRUE;
  }
  for (apos = 0; msa2->aseq[0][apos] != '\0'; apos++) {
    if (esl_abc_CIsResidue(cfg->abc, msa2->aseq[0][apos])) matassign2[apos+1] = TRUE;
  }
  
  if ((tr1 = malloc(sizeof(P7_TRACE *) * msa1->nseq))  == NULL)  esl_fatal("malloc");
  if ((tr2 = malloc(sizeof(P7_TRACE *) * msa2->nseq))  == NULL)  esl_fatal("malloc");
  
  esl_msa_Digitize(cfg->abc, msa1, NULL);
  esl_msa_Digitize(cfg->abc, msa2, NULL);
  p7_trace_FauxFromMSA(msa1, matassign1, optflags, tr1);
  p7_trace_FauxFromMSA(msa2, matassign2, optflags, tr2);
  
  if (msa1->nseq != msa2->nseq) { printf("bad alignments \n"); goto ERROR; }
  if (tr1[0]->M  != tr2[0]->M)  { printf("bad traces \n");     goto ERROR; }
  M    = tr1[0]->M; 
  nseq = msa1->nseq;
  if ((tr = malloc(sizeof(P7_TRACE *) * nseq))  == NULL)  esl_fatal("malloc");
  if ((sq = malloc(sizeof(ESL_SQ   *) * nseq))  == NULL)  esl_fatal("malloc");
  tr[0] = tr1[1];
  tr[1] = tr2[1];
  if (0) {
    p7_trace_Dump(stdout, tr[0]);
    p7_trace_Dump(stdout, tr[1]);
  }
 
  esl_sq_FetchFromMSA(msa1, 1, &sq[0]);
  esl_sq_FetchFromMSA(msa2, 1, &sq[1]);

  p7_tracealign_Seqs(sq, tr, nseq, M, p7_ALL_CONSENSUS_COLS, NULL, &msa);
  esl_msa_MinimGaps(msa, NULL, "-", FALSE);

  if (1||cfg->verbose) {
    fprintf(stdout, "\n B-C alignment\n");
    esl_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM);
  }
  
  *ret_msa = msa;

  free(matassign1);
  free(matassign2);
  free(tr);
  p7_trace_DestroyArray(tr1, msa1->nseq);
  p7_trace_DestroyArray(tr2, msa2->nseq);
  if (sq) {
    esl_sq_Destroy(sq[0]);
    esl_sq_Destroy(sq[1]);
    free(sq);
  }
  return eslOK;

 ERROR:
  if (matassign1) free(matassign1);
  if (matassign1) free(matassign2);
  if (tr ) free(tr);
  if (tr1) p7_trace_DestroyArray(tr1, msa1->nseq);
  if (tr2) p7_trace_DestroyArray(tr2, msa2->nseq);
  if (sq) {
    if (sq[0]) esl_sq_Destroy(sq[0]);
    if (sq[1]) esl_sq_Destroy(sq[1]);
    free(sq);
  }
  return status;
}

static int
run_e2msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *rmsa, ESL_MSA *msa)
{
  ESL_MSA  *e2msa = NULL;
  MSA_STAT *alle2mstat = NULL;
  MSA_STAT *e2mstat = NULL;
  float     sc;
  int       status;

  /* optimal time direct pairwise alignment */
  cfg->fixtime = -1.0;

  esl_stopwatch_Start(cfg->w); 

  /* create the tree */
  esl_msa_Digitize(cfg->abc, msa, NULL);
  status = create_tree(go, cfg, msa);
  if (status != eslOK)  { esl_fatal(cfg->errbuf); }

 /* main function */
  status = e2_msa(cfg->r, cfg->R1[0], NULL, 0, NULL, msa, cfg->msafrq, cfg->T, &e2msa, &sc, cfg->pli, cfg->bg, NULL, cfg->e2ali, cfg->e2optimize, 
		  cfg->mode, cfg->do_viterbi, cfg->tol, cfg->errbuf, cfg->verbose);
  if (status != eslOK)  { esl_fatal(cfg->errbuf); }

  esl_stopwatch_Stop(cfg->w);
  msamanip_CStats(cfg->abc, e2msa, &alle2mstat);
  
  /* The leaves-only alignment */
  msamanip_MSALeaves(&e2msa, FALSE);   
  msamanip_CStats(cfg->abc, e2msa, &e2mstat);
  e2mstat->anclen = alle2mstat->anclen; // transfer the len of the ancestral sequence (cannot be calculated from the leaves-only msa)

  /* score per position */
  sc /= e2mstat->avgsqlen;

  if (cfg->voutput) { /* write output alignment */
    /* print output info */
    fprintf(cfg->outfp, "\ne2msa-%s alignment\n", e1_rate_EvomodelType(cfg->evomodel));
    esl_msafile_Write(cfg->outfp, e2msa, eslMSAFILE_STOCKHOLM);
    msamanip_DumpStats(cfg->outfp, e2msa, e2mstat);
  }

  if ((status = run_benchmark(cfg, cfg->benchfp2, cfg->method, rmsa, cfg->mstat, e2msa, e2mstat, sc)) != eslOK) goto ERROR;
  
  if (e2msa) esl_msa_Destroy(e2msa); 
  if (e2mstat) free(e2mstat);
  if (alle2mstat) free(alle2mstat);
  return eslOK;

 ERROR:
  if (e2msa) esl_msa_Destroy(e2msa);
  return status;
}



static int 
run_benchmark(struct cfg_s *cfg, FILE *benchfp, char *method, ESL_MSA *rmsa, MSA_STAT *mrstat, ESL_MSA *emsa, MSA_STAT *mestat, float sc)
{
  char   *msaname = NULL;
  double  time;
  int     lcu_r;  // if TRUE, consider lower case in reference alignment as unaligned
  int     lcu_e;  // if TRUE, consider lower case in inferred alignment as unaligned
  int     status; 
  
  esl_sprintf(&msaname, "%s.%d", cfg->msaheader, cfg->nmsa); 
 
  time = cfg->w->user + cfg->w->sys;
  
  lcu_r = FALSE; /* it is all in upper case and aligned */
  lcu_e = TRUE; /* the estimated aligment is created with 
		 * p7_trace_FauxFromMSA() such that the residues in the "a" seq are matches
		 * p7_tracealign_Seqs() and includes lower case
		 * for residues that did not align to the more divergent "a" sequence
		 */
  status = FastSP_Benchmark(benchfp, msaname, method, cfg->abc, rmsa, mrstat, emsa, mestat, sc, cfg->treeavgt, time, lcu_r, lcu_e, cfg->errbuf, cfg->verbose);
  if (status != eslOK) goto ERROR;

  free(msaname);
  return eslOK;

 ERROR:
  if (msaname) free(msaname);
  return status;
}

