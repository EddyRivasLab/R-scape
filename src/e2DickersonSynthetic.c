/* e2DickersonSyntheetic -- 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_distance.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2Dickerson.h"
#include "e2.h"
#include "e1_emit.h"
#include "e1_rate.h"
#include "e2_profilesq.h"
#include "branchinference.h"
#include "Dickersonrates.h"
#include "msatree.h"
#include "plot.h"

#define ALIOPTS  "--e2,--e2f,--e2hmmer,--e2fhmmer"   /* Exclusive options for alignment choice */
#define EVOMOPTS "--TKF91,--TKF92,--FID,--GTKF92,--GRTKF92,--ITKF92,--LI,--LR,--AFG,--AFGR,--AFR,--AIF,--GG,--G2"              /* Exclusive options for evolutionary model choice */

static ESL_OPTIONS options[] = {
  /* name           type                                       default   env  range    toggles   reqs   incomp               help                                                                 docgroup*/
  { "-h",           eslARG_NONE,                                FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "show brief help on version and usage",                                     0 },
  { "-v",           eslARG_NONE,                                FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "be verbose?",                                                              0 },
  { "-L",           eslARG_INT,                                 "200",  NULL, "n>0",    NULL,    NULL,  NULL,               "length of ancestral sequence",                                             0 },
  { "-N",           eslARG_INT,                                  "20",  NULL, "n>0",    NULL,    NULL,  NULL,               "number of taxa per tree",                                                  0 },
  { "--ntree",      eslARG_INT,                                  "10",  NULL, "n>0",    NULL,    NULL,  NULL,               "number of trees created",                                                  0 },
  { "--tblmin",     eslARG_REAL,                              "0.001",  NULL, "x>0",    NULL,    NULL,  NULL,               "tree's min mean total branch length",                                      0 },
  { "--tblmax",     eslARG_REAL,                               "10.0",  NULL, "x>0",    NULL,    NULL,  NULL,               "tree's max mean total branch length",                                      0 },
  { "--tblinc",     eslARG_REAL,                               "0.01",  NULL, "x>0",    NULL,    NULL,  NULL,               "tree's average branch length increments",                                  0 },
  { "--nr",         eslARG_INT,                                   "1",  NULL, "n>0",    NULL,    NULL,  NULL,               "number of different rate matrices",                                        0 },
   /* options for msa */
  { "--e2",         eslARG_NONE,                                 TRUE,   NULL, NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, realign pairs of sequences",                                          0 }, 
  { "--e2f",        eslARG_NONE,                                FALSE,   NULL, NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, use alignment as is, do not realign",                                 0 }, 
  { "--e2hmmer",    eslARG_STRING,                              FALSE,   NULL, NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, realign pairs of sequences, use position-specific evomodel",          0 }, 
  { "--e2fhmmer",   eslARG_STRING,                              FALSE,   NULL, NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, use alignment as is, do not realign, use position-specific evomodel", 0 }, 
  /* Control of output */
  { "-o",            eslARG_OUTFILE,                            FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                      0 },
  { "-w",            eslARG_NONE,                               FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "TRUE to display plots",                                                    0 },
  { "--bintime",     eslARG_REAL,                              "0.05",  NULL,"x>0.0",   NULL,    NULL,  NULL,               "time interval for binning data",                                           0 },
  { "--tlinear",     eslARG_REAL,                               "2.0",  NULL, "x>=0",   NULL,    NULL,  NULL,               "max time for a linear regression fit, If x=0 do all points",               0 },
  { "--xlabel",      eslARG_STRING,                  "time parameter",  NULL,  NULL,    NULL,    NULL,  NULL,               "xlabel",                                                                   0 },
  { "--ylabel",      eslARG_STRING,      "rate of change per residue",  NULL,  NULL,    NULL,    NULL,  NULL,               "ylabel",                                                                   0 },
   /* Control of scoring system - substitutions */ 
  { "--mx",         eslARG_STRING,                         "BLOSUM62",  NULL,  NULL,    NULL,    NULL,  "--mxfile",         "substitution rate matrix choice (of some built-in matrices)",              0 },
  { "--mxfile",     eslARG_INFILE,                               NULL,  NULL,  NULL,    NULL,    NULL,  "--mx",             "read substitution rate matrix from file <f>",                              0 },
  /* Control of scoring system - indels */ 
  { "--evomodel",   eslARG_STRING,    "AIF",   NULL,       NULL,   NULL,    NULL,  EVOMOPTS,                                "evolutionary model used",                                                  0 },
  { "--muA",        eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                                    "rate of deletion of ancestral sequences",                                  0 },
  { "--muEM",       eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                                    "rate of deletion of inserted residues",                                    0 },
  { "--ldEM",       eslARG_REAL,    "0.0010",  NULL,     "x>=0",   NULL,    NULL,  NULL,                                    "rate of adding a res to an insert",                                        0 },
  { "--muED",       eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                                    "rate of deletion of inserted residues",                                    0 },
  { "--ldED",       eslARG_REAL,    "0.0010",  NULL,     "x>=0",   NULL,    NULL,  NULL,                                    "rate of adding a res to an insert",                                        0 },
  { "--muI",        eslARG_REAL,    "0.000",   NULL,     "x>=0",   NULL,    NULL,  NULL,                                    "rate of deletion of inserted sequences",                                   0 },
  { "--ldI",        eslARG_REAL,    "0.000",   NULL,     "x>=0",   NULL,    NULL,  NULL,                                    "rate of adding a whole insert",                                            0 },
  { "--rI",         eslARG_REAL,    "0.350",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                                    "fragment parameter for I",                                                 0 },
  { "--rM",         eslARG_REAL,    "0.350",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                                    "fragment parameter for M",                                                 0 },
  { "--rD",         eslARG_REAL,    "0.350",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                                    "fragment parameter for D",                                                 0 },
  /* other options */
  { "--tol",        eslARG_REAL,                             "0.0001",   NULL,  NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                0 },
  { "--seed",        eslARG_INT,                                  "0",   NULL, "n>=0",  NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                      5 },
{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
};
static char usage[]  = "[-options]  ";
static char banner[] = "calculate pairwise rates of evolution using E2";

static int msa_generate(ESL_RANDOMNESS  *r, ESL_TREE *T, int nr, E1_RATE **R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7, E2_ALI e2ali, int L, ESL_MSA **ret_msa, int incnode, double tol, char *errbuf, int verbose);
static int hmm_GenerateAlignment(ESL_RANDOMNESS *r, ESL_TREE *T, P7_RATE *R7, P7_BG *bg7, int L, ESL_MSA **ret_msa, double tol, char *errbuf, int verbose);
static int outfile_header(int nR, int L, int N, char *subsmx, float muA, float muE, float ldE, float muI, float ldI, float rI, char **ret_outheader);

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
  cfg.pfamname = NULL;
  if (esl_opt_ArgNumber(go) != 0)    { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));

  /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  cfg.abc = esl_alphabet_Create(eslAMINO);

  if ( esl_opt_IsOn(go, "-o") ) {
    if ((cfg.outfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else cfg.outfp = stdout;

  /* other options */
  cfg.tol        = esl_opt_GetReal   (go, "--tol");
  cfg.verbose    = esl_opt_GetBoolean(go, "-v");
  cfg.view       = esl_opt_GetBoolean(go, "-w");

  cfg.subsmx     = esl_opt_GetString (go, "--mx");
  cfg.nr         = esl_opt_GetInteger(go, "--nr");
  cfg.evomodel   = e1_rate_Evomodel(esl_opt_GetString(go, "--evomodel"));
  cfg.rateparam.muAM  = esl_opt_GetReal  (go, "--muA");
  cfg.rateparam.muI  = esl_opt_GetReal  (go, "--muI");
  cfg.rateparam.ldI  = esl_opt_GetReal  (go, "--ldI");
  cfg.rateparam.muEM = esl_opt_GetReal  (go, "--muEM");
  cfg.rateparam.ldEM = esl_opt_GetReal  (go, "--ldEM");
  cfg.rateparam.muED = esl_opt_GetReal  (go, "--muED");
  cfg.rateparam.ldED = esl_opt_GetReal  (go, "--ldED");
  cfg.rateparam.rI   = esl_opt_GetReal  (go, "--rI");
  cfg.rateparam.rM   = esl_opt_GetReal  (go, "--rM");
  cfg.rateparam.rD   = esl_opt_GetReal  (go, "--rD");

  cfg.bintime    = esl_opt_GetReal   (go, "--bintime");
  cfg.tlinear    = esl_opt_GetReal   (go, "--tlinear");

  cfg.xlabel     = esl_opt_GetString (go, "--xlabel");
  cfg.ylabel     = esl_opt_GetString( go, "--ylabel");

  if      (esl_opt_GetBoolean(go, "--e2"))      cfg.e2ali = E2;
  else if (esl_opt_GetBoolean(go, "--e2f"))     cfg.e2ali = E2F;
  else if (esl_opt_GetString(go, "--e2hmmer"))  cfg.e2ali = E2HMMER;
  else if (esl_opt_GetString(go, "--e2fhmmer")) cfg.e2ali = E2FHMMER;
  else                                          cfg.e2ali = E2NONE;
 
  /* if e2hmmer or e2fhmmer, read the input HMM */
  if (cfg.e2ali == E2FHMMER || cfg.e2ali == E2HMMER)
  if ((cfg.hmmfile = esl_opt_GetString(go, "--e2hmmer")) == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  if (cfg.e2ali == E2 || cfg.e2ali == E2F) cfg.bg  = e1_bg_Create(cfg.abc);
  else                                     cfg.bg7 = p7_bg_Create(cfg.abc);

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
  ESL_GETOPTS     *go = NULL;	          /* command line processing  */
  struct cfg_s     cfg;
  char            *msg = "e2_DickersonSynthetic() failed";
  char            *histofile = NULL;
  char            *ratefile = NULL;
  char            *binratefile = NULL;
  char            *histops = NULL;
  char            *rateps = NULL;
  char            *binrateps = NULL;
  char            *hmmfile;            /* for e2hmmer case: open input HMM file */
  ESL_DSQ         *dsq = NULL;
  ESL_MSA         *msa = NULL;
  float           *msafrq = NULL;
  ESL_TREE        *T = NULL;
  SPECIES         *SPE = NULL;

  P7_HMMFILE      *hfp = NULL;            /* input HMM file */
  P7_HMM          *hmm = NULL;            /* HMM            */
  P7_RATE         *R7 = NULL;             /* HMM rate       */
  EMRATE          *emR = NULL;

  double           tblmin, tblmax, tblinc;
  double           tbl;
  int              L;                /* length of ancestral sequence */
  int              N;                /* number of taxa */
  int              ntree;            /* number of trees */
  int              t;
  int              r;
  int              status = eslOK;   
 
  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    
  L      = esl_opt_GetInteger(go, "-L");
  N      = esl_opt_GetInteger(go, "-N");
  ntree  = esl_opt_GetInteger(go, "--ntree");
  tblmin = esl_opt_GetReal   (go, "--tblmin");
  tblmax = esl_opt_GetReal   (go, "--tblmax");
  tblinc = esl_opt_GetReal   (go, "--tblinc");
  
  /* Create the evolutionary rate model */
  cfg.nr       = 1;

  if (cfg.e2ali == E2 || cfg.e2ali == E2F) {
    cfg.R1 = malloc(sizeof(E1_RATE) * cfg.nr);
    for (r = 0; r < cfg.nr; r ++) {
      cfg.R1[r] = e1_rate_CreateWithValues(cfg.abc, cfg.evomodel, cfg.rateparam, cfg.subsmx, NULL, NULL, TRUE, cfg.tol, cfg.errbuf, cfg.verbose);
      cfg.R1[r]->evomodel = cfg.evomodel; 
      
      if (cfg.R1[r] == NULL) { printf("%s. bad rate model\n", cfg.errbuf); esl_fatal(msg); }
      e1_rate_SaturationTime(cfg.R1[r], 0.01, cfg.errbuf, cfg.verbose);
    }
  }
  else {  
    /* Open the query profile HMM file */
    status = p7_hmmfile_OpenE(cfg.hmmfile, NULL, &hfp, cfg.errbuf);
    if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg.hmmfile, cfg.errbuf);
    else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg.hmmfile, cfg.errbuf);
    else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg.hmmfile, cfg.errbuf);  
    
    /* <abc> is not known 'til  HMM is read. */
    status = p7_hmmfile_Read(hfp, &(cfg.abc), &hmm);
    if (status != eslOK) p7_Fail("error reading HMM file %s.\n%s\n",  cfg.hmmfile, cfg.errbuf);  
    
    /* Calculate the hmm rate */
    emR = ratematrix_emrate_Create(cfg.abc, 1);
    ratematrix_emrate_Set(cfg.subsmx, NULL, cfg.bg->f, emR, TRUE, cfg.tol, cfg.errbuf, FALSE);

    int betainf;
    int etainf;
    betainf = (cfg.rateparam.muEM > 0.)? ((cfg.rateparam.ldEM < cfg.rateparam.muEM)? cfg.rateparam.ldEM/cfg.rateparam.muEM : 1.0 ) : 0.0;
    etainf  = (cfg.rateparam.muI  > 0.)? ((cfg.rateparam.ldI  < cfg.rateparam.muI)?  cfg.rateparam.ldI/cfg.rateparam.muI : 1.0 ) : 0.0;
    if (p7_RateCalculate(stdout, hmm, cfg.bg7, emR, NULL, &R7, cfg.evomodel, betainf, -1.0, 0.001, cfg.errbuf, FALSE) != eslOK)  
      esl_fatal("%s", cfg.errbuf);
  }

  /* Initializations */
  cfg.pli = e2_pipeline_Create(go, (hmm)? hmm->M:1, 100, 100);
  if (cfg.e2ali == E2 || cfg.e2ali == E2F) cfg.bg  = e1_bg_Create(cfg.abc);
  else                                     cfg.bg7 = p7_bg_Create(cfg.abc);
  
  /* files for writing the plots, joinly for all clusters of orthologs */
  outfile_header(cfg.nr, L, N, cfg.subsmx, cfg.rateparam.muAM, cfg.rateparam.muEM, cfg.rateparam.ldEM, cfg.rateparam.muI, cfg.rateparam.ldI, cfg.rateparam.rI, &cfg.outheader);
  esl_sprintf(&histofile,   "%s.mya.histo",    cfg.outheader);
  esl_sprintf(&ratefile,    "%s.rate.plot",    cfg.outheader);
  esl_sprintf(&binratefile, "%s.binrate.plot", cfg.outheader);
  remove(histofile);
  remove(ratefile);
  remove(binratefile);
  
  /* these simulations don't use species info */
  SPE = species_Create(0); 

  for (t = 0; t < ntree; t ++) {
    /* generate random ultrametric tree T */
    if (esl_tree_Simulate(cfg.r, N, &T) != eslOK) esl_fatal(msg);
    if (esl_tree_SetTaxaParents(T)      != eslOK) esl_fatal(msg);
    if (esl_tree_Validate(T, NULL)      != eslOK) esl_fatal(msg);
    
    tbl = tblmin;
    while (tbl < tblmax) {
      /* scale the tree T */
      if (esl_tree_er_RescaleAverageTotalBL(tbl, &T, cfg.tol, cfg.errbuf, cfg.verbose) != eslOK) esl_fatal(msg);;
      if (1||cfg.verbose) {
	printf("T[%d] tbl=%f\n", t, tbl);
	Tree_Dump(stdout, T, cfg.outheader);
      }
      
      /* emit sequences along T according to evolutionary model */
      if (msa_generate(cfg.r, T, cfg.nr, cfg.R1, cfg.R7, cfg.bg, cfg.bg7, cfg.e2ali, L, &msa, FALSE, cfg.tol, cfg.errbuf, cfg.verbose)!= eslOK) esl_fatal(msg);    
      
      /* Dickerson calculation with R[0] */
      if (e2BranchInference(cfg.outfp, cfg.gnuplot, NULL, histofile, ratefile, binratefile, cfg.r, msa, msafrq, T, SPE, cfg.pli, cfg.bg, cfg.bg7, cfg.e2ali,
			    cfg.evalcutoff, cfg.R1[0], cfg.R7, cfg.noopt, cfg.do_viterbi, cfg.bintime, cfg.tlinear, FALSE, cfg.tol, cfg.view, cfg.errbuf, cfg.verbose) != eslOK) 
	{ printf("%s\n", cfg.errbuf); esl_fatal(msg); }
           
     tbl += tblinc;
     esl_msa_Destroy(msa); msa = NULL;
     free(msafrq); msafrq = NULL;
    }

    esl_tree_Destroy(T);  T = NULL;
  }

  /* plot results */ 
  esl_sprintf(&histops, "%s.ps", histofile);  
  if (plot_gplot_Histogram(cfg.gnuplot, histops, 1, &histofile, "tree distance", cfg.errbuf) == eslOK) { 
    if (cfg.view) plot_file_Display(histops);
  } else esl_fatal("%s\nbad histogram", cfg.errbuf);
  
  esl_sprintf(&rateps,    "%s.ps", ratefile);
  esl_sprintf(&binrateps, "%s.ps", binratefile);
  if (plot_gplot_JointRatesWithLinearRegression(cfg.gnuplot, rateps, 1, &ratefile, cfg.evalcutoff, cfg.tlinear, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
    if (cfg.view) plot_file_Display(rateps);
  } else esl_fatal("%s\nbad rates", cfg.errbuf);
  if (plot_gplot_JointBinnedRatesWithLinearRegression(cfg.gnuplot, binrateps, 1, &ratefile, cfg.bintime, cfg.tlinear, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
    if (cfg.view) plot_file_Display(binrateps);
  } else esl_fatal("%s\nbad binned rates", cfg.errbuf);
      
  /* cleanup */
  if (cfg.R1) {
    for (r = 0; r < cfg.nr; r ++) e1_rate_Destroy(cfg.R1[r]);
    free(cfg.R1);
  }
  if (cfg.R7) p7_RateDestroy(cfg.R7);
  free(dsq);
  fclose(cfg.outfp);
  free(cfg.gnuplot);
  esl_randomness_Destroy(cfg.r);
  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(cfg.abc);
  e2_pipeline_Destroy(cfg.pli);
  if (cfg.bg)  e1_bg_Destroy(cfg.bg);     
  if (cfg.bg7) p7_bg_Destroy(cfg.bg7);     
  if (hmm) p7_hmm_Destroy(hmm); 
  if (hfp) p7_hmmfile_Close(hfp);
  species_Destroy(SPE);
  free(histofile); 
  free(ratefile); 
  free(binratefile); 
  free(histops); 
  if (rateps)free(rateps); 
  if (binrateps) free(binrateps); 
  
  return status;
}


static int
msa_generate(ESL_RANDOMNESS *r, ESL_TREE *T, int nr, E1_RATE **R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7, E2_ALI e2ali, int L, ESL_MSA **ret_msa, int incnode, double tol, char *errbuf, int verbose)
{
  ESL_MSA *msa = NULL;     /* alignment of leaves sequences */
  ESL_MSA *msafull = NULL; /* alignment of leaves and internal node sequences */
  int     *useme = NULL;
  int      i;

  /* Generate the alignment */
  if (e2ali == E2 || e2ali == E2F) {
    if (e1_GenerateAlignment(r, T, nr, R, bg, L, &msafull, tol, errbuf, verbose) != eslOK)
      esl_fatal("%s\nfailed to generate the alignment", errbuf);
  }
  else {
    if (hmm_GenerateAlignment(r, T, R7, bg7, L, &msafull, tol, errbuf, verbose) != eslOK)
      esl_fatal("%s\nfailed to generate the alignment", errbuf);
  }

  if (verbose) 
    eslx_msafile_Write(stdout, msafull, eslMSAFILE_STOCKHOLM);
  
  /* The leaves-only alignment */
  useme = malloc(msafull->nseq * sizeof(int));
  for (i = 0; i < msafull->nseq; i++) {
    if (!incnode && strncmp(msafull->sqname[i], "v", 1) == 0) useme[i] = FALSE; 
    else                                                      useme[i] = TRUE;
  }
  if (esl_msa_SequenceSubset(msafull, useme, &msa) != eslOK)
    esl_fatal("failed to generate leaf alignment");
  if (esl_msa_MinimGaps(msa, NULL, "-", FALSE) != eslOK) 
    esl_fatal("failed to generate leaf alignment");
  
  if (esl_msa_Digitize(bg->abc, msa, NULL)!= eslOK) 
    esl_fatal("failed to digitize alignment");

  *ret_msa = msa;

  free(useme);
  esl_msa_Destroy(msafull);

  return eslOK;
}

static int
hmm_GenerateAlignment(ESL_RANDOMNESS *r, ESL_TREE *T, P7_RATE *R7, P7_BG *bg7, int L, ESL_MSA **ret_msa, double tol, char *errbuf, int verbose)
{
  ESL_MSA *msa = NULL;     /* alignment */
  
  *ret_msa = msa;
  
  return eslOK;
}


static int
outfile_header(int nR, int L, int N, char *subsmx, float muA, float muE, float ldE, float muI, float ldI, float rI, char **ret_outheader)
{      
  char *outheader = NULL;

  esl_sprintf(&outheader, "synthetic_nR%d_L%d_N%d_mu%.4f-%.4f-%.4f_ld%.4f-%.4f_%s", nR, L, N, muA, muE, muI, ldE, ldI, subsmx);
 
  *ret_outheader = outheader;

  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/
