/* pevomark -- 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esl_composition.h"
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
#include "fastsp.h"
#include "msatree.h"
#include "msamanip.h"
#include "minimize.h"
#include "miscellaneous.h"
#include "phmmer.h"

#define ALPHOPTS  "--amino,--dna,--rna"                          /* Exclusive options for alphabet choice */
#define EFFOPTS "--eent,--eclust,--eset,--enone"                 /* Exclusive options for effective sequence number calculation */
#define METHOPTS  "--ehmmsearch,--hmmsearch,--hmmsearch3,--ephmmer,--phmmer,--phmmer3,--ncbiblast"     /* Exclusive options for search method */

/* Exclusive search method */
enum method_e { ehmmsearch,      //ehmmsearch
		hmmsearch,       //hmmsearch trunck
		hmmsearch3,      //hmmsearch 3.1b1
		ephmmer,         //ephmmer
		phmmer,          //phmmer
		phmmer3,         //phmmer 3.1b1
		ncbiblast,       //ncbiblast+
 		undetermined
              };  

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

  char            *msafileom_phmm;
  char            *seqfile;
  char            *msafile;

  char            *outfile;
  FILE            *outfp;
  
  int              voutput;             /* TRUE for verbose output */

  enum method_e    method;

  float            fixpid;

  float            fixtime;
  float            mintime;
  float            maxtime;
  float            inctime;

  float            Eval;
  float            incdomE;
  float            incE;

  int              EmL;
  int              EmN;
  int              EvL;
  int              EvN;
  int              EfL;
  int              EfN;
  float            Eft;

  float            popen;
  float            pextend;

  char            *mx;

  int              infmt;
  int              nmsa;

  char            *gnuplot;
  char            *outheader;          /* header for all output files */
  char            *msaheader;          /* header for all msa-specific output files */

  P7_PIPELINE     *pli;
  P7_BG           *bg;	   	       /* null model (copies made of this into threads) */

  char            *hmmfile;            /* for e2hmmer case: open input HMM file */
  
  double           treeavgt;
  ESL_TREE        *T;

  int              max;                /* for running --max with phmmer or hmmsearch */

  P7_HMMFILE      *hfp;                /* input HMM file */
  P7_HMM          *hmm;                /* HMM            */

  int              enone;
 
  float            tol;
  int              verbose;
};


struct hmmpid_data {
  char    *hmmfile;
  int        N;      // number of samples in hmmemit
  double   time;
  double   pid_target;
  double   pid_obs;
  double   firststep;
  double   tol;
  char    *errbuf;
  int      verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--dout",         eslARG_STRING,    FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "output directory",                                                                          0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
  /* Selecting the alphabet rather than autoguessing it */
  { "--amino",        eslARG_NONE,    "TRUE",    NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is protein sequence data",                                                  2 },
  { "--dna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is DNA sequence data",                                                      2 },
  { "--rna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is RNA sequence data",                                                      2 },
  /* Selecting method */
  { "--ehmmsearch",   eslARG_NONE,    "TRUE",    NULL,       NULL, METHOPTS,  NULL,  NULL,               "search method",                                                                             2 },
  { "--hmmsearch",    eslARG_NONE,      FALSE,   NULL,       NULL, METHOPTS,  NULL,  NULL,               "search method",                                                                             2 },
  { "--hmmsearch3",   eslARG_NONE,      FALSE,   NULL,       NULL, METHOPTS,  NULL,  NULL,               "search method",                                                                             2 },
  { "--ephmmer",      eslARG_NONE,      FALSE,   NULL,       NULL, METHOPTS,  NULL,  NULL,               "search method",                                                                             2 },
  { "--phmmer",       eslARG_NONE,      FALSE,   NULL,       NULL, METHOPTS,  NULL,  NULL,               "search method",                                                                             2 },
  { "--phmmer3",      eslARG_NONE,      FALSE,   NULL,       NULL, METHOPTS,  NULL,  NULL,               "search method",                                                                             2 },
  { "--ncbiblast",    eslARG_NONE,      FALSE,   NULL,       NULL, METHOPTS,  NULL,  NULL,               "search method",                                                                             2 },
  /* Control phmmer/hmmsearch */
  { "--max",          eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "use --max option for phmmer or hmmsearch",                                                  0 },
  /* msa format */
  { "--informat",   eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            0 },
  /* Alternative effective sequence weighting strategies */
  { "--eent",    eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to achieve relative entropy target",  5 },
  { "--eclust",  eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "eff seq # is # of single linkage clusters",            5 },
  { "--enone",   eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "no effective seq # weighting: just use nseq",          5 },
  { "--eset",    eslARG_REAL,   NULL,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "set eff seq # for all models to <x>",                  5 },
  { "--ere",     eslARG_REAL,   NULL,  NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set minimum rel entropy/position to <x>",  5 },
  { "--esigma",  eslARG_REAL, "45.0",  NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set sigma param to <x>",                   5 },
  { "--eid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--eclust",    NULL, "for --eclust: set fractional identity cutoff to <x>",  5 },
  /* other options */
  { "--fixpid",       eslARG_REAL,       NULL,   NULL,"0<=x<=100",   NULL,    NULL,  NULL,               "TRUE for using a fix pid for the ehmm model",           0 },
  { "--fixtime",      eslARG_REAL,       NULL,   NULL,     "x>=0",   NULL,    NULL,"--mintime,--maxtime","TRUE for using a fix time for the ehmm model",          0 },
  { "--mintime",      eslARG_REAL,       NULL,   NULL,     "x>=0",   NULL,    NULL, "--fixtime",         "TRUE for using a range of ehmm models",                 0 },
  { "--maxtime",      eslARG_REAL,     "10.0",   NULL,     "x>=0",   NULL,    NULL, "--fixtime",         "maxtime for using a range of ehmms",                    0 },
  { "--inctime",      eslARG_REAL,     "0.01",   NULL,     "x>=0",   NULL,    NULL, "--fixtime",         "time increment for the ehmm",                           0 },
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",         eslARG_INT,        "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
 };
static char usage[]  = "[-options] <trainmsa> <testseqs>";
static char banner[] = "pevomark";

static int    run_hmmbuild          (struct cfg_s *cfg, ESL_MSA *msa, P7_HMM **ret_hmm);
static int    run_search            (struct cfg_s *cfg, FILE *outfile, ESL_MSA *msa, P7_HMM *hmm);
static int    run_hmmsearch         (struct cfg_s *cfg, FILE *outfp, P7_HMM *hmm);
static int    run_hmmsearch_trunk   (struct cfg_s *cfg, FILE *outfp, char *hmmfile);
static int    run_hmmsearch_3p1b1   (struct cfg_s *cfg, FILE *outfp, char *hmmfile);
static int    run_phmmer            (struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa, int trunk);
static int    run_ehmmsearch        (struct cfg_s *cfg, FILE *outfile, P7_HMM *hmm);
static int    run_ehmmsearch_onetime(struct cfg_s *cfg, FILE *outfp, char *hmmfile, float time);
static int    hmm_from_phmmer(struct cfg_s *cfg, FILE *fp, ESL_SQ *sq, P7_HMM **ret_hmm);
static int    run_ephmmer           (struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa);
static int    run_ephmmer_onetime   (struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa, float time, char **ret_hmmfile);
static int    run_ncbiblast         (struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa);

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
  cfg.seqfile = NULL;
  if (esl_opt_ArgNumber(go) != 2) { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.msafile  = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.seqfile  = esl_opt_GetArg(go, 2)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));

  cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));

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
  cfg.bg = p7_bg_Create(cfg.abc);
  
  cfg.w = esl_stopwatch_Create(); 
  
  /* method */
  if      (esl_opt_GetBoolean(go, "--ehmmsearch"))   cfg.method = ehmmsearch;
  else if (esl_opt_GetBoolean(go, "--hmmsearch"))    cfg.method = hmmsearch;
  else if (esl_opt_GetBoolean(go, "--hmmsearch3"))   cfg.method = hmmsearch3;
  else if (esl_opt_GetBoolean(go, "--ephmmer"))      cfg.method = ephmmer;
  else if (esl_opt_GetBoolean(go, "--phmmer"))       cfg.method = phmmer;
  else if (esl_opt_GetBoolean(go, "--phmmer3"))      cfg.method = phmmer3;
  else if (esl_opt_GetBoolean(go, "--ncbiblast"))    cfg.method = ncbiblast;
  else                                               cfg.method = undetermined;

   /*  output file */
  if ( esl_opt_IsOn(go, "-o") ) {
    if ((cfg.outfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else cfg.outfp = stdout;
  
  /* other options */
  cfg.tol        = esl_opt_GetReal   (go, "--tol");
  cfg.verbose    = esl_opt_GetBoolean(go, "-v");
  cfg.voutput    = esl_opt_GetBoolean(go, "--voutput");
  cfg.fixpid     = esl_opt_IsOn(go, "--fixpid")?  esl_opt_GetReal(go, "--fixpid") : -1.0; 
  cfg.fixtime    = esl_opt_IsOn(go, "--fixtime")? esl_opt_GetReal(go, "--fixtime") : -1.0; 
  cfg.mintime    = esl_opt_IsOn(go, "--mintime")? esl_opt_GetReal(go, "--mintime") : -1.0; 
  cfg.maxtime    = esl_opt_GetReal(go, "--maxtime"); 
  cfg.inctime    = esl_opt_GetReal(go, "--inctime"); 

  if (esl_opt_IsOn(go, "--mintime")) {
    if (cfg.maxtime <= cfg.mintime)              esl_fatal("bad range of times");
    if (cfg.mintime + cfg.inctime > cfg.maxtime) esl_fatal("bad range of times");
  }

  cfg.Eval    = 200.;
  cfg.incdomE = 200.;
  cfg.incE    = 200.;

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

  cfg.enone = FALSE;
  if (esl_opt_IsOn(go, "--enone")) cfg.enone = TRUE;

  /* hmmer/phmmer */
  if (esl_opt_IsOn(go, "--max")) cfg.max = TRUE;
  else                           cfg.max = FALSE;

  /* open outfile with benchmarking results */
  cfg.outfile = NULL;
  cfg.outfp   = NULL;
  if (esl_opt_IsOn(go, "--dout")) esl_sprintf(&cfg.outfile, "%s/%s.out", esl_opt_GetString(go, "--dout"), cfg.outheader); 
  else                            esl_sprintf(&cfg.outfile, "%s.out", cfg.outheader); 
  if ((cfg.outfp = fopen(cfg.outfile, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.outfile);
  fprintf(cfg.outfp, "#ehmm <hmm_name> <hmm_time> <hmm_alen> <hmm_avgsqlen> <hmm_pid> <hmm_me> <hmm_mre>\n");
  fprintf(cfg.outfp, "#<iE> <isc> <E> <sc> <sqfrom> <sqto> <targetname> <queryname>\n");

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
  ESL_MSAFILE   *afp    = NULL;
  ESL_MSA        *msa    = NULL;       /* the query alignment             */
  P7_HMM         *hmm    = NULL;
  int             status = eslOK;
  int             hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* Open the query MSA file */
  status = esl_msafile_Open(&(cfg.abc), cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  /* read the query MSA */
  while ((hstatus = esl_msafile_Read(afp, &msa)) != eslEOF) {
    if (hstatus != eslOK) esl_msafile_ReadFailure(afp, status);
    cfg.nmsa ++;
    
    /* outheader for all msa-output files */
    msamanip_OutfileHeader((msa->acc)?msa->name:cfg.msafile, &cfg.msaheader); 
    
    esl_msa_ConvertDegen2X(msa); 
    esl_msa_Hash(msa);
    
    /* create the query hmm with hmmbuild */
    status = run_hmmbuild(&cfg, msa, &hmm);
    if (status != eslOK) esl_fatal("Failed to create hmm form msa %s", msa->name);
    
    /* run search */
    status = run_search(&cfg, cfg.outfp, msa, hmm);
    if (status != eslOK) esl_fatal("%s\nFailed to run search", cfg.errbuf);
    
    esl_msa_Destroy(msa);         msa = NULL;
    if (hmm) p7_hmm_Destroy(hmm); hmm  = NULL;
  }
  
  esl_stopwatch_Destroy(cfg.w); 
  p7_bg_Destroy(cfg.bg);
   
  fclose(cfg.outfp);
  free(cfg.outheader);
  free(cfg.msaheader);
  free(cfg.outfile);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  esl_msafile_Close(afp);
  if (hmm) p7_hmm_Destroy(hmm);
  if (msa) esl_msa_Destroy(msa); 
 return 0;
}

static int
run_hmmbuild(struct cfg_s *cfg, ESL_MSA *msa, P7_HMM **ret_hmm)
{
  char          tmphmmfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpmsafile[16] = "esltmpXXXXXX"; /* tmpfile template */
  FILE         *fp  = NULL;
  P7_HMMFILE   *hfp = NULL;                      /* open HMM file    */
  P7_HMM       *hmm = NULL;
  char         *args = NULL;
  char         *opts = NULL;
  char         *s = NULL;
  int           status;
  
  if (cfg->method == ephmmer || cfg->method == phmmer || cfg->method == phmmer3 || cfg->method == ncbiblast) 
    {
      *ret_hmm = hmm;
      return eslOK;
    }

  /* MSA input in stockholm format */
  if ((status = esl_tmpfile_named(tmpmsafile,  &fp))                          != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create ms file");
  if ((status = esl_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_STOCKHOLM)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "Failed to write STOCKHOLM file\n");
  fclose(fp);
  
  /* run hmmbuild */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))                     != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hmm file");
  
  esl_sprintf(&opts, " ");
  if (cfg->enone) esl_sprintf(&opts, "--enone "); 

  if      (cfg->method == ehmmsearch)
    esl_sprintf(&args, "%s/lib/hmmer4/src/programs/hmmbuild --amino %s %s %s>/dev/null", s, opts, tmphmmfile, tmpmsafile);
  else if (cfg->method == hmmsearch)
    esl_sprintf(&args, "%s/lib/hmmer4/src/programs/hmmbuild --amino %s %s %s>/dev/null", s, opts, tmphmmfile, tmpmsafile);
  else if (cfg->method == hmmsearch3)
    esl_sprintf(&args, "%s/lib/hmmer-3.1b1/src/hmmbuild --amino %s %s %s>/dev/null",     s, opts, tmphmmfile, tmpmsafile);
  else ESL_XFAIL(status, cfg->errbuf, "error in run_hmmbuild()");
  system(args);
  //printf("%s\n", args);
  fclose(fp);
  
  /* read the hmm */
  status = p7_hmmfile_OpenE(tmphmmfile, NULL, &hfp, cfg->errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", tmphmmfile, cfg->errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                tmphmmfile, cfg->errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, tmphmmfile, cfg->errbuf);  

  status = p7_hmmfile_Read(hfp, &cfg->abc, &hmm);
  if      (status != eslOK)        p7_Fail("Unexpected error %d reading HMM from file %s.\n%s\n",             status, tmphmmfile, cfg->errbuf);  
  
  *ret_hmm = hmm;

  remove(tmphmmfile);
  remove(tmpmsafile);
  
  p7_hmmfile_Close(hfp);
  if (opts != NULL) free(opts);
  if (args != NULL) free(args);
  return eslOK;
  
 ERROR:
  remove(tmphmmfile);
  remove(tmpmsafile);
  if (hfp) p7_hmmfile_Close(hfp);
  if (opts) free(opts);
  if (args) free(args);
  if (hmm) p7_hmm_Destroy(hmm);
  return status;
}

static int
run_search(struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa, P7_HMM *hmm)
{
  int status;

  switch(cfg->method) {
  case ehmmsearch:
    status = run_ehmmsearch(cfg, outfp, hmm);
    break;
  case hmmsearch:
    status = run_hmmsearch(cfg, outfp, hmm);
    break;
  case hmmsearch3:
    status = run_hmmsearch(cfg, outfp, hmm);
    break;
  case ephmmer:
    status = run_ephmmer(cfg, outfp, msa);
    break;
  case phmmer:
    status = run_phmmer(cfg, outfp, msa, TRUE);
    break;
  case phmmer3:
    status = run_phmmer(cfg, outfp, msa, FALSE);
    break;
  case ncbiblast:
    status = run_ncbiblast(cfg, outfp, msa);
    break;
  default:
    ESL_XFAIL(eslFAIL, cfg->errbuf, "no such methods");
  }
  
  return eslOK;

 ERROR:
  return status;
}

static int
run_hmmsearch(struct cfg_s *cfg, FILE *outfp, P7_HMM *hmm)
{
  char             tmphmmfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  FILE            *fp  = NULL;
  int              status;

  /* the hmm tmpfile */
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))  != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tmphmm file");
  if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm))  != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to write hmm file");
  fclose(fp);
  
  if      (cfg->method == hmmsearch)  status = run_hmmsearch_trunk(cfg, outfp, tmphmmfile);
  else if (cfg->method == hmmsearch3) status = run_hmmsearch_3p1b1(cfg, outfp, tmphmmfile);
  else ESL_XFAIL(status, cfg->errbuf, "error in run_hmmserach()");
   
  remove(tmphmmfile);
  return eslOK;

 ERROR:
  remove(tmphmmfile);
  return status;
}

static int
run_hmmsearch_trunk(struct cfg_s *cfg, FILE *outfp, char *hmmfile)
{
  char             tmptblfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmphitfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  FILE            *fp  = NULL;
  struct domain_s *dom = NULL;
  int              ndom;
  char            *args = NULL;
  char            *s = NULL;
  char            *hmmname = NULL;
  int              hmmALEN   = -1;
  float            hmmSQLEN  = -1.0;
  float            hmmPID    = -1.0;
  float            hmmPMATCH = -1.0;
  float            hmmME     = -1.0;
  float            hmmMRE    = -1.0;
  int              d;
  int              status;

  /* run ehmmsearch */
  if ((status = esl_tmpfile_named(tmptblfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tbl file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmphitfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hit file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
  fclose(fp);

  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if (cfg->max) 
    esl_sprintf(&args, "%s/lib/hmmer4/src/programs/hmmsearch --cpu 1 --max -E %f --incdomE %f --incE %f --domtbl %s -A %s %s %s > %s", 
		s, cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile, hmmfile, cfg->seqfile, tmpoutfile);  
  else
    esl_sprintf(&args, "%s/lib/hmmer4/src/programs/hmmsearch --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s %s %s > %s", 
		s, cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile, hmmfile, cfg->seqfile, tmpoutfile);  
  system(args);
  printf("%s\n", args);

  status = misc_ParseHMMERdomtbl(tmptblfile, &dom, &ndom, cfg->errbuf); if (status != eslOK) goto ERROR;
  if (cfg->voutput) misc_DomainDump(stdout, dom, ndom);

  for (d = 0; d < ndom; d++) {
    /* ehmm stats using ehmmemeit */
    status = HMMER_EmitStats(hmmfile, cfg->abc, 1.0, &hmmname, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, &hmmME, &hmmMRE, cfg->errbuf);
    if (status != eslOK) goto ERROR;
 
    misc_DomainWrite(outfp, &dom[d], hmmALEN, hmmSQLEN, hmmPID, hmmPMATCH, hmmME, hmmMRE);
  }

  if (dom) free(dom);

  remove(tmptblfile);
  remove(tmphitfile);
  remove(tmpoutfile);
  free(args);
  return eslOK;

 ERROR:
  remove(tmptblfile);
  remove(tmphitfile);
  remove(tmpoutfile);

  if (args) free(args);
  if (dom) free(dom);

  return status;
}

static int
run_hmmsearch_3p1b1(struct cfg_s *cfg, FILE *outfp, char *hmmfile)
{
  char             tmptblfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmphitfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  FILE            *fp  = NULL;
  struct domain_s *dom = NULL;
  int              ndom;
  char            *args = NULL;
  char            *s = NULL;
  char            *hmmname = NULL;
  int              hmmALEN   = -1;
  float            hmmSQLEN  = -1.0;
  float            hmmPID    = -1.0;
  float            hmmPMATCH = -1.0;
  float            hmmME     = -1.0;
  float            hmmMRE    = -1.0;
  int              d;
  int              status;

  /* run ehmmsearch */
  if ((status = esl_tmpfile_named(tmptblfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tbl file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmphitfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hit file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
  fclose(fp);

  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if (cfg->max)
    esl_sprintf(&args, "%s/lib/hmmer-3.1b1/src/hmmsearch --cpu 1 --max -E %f --incdomE %f --incE %f --domtbl %s -A %s %s %s > %s", 
		s, cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile, hmmfile, cfg->seqfile, tmpoutfile);  
  else
    esl_sprintf(&args, "%s/lib/hmmer-3.1b1/src/hmmsearch --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s %s %s > %s", 
		s, cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile, hmmfile, cfg->seqfile, tmpoutfile);  
  system(args);
  
  status = misc_ParseHMMERdomtbl(tmptblfile, &dom, &ndom, cfg->errbuf); if (status != eslOK) goto ERROR;
  if (cfg->voutput) misc_DomainDump(stdout, dom, ndom);

  for (d = 0; d < ndom; d++) {
    /* ehmm stats using ehmmemeit */
    status = HMMER_EmitStats(hmmfile, cfg->abc, 1.0, &hmmname, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, &hmmME, &hmmMRE, cfg->errbuf);
    if (status != eslOK) goto ERROR;
    
    misc_DomainWrite(outfp, &dom[d], hmmALEN, hmmSQLEN, hmmPID, hmmPMATCH, hmmME, hmmMRE);
   }

  if (dom) free(dom);

  remove(tmptblfile);
  remove(tmphitfile);
  remove(tmpoutfile);
  free(args);
  return eslOK;

 ERROR:
  remove(tmptblfile);
  remove(tmphitfile);
  remove(tmpoutfile);

  if (args) free(args);
  if (dom) free(dom);

  return status;
}


static int
run_ehmmsearch(struct cfg_s *cfg, FILE *outfp, P7_HMM *hmm)
{
  char             tmphmmfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  FILE            *fp  = NULL;
  float            time;
  float            inctime;
  int              status;

  /* the hmm tmpfile */
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))  != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tmphmm file");
  if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm))  != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to write hmm file");
  fclose(fp);
  
  if (cfg->mintime < 0.) {
    if (cfg->fixpid > 0.) EPHMMER_pid2time(hmm, cfg->fixpid, &time, cfg->tol, cfg->errbuf, cfg->verbose); 
    else                  time = cfg->fixtime;
    status = run_ehmmsearch_onetime(cfg, outfp, tmphmmfile, time);
    if (status != eslOK) goto ERROR;
  }
  else {
    inctime = cfg->inctime;
    for (time = cfg->mintime; time <= cfg->maxtime; time += inctime) {
      inctime = (time < 1.0)? cfg->inctime : (time < 10.0)? 10.0*cfg->inctime : 100.0*cfg->inctime;
      status = run_ehmmsearch_onetime(cfg, outfp, tmphmmfile, time);
      if (status != eslOK) esl_fatal("%s\nFailed to run ehmmsearch at time %f", cfg->errbuf, time);
    } 
  }
  
  remove(tmphmmfile);
  return eslOK;

 ERROR:
  remove(tmphmmfile);
  return status;
}


static int
run_ehmmsearch_onetime(struct cfg_s *cfg, FILE *outfp, char *hmmfile, float time)
{
  char             tmptblfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmphitfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  FILE            *fp  = NULL;
  struct domain_s *dom = NULL;
  int              ndom;
  char            *opts = NULL;
  char            *args = NULL;
  char            *s = NULL;
  char            *hmmname = NULL;
  int              hmmALEN   = -1;
  float            hmmSQLEN  = -1.0;
  float            hmmPID    = -1.0;
  float            hmmPMATCH = -1.0;
  float            hmmME     = -1.0;
  float            hmmMRE    = -1.0;
  int              d;
  int              status;

  /* run ehmmsearch */
  if ((status = esl_tmpfile_named(tmptblfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tbl file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmphitfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hit file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
  fclose(fp);

  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if (cfg->max) esl_sprintf(&opts, "--max --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s ", cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile);
  else          esl_sprintf(&opts, "      --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s ", cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile);
  if (time < 0.) esl_sprintf(&args, "%s/src/programs/ehmmsearch %s %s %s > %s", 
			     s,       opts, hmmfile, cfg->seqfile, tmpoutfile);  
  else           esl_sprintf(&args, "%s/src/programs/ehmmsearch %s --fixtime %f %s %s > %s", 
			     s,       opts, time, hmmfile, cfg->seqfile, tmpoutfile);
  system(args);
  //printf("%s\n", args);

  status = misc_ParseHMMERedomtbl(tmptblfile, &dom, &ndom, cfg->errbuf); if (status != eslOK) goto ERROR;
  if (cfg->voutput) misc_DomainDump(stdout, dom, ndom);

  for (d = 0; d < ndom; d++) {
    /* ehmm stats using ehmmemeit */
    status = HMMER_EmitStats(hmmfile, cfg->abc, dom[d].time, &hmmname, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, &hmmME, &hmmMRE, cfg->errbuf);
    if (status != eslOK) goto ERROR;

    misc_DomainWrite(outfp, &dom[d], hmmALEN, hmmSQLEN, hmmPID, hmmPMATCH, hmmME, hmmMRE);
  }

  if (dom)  free(dom);
  if (opts) free(opts);
  remove(tmptblfile);
  remove(tmphitfile);
  remove(tmpoutfile);
  free(args);
  return eslOK;

 ERROR:
  remove(tmptblfile);
  remove(tmphitfile);
  remove(tmpoutfile);

  if (args) free(args);
  if (dom ) free(dom);
  if (opts) free(opts);
  return status;
}

static int
run_phmmer(struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa, int trunk)
{
  char            *besthmmfile = NULL;
  ESL_SQ          *sq = NULL;
  FILE            *fp   = NULL;
  char            *opts = NULL;
  char            *args = NULL;
  char            *s    = NULL;
  struct domain_s *dom  = NULL;
  struct domain_s *dombest = NULL;
  int              ndom;
  int              ndombest = 0;
  float            Ebest = eslINFINITY;
  int              n;
  char            *hmmname = NULL;
  int              hmmALEN   = -1;
  float            hmmSQLEN  = -1.0;
  float            hmmPID    = -1.0;
  float            hmmPMATCH = -1.0;
  float            hmmME     = -1.0;
  float            hmmMRE    = -1.0;
  int              d;
  int              status;


  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;

  // one sequence at the time
  for (n = 0; n < msa->nseq; n ++) {

    esl_sq_FetchFromMSA(msa, n, &sq);
    
    char             tmptargetsqfile[16] = "esltmpXXXXXX"; /* tmpfile template */
    char             tmphmmfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    char             tmptblfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    char             tmphitfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    char             tmpoutfile[16]      = "esltmpXXXXXX"; /* tmpfile template */

    if ((status = esl_tmpfile_named(tmptargetsqfile, &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create targetsq file");
    esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE);
    fclose(fp);

    /* run phmmer */
    if ((status = esl_tmpfile_named(tmptblfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tbl file");
    fclose(fp);
    if ((status = esl_tmpfile_named(tmphitfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hit file");
    fclose(fp);
    if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
    fclose(fp);
    
    if ((status = esl_tmpfile_named(tmphmmfile, &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hmm file");
    if ((status = hmm_from_phmmer(cfg, fp, sq, NULL)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create targetsq file");
    fclose(fp);

    if (cfg->max) esl_sprintf(&opts, "--max --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s ", cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile);
    else          esl_sprintf(&opts, "      --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s ", cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile);

    if (trunk)
      esl_sprintf(&args, "%s/lib/hmmer4/src/programs/phmmer    %s %s %s > %s", 
		  s,       opts, tmptargetsqfile, cfg->seqfile, tmpoutfile);  
    else 
     esl_sprintf(&args, "%s/lib/hmmer-3.1b1/src/phmmer         %s %s %s > %s", 
		  s,       opts, tmptargetsqfile, cfg->seqfile, tmpoutfile);  
    system(args);
    printf("%s\n", args);


    status = misc_ParseHMMERdomtbl(tmptblfile, &dom, &ndom, cfg->errbuf); if (status != eslOK) goto ERROR;
    
    if (ndom > 0 && dom[0].E < Ebest) { 
      Ebest = dom[0].E; 
      misc_DomainClone(ndom, dom, &ndombest, &dombest); 
      
      remove(besthmmfile);
      esl_strdup(tmphmmfile, -1, &besthmmfile);
      if (cfg->voutput) misc_DomainDump(stdout, dom, ndom);
    }
    else remove(tmphmmfile);

    remove(tmptargetsqfile);
    remove(tmptblfile);
    remove(tmphitfile);
    remove(tmpoutfile);

    free(dom);  dom  = NULL;
    free(opts); opts = NULL;
    esl_sq_Destroy(sq); sq = NULL;
  }
  
  for (d = 0; d < ndombest; d++) {
    esl_strdup(msa->name, -1, &(dombest[d].queryname)); // the "query" is the msa
    
    /* ehmm stats using ehmmemeit */
    status = HMMER_EmitStats(besthmmfile, cfg->abc, dombest[d].time, &hmmname, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, &hmmME, &hmmMRE, cfg->errbuf);
    if (status != eslOK) goto ERROR;
    
    misc_DomainWrite(outfp, &dombest[d], hmmALEN, hmmSQLEN, hmmPID, hmmPMATCH, hmmME, hmmMRE);
  }
  
  remove(besthmmfile);

  if (dom)  free(dom);
  if (opts) free(opts);
  if (sq) esl_sq_Destroy(sq);
  return eslOK;
  
 ERROR:
  if (dom)  free(dom);
  if (opts) free(opts);
  if (sq) esl_sq_Destroy(sq);
  return status;
}

static int
hmm_from_phmmer(struct cfg_s *cfg, FILE *fp, ESL_SQ *sq, P7_HMM **ret_hmm)
{
  P7_HMM  *hmm  = NULL;    
  int      status;

  status = PHMMER_gethmm(sq, &hmm, cfg->bg, cfg->EmL, cfg->EmN, cfg->EvL, cfg->EvN, cfg->EfL, cfg->EfN, cfg->Eft, cfg->popen, cfg->pextend, cfg->mx, cfg->errbuf);
  if (status != eslOK) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run PHMMER_gethmm\n%s\n", cfg->errbuf);

  if (fp) {
    if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "HMM save failed");
  }

  if (ret_hmm) *ret_hmm = hmm;
  else         p7_hmm_Destroy(hmm);
  
  return eslOK;
  
 ERROR:
  if (hmm) p7_hmm_Destroy(hmm);
  return status;
}

static int
run_ephmmer(struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa)
{
  float   time;
  float   inctime;
  int     status;
 
  if (cfg->mintime < 0.) {
    status = run_ephmmer_onetime(cfg, outfp, msa, -1.0, NULL);
    if (status != eslOK) goto ERROR;
  }
  else {
    inctime = cfg->inctime;
    for (time = cfg->mintime; time <= cfg->maxtime; time += inctime) {
      inctime = (time < 1.0)? cfg->inctime : (time < 10.0)? 10.0*cfg->inctime : 100.0*cfg->inctime;
      status = run_ephmmer_onetime(cfg, outfp, msa, time, NULL);
      if (status != eslOK) goto ERROR;
    } 
  }
  
   return eslOK;

 ERROR:
   return status;
}

static int
run_ephmmer_onetime(struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa, float time, char **ret_hmmfile)
{
  P7_HMM          *hmm = NULL;
  char            *besthmmfile = NULL;
  ESL_SQ          *sq = NULL;
  FILE            *fp   = NULL;
  char            *opts = NULL;
  char            *args = NULL;
  char            *s    = NULL;
  struct domain_s *dom  = NULL;
  struct domain_s *dombest = NULL;
  float            usetime;
  int              ndom;
  int              ndombest = 0;
  float            Ebest = eslINFINITY;
  int              n;
  char            *hmmname = NULL;
  int              hmmALEN  = -1;
  float            hmmSQLEN   = -1.0;
  float            hmmPID    = -1.0;
  float            hmmPMATCH = -1.0;
  float            hmmME     = -1.0;
  float            hmmMRE    = -1.0;
  int              d;
  int              status;

   if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;

  // one sequence at the time
  for (n = 0; n < msa->nseq; n ++) {
    
    char             tmptargetsqfile[16] = "esltmpXXXXXX"; /* tmpfile template */
    char             tmphmmfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    char             tmptblfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    char             tmphitfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    char             tmpoutfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    
    esl_sq_FetchFromMSA(msa, n, &sq);
 
    if (cfg->fixpid > 0.) {
      hmm_from_phmmer(cfg, NULL, sq, &hmm);
      EPHMMER_pid2time(hmm, cfg->fixpid, &usetime, cfg->tol, cfg->errbuf, cfg->verbose);
    }
    else {
      if (time < 0) usetime = cfg->fixtime;
      else          usetime = time;
    }

   /* run ephmmer */
    if ((status = esl_tmpfile_named(tmptblfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tbl file");
    fclose(fp);
    if ((status = esl_tmpfile_named(tmphitfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hit file");
    fclose(fp);
    if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
    fclose(fp);
    if ((status = esl_tmpfile_named(tmphmmfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hmm file");
    if ((status = hmm_from_phmmer(cfg, fp, sq, NULL)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create targetsq file");
    fclose(fp);
    if ((status = esl_tmpfile_named(tmptargetsqfile, &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create targetsq file");
    esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE);
    fclose(fp);
    
    if (cfg->max) esl_sprintf(&opts, "--max --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s ", cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile);
    else          esl_sprintf(&opts, "      --cpu 1 -E %f --incdomE %f --incE %f --domtbl %s -A %s ", cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmphitfile);


    if (usetime < 0.) esl_sprintf(&args, "%s/src/programs/ephmmer             %s %s %s > %s", 
				  s,          opts, tmptargetsqfile, cfg->seqfile, tmpoutfile);  
    else             esl_sprintf(&args, "%s/src/programs/ephmmer --fixtime %f %s %s %s > %s", 
				 s, usetime, opts, tmptargetsqfile, cfg->seqfile, tmpoutfile);
    system(args);
    printf("%s\n", args);

    status = misc_ParseHMMERedomtbl(tmptblfile, &dom, &ndom, cfg->errbuf); if (status != eslOK) goto ERROR;
    
    if (ndom > 0 && dom[0].E < Ebest) { 
      Ebest = dom[0].E; 
      misc_DomainClone(ndom, dom, &ndombest, &dombest); 
      
      remove(besthmmfile);
      esl_strdup(tmphmmfile, -1, &besthmmfile);
      if (cfg->voutput) misc_DomainDump(stdout, dom, ndom);
    }
    else remove(tmphmmfile);
    
    remove(tmptargetsqfile);
    remove(tmptblfile);
    remove(tmphitfile);
    remove(tmpoutfile);
    
    free(dom);  dom  = NULL;
    free(opts); opts = NULL;
    esl_sq_Destroy(sq); sq = NULL;
  }

  for (d = 0; d < ndombest; d++) {
    esl_strdup(msa->name, -1, &(dombest[d].queryname)); // the "query" is the msa

    /* ehmm stats using ehmmemeit */
    status = HMMER_EmitStats(besthmmfile, cfg->abc, dombest[d].time, &hmmname, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmPMATCH, &hmmME, &hmmMRE, cfg->errbuf);
    if (status != eslOK) goto ERROR;

    if (outfp) misc_DomainWrite(outfp, &dombest[d], hmmALEN, hmmSQLEN, hmmPID, hmmPMATCH, hmmME, hmmMRE);
  }
  
  if (ret_hmmfile) {
    *ret_hmmfile = besthmmfile;
  }
  else  {
    remove(besthmmfile);    
    if (besthmmfile) free(besthmmfile);
  }

  if (dom)  free(dom);  
  if (opts) free(opts); 
  if (hmm) p7_hmm_Destroy(hmm);
  if (sq) esl_sq_Destroy(sq); 
  return eslOK;
  
 ERROR:  
  if (dom)  free(dom);  
  if (opts) free(opts); 
  if (hmm) p7_hmm_Destroy(hmm);
  if (sq) esl_sq_Destroy(sq);
  return status;
}

static int
run_ncbiblast(struct cfg_s *cfg, FILE *outfp, ESL_MSA *msa)
{
  ESL_SQ          *sq = NULL;
  FILE            *fp   = NULL;
  char            *args = NULL;
  char            *s    = NULL;
  struct domain_s *dom  = NULL;
  struct domain_s *dombest = NULL;
  int              ndom;
  int              ndombest = 0;
  float            Ebest = eslINFINITY;
  int              n;
  int              d;
  int              status;


  if ("BLASTDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("BLASTDIR")) == NULL) return eslENOTFOUND;

  // one sequence at the time
  for (n = 0; n < msa->nseq; n ++) {

    esl_sq_FetchFromMSA(msa, n, &sq);
    
    char             tmptargetsqfile[16] = "esltmpXXXXXX"; /* tmpfile template */
    char             tmpoutfile[16]      = "esltmpXXXXXX"; /* tmpfile template */
    
    if ((status = esl_tmpfile_named(tmptargetsqfile, &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create targetsq file");
    esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE);
    fclose(fp);

    /* run ncbiblast+ */
    if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
    fclose(fp);
    
    esl_sprintf(&args, "%s/bin/blastp  --db %s --query %s -num_threads 1 -num_descriptions 9999 > %s", 
		s,     cfg->seqfile, tmptargetsqfile, tmpoutfile);  
    system(args);
     
    status = misc_ParseBLASTout(tmpoutfile, &dom, &ndom, cfg->errbuf); if (status != eslOK) goto ERROR;
    if (cfg->voutput) misc_DomainDump(stdout, dom, ndom);
    
    if (ndom > 0 && dom[0].E < Ebest) { 
      Ebest = dom[0].E; 
      misc_DomainClone(ndom, dom, &ndombest, &dombest); 
    }

    remove(tmptargetsqfile);
    remove(tmpoutfile);

    free(dom);
    esl_sq_Destroy(sq); sq = NULL;
  }
  
  for (d = 0; d < ndombest; d++) {
    esl_strdup(msa->name, -1, &(dombest[d].queryname)); // the "query" is the msa
    misc_DomainWrite(outfp, &dombest[d], -1, -1.0, -1.0, -1.0, -1.0, -1.0);
  }
    
  if (dom) free(dom);
  if (sq) esl_sq_Destroy(sq);
  return eslOK;
  
 ERROR:
  if (dom) free(dom);
  if (sq) esl_sq_Destroy(sq);
  return status;
}





