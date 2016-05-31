/* ptransmark -- align two seqs indirectly through the hmm
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

#define ALPHOPTS "--amino,--dna,--rna"                          /* Exclusive options for alphabet choice        */
#define SHUF_OPTS "--mono,--di,--markov0,--markov1,--reverse"   /* toggle group, seq shuffling options          */

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

  char            *dbfile;	       /* name of seq db file             */
  ESL_SQFILE      *dbfp;   	       /* source database for negatives   */
  int              dbfmt;	       /* format code for dbfile          */
  int              db_nseq;  	       /* # of sequences in the db        */
  int              db_maxL;	       /* maximum seq length in db_lens   */
  double           fq[20];  	       /* background frequency distribution, if we're making iid negatives */

  char            *benchfile1;
  char            *benchfile2;
  FILE            *benchfp1;
  FILE            *benchfp2;
  
  char            *paramfile;

  int              voutput;             /* TRUE for verbose output */

  float            fixtime;
  float            mintime;
  float            maxtime;
  float            inctime;

  float            Eval;
  float            incdomE;
  float            incE;

  int              infmt;
  int              nqmsa;
  int              ntmsa;
  char            *tblfile;
  char            *qmsafile;
  char            *tmsafile;
  char            *gnuplot;
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */
  char            *msaheader;          /* header for all msa-specific output files */

  int              truncate;
  int              partial;

  P7_PIPELINE     *pli;
  P7_BG           *bg;	   	       /* null model (copies made of this into threads) */

  char            *hmmfile;            /* for e2hmmer case: open input HMM file */
  
  double           treeavgt;
  ESL_TREE        *T;

  P7_HMMFILE      *hfp;                /* input HMM file */
  P7_HMM          *hmm;                /* HMM            */
 
  MSA_STAT        *trmstat;            /* statistics of the input alignment */
  float           *trmsafrq;

  int              domax;

  float            tol;
  int              verbose;
};

struct domain_s {
  char   *name;
  int     len;
  ESL_SQ *sql;
  ESL_SQ *sqr;

  float   E;
  float   cE;
  float   sc;
  float   csc;

  char   *sqname;
  int     sqlen;
  int     sqidx;
  int     sqfrom; // [0..sqlen-1]
  int     sqto;   // [0..sqlen-1]

  float   time;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  /* options for msa */
  { "--domax",        eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to use the max domain",                                                        0 }, 
  { "--truncate",     eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to truncate the b,c sequences",                                                        0 }, 
  { "--partial",      eslARG_STRING,     NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "TRUE to truncate+add other random nts to the b,c sequences",                                0 }, 
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
 /* Selecting the alphabet rather than autoguessing it */
  { "--amino",        eslARG_NONE,      TRUE,    NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is protein sequence data",                                                  2 },
  { "--dna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is DNA sequence data",                                                      2 },
  { "--rna",          eslARG_NONE,      FALSE,   NULL,       NULL, ALPHOPTS,  NULL,  NULL,               "input alignment is RNA sequence data",                                                      2 },
  /* msa format */
  { "--informat",  eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                             0 },
/* Options controlling negative segment randomization method  */
  { "--mono",        eslARG_NONE,   "default",   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "shuffle preserving monoresidue composition",                                                2 },
  { "--di",          eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "shuffle preserving mono- and di-residue composition",                                       2 },
  { "--markov0",     eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "generate with 0th order Markov properties per input",                                       2 },
  { "--markov1",     eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "generate with 1st order Markov properties per input",                                       2 },
  { "--reverse",     eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "reverse each input",                                                                        2 },
  { "--iid",         eslARG_NONE,       FALSE,   NULL,      NULL, SHUF_OPTS,  NULL, NULL,                "generate random iid sequence for negatives",                                                2 },
  /* other options */
  { "--fixtime",      eslARG_REAL,       NULL,   NULL,     "x>=0",   NULL,    NULL,"--mintime,--maxtime","TRUE for using a fix time for the ehmm model",          0 },
  { "--mintime",      eslARG_REAL,       NULL,   NULL,     "x>=0",   NULL,    NULL, "--fixtime",         "TRUE for using a range of ehmm models",                 0 },
  { "--maxtime",      eslARG_REAL,     "10.0",   NULL,     "x>=0",   NULL,    NULL, "--fixtime",         "maxtime for using a range of ehmms",                    0 },
  { "--inctime",      eslARG_REAL,     "0.01",   NULL,     "x>=0",   NULL,    NULL, "--fixtime",         "time increment for the ehmm",                           0 },
  { "--minDPL",       eslARG_INT,       "100",   NULL,      NULL,      NULL,  NULL,  NULL,               "minimum segment length for DP shuffling",                                                   4 },
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,       "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <tbl> <trainmsa> <test_msa>";
static char banner[] = "proftransmark";

static int process_dbfile(struct cfg_s *cfg);
static int run_hmmbuild(struct cfg_s *cfg, ESL_MSA *msa, P7_HMM **ret_hmm);
static int ntestmsa(char *tmsafile, char *qmsaname, char *errbuf);
static int read_testmsa(struct cfg_s *cfg, char *tmsafile, char *msaname, ESL_MSA **ret_msa);
static int doctor_testmsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa, ESL_MSA **ret_rmsa);
static int select_domain(struct cfg_s *cfg, ESL_MSA *msa, int *ret_ldom, int *ret_L1, int *ret_L2);
static int domain_notallgaps(ESL_ALPHABET *abc, ESL_DSQ *dsq, int L);
static int truncate_testmsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, int *ret_L1, int *ret_L2);
static int expand_partialmsas(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa1, ESL_MSA **omsa2, int L1, int L2);
static int set_gaps(ESL_ALPHABET *abc, ESL_DSQ *dsq, int L);
static int set_random_segment(ESL_GETOPTS *go, struct cfg_s *cfg, FILE *logfp, ESL_DSQ *dsq, int L);
static int run_ehmmsearch(struct cfg_s *cfg, FILE *benchfp, float time, P7_HMM *hmm, ESL_MSA *rmsa, MSA_STAT *rmsastat, int domax);
static int domain_dump(struct domain_s *dom, int ndom);
static int parse_edomtbl(char *tblfile, ESL_MSA *msa, struct domain_s **ret_dom, int *ret_ndom, char *errbuf);
static int domain2domain_msa(struct cfg_s *cfg, char *hmmfile, FILE *benchfp, struct domain_s *dom, int ndom, ESL_MSA **omsa, ESL_MSA *rmsa, MSA_STAT *rmsastat, 
			     int domax, char *errbuf);
static int expand_msa(ESL_MSA **omsa, struct domain_s *dom, int alen);
static int swap_msa(ESL_MSA **omsa);
static int ehmm_stats(struct cfg_s *cfg, char *hmmfile, float time, char **ret_hmmname, int *ret_hmmALEN, float *ret_hmmSQLEN, float *ret_hmmPID, 
		      float *ret_hmmME, float *ret_hmmMRE);
static int parse_ehmmemit(char *file, char **ret_hmmname, float *ret_hmmME, float *ret_hmmMRE, char *errbuf);
static int parse_alistat(char *file, int *ret_hmmALEN, float *ret_hmmSQLEN, float *ret_hmmPID, char *errbuf);
static int run_benchmark(struct cfg_s *cfg, char *hmmfile, FILE *benchfp, ESL_MSA *rmsa, MSA_STAT *mrstat, ESL_MSA *emsa, MSA_STAT *mestat, float avgsc, float avgtime);

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

  cfg.qmsafile = NULL;
  cfg.tmsafile = NULL;
  if (esl_opt_ArgNumber(go) != 3) { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  if ((cfg.tblfile  = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));
  if ((cfg.qmsafile  = esl_opt_GetArg(go, 2)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.tmsafile  = esl_opt_GetArg(go, 3)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));

  /* outheader for all output files */
  cfg.outheader = NULL;
  msamanip_OutfileHeader(cfg.qmsafile, &cfg.outheader); 
  
   /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nqmsa = 0;
  cfg.ntmsa = 0;

  /* alphabet */
  if      (esl_opt_GetBoolean(go, "--amino"))   cfg.abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg.abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg.abc = esl_alphabet_Create(eslRNA);
  else                                          cfg.abc = NULL;
 cfg.bg = p7_bg_Create(cfg.abc);

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
  cfg.partial    = esl_opt_IsOn(go, "--partial")? TRUE : FALSE;
  cfg.domax      = esl_opt_GetBoolean(go, "--domax");
  cfg.fixtime    = esl_opt_IsOn(go, "--fixtime")? esl_opt_GetReal(go, "--fixtime") : -1.0; 
  cfg.mintime    = esl_opt_IsOn(go, "--mintime")? esl_opt_GetReal(go, "--mintime") : -1.0; 
  cfg.maxtime    = esl_opt_GetReal(go, "--maxtime"); 
  cfg.inctime    = esl_opt_GetReal(go, "--inctime"); 
  cfg.dbfile     = NULL;
  cfg.dbfp       = NULL;
  cfg.dbfmt      = eslSQFILE_FASTA;

  if (esl_opt_IsOn(go, "--mintime")) {
    if (cfg.maxtime <= cfg.mintime)              esl_fatal("bad range of times");
    if (cfg.mintime + cfg.inctime > cfg.maxtime) esl_fatal("bad range of times");
  }

  if (esl_opt_IsOn(go, "--partial")) {
    cfg.dbfile = esl_opt_GetString(go, "--partial");
    
    if (cfg.abc->type == eslAMINO) esl_composition_SW34(cfg.fq);
    else                           esl_vec_DSet(cfg.fq, cfg.abc->K, 1.0 / (double) cfg.abc->K);
    
    /* Open and process the dbfile; make sure it's in the same alphabet */
    process_dbfile(&cfg);
  }

  cfg.Eval    = 2000.;
  cfg.incdomE = 2000.;
  cfg.incE    = 2000.;

  /* open outfile with benchmarking results */
  cfg.benchfile1 = NULL;
  cfg.benchfile2 = NULL;
  cfg.benchfp1   = NULL;
  cfg.benchfp2   = NULL;
  esl_sprintf(&cfg.benchfile1, "%s.bench1", cfg.outheader); 
  esl_sprintf(&cfg.benchfile2, "%s.bench2", cfg.outheader); 
  if ((cfg.benchfp1 = fopen(cfg.benchfile1, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.benchfile1);
  if ((cfg.benchfp2 = fopen(cfg.benchfile2, "w")) == NULL) esl_fatal("Failed to open output file %s", cfg.benchfile2);
  fprintf(cfg.benchfp1, "#MSA_NAME           TP_hom   TRUE_hom  FOUND_hom      \n");
  fprintf(cfg.benchfp2, "#MSA_NAME           TP_hom   TRUE_hom  FOUND_hom      \n");

    
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
  char           *tmsaname = NULL;
  ESL_MSAFILE   *afp      = NULL;
  ESL_MSA        *qmsa     = NULL;       /* the query alignment             */
  ESL_MSA        *tmsa     = NULL;       /* the given pairwise testmsa      */
  ESL_MSA        *trmsa    = NULL;       /* the referenc pairwise testmsa   */
  P7_HMM         *hmm      = NULL;
  float           time;
  float           inctime;
  int             ntmsa;
  int             x;
  int             status = eslOK;
  int             hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* Open the MSA file */
  status = esl_msafile_Open(&(cfg.abc), cfg.qmsafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  /* read the query MSA */
  while ((hstatus = esl_msafile_Read(afp, &qmsa)) != eslEOF) {
    if (hstatus != eslOK) esl_msafile_ReadFailure(afp, status);
    cfg.nqmsa ++;
 
    /* number of tmsas associated to this qmsa */
    ntmsa = ntestmsa(cfg.tblfile, qmsa->name, cfg.errbuf);
    if (ntmsa < 0) esl_fatal("%s Failed to find any test sequences", cfg.errbuf);
    printf("%d>qmsa %s ntest %d\n", cfg.nqmsa, qmsa->name, ntmsa);
    
    /* outheader for all msa-output files */
    msamanip_OutfileHeader((qmsa->acc)?qmsa->name:cfg.qmsafile, &cfg.msaheader); 
    
    esl_msa_ConvertDegen2X(qmsa); 
    esl_msa_Hash(qmsa);
    
    /* create the query hmm with hmmbuild */
    status = run_hmmbuild(&cfg, qmsa, &hmm);
    if (status != eslOK || hmm == NULL) esl_fatal("Failed to create hmm form msa %s", qmsa->name);
    
    /* read the test msas */
    for (x = 0; x < ntmsa; x ++) {

      esl_sprintf(&tmsaname, "%s-%d", qmsa->name, x+1);
      read_testmsa(&cfg, cfg.tmsafile, tmsaname, &tmsa); 
      if (tmsa == NULL) esl_fatal("Failed to find test msa %s",tmsaname);
      if (cfg.voutput) {/* the testmsa */
	fprintf(stdout, "\nORIGINAL TEST alignment\n");
	esl_msafile_Write(stdout, tmsa, eslMSAFILE_STOCKHOLM);
      }

      status = doctor_testmsa(go, &cfg, &tmsa, &trmsa);
      if (status != eslOK) esl_fatal("%s\nFailed to doctor testmsa", cfg.errbuf);
 
      if (cfg.voutput) {/* the testmsa */
	fprintf(stdout, "\nTEST alignment\n");
	esl_msafile_Write(stdout, tmsa, eslMSAFILE_STOCKHOLM);
      }
      if (cfg.voutput) {/* the MSAProbs alignment */
	fprintf(stdout, "\nMSAProbs alignment\n");
	esl_msafile_Write(stdout, trmsa, eslMSAFILE_STOCKHOLM);
      }
      
      /* fixed-time  transitive pairwise alignment of test sequences */
      if (cfg.fixtime > 0.) {
	status = run_ehmmsearch(&cfg, cfg.benchfp1, cfg.fixtime, hmm, trmsa, cfg.trmstat, cfg.domax);
	if (status != eslOK) esl_fatal("%s\nFailed to run fixtime ehmmsearch", cfg.errbuf);
      }
      else if (cfg.mintime > 0.) {
	inctime = cfg.inctime;
	for (time = cfg.mintime; time <= cfg.maxtime; time += inctime) {
	  inctime = (time < 1.0)? cfg.inctime : (time < 10.0)? 10.0*cfg.inctime : 100.0*cfg.inctime;
	  status = run_ehmmsearch(&cfg, cfg.benchfp1, time, hmm, trmsa, cfg.trmstat, cfg.domax);
	  if (status != eslOK) esl_fatal("%s\nFailed to run ehmmsearch at time %f", cfg.errbuf, time);
	}
      }

      if (cfg.fixtime > 0. || cfg.mintime > 0.) {
	/* optimal-time transitive pairwise alignment of test sequences */
	status = run_ehmmsearch(&cfg, cfg.benchfp2, -1.0, hmm, trmsa, cfg.trmstat, cfg.domax);
	if (status != eslOK) esl_fatal("%s\nFailed to run optimal ehmmsearch", cfg.errbuf);
      }
      
      esl_msa_Destroy(tmsa);     tmsa  = NULL;
      esl_msa_Destroy(trmsa);    trmsa = NULL;
      free(cfg.trmsafrq); cfg.trmsafrq = NULL;
    }

    esl_msa_Destroy(qmsa);  qmsa = NULL;
    p7_hmm_Destroy(hmm);    hmm  = NULL;
  }
  
  esl_stopwatch_Destroy(cfg.w); 
  p7_bg_Destroy(cfg.bg);
   
  if (cfg.dbfp) esl_sqfile_Close(cfg.dbfp); 
  fclose(cfg.outfp);
  free(cfg.outheader);
  free(cfg.msaheader);
  fclose(cfg.benchfp1);
  fclose(cfg.benchfp2);
  free(cfg.benchfile1);
  free(cfg.benchfile2);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  esl_msafile_Close(afp);
  if (hmm) p7_hmm_Destroy(hmm);
  if (qmsa) esl_msa_Destroy(qmsa); 
 return 0;
}

/* Open the source sequence database for negative subseqs;
 * upon return, cfg->dbfp is open (digital, SSI indexed);
 * cfg->db_maxL and cfg->db_nseq are set.
 */
static int
process_dbfile(struct cfg_s *cfg)
{
  ESL_SQ     *sq = esl_sq_CreateDigital(cfg->abc);
  int         status;
  
  /* Open the sequence file in digital mode */
  status = esl_sqfile_OpenDigital(cfg->abc, cfg->dbfile, cfg->dbfmt, NULL, &(cfg->dbfp));
  if      (status == eslENOTFOUND) esl_fatal("No such file %s", cfg->dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of seqfile %s unrecognized.", cfg->dbfile);
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
ntestmsa(char *tblfile, char *qmsaname, char *errbuf)
{
  ESL_FILEPARSER      *tblfp = NULL;
  char                *tok1;
  char                *tok2;
  int                  ntest = -1;
  int                  status;
  
  if (esl_fileparser_Open(tblfile, NULL, &tblfp) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to open file %s", tblfile);
  esl_fileparser_SetCommentChar(tblfp, '#');
  
  while (esl_fileparser_NextLine(tblfp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok1, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse name from file %s", tblfile);

      if (strcmp(tok1, qmsaname) == 0) {
	if (esl_fileparser_GetTokenOnLine(tblfp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse pid from file %s", tblfile);	       
	if (esl_fileparser_GetTokenOnLine(tblfp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse pid from file %s", tblfile);	       
	if (esl_fileparser_GetTokenOnLine(tblfp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse pid from file %s", tblfile);	       
	if (esl_fileparser_GetTokenOnLine(tblfp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse pid from file %s", tblfile);	       
	if (esl_fileparser_GetTokenOnLine(tblfp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse pid from file %s", tblfile);
	if (esl_fileparser_GetTokenOnLine(tblfp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse pid from file %s", tblfile);
	ntest = atoi(tok2);
	if (esl_fileparser_GetTokenOnLine(tblfp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse pid from file %s", tblfile);	       
      } 
    }    
  esl_fileparser_Close(tblfp); tblfp = NULL;
  
  return ntest;

 ERROR:
  return -1;
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
  char         *s = NULL;
  int           status;
  
  /* MSA input in stockholm format */
  if ((status = esl_tmpfile_named(tmpmsafile,  &fp))                          != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create ms file");
  if ((status = esl_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_STOCKHOLM)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "Failed to write STOCKHOLM file\n");
  fclose(fp);
  
  /* run hmmbuild */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))                     != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create hmm file");
  esl_sprintf(&args, "%s/lib/hmmer4/src/programs/hmmbuild %s %s>/dev/null", s, tmphmmfile, tmpmsafile);
  system(args);
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
  if (args != NULL) free(args);
  return eslOK;
  
 ERROR:
  remove(tmphmmfile);
  remove(tmpmsafile);
  if (hfp) p7_hmmfile_Close(hfp);
  if (args) free(args);
  if (hmm) p7_hmm_Destroy(hmm);
  return status;
}

static int
read_testmsa(struct cfg_s *cfg, char *tmsafile, char *name, ESL_MSA **ret_msa)
{
  char          tmpmsafile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpnamefile[16] = "esltmpXXXXXX"; /* tmpfile template */
  ESL_MSAFILE *afp = NULL;
  FILE         *fp  = NULL;
  ESL_MSA      *msa = NULL;
  char         *args = NULL;
  char         *s = NULL;
  int           status;
  
  /* name file */
  if ((status = esl_tmpfile_named(tmpnamefile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create namemsa file");
  fprintf(fp, "%s\n", name);
  fclose(fp);

  /* run esl-afetch */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpmsafile, &fp))                     != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create msa file");
  esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-afetch -f %s %s > %s", s, tmsafile, tmpnamefile, tmpmsafile);
  system(args);
  fclose(fp);
  
  /* read MSA from tempmsafile */
  if (esl_msafile_Open(NULL, tmpmsafile, NULL, eslMSAFILE_STOCKHOLM, NULL, &afp) != eslOK) esl_msafile_OpenFailure(afp, status);
  afp->format = eslMSAFILE_STOCKHOLM;
  if (esl_msafile_Read(afp, &msa) != eslOK) esl_msafile_ReadFailure(afp, status);
  esl_msafile_Close(afp);

  *ret_msa = msa;
  
  remove(tmpmsafile); 
  remove(tmpnamefile); 
  
  free(args);
  return eslOK;

 ERROR:
  remove(tmpmsafile);
  remove(tmpnamefile); 

  if (afp) esl_msafile_Close(afp);
  if (args) free(args);
  if (msa) esl_msa_Destroy(msa);
  return status;
}

static int
doctor_testmsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **otmsa, ESL_MSA **ret_trmsa)
{
  ESL_MSA *trmsa = NULL;
  ESL_MSA *tmsa;
  int      L1, L2;
  int      status;

  tmsa = *otmsa;

  /* first truncate and get the reference alignment and stats */
  if (cfg->truncate || cfg->partial) {
    status = truncate_testmsa(go, cfg, tmsa, &L1, &L2);
    if (status != eslOK) goto ERROR;
  }

  /* The msaprobs alignment of the test sequences (before add any non-homologous regions) */
  status = MSAProbs_Align(tmsa, &trmsa, cfg->errbuf, cfg->verbose);
  esl_msa_SetName(trmsa, tmsa->name, -1);
  if (status != eslOK || trmsa == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "Failed to run msaprobs");
    
  /* then (if partial) add the non-homologous regions both to the tmsa and rmsa */
  if (cfg->partial) {
    expand_partialmsas(go, cfg, &tmsa, &trmsa, L1, L2);
  }
  
  /* reference msa aveid and avematch */
  msamanip_CStats(cfg->abc, trmsa, &cfg->trmstat);
  msamanip_CBaseComp(cfg->abc, trmsa, cfg->bg->f, &cfg->trmsafrq);

  *otmsa = tmsa;
  *ret_trmsa = trmsa;
  return eslOK;
  
 ERROR:
  if (trmsa) esl_msa_Destroy(trmsa);
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
truncate_testmsa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, int *ret_L1, int *ret_L2)
{
  ESL_DSQ *dsq;
  int      ldom;
  int      L1, L2;
  int      n;
  int      status;

 esl_msa_Digitize(cfg->abc, msa, NULL);
 
 select_domain(cfg, msa, &ldom, &L1, &L2);
 //printf("ldom %d L1 %d L2 %d alen %"PRId64"\n", ldom, L1, L2, msa->alen);
 if (ldom < 0 || L1 < 0 || L2 < 0) ESL_XFAIL(eslFAIL, cfg->errbuf, "select_domain() failed");

  /* for all the seq in the alignment */
  for (n = 0; n < msa->nseq; n ++) {
    dsq = msa->ax[n];
    set_gaps(cfg->abc, dsq+1,         L1);
    set_gaps(cfg->abc, dsq+L1+ldom+1, L2);
  }
  esl_msa_MinimGaps(msa, NULL, "-", FALSE);

  if (ret_L1) *ret_L1 = L1;
  if (ret_L2) *ret_L2 = L2;
  return eslOK;

 ERROR:
  return status;
}

static int
expand_partialmsas(ESL_GETOPTS *go, struct cfg_s *cfg,  ESL_MSA **omsa1, ESL_MSA **omsa2, int L1, int L2)
{
  ESL_MSA  *new1 = NULL;
  ESL_MSA  *new2 = NULL;
  ESL_MSA  *msa1;
  ESL_MSA  *msa2;
  ESL_SQ  **sql = NULL;
  ESL_SQ  **sqr = NULL;
  int       nseq;
  int       n;
  int       x;
  int       i;
  int       status;

  msa1 = *omsa1;
  msa2 = *omsa2;
  if (msa1->nseq != msa2->nseq) { status = eslFAIL; goto ERROR; }
  nseq = msa1->nseq;

  if (msa1->aseq == NULL) esl_msa_Textize(msa1);
  if (msa2->aseq == NULL) esl_msa_Textize(msa2);

  // create the flanking regions
  ESL_ALLOC(sql, sizeof(ESL_SQ) * nseq);
  ESL_ALLOC(sqr, sizeof(ESL_SQ) * nseq);
  for (n = 0; n < nseq; n ++) {
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
  
  // then add the flanking regions 
  for (x = 0; x < nseq; x ++) {
    new1 = esl_msa_Create(nseq, msa1->alen+L1+L2);  
    esl_msa_SetName(new1, msa1->name, -1);

    for (n = 0; n < nseq; n ++) {
      esl_msa_SetSeqName(new1, n, msa1->sqname[n], -1);
     
      for (i = 0; i < new1->alen; i ++) {	  
       if (n==x) {
	 if      (i <  L1)            new1->aseq[x][i] = tolower(sql[x]->seq[i]);
	 else if (i >= L1+msa1->alen) new1->aseq[x][i] = tolower(sqr[x]->seq[i-(L1+msa1->alen)]);
	 else                         new1->aseq[x][i] = msa1->aseq[x][i-L1];
       }
       else {
	 if (i < L1 || i >= L1+msa1->alen) new1->aseq[n][i] = '-';
	 else                              new1->aseq[n][i] = msa1->aseq[n][i-L1];
       }
      }
   }
    
    esl_msa_Destroy(msa1); msa1 = NULL;
    msa1 = new1;
  }
  
  for (x = 0; x < nseq; x ++) {
   new2 = esl_msa_Create(nseq, msa2->alen+L1+L2);  
   esl_msa_SetName(new2, msa2->name, -1);

    for (n = 0; n < nseq; n ++) {
      esl_msa_SetSeqName(new2, n, msa2->sqname[n], -1);

      for (i = 0; i < new2->alen; i ++) {	  
	if (n==x) {
	  if      (i <  L1)            new2->aseq[x][i] = tolower(sql[x]->seq[i]);
	  else if (i >= L1+msa2->alen) new2->aseq[x][i] = tolower(sqr[x]->seq[i-(L1+msa2->alen)]);
	  else                         new2->aseq[x][i] = msa2->aseq[x][i-L1];
	}
	else {
	  if (i < L1 || i >= L1+msa2->alen) new2->aseq[n][i] = '-';
	  else                              new2->aseq[n][i] = msa2->aseq[n][i-L1];
	}
      }
    }

    esl_msa_Destroy(msa2); msa2 = NULL;
    msa2 = new2;
   }
  
  if (sql) {
    for (n = 0; n < nseq; n ++) 
      if (sql[n]) esl_sq_Destroy(sql[n]);
    free(sql);
  }
  if (sqr) {
    for (n = 0; n < nseq; n ++) 
      if (sqr[n]) esl_sq_Destroy(sqr[n]);
    free(sqr);
  }

  *omsa1 = msa1;
  *omsa2 = msa2;
  
  //printf("O1 %d\n", x);
  //esl_msafile_Write(stdout, *omsa1, eslMSAFILE_STOCKHOLM);
  //printf("O2 %d\n", x);
  //esl_msafile_Write(stdout, *omsa2, eslMSAFILE_STOCKHOLM);
  
  return eslOK;

 ERROR:
  if (sql) {
    for (n = 0; n < nseq; n ++) 
      if (sql[n]) esl_sq_Destroy(sql[n]);
    free(sql);
  }
  if (sqr) {
    for (n = 0; n < nseq; n ++) 
      if (sqr[n]) esl_sq_Destroy(sqr[n]);
    free(sqr);
  }
  if (new1) esl_msa_Destroy(new1);
  if (new2) esl_msa_Destroy(new2);
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
run_ehmmsearch(struct cfg_s *cfg, FILE *benchfp, float time, P7_HMM *hmm, ESL_MSA *rmsa, MSA_STAT *rmsastat, int domax)
{
  char             tmprmsafile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpmsafile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpfastafile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char             tmptblfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmphmmfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  ESL_MSAFILE    *afp = NULL;
  FILE            *fp  = NULL;
  ESL_MSA         *msa = NULL;
  struct domain_s *dom;
  int              ndom;
  char            *args = NULL;
  char            *s = NULL;
  int              status;

  /* the rmsa converted to tmpfile */
  if ((status = esl_tmpfile_named(tmprmsafile,  &fp))                          != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create rmsa file");
  if ((status = esl_msafile_Write(fp, (ESL_MSA *)rmsa, eslMSAFILE_STOCKHOLM)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "Failed to write STOCKHOLM file\n");
  fclose(fp);

  /* run esl-reformat */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpfastafile, &fp))                     != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create msa file");
  esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-reformat fasta %s > %s", s, tmprmsafile, tmpfastafile);
  system(args);
  fclose(fp);

  /* the hmm tmpfile */
  if ((status = esl_tmpfile_named(tmphmmfile, &fp))  != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tmphmm file");
  if ((status = p7_hmmfile_WriteASCII(fp, -1, hmm))  != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to write hmm file");
  fclose(fp);
  
  /* run ehmmsearch */
  if ((status = esl_tmpfile_named(tmpmsafile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create msa file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmptblfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tbl file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
  fclose(fp);
  if (time < 0.) 
    esl_sprintf(&args, "%s/src/programs/ehmmsearch              --cpu 1 -E %f --incdomE %f --incE %f --max --domtbl %s -A %s %s %s > %s", 
		s,       cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmpmsafile, tmphmmfile, tmpfastafile, tmpoutfile);  
  else
    esl_sprintf(&args, "%s/src/programs/ehmmsearch --fixtime %f --cpu 1 -E %f --incdomE %f --incE %f --max --domtbl %s -A %s %s %s > %s", 
		s, time, cfg->Eval, cfg->incdomE, cfg->incE, tmptblfile, tmpmsafile, tmphmmfile, tmpfastafile, tmpoutfile);
  system(args);

  /* the msa */
  if (esl_msafile_Open(NULL, tmpmsafile, NULL, eslMSAFILE_STOCKHOLM, NULL, &afp) != eslOK) esl_msafile_OpenFailure(afp, status);
  afp->format = eslMSAFILE_STOCKHOLM;
  if (esl_msafile_Read(afp, &msa) != eslOK) esl_msafile_ReadFailure(afp, status);
  esl_msafile_Close(afp); afp = NULL;
  if (msa == NULL) ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to create the msa");

  esl_msa_MinimGaps(msa, NULL, "-.~", FALSE);
  if (cfg->voutput) esl_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM);

  status = parse_edomtbl(tmptblfile, rmsa, &dom, &ndom, cfg->errbuf); if (status != eslOK) goto ERROR;
  if (cfg->voutput) domain_dump(dom, ndom);

  status = domain2domain_msa(cfg, tmphmmfile, benchfp, dom, ndom, &msa, rmsa, rmsastat, domax, cfg->errbuf); if (status != eslOK) goto ERROR;
  
  if (dom) free(dom);
  if (msa) esl_msa_Destroy(msa);

  remove(tmphmmfile);
  remove(tmpfastafile);
  remove(tmprmsafile);
  remove(tmpmsafile);
  remove(tmptblfile);
  remove(tmpoutfile);
  free(args);
  return eslOK;

 ERROR:
  remove(tmphmmfile);
  remove(tmpfastafile);
  remove(tmprmsafile);
  remove(tmpmsafile);
  remove(tmptblfile);
  remove(tmpoutfile);

  if (args) free(args);
  if (msa)  esl_msa_Destroy(msa);  
  if (afp) esl_msafile_Close(afp);
  if (dom) free(dom);

  return status;
}

static int
parse_edomtbl(char *tblfile, ESL_MSA *msa, struct domain_s **ret_dom, int *ret_ndom, char *errbuf)
{
  ESL_FILEPARSER  *tblfp  = NULL;
  struct domain_s *domain = NULL;
  struct domain_s *dom;
  ESL_SQ          *sq     = NULL;
  char            *tok;
  char            *name, *prvname;
  int              ndom = 0;
  int              idx = 0;
  int              L1, L2;
  int              n;
  int              i;
  int              status;

  if (esl_fileparser_Open(tblfile, NULL, &tblfp) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to open file %s", tblfile);
  esl_fileparser_SetCommentChar(tblfp, '#');

  ESL_ALLOC(domain, sizeof(struct domain_s) * (ndom+1));

  while (esl_fileparser_NextLine(tblfp) == eslOK)
    {
      dom = &domain[ndom];
      dom->time = -1.0;

      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse sqname    from file %s", tblfile);
      esl_strdup(tok, -1, &dom->sqname);
      esl_strdup(tok, -1, &name);

      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse accesion  from file %s", tblfile);	  
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse tlen      from file %s", tblfile);
      dom->sqlen = atoi(tok);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse queryname from file %s", tblfile);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse accesion  from file %s", tblfile);	  
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse qlen      from file %s", tblfile);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse Evalue    from file %s", tblfile);
      dom->E = atof(tok);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse score     from file %s", tblfile);
      dom->sc = atof(tok);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse bias      from file %s", tblfile);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse time      from file %s", tblfile);
      dom->time = atof(tok);

      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse #         from file %s", tblfile);	       
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse of        from file %s", tblfile);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse cEvalue   from file %s", tblfile);	       
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse iEvalue   from file %s", tblfile);
      dom->cE = atof(tok);	       
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse score     from file %s", tblfile);
      dom->csc = atof(tok);	       
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse bias      from file %s", tblfile);
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse hmmfrom   from file %s", tblfile);	       
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse hmmto     from file %s", tblfile);	       
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse alifrom   from file %s", tblfile);
      dom->sqfrom = atoi(tok) - 1;
      if (esl_fileparser_GetTokenOnLine(tblfp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse alito     from file %s", tblfile);	     
      dom->sqto = atoi(tok) - 1;
      esl_sprintf(&dom->name, "%s/%d-%d", dom->sqname, dom->sqfrom+1, dom->sqto+1);     

      dom->sql = NULL;
      dom->sqr = NULL;
      dom->len  = dom->sqto - dom->sqfrom + 1;
      if (strcmp(name,prvname) != 0) idx ++;
      dom->sqidx = idx;

      /* identify the unaligned flanquing regions */
      for (n = 0; n < msa->nseq; n ++) {
	if (esl_strcmp(msa->sqname[n], dom->sqname) == 0) {
	  //printf("found sq %s from %d to %d | sqlen %d\n", dom->sqname, dom->sqfrom, dom->sqto, dom->sqlen); 
	  break;
	}
      }
      if (n == msa->nseq) ESL_XFAIL(eslFAIL, errbuf, "could not find sq %s\n", dom->sqname);
      
      status = esl_sq_FetchFromMSA(msa, n, &sq); /* extract the seq from the rmsa */
      if (status != eslOK) ESL_XFAIL(status, errbuf, "could not extract sq %s\n", dom->sqname);
      if (dom->sqlen != sq->n) ESL_XFAIL(eslFAIL, errbuf, "wrong lenght for sq %s is %d sqlen should be %" PRId64"\n", dom->sqname, dom->sqlen, sq->n); 
      
      L1 = dom->sqfrom;
      L2 = dom->sqlen - 1 - dom->sqto;
      dom->sql = esl_sq_Create();
      dom->sqr = esl_sq_Create();
      esl_sq_GrowTo(dom->sql, L1+1);
      esl_sq_GrowTo(dom->sqr, L2+1);
      dom->sql->n = L1;
      dom->sqr->n = L2;
      for (i = 0; i < L1; i ++) 
	dom->sql->seq[i] = tolower(sq->seq[i]);
      dom->sql->seq[L1] = '\0';
      for (i = 0; i < L2; i ++) 
	dom->sqr->seq[i] = tolower(sq->seq[dom->sqto+i+1]);
      dom->sqr->seq[L2] = '\0';
      prvname = name;
      ndom ++;
      if (domain) ESL_REALLOC(domain, sizeof(struct domain_s) * (ndom+1));

    }
  esl_fileparser_Close(tblfp); tblfp = NULL;

  *ret_dom  = domain;
  *ret_ndom = ndom;
  return eslOK;

  ERROR:
  if (domain) free(domain);
  return status;
}

static int
domain_dump(struct domain_s *dom, int ndom) 
{
  int d;

  printf("ndom %d\n", ndom);
  for (d = 0; d < ndom; d ++) {
    printf("dom[%d]\t%d-%d\tlen %d\tcsc %.2f\tcE %f\ttime %f\t", d, dom[d].sqfrom, dom[d].sqto, dom[d].len, dom[d].csc, dom[d].cE, dom[d].time);
    printf("sq %s len %d\n", dom[d].sqname, dom[d].sqlen); 
    printf("sql[%"PRId64"]=%s\n", dom[d].sql->n, dom[d].sql->seq);
    printf("sqr[%"PRId64"]=%s\n", dom[d].sqr->n, dom[d].sqr->seq);
  }

  return eslOK;
}

static int
domain2domain_msa(struct cfg_s *cfg, char *hmmfile, FILE *benchfp, struct domain_s *domain, int ndom, ESL_MSA **omsa, ESL_MSA *rmsa, MSA_STAT *rmsastat, 
		  int domax, char *errbuf)
{
  struct domain_s *dom1, *dom2;
  ESL_MSA         *dmsa = NULL;
  ESL_MSA         *msa;
  MSA_STAT        *dmsastat = NULL;
  float            avgsc;
  float            avgtime;
  float            sc1max = -eslINFINITY;
  float            sc2max = -eslINFINITY;
  int             *useme = NULL;
  int              d, d1, d2;
  int              d1max = -1;
  int              d2max = -1;
  int              n;
  int              status;
  
  msa = *omsa;

  ESL_ALLOC(useme, sizeof(int)*ndom);
  
  if (domax) {
    for (d = 0; d < ndom; d ++) {
      if      (domain[d].sqidx == 1 && domain[d].csc > sc1max) { sc1max = domain[d].csc; d1max = d; }
      else if (domain[d].sqidx == 2 && domain[d].csc > sc2max) { sc2max = domain[d].csc; d2max = d; }
    }
  }
  
  for (d1 = 0; d1 < ndom; d1 ++) {
    dom1 = &domain[d1];
    if (domax && d1 != d1max) continue;

    for (d2 = d1+1; d2 < ndom; d2 ++) {
      dom2 = &domain[d2];
      
      if (dom1->sqidx == dom2->sqidx) continue;
      if (domax && d2 != d2max) continue;
    
      avgsc   = 0.5 * (dom1->sc + dom2->sc);
      avgtime = 0.5 * (dom1->time + dom2->time);

      esl_vec_ISet(useme, ndom, FALSE);
      useme[d1] = TRUE;
      useme[d2] = TRUE;
      if ((status = esl_msa_SequenceSubset(msa, useme, &dmsa)) != eslOK) goto ERROR;
      esl_msa_MinimGaps(dmsa, errbuf, "-_.~", FALSE);
      for (n = 0; n < dmsa->nseq; n ++) {
	if      (strcmp(dmsa->sqname[n], dom1->name) == 0) esl_msa_SetSeqName(dmsa, n, dom1->sqname, -1);
	else if (strcmp(dmsa->sqname[n], dom2->name) == 0) esl_msa_SetSeqName(dmsa, n, dom2->sqname, -1);
	else ESL_XFAIL(eslFAIL, cfg->errbuf, "failed to find domain for seq %s", dmsa->sqname[n]);
      }

      expand_msa(&dmsa, dom1, dmsa->alen);
      expand_msa(&dmsa, dom2, dmsa->alen);
      if (cfg->voutput) {
	printf("\nalignment dom %d-%d\n", d1, d2);
	esl_msafile_Write(stdout, dmsa, eslMSAFILE_STOCKHOLM);
      }
   
      msamanip_CStats(cfg->abc, dmsa, &dmsastat);
      
      if (strcmp(rmsa->sqname[0], dmsa->sqname[0]) != 0) {
	//printf("\n need to swap sqs in the dmsa\n");
	status = swap_msa(&dmsa);
	if (status != eslOK) goto ERROR;
	if (strcmp(rmsa->sqname[0], dmsa->sqname[0]) != 0) { status = eslFAIL; goto ERROR; }
      }
     
      status = run_benchmark(cfg, hmmfile, benchfp, rmsa, rmsastat, dmsa, dmsastat, avgsc, avgtime);
      if (status != eslOK) goto ERROR;

      esl_msa_Destroy(dmsa); dmsa = NULL;
    }
  }
 
  if (useme) free(useme);
  return eslOK;
  
 ERROR:
  if (useme) free(useme);
  
  return status;
}

static int
swap_msa(ESL_MSA **omsa)
{
  ESL_MSA *new = NULL;
  ESL_MSA *msa;
  int      n;
  int      nn;
  int      i;
  int      status;

  msa = *omsa;
  if (msa->nseq != 2) { status = eslFAIL; goto ERROR; }

  new = esl_msa_Create(msa->nseq, msa->alen);  
  esl_msa_SetName(new, msa->name, -1);

  for (n = 0; n < new->nseq; n ++) {
    nn = (n == 0)? 1 : 0;
    esl_msa_SetSeqName(new, n, msa->sqname[nn], -1);

    for (i = 0; i < new->alen; i ++) {	
      new->aseq[n][i] = msa->aseq[nn][i];  
    }
  }
  esl_msa_Destroy(msa);
  *omsa = new;

  return eslOK;

 ERROR:
  if (new) esl_msa_Destroy(new);
  return status;

}

static int
expand_msa(ESL_MSA **omsa, struct domain_s *dom, int alen)
{
  ESL_MSA *new = NULL;
  ESL_MSA *msa;
  int      L1, L2;
  int      n;
  int      i;
  int      thisone;
  
  msa = *omsa;

  L1 = dom->sql->n;
  L2 = dom->sqr->n;

  new = esl_msa_Create(msa->nseq, msa->alen+L1+L2);  
  esl_msa_SetName(new, msa->name, -1);

  for (n = 0; n < new->nseq; n ++) {
    esl_msa_SetSeqName(new, n, msa->sqname[n], -1);
    thisone = (esl_strcmp(msa->sqname[n], dom->sqname) == 0)? TRUE : FALSE;

    for (i = 0; i < new->alen; i ++) {	  
      if (thisone) {
	if      (i <  L1)      new->aseq[n][i] = dom->sql->seq[i];
	else if (i >= L1+alen) new->aseq[n][i] = dom->sqr->seq[i-(L1+alen)];
	else                   new->aseq[n][i] = msa->aseq[n][i-L1];
      }
      else {
	if (i < L1 || i >= L1+alen) new->aseq[n][i] = '-';
	else                        new->aseq[n][i] = msa->aseq[n][i-L1];
      }
    }
  }
  //esl_msafile_Write(stdout, new, eslMSAFILE_STOCKHOLM);
 
  esl_msa_Destroy(msa);
  *omsa = new;

  return eslOK;
}

static int
ehmm_stats(struct cfg_s *cfg, char *hmmfile, float time, char **ret_hmmname, int *ret_hmmALEN, float *ret_hmmSQLEN, float *ret_hmmPID, 
	   float *ret_hmmME, float *ret_hmmMRE)
{
  char   tmpmsafile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char   tmpout1file[16] = "esltmpXXXXXX"; /* tmpfile template */
  char   tmpout2file[16] = "esltmpXXXXXX"; /* tmpfile template */
  char  *args = NULL;
  char  *hmmname = NULL;
  char  *s = NULL;
  FILE  *fp = NULL;
  int    hmmALEN;
  float  hmmSQLEN;
  float  hmmPID;
  float  hmmME;
  float  hmmMRE;
  int    N = 500;
  int    status;

  /* the tmpfiles */
  if ((status = esl_tmpfile_named(tmpmsafile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create msa file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpout1file,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
  fclose(fp);
  if ((status = esl_tmpfile_named(tmpout2file,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
  fclose(fp);

  /* run ehmmemit and esl-alistat to get stats on the actual hmm */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&args, "%s/src/programs/ehmmemit --time %f -N %d -a -o %s --statfile %s %s", s, time, N, tmpmsafile, tmpout1file, hmmfile);  
  system(args);

  /* run esl-alistat */
  esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-alistat %s > %s", s, tmpmsafile, tmpout2file);
  system(args);
  fclose(fp);

  /* parse hmmemit output file to get hmmME and hmmMRE */
  parse_ehmmemit(tmpout1file, &hmmname, &hmmME, &hmmMRE, cfg->errbuf);

  /* parse esl-alistat output file to get hmmPID */
  parse_alistat(tmpout2file, &hmmALEN, &hmmSQLEN, &hmmPID, cfg->errbuf);

  if (ret_hmmname)  *ret_hmmname  = hmmname;
  if (ret_hmmME)    *ret_hmmME    = hmmME;
  if (ret_hmmMRE)   *ret_hmmMRE   = hmmMRE;
  if (ret_hmmALEN)  *ret_hmmALEN  = hmmALEN;
  if (ret_hmmSQLEN) *ret_hmmSQLEN = hmmSQLEN;
  if (ret_hmmPID)   *ret_hmmPID   = hmmPID;

  remove(tmpmsafile);
  remove(tmpout1file);
  remove(tmpout2file);
  if (args) free(args);
  return eslOK;

 ERROR:
  remove(tmpmsafile);
  remove(tmpout1file);
  remove(tmpout2file);
  if (args) free(args);
  return status;
}

static int
parse_ehmmemit(char *file, char **ret_hmmname, float *ret_hmmME, float *ret_hmmMRE, char *errbuf)
{
  ESL_FILEPARSER *fp  = NULL;
  char           *hmmname = NULL;
  float           ME  = -1.0;
  float           MRE = -1.0;
  char           *tok;
  float           time;
  int             status;

  if (esl_fileparser_Open(file, NULL, &fp) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to open file %s", file);
  esl_fileparser_SetCommentChar(fp, '#');
  
  while (esl_fileparser_NextLine(fp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      esl_strdup(tok, -1, &hmmname);
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      time = atof(tok);
      if (time == 0.0 || time == 1.0 || time == eslINFINITY) continue;
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      ME = atof(tok);
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
      MRE = atof(tok);
      if (ME < 0. || MRE < 0.) ESL_XFAIL(eslFAIL, errbuf, "error extracting ME or MRE from file %s", file);
    }
  
  *ret_hmmname = hmmname;
  *ret_hmmME   = ME;
  *ret_hmmMRE  = MRE;
  esl_fileparser_Close(fp);
  return eslOK;

 ERROR:
  if (fp) esl_fileparser_Close(fp);
  return status;
}

static int
parse_alistat(char *file, int *ret_hmmALEN, float *ret_hmmSQLEN, float *ret_hmmPID, char *errbuf)
{
  ESL_FILEPARSER *fp  = NULL;
  int             alen     = -1;
  float           avgpid   = -1.0;
  float           avgsqlen = -1.0;
  char           *tok;
  int             nl = 0;
  int             status;
  
  if (esl_fileparser_Open(file, NULL, &fp) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to open file %s", file);
  esl_fileparser_SetCommentChar(fp, '#');
  
  while (esl_fileparser_NextLine(fp) == eslOK)
    {
      nl ++;
      if (nl == 4) {
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);	
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);	
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
	alen = atoi(tok);
      }
      else if (nl == 8) {
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);	
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);	
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
	avgsqlen = atof(tok);
      }
      else if (nl == 9) {
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);	
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);	
	if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse from file %s", file);
	avgpid = atof(tok);
      }
      else continue;
    }
  
  if (alen < 0 || avgsqlen < 0.|| avgpid < 0.) ESL_XFAIL(eslFAIL, errbuf, "failed to extract parameters from file %s\n", file);
  
  if (ret_hmmALEN)  *ret_hmmALEN  = alen;
  if (ret_hmmSQLEN) *ret_hmmSQLEN = avgsqlen;
  if (ret_hmmPID)   *ret_hmmPID   = avgpid;
  esl_fileparser_Close(fp);
  return eslOK;

 ERROR:
  if (fp) esl_fileparser_Close(fp);
  return status;
}

static int 
run_benchmark(struct cfg_s *cfg, char *hmmfile, FILE *benchfp, ESL_MSA *rmsa, MSA_STAT *mrstat, ESL_MSA *emsa, MSA_STAT *mestat, float avgsc, float avgtime)
{
  char   *msaname = NULL;
  double  cputime;
  char   *hmmname = NULL;
  int     hmmALEN  = -1;
  float   hmmSQLEN = -1.0;
  float   hmmPID   = -1.0;
  float   hmmME    = -1.0;
  float   hmmMRE   = -1.0;
  int     lcu_r;  // if TRUE, consider lower case in reference alignment as unaligned
  int     lcu_e;  // if TRUE, consider lower case in inferred alignment as unaligned
  int     status; 
  
  esl_sprintf(&msaname, "%s", rmsa->name); 
 
  /* ehmm stats using ehmmemeit */
  status = ehmm_stats(cfg, hmmfile, avgtime, &hmmname, &hmmALEN, &hmmSQLEN, &hmmPID, &hmmME, &hmmMRE);
  if (status != eslOK) goto ERROR;
  printf("\nehmm %s %f %d %f %f %f %f\n", hmmname, avgtime, hmmALEN, hmmSQLEN, hmmPID, hmmME, hmmMRE);
  fprintf(benchfp, "ehmm %s %f %d %f %f %f %f\n", hmmname, avgtime, hmmALEN, hmmSQLEN, hmmPID, hmmME, hmmMRE);
  
  cputime = cfg->w->user + cfg->w->sys;
  
  lcu_r = FALSE; /* it is all in upper case and aligned */
  lcu_e = TRUE; /* the estimated aligment is created with 
		 * p7_trace_FauxFromMSA() such that the residues in the "a" seq are matches
		 * p7_tracealign_Seqs() and includes lower case
		 * for residues that did not align to the more divergent "a" sequence
		 */
  status = FastSP_Benchmark(benchfp, msaname, NULL, cfg->abc, rmsa, mrstat, emsa, mestat, avgsc, avgtime, cputime, lcu_r, lcu_e, cfg->errbuf, cfg->verbose);
  if (status != eslOK) goto ERROR;
  
  free(hmmname);
  free(msaname);
  return eslOK;

 ERROR:
  if (hmmname) free(hmmname);
  if (msaname) free(msaname);
  return status;
}

