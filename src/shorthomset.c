/* shorthomset -- from a db of pairwise alignments, extract short regions of homology
 *                and add flaking random sequences.
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
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2.h"
#include "e2_msa.h"
#include "e2_train.h"
#include "e2_tree.h"
#include "msamanip.h"
#include "msatree.h"


/* Exclusive options for evolutionary model choice */

/* struct cfg_s : "Global" application configuration shared by all threads/processes.
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  /* Shared configuration in masters & workers */
  int                  argc;
  char               **argv;
  
  char                 errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS      *r;	               /* random numbers for stochastic sampling (grm-emit) */
  ESL_ALPHABET        *abc;                    /* the alphabet */
  double               fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double               idthresh;	       /* fractional identity threshold for selecting subset of sequences */
  
  char                *bmarkname;
  
  char                 summfile[256];	       /* name of summary file          */
  char                *setfile;
  FILE                *setfp;
  FILE                *summfp;                 /* output stream: summary table of set */

  int                  fragone;                 /* fragment only one of the two sequences */
  int                  ndomains;
  double               mu;                      /* mu of gaussian to sample fragments */
  double               sigma;                   /* sigma of gaussian to sample fragments */
  int                  voutput;                 /* TRUE for verbose output */
  char                *outheader;               /* header for all output files */
  char                *msaheader;               /* header for all msa-specific output files */
 
  int                  nmsa;
  char                *msafile;
  ESL_MSA             *msa;
  MSA_STAT            *mstat;                   /* statistics of the individual alignment */

  char                *dbfile;
  int                  dbfmt;
  ESL_SQFILE          *dbfp;   	                /* source database for negatives                           */
  int                  db_nseq;	                /* # of sequences in the db                                */
  int                  db_maxL;          	/* maximum seq length in db_lens                           */

  int                  infmt;
  float                tol;
  int                  verbose;
};

 static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      0 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                0 },
  /* options for msa */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "-I",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < x id",                                                               0 },
  { "--mu",           eslARG_REAL,     "50.",    NULL,      "x>0",   NULL,    NULL,  NULL,               "length of homology extracted",                                                              0 },
  { "--sigma",        eslARG_REAL,     "10.",    NULL,      "x>0",   NULL,    NULL,  NULL,               "length of homology extracted",                                                              0 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE,   FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                                       0 },
  { "--voutput",      eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "verbose output",                                                                            0 },
  /* msa format */ 
  { "--informat",  eslARG_STRING,        NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            0 },
  /* other options */
  { "--single",       eslARG_NONE,     FALSE,    NULL,       NULL,   NULL,    NULL,  NULL,               "embed one, not two domains in each positive",                                               4 },
  { "--fragone",      eslARG_NONE,     FALSE,    NULL,       NULL,   NULL,    NULL,  NULL,               "embed one only in one of the two seqs",                                                     4 },
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",         eslARG_INT,        "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "shorthomset [-options] <benchmark_prefix> <Stockholm MSA file> <FASTA sequence database> ";
static char banner[] = "short homologies flaked by random sequence";

static double avg_pid_msa(ESL_MSA *msa);
static int    length_dealigned(ESL_ALPHABET *abc, ESL_DSQ *dsq, int alen, int idx, int *l);
static int    process_dbfile(struct cfg_s *cfg, char *dbfile, int dbfmt);
static int    set_random_segment(struct cfg_s *cfg, ESL_DSQ *dsq, int L, int idx, int *l);
static int    short_homology(struct cfg_s *cfg, int *ret_nhom);
static int    short_homology_embedding(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int *ret_L1, int *ret_L2, int *ret_L3, int *dd1n, int *dd2n, int *l1, int *l2, int *l3);
static int    msa_annotate_dompid(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int L1, int L2, int L3, double *ret_d1pid, double *ret_d2pid);
static int    msa_annotate_rf(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int L1, int L2, int L3);
static int    msa_annotate_nonhomologs(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int L1, int L2, int L3);
static int    write_msa(FILE *fp, ESL_MSA *msa, int verbose);
static int    write_summary(struct cfg_s *cfg, ESL_MSA *msa, int d1n, double d1pid, int d2n, double d2pid, int L1, int L2, int L3, int *dd1n, int *dd2n, int *l1, int *l2, int *l3);

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
  if (esl_opt_ArgNumber(go) != 3) { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  if ((cfg.bmarkname  = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <bmkname> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.msafile  = esl_opt_GetArg(go, 2)) == NULL) { 
    if (puts("Failed to get <slifile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.dbfile   = esl_opt_GetArg(go, 3)) == NULL) { 
    if (puts("Failed to get <dbfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  cfg.dbfmt    = eslSQFILE_FASTA;
  
   cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));  
  
  /* outheader for all output files */
  cfg.outheader = NULL;
  esl_sprintf(&cfg.outheader, cfg.bmarkname);

  if (snprintf(cfg.summfile, 256, "%s.summary", cfg.bmarkname) >= 256)  esl_fatal("Failed to construct summfile");
   if ((cfg.summfp = fopen(cfg.summfile, "w"))      == NULL) esl_fatal("Failed to open set summary file\n");
   
  /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nmsa = 0;
  cfg.msa  = NULL;
  
  /* alphabet */
  cfg.abc = esl_alphabet_Create(eslAMINO);
  
  /* other options */
  cfg.mu       = esl_opt_GetReal(go, "--mu");
  cfg.sigma    = esl_opt_GetReal(go, "--sigma");
  cfg.fragfrac = esl_opt_IsOn(go, "-F")? esl_opt_GetReal(go, "-F") : -1.0;
  cfg.idthresh = esl_opt_IsOn(go, "-I")? esl_opt_GetReal(go, "-I") : -1.0;
  cfg.tol      = esl_opt_GetReal   (go, "--tol");
  cfg.verbose  = esl_opt_GetBoolean(go, "-v");
  cfg.voutput  = esl_opt_GetBoolean(go, "--voutput");
  
  cfg.ndomains = (esl_opt_GetBoolean(go, "--single")  == TRUE) ? 1 : 2;
  cfg.fragone  = (esl_opt_GetBoolean(go, "--fragone") == TRUE) ? TRUE : FALSE;

  /* open outfile */
  cfg.setfile = NULL;
  cfg.setfp   = NULL;
  esl_sprintf(&cfg.setfile, "%s.set.sto", cfg.outheader);
  if ((cfg.setfp = fopen(cfg.setfile, "w")) == NULL) esl_fatal("Failed to open set file %s", cfg.setfile);
  
  *ret_go  = go;
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
  char           *msg = "shorthomset failed";
  ESL_GETOPTS    *go;
  struct cfg_s    cfg;
  ESL_MSAFILE   *afp = NULL;
  int             seq_cons_len = 0;
  int             nfrags = 0;	  	  /* # of fragments removed */
  int             nremoved = 0;	          /* # of identical sequences removed */
  int             status = eslOK;
  int             hstatus = eslOK;
  int             nhom = 0;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  /* outheader for all msa-output files */
  msamanip_OutfileHeader(cfg.msafile, &cfg.msaheader); 
  
  /* Open and process the dbfile; make sure it's in the same alphabet */
  process_dbfile(&cfg, cfg.dbfile, cfg.dbfmt);
  
  /* Open the MSA file */
  status = esl_msafile_Open(&(cfg.abc), cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  /* read the training MSAs */
  while ((hstatus = esl_msafile_Read(afp, &cfg.msa)) != eslEOF) {
    if (hstatus != eslOK) esl_msafile_ReadFailure(afp, status);
    
    esl_msa_ConvertDegen2X(cfg.msa); 
    esl_msa_Hash(cfg.msa);
    
    if (esl_opt_IsOn(go, "-F") && msamanip_RemoveFragments(cfg.fragfrac, &cfg.msa, &nfrags, &seq_cons_len) != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-I")) msamanip_SelectSubsetBymaxID(cfg.r, &cfg.msa, cfg.idthresh, &nremoved);
    
    /* print some info */
    if (cfg.voutput) {
      msamanip_XStats(cfg.msa, &cfg.mstat); //  msa aveid and avematch 
      fprintf(stdout, "%s %6d sqs \n", cfg.msa->acc, cfg.msa->nseq);
      msamanip_DumpStats(stdout, cfg.msa, cfg.mstat); 
    }       

    cfg.nmsa ++;
    /* extract a pairwise alignment at random */
    short_homology(&cfg, &nhom); 
 
    esl_msa_Destroy(cfg.msa); cfg.msa = NULL;
  }
  esl_msafile_Close(afp);
  
  printf("found %d/%d homologies\n", nhom, cfg.nmsa);

  /* cleanup */          
  fclose(cfg.setfp);
  fclose(cfg.summfp);
  free(cfg.outheader);
  free(cfg.setfile);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);

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
  printf("dbfile: %s nseq %d maxL %d\n", dbfile, cfg->db_maxL, cfg->db_nseq);

  /* Open SSI index */
  if (esl_sqfile_OpenSSI(cfg->dbfp, NULL) != eslOK) esl_fatal("Failed to open SSI index file");
  if (cfg->dbfp->data.ascii.ssi->nprimary != cfg->db_nseq)     esl_fatal("oops, nprimary != nseq");

  esl_sq_Destroy(sq);
  return eslOK;
}

static int
short_homology(struct cfg_s *cfg, int *ret_nhom)
{
  char    *msaname = NULL;
  ESL_MSA *shortmsa = NULL;
  double   d1pid, d2pid;        /* percentate id of domains */
  double   mu, sigma;
  int      N = cfg->msa->nseq;
  int      L = cfg->msa->alen;
  int      d1n, d2n;		/* lengths of two domains             */
  int      L1,L2,L3;		/* lengths of three random regions    */
  int     *dd1n = NULL;		/* lengths of two domains (unaligned) */
  int     *dd2n = NULL;		/* lengths of two domains (unaligned) */
  int     *l1 = NULL;	  	/* lengths of three unaligned random regions    */
  int     *l2 = NULL;	   	/* lengths of three unaligned random regions    */
  int     *l3 = NULL;	     	/* lengths of three unaligned random regions    */
  int      MINHOM;
  int      MAXNIT = 1000;
  int      minhom = 0;
  int      nit = 0;
  int      status;
  
  shortmsa = esl_msa_Clone(cfg->msa); 

  esl_sprintf(&msaname, "%s.%d", cfg->outheader, cfg->nmsa);
  esl_msa_SetName(shortmsa, msaname, -1);

  ESL_ALLOC(dd1n, sizeof(int) * N);
  ESL_ALLOC(l1,   sizeof(int) * N);
  ESL_ALLOC(l2,   sizeof(int) * N);

  if (cfg->ndomains == 2) {
    ESL_ALLOC(dd2n, sizeof(int) * N);
    ESL_ALLOC(l3,   sizeof(int) * N);
  }

  L1 = L2 = L3 = 0;
  esl_vec_ISet(dd1n,   N, 0.0);
  esl_vec_ISet(l1,     N, 0.0);
  esl_vec_ISet(l2,     N, 0.0);
  if (cfg->ndomains == 2) {
    esl_vec_ISet(dd2n, N, 0.0);
    esl_vec_ISet(l3,   N, 0.0);    
  }

  /* select the length of the domains from a geometric of paramter mu
   */
  mu    = cfg->mu;
  sigma = cfg->sigma;

  d1n = d2n = 0;
  do {
    d1n   = (int)fabs(esl_rnd_Gaussian(cfg->r, mu, sigma));
    if (cfg->ndomains == 2 ) 
      d2n = (int)fabs(esl_rnd_Gaussian(cfg->r, mu, sigma));
    nit++;
  } while (nit < MAXNIT && ( d1n <= 0 || (cfg->ndomains==2 && d2n <= 0) || L - (d1n+d2n) < (int)((float)L/5.) ) );

  if (nit < MAXNIT) {
    MINHOM = (int)((float)d1n/3.);
    
    nit = 0;
    while (minhom <= MINHOM && nit < MAXNIT) {
      nit ++;
      status = short_homology_embedding(cfg, shortmsa, d1n, d2n, &L1, &L2, &L3, dd1n, dd2n, l1, l2, l3);
      if (status != eslOK) goto ERROR;
      minhom = esl_vec_IMin(dd1n, N);
      if (cfg->ndomains == 2 && esl_vec_IMin(dd2n, N) < minhom) minhom =  esl_vec_IMin(dd2n, N);
    }
  }

  if (nit == MAXNIT) printf("could not find short homology\n");
  else  {
    (*ret_nhom) ++;    
    msa_annotate_dompid(cfg, shortmsa, d1n, d2n, L1, L2, L3, &d1pid, &d2pid); /* annnotate homologies pid */
    write_summary(cfg, shortmsa, d1n, d1pid, d2n, d2pid, L1, L2, L3, dd1n, dd2n, l1, l2, l3);

    msa_annotate_rf(cfg, shortmsa, d1n, d2n, L1, L2, L3);                     /* annnotate homologies rf*/
    msa_annotate_nonhomologs(cfg, shortmsa, d1n, d2n, L1, L2, L3);            /* annotate lowercase for nonhomologous positions */
    

    //write_msa(stdout, cfg->msa, FALSE); 
    write_msa(stdout, shortmsa, FALSE); 
    write_msa(cfg->setfp, shortmsa, FALSE); 
  }
  
  esl_msa_Destroy(shortmsa); 
  free(msaname);
  free(dd1n);
  free(l1);
  free(l2);
  if (cfg->ndomains == 2) {
    free(dd2n);
    free(l3);
  }

  return eslOK;

 ERROR:

  if (shortmsa) esl_msa_Destroy(shortmsa); 
  if (msaname) free(msaname);
  if (dd1n) free(dd1n);
  if (dd2n) free(dd2n);
  if (l1) free(l1);
  if (l2) free(l2);
  if (l3)free(l3);
 return eslOK;
}


static int
short_homology_embedding(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int *ret_L1, int *ret_L2, int *ret_L3, int *dd1n, int *dd2n, int *l1, int *l2, int *l3)
{
  ESL_DSQ *dsq;
  int      N = msa->nseq;
  int      L = msa->alen;
  int      L1,L2,L3;		/* lengths of three random regions    */
  int      i, j;
  int      x;
  
  L1 = L2 = L3 = 0;
  esl_vec_ISet(dd1n,   N, 0.0);
  esl_vec_ISet(l1,     N, 0.0);
  esl_vec_ISet(l2,     N, 0.0);
  if (cfg->ndomains == 2) {
    esl_vec_ISet(dd2n, N, 0.0);
    esl_vec_ISet(l3,   N, 0.0);    
  }

  if (cfg->ndomains == 2) 
    {
      /* Select random lengths of three flanking domains;
       * Imagine picking two "insert after" points i,j in sequence 1..L', for
       * L' = L-d1n-d2n (the total length of nonhomologous test seq)
       */
      do {
	i = esl_rnd_Roll(cfg->r, L - d1n - d2n + 1 ); /* i = 0..L' */
	j = esl_rnd_Roll(cfg->r, L - d1n - d2n + 1 ); /* j = 0..L' */
      } while (i > j);
      
      /* now 1           .. i         = random region 1 (if i==0, there's none); 
       *     i+1         .. i+d1n     = domain 1
       *     i+d1n+1     .. j+d1n     = random region 2 (if i==j, there's none);
       *     j+d1n+1     .. j+d1n+d2n = domain 2
       *     j+d1n+d2n+1 .. L         = random region 3 (if j == L-d1n-d2n, there's none);
       */
      L1 = i;			
      L2 = j-i;
      L3 = L - d1n - d2n - j;
    }
  else 
    { /* embedding one domain */
      i = esl_rnd_Roll(cfg->r, L - d1n + 1 ); /* i = 0..L' */
      /* now 1           .. i         = random region 1 (if i==0, there's none); 
       *     i+1         .. i+d1n     = domain 1
       *     i+d1n+1     .. L         = random region 2 (if i==j, there's none);
       */
      L1 = i;			
      L2 = L - d1n - L1;
      L3 = 0;
    }
  
  for (x = 0; x < msa->nseq; x ++) {
    dsq = msa->ax[x];
 
    length_dealigned(cfg->abc, dsq+i+1, d1n, x, dd1n);   
    if (cfg->ndomains == 2) 
      length_dealigned(cfg->abc, dsq+j+d1n+1, d2n, x, dd2n);
    
    if (!cfg->fragone || x == 1) {
      set_random_segment(cfg, dsq+1,       L1, x, l1); if (l1[x] > L1)  esl_fatal("bad domain L1 %d l1 %d\n", L1, l1[x]);
      set_random_segment(cfg, dsq+i+d1n+1, L2, x, l2); if (l2[x] > L2)  esl_fatal("bad domain L2 %d l2 %d\n", L2, l2[x]);      
      if (cfg->ndomains == 2) {
	set_random_segment(cfg, dsq+j+d1n+d2n+1, L3, x, l3); if (l3[x] > L3)  esl_fatal("bad domain L3 %d l3 %d\n", L3, l3[x]);
      }
    }

  }
  
  *ret_L1    = L1;
  *ret_L2    = L2;
  *ret_L3    = L3;

  return eslOK;
}


static int
set_random_segment(struct cfg_s *cfg, ESL_DSQ *dsq, int L, int idx, int *l)
{
  ESL_SQ  *sq           = esl_sq_CreateDigital(cfg->abc);
  char    *pkey         = NULL;
  int      start, end;
  int64_t  Lseq;
  int      len = 0;
  int      x;
  int      status;

  if (L==0) { l[idx] = 0; return eslOK; }

  if (L > cfg->db_maxL) esl_fatal("can't fetch a segment of length %d; database max is %d\n", L, cfg->db_maxL);

  /* fetch a random subseq from the source database */
  esl_sq_GrowTo(sq, L);
  do {                                                     
    if (pkey != NULL) free(pkey);
    if (esl_ssi_FindNumber(cfg->dbfp->data.ascii.ssi, esl_rnd_Roll(cfg->r, cfg->db_nseq), NULL, NULL, NULL, &Lseq, &pkey) != eslOK)
      esl_fatal("failed to look up a random seq");
  } while (Lseq < L);
  
  start = 1 + esl_rnd_Roll(cfg->r, Lseq-L);              
  end   = start + L - 1;
  if (esl_sqio_FetchSubseq(cfg->dbfp, pkey, start, end, sq) != eslOK) esl_fatal("failed to fetch subseq");
  esl_sq_ConvertDegen2X(sq);
  
  /* Now apply the appropriate randomization algorithm */
  status = esl_rsq_XShuffle  (cfg->r, sq->dsq, L, sq->dsq);
  if (status != eslOK) esl_fatal("esl's shuffling failed");

  for (x = 0; x < L; x ++) 
    if (! esl_abc_XIsGap(cfg->abc, dsq[x]) && ! esl_abc_XIsMissing(cfg->abc, dsq[x]) ) {
      memcpy(dsq+x, sq->dsq+1+len, sizeof(ESL_DSQ));
      len ++;
    }
  
  l[idx] = len;

  esl_sq_Destroy(sq);
  free(pkey);
  return eslOK;
}
 
static int
length_dealigned(ESL_ALPHABET *abc, ESL_DSQ *dsq, int alen, int idx, int *l)
{
  int len = 0;
  int x;

  for (x = 0; x < alen; x ++) 
    if (! esl_abc_XIsGap(abc, dsq[x]) && ! esl_abc_XIsMissing(abc, dsq[x]) ) {
      len ++;
    }

  l[idx] = len;

  return eslOK;
} 

static int
write_msa(FILE *fp, ESL_MSA *msa, int verbose)
{
  MSA_STAT *msastat = NULL;

  if (esl_msafile_Write(fp, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write train msa to file"); 
  if (verbose) {
    msamanip_XStats(msa, &msastat); //  msa aveid and avematch 
    msamanip_DumpStats(stdout, msa, msastat); 
  }
  
  if (msastat) free(msastat);
  return eslOK;
}

static int
msa_annotate_dompid(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int L1, int L2, int L3, double *ret_d1pid, double *ret_d2pid)
{
  char             tmpmsafile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  ESL_MSAFILE   *afp  = NULL;
  FILE            *fp  = NULL;
  char            *args = NULL;
  char            *s = NULL;
  ESL_MSA         *msadom = NULL;
  double           d1pid = -1.0;
  double           d2pid = -1.0;
  int              i, j;
  int               L;
  int              status;

  /* MSA input in stockholm format */
  if ((status = esl_tmpfile_named(tmpmsafile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create msa file");
  if ((status = esl_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_STOCKHOLM)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "Failed to write STOCKHOLM file\n");
  fclose(fp);

  /*     1           .. i         = random region 1 (if i==0, there's none); 
   *     i+1         .. i+d1n     = domain 1
   *     i+d1n+1     .. j+d1n     = random region 2 (if i==j, there's none);
   *     j+d1n+1     .. j+d1n+d2n = domain 2
   *     j+d1n+d2n+1 .. L         = random region 3 (if j == L-d1n-d2n, there's none);
   */
  /*     1           .. i         = random region 1 (if i==0, there's none); 
   *     i+1         .. i+d1n     = domain 1
   *     i+d1n+1     .. L         = random region 2 (if i==j, there's none);
   */
  i = L1;
  j = i + L2;
  L = L1 + L2 + L3 + d1n + d2n;
  
 
  /* the truncated msa file  */
  if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
  fclose(fp);  
  /* run alimask */
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-alimask -t %s %d-%d > %s", s, tmpmsafile, i+1, i+d1n, tmpoutfile);  
  system(args);
  /* the truncated msa  */
  status = esl_msafile_Open(&(cfg->abc), tmpoutfile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  status = esl_msafile_Read(afp, &msadom);
  esl_msafile_Close(afp);
  d1pid = 100. * avg_pid_msa(msadom);

  if (cfg->ndomains == 2) {
    esl_msa_Destroy(msadom); msadom = NULL;
    remove(tmpoutfile);
    
    /* the truncated msa file  */
    if ((status = esl_tmpfile_named(tmpoutfile,  &fp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create out file");
    fclose(fp);  
    /* run alimask */
    if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
    if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
    esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-alimask -t %s %d-%d > %s", s, tmpmsafile, j+d1n+1, j+d1n+d2n, tmpoutfile);  
    system(args);
    /* the truncated msa  */
    status = esl_msafile_Open(&(cfg->abc), tmpoutfile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
    if (status != eslOK) esl_msafile_OpenFailure(afp, status);
    status = esl_msafile_Read(afp, &msadom);
    esl_msafile_Close(afp);
    d2pid = 100. * avg_pid_msa(msadom);
  }
   
  *ret_d1pid = d1pid;
  *ret_d2pid = d2pid;

  esl_msa_Destroy(msadom);
  remove(tmpmsafile);
  remove(tmpoutfile);
  return eslOK;

 ERROR:
  if (msadom) esl_msa_Destroy(msadom);
  remove(tmpmsafile);
  remove(tmpoutfile);
  return status;
}

static int
msa_annotate_rf(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int L1, int L2, int L3)
{
  int apos;
  int i, j;
  int L;
  int status;

  ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
  for (apos = 0; apos < msa->alen; apos++) 
    msa->rf[apos] = '.';
  msa->rf[msa->alen] = '\0';

  /*     1           .. i         = random region 1 (if i==0, there's none); 
   *     i+1         .. i+d1n     = domain 1
   *     i+d1n+1     .. j+d1n     = random region 2 (if i==j, there's none);
   *     j+d1n+1     .. j+d1n+d2n = domain 2
   *     j+d1n+d2n+1 .. L         = random region 3 (if j == L-d1n-d2n, there's none);
   */
  /*     1           .. i         = random region 1 (if i==0, there's none); 
   *     i+1         .. i+d1n     = domain 1
   *     i+d1n+1     .. L         = random region 2 (if i==j, there's none);
   */
  i = L1;
  j = i + L2;
  L = L1 + L2 + L3 + d1n + d2n;
  
  for (apos = 0; apos < msa->alen; apos++) { 
    if (apos >= i     && apos < i+d1n)     msa->rf[apos] = 'x'; /* domain 1: watch off by one: rf[0..alen-1]; */
    if (apos >= j+d1n && apos < j+d1n+d2n) msa->rf[apos] = 'x'; /* domain 2: watch off by one: rf[0..alen-1]; */
  }
  
  return eslOK;

 ERROR:
  return status;
}
  
static int
msa_annotate_nonhomologs(struct cfg_s *cfg, ESL_MSA *msa, int d1n, int d2n, int L1, int L2, int L3)
{
  int idx;
  int apos;
  int i, j;
  int L;

  esl_msa_Textize(msa);
  
  /*     1           .. i         = random region 1 (if i==0, there's none); 
   *     i+1         .. i+d1n     = domain 1
   *     i+d1n+1     .. j+d1n     = random region 2 (if i==j, there's none);
   *     j+d1n+1     .. j+d1n+d2n = domain 2
   *     j+d1n+d2n+1 .. L         = random region 3 (if j == L-d1n-d2n, there's none);
   */
  /*     1           .. i         = random region 1 (if i==0, there's none); 
   *     i+1         .. i+d1n     = domain 1
   *     i+d1n+1     .. L         = random region 2 (if i==j, there's none);
   */
  i = L1;
  j = i + L2;
  L = L1 + L2 + L3 + d1n + d2n;
  
  for (idx = 0; idx < msa->nseq; idx++)
    {
      for (apos = 0; apos < msa->alen; apos++) {
	if (                     apos < i)     msa->aseq[idx][apos] = tolower(msa->aseq[idx][apos]);    /* ramdom region 1 */
	if (apos >= i+d1n     && apos < j+d1n) msa->aseq[idx][apos] = tolower(msa->aseq[idx][apos]);    /* ramdom region 2 */
 	if (apos >= j+d1n+d2n && apos < L)     msa->aseq[idx][apos] = tolower(msa->aseq[idx][apos]);    /* ramdom region 3 */
      }
    } 

  return eslOK;
}
  

static int
write_summary(struct cfg_s *cfg, ESL_MSA *msa, int d1n, double d1pid, int d2n, double d2pid, int L1, int L2, int L3, int *dd1n, int *dd2n, int *l1, int *l2, int *l3)
{
  double avgpid = 100.*avg_pid_msa(msa);
  int    L      = cfg->msa->alen;
  int    x;
  int    l;

  if (cfg->ndomains == 2) 
    {
      fprintf(cfg->summfp, "%-35s %5d %.2f %5d %5d %.2f %5d %5d %.2f %5d", msa->name, L, avgpid, L1, d1n, d1pid, L2, d2n, d2pid, L3);
      printf(              "%-35s %5d %.2f %5d %5d %.2f %5d %5d %.2f %5d", msa->name, L, avgpid, L1, d1n, d1pid, L2, d2n, d2pid, L3);
    }
  else 
    { 
      fprintf(cfg->summfp, "%-35s %5d %.2f %5d %5d %.2f %5d", msa->name, L, avgpid, L1, d1n, d1pid, L2);
      printf(              "%-35s %5d %.2f %5d %5d %.2f %5d", msa->name, L, avgpid, L1, d1n, d1pid, L2);
    }
  
  for (x = 0; x < msa->nseq; x ++) {    
    if (cfg->ndomains == 2) {
      l = l1[x] + l2[x] + dd1n[x] + l3[x] + dd2n[x];
      printf(              " %-35s  %5d %5d %5d %5d %5d %5d", msa->sqname[x], l, l1[x], dd1n[x], l2[x], dd2n[x], l3[x]);
      fprintf(cfg->summfp, " %-35s  %5d %5d %5d %5d %5d %5d", msa->sqname[x], l, l1[x], dd1n[x], l2[x], dd2n[x], l3[x]);
    }
    else {
      l = l1[x] + l2[x] + dd1n[x];      
      printf(              " %-35s  %5d %5d %5d %5d", msa->sqname[x], l, l1[x], dd1n[x], l2[x]);
      fprintf(cfg->summfp, " %-35s  %5d %5d %5d %5d", msa->sqname[x], l, l1[x], dd1n[x], l2[x]);
    }
  }
  fprintf(cfg->summfp, "\n");
  fprintf(stdout,      "\n");

  return eslOK;
}

static double
avg_pid_msa(ESL_MSA *msa)
{
  double avg_pid = 0.0;
  double pid;
  int    n = 0;
  int    i, j;
  
  if (msa->nseq == 1) return -1.0;

  for (i = 0; i < msa->nseq; i++)
    for (j = i+1; j < msa->nseq; j++) /* traverse test seq stack without destroying/popping */
      {
	esl_dst_XPairId(msa->abc, msa->ax[i], msa->ax[j], &pid, NULL, NULL);	
	avg_pid += pid;
	n ++;
      }  
  if (n > 0) avg_pid /= (double)n;
  
  return avg_pid;
}
