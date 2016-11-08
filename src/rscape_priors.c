/* rscape_priors.c
 *
 * use easel code to calculate dirichlet priors for the
 * R-scape joint probabilities pij (from the ssu-lsu cwr alignments)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "esl_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dirichlet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msaweight.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "msamanip.h"
#include "msatree.h"

#define IDX(i,j,L)  ( (i) * (L) + (j) )

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
  
  int              nmsa;
  char            *msafile;
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */
  char            *msaheader;          /* header for all msa-specific output files */
  int              infmt;
  char            *name;
  
  int              submsa;              /* set to the number of random seqs taken from original msa.
					 * Set to 0 if we are taking all */
  MSA_STAT        *mstat;               /* statistics of the input alignment */
  float           *msafrq;

  float            thresh1;            /* cluster threshold */
  float            thresh2;            /* percentage id threshold */

  int              N;                  /* number of mixture diritchlet */
  float            tol;
  int              verbose;
};

static ESL_OPTIONS options[] = {
  /* name           type              default   env  range   toggles        reqs   incomp               help                                       docgroup*/
  { "-h",           eslARG_NONE,        FALSE,   NULL, NULL, NULL,   NULL,    NULL,               "show brief help on version and usage",       
                            0 },
  { "-v",           eslARG_NONE,        FALSE,   NULL, NULL, NULL,   NULL,    NULL,               "be verbose",                                                                                       0 },  
  { "-N",            eslARG_INT,          "1",   NULL,"n>0", NULL,   NULL,    NULL,               "number of mixture dirichlets",                                                                     0 },
  /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "-I",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < x id",                                                               0 },
  { "--submsa",       eslARG_INT,       FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   0 },
  { "--minid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      0 },
  { "--maxid",        eslARG_REAL,      FALSE,   NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",  
                             0},
  /* options for ribosum thresholds */
  { "-1",             eslARG_REAL,    "0.80",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "clustering threshold",                                                                      0 },
  { "-2",             eslARG_REAL,    "0.65",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "percentage id threshold",                                                                   0 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE,      NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "direct output to file <f>, not stdout",                                                     0 },
  /* msa format */
  { "--informat",    eslARG_STRING,      NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            0 },
  /* other options */
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,       "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                                       5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "given an msa annotated with RNA secondary structure, calculate dirichlet priors for the joints p_ij";

static int fit_dirichlet(FILE *fp, int N, int K, double **c, int verbose);
static int add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, double *xrnaJ, double *prnaJ, double *urnaJ, double *xrnaM, double *prnaM, double *urnaM, char *errbuf);
static int xrnaJ_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa,           double *xrnaJ, char *errbuf);
static int prnaJ_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *prnaJ, char *errbuf);
static int urnaJ_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *urnaJ, char *errbuf);
static int xrnaM_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa,           double *xrnaJ, char *errbuf);
static int prnaM_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *prnaJ, char *errbuf);
static int urnaM_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *urnaJ, char *errbuf);

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, struct cfg_s *ret_cfg)
{
  ESL_GETOPTS  *go = esl_getopts_Create(options);
  struct cfg_s  cfg;
  int           status;
  
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\noptions:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/
      exit(0);
    }
  
  /* the msa with RNA secondary structure */
  cfg.msafile = NULL;
  if (esl_opt_ArgNumber(go) != 1) { if (puts("Incorrect number of command line arguments.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.msafile  = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
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
  cfg.abc = esl_alphabet_Create(eslRNA);
  
  cfg.w = esl_stopwatch_Create(); 
  
  /*  output file */
  if ( esl_opt_IsOn(go, "-o") ) {
    if ((cfg.outfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else cfg.outfp = stdout;
  
  /* options for mdifying the input msa */
  cfg.fragfrac   = esl_opt_IsOn(go, "-F")? esl_opt_GetReal(go, "-F") : -1.0;
  cfg.idthresh   = esl_opt_IsOn(go, "-I")? esl_opt_GetReal(go, "-I") : -1.0;
  if (esl_opt_IsOn(go, "--submsa")) cfg.submsa = esl_opt_GetInteger(go, "--submsa");
  else                              cfg.submsa = 0;
  
  /* options for mdifying the input msa */
  cfg.thresh1 = esl_opt_GetReal(go, "-1");
  cfg.thresh2 = esl_opt_GetReal(go, "-2");

  cfg.N = esl_opt_GetInteger(go, "-N");

  /* other options */
  cfg.tol     = esl_opt_GetReal   (go, "--tol");
  cfg.verbose = esl_opt_GetBoolean(go, "-v");
  
  esl_sprintf(&cfg.name, "%s%d-%d", "DIRICHLET", (int)(cfg.thresh1*100), (int)(cfg.thresh2*100));
  
  *ret_go  = go;
  *ret_cfg = cfg;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere options are:")                                       < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 0, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  
  
 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
output_header(FILE *ofp, ESL_GETOPTS *go)
{
  esl_banner(ofp, go->argv[0], banner);
  
  if (esl_opt_IsUsed(go, "-o")          && fprintf(ofp, "# output directed to file:         %s\n",             esl_opt_GetString(go, "-o"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}


int
main(int argc, char **argv)
{
  char                *msg = "rscape-priors failed";
  ESL_GETOPTS         *go;
  struct cfg_s         cfg;
  ESL_MSAFILE        *afp = NULL;
  ESL_MSA             *msa = NULL;                /* the input alignment  */
  double             **xrnaJ = NULL;              // all      positions joint     16 array:     xrnaJ[0..nc-1][IDX(a,b)] = xrnaJ[0..nc-1][IDX(b,a)]
  double             **prnaJ = NULL;              // paired   positions joint     16 array:     prnaJ[0..nc-1][IDX(a,b)] = prnaJ[0..nc-1][IDX(b,a)]
  double             **urnaJ = NULL;              // unpaired positions joint     16 array:     urnaJ[0..nc-1][IDX(a,b)] = urnaJ[0..nc-1][IDX(b,a)]
  double             **xrnaM = NULL;              // all      positions marginals  4 array:     xrnaM[0..nc-1][a] 
  double             **prnaM = NULL;              // paired   positions marginals  4 array:     prnaM[0..nc-1][a]
  double             **urnaM = NULL;              // unpaired positions marginals  4 array:     urnaM[0..nc-1][a]
  int                  K, K2;
  int                  nalloc = 10;
  int                  seq_cons_len = 0;
  int                  nfrags = 0;	  	  /* # of fragments removed */
  int                  nremoved = 0;	          /* # of identical sequences removed */
  int                  n;
  int                  idx;
  int                  i;
  int                  status = eslOK;
  int                  hstatus = eslOK;

  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    
  
  K  = cfg.abc->K;
  K2 = K * K;
  ESL_ALLOC(xrnaJ, sizeof(double *) * nalloc);
  ESL_ALLOC(prnaJ, sizeof(double *) * nalloc);
  ESL_ALLOC(urnaJ, sizeof(double *) * nalloc);
  ESL_ALLOC(xrnaM, sizeof(double *) * nalloc);
  ESL_ALLOC(prnaM, sizeof(double *) * nalloc);
  ESL_ALLOC(urnaM, sizeof(double *) * nalloc);
  ESL_ALLOC(xrnaJ[0], sizeof(double) * K2 * nalloc);
  ESL_ALLOC(prnaJ[0], sizeof(double) * K2 * nalloc); 
  ESL_ALLOC(urnaJ[0], sizeof(double) * K2 * nalloc); 
  ESL_ALLOC(xrnaM[0], sizeof(double) *  K * nalloc); 
  ESL_ALLOC(prnaM[0], sizeof(double) *  K * nalloc); 
  ESL_ALLOC(urnaM[0], sizeof(double) *  K * nalloc); 

  for (n = 1; n < nalloc; n ++) xrnaJ[n] = xrnaJ[0] + n * K2;
  for (n = 1; n < nalloc; n ++) prnaJ[n] = prnaJ[0] + n * K2;
  for (n = 1; n < nalloc; n ++) urnaJ[n] = urnaJ[0] + n * K2;
  for (n = 1; n < nalloc; n ++) xrnaM[n] = xrnaM[0] + n *  K;
  for (n = 1; n < nalloc; n ++) prnaM[n] = prnaM[0] + n *  K;
  for (n = 1; n < nalloc; n ++) urnaM[n] = urnaM[0] + n *  K;

  /* Open the MSA file */
  status = esl_msafile_Open(NULL, cfg.msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  
  /* read the MSA */
  while ((hstatus = esl_msafile_Read(afp, &msa)) != eslEOF) {
    if (hstatus != eslOK) esl_msafile_ReadFailure(afp, status);
    cfg.nmsa ++;
    if (cfg.nmsa == nalloc) {
      nalloc += 10;
      ESL_REALLOC(xrnaJ,    sizeof(double *) * nalloc);
      ESL_REALLOC(prnaJ,    sizeof(double *) * nalloc);
      ESL_REALLOC(urnaJ,    sizeof(double *) * nalloc);
      ESL_REALLOC(xrnaM,    sizeof(double *) * nalloc);
      ESL_REALLOC(prnaM,    sizeof(double *) * nalloc);
      ESL_REALLOC(urnaM,    sizeof(double *) * nalloc);
      ESL_REALLOC(xrnaJ[0], sizeof(double *) * nalloc * K2);
      ESL_REALLOC(prnaJ[0], sizeof(double *) * nalloc * K2);
      ESL_REALLOC(urnaJ[0], sizeof(double *) * nalloc * K2);
      ESL_REALLOC(xrnaM[0], sizeof(double *) * nalloc *  K);
      ESL_REALLOC(prnaM[0], sizeof(double *) * nalloc *  K);
      ESL_REALLOC(urnaM[0], sizeof(double *) * nalloc *  K);
      for (n = cfg.nmsa; n < nalloc; n ++) xrnaJ[n] = xrnaJ[0] + n * K2;
      for (n = cfg.nmsa; n < nalloc; n ++) prnaJ[n] = prnaJ[0] + n * K2;
      for (n = cfg.nmsa; n < nalloc; n ++) urnaJ[n] = urnaJ[0] + n * K2;
      for (n = cfg.nmsa; n < nalloc; n ++) xrnaM[n] = xrnaM[0] + n *  K;
      for (n = cfg.nmsa; n < nalloc; n ++) prnaM[n] = prnaM[0] + n *  K;
      for (n = cfg.nmsa; n < nalloc; n ++) urnaM[n] = urnaM[0] + n *  K;
   }
   
    /* select submsa */
    if (cfg.submsa) {
      if (msamanip_SelectSubset(cfg.r, cfg.submsa, &msa, NULL, cfg.errbuf, cfg.verbose) != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
    }
    
    /* outheader for all msa-output files */
    msamanip_OutfileHeader((msa->acc)?msa->acc:cfg.msafile, &cfg.msaheader); 
    
    esl_msa_Hash(msa);
    
    if (esl_opt_IsOn(go, "-F") && msamanip_RemoveFragments(cfg.fragfrac, &msa, &nfrags, &seq_cons_len) != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
    if (esl_opt_IsOn(go, "-I"))   msamanip_SelectSubsetBymaxID(cfg.r, &msa, cfg.idthresh, TRUE, &nremoved);
    
    /* given msa aveid and avematch */
    msamanip_CStats(cfg.abc, msa, &cfg.mstat);
    cfg.msafrq = NULL;
    //msamanip_CBaseComp(cfg.abc, msa, cfg.bg->f, &cfg.msafrq);
    
    if (esl_opt_IsOn(go, "--minid") && cfg.mstat->avgid < 100.*esl_opt_GetReal(go, "--minid")) continue;
    if (esl_opt_IsOn(go, "--maxid") && cfg.mstat->avgid > 100.*esl_opt_GetReal(go, "--maxid")) continue;
    
    /* Make all seqs upper case */
    for (idx = 0; idx < msa->nseq; idx++) 
      for (i = 0; i < msa->alen; i ++)
      msa->aseq[idx][i] = toupper(msa->aseq[idx][i]);
    
    /* Print some stuff about what we're about to do.
     */
    if (msa->name != NULL) printf("Alignment:           %s\n",  msa->name);
    else                   printf("Alignment:           #%d\n", cfg.nmsa);
    printf                       ("Number of sequences: %d\n",  msa->nseq);
    printf                       ("Number of columns:   %"PRId64"\n",  msa->alen);
    puts("");
    fflush(stdout);
 
    status = esl_msaweight_GSC(msa, NULL);
    if (status != eslOK) ESL_XFAIL(eslFAIL, cfg.errbuf, "failed msaweight_GSC");
 
    add_counts(cfg.abc, msa, xrnaJ[cfg.nmsa-1], prnaJ[cfg.nmsa-1], urnaJ[cfg.nmsa-1], xrnaM[cfg.nmsa-1], prnaM[cfg.nmsa-1], urnaM[cfg.nmsa-1], cfg.errbuf);

    esl_msa_Destroy(msa); msa = NULL;
    if (cfg.msafrq) free(cfg.msafrq); cfg.msafrq = NULL;
    free(cfg.msaheader); cfg.msaheader = NULL;
  }

  status = fit_dirichlet(cfg.outfp, cfg.nmsa, K2, xrnaJ, cfg.verbose);
  status = fit_dirichlet(cfg.outfp, cfg.nmsa, K2, prnaJ, cfg.verbose);
  status = fit_dirichlet(cfg.outfp, cfg.nmsa, K2, urnaJ, cfg.verbose);
  status = fit_dirichlet(cfg.outfp, cfg.nmsa,  K, xrnaM, cfg.verbose);
  status = fit_dirichlet(cfg.outfp, cfg.nmsa,  K, prnaM, cfg.verbose);
  status = fit_dirichlet(cfg.outfp, cfg.nmsa,  K, urnaM, cfg.verbose);
   
  if (status != eslOK) esl_fatal(msg);
 
  /* cleanup */
  esl_stopwatch_Destroy(cfg.w);
  fclose(cfg.outfp);
  free(cfg.outheader);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  esl_msafile_Close(afp);
  free(xrnaJ[0]);
  free(prnaJ[0]);
  free(urnaJ[0]);
  free(xrnaM[0]);
  free(prnaM[0]);
  free(urnaM[0]);
  free(xrnaJ);
  free(prnaJ);
  free(urnaJ);
  free(xrnaM);
  free(prnaM);
  free(urnaM);
  return status;

 ERROR:
  exit(1);
}


static int
fit_dirichlet(FILE *fp, int nc, int K, double **c, int verbose)
{
  ESL_MIXDCHLET       *d = NULL;
  double               guess;
  int                  status;
 
  guess = esl_vec_DArgMax(c[0], K); 
  guess = 1.;

  d = esl_mixdchlet_Create(1, K);
  if (status != eslOK) goto ERROR;
  esl_vec_DSet(d->pq,       1, 1.0);
  esl_vec_DSet(d->alpha[0], K, guess);

  status = esl_mixdchlet_Fit(c, nc, d, TRUE);
  if (status != eslOK) goto ERROR;

  /* write the dirichlets */
  esl_mixdchlet_Write(fp, d);
 
  esl_mixdchlet_Destroy(d);

  return eslOK;
  
 ERROR:
  return status;
}
 

static int                  
add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, double *xrnaJ, double *prnaJ, double *urnaJ, double *xrnaM, double *prnaM, double *urnaM, char *errbuf)
{
  int    **ct = NULL;   // ct array for each seq
  int      i;
  int      status;
  
  ESL_ALLOC(ct, sizeof(int *) * msa->nseq);
  ct[0] = NULL;

  ESL_ALLOC(ct[0], sizeof(int) * (msa->alen+1) * msa->nseq);
  for (i = 1; i < msa->nseq; i ++) ct[i] = ct[0] + i * (msa->alen+1);

  for (i = 0; i < msa->nseq; i ++) {
    if      (msa->ss_cons) esl_wuss2ct(msa->ss_cons, msa->alen, ct[i]);
    else if (msa->ss[i])   esl_wuss2ct(msa->ss[i],   msa->alen, ct[i]);
    else                   ESL_XFAIL(eslFAIL, "no ss for sequence %d\n", errbuf);
  }  

  xrnaJ_add_counts(abc, msa,     xrnaJ, errbuf);
  prnaJ_add_counts(abc, msa, ct, prnaJ, errbuf);
  urnaJ_add_counts(abc, msa, ct, urnaJ, errbuf);
  xrnaM_add_counts(abc, msa,     xrnaM, errbuf);
  prnaM_add_counts(abc, msa, ct, prnaM, errbuf);
  urnaM_add_counts(abc, msa, ct, urnaM, errbuf);
 
  free(ct[0]);
  free(ct);
  
  return eslOK;

 ERROR:
  if (ct[0]) free(ct[0]);
  if (ct) free(ct);
  
  return status;
}

static int                  
xrnaJ_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, double *xrnaJ, char *errbuf)
{
  double wgt;
  int    K2 = abc->K * abc->K;
  int    i;
  int    posa, posb;
  int    c, cc;
  int    pair;

  esl_vec_DSet(xrnaJ, K2, 0.0);

  for (i = 0; i < msa->nseq; i++) {
    wgt = msa->wgt[i]; 
    
    for (posa = 0; posa < msa->alen; posa++) {
      c = msa->aseq[i][posa];
      
      if (esl_abc_CIsCanonical(abc, c)) 
	{
	  for (posb = 0; posb < msa->alen; posb++) {
	    if (posb == posa) continue;
	    
	    cc = msa->aseq[i][posb];
	    if (esl_abc_CIsCanonical(abc, cc)) 
	      {
		pair = IDX(esl_abc_DigitizeSymbol(abc,  c), esl_abc_DigitizeSymbol(abc, cc), abc->K);
		xrnaJ[pair] += wgt;
	      }
	  }
	}
    }
  }
  return eslOK;
}

static int                  
prnaJ_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *prnaJ, char *errbuf)
{
  double wgt;
  int    K2 = abc->K * abc->K;
  int    i;
  int    c, cc;
  int    posa;
  int    cposa;
  int    pair;

  esl_vec_DSet(prnaJ, K2, 0.0);

  for (i = 0; i < msa->nseq; i++) {
      wgt = msa->wgt[i]; 
      
      for (posa = 0; posa < msa->alen; posa++) {
	c = msa->aseq[i][posa];
	
	if (esl_abc_CIsCanonical(abc, c)) 
	  {
	    cposa = ct[i][posa+1] - 1;
	    if (cposa == -1) continue;

	    cc = msa->aseq[i][cposa];	    	    
	    if (esl_abc_CIsCanonical(abc, cc)) 
	      {
		pair = IDX(esl_abc_DigitizeSymbol(abc, c), esl_abc_DigitizeSymbol(abc, cc), abc->K);
		prnaJ[pair] += wgt;
	      } 
	  }
      }
  }

  return eslOK;
}

static int                  
urnaJ_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *urnaJ, char *errbuf)
{
  double wgt;
  int    K2 = abc->K * abc->K;
  int    i;
  int    c, cc;
  int    posa, posb;
  int    cposa;
  int    pair;

  esl_vec_DSet(urnaJ, K2, 0.0);

  for (i = 0; i < msa->nseq; i++) {
    wgt = msa->wgt[i]; 
    
    for (posa = 0; posa < msa->alen; posa++) {
      c = msa->aseq[i][posa];
      
      if (esl_abc_CIsCanonical(abc, c)) 
	{
	  cposa = ct[i][posa+1] - 1;
	  if (cposa >= 0) continue;
	  
	  for (posb = 0; posb < msa->alen; posb++) {
	    if (posb == posa) continue;
	    
	    cc = msa->aseq[i][posb];	    	    
	    if (esl_abc_CIsCanonical(abc, cc)) 
	      {
		pair = IDX(esl_abc_DigitizeSymbol(abc, c), esl_abc_DigitizeSymbol(abc, cc), abc->K);
		urnaJ[pair] += wgt;
	      } 
	  }
	}
    }
  }
  
   return eslOK;
}
static int                  
xrnaM_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, double *xrnaM, char *errbuf)
{
  double wgt;
  int    K = abc->K;
  int    i;
  int    pos;
  int    c;

  esl_vec_DSet(xrnaM, K, 0.0);

  for (i = 0; i < msa->nseq; i++) {
    wgt = msa->wgt[i]; 
    
    for (pos = 0; pos < msa->alen; pos++) {
      c = msa->aseq[i][pos];
      
      if (esl_abc_CIsCanonical(abc, c)) xrnaM[esl_abc_DigitizeSymbol(abc,c)] += wgt;
	      
    }
  }
  return eslOK;
}

static int                  
prnaM_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *prnaM, char *errbuf)
{
  double wgt;
  int    K = abc->K;
  int    i;
  int    c, cc;
  int    pos;
  int    cpos;

  esl_vec_DSet(prnaM, 2, 0.0);

  for (i = 0; i < msa->nseq; i++) {
      wgt = msa->wgt[i]; 
      
      for (pos = 0; pos < msa->alen; pos++) {
	c = msa->aseq[i][pos];
	
	if (esl_abc_CIsCanonical(abc, c)) 
	  {
	    cpos = ct[i][pos+1] - 1;
	    if (cpos == -1) continue;

	    cc = msa->aseq[i][cpos];	    	    
	    if (esl_abc_CIsCanonical(abc, cc)) prnaM[esl_abc_DigitizeSymbol(abc, c)] += wgt;
	  }
      }
  } 
  return eslOK;
}

static int                  
urnaM_add_counts(ESL_ALPHABET *abc, ESL_MSA *msa, int **ct, double *urnaM, char *errbuf)
{
  double wgt;
  int    K = abc->K;
  int    i;
  int    c, cc;
  int    pos;
  int    cpos;

  esl_vec_DSet(urnaM, K, 0.0);

  for (i = 0; i < msa->nseq; i++) {
    wgt = msa->wgt[i]; 
    
    for (pos = 0; pos < msa->alen; pos++) {
      c = msa->aseq[i][pos];
      
      if (esl_abc_CIsCanonical(abc, c)) 
	{
	  cpos = ct[i][pos+1] - 1;
	  if (cpos >= 0) continue;
	  
	  cc = msa->aseq[i][cpos];	    	    
	  if (esl_abc_CIsCanonical(abc, cc)) urnaM[esl_abc_DigitizeSymbol(abc,c)] += wgt;
	}
    }
  }

 esl_vec_DDump(stdout, urnaM, K, NULL);

  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/

 
