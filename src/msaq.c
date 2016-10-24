/* msaq
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
#include "esl_stopwatch.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "rscape_config.h"

#include "msamanip.h"
#include "msatree.h"

#define ALPHOPTS     "--amino,--dna,--rna"                      /* Exclusive options for alphabet choice */
#define MODOPTS      "--opt,--hmm,--cm"                         /* Exclusive options for model */

typedef enum{
  OPT = 0,
  HMM = 1,
  CM  = 2,
} MODEL;

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

  char            *gnuplot;
  
  double           fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double           idthresh;	       /* fractional identity threshold for selecting subset of sequences */
  double           gapthresh;          /* only keep columns with <= gapthresh fraction of gaps in them */
  int              singlelink;         /* if TRUE use single linkage clustering to remove highly identical sequences */
  double           maxnowc;            /* max fraction of non wc basepairs allowd per column */   
  
  double           minidthresh;	       /* fractional minimal identity threshold for selecting subset of sequences */
    
  int              onemsa;
  int              nmsa;
  char            *msafile;
  char            *filename;
  char            *msaname;
  int             *msamap;
  int             *msarevmap;
  int              firstpos;
  int              maxsq_gsc;        /* maximum number of seq to switch from GSC to PG sequence weighting */
  int              tstart;
  int              tend;

  int              onbpairs;
  int              nbpairs;

  MODEL            model;              /* cm or hmm */
  
  char            *outdir;
  char            *outheader;          /* header for all output files */
  int              infmt;
 
  char            *msaqfile;
  FILE            *msaqfp;
  char            *seqfile;
  
  double          *ft;
  double          *fbp;
  double          *fnbp;

  char            *treefile;
  FILE            *treefp;
  ESL_TREE        *T;
  double           treeavgt;
 
  int              consensus;           /* if TRUE, analyze only consensus positions in the alignment */

  int              nseqmin;
  int              submsa;              /* set to the number of random seqs taken from original msa.
					 * Set to 0 if we are taking all */
  MSA_STAT        *omstat;              /* statistics of the original alignment */
  MSA_STAT        *mstat;               /* statistics of the analyzed alignment */
  float           *msafrq;

  char            *plotfile;
  int              makeplot;
  
  float            tol;
  int              verbose;
};


static ESL_OPTIONS options[] = {
  /* name             type             default   env         range   toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,     FALSE,    NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      1 },
  { "--outdir",     eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify a directory for all output files",                                                  1 },
  { "-v",             eslARG_NONE,     FALSE,    NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                1 },
  { "--onemsa",       eslARG_NONE,     FALSE,    NULL,       NULL,   NULL,    NULL,  NULL,               "if file has more than one msa, analyze only the first one",                                 1 },
  /* options */
  { "--maxnowc",      eslARG_REAL,      NULL,    NULL,  "0<=x<=1",   NULL,    NULL,  NULL,               "max fraction of non-WC allowed",                                                            1 },
  { "--opt",          eslARG_NONE,    "TRUE",    NULL,       NULL,   MODOPTS, NULL,  NULL,               "build a CM if structure is given, otherwise an HMM  ",                                      1 },
  { "--hmm",          eslARG_NONE,     FALSE,    NULL,       NULL,   MODOPTS, NULL,  NULL,               "build an HMM ",                                                                             1 },
  { "--cm",           eslARG_NONE,     FALSE,    NULL,       NULL,   MODOPTS, NULL,  NULL,               "build an infernal cm ",                                                                             1 },
 /* options for input msa (if seqs are given as a reference msa) */
  { "-F",             eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      1 },
  { "-I",             eslARG_REAL,     "1.0",    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "require seqs to have < <x> id",                                                             1 },
  { "-i",             eslARG_REAL,      NULL,    NULL, "0<=x<1.0",   NULL,    NULL,  NULL,               "require seqs to have >= <x> id",                                                            1 },
  { "--tstart",        eslARG_INT,      FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "min alignment position to analyze [1..alen]",                                               1 },
  { "--tend",          eslARG_INT,      FALSE,   NULL,      "n>0",   NULL,    NULL,  NULL,               "max alignment position to analyze [1..alen]",                                               1 },
  { "--consensus",    eslARG_NONE,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "analyze only consensus (seq_cons) positions",                                               1 },
  { "--submsa",       eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "take n random sequences from the alignment, all if NULL",                                   1 },
  { "--nseqmin",      eslARG_INT,       NULL,    NULL,      "n>0",   NULL,    NULL,  NULL,               "minimum number of sequences in the alignment",                                              1 },  
  { "--gapthresh",    eslARG_REAL,     "1.0",    NULL,  "0<=x<=1",   NULL,    NULL,  NULL,               "keep columns with < <x> fraction of gaps",                                                  1 },
  { "--minid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "minimum avgid of the given alignment",                                                      1 },
  { "--maxid",        eslARG_REAL,      NULL,    NULL, "0<x<=1.0",   NULL,    NULL,  NULL,               "maximum avgid of the given alignment",                                                      1 },
  { "--singlelink",   eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "true to use single linkage clustering (default is esl_msaweight_IDFilter)",                 0},
  /* alphabet type */
  { "--dna",          eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "use DNA alphabet",                                                                          0 },
  { "--rna",          eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "use RNA alphabet",                                                                          0 },
  { "--amino",        eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "use protein alphabet",                                                                      0 },
  /* msa format */
  { "--informat",   eslARG_STRING,      NULL,    NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            1 },
 /* other options */  
  { "--tol",          eslARG_REAL,    "1e-3",    NULL,       NULL,   NULL,    NULL,  NULL,               "tolerance",                                                                                 0 },
  { "--seed",          eslARG_INT,      "42",    NULL,     "n>=0",   NULL,    NULL,  NULL,               "set RNG seed to <n>. Use 0 for a random seed.",                                             1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile> <seqfile> <msaqfile>";
static char banner[] = "realigns sequences by Infernal or HMM";

static int         msaq_banner(FILE *fp, char *progname, char *banner);
static int         MSA_banner(FILE *fp, char *msaname, MSA_STAT *mstat, MSA_STAT *omstat, int nbpairs, int onbpairs);
static int         get_msaname(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int         msa_manipulate(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **msa);
static int         msaq(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int         model_build(struct cfg_s *cfg, char *msafile, char *modfile);
static int         model_align(struct cfg_s *cfg, char *modfile, char *msafile, char *msaqfile);

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
 
  /* banner */
  msaq_banner(stdout, cfg.argv[0], banner);

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_usage(stdout,  cfg.argv[0], usage);
      if (puts("\noptions:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
     exit(0);
    }

  cfg.msafile  = NULL;  
  cfg.msaqfile = NULL;  
  if (esl_opt_ArgNumber(go) != 3) { 
    if (puts("Incorrect number of command line arguments.")       < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
 
  if ((cfg.msafile   = esl_opt_GetArg(go, 1)) == NULL) { 
    if (puts("Failed to get <msafile> argument on command line")  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.seqfile   = esl_opt_GetArg(go, 2)) == NULL) { 
    if (puts("Failed to get <seqqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((cfg.msaqfile  = esl_opt_GetArg(go, 3)) == NULL) { 
    if (puts("Failed to get <msaqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

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
  cfg.msarevmap = NULL;

  /* alphabet */
  cfg.abc = NULL;
  cfg.abcisRNA = FALSE;
  if      (esl_opt_GetBoolean(go, "--rna"))   { cfg.abc = esl_alphabet_Create(eslRNA);   cfg.abcisRNA = TRUE;  }
  else if (esl_opt_GetBoolean(go, "--dna"))   { cfg.abc = esl_alphabet_Create(eslDNA);   cfg.abcisRNA = TRUE;  }
  else if (esl_opt_GetBoolean(go, "--amino")) { cfg.abc = esl_alphabet_Create(eslAMINO);                       }
  
  cfg.watch = esl_stopwatch_Create(); 
  
  cfg.outdir = NULL;
  if (esl_opt_IsOn(go, "--outdir")) esl_sprintf( &cfg.outdir, "%s", esl_opt_GetString(go, "--outdir"));

  /* options */
  cfg.maxnowc     = esl_opt_IsOn(go, "--maxnowc")?    esl_opt_GetReal   (go, "--maxnowc") : -1.0;
  if      (esl_opt_GetBoolean(go, "--opt")) cfg.model = OPT;
  else if (esl_opt_GetBoolean(go, "--hmm")) cfg.model = HMM;
  else if (esl_opt_GetBoolean(go, "--cm"))  cfg.model = CM;

  /* other options */
  cfg.consensus   = esl_opt_IsOn(go, "--consensus")?                                   TRUE : FALSE;
  cfg.maxsq_gsc   = 1000;
  cfg.fragfrac    = esl_opt_IsOn(go, "-F")?           esl_opt_GetReal   (go, "-F")          : 0.0;
  cfg.idthresh    = esl_opt_IsOn(go, "-I")?           esl_opt_GetReal   (go, "-I")          : -1.0;
  cfg.minidthresh = esl_opt_IsOn(go, "-i")?           esl_opt_GetReal   (go, "-i")          : -1.0;
  cfg.nseqmin     = esl_opt_IsOn(go, "--nseqmin")?    esl_opt_GetInteger(go, "--nseqmin")   : -1;
  cfg.gapthresh   = esl_opt_GetReal(go, "--gapthresh");
  cfg.tstart      = esl_opt_IsOn(go, "--tstart")?     esl_opt_GetInteger(go, "--tstart")    : 0;
  cfg.tend        = esl_opt_IsOn(go, "--tend")?       esl_opt_GetInteger(go, "--tend")      : 0;
  cfg.tol         = esl_opt_GetReal   (go, "--tol");
  cfg.verbose     = esl_opt_GetBoolean(go, "-v");
  cfg.onemsa      = esl_opt_IsOn(go, "--onemsa")?     esl_opt_GetBoolean(go, "--onemsa")    : FALSE;
  cfg.singlelink  = esl_opt_GetBoolean(go, "--singlelink");

  if (cfg.minidthresh > cfg. idthresh) esl_fatal("minidthesh has to be smaller than idthresh");

  esl_FileTail(cfg.msafile, TRUE, &cfg.filename);
  if (cfg.outdir) esl_sprintf( &cfg.outheader, "%s/%s.filter", cfg.outdir, cfg.filename);
  else            esl_sprintf( &cfg.outheader, "%s.filter", cfg.outheader);
  esl_sprintf( &cfg.outheader, "%s.I%.2f.G%.2f", cfg.outheader, cfg.idthresh, cfg.gapthresh);
  if (cfg.fragfrac > 0.) esl_sprintf( &cfg.outheader, "%s.F%.2f", cfg.outheader, cfg.fragfrac);

  if (esl_opt_IsOn(go, "--submsa")) { 
    cfg.submsa = esl_opt_GetInteger(go, "--submsa"); 
    esl_sprintf(&outname, "%s.select%d", cfg.outheader, cfg.submsa);  
    free(cfg.outheader); cfg.outheader = NULL;
    esl_sprintf(&cfg.outheader, "%s", outname); 
  }
  else { cfg.submsa = 0; }
  
  cfg.mstat  = NULL;
  cfg.omstat = NULL;

  cfg.plotfile = NULL;

  cfg.treefile  = NULL;
  cfg.treefp    = NULL;
  cfg.T         = NULL;

  cfg.ft   = NULL;
  cfg.fbp  = NULL;
  cfg.fnbp = NULL;

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
msaq_banner(FILE *fp, char *progname, char *banner)
{
  char *name = NULL;
  int   status;

  if ((status = esl_FileTail(progname, FALSE, &name)) != eslOK) return status;

  if (fprintf(fp, "# %s :: %s\n", name, banner)                                               < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");

  if (name) free(name);
  return eslOK;

 ERROR:
  if (name) free(name);
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
    
    status = msa_manipulate(go, &cfg, &msa);
    if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to manipulate alignment"); }
    if (msa == NULL) continue;
    if (msa->alen == 0) {
      MSA_banner(stdout, cfg.msaname, cfg.mstat, cfg.omstat, cfg.nbpairs, cfg.onbpairs);
      continue;
    }
    
    cfg.firstpos = 1;
    if (cfg.makeplot) {
      if (cfg.outdir) esl_sprintf( &cfg.plotfile, "%s/%s", cfg.outdir, cfg.msaname);
      else            esl_sprintf( &cfg.plotfile, "%s", cfg.msaname);
    }
    status = msaq(go, &cfg, msa);
    if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to calculate msaq"); }
    
    if (omsaname) free(omsaname); omsaname = NULL;
    if (useme) free(useme); useme = NULL;
    if (msa) esl_msa_Destroy(msa); msa = NULL;
    if (cfg.msaname) free(cfg.msaname); cfg.msaname = NULL;
    if (cfg.msamap) free(cfg.msamap); cfg.msamap = NULL;
    if (cfg.msarevmap) free(cfg.msarevmap); cfg.msarevmap = NULL;
    if (cfg.omstat) free(cfg.omstat); cfg.omstat = NULL;
    if (cfg.mstat) free(cfg.mstat); cfg.mstat = NULL;
  }

  /* cleanup */
  free(cfg.gnuplot);
  free(cfg.filename);
  esl_stopwatch_Destroy(cfg.watch);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  esl_msafile_Close(afp);
  if (cfg.msaname) free(cfg.msaname);
  if (cfg.treefile) free(cfg.treefile);
  if (cfg.outdir) free(cfg.outdir);
  free(cfg.outheader);
  if (cfg.plotfile)   free(cfg.plotfile);
  if (cfg.ft) free(cfg.ft);
  if (cfg.fbp) free(cfg.fbp);
  if (cfg.fnbp) free(cfg.fnbp);
  if (useme) free(useme);
  if (cfg.msamap) free(cfg.msamap); 
  if (cfg.msarevmap) free(cfg.msarevmap); 
  if (cfg.omstat) free(cfg.omstat);
  if (cfg.mstat) free(cfg.mstat);

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
  else if (msa->name)             esl_sprintf(&cfg->msaname, "%s",      msa->name, type);
  else if (msa->acc)              esl_sprintf(&cfg->msaname, "%s",      msa->acc);
  else if (cfg->onemsa)           esl_sprintf(&cfg->msaname, "%s",      cfg->filename);
  else                            esl_sprintf(&cfg->msaname, "%s_%d", cfg->filename, cfg->nmsa);
  if (esl_opt_IsOn(go, "--submsa")) {					       
    esl_sprintf(&submsaname, "%s.select%d", cfg->msaname, esl_opt_GetInteger(go, "--submsa"));
    free(cfg->msaname); cfg->msaname = NULL;
    esl_sprintf(&cfg->msaname, "%s", submsaname);
  }

  /* add tstat..tend information */
  if (cfg->tstart > 0 || cfg->tend > 0) 
    esl_sprintf(&cfg->msaname, "%s_%d-%d", cfg->msaname, cfg->tstart, cfg->tend);
  
  if (tok) free(tok);
  if (type) free(type);
  if (submsaname) free(submsaname);
  return eslOK;
}

static int
msa_manipulate(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA **omsa)
{
  ESL_MSA *msa = *omsa;
  int     *useme = NULL;
  char    *msg = "original_msa_manipulate failed";
  char    *type = NULL;
  char    *tok = NULL;
  char    *submsaname = NULL;
  int64_t  alen = (*omsa)->alen;
  int64_t  startpos, endpos;
  int      seq_cons_len = 0;
  int      nremoved = 0;	  /* # of identical sequences removed */
  int      nfrags = 0;	          /* # of fragments removed */
  int      pos;

  printf("MSA          %s (alen = %" PRIu64 " nseq = %d)\n", cfg->filename, alen, (*omsa)->nseq);
  printf("ID  %.4f G   %.4f ", cfg->idthresh, cfg->gapthresh);
  if (cfg->fragfrac > 0) printf("F   %.2f\n", cfg->fragfrac); else printf("\n");

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
  
  /* apply msa filters and than select submsa
   */
  if (cfg->fragfrac >= 0.    && msamanip_RemoveFragments(cfg->fragfrac, omsa, &nfrags, &seq_cons_len)                 != eslOK) {
    printf("%s\nremove_fragments failed\n", cfg->errbuf);                 esl_fatal(msg); }
  if (esl_opt_IsOn(go, "-I") && msamanip_SelectSubsetBymaxID(cfg->r, omsa, cfg->idthresh, cfg->singlelink, &nremoved) != eslOK) {
    printf("%s\n", cfg->errbuf); printf("select_subsetBymaxID failed\n"); esl_fatal(msg); }
  if (esl_opt_IsOn(go, "-i") && msamanip_SelectSubsetByminID(cfg->r, omsa, cfg->minidthresh, &nremoved)               != eslOK) {
    printf("%s\n", cfg->errbuf); printf("select_subsetByminID failed\n"); esl_fatal(msg); }
  if (cfg->submsa            && msamanip_SelectSubset(cfg->r, cfg->submsa, omsa, NULL, cfg->errbuf, cfg->verbose)     != eslOK) {
    printf("%s\n", cfg->errbuf);                                          esl_fatal(msg); }

  /* What is left after all the filtering?
   */
  msa = *omsa;  
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

  /* Now apply [tstart,tend] restriction if given
   */
  if (msamanip_Truncate(msa, cfg->tstart, cfg->tend, &startpos, &endpos, cfg->errbuf) != eslOK) {
    printf("%s\nTruncate failed\n", cfg->errbuf); esl_fatal(msg); }

  /* remove columns with gaps.
   * Important: the mapping is done here; cannot remove any other columns beyond this point.
   */
  if (cfg->consensus) {
    if (msamanip_SelectConsensus(msa, &useme, cfg->verbose) != eslOK) {
      printf("%s\nconsensus selection fails\n", cfg->errbuf); esl_fatal(msg);
    }
  }
  if (msamanip_RemoveGapColumns(cfg->gapthresh, msa, startpos, endpos, alen, &cfg->msamap, &cfg->msarevmap, useme, cfg->errbuf, cfg->verbose) != eslOK) {
    printf("%s\n", cfg->errbuf); esl_fatal(msg);
  }
  msamanip_ConvertDegen2N(msa);

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
msaq(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  char     tmpmsafile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char     tmpmodfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  FILE    *msafp = NULL;
  int     *ct = NULL;
  int      K = msa->abc->K;
  int      status;
  
  /* the ct vector  */
  status = msamanip_CalculateCT(msa, &ct, &cfg->nbpairs, cfg->maxnowc, cfg->errbuf);

  /* Print some alignment information */
  MSA_banner(stdout, cfg->msaname, cfg->mstat, cfg->omstat, cfg->nbpairs, cfg->onbpairs);
  
  /* write the msa */
  if ((status = esl_tmpfile_named(tmpmsafile, &msafp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create tmpmsafile");
  esl_msafile_Write(msafp, msa, eslMSAFILE_STOCKHOLM);
  fclose(msafp);

  /* make cm/hmm */
  model_build(cfg, tmpmsafile, tmpmodfile);
  
  /* align sequences to cm/hmm */
  model_align(cfg, tmpmodfile, cfg->seqfile, cfg->msaqfile);

  /* free and return */
  remove(tmpmsafile);
  remove(tmpmodfile);
  free(ct);
  return eslOK;
  
 ERROR:
  if (ct) free(ct);
  if (msa) esl_msa_Destroy(msa);
  remove(tmpmsafile);
  remove(tmpmodfile);
  return status;
}

static int
model_build(struct cfg_s *cfg, char *msafile, char *modfile)
{
  FILE *modfp = NULL;
  char *cmd   = NULL;
  char *args  = NULL;
  char *s     = NULL;
  int   status;

  if      (cfg->model == HMM || (cfg->model == OPT && cfg->nbpairs == 0) ) { // build an HMM
    if (RSCAPE_BIN) esl_sprintf(&cmd, "%s/hmmer/bin/hmmbuild", RSCAPE_BIN);  
    else            ESL_XFAIL(status, cfg->errbuf, "Failed to find hmmbuild executable\n");
  }
  else if (cfg->model == CM  || (cfg->model == OPT && cfg->nbpairs >  0) ) { // build an CM model
    if (RSCAPE_BIN) esl_sprintf(&cmd, "%s/infernal-1.1.2/bin/cmbuild -F", RSCAPE_BIN);  
    else            ESL_XFAIL(status, cfg->errbuf, "Failed to find cmbuild executable\n");
  }
  
  if ((status = esl_tmpfile_named(modfile, &modfp)) != eslOK) ESL_XFAIL(status, cfg->errbuf, "failed to create modfile");
  fclose(modfp);
  esl_sprintf(&args, "%s %s %s", cmd, modfile, msafile);
  if (cfg->verbose) { printf("%s\n", args); }
  system(args);
    
  if (cmd  != NULL) free(cmd);
  if (args != NULL) free(args);
  return eslOK;

 ERROR:
  if (cmd  != NULL) free(cmd);
  if (args != NULL) free(args);
  return status;
}

static int
model_align(struct cfg_s *cfg, char *modfile, char *msafile, char *msaqfile)
{
  char *sfile = NULL;
  char *cmd   = NULL;
  char *args  = NULL;
  char *s     = NULL;
  int   status;
  
  if      (cfg->model == HMM || (cfg->model == OPT && cfg->nbpairs == 0) ) { // it is an HMM
    if (RSCAPE_BIN) esl_sprintf(&cmd, "%s/hmmer/bin/hmmalign", RSCAPE_BIN);  
    else            ESL_XFAIL(status, cfg->errbuf, "Failed to find hmmalign executable\n");
    esl_sprintf(&args, "%s %s %s > %s", cmd, modfile, msafile, msaqfile);
  }
  else if (cfg->model == CM  || (cfg->model == OPT && cfg->nbpairs >  0) ) { // it is a infernal CM model
    esl_sprintf(&sfile, "%s.sfile", msaqfile);
    if (RSCAPE_BIN) esl_sprintf(&cmd, "%s/infernal-1.1.2/bin/cmalign --sfile %s -o %s", RSCAPE_BIN, sfile, msaqfile);  
    else            ESL_XFAIL(status, cfg->errbuf, "Failed to find cmbuild executable\n");
    esl_sprintf(&args, "%s %s %s > /tmp/null", cmd, modfile, msafile);
  }

  if (1||cfg->verbose) { printf("%s\n", args); }
  system(args);
  
  if (cmd  != NULL) free(cmd);
  if (args != NULL) free(args);
  return eslOK;

 ERROR:
  if (cmd  != NULL) free(cmd);
  if (args != NULL) free(args);
  return status;
}
