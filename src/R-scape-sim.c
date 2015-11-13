/* R-scape-sim -- simulate alignments to test R-scape
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

#include "e1_rate.h"

#include "msamanip.h"
#include "msatree.h"
#include "covariation.h"
#include "covgrammars.h"
#include "ribosum_matrix.h"
#include "cov_simulate.h"

#define TREEOPTS   "--star,--given,--sim"                                          

/* Exclusive options for evolutionary model choice */

/* struct cfg_s : "Global" application configuration shared by all threads/processes.
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s { /* Shared configuration in masters & workers */
  int                  argc;
  char               **argv;
  ESL_STOPWATCH       *watch;
  char                 errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS      *r;	               /* random numbers */
  ESL_ALPHABET        *abc;                    /* the alphabet */
  
  int                  onemsa;
  int                  nmsa;
  char                *msafile;
  char                *msaname;
  
  char                *outdir;
  char                *filename;
  char                *outheader;              // header for all output files 
  char                *simsafile;
  FILE                *simsafp;
  int                  infmt;
  
  int                  N;                      // number of sequences in the alignment
  double               target_abl;             // average branch length (in number of changes per site)
  double               target_atbl;            // average total branch length (in number of changes per site)
  double               abl;                    // average branch length (in number of changes per site)
  double               atbl;                   // average total branch length (in number of changes per site)
  TREETYPE             treetype;               // star, given or simulated tree topology
  ESL_TREE            *T;
  
  int                  noss;                   // asume unstructured, do not use the given secondary structure if any
  int                  noindels;               // make ungapped alignmenst
  
  MSA_STAT            *mstat;                  // statistics of the given alignment 
  MSA_STAT            *simstat;                // statistics of the simulated alignment 
  int                 *ct;
  int                  nbpairs;
  int                  simnbpairs;
  int                  usesq;                   // use a specific seq in the original msa as root

  EVOM                 evomodel;
  char                *paramfile;
  struct rateparam_s   rateparam;
  E1_RATE             *e1rate;
  E1_RATE             *e1rateB;                 // only difference with e1rate: rate of deletion of ancestral residues is lower.
                                                // used to delete basepaired residues
  ESL_DMATRIX         *R1;                      // 4x4 rate matrix
  double              *bg;                      // background frequencies

  char                *ribofile;
  struct ribomatrix_s *ribosum;
  float                tol;
  int                  verbose;
};

static ESL_OPTIONS options[] = {
  /* name             type              default  env        range    toggles  reqs   incomp              help                                                                                  docgroup*/
  { "-h",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "show brief help on version and usage",                                                      1 },
  { "-v",             eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "be verbose",                                                                                1 },
  /* parameters to control the simulation */
  { "-N",              eslARG_INT,       "40",   NULL,      "n>1",   NULL,    NULL,  NULL,               "number of sequences in the simulated msa",                                                  0 }, 
  { "--abl",           eslARG_REAL,      NULL,   NULL,      "x>0",   NULL,    NULL,"--atbl",             "tree average branch length in number of changes per site",                                  0 }, 
  { "--atbl",          eslARG_REAL,    "0.60",   NULL,      "x>0",   NULL,    NULL,"--abl",              "tree average total branch length in number of changes per site",                            0 }, 
  { "--noss",          eslARG_NONE,     FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "assume unstructured, even if msa has a given ss_cons",                                      0 }, 
  { "--noindels",      eslARG_NONE,     FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,               "produces ungapped alignments",                                                              0 }, 
  { "--star",          eslARG_NONE,     FALSE,   NULL,       NULL,  TREEOPTS, NULL,  NULL,               "star topology",                                                                             0 },
  { "--rand",          eslARG_NONE,     FALSE,   NULL,       NULL,  TREEOPTS, NULL,  NULL,               "independent sequences",                                                                     0 },
  { "--given",         eslARG_NONE,     FALSE,   NULL,       NULL,  TREEOPTS, NULL,  NULL,               "given msa topology",                                                                        0 },
  { "--sim",           eslARG_NONE,      TRUE,   NULL,       NULL,  TREEOPTS, NULL,  NULL,               "simulated topology",                                                                        0 },
  { "--usesq",          eslARG_INT,      NULL,   NULL,      "n>=1",  NULL,    NULL,  NULL,               "sq from the origional msa used as root (default random)",                                   0 }, 
 /* options for input msa (if seqs are given as a reference msa) */
  { "--informat",   eslARG_STRING,       NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "specify format",                                                                            1 },
   /* Control of scoring system - substitutions */ 
  { "--mxfile",     eslARG_INFILE,       NULL,   NULL,       NULL,   NULL,    NULL,  NULL,               "read substitution rate matrix from file <f>",                                              0 },
  /* Control of scoring system - indels */ 
  { "--evomodel",     eslARG_STRING,    "AIF",   NULL,       NULL,   NULL,    NULL,  NULL,                "evolutionary model used",                                                                 0 },
  /* Control of scoring system - ribosum */
  { "--ribofile",     eslARG_INFILE,     NULL,   NULL,       NULL,   NULL,    NULL,  NULL,                "read ribosum structure from file <f>",                                                    0 },
  /* Control of output */
  { "--outdir",     eslARG_STRING,       NULL,   NULL,       NULL,   NULL,    NULL,  NULL,                "specify a directory for all output files",                                                1 },
  { "-o",          eslARG_OUTFILE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,                "send output to file <f>, not stdout",                                                     1 },
  /* other options */  
  { "--onemsa",       eslARG_NONE,      FALSE,   NULL,       NULL,   NULL,    NULL,  NULL,                "if file has more than one msa, analyze only the first one",                               1 },
  { "--tol",          eslARG_REAL,     "1e-3",   NULL,       NULL,   NULL,    NULL,  NULL,                "tolerance",                                                                               0 },
  { "--seed",          eslARG_INT,       "0",    NULL,     "n>=0",   NULL,    NULL,  NULL,                "set RNG seed to <n>",                                                                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msa>";
static char banner[] = "R-scape-sim - synthetic alignments to test R-scape";

static int MSA_banner(FILE *fp, char *msaname, MSA_STAT *mstat, MSA_STAT *omstat, int nbpairs, int onbpairs);
static int get_msaname(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);
static int simulate_msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA **simmsa);
static int create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa);

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static int process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, struct cfg_s *ret_cfg)
{
  ESL_GETOPTS  *go = esl_getopts_Create(options);
  struct cfg_s  cfg;
  int           status;

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
  
  /* outheader for all output files */
  cfg.outheader = NULL;
  msamanip_OutfileHeader(cfg.msafile, &cfg.outheader); 
  
  /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = eslx_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  cfg.nmsa = 0;
  cfg.msaname = NULL;
  
  /* alphabet */
  cfg.abc = esl_alphabet_Create(eslRNA);
  esl_alphabet_SetEquiv(cfg.abc, '=', '-');     /* allow = as a gap character too */
  esl_alphabet_SetEquiv(cfg.abc, '.', '-');     /* allow . as a gap character too */
  
  cfg.watch = esl_stopwatch_Create(); 
  
  cfg.outdir = NULL;
  if (esl_opt_IsOn(go, "--outdir")) esl_sprintf( &cfg.outdir, "%s", esl_opt_GetString(go, "--outdir"));
  
  esl_FileTail(cfg.msafile, TRUE, &cfg.filename);
  if ( cfg.outdir ) esl_sprintf( &cfg.outheader, "%s/%s", cfg.outdir, cfg.filename);
  
  /* parameters of the simulation */
  cfg.N        = esl_opt_GetInteger(go, "-N");
  cfg.noss     = esl_opt_GetBoolean(go, "--noss");
  cfg.noindels = esl_opt_GetBoolean(go, "--noindels");
  if      (esl_opt_GetBoolean  (go, "--star"))  cfg.treetype = STAR;
  else if (esl_opt_GetBoolean  (go, "--given")) cfg.treetype = GIVEN;
  else if (esl_opt_GetBoolean  (go, "--sim"))   cfg.treetype = SIM; 
  else if (esl_opt_GetBoolean  (go, "--rand"))  cfg.treetype = RAND; 

  cfg.usesq       = esl_opt_IsOn(go, "--usesq")? esl_opt_GetInteger(go, "--usesq") : -1.0;
  cfg.target_atbl = esl_opt_GetReal(go, "--atbl");
  cfg.target_abl  = esl_opt_IsOn(go, "--abl")?  esl_opt_GetReal(go, "--abl")  : -1.0;
  cfg.abl  = -1.0;
  cfg.atbl = -1.0;
  
  /* other options */
  cfg.onemsa  = esl_opt_IsOn(go, "--onemsa")?     esl_opt_GetBoolean(go, "--onemsa")    : FALSE;
  cfg.tol     = esl_opt_GetReal   (go, "--tol");
  cfg.verbose = esl_opt_GetBoolean(go, "-v");

  /* file with the simulated msa */
  cfg.simsafile = NULL;
  cfg.simsafp   = NULL;
  if ( esl_opt_IsOn(go, "-o") ) 
    esl_sprintf(&cfg.simsafile, "%s", esl_opt_GetString(go, "-o"));
  else
    esl_sprintf(&cfg.simsafile, "%s_synthetic.sto", cfg.filename);

  if ((cfg.simsafp = fopen(cfg.simsafile, "w")) == NULL) esl_fatal("Failed to open simsa file %s", cfg.simsafile);
  
  cfg.T  = NULL;
  cfg.ct = NULL;
  cfg.mstat   = NULL;
  cfg.simstat = NULL;
  cfg.nbpairs    = 0;
  cfg.simnbpairs = 0;

  /* the evolutionary model */
  cfg.evomodel = e1_rate_Evomodel(esl_opt_GetString(go, "--evomodel"));
  
  /* the paramfile  */
  cfg.paramfile = NULL;
  if (cfg.evomodel == AIF)
    esl_sprintf(&cfg.paramfile, "lib/evoparam/Pfam.seed.S1000.trainGD.AIF.param");
  else if (cfg.evomodel == AFG)
    esl_sprintf(&cfg.paramfile, "lib/evoparam/Pfam.seed.S1000.trainGD.AFG.param");
  else esl_fatal("could not identify evomodel");
  status = e1_rate_ReadParamfile(cfg.paramfile, &cfg.rateparam, &cfg.evomodel, cfg.errbuf, cfg.verbose);
  if (status != eslOK) esl_fatal("Failed to read paramfile %s\n%s", cfg.paramfile, cfg.errbuf);

  /* the ribosum matrices */
  cfg.ribofile = NULL;
  cfg.ribosum  = NULL;
  if ( esl_opt_IsOn(go, "--ribofile") ) { cfg.ribofile = esl_opt_GetString(go, "--ribofile"); }
  else esl_sprintf(&cfg.ribofile, "lib/ribosum/ssu-lsu.final.er.ribosum");
  cfg.ribosum = Ribosum_matrix_Read(cfg.ribofile, cfg.abc, FALSE, cfg.errbuf);
  if (cfg.ribosum == NULL) esl_fatal("%s\nfailed to create ribosum matrices from file %s\n", cfg.errbuf, cfg.ribofile);
  if (cfg.verbose) Ribosum_matrix_Write(stdout, cfg.ribosum);

  /* the e1_rate  */
  cfg.e1rate = NULL;
  cfg.e1rateB = NULL;
  cfg.R1 = cfg.ribosum->xrnaQ;
  cfg.bg = cfg.ribosum->bg;

  cfg.e1rate = e1_rate_CreateWithValues(cfg.abc, cfg.evomodel, cfg.rateparam, NULL, cfg.R1, cfg.bg, TRUE, cfg.tol, cfg.errbuf, cfg.verbose);
  cfg.e1rate->evomodel = cfg.evomodel; 
  if (cfg.e1rate == NULL) { printf("%s. bad rate model\n", cfg.errbuf); esl_fatal("Failed to create e1rate"); }
  if (cfg.verbose) {
    e1_rate_Dump(stdout, cfg.e1rate);
    esl_dmatrix_Dump(stdout, cfg.e1rate->em->Qstar, "ACGU", "ACGU");
    esl_dmatrix_Dump(stdout, cfg.e1rate->em->E, "ACGU", "ACGU");
    esl_vec_DDump(stdout, cfg.e1rate->em->f, cfg.abc->K, "ACGU");
  }
  // for e1rateB, lower the rate of ancestral deletions
  cfg.rateparam.muAM /= 2.0;
  cfg.e1rateB = e1_rate_CreateWithValues(cfg.abc, cfg.evomodel, cfg.rateparam, NULL, cfg.R1, cfg.bg, TRUE, cfg.tol, cfg.errbuf, cfg.verbose);
  cfg.e1rateB->evomodel = cfg.evomodel; 
  if (cfg.e1rateB == NULL) { printf("%s. bad rate model\n", cfg.errbuf); esl_fatal("Failed to create e1rateB"); }
  if (cfg.verbose) {
    e1_rate_Dump(stdout, cfg.e1rateB);
    esl_dmatrix_Dump(stdout, cfg.e1rateB->em->Qstar, "ACGU", "ACGU");
    esl_dmatrix_Dump(stdout, cfg.e1rateB->em->E, "ACGU", "ACGU");
    esl_vec_DDump(stdout, cfg.e1rateB->em->f, cfg.abc->K, "ACGU");
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
MSA_banner (FILE *fp, char *msaname, MSA_STAT *simstat, MSA_STAT *mstat, int simnbpairs, int nbpairs)
{
  if (simstat) 
    fprintf(fp, "# simMSA %s nseq %d (%d) alen %" PRId64 " (%" PRId64 ") avgid %.2f (%.2f) nbpairs %d (%d)\n", 
	    msaname, simstat->nseq, mstat->nseq, simstat->alen, mstat->alen, 
	    simstat->avgid, mstat->avgid, simnbpairs, nbpairs);
  else 
    fprintf(fp, "# givenMSA %s nseq %d alen %" PRId64 " avgid %.2f nbpairs %d\n", 
	    msaname, mstat->nseq, mstat->alen, mstat->avgid, nbpairs);

  return eslOK;
}

int
main(int argc, char **argv)
{ 
  ESL_GETOPTS     *go;
  struct cfg_s     cfg;
  ESLX_MSAFILE    *afp = NULL;
  ESL_MSA         *msa = NULL;            /* the input alignment    */
  ESL_MSA         *simsa = NULL;          /* the simulated alignment    */
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
    
    if (cfg.N < msa->nseq && cfg.treetype == GIVEN) {
      status = msamanip_SelectSubset(cfg.r, cfg.N, &msa, NULL, cfg.errbuf, cfg.verbose);
      if (status != eslOK) {
	printf("%s\n", cfg.errbuf);              
	esl_fatal("Failed to manipulate original alignment"); 
      }
      if (cfg.verbose) eslx_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM);
    }

    /* remove degenacies */
    msamanip_ConvertDegen2RandomCanonical(cfg.r, msa);

    /* the msaname */
    status = get_msaname(go, &cfg, msa);
    if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to manipulate original alignment"); }
    
    /* stats of the given alignment */
    msamanip_XStats(msa, &cfg.mstat);
    msamanip_CalculateCT(msa, NULL, &cfg.nbpairs, cfg.errbuf);
 
    status = simulate_msa(go, &cfg, msa, &simsa);
    if (status != eslOK)  { printf("%s\n", cfg.errbuf); esl_fatal("Failed to simulate msa"); }
    if (simsa == NULL)    { printf("%s\n", cfg.errbuf); esl_fatal("Failed to create msa"); }
    
    /* stats of the simulated alignment */
    msamanip_XStats(simsa, &cfg.simstat);
    msamanip_CalculateCT(simsa, NULL, &cfg.simnbpairs, cfg.errbuf);
 
    /* write the simulated msa to file */
    if (cfg.simsafp && simsa) eslx_msafile_Write(cfg.simsafp, simsa, eslMSAFILE_STOCKHOLM);
    if (1||cfg.verbose) 
      MSA_banner(stdout, cfg.msaname, cfg.simstat, cfg.mstat, cfg.simnbpairs, cfg.nbpairs);
    if (cfg.verbose) 
      eslx_msafile_Write(stdout, simsa, eslMSAFILE_STOCKHOLM);
  
    if (msa) esl_msa_Destroy(msa); msa = NULL;
    if (simsa) esl_msa_Destroy(simsa); simsa = NULL;
    if (cfg.msaname) free(cfg.msaname); cfg.msaname = NULL;
    if (cfg.mstat) free(cfg.mstat); cfg.mstat = NULL;
    if (cfg.simstat) free(cfg.simstat); cfg.simstat = NULL;
    if (cfg.T) esl_tree_Destroy(cfg.T); cfg.T = NULL;
  }

  /* cleanup */
  if (cfg.paramfile) free(cfg.paramfile);
  if (cfg.filename) free(cfg.filename);
  if (cfg.simsafp) fclose(cfg.simsafp);
  esl_stopwatch_Destroy(cfg.watch);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(cfg.r);
  if (cfg.ct) free(cfg.ct);
  eslx_msafile_Close(afp);
  if (cfg.msaname) free(cfg.msaname);
  free(cfg.outheader);
  if (cfg.e1rate)  e1_rate_Destroy(cfg.e1rate);
  if (cfg.e1rateB) e1_rate_Destroy(cfg.e1rateB);
  if (cfg.ribofile) free(cfg.ribofile);
  if (cfg.ribosum) Ribosum_matrix_Destroy(cfg.ribosum);
  if (cfg.simsafile) free(cfg.simsafile);
  if (cfg.mstat) free(cfg.mstat); 
  if (cfg.simstat) free(cfg.simstat); 

  return 0;
}

static int
get_msaname(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  char *msg = "get_msaname failed";
  char *type = NULL;
  char *tp;
  char *tok1;
  char *tok2;
  char *tok = NULL;
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
  if      (msa->acc  && msa->name && type) esl_sprintf(&cfg->msaname, "%s_%s%s", msa->acc, msa->name, type);
  else if (msa->acc  && msa->name)         esl_sprintf(&cfg->msaname, "%s_%s",   msa->acc, msa->name);
  else if (msa->name && type)              esl_sprintf(&cfg->msaname, "%s%s",    msa->name, type);
  else if (msa->acc)                       esl_sprintf(&cfg->msaname, "%s",      msa->acc);
  else if (msa->name)                      esl_sprintf(&cfg->msaname, "%s",      msa->name);
  else if (cfg->onemsa)                    esl_sprintf(&cfg->msaname, "%s",      cfg->filename);
  else                                     esl_sprintf(&cfg->msaname, "%s_%d",   cfg->filename, cfg->nmsa);

  if (tok) free(tok);
  if (type) free(type);
  return eslOK;
}


static int
simulate_msa(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA **ret_simsa)
{
  ESL_MSA *root    = NULL;    /* the ancestral sq and structure */
  ESL_MSA *msafull = NULL;    /* alignment of leaves and internal node sequences */
  ESL_MSA *simsa   = NULL;    /* alignment of leave sequences */
  char    *rootname = NULL;
  char    *rootdesc = NULL;
  int     *useme = NULL;
  int      root_bp;
  int      usesq;
  int      i;
  
  if (msa == NULL) return eslOK;
  if (cfg->treetype == GIVEN) cfg->N = msa->nseq;
  if (cfg->N <= 1) return eslOK;

  // create the tree with target abl
  create_tree(go, cfg, msa);
  
  // sample the ancestral sequence from msa, adjust the secondary structure
  useme = malloc(msa->nseq * sizeof(int));
  esl_vec_ISet(useme, msa->nseq, FALSE);
  usesq = (cfg->usesq > 0 && cfg->usesq < msa->nseq)? cfg->usesq-1 : (int)(esl_random(cfg->r)*msa->nseq);
  useme[usesq] = TRUE;

  esl_msa_SequenceSubset(msa, useme, &root);
  if (esl_msa_MinimGaps(root, NULL, "-", FALSE) != eslOK) 
    esl_fatal("failed to generate the root sequence");
  esl_sprintf(&rootname, "%s", cfg->filename);
  esl_sprintf(&rootdesc, "%s-%s", cfg->filename, root->sqname[0]);
  esl_msa_SetName(root, rootname, -1);
  esl_msa_SetDesc(root, rootdesc, -1);
  free(useme); useme = NULL;
 
  if (cfg->verbose) {
    msamanip_CalculateCT(root, NULL, &root_bp, cfg->errbuf);
    printf("root sequence L = %d nbps = %d \n", (int)root->alen, root_bp);
    eslx_msafile_Write(stdout, root, eslMSAFILE_STOCKHOLM);
 }

  /* Generate the simulated alignment */
  if (cov_GenerateAlignment(cfg->r, cfg->treetype, cfg->N, cfg->atbl, cfg->T, root, cfg->e1rate, cfg->e1rateB, cfg->ribosum, &msafull, 
			    cfg->noss, cfg->noindels, cfg->tol, cfg->errbuf, cfg->verbose) != eslOK)
    esl_fatal("%s\nfailed to generate the simulated alignment", cfg->errbuf);
  if (msafull == NULL) esl_fatal("%s\nfailed to generate the simulated alignment", cfg->errbuf);

  /* The leaves-only alignment */
  useme = malloc(msafull->nseq * sizeof(int));
  for (i = 0; i < msafull->nseq; i++) {
    if (strncmp(msafull->sqname[i], "v", 1) == 0) useme[i] = FALSE; 
    else                                          useme[i] = TRUE;
  }
  if (esl_msa_SequenceSubset(msafull, useme, &simsa) != eslOK)
    esl_fatal("failed to generate leaf alignment");
  if (esl_msa_MinimGaps(simsa, NULL, "-", FALSE) != eslOK) 
    esl_fatal("failed to remove gaps alignment");
  
 if (1||cfg->verbose) 
    eslx_msafile_Write(stdout, simsa, eslMSAFILE_STOCKHOLM);

  *ret_simsa = simsa;

  free(useme);
  free(rootname);
  free(rootdesc);
  esl_msa_Destroy(root);
  if (msafull) esl_msa_Destroy(msafull);
  return eslOK;
}

static int
create_tree(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_MSA *msa)
{
  char *msg = "bad tree";
  int   status;
  
  /* the TREE */
  if (cfg->treetype == RAND || cfg->treetype == STAR) { // sequences are independent or sequences are independent but derived form a common ancestor
    cfg->T = NULL;
    cfg->atbl = cfg->abl = cfg->target_atbl;

  }
  else if (cfg->treetype == GIVEN) {   // use the tree determined by the given alignment
    status = Tree_CalculateExtFromMSA(msa, &cfg->T, TRUE, cfg->errbuf, cfg->verbose);
    if (status != eslOK) { printf("%s\n", cfg->errbuf); esl_fatal(cfg->errbuf); }
  }
  else if (cfg->treetype == SIM) {  // generate random ultrametric tree T
    if (esl_tree_Simulate(cfg->r, cfg->N, &cfg->T) != eslOK) esl_fatal(msg);
    if (esl_tree_SetTaxaParents(cfg->T)            != eslOK) esl_fatal(msg);
    if (esl_tree_Validate(cfg->T, NULL)            != eslOK) esl_fatal(msg);
  }

  /* scale the tree */
  if (cfg->T) {
    if (cfg->target_abl > 0) {
      if (esl_tree_er_RescaleAverageBL(cfg->target_abl, cfg->T, cfg->tol, cfg->errbuf, cfg->verbose) != eslOK) esl_fatal(msg);
    }
    else if (cfg->target_atbl > 0) {
      if (esl_tree_er_RescaleAverageTotalBL(cfg->target_atbl, cfg->T, cfg->tol, cfg->errbuf, cfg->verbose) != eslOK) esl_fatal(msg);
    }
    
    cfg->abl = esl_tree_er_AverageBL(cfg->T);
    Tree_GetNodeTime(0, cfg->T, &cfg->atbl, NULL, NULL, cfg->errbuf, cfg->verbose);
    if (cfg->verbose) printf("# average leave-to-root length: %f average branch length: %f\n", cfg->atbl, cfg->abl);
    if (cfg->verbose) Tree_Dump(stdout, cfg->T, "Tree");
  }
  else if (cfg->verbose) printf("# average leave-to-root length: %f \n", cfg->atbl);
  
  return eslOK;
}


