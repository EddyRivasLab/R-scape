/* e2Dickerson -- given an msa (or a collection of sequences and a tree)
 *                estimate the pairwise rates of molecular evolution
 *                using either a simple counting method or the e2 method
 *                based on an E2 affine evolutionary mode.
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
#include "esl_msacluster.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2Dickerson.h"
#include "e2.h"
#include "e2_profilesq.h"
#include "e1_rate.h"
#include "branchinference.h"
#include "fetchpfamdb.h"
#include "msatree.h"
#include "Dickersonrates.h"
#include "mya.h"
#include "msamanip.h"
#include "orthologs.h"
#include "plot.h"
#include "ssifile.h"

#define ALPHOPTS "--amino,--dna,--rna"                                                                                          /* Exclusive options for alphabet choice */
#define TAXOPTS  "--superkingdom,--kingdom,--phylum,--class,--order,--family,--genus"                                           /* Exclusive options for taxonomy choice to group sequences */                     
#define MYAOPTS  "--domya,--fromfile,--onthespot"                                                                               /* Exclusive options for getting mya for species */                     
#define PLOTOPTS "--serial,--together,--joint,--plotall"                                                                        /* Exclusive options for plotting rates oortholog clusters from a given family */                     
#define ALIOPTS  "--e2,--e2f,--e2hmmer,--e2fhmmer"                                                                              /* Exclusive options for alignment choice */
#define EVOMOPTS "--TKF91,--TKF92,--FID,--GTKF92,--GRTKF92,--ITKF92,--LI,--LR,--AFG,--AFGR,--AFR,--AIF,--GG,--G2"              /* Exclusive options for evolutionary model choice */

static ESL_OPTIONS options[] = {
  /* name           type                                      default   env  range    toggles   reqs   incomp               help                                                                 docgroup*/
  { "-h",           eslARG_NONE,                               FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "show brief help on version and usage",                                     0 },
  { "-v",           eslARG_NONE,                               FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "be verbose?",                                                              0 },
  { "--makeT",      eslARG_NONE,                               FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "make the tree from the msa using fasttree",                                0 },                    
  { "--group",      eslARG_INT,                                 "10",  NULL, "n>0",    NULL,    NULL,  NULL,               "the approximate number of groups that the species should be grouped into", 0 }, { "--noopt",      eslARG_NONE,                               FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "TRUE for no optimization of parameters",                                   0 }, 
   /* options for msa */
  { "-F",           eslARG_REAL,                              "0.80",  NULL,"0<x<=1.0",NULL,    NULL,  NULL,               "filter out seqs <x*seq_cons residues",                                                      0 },
  { "-I",           eslARG_REAL,                             "0.999",  NULL,"0<x<=1.0",NULL,    NULL,  NULL,               "require seqs to have < x id",                                                               0 },
  { "--full",       eslARG_NONE,                               FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "TRUE if a full alignment, default is a seed alignment",                                     0 },
  { "--submsa",     eslARG_INT,                                FALSE,  NULL, "n>0",    NULL,    NULL,  NULL,               "take n random sequences from the alignment to continue the analysis, all fif NULL",         0 },
  { "--msafile",    eslARG_NONE,                               FALSE,  NULL,  NULL,    NULL,    NULL,  NULL,               "main argument is an msa instead of a pfam name",                                            0 },
  { "--local",      eslARG_NONE,                               FALSE,  NULL,   NULL,   NULL,    NULL,  "--e2f,--e2fhmmer", "TRUE for a local alignment, default is global",                                             0 }, 
  { "--e2",         eslARG_NONE,                                TRUE,  NULL,  NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, realign pairs of sequences",                                          0 }, 
  { "--e2f",        eslARG_NONE,                               FALSE,  NULL,  NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, use alignment as is, do not realign",                                 0 }, 
  { "--e2hmmer",    eslARG_STRING,                             FALSE,  NULL,  NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, realign pairs of sequences, use position-specific evomodel",          0 }, 
  { "--e2fhmmer",   eslARG_STRING,                             FALSE,  NULL,  NULL,    ALIOPTS, NULL,  NULL,               "for the e2 algorithm, use alignment as is, do not realign, use position-specific evomodel", 0 }, 
  /* options for mya */
  { "--domya",      eslARG_NONE,                               FALSE,  NULL,  NULL, MYAOPTS,    NULL,  NULL,               "calculate time in mya for all comparisons, store in distfile",             0 },
  { "--fromfile",   eslARG_NONE,                           "default",  NULL,  NULL, MYAOPTS,    NULL,  NULL,               "distfile exist",                                                           0 },
  { "--onthespot",  eslARG_NONE,                               FALSE,  NULL,  NULL, MYAOPTS,    NULL,  NULL,               "calculate time in mya on the spot for each comparison",                    0 },
  /* options for orthologs */  
  { "--nboot",      eslARG_INT,                                "100",   NULL, "n>0",    NULL,    NULL,  NULL,               "number of boostrapped genes trees",                                        0 },
  { "--reboot",     eslARG_NONE,                               FALSE,   NULL,  NULL,    NULL,    NULL,  NULL,               "refresh the treeboot file even if it exists",                              0 },
  { "--clsize",     eslARG_INT,                                  "2",   NULL, "n>0",    NULL,    NULL,  NULL,               "minimum size of ortholog clusters to analyze",                             0 },
  { "--clmax",      eslARG_NONE,                               FALSE,   NULL,  NULL,    NULL,    NULL,  "--clsize",         "use the maximal cluster of orthologs (superseeds --clsize)",               0 },
  /* alternative options for taxonomy type. Supersedes group */
  { "--superkingdom", eslARG_NONE,                              FALSE,  NULL,  NULL, TAXOPTS,    NULL,  NULL,               "taxonmy type to group sequences by",                                       0 },
  { "--kingdom",      eslARG_NONE,                              FALSE,  NULL,  NULL, TAXOPTS,    NULL,  NULL,               "taxonmy type to group sequences by",                                       0 },
  { "--phylum",       eslARG_NONE,                              FALSE,  NULL,  NULL, TAXOPTS,    NULL,  NULL,               "taxonmy type to group sequences by",                                       0 },
  { "--class",        eslARG_NONE,                              FALSE,  NULL,  NULL, TAXOPTS,    NULL,  NULL,               "taxonmy type to group sequences by",                                       0 },
  { "--order",        eslARG_NONE,                              FALSE,  NULL,  NULL, TAXOPTS,    NULL,  NULL,               "taxonmy type to group sequences by",                                       0 },
  { "--family",       eslARG_NONE,                              FALSE,  NULL,  NULL, TAXOPTS,    NULL,  NULL,               "taxonmy type to group sequences by",                                       0 },  
  { "--genus",        eslARG_NONE,                              FALSE,  NULL,  NULL, TAXOPTS,    NULL,  NULL,               "taxonmy type to group sequences by",                                       0 },
  /* Control of output */
  { "-o",            eslARG_OUTFILE,                            FALSE,  NULL,   NULL,   NULL,    NULL,  NULL,               "send output to file <f>, not stdout",                                      0 },
  { "-w",            eslARG_NONE,                               FALSE,  NULL,   NULL,   NULL,    NULL,  NULL,               "TRUE to display plots",                                                    0 },
  { "--bintime",     eslARG_REAL,                               "1.0",  NULL,"x>0.0",   NULL,    NULL,  NULL,               "time interval for binning data",                                           0 },
  { "--xlabel",      eslARG_STRING,                      "time (MYA)",  NULL,    NULL,  NULL,    NULL,  NULL,               "xlabel",                                                                   0 },
  { "--ylabel",      eslARG_STRING,      "rate of change per residue",  NULL,    NULL,  NULL,    NULL,  NULL,               "ylabel",                                                                   0 },
  /* msa format */
  { "--informat",   eslARG_STRING,                               NULL,  NULL,    NULL,   NULL,   NULL,  NULL,               "specify format",                                                           0 },
  /* Selecting the alphabet rather than autoguessing it */
  { "--amino",      eslARG_NONE,                                 FALSE, NULL, NULL,  ALPHOPTS,   NULL,  NULL,               "input alignment is protein sequence data",                                 2 },
  { "--dna",        eslARG_NONE,                                 FALSE, NULL, NULL,  ALPHOPTS,   NULL,  NULL,               "input alignment is DNA sequence data",                                     2 },
  { "--rna",        eslARG_NONE,                                 FALSE, NULL, NULL,  ALPHOPTS,   NULL,  NULL,               "input alignment is RNA sequence data",                                     2 },
  /* Control of orthologs */
  { "--orthcutoff",  eslARG_REAL,                               "0.70",  NULL,  NULL,   NULL,    NULL,  NULL,               "min rio boostrap value for orthologs single linkage",                      0 },
  /* Control of homolgy */
  { "--evalcutoff",  eslARG_REAL,                             "0.0001",  NULL,  NULL,   NULL,    NULL,  NULL,               "evalue cutoff for homology searches",                                      0 },
  /* Control of scoring system - substitutions */
  { "--mx",         eslARG_STRING,                          "BLOSUM62",  NULL,  NULL,   NULL,    NULL,  "--mxfile",         "substitution rate matrix choice (of some built-in matrices)",              0 },
  { "--mxfile",     eslARG_INFILE,                                NULL,  NULL,  NULL,   NULL,    NULL,  "--mx",             "read substitution rate matrix from file <f>",                              0 },
  /* Control of scoring system - indels */ 
  { "--evomodel",   eslARG_STRING,   "AFGC",   NULL,       NULL,   NULL,    NULL,  EVOMOPTS,                                "evolutionary model used",                                                  0 },
  { "--paramfile",  eslARG_STRING,     NULL,   NULL,       NULL,   NULL,    NULL,  NULL,                                    "file with rate parameters (overrides the individual options below)",       0 },
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
  /* alternative options for plotting clusters of ortholgos from a given family */
  { "--serial",     eslARG_NONE,                                 FALSE,   NULL,  NULL, PLOTOPTS,  NULL,  NULL,               "plot clusters of ortholgs in serial figures",                             0 },
  { "--together",   eslARG_NONE,                                 FALSE,   NULL,  NULL, PLOTOPTS,  NULL,  NULL,               "plot clusters of ortholgs all in the same figure",                        0 },
  { "--joint",      eslARG_NONE,                                 FALSE,   NULL,  NULL, PLOTOPTS,  NULL,  NULL,               "plot clusters of ortholgs usign joint linear regressions",                0 },
  { "--plotall",    eslARG_NONE,                             "default",   NULL,  NULL, PLOTOPTS,  NULL,  NULL,               "plot all the above",                                                      0 },
  /* other options */
  { "--tol",        eslARG_REAL,                              "0.0001",   NULL,  NULL,   NULL,    NULL,  NULL,               "tolerance",                                                               0 },
  { "--seed",        eslARG_INT,                                   "0",   NULL, "n>=0",  NULL,    NULL,  NULL,               "set RNG seed to <n>",                                                     5 },
{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
};
static char usage[]  = "[-options] <pfamname> ";
static char banner[] = "calculate pairwise rates of evolution using E2";

static int e2Dickerson_onecluster(FILE *outfp, char *gnuplot, char *spfile, char *distfile, char *histofile, char *ratefile, char *binratefile,
				  ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, ORTHO *ortho, E2_PIPELINE *pli, E1_BG *bg, P7_BG *bg7, E2_ALI e2ali, 
				  float evalcutoff, E1_RATE *R, P7_RATE *R7, int mode, int do_viterbi, float bintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose);
static int remove_fragments(struct cfg_s *cfg, ESL_MSA **msa, int *ret_nfrags, int *ret_seq_cons_len);
static int select_subset(struct cfg_s *cfg, ESL_MSA **msa, int *ret_nremoved);
static int outfile_header(char *acc, int full, char *subsmx, float muA, float muE, float ldE, float muI, float ldI, float rI, char **ret_outheader);

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
  if (esl_opt_ArgNumber(go) != 1)    { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_IsOn(go, "--msafile")) {
    if ((cfg.msafile  = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  }
  else {
      if ((cfg.pfamname  = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  }
  
  cfg.r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));
  
  esl_sprintf(&cfg.gnuplot, "%s -persist", getenv("GNUPLOT"));

  /* If you know the MSA file format, set it (<infmt>, here). */
  cfg.infmt = eslMSAFILE_UNKNOWN;
  if (esl_opt_IsOn(go, "--informat") &&
      (cfg.infmt = eslx_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));

  if      (esl_opt_GetBoolean(go, "--amino"))   cfg.abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg.abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg.abc = esl_alphabet_Create(eslRNA);
  else                                          cfg.abc = NULL;

  if ( esl_opt_IsOn(go, "-o") ) {
    if ((cfg.outfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else cfg.outfp = stdout;

  /* other options */
  cfg.fragfrac   = esl_opt_GetReal(go, "-F");
  cfg.idthresh   = esl_opt_GetReal(go, "-I");
  cfg.tol        = esl_opt_GetReal(go, "--tol");
  cfg.verbose    = esl_opt_GetBoolean(go, "-v");
  cfg.view       = esl_opt_GetBoolean(go, "-w");
  cfg.makeT      = esl_opt_GetBoolean(go, "--makeT");
  cfg.full       = esl_opt_GetBoolean(go, "--full");
  cfg.group      = esl_opt_GetInteger(go, "--group");
  cfg.orthcutoff = esl_opt_GetReal(go, "--orthcutoff");
  cfg.evalcutoff = esl_opt_GetReal(go, "--evalcutoff");

  cfg.evomodel       = e1_rate_Evomodel(esl_opt_GetString(go, "--evomodel"));
  cfg.rateparam.muAM = esl_opt_GetReal  (go, "--muA");
  cfg.rateparam.muI  = esl_opt_GetReal  (go, "--muI");
  cfg.rateparam.ldI  = esl_opt_GetReal  (go, "--ldI");
  cfg.rateparam.muEM = esl_opt_GetReal  (go, "--muEM");
  cfg.rateparam.ldEM = esl_opt_GetReal  (go, "--ldEM");
  cfg.rateparam.muED = esl_opt_GetReal  (go, "--muED");
  cfg.rateparam.ldED = esl_opt_GetReal  (go, "--ldED");
  cfg.rateparam.rI   = esl_opt_GetReal  (go, "--rI");
  cfg.rateparam.rM   = esl_opt_GetReal  (go, "--rM");
  cfg.rateparam.rD   = esl_opt_GetReal  (go, "--rD");

  /* the paramfile takes precedent over individually feed values above */
  cfg.paramfile = NULL;
  if (esl_opt_IsOn(go, "--paramfile")) {
    cfg.paramfile = esl_opt_GetString(go, "--paramfile");
    status = e1_rate_ReadParamfile(cfg.paramfile, &cfg.rateparam, &cfg.evomodel, cfg.errbuf, cfg.verbose);
    if (status != eslOK) esl_fatal("Failed to read paramfile %s\n%s", cfg.paramfile, cfg.errbuf);
  }

  cfg.subsmx     = esl_opt_GetString(go, "--mx");
  cfg.mode       = esl_opt_IsOn(go, "--local")? e2_LOCAL : e2_GLOBAL;
  cfg.do_viterbi = FALSE;

  cfg.bintime    = esl_opt_GetReal  (go, "--bintime");
  cfg.tlinear    = 0.0;

  cfg.xlabel     = esl_opt_GetString(go, "--xlabel");
  cfg.ylabel     = esl_opt_GetString(go, "--ylabel");

  cfg.nboot      = esl_opt_GetInteger(go, "--nboot");
  cfg.reboot     = esl_opt_GetBoolean(go, "--reboot");
  cfg.clsize     = esl_opt_GetInteger(go, "--clsize");
  cfg.clmax      = esl_opt_GetBoolean(go, "--clmax");

  if (esl_opt_IsOn(go, "--submsa")) cfg.submsa = esl_opt_GetInteger(go, "--submsa");
  else                              cfg.submsa = 0;
 
  if      (esl_opt_GetBoolean(go, "--e2"))      cfg.e2ali = E2;
  else if (esl_opt_GetBoolean(go, "--e2f"))     cfg.e2ali = E2F;
  else if (esl_opt_GetString(go, "--e2hmmer"))  cfg.e2ali = E2HMMER;
  else if (esl_opt_GetString(go, "--e2fhmmer")) cfg.e2ali = E2FHMMER;
  else                                          cfg.e2ali = E2NONE;
  
  cfg.R1  = NULL;
  cfg.R7  = NULL;
  cfg.bg  = NULL;
  cfg.bg7 = NULL;
  if (cfg.e2ali == E2 || cfg.e2ali == E2F) cfg.bg  = e1_bg_Create(cfg.abc);
  else                                     cfg.bg7 = p7_bg_Create(cfg.abc);

  /* if e2hmmer or e2fhmmer, read the input HMM */
  if (cfg.e2ali == E2FHMMER || cfg.e2ali == E2HMMER)
    if ((cfg.hmmfile = esl_opt_GetString(go, "--e2hmmer")) == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  if      (esl_opt_GetBoolean(go, "--superkingdom"))   cfg.taxotype = SUPERKINGDOM;
  else if (esl_opt_GetBoolean(go, "--kingdom"))        cfg.taxotype = KINGDOM;
  else if (esl_opt_GetBoolean(go, "--phylum"))         cfg.taxotype = PHYLUM;
  else if (esl_opt_GetBoolean(go, "--class"))          cfg.taxotype = CLASS;
  else if (esl_opt_GetBoolean(go, "--order"))          cfg.taxotype = ORDER;
  else if (esl_opt_GetBoolean(go, "--family"))         cfg.taxotype = FAMILY;
  else if (esl_opt_GetBoolean(go, "--genus"))          cfg.taxotype = GENUS;
  else                                                 cfg.taxotype = UNKNOWN;
  
  if      (esl_opt_GetBoolean(go, "--domya"))          cfg.myatype  = DOMYA;
  else if (esl_opt_GetBoolean(go, "--fromfile"))       cfg.myatype  = FROMFILE;
  else if (esl_opt_GetBoolean(go, "--onthespot"))      cfg.myatype  = ONTHESPOT;

  if      (esl_opt_GetBoolean(go, "--serial"))         cfg.plottype = SERIAL;
  else if (esl_opt_GetBoolean(go, "--together"))       cfg.plottype = TOGETHER;
  else if (esl_opt_GetBoolean(go, "--joint"))          cfg.plottype = JOINT;
  else if (esl_opt_GetBoolean(go, "--plotall"))        cfg.plottype = PLOTALL;
 
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

static int
output_header(FILE *ofp, ESL_GETOPTS *go)
{
  p7_banner(ofp, go->argv[0], banner);
  
  if (esl_opt_IsUsed(go, "-o")          && fprintf(ofp, "# output directed to file:        %s\n",             esl_opt_GetString(go, "-o"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--mx")        && fprintf(ofp, "# subst rate matrix (built-in):   %s\n",             esl_opt_GetString(go, "--mx"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--mxfile")    && fprintf(ofp, "# subst rate matrix (file):       %s\n",             esl_opt_GetString(go, "--mxfile"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--amino")     && fprintf(ofp, "# input alignment is asserted as:   protein\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--dna")       && fprintf(ofp, "# input alignment is asserted as:   DNA\n")                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--rna")       && fprintf(ofp, "# input alignment is asserted as:   RNA\n")                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;	          /* command line processing  */
  struct cfg_s     cfg;
  char            *msg = "e2_Dickerson() failed";
  char            *famtreefile = NULL;
  char            *famspfile = NULL;
  char            *famsptreefile = NULL;
  char            *famdistfile = NULL;
  char            *histofile = NULL;
  char            *ratefile = NULL;
  char            *binratefile = NULL;
  char            *histopdf = NULL;
  char            *ratepdf = NULL;
  char            *binratepdf = NULL;
  FILE            *treefp = NULL; 
  FILE            *spfp = NULL; 
  ESLX_MSAFILE    *msafp = NULL;
  ESL_MSA         *msa = NULL;
  ORTHO           *cortho = NULL;         /* clusters of ortholgos */
  ORTHO           *ortho;                 /* pointer to a particular orthlog cluster */
  ESL_TREE        *T = NULL;              /* gene tree */

  P7_HMMFILE      *hfp = NULL;            /* input HMM file */
  P7_HMM          *hmm = NULL;            /* HMM            */
  EMRATE          *emR = NULL;
  double           avgid;
  int              seq_cons_len;
  int              nfrags = 0;	  	  /* # of fragments removed */
  int              nremoved = 0;	  /* # of identical sequences removed */
  int              ncl;                   /* total number of ortholog clusters per family to analyze */
  int              nc;
  int              c;                     /* index for orthlogs clusters */
  int              cmax;
  int              cidx;
  int              r;
  int              hstatus = eslOK;   
  int              status = eslOK;   
  
  /* Initializations */
  process_commandline(argc, argv, &go, &cfg);    

  if (!cfg.msafile) { /* take the .msa, .tree and .species files from the pfamdb */
    if (fetch_MSATreeAndSpecies(&cfg.msafile, &famtreefile, &famspfile, &famsptreefile, &famdistfile, cfg.pfamname, 
				cfg.full, cfg.group, cfg.taxotype, cfg.myatype, cfg.errbuf) != eslOK) 
      esl_fatal("Failed to access pfamdb. %s", cfg.errbuf);
  }
  
  /* Open the MSA file */
  status = eslx_msafile_Open(&cfg.abc, cfg.msafile, NULL, cfg.infmt, NULL, &msafp);
  if (status != eslOK) eslx_msafile_OpenFailure(msafp, status);
  
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
    betainf = (cfg.rateparam.muEM > 0.)? ((cfg.rateparam.ldEM < cfg.rateparam.muEM)? cfg.rateparam.ldEM/cfg.rateparam.muEM : 1.0 ) : -1.;
    etainf  = (cfg.rateparam.muI  > 0.)? ((cfg.rateparam.ldI  < cfg.rateparam.muI)?  cfg.rateparam.ldI/cfg.rateparam.muI : 1.0 ) : -1.;
    if (cfg.evomodel != GG) cfg.rateparam.rI = -1.0; /* does not have to be provided */
    if (cfg.evomodel != GG) etainf           = -1.0; /* does not have to be provided */
    if (p7_RateCalculate(stdout, hmm, cfg.bg7, emR, NULL, &(cfg.R7), cfg.evomodel, betainf, -1.0, 0.001, cfg.errbuf, FALSE) != eslOK)  
      esl_fatal("%s", cfg.errbuf);
  }

  /* Initializations */
  cfg.pli = e2_pipeline_Create(go, (hmm)? hmm->M:1, 100, 100);
  if (cfg.e2ali == E2 || cfg.e2ali == E2F) cfg.bg  = e1_bg_Create(cfg.abc);
  else                                     cfg.bg7 = p7_bg_Create(cfg.abc);
  
  output_header(cfg.outfp, go);   
 
  while ((hstatus = eslx_msafile_Read(msafp, &msa)) == eslOK) 
    {
      esl_msa_ConvertDegen2X(msa); 
      esl_msa_Hash(msa);

      if (remove_fragments(&cfg, &msa, &nfrags, &seq_cons_len) != eslOK) { printf("remove_fragments failed\n"); esl_fatal(msg); }
      select_subset(&cfg, &msa, &nremoved);
      esl_dst_XAverageId(cfg.abc, msa->ax, msa->nseq, 10000, &avgid); /* 10000 is max_comparisons, before sampling kicks in */

      /* print some info */
      fprintf(cfg.outfp, "name                  %%id     alen   cons_len nseq    nfrag    nid     used         msa\n");
      fprintf(cfg.outfp, "%-20s  %3.0f%% %6d %6d   %6d %6d   %6d   ", msa->name, 100.*avgid, (int)msa->alen, seq_cons_len, msa->nseq, nfrags, nremoved);
      
      /* select submsa (RIO has space limitations) */
      if (cfg.submsa) {
	if (msamanip_SelectSubset(cfg.r, cfg.submsa, &msa, &cfg.msafile, cfg.errbuf, cfg.verbose) != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
      }
      fprintf(cfg.outfp, "%6d          %s\n", msa->nseq, cfg.msafile);
      if (cfg.verbose) { if (eslx_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal(msg); }
      /* outheader for all output files. PFxxxxxx.{seed,full} || PFxxxxxx.{seed,full}.random<n> */
      outfile_header(msa->acc, cfg.full, cfg.subsmx, cfg.rateparam.muAM, cfg.rateparam.muEM, cfg.rateparam.ldEM, cfg.rateparam.muI, cfg.rateparam.ldI, cfg.rateparam.rI, &cfg.outheader); 
      
       /* The genes tree T 
	*/
      if (cfg.makeT) { /* If requested, calculate the tree using FastTree */
	if (Tree_CalculateExtFromMSA(msa, &T, TRUE, cfg.errbuf, cfg.verbose) != eslOK) { printf("%s\n", cfg.errbuf); esl_fatal(msg); }
      }
      else { /* read the pfam tree */
	if ((treefp = fopen(famtreefile, "r")) == NULL) { printf("failed to open %s\n", famtreefile); esl_fatal(msg); }
	if (esl_tree_ReadNewick(treefp, cfg.errbuf, &T) != eslOK) esl_fatal("Failed to read tree %s: %s", famtreefile, cfg.errbuf);
	fclose(treefp); treefp = NULL;
      }
      if (cfg.verbose) Tree_Dump(cfg.outfp, T, "original Tree");

      /* find orthologs */
      if ((status = ortho_FindClusters(&famsptreefile, cfg.outfp, cfg.r, msa, cfg.outheader, &cortho, &nc, cfg.nboot, cfg.reboot,
				       cfg.orthcutoff, cfg.errbuf, cfg.verbose)) != eslOK) 
	{ printf("%s\n", cfg.errbuf); esl_fatal(msg); } 
      if (cortho == NULL) { printf("Failed to produce orthologs cluster\n"); esl_fatal(msg); } 
      
      /* files for writing the plots, joinly for all clusters of orthologs */
      esl_sprintf(&histofile,   "%s.mya.histo",    cfg.outheader);
      esl_sprintf(&ratefile,    "%s.rate.plot",    cfg.outheader);
      esl_sprintf(&binratefile, "%s.binrate.plot", cfg.outheader);
      remove(histofile);
      remove(ratefile);
      remove(binratefile);

      if (cfg.clmax) { /* for the maximal cluster run the whole method */
	ncl = 1;  
	cmax = 0;
	for (c = 0; c < nc; c ++) { if (cortho[c].no > cmax) { cmax = cortho[c].no; cidx = c; } }
	ortho = &(cortho[cidx]);
	fprintf(cfg.outfp, "\nAnalysing maximal cluster with %d sequences\n", ortho->no);
	if (e2Dickerson_onecluster(cfg.outfp, cfg.gnuplot, famspfile, famdistfile, histofile, ratefile, binratefile, cfg.r, msa, T, ortho, cfg.pli, cfg.bg, cfg.bg7, cfg.e2ali, cfg.evalcutoff, 
				   cfg.R1[0], cfg.R7, cfg.mode, cfg.do_viterbi, cfg.bintime, cfg.tlinear, cfg.full, cfg.tol, cfg.view, cfg.errbuf, cfg.verbose)
	    != eslOK) esl_fatal("%s", cfg.errbuf);
      }
      else {      /* for each cluster with more at least 'clsize' species, run the whole method */
	ncl = 0;
	for (c = 0; c < nc; c ++) {
	  ortho =  &(cortho[c]);
	  if (ortho->no >= cfg.clsize) {
	    fprintf(cfg.outfp, "\nAnalysing cluster[%d] with %d sequences\n", ncl, ortho->no);
	    if (e2Dickerson_onecluster(cfg.outfp, cfg.gnuplot, famspfile, famdistfile, histofile, ratefile, binratefile, cfg.r, msa, T, ortho, cfg.pli, cfg.bg, cfg.bg7, cfg.e2ali, cfg.evalcutoff, 
				       cfg.R1[0], cfg.R7, cfg.mode, cfg.do_viterbi, cfg.bintime, cfg.tlinear, cfg.full, cfg.tol, cfg.view, cfg.errbuf, cfg.verbose)
		!= eslOK) esl_fatal("%s", cfg.errbuf);
	    ncl ++;
	  }
	}
      }
      if (ncl == 0) { printf("found no group of ortholgs with at least %d sequences to analyze\n", cfg.clsize); return status; }

      /* plot results */
      esl_sprintf(&histopdf,   "%s.pdf", histofile);
      remove(histopdf);
      
      if (plot_gplot_Histogram(cfg.gnuplot, histopdf, 1, &histofile, "MYA", cfg.errbuf) == eslOK) { 
	if (cfg.view) plot_file_Display(histopdf);
      } else esl_fatal("%s\nbad histogram", cfg.errbuf);

      switch(cfg.plottype) {
      case SERIAL:
	esl_sprintf(&ratepdf,    "%s.serial.pdf", ratefile);
	esl_sprintf(&binratepdf, "%s.serial.pdf", binratefile);
	remove(ratepdf);
	remove(binratepdf);
	if (plot_gplot_SerialRatesWithLinearRegression(cfg.gnuplot, ratepdf, 1, &ratefile, cfg.evalcutoff, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(ratepdf);
	} else esl_fatal("%s\nbad rates", cfg.errbuf);
	if (plot_gplot_SerialBinnedRatesWithLinearRegression(cfg.gnuplot, binratepdf, 1, &binratefile, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(binratepdf);
	} else esl_fatal("%s\nbad binned rates", cfg.errbuf);
	break;
      case TOGETHER:
	esl_sprintf(&ratepdf,    "%s.together.pdf", ratefile);
	esl_sprintf(&binratepdf, "%s.together.pdf", binratefile);
	remove(ratepdf);
	remove(binratepdf);
	if (plot_gplot_TogetherRatesWithLinearRegression(cfg.gnuplot, ratepdf, 1, &ratefile, cfg.evalcutoff, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(ratepdf);
	} else esl_fatal("%s\nbad rates", cfg.errbuf);
	if (plot_gplot_TogetherBinnedRatesWithLinearRegression(cfg.gnuplot, binratepdf, 1, &binratefile, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(binratepdf);
	} else esl_fatal("%s\nbad binned rates", cfg.errbuf);
	break;
      case JOINT:
	esl_sprintf(&ratepdf,    "%s.joint.pdf", ratefile);
	esl_sprintf(&binratepdf, "%s.joint.pdf", binratefile);
	remove(ratepdf);
	remove(binratepdf);
	if (plot_gplot_JointRatesWithLinearRegression(cfg.gnuplot, ratepdf, 1, &ratefile, cfg.evalcutoff, cfg.tlinear, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(ratepdf);
	} else esl_fatal("%s\nbad rates", cfg.errbuf);
	if (plot_gplot_JointBinnedRatesWithLinearRegression(cfg.gnuplot, binratepdf, 1, &ratefile, cfg.bintime, cfg.tlinear, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(binratepdf);
	} else esl_fatal("%s\nbad binned rates", cfg.errbuf);
	break;
      case PLOTALL:
	esl_sprintf(&ratepdf,    "%s.serial.pdf", ratefile);
	esl_sprintf(&binratepdf, "%s.serial.pdf", binratefile);
	remove(ratepdf);
	remove(binratepdf);
	if (plot_gplot_SerialRatesWithLinearRegression(cfg.gnuplot, ratepdf, 1, &ratefile, cfg.evalcutoff, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(ratepdf);
	} else esl_fatal("%s\nbad rates", cfg.errbuf);
	if (plot_gplot_SerialBinnedRatesWithLinearRegression(cfg.gnuplot, binratepdf, 1, &binratefile, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(binratepdf);
	} else esl_fatal("%s\nbad binned rates", cfg.errbuf);
	free(ratepdf);    ratepdf = NULL;
	free(binratepdf); binratepdf = NULL;

	esl_sprintf(&ratepdf,    "%s.together.pdf", ratefile);
	esl_sprintf(&binratepdf, "%s.together.pdf", binratefile);
	remove(ratepdf);
	remove(binratepdf);
	if (plot_gplot_TogetherRatesWithLinearRegression(cfg.gnuplot, ratepdf, 1, &ratefile, cfg.evalcutoff, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(ratepdf);
	} else esl_fatal("%s\nbad rates", cfg.errbuf);
	if (plot_gplot_TogetherBinnedRatesWithLinearRegression(cfg.gnuplot, binratepdf, 1, &binratefile, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(binratepdf);
	} else esl_fatal("%s\nbad binned rates", cfg.errbuf);
	free(ratepdf);    ratepdf = NULL;
	free(binratepdf); binratepdf = NULL;

	esl_sprintf(&ratepdf,    "%s.joint.pdf", ratefile);
	esl_sprintf(&binratepdf, "%s.joint.pdf", binratefile);
	remove(ratepdf);
	remove(binratepdf);
	if (plot_gplot_JointRatesWithLinearRegression(cfg.gnuplot, ratepdf, 1, &ratefile, cfg.evalcutoff, cfg.tlinear, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(ratepdf);
	} else esl_fatal("%s\nbad rates", cfg.errbuf);
	if (plot_gplot_JointBinnedRatesWithLinearRegression(cfg.gnuplot, binratepdf, 1, &ratefile, cfg.bintime, cfg.tlinear, cfg.xlabel, cfg.ylabel, cfg.errbuf) == eslOK) { 
	  if (cfg.view) plot_file_Display(binratepdf);
	} else esl_fatal("%s\nbad binned rates", cfg.errbuf);
	break;
       default: printf("unknown plotting type"); return eslFAIL;
      }
  
      /* cleanup */
      free(histofile); histofile = NULL;
      free(ratefile); ratefile = NULL;
      free(binratefile); binratefile = NULL;
      free(histopdf); histopdf = NULL;
      if (ratepdf)free(ratepdf); ratepdf = NULL;
      if (binratepdf) free(binratepdf); binratepdf = NULL;
      esl_msa_Destroy(msa); msa = NULL;
      ortho_Destroy(nc, cortho); cortho = NULL;
      esl_tree_Destroy(T); T = NULL;
    }
  if (hstatus != eslEOF) eslx_msafile_ReadFailure(msafp, hstatus);

  fclose(cfg.outfp);
  if (treefp) fclose(treefp);
  if (spfp) fclose(spfp);      
  eslx_msafile_Close(msafp); 

  /* cleanup */
  if (cfg.R1) {
    for (r = 0; r < cfg.nr; r ++) e1_rate_Destroy(cfg.R1[r]);
    free(cfg.R1);
  }
  if (cfg.R7) p7_RateDestroy(cfg.R7);
  free(cfg.msafile);
  free(famtreefile);	
  free(famspfile); 
  free(famsptreefile); 
  free(famdistfile); 
  free(cfg.gnuplot);
  esl_randomness_Destroy(cfg.r);
  if (msa) esl_msa_Destroy(msa);
  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(cfg.abc);
  if (T) esl_tree_Destroy(T);
  e2_pipeline_Destroy(cfg.pli);
  if (cfg.bg)  e1_bg_Destroy(cfg.bg);     
  if (cfg.bg7) p7_bg_Destroy(cfg.bg7);     
  if (hmm) p7_hmm_Destroy(hmm); 
  if (hfp) p7_hmmfile_Close(hfp);

  return status;
}

static int
e2Dickerson_onecluster(FILE *outfp, char *gnuplot, char *spfile, char *distfile, char *histofile, char *ratefile, char *binratefile,
		       ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_TREE *T, ORTHO *ortho, E2_PIPELINE *pli, E1_BG *bg, P7_BG *bg7, E2_ALI e2ali,  
		       float evalcutoff, E1_RATE *R, P7_RATE *R7, int mode, int do_viterbi, float bintime, float tlinear, int full, float tol, int view, char *errbuf, int verbose)
{
  ESL_MSA         *omsa = NULL;
  ESL_MSA         *msa2 = NULL;           /* msa for pair of orthogonal sequences*/
  float           *msa2frq = NULL;
  ESL_TREE        *oT = NULL;             /* orthogonal gene tree */
  ESL_TREE        *T2 = NULL;             /* gene tree for pair of orthogonal sequences*/
  SPECIES         *SPE = NULL;            /* species */
  SPECIES         *SPE2 = NULL;           /* species */
  int             *usepair = NULL;
  int              no = 0;
  int              i;
  int              n, m;
  int              status;
 
  if ((status = esl_msa_SequenceSubset(msa, ortho->useme, &omsa)) != eslOK) { printf("msa subset failed\n"); status = eslFAIL; goto ERROR; }

  /* allocate */
  ESL_ALLOC(usepair, sizeof(int) * omsa->nseq);

#if 0
  /* collapse T to the orhtolog sequences */
  oT = Tree_Collapse(T, msa, ortho->useme, errbuf, verbose);
  if (oT == NULL) esl_fatal(msg); 
  if ((status = Tree_ReorderTaxaAccordingMSA(omsa, oT, errbuf, verbose)) != eslOK) esl_fatal(msg);
  if ((status = Tree_RootAtMidPoint(&oT, NULL, errbuf, verbose))         != eslOK) esl_fatal(msg);
#else
  /* make the T for the orthogonal sequences with FastTree */
  if (Tree_CalculateExtFromMSA(omsa, &oT, TRUE, errbuf, verbose) != eslOK) { printf("%s\n", errbuf); status = eslFAIL; goto ERROR; }
#endif
  
  if (verbose) Tree_Dump(outfp, oT, "orthologs Tree");
  
  if (verbose) {
    for (i = 0; i < ortho->ntaxa; i ++) {
      if (ortho->useme[i]) {
	fprintf(outfp, "oseq[%d]: %s\n", no++, ortho->taxalist[i]);
      }
    }
  }
  
  /* From the species file, find species and taxonomy for each sequence in the alignment */
  if (fetch_ReadSpeciesFile(outfp, spfile, errbuf, &SPE, omsa, TRUE) != eslOK) 
    esl_fatal("Failed to read species file %s: %s", spfile, errbuf);
  esl_msa_Hash(omsa);

  /* do the actual calculations */  
  for (n = 0; n < omsa->nseq-1; n ++) {
    for (m = n+1; m < omsa->nseq; m ++) {  

      esl_vec_ISet(usepair, omsa->nseq, FALSE);
      usepair[n] = TRUE;
      usepair[m] = TRUE;

      if ((status = esl_msa_SequenceSubset(omsa, usepair, &msa2)) != eslOK) { printf("msa subset failed\n"); status = eslFAIL; goto ERROR; }
      msamanip_XBaseComp(msa2, bg->f, &msa2frq);

      /* From the species file, find species and taxonomy for each sequence in the alignment */
      if (fetch_ReadSpeciesFile(outfp, spfile, errbuf, &SPE2, msa2, FALSE) != eslOK) 
	esl_fatal("Failed to read species file %s: %s", spfile, errbuf);

      if (Tree_CalculateExtFromMSA(msa2, &T2, TRUE, errbuf, verbose) != eslOK) { printf("%s\n", errbuf); status = eslFAIL; goto ERROR; }
      if (verbose) {
	Tree_Dump(stdout, T2, "T2");
      }
      
      if (e2BranchInference(outfp, gnuplot, distfile, histofile, ratefile, binratefile, 
			    r, msa2, msa2frq, T2, SPE2, pli, bg, bg7, e2ali, evalcutoff, R, R7, mode, do_viterbi, bintime, tlinear, full, tol, view, errbuf, verbose) != eslOK) 
	{ printf("%s\n", errbuf); status = eslFAIL; goto ERROR; }
      
      esl_msa_Destroy(msa2); msa2 = NULL;
      esl_tree_Destroy(T2); T2 = NULL;
      species_Destroy(SPE2); SPE2 = NULL;
    }
  }
  
  free(usepair);
  esl_msa_Destroy(omsa);
  free(msa2frq);
  esl_tree_Destroy(oT); 
  species_Destroy(SPE); 
  return eslOK;
  
 ERROR:
  if (usepair) free(usepair);
  if (omsa)    esl_msa_Destroy(omsa);
  if (msa2)    esl_msa_Destroy(msa2);
  if (msa2frq) free(msa2frq);
  if (oT)      esl_tree_Destroy(oT); 
  if (T2)      esl_tree_Destroy(T2); 
  if (SPE)     species_Destroy(SPE); 
  if (SPE2)    species_Destroy(SPE2); 
  return status;
 }



/* Step 1. Label all sequence fragments < fragfrac of average raw length 
   msa is changed in place
*/
static int
remove_fragments(struct cfg_s *cfg, ESL_MSA **msa, int *ret_nfrags, int *ret_seq_cons_len)
{
  ESL_MSA *omsa;
  ESL_MSA *new = NULL;
  ESL_DSQ *dsq = NULL;
  int     *useme = NULL;
  double   flen;
  int64_t  clen;           /* seq_cons length */
  int64_t  alen;
  int      len;
  int      n;
  int      i;
  int      x;
  int      nfrags;
  int      status;

  omsa = *msa; /* pointer to original msa */

  for (n = 0; n < omsa->ngc; n++) {
    if (strcmp(omsa->gc_tag[n], "seq_cons") == 0) {
      esl_abc_CreateDsq(omsa->abc, omsa->gc[n], &dsq);
      clen = esl_abc_dsqrlen(omsa->abc, dsq);
      alen = esl_abc_dsqlen(dsq);
      printf("%s: len %d alen %d\n%s\n\n", omsa->gc_tag[n], (int)clen, (int)alen, omsa->gc[n]);
    }
  }
  if (clen == 0.) { printf("couldn't find 'seq_cons'\n"); status = eslFAIL; goto ERROR; } 
  flen = (double)clen * cfg->fragfrac;

  /* for each seq in msa, calculate number of residues that match the seq_cons */
  ESL_ALLOC(useme, sizeof(int) * omsa->nseq);
  for (nfrags = 0, i = 0; i < omsa->nseq; i++)  {
    len = 0;
    for (x = 1; dsq[x] != eslDSQ_SENTINEL; x++) { 
      if (esl_abc_XIsResidue(omsa->abc, dsq[x]) && esl_abc_XIsResidue(omsa->abc, omsa->ax[i][x])) len ++;
    }
    useme[i] = (len < flen) ? 0 : 1;
  }
  if ((status = esl_msa_SequenceSubset(omsa, useme, &new)) != eslOK) goto ERROR;

  *ret_seq_cons_len = clen;
  *ret_nfrags = omsa->nseq - esl_vec_ISum(useme, omsa->nseq);

 /* replace msa */
  esl_msa_Destroy(omsa);
  *msa = new;

  free(dsq);
  free(useme);
  return eslOK;

 ERROR:
  if (new) esl_msa_Destroy(new);
  if (dsq   != NULL) free(dsq);
  if (useme != NULL) free(useme); 
  return status;
}

/* Step 2. Extract subset with non-identical sequences
 */
static int
select_subset(struct cfg_s *cfg, ESL_MSA **msa, int *ret_nremoved)
{      
  ESL_MSA   *omsa  = NULL;
  ESL_MSA   *new  = NULL;
  int       *assignment = NULL;
  int       *nin        = NULL;
  int       *useme      = NULL;
  int        nused      = 0;
  int        nc         = 0;
  int        c;
  int        nskip;
  int        i;
  int        status;

  omsa = *msa;

  ESL_ALLOC(useme, sizeof(int) * omsa->nseq);
  esl_vec_ISet(useme, omsa->nseq, 0);

  if ((status = esl_msacluster_SingleLinkage(omsa, cfg->idthresh, &assignment, &nin, &nc)) != eslOK) goto ERROR;
 
  for (c = 0; c < nc; c++)
    {
      nskip = esl_rnd_Roll(cfg->r, nin[c]); /* pick a random seq in this cluster to be the test. */
      for (i = 0; i < omsa->nseq; i++)
	if (assignment[i] == c) {
	  if (nskip == 0) {
	    nused ++;
	    useme[i] = 1;
	    break;
	  } else nskip--; 
	}
    }
  *ret_nremoved = omsa->nseq - nused;

  if ((status = esl_msa_SequenceSubset(omsa, useme, &new)) != eslOK) goto ERROR;


 /* replace msa */
  esl_msa_Destroy(omsa);
  *msa = new;

  free(useme);
  free(nin);
  free(assignment);

 return eslOK;

 ERROR:
  if (useme      != NULL) free(useme);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  if (new        != NULL) esl_msa_Destroy(new);
 return status;
}

static int
outfile_header(char *acc, int full, char *subsmx, float muA, float muE, float ldE, float muI, float ldI, float rI, char **ret_outheader)
{      
  char *outheader = NULL;
  char *st;
  char *name;
  char *nametail;

  st = acc; /* remove the version from the accession, but keep the possible 'random<n>' tail present in msa->acc if option --msasub has been used */
  esl_strtok(&st, ".", &name);
  esl_strtok(&st, ".", &nametail);
  if (full) (strcmp(st, "\0")==0)? esl_sprintf(&outheader, "%s.full", name) : 
	      esl_sprintf(&outheader, "%s.full.%s", name, nametail);
  else      (strcmp(st, "\0")==0)? esl_sprintf(&outheader, "%s.seed", name) : 
	      esl_sprintf(&outheader, "%s.seed.%s", name, nametail);
 
  *ret_outheader = outheader;

  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/
