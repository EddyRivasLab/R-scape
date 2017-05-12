/* pottsim -- 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_stats.h"
#include "esl_vectorops.h"

#include "pottsbuild.h"
#include "pottsscore.h"


#define ALPHOPTS     "--amino,--dna,--rna"                           /* Exclusive options for alphabet choice */

static ESL_OPTIONS options[] = {
  /* name           type              default  env       range      toggles     reqs        incomp        help                                                    docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "show brief help on version and usage",                        1 },
  { "-v",           eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "be verbose",                                                  1 },
  /* parameters to control the simulation */
  { "-L",            eslARG_INT,       "100",  NULL,     "n>0",     NULL,       NULL,       NULL,         "length of ancestral sequence",                                1 },
  { "-N",            eslARG_INT,     "10000",  NULL,     "n>0",     NULL,       NULL,       NULL,         "number of sequences for a given set of parameters and time",  1 },
  /* alphabet type */
  { "--dna",          eslARG_NONE,     FALSE,  NULL,      NULL, ALPHOPTS,       NULL,       NULL,         "use DNA alphabet",                                            1 },
  { "--rna",          eslARG_NONE,    "TRUE",  NULL,      NULL, ALPHOPTS,       NULL,       NULL,         "use RNA alphabet",                                            1 },
  { "--amino",        eslARG_NONE,     FALSE,  NULL,      NULL, ALPHOPTS,       NULL,       NULL,         "use protein alphabet",                                        1 },
  /* parameters of the potts model */
  { "--pottsfile",  eslARG_INFILE,    NULL,    NULL,      NULL,     NULL,       NULL,       NULL,         "read potts params from file <f>",                             1 },
  { "--sigma",        eslARG_REAL,   "0.1",    NULL,    "x>=0",     NULL,       NULL,       NULL,         "if sampling param from a N(0,sigma)    ",                     1 },
  /* other options */  
  { "--tol",        eslARG_REAL,     "0.0001", NULL,      NULL,     NULL,       NULL,       NULL,         "tolerance",                                                   1 },
  { "--seed",        eslARG_INT,          "0", NULL,    "n>=0",     NULL,       NULL,       NULL,         "set RNG seed to <n>",                                         1 },
   {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msaoutfile>";
static char banner[] = "pottsim";

static int pottsim_banner(FILE *fp, char *progname, char *banner);
static int pottsim(ESL_GETOPTS *go, ESL_RANDOMNESS *r, FILE *fp, int N, PT *pt, double tol, char *errbuf, int verbose);

int
main(int argc, char **argv)
{
  ESL_GETOPTS        *go  = esl_getopts_Create(options);
  ESL_RANDOMNESS     *r = NULL;
  ESL_ALPHABET       *abc = NULL;
  char                errbuf[eslERRBUFSIZE];
  char               *pottsfile = NULL;
  char               *msaoutfile;
  FILE               *msafp;
  PT                 *pt = NULL;
  double              tol;
  double              sigma;
  int                 N;
  int                 L;
  int                 verbose;
  int                 status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* R-scape banner */
  pottsim_banner(stdout, argv[0], banner);

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_usage(stdout,  argv[0], usage);
      if (puts("\noptions:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
     exit(0);
    }

  msaoutfile = NULL;  
  if (esl_opt_ArgNumber(go) != 1)                    { puts("Incorrect number of command line arguments");          exit(1); }
  if ((msaoutfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <msaoutfile> argument on command line"); exit(1); }

  /* Initializations */
  r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  
 /* Open msa output file */
  msafp = fopen(msaoutfile, "w");
  if (msafp == NULL) esl_fail("Failed to open file %s for writing", msaoutfile);

  /* Options */
  L         = esl_opt_GetInteger(go, "-L");
  N         = esl_opt_GetInteger(go, "-N");
  verbose   = esl_opt_GetBoolean(go, "-v");
  sigma     = esl_opt_GetReal(go, "--sigma");
  tol       = esl_opt_GetReal(go, "--tol");

   /* alphabet */
  abc = NULL;
  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);   
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);  
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO);
  
  if ( esl_opt_IsOn(go, "--pottsfile") ) {
    pottsfile = esl_opt_GetString(go, "--pottsfile");
    if (!esl_FileExists(pottsfile))  esl_fatal("pottsfile %s does not seem to exist\n", pottsfile);    
    pt = potts_Read(pottsfile, abc, errbuf);
    if (pt == NULL) esl_fatal("error creating pt structure from %s\n", pottsfile);  
  }
  else {
    pt = potts_Create(L, abc->K, abc, 0.0, NONE, SCNONE);
    potts_InitGaussian(r, pt, 0.0, sigma);
  }

  if (1||verbose) {
    fprintf(stdout, "# L         = %d\n", L);
    fprintf(stdout, "# N         = %d\n", N);
    fprintf(stdout, "# K         = %d\n", abc->K);
    potts_Write(stdout, pt);
  }
  
  pottsim(go, r, msafp, N, pt, tol, errbuf, verbose);

  potts_Destroy(pt);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  fclose(msafp);
  return 0;

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
pottsim_banner(FILE *fp, char *progname, char *banner)
{
  char *appname = NULL;
  int   status;

  if ((status = esl_FileTail(progname, FALSE, &appname)) != eslOK) return status;

  if (fprintf(fp, "# %s :: %s\n", appname, banner)                                               < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");

  if (appname) free(appname);
  return eslOK;

 ERROR:
  if (appname) free(appname);
  return status;
}

static int
pottsim(ESL_GETOPTS  *go, ESL_RANDOMNESS *r, FILE *fp, int N, PT *pt, double tol, char *errbuf, int verbose)
{
  char           *msg = "pottsim failed";
  ESL_MSA        *msa = NULL;
  ESL_HISTOGRAM   *h = NULL;
  int              L = pt->L;
  int              n;
  int              status;

  msa = esl_msa_Create(N, L);
  
  //esl_msafile_Write(fp, msa, eslMSAFILE_STOCKHOLM);

  esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
   if (msa) esl_msa_Destroy(msa);
   return status;
}




/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/
