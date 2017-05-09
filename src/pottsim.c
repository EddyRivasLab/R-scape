/* pottsim -- 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_stats.h"
#include "esl_vectorops.h"

#include "pottsbuild.h"
#include "pottsscore.h"


#define POTTSOPTS    "LI,LR,AFG,AFGR,AFR,AIF,GG,AG,AGA"              /* Exclusive options for evolutionary model choice */
#define ALPHOPTS     "--amino,--dna,--rna"                           /* Exclusive options for alphabet choice */

static ESL_OPTIONS options[] = {
  /* name           type              default  env       range      toggles     reqs        incomp        help                                                    docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "show brief help on version and usage",                        0 },
  { "-v",           eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "be verbose",                                                  0 },
  /* parameters to control the simulation */
  { "-L",            eslARG_INT,      "1000",  NULL,     "n>0",     NULL,       NULL,       NULL,         "length of ancestral sequence",                                0 },
  { "-N",            eslARG_INT,     "10000",  NULL,     "n>0",     NULL,       NULL,       NULL,         "number of sequences for a given set of parameters and time",  0 },
  /* alphabet type */
  { "--dna",          eslARG_NONE,     FALSE,  NULL,      NULL,     NULL,       NULL,       NULL,         "use DNA alphabet",                                            0 },
  { "--rna",          eslARG_NONE,    "TRUE",  NULL,      NULL,     NULL,       NULL,       NULL,         "use RNA alphabet",                                            0 },
  { "--amino",        eslARG_NONE,     FALSE,  NULL,      NULL,     NULL,       NULL,       NULL,         "use protein alphabet",                                        0 },
  /* parameters of the potts model */
  { "--evomodel",   eslARG_STRING,    "AIF",   NULL,      NULL,     NULL,       NULL,       NULL,         "evolutionary model used",                                     0 },
  { "--mu",         eslARG_REAL,    "0.0004",  NULL,    "x>=0",     NULL,       NULL,       NULL,         "rate of deletion of ancestral sequences",                     0 },
  /* other options */  
  { "--tol",        eslARG_REAL,     "0.0001", NULL,      NULL,     NULL,       NULL,       NULL,         "tolerance",                                                   0 },
  { "--seed",        eslARG_INT,          "0", NULL,    "n>=0",     NULL,       NULL,       NULL,         "set RNG seed to <n>",                                         5 },
   {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <outfile>";
static char banner[] = "pottsim";


static int pottsim(ESL_GETOPTS *go, ESL_RANDOMNESS *r, FILE *fp, FILE *fitfp, int N, PT *pt,
		   double tol, char *errbuf, int verbose);

int
main(int argc, char **argv)
{
  ESL_GETOPTS        *go  = esl_getopts_Create(options);
  ESL_RANDOMNESS     *r = NULL;
  ESL_ALPHABET       *abc = NULL;
  char                errbuf[eslERRBUFSIZE];
  char               *outfile;
  char               *fitfile;
  FILE               *fp = NULL;                /* output file handle  */
  FILE               *fitfp = NULL;             /* goodnes of fit file handle  */
  PT                 *pt = NULL;
  double              tol;
  int                 N;
  int                 L;
  int                 verbose;

  if (esl_opt_ArgNumber(go) != 1)                 { puts("Incorrect number of command line arguments");          exit(1); }
  if ((outfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <evohmmfile> argument on command line"); exit(1); }

  /* Initializations */
  r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  
 /* Open evom output file */
  fp = fopen(outfile, "w");
  if (fp == NULL) esl_fail("Failed to open file %s for writing", outfile);

  /* Options */
  L         = esl_opt_GetInteger(go, "-L");
  N         = esl_opt_GetInteger(go, "-N");
  verbose   = esl_opt_GetBoolean(go, "-v");
  tol       = esl_opt_GetReal(go, "--tol");

   /* alphabet */
  abc = NULL;
  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);   
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);  
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO);
  
  pt = potts_Create(L, abc, 0.0, NONE);
 
  /* geometric goodness of fit */
  if (esl_opt_IsUsed(go, "--geofit")) {
    fitfile = esl_opt_GetString(go, "--geofit");
    fitfp = fopen(fitfile, "w");
    if (fitfp == NULL) esl_fail("Failed to open file %s for writing", fitfile);
  }
  else printf("Option --geofit:  (not set)\n");
  
  fprintf(stdout, "# L         = %d\n", L);
  fprintf(stdout, "# N         = %d\n", N);
 
  fprintf(fp, "# L         = %d\n", L);
  fprintf(fp, "# N         = %d\n", N);
  
  pottsim(go, r, fp, fitfp, N, pt, tol, errbuf, verbose);

  potts_Destroy(pt);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  fclose(fp);
  if (fitfp) fclose(fitfp);
  return 0;
}


static int
pottsim(ESL_GETOPTS  *go, ESL_RANDOMNESS *r, FILE *fp, FILE *fitfp, int N, PT *pt,
      double tol, char *errbuf, int verbose)
{
  char           *msg = "pottsim failed";
  ESL_MSA        *msa = NULL;
  ESL_HISTOGRAM   *h = NULL;
  int              L = pt->L;
  int              n;
  int              status;


 
  esl_msa_Destroy(msa);
  if (h) esl_histogram_Destroy(h);
   return eslOK;

 ERROR:
   if (msa) esl_msa_Destroy(msa);
  if (h) esl_histogram_Destroy(h);
   return status;
}



/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/
