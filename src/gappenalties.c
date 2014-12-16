/* gappenalties -- calculate rate matrices from a substitution scoring matrix
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_rootfinder.h"
#include "esl_scorematrix.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2.h"
#include "e1_rate.h"
#include "e1_model.h"
#include "evohmmer.h"
#include "p7_evopipeline.h"
#include "ratematrix.h"
#include "ratebuilder.h"


static int e1model(ESL_GETOPTS  *go, FILE *fp, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double gapo, double gape, EVOM evomodel, double betainf, double etainf, double rI, int N, 
		   double tol, char *errbuf, int verbose);
static int e1model_calculate_gapcosts(E1_RATE *R1, double time, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double *ret_gape, double *ret_gapo, double *ret_gapO, double *ret_logbeta, 
				      double tol, char *errbuf, int verbose);
static int e1model_calculate_subrate(ESL_GETOPTS *go, ESL_ALPHABET *abc, P7_BG *bg, RATEBUILDER **ret_ratebld, double tol, char *errbuf, int verbose);
static int e1model_calculate_insrate(FILE *fp, double *ret_rI, double *ret_muE, double *ret_muI, double *ret_ld, double *ret_ldI, double gapo, double gape, EVOM evomodel, double betainf, 
				     double etainf, double rI, double tol, char *errbuf, int verbose);

static ESL_OPTIONS options[] = {
  /* name           type              default  env  range       toggles     reqs        incomp             help                                                                  docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL, NULL,      NULL,       NULL,       NULL,              "show brief help on version and usage",                                               0 },
  { "--gapo",       eslARG_REAL,      "-11.0", NULL, "x<=0",    NULL,       NULL,       NULL,              "gap open",                                                                           0 },
  { "--gape",       eslARG_REAL,      "-1.0",  NULL, "x<=0",    NULL,       NULL,       NULL,              "gap exted",                                                                          0 },
  { "--evomodel",   eslARG_STRING,      "AIF", NULL, NULL,      NULL,       NULL,       NULL,             "evolutionary model used",                                      0 },
  { "--betainf",    eslARG_REAL,       "0.69", NULL, "x>=0",    NULL,       NULL,       NULL,              "betainf = ldf/muf = beta at time infinity (if ldf<muf)",       0 },
  { "--etainf",     eslARG_REAL,       "0.69", NULL, "x>=0",    NULL,       NULL,       NULL,              "etainf = ldI/muI = eta at time infinity (if ldI<muI)",         0 },
  { "--rI",       eslARG_REAL,       "0.40", NULL, "x>=0",    NULL,       NULL,       NULL,              "rI = eta at time zero",                                      0 },
  { "-v",           eslARG_NONE,        FALSE, NULL, NULL,      NULL,       NULL,       NULL,              "be verbose",                                                                         0 },

/* Control of scoring system */
  { "--mx",       eslARG_STRING, "BLOSUM62",   NULL, NULL,      NULL,       NULL, "--mxfile",              "substitution score matrix choice (of some built-in matrices)",                       0 },
  { "--mxfile",   eslARG_INFILE,       NULL,   NULL, NULL,      NULL,       NULL,     "--mx",              "read substitution score matrix from file <f>",                                       0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <evomfile>";
static char banner[] = "evolve a selection of substitution matrix and gap penalties";


int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go  = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char          errbuf[eslERRBUFSIZE];
  char         *evomfile;
  FILE         *fp = NULL;                /* output file handle  */
  ESL_ALPHABET *abc = NULL;               /* digital alphabet    */
  P7_BG        *bg  = NULL;		  /* null model (copies made of this into threads)    */
  double        gapo;
  double        gape;
  EVOM          evomodel;
  double        betainf;
  double        etainf;
  double        rI;
  double        tol = 0.01;
  int           mode = e2_GLOBAL;
  int           N = 1000;
  int           L = 1e+20;
  int           verbose;

  if (esl_opt_ArgNumber(go) != 1)                  { puts("Incorrect number of command line arguments");          exit(1); }
  if ((evomfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <evohmmfile> argument on command line"); exit(1); }

 /* Initializations */
  abc = esl_alphabet_Create(eslAMINO);
  bg  = p7_bg_Create(abc);

 /* Open evom output file */
  fp = fopen(evomfile, "w");
  if (fp == NULL) p7_Fail("Failed to open file %s for writing", evomfile);

  /* Options */
  gapo     = esl_opt_GetReal(go, "--gapo");
  gape     = esl_opt_GetReal(go, "--gape");
  evomodel = e1_rate_Evomodel(esl_opt_GetString(go, "--evomodel"));
  betainf  = esl_opt_GetReal(go, "--betainf");
  etainf   = esl_opt_GetReal(go, "--etainf");
  rI     = esl_opt_GetReal(go, "--rI");
  verbose  = esl_opt_GetBoolean(go, "-v");
  fprintf(stdout, "# gap open:               = %g\n", gapo);
  fprintf(stdout, "# gap extend:             = %g\n", gape);
  fprintf(stdout, "# beta at time infinity:  = %g\n", betainf);
  
  fprintf(fp, "# gap open:               = %g\n", gapo);
  fprintf(fp, "# gap extend:             = %g\n", gape);
  fprintf(fp, "# beta at time infinity:  = %g\n", betainf);
 
  e1model(go, fp, abc, bg, mode, L, gapo, gape, evomodel, betainf, etainf, rI, N, tol, errbuf, verbose);

  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  fclose(fp);
  return 0;
}


static int
e1model(ESL_GETOPTS  *go, FILE *fp, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double gapo, double gape, EVOM evomodel, double betainf, double etainf, double rI, int N, double tol, char *errbuf, int verbose)
{
  char               *msg = "e1model unit test failed";
  RATEBUILDER        *ratebld = NULL;           /* construction configuration                       */
  E1_RATE            *R1   = NULL; 
  struct rateparam_s  rateparam;
  double              ttz = 1e-10;
  double              tt;
  double              gapet, gapot, gapOt;
  double              logbetat;
  double              subsite;
  int                 x;
  int                 status;

  status = e1model_calculate_subrate(go, abc, bg, &ratebld, tol, errbuf, FALSE);
  if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  subsite = ratematrix_SubsPerSite(ratebld->Q, ratebld->p);

  rateparam.rI = rI;
  rateparam.muAM  = 0.0; /*  bogus, not needed */
  status = e1model_calculate_insrate(fp, &rateparam.rI, &rateparam.muEM, &rateparam.muI, &rateparam.ldEM, &rateparam.ldI, gapo, gape, evomodel, betainf, etainf, rI, tol, errbuf, FALSE);
  if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }

  R1 = e1_rate_CreateWithValues(abc, evomodel, rateparam, NULL, ratebld->Q, tol, errbuf, FALSE);
  if (R1 == NULL) { printf("%s\n", errbuf); esl_fatal(msg); }

  for (x = 0; x <= 2*N; x ++) {
    tt = (x<=N)? ttz + x*(1.0-ttz)/N : ttz + x*(1.0-ttz)/N;

    status = e1model_calculate_gapcosts(R1, tt, abc, bg, mode, L, &gapet, &gapot, &gapOt, &logbetat, tol, errbuf, FALSE);
    if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }

    if (verbose) fprintf(stdout, "%f %f %f %f %f %f\n", tt, subsite*tt, gapet, logbetat, gapot, gapOt);
    fprintf(fp,     "%f %f %f %f %f %f\n", tt, subsite*tt, gapet, logbetat, gapot, gapOt);
  }
  /* at infinity */
  tt = 1e+8;
  status = e1model_calculate_gapcosts(R1, tt, abc, bg, mode, L, &gapet, &gapot, &gapOt, &logbetat, tol, errbuf, FALSE); 
  if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  fprintf(stdout, "# etainf = exp(%f) = %f should be (%f) betainf = exp(%f) = %f \n", gapet, exp(gapet), (betainf<=1.)?log(betainf):1.0, logbetat, exp(logbetat));
  fprintf(fp,     "# etainf = exp(%f) = %f should be (%f) betainf = exp(%f) = %f \n", gapet, exp(gapet), (betainf<=1.)?log(betainf):1.0, logbetat, exp(logbetat));
   
  e1_rate_Destroy(R1);
  ratebuilder_Destroy(ratebld);
  return eslOK;
}


static int
e1model_calculate_gapcosts(E1_RATE *R1, double time, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double *ret_gape, double *ret_gapo, double *ret_gapO, double *ret_logbeta, 
			   double tol, char *errbuf, int verbose)
{
  E1_MODEL *evom = NULL;
  double    gape;    /* log eta_t                                                  */               
  double    gapo;    /* log beta_t*(1-eta_t)/eta_t ||  gapcost = gapo + (n)  *gape */
  double    gapO;    /* log beta_t*(1-eta_t)       ||  gapcost = gapO + (n-1)*gape */
  double    logbeta;
  int       status;

  evom = e1_model_Create(R1, time, NULL, bg->f, mode, L, abc, tol, errbuf, verbose);
  if (evom == NULL) { status = eslFAIL; goto ERROR; }

  if (verbose) e1_model_Dump(stdout, evom); 
  
  gape    = log(evom->t[e1H_II]);
  logbeta = log(evom->t[e1H_SI]);
  gapO    = logbeta + log(1. - evom->t[e1H_II]);
  gapo    = gapO - gape;
  
  e1_model_Destroy(evom); evom = NULL;
  
  *ret_gape    = gape;
  *ret_gapo    = gapo;
  *ret_gapO    = gapO;
  *ret_logbeta = logbeta;

  return eslOK;

 ERROR:
  if (evom)  e1_model_Destroy(evom);
  return status;
}

static int
e1model_calculate_subrate(ESL_GETOPTS *go, ESL_ALPHABET *abc, P7_BG *bg, RATEBUILDER **ret_ratebld, double tol, char *errbuf, int verbose)
{
  RATEBUILDER  *ratebld = NULL;           /* construction configuration   */
  int           status;

  /* Initialize a default ratebuilder configuration,
   * then set only the options we need for single sequence search
   */
  ratebld = ratebuilder_Create(abc);

   /* Default is stored in the --mx option, so it's always IsOn(). Check --mxfile first; then go to the --mx option and the default. */
  if (esl_opt_IsOn(go, "--mxfile")) status = ratebuilder_SetScoreSystem (ratebld, esl_opt_GetString(go, "--mxfile"), NULL, bg);
  else                              status = ratebuilder_LoadScoreSystem(ratebld, esl_opt_GetString(go, "--mx"),           bg); 
  if (status != eslOK) p7_Fail("Failed to set single query seq score system:\n%s\n", ratebld->errbuf);
  
  status = ratematrix_CreateFromConditionals(ratebld->P, ratebld->p, &(ratebld->Q), &(ratebld->E), tol, errbuf, verbose);
  if (status != eslOK) p7_Fail("Failed to calculate rate matrix\n%s\n", errbuf);
  if (verbose) {
    esl_scorematrix_Write(stdout, ratebld->S);
    esl_dmatrix_Dump(stdout, ratebld->P, NULL, NULL);
    esl_dmatrix_Dump(stdout, ratebld->Q, NULL, NULL);
  }
  
  *ret_ratebld = ratebld;
  return eslOK;
 }

static int
e1model_calculate_insrate(FILE *fp, double *ret_rI, double *ret_muE, double *ret_muI, double *ret_ld, double *ret_ldI, double gapo, double gape, EVOM evomodel, double betainf, double etainf, double rI,
			  double tol, char *errbuf, int verbose)
{
  struct rateparam_s rateparam;
  double eta_star;
  double beta_star;
  double muE;
  double muI;
  double ld;
  double ldI;
  int    status;

  /* the insertion/deletion rates */
  eta_star  = exp(gape);
  beta_star = exp(gapo + gape - log(1. - eta_star));

  status = e1_rate_CalculateInsertRates(evomodel, &rateparam, betainf, betainf, eta_star, beta_star, 0.0, tol, errbuf, verbose);
  if (status != eslOK) goto ERROR;
  
  fprintf(stdout, "# gap open/extend: %f %f | eta %.10f (%f) beta %.10f (%f)\n", gapo, gape, eta_star, log(eta_star), beta_star, log(beta_star));
  fprintf(stdout, "# rI = %g\n", rI);
  fprintf(stdout, "# ldf  = %g\n", ld);
  fprintf(stdout, "# muf  = %g\n", muE);
  fprintf(stdout, "# ldI  = %g\n", ldI);
  fprintf(stdout, "# muI  = %g\n", muI);
  
  fprintf(fp, "# gap open/extend: %f %f | eta %.10f (%f) beta %.10f (%f)\n", gapo, gape, eta_star, log(eta_star), beta_star, log(beta_star));
  fprintf(fp, "# rI = %g\n", rI);
  fprintf(fp, "# ldf  = %g\n", ld);
  fprintf(fp, "# muf  = %g\n", muE);
  fprintf(fp, "# ldI = %g\n", ldI);
  fprintf(fp, "# muI = %g\n", muI);

  *ret_rI = rI;
  *ret_muE  = muE;
  *ret_muI  = muI;
  *ret_ld   = ld;
  *ret_ldI  = ldI;
  
  return eslOK;

 ERROR:
  return status;
 }




/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/
