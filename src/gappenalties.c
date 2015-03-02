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


static int AFRmodel(ESL_GETOPTS  *go, FILE *fp, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double gapscale, double gapo, double gape, EVOM evomodel, double betainf, double rM, int scaledrate, int N, 
		    double tol, char *errbuf, int verbose);
static int AFR_calculate_insrate(FILE *fp, double *ret_rI, double *ret_muA, double *ret_ld, double gaposc, double gapesc, double betainf, double rM, double tol, char *errbuf, int verbose);
static int AFR_calculate_gapcosts(E1_RATE *R1, double time, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double gapscale, double *ret_gape, double *ret_gapo, double tol, char *errbuf, int verbose);
static int e1model_calculate_subrate(ESL_GETOPTS *go, ESL_ALPHABET *abc, P7_BG *bg, RATEBUILDER **ret_ratebld, int scaled, double tol, char *errbuf, int verbose);

static ESL_OPTIONS options[] = {
  /* name           type              default  env  range        toggles     reqs        incomp             help                                                                  docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL, NULL,       NULL,       NULL,       NULL,              "show brief help on version and usage",                                               0 },
  { "--gapo",       eslARG_REAL,      "-11.0", NULL, "x<0",      NULL,       NULL,       NULL,              "gap open",                                                                           0 },
  { "--gape",       eslARG_REAL,       "-1.0", NULL, "x<0",      NULL,       NULL,       NULL,              "gap exted",                                                                          0 },
  { "--gapsc",      eslARG_REAL,       "2.0",  NULL, "x>0",      NULL,       NULL,       NULL,              "gap scale",                                                                          0 },
  { "--betainf",    eslARG_REAL,     "0.499",  NULL, "x>=0.",    NULL,       NULL,       NULL,              "betainf = beta at time infinity (if ld<muA)",                                        0 },
  { "--rM",         eslARG_REAL,       "0.00", NULL, "x>=0",     NULL,       NULL,       NULL,              "fragment parameter rM",                                                              0 },
  { "-N",           eslARG_INT,        "100",  NULL, "n>0",      NULL,       NULL,       NULL,              "number of time points",                                                              0 },
  { "-v",           eslARG_NONE,        FALSE, NULL, NULL,       NULL,       NULL,       NULL,              "be verbose",                                                                         0 },
  
  /* Control of scoring system */
  { "--mx",       eslARG_STRING, "BLOSUM62",   NULL, NULL,       NULL,       NULL, "--mxfile",              "substitution score matrix choice (of some built-in matrices)",                       0 },
  { "--mxfile",   eslARG_INFILE,       NULL,   NULL, NULL,       NULL,       NULL,     "--mx",              "read substitution score matrix from file <f>",                                       0 },
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
  double        gaposc;
  double        gapesc;
  double        gapscale;
  EVOM          evomodel;
  double        betainf;
  double        rM;
  double        tol = 0.01;
  int           scaledrate = FALSE;
  int           mode = JOINT;
  int           N;
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
  gapscale = esl_opt_GetReal(go, "--gapsc");
  gaposc   = gapo/gapscale;
  gapesc   = gape/gapscale;
  evomodel = AFR;
  N        = esl_opt_GetInteger(go, "-N");
  betainf  = esl_opt_GetReal(go, "--betainf");
  rM       = esl_opt_GetReal(go, "--rM");
  verbose  = esl_opt_GetBoolean(go, "-v");
  if (betainf >= 0.5)       { printf("betainf should be smaller than 1/2\n"); exit(1); }
  if (rM < 0.0 || rM > 1.0) { printf("rm should be 0 < rM < 1\n");            exit(1); }
  
  fprintf(stdout, "# gap open:               = %g\n", gapo);
  fprintf(stdout, "# gap extend:             = %g\n", gape);
  fprintf(stdout, "# rM:                     = %g\n", rM);
  fprintf(stdout, "# beta at time infinity:  = %g\n", betainf);
  
  fprintf(fp, "# gap open:               = %g\n", gapo);
  fprintf(fp, "# gap extend:             = %g\n", gape);
  fprintf(fp, "# rM:                     = %g\n", rM);
  fprintf(fp, "# beta at time infinity:  = %g\n", betainf);
  
  AFRmodel(go, fp, abc, bg, mode, L, gapscale, gaposc, gapesc, evomodel, betainf, rM, scaledrate, N, tol, errbuf, verbose);
  
  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  fclose(fp);
  return 0;
}


static int
AFRmodel(ESL_GETOPTS  *go, FILE *fp, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double gapscale, double gapo, double gape, EVOM evomodel, 
	 double betainf, double rM, int scaledrate, int N, double tol, char *errbuf, int verbose)
{
  char               *msg = "e1model unit test failed";
  RATEBUILDER        *ratebld = NULL;           /* construction configuration */
  ESL_DMATRIX        *P = NULL;
  E1_RATE            *R1   = NULL; 
  struct rateparam_s  rateparam;
  double              rX;
  double              ld;
  double              muA;
  double              ttz = 1e-10;
  double              tt;
  double              gapet, gapot;
  double              subsite;
  double              expsc;
  double              pid;
  double              totalt = 50.0;
  int                 x;
  int                 status;
  
  status = e1model_calculate_subrate(go, abc, bg, &ratebld, scaledrate, tol, errbuf, FALSE);
  if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }

  status = AFR_calculate_insrate(fp, &rX, &muA, &ld, gapo, gape, betainf, rM, tol, errbuf, FALSE);
  if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }

  /* pass the values to structure rateparam */
  rateparam.rI   = rX;
  rateparam.rD   = rX;
  rateparam.rM   = rM;
  rateparam.ldI  = ld;
  rateparam.muAM = muA;

  R1 = e1_rate_CreateWithValues(abc, evomodel, rateparam, NULL, ratebld->Q, gapscale, tol, errbuf, FALSE);
  if (R1 == NULL) { printf("%s\n", errbuf); esl_fatal(msg); }

  for (x = 0; x <= totalt*N; x ++) {
    tt = (x<=N)? ttz + x*(1.0-ttz)/N : ttz + x*(1.0-ttz)/N;

    P = ratematrix_ConditionalsFromRate(tt, ratebld->Q, tol, errbuf, verbose);
    subsite = ratematrix_DFreqSubsPerSite(P, ratebld->p);
    pid = (1.0 - subsite) * 100.;
    expsc = 2.0 * ratematrix_ExpScore(P, ratebld->p) / eslCONST_LOG2; // in half bits

    status = AFR_calculate_gapcosts(R1, tt, abc, bg, mode, L, gapscale, &gapet, &gapot, tol, errbuf, FALSE);
    if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
    
    if (verbose) fprintf(stdout, "%f %f %f %f %f %f\n", tt, pid, 100*subsite, expsc, gapet, gapot);
    fprintf(fp,     "%f %f %f %f %f %f\n", tt, pid, 100*subsite, expsc, gapet, gapot);

    esl_dmatrix_Destroy(P); P = NULL;
  }
  /* at infinity */
  tt = 1e+8;
  status = AFR_calculate_gapcosts(R1, tt, abc, bg, mode, L, gapscale, &gapet, &gapot, tol, errbuf, FALSE); 
  if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  
  e1_rate_Destroy(R1);
  ratebuilder_Destroy(ratebld);
  return eslOK;
}

static int
AFR_calculate_insrate(FILE *fp, double *ret_rI, double *ret_muA, double *ret_ld, double gaposc, double gapesc, double betainf, double rM, 
		      double tol, char *errbuf, int verbose)
{
  double Tstar_MX;
  double Tstar_XX;
  double betastar;
  double newbetastar;
  double muA;
  double ld;
  double rI;
  int    status;

  /* the insertion/deletion rates */
  Tstar_XX = exp(gapesc * eslCONST_LOG2);
  Tstar_MX = exp(gaposc * eslCONST_LOG2) *  Tstar_XX / (1.0 - Tstar_XX);
  betastar = Tstar_MX/(1.0 - rM);

  muA = log(betainf) + log(1.0 - betastar) - log(betainf - betastar);
  ld  = muA * betainf / (1.0 - betainf);
  rI  = (Tstar_XX - betastar) / (1.0 - betastar);
  
  newbetastar = betainf * (1.0 - exp(-muA)) / (1.0 - betainf*exp(-muA));
  if (muA < 0.0) { status = eslFAIL; goto ERROR; }
  if (ld  < 0.0){  status = eslFAIL; goto ERROR; }
  if (rI >= 1.0 || rI <= 0.0) { status = eslFAIL; goto ERROR; }

  fprintf(fp, "# gap open/extend: %f %f | TstarXX %f TstarMX %f betastar %.10f (%f) new %f\n", 
	  gaposc, gapesc, Tstar_XX, Tstar_MX, betastar, log(betastar), newbetastar);
  fprintf(fp, "# rI  = %f\n", rI);
  fprintf(fp, "# ld  = %f\n", ld);
  fprintf(fp, "# muA = %f\n", muA);

  fprintf(stdout, "# gap open/extend: %f %f | TstarXX %f TstarMX %f betastar %.10f (%f) new %f\n", 
	  gaposc, gapesc, Tstar_XX, Tstar_MX, betastar, log(betastar), newbetastar);
  fprintf(stdout, "# rI  = %f\n", rI);
  fprintf(stdout, "# ld  = %f\n", ld);
  fprintf(stdout, "# muA = %f\n", muA);

  *ret_rI   = rI;
  *ret_muA  = muA;
  *ret_ld   = ld;
  
  return eslOK;

 ERROR:
  return status;
}

static int
AFR_calculate_gapcosts(E1_RATE *R1, double time, ESL_ALPHABET *abc, P7_BG *bg, int mode, int L, double gapscale, double *ret_gapet, double *ret_gapot,
		       double tol, char *errbuf, int verbose)
{
  E1_MODEL *evom = NULL;
  double    TMX;
  double    TXX;
  double    TXM;
  double    TXY;
  double    TXE;
  double    gapet;    
  double    gapot;   
  int       status;

  evom = e1_model_Create(R1, time, NULL, bg->f, mode, L, abc, tol, errbuf, verbose);
  if (evom == NULL) { status = eslFAIL; goto ERROR; }

  if (verbose) e1_model_Dump(stdout, evom); 
  
  TMX = evom->t[e1H_SI];
  TXX = evom->t[e1H_II];
  TXM = evom->t[e1H_IS];
  TXY = evom->t[e1H_ID];
  TXE = evom->t[e1H_IE];

  printf("time %f MX-XM %f XY %f MX-XE %f MX(1-TXX) %f\n", time, log(TMX) + log(TXM), log(TXY), log(TMX) + log(TXE), log(TMX) + log(1-TXX));
  gapet = log(TXX);
  gapot = log(TMX) + log(1.0 - TXX) - log(TXX);
 
  gapet /= eslCONST_LOG2;
  gapot /= eslCONST_LOG2;
		     
  gapet *= gapscale;
  gapot *= gapscale;

  e1_model_Destroy(evom); evom = NULL;
  
  *ret_gapet = gapet;
  *ret_gapot = gapot;

  return eslOK;

 ERROR:
  if (evom)  e1_model_Destroy(evom);
  return status;
}



static int
e1model_calculate_subrate(ESL_GETOPTS *go, ESL_ALPHABET *abc, P7_BG *bg, RATEBUILDER **ret_ratebld, int scaledrate, double tol, char *errbuf, int verbose)
{
  RATEBUILDER  *ratebld = NULL;           /* construction configuration   */
  int           status;

  /* Initialize a default ratebuilder configuration,
   * then set only the options we need for single sequence search
   */
  ratebld = ratebuilder_Create(abc);

   /* Default is stored in the --mx option, so it's always IsOn(). Check --mxfile first; then go to the --mx option and the default. */
  if (esl_opt_IsOn(go, "--mxfile")) status = ratebuilder_SetScoreSystem (ratebld, esl_opt_GetString(go, "--mxfile"), NULL, bg);
  else                              status = ratebuilder_LoadScoreSystem(ratebld, esl_opt_GetString(go, "--mx"),           bg, scaledrate); 
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



/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/
