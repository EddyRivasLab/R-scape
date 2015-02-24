/* e1sim -- evolve according to the evolutionary model.
 *          both using the finite-form model, and the rates at infinitesimal time to compare.
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

#include "hmmer.h"

#include "e2.h"
#include "e1_simulate.h"
#include "e1_goodness.h"

#define EVOMOPTS "LI,LR,AFG,AFGR,AFR,AIF,GG,AG,AGA"              /* Exclusive options for evolutionary model choice */

static ESL_OPTIONS options[] = {
  /* name           type              default  env       range      toggles     reqs        incomp        help                                                    docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "show brief help on version and usage",                        0 },
  { "--geofit",   eslARG_STRING,         NULL, NULL,      NULL,     NULL,       NULL,       NULL,         "geometric goodness of fit",                                   0 },
  /* parameters to control the simulation */
  { "-L",            eslARG_INT,      "1000",  NULL,     "n>0",     NULL,       NULL,       NULL,         "length of ancestral sequence",                                0 },
  { "-N",            eslARG_INT,     "10000",  NULL,     "n>0",     NULL,       NULL,       NULL,         "number of sequences for a given set of parameters and time",  0 },
  { "--joint",      eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "generate ancestral sqs with geom L/(L+1)",                    0 },
  { "--rev",        eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "impose the condition L = ldI/(muA-ldI)  ",                    0 },
  { "--tmin",       eslARG_REAL,        "0.0", NULL,    "x>=0",     NULL,       NULL,       NULL,         "min time to evolve sequence",                                 0 },
  { "--tmax",       eslARG_REAL,        "3.0", NULL,     "x>0",     NULL,       NULL,       NULL,         "max time to evolve sequence",                                 0 },
  { "--tinterval",  eslARG_REAL,       "0.20", NULL,     "x>0",     NULL,       NULL,       NULL,         "time interval to perform a given infinitesimal calculation",  0 },
  { "--tinc",       eslARG_REAL,       "0.01", NULL,     "x>0",     NULL,       NULL,       NULL,         "time increments to calculate statistics",                     0 },
  { "--teps",       eslARG_REAL,       "1e-6", NULL,     "x>0",     NULL,       NULL,       NULL,         "small time increments for the infinitesimal simulation",      0 },
  /* parameters of the evolutionary model */
  { "--evomodel",   eslARG_STRING,    "AIF",   NULL,       NULL,   NULL,    NULL,  NULL,                 "evolutionary model used",                                      0 },
  { "--muAM",       eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of deletion of ancestral sequences",                     0 },
  { "--muAD",       eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of deletion of ancestral sequences",                     0 },
  { "--muAI",       eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of deletion of ancestral sequences",                     0 },
  { "--muEM",       eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of deletion of inserted residues",                       0 },
  { "--ldEM",       eslARG_REAL,    "0.0010",  NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of adding a res to an insert",                           0 },
  { "--muED",       eslARG_REAL,    "0.0004",  NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of deletion of inserted residues",                       0 },
  { "--ldED",       eslARG_REAL,    "0.0010",  NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of adding a res to an insert",                           0 },
  { "--muI",        eslARG_REAL,    "0.000",   NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of deletion of inserted sequences",                      0 },
  { "--ldI",        eslARG_REAL,    "0.000",   NULL,     "x>=0",   NULL,    NULL,  NULL,                  "rate of adding a whole insert",                               0 },
  { "--rI",         eslARG_REAL,     "0.89",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                  "fragment parameter for I",                                    0 },
  { "--rM",         eslARG_REAL,     "0.90",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                  "fragment parameter for M",                                    0 },
  { "--rD",         eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                  "fragment parameter for D",                                    0 },
  { "--sI",         eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                  "mic: adding an insert",                                       0 },
  { "--sD",         eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                  "mic: removing an insert",                                     0 },
  { "--vI",         eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                  "mic: adding to an existing insert",                           0 },
  { "--vD",         eslARG_REAL,     "0.55",   NULL,"0<=x<=1.0",   NULL,    NULL,  NULL,                  "mic: removing from an existing an insert",                    0 },
  { "-v",           eslARG_NONE,        FALSE, NULL,      NULL,    NULL,    NULL,  NULL,                  "be verbose",                                                  0 },
  /* other options */
  { "--tol",        eslARG_REAL,     "0.0001", NULL,      NULL,    NULL,    NULL,  NULL,                  "tolerance",                                                   0 },
  { "--seed",        eslARG_INT,          "0", NULL,    "n>=0",    NULL,    NULL,  NULL,                  "set RNG seed to <n>",                                         5 },
   {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <outfile>";
static char banner[] = "e1sim";


static int e1sim(ESL_GETOPTS *go, ESL_RANDOMNESS *r, FILE *fp, FILE *fitfp, int L, int N, 
		 double tmin, double tmax, double tinterval, double tinc, double teps, 
		 EVOM evomodel, struct rateparam_s  rateparam, int joint, int rev,
		 double tol, char *errbuf, int verbose);
static int  generate_ancestral(ESL_RANDOMNESS *r, E1_RATE *R,  ESQ ***ret_esq, int *L, int N, int joint, int rev, char *errbuf, int verbose); 

int
main(int argc, char **argv)
{
  ESL_GETOPTS        *go  = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS     *r = NULL;	  
  char                errbuf[eslERRBUFSIZE];
  char               *outfile;
  char               *fitfile;
  FILE               *fp = NULL;                /* output file handle  */
  FILE               *fitfp = NULL;             /* goodnes of fit file handle  */
  EVOM                evomodel;
  struct rateparam_s  rateparam;
  double              tmin;
  double              tmax;
  double              tinterval;
  double              tinc;
  double              teps;      /* tespsilon for infinitesimal simulations */
  double              tol;
  int                 N;
  int                 L;
  int                 joint;
  int                 rev;
  int                 verbose;

  if (esl_opt_ArgNumber(go) != 1)                 { puts("Incorrect number of command line arguments");          exit(1); }
  if ((outfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <evohmmfile> argument on command line"); exit(1); }

  /* Initializations */
  r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  
 /* Open evom output file */
  fp = fopen(outfile, "w");
  if (fp == NULL) p7_Fail("Failed to open file %s for writing", outfile);

  /* Options */
  L         = esl_opt_GetInteger(go, "-L");
  N         = esl_opt_GetInteger(go, "-N");
  joint     = esl_opt_GetBoolean(go, "--joint");
  rev       = esl_opt_GetBoolean(go, "--rev");
  tmin      = esl_opt_GetReal(go, "--tmin");
  tmax      = esl_opt_GetReal(go, "--tmax");
  tinterval = esl_opt_GetReal(go, "--tinterval");
  tinc      = esl_opt_GetReal(go, "--tinc");
  teps      = esl_opt_GetReal(go, "--teps");
  verbose   = esl_opt_GetBoolean(go, "-v");
  tol       = esl_opt_GetReal(go, "--tol");
 
  evomodel       = e1_rate_Evomodel(esl_opt_GetString(go, "--evomodel"));
  rateparam.muAM = esl_opt_GetReal  (go, "--muAM");
  rateparam.muAD = esl_opt_GetReal  (go, "--muAD");
  rateparam.muAI = esl_opt_GetReal  (go, "--muAI");
  rateparam.muI  = esl_opt_GetReal  (go, "--muI");
  rateparam.ldI  = esl_opt_GetReal  (go, "--ldI");
  rateparam.muEM = esl_opt_GetReal  (go, "--muEM");
  rateparam.ldEM = esl_opt_GetReal  (go, "--ldEM");
  rateparam.muED = esl_opt_GetReal  (go, "--muED");
  rateparam.ldED = esl_opt_GetReal  (go, "--ldED");
  rateparam.rI   = esl_opt_GetReal  (go, "--rI");
  rateparam.rM   = esl_opt_GetReal  (go, "--rM");
  rateparam.rD   = esl_opt_GetReal  (go, "--rD");
  rateparam.sI   = esl_opt_GetReal  (go, "--sI");
  rateparam.sD   = esl_opt_GetReal  (go, "--sD");
  rateparam.vI   = esl_opt_GetReal  (go, "--vI");
  rateparam.vD   = esl_opt_GetReal  (go, "--vD");

 
  /* geometric goodness of fit */
  if (esl_opt_IsUsed(go, "--geofit")) {
    fitfile = esl_opt_GetString(go, "--geofit");
    fitfp = fopen(fitfile, "w");
    if (fitfp == NULL) p7_Fail("Failed to open file %s for writing", fitfile);
  }
  else printf("Option --geofit:  (not set)\n");
  
  if (!joint) fprintf(stdout, "# L         = %d\n", L);
  else        fprintf(stdout, "# p         = %f\n", (double)L/((double)L+1.));
  fprintf(stdout, "# N         = %d\n", N);
  fprintf(stdout, "# tmin      = %g\n", tmin);
  fprintf(stdout, "# tmax      = %g\n", tmax);
  fprintf(stdout, "# tinterval = %g\n", tinterval);
  fprintf(stdout, "# tinc      = %g\n", tinc);
  fprintf(stdout, "# teps      = %g\n", teps);
  fprintf(stdout, "# \n");
  fprintf(stdout, "# rI:  = %g\n", rateparam.rI);
  fprintf(stdout, "# rM:  = %g\n", rateparam.rM);
  fprintf(stdout, "# rD:  = %g\n", rateparam.rD);
  fprintf(stdout, "# sI:  = %g\n", rateparam.sI);
  fprintf(stdout, "# sD:  = %g\n", rateparam.sD);
  fprintf(stdout, "# vI:  = %g\n", rateparam.vI);
  fprintf(stdout, "# vD:  = %g\n", rateparam.vD);
  fprintf(stdout, "# muAM = %g\n", rateparam.muAM);
  fprintf(stdout, "# muAD = %g\n", rateparam.muAD);
  fprintf(stdout, "# muAI = %g\n", rateparam.muAI);
  fprintf(stdout, "# ldI  = %g\n", rateparam.ldI);
  fprintf(stdout, "# muI  = %g\n", rateparam.muI);
  fprintf(stdout, "# ldEM = %g\n", rateparam.ldEM);
  fprintf(stdout, "# muEM = %g\n", rateparam.muEM);
  fprintf(stdout, "# ldED = %g\n", rateparam.ldED);
  fprintf(stdout, "# muED = %g\n", rateparam.muED);

  if (!joint) fprintf(fp, "# L         = %d\n", L);
  else        fprintf(fp, "# p         = %f\n", (double)L/((double)L+1.));
  fprintf(fp, "# N         = %d\n", N);
  fprintf(fp, "# tmin      = %g\n", tmin);
  fprintf(fp, "# tmax      = %g\n", tmax);
  fprintf(fp, "# tinterval = %g\n", tinterval);
  fprintf(fp, "# tinc      = %g\n", tinc);
  fprintf(fp, "# teps      = %g\n", teps);
  fprintf(fp, "# \n");
  fprintf(fp, "# eta at time zero: = %g\n", rateparam.rI);
  fprintf(fp, "# sI:  = %g\n", rateparam.sI);
  fprintf(fp, "# sD:  = %g\n", rateparam.sD);
  fprintf(fp, "# vI:  = %g\n", rateparam.vI);
  fprintf(fp, "# vD:  = %g\n", rateparam.vD);
  fprintf(fp, "# muAM = %g\n", rateparam.muAM);
  fprintf(fp, "# muAD = %g\n", rateparam.muAD);
  fprintf(fp, "# muAI = %g\n", rateparam.muAI);
  fprintf(fp, "# ldI  = %g\n", rateparam.ldI);
  fprintf(fp, "# muI  = %g\n", rateparam.muI);
  fprintf(fp, "# ldEM = %g\n", rateparam.ldEM);
  fprintf(fp, "# muEM = %g\n", rateparam.muEM);
  fprintf(fp, "# ldED = %g\n", rateparam.ldED);
  fprintf(fp, "# muED = %g\n", rateparam.muED);
  
  e1sim(go, r, fp, fitfp, L, N, tmin, tmax, tinterval, tinc, teps, evomodel, rateparam, joint, rev, tol, errbuf, verbose);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  fclose(fp);
  if (fitfp) fclose(fitfp);
  return 0;
}


static int
e1sim(ESL_GETOPTS  *go, ESL_RANDOMNESS *r, FILE *fp, FILE *fitfp, int L, int N, double tmin, double tmax, double tinterval, double tinc, double teps,
      EVOM evomodel, struct rateparam_s rateparam, int joint, int rev, double tol, char *errbuf, int verbose)
{
  char           *msg  = "e1sim unit test failed";
  E1_RATE         *R   = NULL;
  ESQ           **esq  = NULL;
  ESL_HISTOGRAM   *h = NULL;
  double           tbeg;               /* time beging */
  double           tend;               /* time end    */
  double           time;
  double           sE_the;             /*          number       of surviving ancestral residues - theorerical */
  double           sE_fin, sV_fin;     /* expected number (std) of surviving ancestral residues - sampling from the finite time model */
  double           sE_inf, sV_inf;     /* expected number (std) of surviving ancestral residues - sampling from the infinitesimal model */
  double           nE_the;             /*          number       of inserted bloks (inserts) - theorerical */
  double           nE_fin, nV_fin;     /* expected number (std) of inserted bloks (inserts) - sampling from the finite time model */
  double           nE_inf, nV_inf;     /* expected number (std) of inserted bloks (inserts) - sampling from the infinitesimal model */
  double           eE_the;             /*          number       of inserted residues - theorerical */
  double           eE_fin, eV_fin;     /* expected number (std) of inserted residues - sampling from the finite time model */
  double           eE_inf, eV_inf;     /* expected number (std) of inserted residues - sampling from the infinitesimal model */
  double           lE_the;             /*          number       of total    residues - theorerical */
  double           lE_fin, lV_fin;     /* expected number (std) of total    residues - sampling from the finite time model */
  double           lE_inf, lV_inf;     /* expected number (std) of total residues - sampling from the infinitesimal model */
  double           gamma_the;
  double           gammaE_fin, gammaV_fin;
  double           gammaE_inf, gammaV_inf;
  double           eta_the;
  double           etaE_fin, etaV_fin;
  double           etaE_inf, etaV_inf;
  double           beta_the;
  double           betaE_fin, betaV_fin;
  double           betaE_inf, betaV_inf;
  double          *gammaS_fin = NULL;
  double          *gammaS_inf = NULL;
  double          *betaS_fin = NULL;
  double          *betaS_inf = NULL;
  double          *etaS_fin = NULL;
  double          *etaS_inf = NULL;
  double          *sS_fin = NULL;
  double          *nS_fin = NULL;
  double          *eS_fin = NULL;
  double          *lS_fin = NULL;
  double          *sS_inf = NULL;
  double          *nS_inf = NULL;
  double          *eS_inf = NULL;
  double          *lS_inf = NULL;
  double           alen_the;
  double           alen_fin;
  double           alen_inf;
  int              ntimes = 0;
  int              nbad;
  int              n;
  int              status;

  if (rev) {
    if (rateparam.ldEM != 0. || rateparam.muEM != 0. || rateparam.rI != 0.) { 
      printf("not a linear model rI or ld or muE are not zero\n"); 
      esl_fatal(msg);
    }
    if (fabs(rateparam.muI-(rateparam.ldI+rateparam.muAM)) > 1e-6) { 
      printf("not a reversible model muI %f should be equal to ldI+muA = %f\n", rateparam.muI, rateparam.ldI+rateparam.muAM); 
      esl_fatal(msg); 
    }
  }

  ESL_ALLOC(sS_fin,     sizeof(double *) * N); if (sS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of sS_fin failed");
  ESL_ALLOC(nS_fin,     sizeof(double *) * N); if (nS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of nS_fin failed");
  ESL_ALLOC(eS_fin,     sizeof(double *) * N); if (eS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of eS_fin failed");
  ESL_ALLOC(lS_fin,     sizeof(double *) * N); if (lS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of lS_fin failed");
  ESL_ALLOC(sS_inf,     sizeof(double *) * N); if (sS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of sS_inf failed");
  ESL_ALLOC(nS_inf,     sizeof(double *) * N); if (nS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of nS_inf failed");
  ESL_ALLOC(eS_inf,     sizeof(double *) * N); if (eS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of eS_inf failed");
  ESL_ALLOC(lS_inf,     sizeof(double *) * N); if (lS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of lS_inf failed");
  ESL_ALLOC(gammaS_fin, sizeof(double *) * N); if (gammaS_fin == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of gammaS_fin failed");
  ESL_ALLOC(gammaS_inf, sizeof(double *) * N); if (gammaS_inf == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of gammaS_inf failed");
  ESL_ALLOC(betaS_fin,  sizeof(double *) * N); if (betaS_fin  == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of beatS_fin failed");
  ESL_ALLOC(betaS_inf,  sizeof(double *) * N); if (betaS_inf  == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of betaS_inf failed");
  ESL_ALLOC(etaS_fin,   sizeof(double *) * N); if (etaS_fin   == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of etaS_fin failed");
  ESL_ALLOC(etaS_inf,   sizeof(double *) * N); if (etaS_inf   == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of etaS_inf failed");

  /* initialize */
  esl_vec_DSet(sS_fin,  N, 0.0);
  esl_vec_DSet(nS_fin,  N, 0.0);
  esl_vec_DSet(eS_fin,  N, 0.0);
  esl_vec_DSet(lS_fin,  N, 0.0);

  esl_vec_DSet(sS_inf,  N, 0.0);
  esl_vec_DSet(nS_inf,  N, 0.0);
  esl_vec_DSet(eS_inf,  N, 0.0);
  esl_vec_DSet(lS_inf,  N, 0.0);

  R = e1_rate_CreateWithValues(NULL, evomodel, rateparam, NULL, NULL, TRUE, tol, errbuf, FALSE);
  if (R == NULL) { printf("%s\n", errbuf); esl_fatal(msg); }
  printf("Simulations with model: %s\n", e1_rate_EvomodelType(R->evomodel));

  /* generate the ancestral sequences */
  status = generate_ancestral(r, R, &esq, &L, N, joint, rev, errbuf, verbose);
  if (status != eslOK) { printf("generate_ancestral() %s\n", errbuf); esl_fatal(msg); }
  
 /* first time interval */
  tbeg = tmin; 
  tend = tbeg + tinterval;
  
  while (tend < tmax+tinc) {
    
    /* create sequences at tbeg to start the inf evolutionary process in region (tbeg,tend) */
    if (evomodel != GG && evomodel != AG) {
      status = e1_sim_FiniteTime(r, R, tbeg, N, esq, esq, NULL, NULL, NULL, NULL, tol, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
    }
 
    nbad = 0;
    time = tbeg + tinc;
    while (time <= tend) {    
      
      ntimes ++;
      h = esl_histogram_CreateFull(-1.0, (double)L, 1.0);
      
      status = e1_sim_Theoretical(R, time, L, &sE_the, &nE_the, &eE_the, &lE_the, &gamma_the, &eta_the, &beta_the, tol, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
      status = e1_sim_FiniteTime(r, R, time, N, esq, NULL, sS_fin, nS_fin, eS_fin, lS_fin, tol, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
      status = e1_sim_Infinitesimal(r, R,  N, esq, time-tinc, tinc, teps, sS_inf, nS_inf, eS_inf, lS_inf, &nbad, h, tol, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
     
      esl_stats_DMean(sS_fin, N, &sE_fin, &sV_fin); sV_fin = sqrt(sV_fin);
      esl_stats_DMean(nS_fin, N, &nE_fin, &nV_fin); nV_fin = sqrt(nV_fin);
      esl_stats_DMean(eS_fin, N, &eE_fin, &eV_fin); eV_fin = sqrt(eV_fin);
      esl_stats_DMean(lS_fin, N, &lE_fin, &lV_fin); lV_fin = sqrt(lV_fin);
       
      esl_stats_DMean(sS_inf, N, &sE_inf, &sV_inf); sV_inf = sqrt(sV_inf);
      esl_stats_DMean(nS_inf, N, &nE_inf, &nV_inf); nV_inf = sqrt(nV_inf);
      esl_stats_DMean(eS_inf, N, &eE_inf, &eV_inf); eV_inf = sqrt(eV_inf);
      esl_stats_DMean(lS_inf, N, &lE_inf, &lV_inf); lV_inf = sqrt(lV_inf);
       
      for (n = 0; n < N; n ++) gammaS_fin[n] = (esq[n]->L > 0.)? 1.0 - sS_fin[n]/esq[n]->L : 0.0;
      for (n = 0; n < N; n ++) gammaS_inf[n] = (esq[n]->L > 0.)? 1.0 - sS_inf[n]/esq[n]->L : 0.0;
      esl_stats_DMean(gammaS_fin, N, &gammaE_fin, &gammaV_fin); gammaV_fin = sqrt(gammaV_fin);
      esl_stats_DMean(gammaS_inf, N, &gammaE_inf, &gammaV_inf); gammaV_inf = sqrt(gammaV_inf);
      
      for (n = 0; n < N; n ++) betaS_fin[n] = nS_fin[n]/(esq[n]->L + 1.);
      for (n = 0; n < N; n ++) betaS_inf[n] = nS_inf[n]/(esq[n]->L + 1.);
      esl_stats_DMean(betaS_fin, N, &betaE_fin, &betaV_fin); betaV_fin = sqrt(betaV_fin);
      esl_stats_DMean(betaS_inf, N, &betaE_inf, &betaV_inf); betaV_inf = sqrt(betaV_inf);
      
      for (n = 0; n < N; n ++) etaS_fin[n] = (eS_fin[n] > 0.)? 1.0 - nS_fin[n]/eS_fin[n] : 0.0;
      for (n = 0; n < N; n ++) etaS_inf[n] = (eS_inf[n] > 0.)? 1.0 - nS_inf[n]/eS_inf[n] : 0.0;
      esl_stats_DMean(etaS_fin, N, &etaE_fin, &etaV_fin); etaV_fin = sqrt(etaV_fin);
      esl_stats_DMean(etaS_inf, N, &etaE_inf, &etaV_inf); etaV_inf = sqrt(etaV_inf);
      
      if (joint) {
	alen_the  = L + eE_the;
	sE_the   /= L;
	nE_the   /= (L+1);
	eE_the   /= alen_the;
	lE_the   /= L;

	for (n = 0; n < N; n ++) {
	  alen_fin   = (double)esq[n]->L + eS_fin[n];
	  sS_fin[n] /= (double)esq[n]->L;
	  nS_fin[n] /= ((double)esq[n]->L + 1.);
	  eS_fin[n] /= alen_fin;
	  lS_fin[n] /= (double)esq[n]->L;

	  alen_inf   = (double)esq[n]->L + eS_inf[n];
	  sS_inf[n] /= (double)esq[n]->L;
	  nS_inf[n] /= ((double)esq[n]->L + 1.);
	  eS_inf[n] /= alen_inf;
	  lS_inf[n] /= (double)esq[n]->L;
	}
 
	esl_stats_DMean(sS_fin, N, &sE_fin, &sV_fin); sV_fin = sqrt(sV_fin);
	esl_stats_DMean(nS_fin, N, &nE_fin, &nV_fin); nV_fin = sqrt(nV_fin);
	esl_stats_DMean(eS_fin, N, &eE_fin, &eV_fin); eV_fin = sqrt(eV_fin);
	esl_stats_DMean(lS_fin, N, &lE_fin, &lV_fin); lV_fin = sqrt(lV_fin);
	
	esl_stats_DMean(sS_inf, N, &sE_inf, &sV_inf); sV_inf = sqrt(sV_inf);
	esl_stats_DMean(nS_inf, N, &nE_inf, &nV_inf); nV_inf = sqrt(nV_inf);
	esl_stats_DMean(eS_inf, N, &eE_inf, &eV_inf); eV_inf = sqrt(eV_inf);
	esl_stats_DMean(lS_inf, N, &lE_inf, &lV_inf); lV_inf = sqrt(lV_inf);
     }

      if (1||verbose) fprintf(stdout, "\n n=%d\n", ntimes);
      if (1||verbose) fprintf(stdout, "%f %d | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f | %d \n", time, L, 
			      sE_the, sE_fin, sV_fin, sE_inf, sV_inf, 
			      nE_the, nE_fin, nV_fin, nE_inf, nV_inf, 
			      eE_the, eE_fin, eV_fin, eE_inf, eV_inf, nbad);
      if (1||verbose) fprintf(stdout, "%f %d | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f | etaslp %.2f %.2f %.2f \n", time, L, 
			      gamma_the, gammaE_fin, gammaV_fin, gammaE_inf, gammaV_inf, 
			      beta_the,  betaE_fin, betaV_fin,  betaE_inf, betaV_inf,
			      eta_the,   etaE_fin,  etaV_fin,   etaE_inf,  etaV_inf,
			      eta_the/time,  etaE_fin/time,  etaE_inf/time);

      fprintf(fp, "%f %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", time, L, 
	      sE_the, sE_fin, sV_fin, sE_inf, sV_inf, 
	      nE_the, nE_fin, nV_fin, nE_inf, nV_inf, 
	      eE_the, eE_fin, eV_fin, eE_inf, eV_inf, 
	      lE_the, lE_fin, lV_fin, lE_inf, lV_inf,
	      gamma_the, gammaE_fin, gammaV_fin, gammaE_inf, gammaV_inf,  
	      beta_the,  betaE_fin, betaV_fin,  betaE_inf, betaV_inf,
	      eta_the,   etaE_fin,  etaV_fin,   etaE_inf,  etaV_inf);
      
       time += (time <= 100.)? tinc : ( (time <= 1000.)? tinc*10. : tinc*1000. );

      status = GoodnessGeometric(fitfp, h, NULL, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }

      esl_histogram_Destroy(h); h = NULL;
    }

    tbeg  = tend;
    tend += tinterval;
  }

  e1_rate_Destroy(R);
  for (n = 0; n < N; n++) e1_sim_ESQDestroy(esq[n]);
  free(esq);
  if (h)          esl_histogram_Destroy(h);
  if (sS_fin)     free(sS_fin);
  if (nS_fin)     free(nS_fin);
  if (eS_fin)     free(eS_fin);
  if (lS_fin)     free(lS_fin);
  if (sS_inf)     free(sS_inf);
  if (nS_inf)     free(nS_inf);
  if (eS_inf)     free(eS_inf);
  if (lS_inf)     free(lS_inf);
  if (gammaS_fin) free(gammaS_fin);
  if (gammaS_inf) free(gammaS_inf);
  if (betaS_fin)  free(betaS_fin);
  if (betaS_inf)  free(betaS_inf);
  if (etaS_fin)   free(etaS_fin);
  if (etaS_inf)   free(etaS_inf);
  return eslOK;

 ERROR:
  if (R) e1_rate_Destroy(R);
  if (esq) {
    for (n = 0; n < N; n++) if (esq[n]) e1_sim_ESQDestroy(esq[n]);
    free(esq);
  }
  if (h)          esl_histogram_Destroy(h);
  if (sS_fin)     free(sS_fin);
  if (nS_fin)     free(nS_fin);
  if (eS_fin)     free(eS_fin);
  if (sS_inf)     free(sS_inf);
  if (nS_inf)     free(nS_inf);
  if (eS_inf)     free(eS_inf);
  if (gammaS_fin) free(gammaS_fin);
  if (gammaS_inf) free(gammaS_inf);
  if (betaS_fin)  free(betaS_fin);
  if (betaS_inf)  free(betaS_inf);
  if (etaS_fin)   free(etaS_fin);
  if (etaS_inf)   free(etaS_inf);
   return status;
}

static int
generate_ancestral(ESL_RANDOMNESS *r, E1_RATE *R, ESQ ***ret_esq, int *L, int N, int joint, int rev, char *errbuf, int verbose)
{
  ESQ    **esq = NULL;
  double   p;
  double   l = (double)(*L);
  int      len;
  int      n;
  int      status;

  /* override L and use a length that will make the linear
   * model (ie ld=0^muE=0^rI=0) reversible 
   */
  if (rev) l = R->ldE[e1R_I] / (R->muA[e1R_S] - R->ldE[e1R_I]);
  if (l < 0.) ESL_XFAIL(eslFAIL, errbuf, "generate_ancestral(): negative length l=%d. ldI %f muAM=%f", (int)l, R->ldE[e1R_I], R->muA[e1R_S]); 

  p = l/(l+1.0);

  ESL_ALLOC(esq, sizeof(ESQ *) * N);
  for (n = 0; n < N; n ++) esq[n] = NULL;
  
  for (n = 0; n < N; n++) {
    /* fix length sequences or 
     * sample the length from the geometric distribution */     
    if (!joint) len = (int)l; 
    else        len = (p == 0.0)? 1 : (int)floor(log(esl_random(r)) / log(p) + 1.);      

    esq[n] = e1_sim_ESQCreate(len);
    if (esq[n] == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of ESQ failed");
  }

  *L = (int)ceil(l);
  *ret_esq = esq;
  return eslOK;

 ERROR:
  if (esq) free(esq);
  return status;
}



/*****************************************************************
 * @LICENSE@
 *
 *****************************************************************/
