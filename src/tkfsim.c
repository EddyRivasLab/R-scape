/* TKFsim -- evolve according to the TKF evolutionary model.
 *          both using the finite-form model, and the rates at infinitesimal time to compare.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stats.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2.h"
#include "e1_goodness.h"
#include "tkf_rate.h"
#include "tkf_simulate.h"


static ESL_OPTIONS options[] = {
  /* name           type              default  env       range      toggles     reqs        incomp        help                                                    docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "show brief help on version and usage",                        0 },
  { "--geofit",   eslARG_STRING,         NULL, NULL,      NULL,     NULL,       NULL,       NULL,         "geometric goodness of fit",                                   0 },
  /* parameters to control the simulation */
  { "-L",           eslARG_INT,        "1000", NULL,     "n>0",     NULL,       NULL,       NULL,         "length of ancestral sequence",                                0 },
  { "-N",           eslARG_INT,       "10000", NULL,     "n>0",     NULL,       NULL,       NULL,         "number of sequences for a given set of parameters and time",  0 },
  { "--joint",      eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "generate ancestral sqs with geom L/(L+1)",                    0 },
  { "--rev",        eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "impose the condition L = ldI/(muD-ldI)  ",                    0 },
  { "--tmin",       eslARG_REAL,        "0.0", NULL,    "x>=0",     NULL,       NULL,       NULL,         "min time to evolve sequence",                                 0 },
  { "--tmax",       eslARG_REAL,        "3.0", NULL,     "x>0",     NULL,       NULL,       NULL,         "max time to evolve sequence",                                 0 },
  { "--tinterval",  eslARG_REAL,       "0.20", NULL,     "x>0",     NULL,       NULL,       NULL,         "time interval to perform a given infinitesimal calculation",  0 },
  { "--tinc",       eslARG_REAL,       "0.01", NULL,     "x>0",     NULL,       NULL,       NULL,         "time increments to calculate statistics",                     0 },
  { "--teps",       eslARG_REAL,       "1e-6", NULL,     "x>0",     NULL,       NULL,       NULL,         "small time increments for the infinitesimal simulation",      0 },
  /* parameters of the evolutionary model */
  { "--etaz",       eslARG_REAL,       "0.01", NULL, "0<=x<=1",     NULL,       NULL,       NULL,         "eta at time zero",                                            0 },
  { "--mu",         eslARG_REAL,       "0.01", NULL,    "x>=0",     NULL,       NULL,       NULL,         "rate of deletions",                                           0 },
  { "--ld",        eslARG_REAL,       "0.02",  NULL,    "x>=0",     NULL,       NULL,       NULL,         "rate of insertions",                                               0 },
  { "-v",           eslARG_NONE,        FALSE, NULL,      NULL,     NULL,       NULL,       NULL,         "be verbose",                                                  0 },
  /* other options */
  { "--tol",        eslARG_REAL,     "0.0001", NULL,      NULL,     NULL,       NULL,       NULL,         "tolerance",                                                   0 },
  { "--seed",        eslARG_INT,          "0", NULL,    "n>=0",     NULL,       NULL,       NULL,         "set RNG seed to <n>",                                         5 },
   {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <outfile>";
static char banner[] = "e1sim";


static int tkfsim(ESL_GETOPTS *go, ESL_RANDOMNESS *r, FILE *fp, FILE *fitfp, int L, int N, 
		  double tmin, double tmax, double tinterval, double tinc, double teps, 
		  double mu, double ld, double etaz, int joint, int rev, double tol, char *errbuf, int verbose);
static int generate_ancestral(ESL_RANDOMNESS *r, TKF_RATE *R,  ESQ ***ret_esq, int *L, int N, int joint, int rev, char *errbuf, int verbose);

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS *r = NULL;	  
  char            errbuf[eslERRBUFSIZE];
  char           *outfile;
  char           *fitfile;
  FILE           *fp = NULL;                /* output file handle  */
  FILE           *fitfp = NULL;             /* goodnes of fit file handle  */
  double          mu;
  double          ld;
  double          etaz;
  double          tmin;
  double          tmax;
  double          tinterval;
  double          tinc;
  double          teps;      /* tespsilon for infinitesimal simulations */
  double          tol;
  int             N;
  int             L;
  int             joint;
  int             rev;
  int             verbose;

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
  etaz      = esl_opt_GetReal(go, "--etaz");
  mu        = esl_opt_GetReal(go, "--mu");
  ld        = esl_opt_GetReal(go, "--ld");
  verbose   = esl_opt_GetBoolean(go, "-v");
  tol       = esl_opt_GetReal(go, "--tol");
 
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
  fprintf(stdout, "# eta at time zero: = %g\n", etaz);
  fprintf(stdout, "# mu  = %g\n", mu);
  fprintf(stdout, "# ld  = %g\n", ld);

  if (!joint) fprintf(fp, "# L         = %d\n", L);
  else        fprintf(fp, "# p         = %f\n", (double)L/((double)L+1.));
  fprintf(fp, "# N         = %d\n", N);
  fprintf(fp, "# tmin      = %g\n", tmin);
  fprintf(fp, "# tmax      = %g\n", tmax);
  fprintf(fp, "# tinterval = %g\n", tinterval);
  fprintf(fp, "# tinc      = %g\n", tinc);
  fprintf(fp, "# teps      = %g\n", teps);
  fprintf(fp, "# \n");
  fprintf(fp, "# r       = %g\n", etaz);
  fprintf(fp, "# mu_tkf  = %g\n", mu);
  fprintf(fp, "# ld_tkf  = %g\n", ld);
  
  tkfsim(go, r, fp, fitfp, L, N, tmin, tmax, tinterval, tinc, teps, mu, ld, etaz, joint, rev, tol, errbuf, verbose);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  fclose(fp);
  if (fitfp) fclose(fitfp);
  return 0;
}


static int
tkfsim(ESL_GETOPTS  *go, ESL_RANDOMNESS *r, FILE *fp, FILE *fitfp, int L, int N, double tmin, double tmax, double tinterval, double tinc, double teps,
       double mu, double ld, double etaz, int joint, int rev, double tol, char *errbuf, int verbose)
{
  char            *msg  = "e1sim unit test failed";
  TKF_RATE        *R   = NULL;
  ESQ            **esq  = NULL;
  ESL_HISTOGRAM   *h = NULL;
  double           tbeg;               /* time beging */
  double           tend;               /* time end    */
  double           time;
  double           sE_the;             /*          number       of surviving residues - theorerical */
  double           sE_fin, sV_fin;     /* expected number (std) of surviving residues - sampling from the finite time model */
  double           sE_inf, sV_inf;     /* expected number (std) of surviving residues - sampling from the infinitesimal model */
  double           d0E_the;            /*          number       of lone deletions     - theorerical */
  double           d0E_fin, d0V_fin;   /* expected number (std) of lone deletions     - sampling from the finite time model */
  double           d0E_inf, d0V_inf;   /* expected number (std) of lone deletions     - sampling from the infinitesimal model */
  double           d1E_the;            /*          number       of deletions followed by ins     - theorerical */
  double           d1E_fin, d1V_fin;   /* expected number (std) of deletions followed by ins     - sampling from the finite time model */
  double           d1E_inf, d1V_inf;   /* expected number (std) of deletions followed by ins     - sampling from the infinitesimal model */
  double           iE_the;             /*          number       of inserted residues  - theorerical */
  double           iE_fin, iV_fin;     /* expected number (std) of inserted residues  - sampling from the finite time model */
  double           iE_inf, iV_inf;     /* expected number (std) of inserted residues  - sampling from the infinitesimal model */
  double           nE_the;             /*          number       of inserted bloks (inserts) - theorerical */
  double           nE_fin, nV_fin;     /* expected number (std) of inserted bloks (inserts) - sampling from the finite time model */
  double           nE_inf, nV_inf;     /* expected number (std) of inserted bloks (inserts) - sampling from the infinitesimal model */
  double           eE_the;             /*          number       of inserted residues - theorerical */
  double           eE_fin, eV_fin;     /* expected number (std) of inserted residues - sampling from the finite time model */
  double           eE_inf, eV_inf;     /* expected number (std) of inserted residues - sampling from the infinitesimal model */
  double           lE_the;             /*          number       of total    residues  - theorerical */
  double           lE_fin, lV_fin;     /* expected number (std) of total    residues  - sampling from the finite time model */
  double           lE_inf, lV_inf;     /* expected number (std) of total residues - sampling from the infinitesimal model */
  double           gamma_the;
  double           gammaE_fin, gammaV_fin;
  double           gammaE_inf, gammaV_inf;
  double           beta_the;
  double           heta_the;
  double           eta_the;
  double           betaE_fin, betaV_fin;
  double           betaE_inf, betaV_inf;
  double           hetaE_fin, hetaV_fin;
  double           hetaE_inf, hetaV_inf;
  double           etaE_fin, etaV_fin;
  double           etaE_inf, etaV_inf;
  double          *gammaS_fin = NULL;
  double          *gammaS_inf = NULL;
  double          *betaS_fin = NULL;
  double          *betaS_inf = NULL;
  double          *hetaS_fin = NULL;
  double          *hetaS_inf = NULL;
  double          *etaS_fin = NULL;
  double          *etaS_inf = NULL;
  double           alen_the;
  double           alen_fin;
  double           alen_inf;
  double          *sS_fin  = NULL;
  double          *d0S_fin = NULL;
  double          *d1S_fin = NULL;
  double          *iS_fin  = NULL;
  double          *nS_fin  = NULL;
  double          *eS_fin  = NULL;
  double          *lS_fin  = NULL;
  double          *sS_inf  = NULL;
  double          *d0S_inf = NULL;
  double          *d1S_inf = NULL;
  double          *iS_inf  = NULL;
  double          *nS_inf  = NULL;
  double          *eS_inf  = NULL;
  double          *lS_inf  = NULL;
  int              ntimes = 0;
  int              nbad;
  int              n;
  int              status;

  if (rev) {
    if (ld >= mu) { printf("not a reversible model lambda %f should be smaller than mu = %f\n", ld, mu); esl_fatal(msg); }
  }

  R = tkf_rate_CreateWithValues(NULL, mu, ld, etaz, NULL, NULL, TRUE, tol, errbuf, FALSE);
  if (R == NULL) { printf("%s\n", errbuf); esl_fatal(msg); }

  /* generate the ancestral sequences */
  status = generate_ancestral(r, R, &esq, &L, N, joint, rev, errbuf, verbose);
  if (status != eslOK) { printf("generate_ancestral() error\n"); esl_fatal(msg); }

  ESL_ALLOC(sS_fin,     sizeof(double *) * N); if (sS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of sS_fin failed");
  ESL_ALLOC(d0S_fin,    sizeof(double *) * N); if (d0S_fin    == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of d0S_fin failed");
  ESL_ALLOC(d1S_fin,    sizeof(double *) * N); if (d1S_fin    == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of d1S_fin failed");
  ESL_ALLOC(iS_fin,     sizeof(double *) * N); if (iS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of iS_fin failed");
  ESL_ALLOC(nS_fin,     sizeof(double *) * N); if (nS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of nS_fin failed");
  ESL_ALLOC(eS_fin,     sizeof(double *) * N); if (eS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of eS_fin failed");
  ESL_ALLOC(lS_fin,     sizeof(double *) * N); if (lS_fin     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of lS_fin failed");
  ESL_ALLOC(sS_inf,     sizeof(double *) * N); if (sS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of sS_fin failed");
  ESL_ALLOC(d0S_inf,    sizeof(double *) * N); if (d0S_inf    == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of d0S_inf failed");
  ESL_ALLOC(d1S_inf,    sizeof(double *) * N); if (d1S_inf    == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of d1S_inf failed");
  ESL_ALLOC(iS_inf,     sizeof(double *) * N); if (iS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of iS_inf failed");
  ESL_ALLOC(nS_inf,     sizeof(double *) * N); if (nS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of nS_inf failed");
  ESL_ALLOC(eS_inf,     sizeof(double *) * N); if (eS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of eS_inf failed");
  ESL_ALLOC(lS_inf,     sizeof(double *) * N); if (lS_inf     == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of lS_inf failed");
  ESL_ALLOC(gammaS_fin, sizeof(double *) * N); if (gammaS_fin == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of gammaS_fin failed");
  ESL_ALLOC(gammaS_inf, sizeof(double *) * N); if (gammaS_inf == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of gammaS_inf failed");
  ESL_ALLOC(betaS_fin,  sizeof(double *) * N); if (betaS_fin  == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of beatS_fin failed");
  ESL_ALLOC(betaS_inf,  sizeof(double *) * N); if (betaS_inf  == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of betaS_inf failed");
  ESL_ALLOC(hetaS_fin,  sizeof(double *) * N); if (hetaS_fin  == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of heatS_fin failed");
  ESL_ALLOC(hetaS_inf,  sizeof(double *) * N); if (hetaS_inf  == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of hetaS_inf failed");
  ESL_ALLOC(etaS_fin,   sizeof(double *) * N); if (etaS_fin   == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of eatS_fin failed");
  ESL_ALLOC(etaS_inf,   sizeof(double *) * N); if (etaS_inf   == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of etaS_inf failed");
 
  /* initialize */
  esl_vec_DSet(sS_fin,  N, 0.0);
  esl_vec_DSet(d0S_fin, N, 0.0);
  esl_vec_DSet(d1S_fin, N, 0.0);
  esl_vec_DSet(iS_fin,  N, 0.0);
  esl_vec_DSet(nS_fin,  N, 0.0);
  esl_vec_DSet(eS_fin,  N, 0.0);
  esl_vec_DSet(lS_fin,  N, 0.0);

  esl_vec_DSet(sS_inf,  N, 0.0);
  esl_vec_DSet(d0S_inf, N, 0.0);
  esl_vec_DSet(d1S_inf, N, 0.0);
  esl_vec_DSet(iS_inf,  N, 0.0);
  esl_vec_DSet(nS_inf,  N, 0.0);
  esl_vec_DSet(eS_inf,  N, 0.0);
  esl_vec_DSet(lS_inf,  N, 0.0);

   /* first time interval */
  tbeg = tmin; 
  tend = tbeg + tinterval;

   while (tend <= tmax) {
    
    /* create sequences at atbeg to start the inf evolutionary process in region (tbeg,tend) */
     status = tkf_sim_FiniteTime(r, R, tbeg, L, N, esq, esq, NULL, NULL, NULL, NULL, NULL, NULL, NULL, tol, errbuf, FALSE);
    if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }

    nbad = 0;
    time = tbeg + tinc;
    while (time <= tend) {    
      
      ntimes ++;
      h = esl_histogram_CreateFull(-1.0, (double)L, 1.0);
      
      status = tkf_sim_Theoretical(R, time, L, &sE_the, &d0E_the, &d1E_the, &iE_the, &nE_the, &eE_the, &lE_the, &gamma_the, &beta_the, &heta_the, &eta_the, tol, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
      status = tkf_sim_FiniteTime(r, R, time, L, N, esq, NULL, sS_fin, d0S_fin, d1S_fin, iS_fin, nS_fin, eS_fin, lS_fin, tol, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
      status = tkf_sim_Infinitesimal(r, R, N, esq, time-tinc, tinc, teps, sS_inf, d0S_inf, d1S_inf, iS_inf, nS_inf, eS_inf, lS_inf, &nbad, h, tol, errbuf, FALSE);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
      
      esl_stats_DMean(sS_fin,  N, &sE_fin,  &sV_fin);  sV_fin  = sqrt(sV_fin);
      esl_stats_DMean(d0S_fin, N, &d0E_fin, &d0V_fin); d0V_fin = sqrt(d0V_fin);
      esl_stats_DMean(d1S_fin, N, &d1E_fin, &d1V_fin); d1V_fin = sqrt(d1V_fin);
      esl_stats_DMean(iS_fin,  N, &iE_fin,  &iV_fin);  iV_fin  = sqrt(iV_fin);
      esl_stats_DMean(nS_fin,  N, &nE_fin,  &nV_fin);  nV_fin  = sqrt(nV_fin);
      esl_stats_DMean(eS_fin,  N, &eE_fin,  &eV_fin);  eV_fin  = sqrt(eV_fin);
      esl_stats_DMean(lS_fin,  N, &lE_fin,  &lV_fin);  lV_fin  = sqrt(lV_fin);
       
      esl_stats_DMean(sS_inf,  N, &sE_inf,  &sV_inf);  sV_inf  = sqrt(sV_inf);
      esl_stats_DMean(d0S_inf, N, &d0E_inf, &d0V_inf); d0V_inf = sqrt(d0V_inf);
      esl_stats_DMean(d1S_inf, N, &d1E_inf, &d1V_inf); d1V_inf = sqrt(d1V_inf);
      esl_stats_DMean(iS_inf,  N, &iE_inf,  &iV_inf);  iV_inf  = sqrt(iV_inf);
      esl_stats_DMean(nS_inf,  N, &nE_inf,  &nV_inf);  nV_inf  = sqrt(nV_inf);
      esl_stats_DMean(eS_inf,  N, &eE_inf,  &eV_inf);  eV_inf  = sqrt(eV_inf);
      esl_stats_DMean(lS_inf,  N, &lE_inf,  &lV_inf);  lV_inf  = sqrt(lV_inf);
       
      for (n = 0; n < N; n ++) gammaS_fin[n] = (esq[n]->L > 0.)? 1.0 - (double)sS_fin[n]/(double)esq[n]->L : 0.0;
      for (n = 0; n < N; n ++) gammaS_inf[n] = (esq[n]->L > 0.)? 1.0 - (double)sS_inf[n]/(double)esq[n]->L : 0.0;
      esl_stats_DMean(gammaS_fin, N, &gammaE_fin, &gammaV_fin); gammaV_fin = sqrt(gammaV_fin);
      esl_stats_DMean(gammaS_inf, N, &gammaE_inf, &gammaV_inf); gammaV_inf = sqrt(gammaV_inf);
      
      for (n = 0; n < N; n ++) betaS_fin[n] = (esq[n]->L     > 0.)?                       (double)d0S_fin[n] * R->lambda / ( (double)esq[n]->L * R->mu )         : 0.0;
      for (n = 0; n < N; n ++) betaS_inf[n] = (esq[n]->L     > 0.)?                       (double)d0S_inf[n] * R->lambda / ( (double)esq[n]->L * R->mu )         : 0.0;
      for (n = 0; n < N; n ++) hetaS_fin[n] = (esq[n]->L     > 0. && gammaS_fin[n] > 0.)? (double)d1S_fin[n]             / ( (double)esq[n]->L * gammaS_fin[n] ) : 0.0;
      for (n = 0; n < N; n ++) hetaS_inf[n] = (esq[n]->L     > 0. && gammaS_inf[n] > 0.)? (double)d1S_inf[n]             / ( (double)esq[n]->L * gammaS_inf[n] ) : 0.0;
      for (n = 0; n < N; n ++) etaS_fin[n]  = R->etaz + (1.0 - R->etaz)*betaS_fin[n];
      for (n = 0; n < N; n ++) etaS_inf[n]  = R->etaz + (1.0 - R->etaz)*betaS_inf[n];

      esl_stats_DMean(betaS_fin, N, &betaE_fin, &betaV_fin); betaV_fin = sqrt(betaV_fin);
      esl_stats_DMean(betaS_inf, N, &betaE_inf, &betaV_inf); betaV_inf = sqrt(betaV_inf);
      esl_stats_DMean(hetaS_fin, N, &hetaE_fin, &hetaV_fin); hetaV_fin = sqrt(hetaV_fin);
      esl_stats_DMean(hetaS_inf, N, &hetaE_inf, &hetaV_inf); hetaV_inf = sqrt(hetaV_inf);
      esl_stats_DMean(etaS_fin,  N, &etaE_fin,  &etaV_fin);  etaV_fin  = sqrt(etaV_fin);
      esl_stats_DMean(etaS_inf,  N, &etaE_inf,  &etaV_inf);  etaV_inf  = sqrt(etaV_inf);
      
      if (joint) {
 	alen_the  = (double)L + d1E_the + iE_the;
	sE_the   /= (double)L;
	d0E_the  /= (double)L;
	d1E_the  /= (double)L;
	iE_the   /= alen_the;
	lE_the   /= (double)L;

	for (n = 0; n < N; n ++) {
	  alen_fin    = (double)esq[n]->L + (double)d1S_fin[n] + (double)iS_fin[n];
	  sS_fin[n]  /= (double)esq[n]->L;
	  d0S_fin[n] /= (double)esq[n]->L;
	  d1S_fin[n] /= (double)esq[n]->L;
	  iS_fin[n]  /= alen_fin;
	  nS_fin[n]  /= ((double)esq[n]->L + 1.);
	  eS_fin[n]  /= alen_fin;
	  lS_fin[n]  /= (double)esq[n]->L;

	  alen_inf    = (double)esq[n]->L + (double)d1S_inf[n] + (double)iS_inf[n];
	  sS_inf[n]  /= (double)esq[n]->L;
	  d0S_inf[n] /= (double)esq[n]->L;
	  d1S_inf[n] /= (double)esq[n]->L;
	  iS_inf[n]  /= alen_inf;
	  nS_inf[n]  /= ((double)esq[n]->L + 1.);
	  eS_inf[n]  /= alen_inf;
	  lS_inf[n]  /= (double)esq[n]->L;
	}
 
	esl_stats_DMean(sS_fin,  N, &sE_fin,  &sV_fin);  sV_fin  = sqrt(sV_fin);
	esl_stats_DMean(d0S_fin, N, &d0E_fin, &d0V_fin); d0V_fin = sqrt(d0V_fin);
	esl_stats_DMean(d1S_fin, N, &d1E_fin, &d1V_fin); d1V_fin = sqrt(d1V_fin);
	esl_stats_DMean(iS_fin,  N, &iE_fin,  &iV_fin);  iV_fin  = sqrt(iV_fin);
	esl_stats_DMean(nS_fin,  N, &nE_fin,  &nV_fin);  nV_fin  = sqrt(nV_fin);
	esl_stats_DMean(eS_fin,  N, &eE_fin,  &eV_fin);  eV_fin  = sqrt(eV_fin);
 	esl_stats_DMean(lS_fin,  N, &lE_fin,  &lV_fin);  lV_fin  = sqrt(lV_fin);
	
	esl_stats_DMean(sS_inf,  N, &sE_inf,  &sV_inf);  sV_inf  = sqrt(sV_inf);
	esl_stats_DMean(d0S_inf, N, &d0E_inf, &d0V_inf); d0V_inf = sqrt(d0V_inf);
	esl_stats_DMean(d1S_inf, N, &d1E_inf, &d1V_inf); d1V_inf = sqrt(d1V_inf);
	esl_stats_DMean(iS_inf,  N, &iE_inf,  &iV_inf);  iV_inf  = sqrt(iV_inf);
	esl_stats_DMean(nS_inf,  N, &nE_inf,  &nV_inf);  nV_inf  = sqrt(nV_inf);
	esl_stats_DMean(eS_inf,  N, &eE_inf,  &eV_inf);  eV_inf  = sqrt(eV_inf);
 	esl_stats_DMean(lS_inf,  N, &lE_inf,  &lV_inf);  lV_inf  = sqrt(lV_inf);
     }

      if (1||verbose) fprintf(stdout, "\n n=%d\n", ntimes);
      if (1||verbose) fprintf(stdout, "%f %d | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f | %d\n", time, L, 
			      sE_the,  sE_fin,  sV_fin,  sE_inf,  sV_inf, 
			      nE_the,  nE_fin,  nV_fin,  nE_inf,  nV_inf, 
			      eE_the,  eE_fin,  eV_fin,  eE_inf,  eV_inf,  nbad);
      if (1||verbose) fprintf(stdout, "%f %d | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f %.2f\n", time, L, 
			      d1E_the, d1E_fin, d1V_fin, d1E_inf, d1V_inf, 
			      d0E_the, d0E_fin, d0V_fin, d0E_inf, d0V_inf, 
			      iE_the,  iE_fin,  iV_fin,  iE_inf,  iV_inf);
      if (1||verbose) fprintf(stdout, "%f %d | %f %f %f %f %f | %f %f %f %f %f | %f %f %f %f %f \n", time, L, 
			      gamma_the, gammaE_fin, gammaV_fin, gammaE_inf, gammaV_inf, 
			      beta_the,  betaE_fin, betaV_fin,  betaE_inf, betaV_inf, 
			      eta_the,   etaE_fin,  etaV_fin,   etaE_inf,  etaV_inf);

      fprintf(fp, "%f %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", time, L, 
	      sE_the,  sE_fin,  sV_fin,  sE_inf,  sV_inf, 
	      nE_the,  nE_fin,  nV_fin,  nE_inf,  nV_inf, 
	      eE_the,  eE_fin,  eV_fin,  eE_inf,  eV_inf, 
	      lE_the,  lE_fin,  lV_fin,  lE_inf,  lV_inf,
	      gamma_the, gammaE_fin, gammaV_fin, gammaE_inf, gammaV_inf,  
	      beta_the,  betaE_fin, betaV_fin,  betaE_inf, betaV_inf,  
	      eta_the,   etaE_fin,  etaV_fin,   etaE_inf,  etaV_inf);
      
      status = GoodnessGeometric(fitfp, h, NULL, errbuf, verbose);
      if (status != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }

      time += (time <= 10.)? tinc : ( (time <= 100.)? tinc*10. : ( (time <= 1000.)? tinc*100. : ((time <= 10000.)?  tinc*1000.: tinc*5000.)) );
     }

    tbeg  = tend;
    tend += tinterval;
  }

  tkf_rate_Destroy(R);
  for (n = 0; n < N; n++) e1_sim_ESQDestroy(esq[n]);
  free(esq);
  if (h)          esl_histogram_Destroy(h);
  if (sS_fin)     free(sS_fin);
  if (d0S_fin)    free(d0S_fin);
  if (d1S_fin)    free(d1S_fin);
  if (iS_fin)     free(iS_fin);
  if (nS_fin)     free(nS_fin);
  if (eS_fin)     free(eS_fin);
  if (sS_inf)     free(sS_inf);
  if (d0S_inf)    free(d0S_inf);
  if (d1S_inf)    free(d1S_inf);
  if (iS_inf)     free(iS_inf);
  if (nS_inf)     free(nS_inf);
  if (eS_inf)     free(eS_inf);
  if (lS_inf)     free(lS_inf);
  if (gammaS_fin) free(gammaS_fin);
  if (gammaS_inf) free(gammaS_inf);
  if (betaS_fin)  free(betaS_fin);
  if (betaS_inf)  free(betaS_inf);
  if (hetaS_fin)  free(hetaS_fin);
  if (hetaS_inf)  free(hetaS_inf);
  if (etaS_fin)   free(etaS_fin);
  if (etaS_inf)   free(etaS_inf);
  return eslOK;
  
 ERROR:
  if (R) tkf_rate_Destroy(R);
  if (esq) {
    for (n = 0; n < N; n++) if (esq[n]) e1_sim_ESQDestroy(esq[n]);
    free(esq);
  }
  if (h)          esl_histogram_Destroy(h);
  if (sS_fin)     free(sS_fin);
  if (d0S_fin)    free(d0S_fin);
  if (d1S_fin)    free(d1S_fin);
  if (iS_fin)     free(iS_fin);
  if (nS_fin)     free(nS_fin);
  if (eS_fin)     free(eS_fin);
  if (lS_fin)     free(lS_fin);
  if (sS_inf)     free(sS_inf);
  if (d0S_inf)    free(d0S_inf);
  if (d1S_inf)    free(d1S_inf);
  if (iS_inf)     free(iS_inf);
  if (nS_inf)     free(nS_inf);
  if (eS_inf)     free(eS_inf);
  if (lS_inf)     free(lS_inf);
  if (gammaS_fin) free(gammaS_fin);
  if (gammaS_inf) free(gammaS_inf);
  if (betaS_fin)  free(betaS_fin);
  if (betaS_inf)  free(betaS_inf);
  if (hetaS_fin)  free(hetaS_fin);
  if (hetaS_inf)  free(hetaS_inf);
  if (etaS_fin)   free(etaS_fin);
  if (etaS_inf)   free(etaS_inf);
  return status;
}

static int
generate_ancestral(ESL_RANDOMNESS *r, TKF_RATE *R, ESQ ***ret_esq, int *L, int N, int joint, int rev, char *errbuf, int verbose)
{
  ESQ    **esq = NULL;
  double   pzero;
  double   p;
  double   l = (double)(*L);
  int      len;
  int      n;
  int      status;

  pzero = R->lambda / R->mu;
  p     = pzero * (1.0 - R->etaz) + R->etaz;

  /* override L and use a length that will make
   * the model reversible
   */
  if (rev) l = pzero / (1.0 - p);
  if (l < 0) ESL_XFAIL(eslFAIL, errbuf, "TKF generate_ancestral(): negative length"); 

  ESL_ALLOC(esq, sizeof(ESQ *) * N);
  for (n = 0; n < N; n ++) esq[n] = NULL;
  
  for (n = 0; n < N; n++) {
    /* fix length sequences or 
     * sample the length from the geometric distribution 
     *  gamma_0 = 1 - pzero;
     *  gamma_n = pzero * (1 - p) p ^{n-1},  n > 1
     */     
    if (!joint) {
      len = (int)l; 
    }
    else {
      if ((1.0-pzero) < esl_random(r)) {
	len = (p == 0.0)? 1 : (int)floor(log(esl_random(r)) / log(p) + 1.); 
      }
    }     

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
