/* e1_simulate - funtions to evolve sequences either from
 *               the finite-time distribution or from the infinitesimal rate
 * 
 * Contents:
 *
 * ER, Tue Sep 27 13:19:30 2011 [Janelia] 
 * SVN $Id:$
 */

#include "p7_config.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_rootfinder.h"
#include "esl_stats.h"
#include "esl_vectorops.h"
#include "hmmer.h"

#include "e2.h"
#include "tkf_rate.h"
#include "tkf_model.h"
#include "tkf_simulate.h"
#include "e1_simulate.h"
#include "ratematrix.h"

static int tkf_finite_fate_residue(ESL_RANDOMNESS *r, TKF_MODEL *tkf, double tt, ESQ *Asq, ESQ *Dsq);
static int tkf_finite_fate_insert (ESL_RANDOMNESS *r, TKF_MODEL *tkf, double tt, ESQ *Asq, ESQ *Dsq);
static int tkf_inf_fate_residue   (ESL_RANDOMNESS *r, TKF_RATE  *R,   double te, ESQ *Asq, ESQ *Dsq);
static int tkf_inf_fate_insert    (ESL_RANDOMNESS *r, TKF_RATE  *R,   double te, ESQ *Asq, ESQ *Dsq, int *ret_nbad);

int 
tkf_sim_Theoretical(TKF_RATE *R, double time, int L, double *ret_sE, double *ret_d0E, double *ret_d1E, double *ret_iE, double *ret_nE, double *ret_eE, double *ret_lE, 
		    double *ret_gamma, double *ret_beta, double *ret_heta, double *ret_eta, double tol, char *errbuf, int verbose)
{
  TKF_MODEL *tkf = NULL;
  double     nogammat;
  double     etat;
  double     sE;
  double     d0E;
  double     d1E;
  double     iE;
  double     nE;
  double     eE;
  int        status;

  tkf = tkf_model_Create(R, time, 1.0, NULL, L, NULL, tol, verbose);
  if (tkf == NULL) ESL_XFAIL(eslFAIL, errbuf, "error creating tkfmodel");
    
  nogammat = 1.0 - tkf->R->mu * tkf->betat / (tkf->R->lambda * tkf->gammat);
  etat = tkf->R->etaz + (1.0 - tkf->R->etaz) * tkf->betat;

  sE  = (double)L * (1.0-tkf->gammat);
  d1E = (double)L * tkf->gammat * nogammat;
  d0E = (double)L * tkf->gammat * (1.0 - nogammat);
  iE  = tkf->betat / (1.0 - tkf->betat) * ( (double)L + 1.0 - d0E );    
  
  nE = d1E + tkf->betat * (sE + 1.);
  eE = (d1E + iE) / (1.0 - tkf->R->etaz);

  if (ret_sE)  *ret_sE  = sE;
  if (ret_d0E) *ret_d0E = d0E;
  if (ret_d1E) *ret_d1E = d1E;
  if (ret_iE)  *ret_iE  = iE;
  if (ret_nE)  *ret_nE  = nE;
  if (ret_eE)  *ret_eE  = eE;
  if (ret_lE)  *ret_lE  = sE + eE;

  if (ret_gamma) *ret_gamma = tkf->gammat;
  if (ret_beta)  *ret_beta  = tkf->betat;
  if (ret_heta)  *ret_heta  = nogammat;
  if (ret_eta)   *ret_eta   = etat;
  
  tkf_model_Destroy(tkf); 
  return eslOK;

 ERROR:
  if (tkf) tkf_model_Destroy(tkf); 
  return status;
}


int 
tkf_sim_FiniteTime(ESL_RANDOMNESS *r, TKF_RATE *R, double time, int L, int N,  ESQ **esq, ESQ **newesq, double *sS, double *d0S, double *d1S, double *iS, 
		   double *nS, double *eS, double *lS, double tol, char *errbuf, int verbose)
{
  TKF_MODEL *tkf = NULL;
  ESQ       *Asq = NULL;
  ESQ       *Dsq = NULL;
  int        s;                  /* number of surviving residues */
  int        nb;                 /* total number of insert blocks */
  int        e;                  /* total number of inserted fragments */
  int        d0;
  int        d1;
  int        fr;                 /* total number of inserted residues */
  int        n;
  int        status;

  tkf = tkf_model_Create(R, time, 1.0, NULL, L, NULL, tol, verbose);
  if (tkf == NULL) ESL_XFAIL(eslFAIL, errbuf, "error creating tkfodel");

  for (n = 0; n < N; n++) {
     /* initialize 
      * from esq we only use the length. Asq and Dsq below are
      * fresh new sequences with just L alive ancestrals which are
      * evolved as a whole to t=time.
      */
    
    /* allocate */
    Asq = e1_sim_ESQCreate(esq[n]->L);
    Dsq = e1_sim_ESQCreate(esq[n]->L);

    tkf_finite_fate_residue(r, tkf, time, Asq, Dsq);
    tkf_finite_fate_insert (r, tkf, time, Asq, Dsq); 
 
    /* residue counting */
    e1_sim_ESQCount(Dsq, &s, &nb, &e, &d0, &d1, &fr);
    if (tkf->R->etaz == 0.0 && e != fr) { printf("TKF91 e=%d different from fr %d\n", e, fr); exit(1); }
    if (sS)  sS[n]  = (double)s;               /* ancestral residues still alive */
    if (d0S) d0S[n] = (double)d0;              /* number of d0 deletions */
    if (d1S) d1S[n] = (double)d1;              /* number of d1 deletions */
    if (iS)  iS[n]  = (double)e-(double)d1;    /* number of i residues/fragments */
    if (nS)  nS[n]  = (double)nb;              /* number of insert blocks */
    if (eS)  eS[n]  = (double)fr;              /* number of total inserted residues */
    if (lS)  lS[n]  = (double)s+(double)fr;    /* number of total descendent residues  */
        
    if (newesq) { e1_sim_ESQCopy(Dsq, newesq[n]); }
    e1_sim_ESQDestroy(Asq); Asq = NULL;
    e1_sim_ESQDestroy(Dsq); Dsq = NULL;
  }
 
  tkf_model_Destroy(tkf); 

  return eslOK;

 ERROR:  
  if (tkf) tkf_model_Destroy(tkf); 
  if (Asq) e1_sim_ESQDestroy(Asq);
  if (Dsq) e1_sim_ESQDestroy(Dsq);
  return status;
}


int
tkf_sim_Infinitesimal(ESL_RANDOMNESS *r, TKF_RATE *R, int N, ESQ **esq, double time, double tinc, double tepsilon, double *sS, double *d0S, double *d1S, double *iS, 
		      double *nS, double *eS, double *lS, int *ret_nbad, ESL_HISTOGRAM *h, double tol, char *errbuf, int verbose)
{
  ESQ      *Asq = NULL;            /* ancestral esq */
  ESQ      *Dsq = NULL;
  double    tt;
  int       s;                  /* number of surviving residues */
  int       nb;                 /* number of inserts */
  int       e;                  /* total number of inserted fragments */
  int       d0;
  int       d1;
  int       fr;                 /* total number of inserted residues */
  int       n;
  int       nbad = *ret_nbad;
  int       status;

  /* initialize */
  esl_vec_DSet(sS,  N, 0.0);
  esl_vec_DSet(d0S, N, 0.0);
  esl_vec_DSet(d1S, N, 0.0);
  esl_vec_DSet(iS,  N, 0.0);
  esl_vec_DSet(nS,  N, 0.0);
  esl_vec_DSet(eS,  N, 0.0);
  esl_vec_DSet(lS,  N, 0.0);
  
  for (n = 0; n < N; n++) {
    /* initialize */
    tt  = time + tepsilon;
    
    Asq = e1_sim_ESQClone(esq[n]);
    Dsq = e1_sim_ESQClone(esq[n]); 
    if (Asq == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of ancestral ESQ failed");
    if (Dsq == NULL) ESL_XFAIL(eslFAIL, errbuf, "allocation of descendant ESQ failed");
       
    e1_sim_ESQCount(Asq, &s, &nb, &e, &d0, &d1, &fr);
    //printf("^^N[%d]time %f time+tinc %f s %d n %d e %d d0 %d d1 %d fr %d\n", n, time, time+tinc, s, nb, e, d0, d1, fr);
    while (tt < time+tinc+tepsilon) {
      tkf_inf_fate_residue(r, R, tepsilon, Asq, Dsq);
      e1_sim_ESQCopy(Dsq, Asq);
      tkf_inf_fate_insert (r, R, tepsilon, Asq, Dsq, &nbad); 
      e1_sim_ESQCopy(Dsq, Asq);

      tt += tepsilon;
    }    

    /* residue counting */
    e1_sim_ESQCount(Dsq, &s, &nb, &e, &d0, &d1, &fr);
    if (R->etaz == 0.0 && e != fr) { printf("TKF91 e=%d fr =%d\n", e, fr); exit(1); }

    //printf("++N[%d]time %f s %d n %d e %d d0 %d d1 %d fr %d\n\n", n, tt-tepsilon, s, nb, e, d0, d1, fr);
    if (sS)  sS[n]  = (double)s;               /* ancestral residues still alive */
    if (d0S) d0S[n] = (double)d0;              /* number of d0 deletions */
    if (d1S) d1S[n] = (double)d1;              /* number of d1 deletions */
    if (iS)  iS[n]  = (double)e-(double)d1;    /* number of i residues/fragmenst */
    if (nS)  nS[n]  = (double)nb;              /* number of insert blocks */
    if (eS)  eS[n]  = (double)fr;              /* number of total inserted residues */
    if (lS)  lS[n]  = (double)s+(double)fr;    /* number of total descendent residues  */
  
   /* collate insert lengths */
    if (h) e1_sim_ESQInsertLenHisto(h, Dsq);

    /* return the descendant sequence */
    e1_sim_ESQCopy(Dsq, esq[n]);

    e1_sim_ESQDestroy(Asq); Asq = NULL;
    e1_sim_ESQDestroy(Dsq); Dsq = NULL;
  }

  *ret_nbad = nbad;

  return eslOK;

 ERROR:
  if (Dsq) e1_sim_ESQDestroy(Dsq);
  return status;
}


/* Ancestral residue fate:
 *
 *        T{A->0} = muD * te;        // ancestral residue dies
 *        T{A->A} = 1 - T{A->0};     // ancestral residue survives
 */
static int
tkf_finite_fate_residue(ESL_RANDOMNESS *r, TKF_MODEL *tkf, double tt, ESQ *Asq, ESQ *Dsq)
{
  double Tdead = tkf->gammat;
  int    i;

  for (i = 1; i <= Asq->L; i ++) {
    if (Asq->Alive[i]) {
      if (esl_random(r) < Tdead) Dsq->Alive[i] = FALSE;
    }
  }

  return eslOK;
} 
static int
tkf_inf_fate_residue(ESL_RANDOMNESS *r, TKF_RATE *R, double tepsilon, ESQ *Asq, ESQ *Dsq)
{
  double Tdead = R->mu * tepsilon;
  int    i;
  
  for (i = 1; i <= Asq->L; i ++) {    
    if (Asq->Alive[i]) {
      if (esl_random(r) < Tdead) Dsq->Alive[i] = FALSE;
    }
  }

  return eslOK;
} 

static int
tkf_finite_fate_insert(ESL_RANDOMNESS *r, TKF_MODEL *tkf, double tt, ESQ *Asq, ESQ *Dsq)
{
  double  nobeta;
  int     l;
  int     i;
  int     x;
  int     status;

  nobeta = (tkf->gammat > 0.0 && tkf->R->lambda > 0.0)? 1.0 - tkf->R->mu * tkf->betat / ( tkf->gammat * tkf->R->lambda) : 1.0;

  for (i = 0; i < Asq->L; i ++) {
    if (Asq->Alive[i] && Dsq->Alive[i]) { /* after a surviving ancestral residued, 
			  * sample with geometric of parameter beta_t 
			  * bn = 0 (1-beta) 
			  * bn = 1 (1-beta) * beta
			  * bn = 2 (1-beta) * beta^2
			  */
      l = (tkf->betat == 0.0)? 0 : (int)floor(log(esl_random(r)) / log(tkf->betat));
    }
    else if (Asq->Alive[i] && !Dsq->Alive[i]) { /* after a deletion, 
	    * sample with weird geometric:
	    * bn = 0 (1-nobeta)
	    * bn = 1 (1-beta) * nobeta
	    * bn = 2 (1-beta) * nobeta * beta
	    * bn = 3 (1-beta) * nobeta * beta^2
	    */
      l = 0;
      if (1.0 - nobeta < esl_random(r)) {
	l = (tkf->betat == 0.0)? 1 : (int)floor(log(esl_random(r)) / log(tkf->betat) + 1.0);
      }
    } 
    Dsq->Bn[i] = l;
    
    /* if TKF92 (ie T->etaz > 0), for each new l add the fragment */
    ESL_ALLOC(Dsq->Fn[i], sizeof(int) * Dsq->Bn[i]);
    for (x = 0; x < Dsq->Bn[i]; x ++) { 
      Dsq->Fn[i][x] = (tkf->R->etaz == 0.0)? 1 : (int)floor(log(esl_random(r)) / log(tkf->R->etaz) + 1.);
    }
  }
  
  return eslOK;

 ERROR:
  return status;
}

static int
tkf_inf_fate_insert(ESL_RANDOMNESS *r, TKF_RATE *R, double te, ESQ *Asq, ESQ *Dsq, int *ret_nbad)
{
  double  ins = te * R->lambda;
  double  del = te * R->mu;
  int     nbad     = *ret_nbad;
  int     newevent = 0;
  int     i;
  int     x;
  int     y;
  int     status;

  /* an insert can go in s+1+e places */
  for (i = 0; i <= Asq->L; i ++) {
    if (Asq->Alive[i]) {
      if (esl_random(r) < ins) { 
	Dsq->Bn[i] ++;

	/* if TKF92 */
	if (Dsq->Fn[i]) ESL_REALLOC(Dsq->Fn[i], sizeof(int) * Dsq->Bn[i]);
	else            ESL_ALLOC  (Dsq->Fn[i], sizeof(int) * Dsq->Bn[i]);
	Dsq->Fn[i][Dsq->Bn[i]-1] = (R->etaz == 0.0)? 1 : (int)floor(log(esl_random(r)) / log(R->etaz) + 1.);

	newevent ++; if (newevent > 1) nbad ++;
      }
    }
    
    for (x = 1; x <= Asq->Bn[i]; x ++) {
      if (esl_random(r) < ins) { 
	Dsq->Bn[i] ++; 

	/* if TKF92 */
	if (Dsq->Fn[i]) ESL_REALLOC(Dsq->Fn[i], sizeof(int) * Dsq->Bn[i]);
	else            ESL_ALLOC  (Dsq->Fn[i], sizeof(int) * Dsq->Bn[i]);
	Dsq->Fn[i][Dsq->Bn[i]-1] = (R->etaz == 0.0)? 1 : (int)floor(log(esl_random(r)) / log(R->etaz) + 1.);

	newevent ++; if (newevent > 1) nbad ++;
      }
    }
  }
 
   /* all Bn's can be deleted (the whole fragment) */
  for (i = 0; i <= Asq->L; i ++) {
    for (x = 1; x <= Asq->Bn[i]; x ++) {
      if (esl_random(r) < del) { 
	for (y = x + 1; y < Dsq->Bn[i]; y ++) {
	  Dsq->Fn[i][y-1] = Dsq->Fn[i][y];
	}
	
	Dsq->Bn[i] --; 
	newevent ++; if (newevent > 1) nbad ++; 
      }
    }
  }

  *ret_nbad = nbad;
  return eslOK;

 ERROR:
  return status;
}



/*****************************************************************
 * 
 *****************************************************************/
