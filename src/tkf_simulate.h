/* tkf_simulate - funtions to evolve sequences either from
 *               the finite-time distribution or from the infinitesimal rate
 *
 */
#ifndef TKFSIMULATE_INCLUDED
#define TKFSIMULATE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_histogram.h"
#include "esl_ratematrix.h"

#include "e2.h"
#include "e1_simulate.h"
#include "ratematrix.h"


extern int  tkf_sim_Theoretical(TKF_RATE *R, double time, int L, double *ret_sE, double *ret_d0E, double *ret_d1E, double *ret_iE, 
				double *ret_nE, double *ret_eE, double *ret_lE, 
				double *ret_gamma, double *ret_beta, double *ret_heta, double *ret_eta, double tol, char *errbuf, int verbose);
extern int  tkf_sim_FiniteTime(ESL_RANDOMNESS *r, TKF_RATE *R, double time, int L, int N, ESQ **esq, ESQ **newesq, double *sS, double *d0S, double *d1S, double *iS, 
			       double *nS, double *eS, double *lS, double tol, char *errbuf, int verbose);
extern int  tkf_sim_Infinitesimal(ESL_RANDOMNESS *r, TKF_RATE *R1, int N, ESQ **esq, double time, double tinc, double tepsilon, double *sS, double *d0S, double *d1S, 
				  double *iS, double *nS, double *eS, double *lS, int *ret_nbad, ESL_HISTOGRAM *h, double tol, char *errbuf, int verbose);

#endif /*TKFSIMULATE_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
