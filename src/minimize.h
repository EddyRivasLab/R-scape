/* minimize
 *
 */
#ifndef MINIMIZE_INCLUDED
#define MINIMIZE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_dmatrix.h"

#include "lbfgs.h"

#define MAXITERATIONS 100

enum NMtransf_e {
  REFLECT   = 0,
  EXPAND    = 1,
  OCONTRACT = 2,
  ICONTRACT = 3,
  SHRINK    = 4,
  TNONE     = 5
};

extern int min_Bracket(double *x, double *u, long n, double firststep,
		       double (*func)(double *, long, void *),
		       void *prm, double tol, double *wrk, double *ret_fx);
extern int min_NelderMead(double *x, double *u, long n,
			  double (*func)(double *, long, void *),
			  void *prm, double tol, double *wrk, ESL_DMATRIX *simplex, double *ret_fx);
extern int min_LBFGS(int n,
		     lbfgsfloatval_t (evaluate)(void *, const lbfgsfloatval_t *, lbfgsfloatval_t *, const int, const lbfgsfloatval_t),
		     double tol, double *ret_fx);

#endif /*MINIMIZE_INCLUDED*/


/************************************************************
 * @LICENSE@
 ************************************************************/
