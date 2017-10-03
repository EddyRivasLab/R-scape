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

extern int min_ConjugateGradientDescent(double *x, double *u, int n, 
					double (*func)(double *, int, void *),
					double (*bothfunc)(double *, int, void *, double *),
					void *prm, double tol, double *wrk, double *ret_fx);

extern int min_Bracket(double *x, double *u, int n, double firststep,
		       double (*func)(double *, int, void *),
		       void *prm, double tol, double *wrk, double *ret_fx);
extern int min_NelderMead(double *x, double *u, int n,
			  double (*func)(double *, int, void *),
			  void *prm, double tol, double *wrk, ESL_DMATRIX *simplex, double *ret_fx);

#endif /*MINIMIZE_INCLUDED*/


/************************************************************
 * @LICENSE@
 ************************************************************/
