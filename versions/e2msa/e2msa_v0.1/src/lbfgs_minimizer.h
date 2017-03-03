/* ratematrix - funtions to calculate a rate matrix from conditional probabilities
 *
 */
#ifndef LBFGSMINIMIZER_INCLUDED
#define LBFGSMINIMIZER_INCLUDED

#include "evopipeline.h"

struct lbfgs_data {
  struct optimize_data  *data;
  double                *u;
  double               (*func)(double *, int, void *); 
  void                 (*dfunc)(double *, int, void *, double *);
};

extern int min_LBFGS(double *x, double *u, int n, 
		     double (*func)(double *, int, void *),
		     void (*dfunc)(double *, int, void *, double *),
		     void *prm, double tol, double *wrk, double *ret_fx);

#endif /* LBFGSMINIMIZER_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
