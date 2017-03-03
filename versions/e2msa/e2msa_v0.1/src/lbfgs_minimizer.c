/* Multidimensional optimization using 
 * 
 * Can be used even without derivative information; falls back to
 * a numeric gradient if analytic gradient is unavailable.
 */
#include "esl_config.h"

#include <math.h>
#include <float.h>

#include "lbfgs.h"

#include "evopipeline.h"
#include "lbfgs_minimizer.h"

#include "easel.h"
#include "esl_vectorops.h"
#include "esl_minimizer.h"


/* Return the negative gradient at a point, determined 
 * numerically.
 */
static void
numeric_derivative(double *x, double *u, int n, 
		   double (*func)(double *, int, void*),
		   void *prm, double relstep,
		   double *dx)
{
  int    i;
  double delta;
  double f1, f2;
  double tmp;

  for (i = 0; i < n; i++)
    {
      delta = fabs(u[i] * relstep);

      tmp = x[i]; 
      x[i] = tmp + delta;
      f1  = (*func)(x, n, prm);
      x[i] = tmp - delta;
      f2  = (*func)(x, n, prm);
      x[i] = tmp;

      dx[i] = (-0.5 * (f1-f2)) / delta;

      ESL_DASSERT1((! isnan(dx[i])));
    }
}

/* function invoked by lbfgs()
 */
static lbfgsfloatval_t 
evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{ 
  struct lbfgs_data     *prm = (struct lbfgs_data *) instance;
  struct optimize_data  *data = prm->data;
  double               (*func)(double *, int, void *) = prm->func;
  void                 (*dfunc)(double *, int, void *, double *) = prm->dfunc;
  lbfgsfloatval_t        fx;
  int                    i;
  
  fx = (lbfgsfloatval_t)(*func)((double *)x, n, data);	/* the objective function */
 
  if (dfunc != NULL) 
    {
      (*dfunc)((double *)x, n, data, (double *)g);	/* find the current negative gradient, - df(x)/dxi  */
      esl_vec_DScale((double *)g, n, -1.0);
    } 
  else numeric_derivative((double *)x, prm->u, n, func, data, 1e-4, (double *)g); /* resort to brute force */
 
#if 1
  printf("\nEVALUATE step %f\n", step);
  printf("fx = %f\n", fx);
  printf("point   :    ");
  for (i = 0; i < n; i++)
    printf("%f ", x[i]);  
  printf("\ngradient:    ");
  for (i = 0; i < n; i++)
    printf("%f ", g[i]);
  printf("\n");
 #endif

  return fx;
}

/* function invoked by lbfgs()
 */
static int
progress(void *instance, 
	 const lbfgsfloatval_t *x,    const lbfgsfloatval_t *g,    const lbfgsfloatval_t fx, 
	 const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
	 int n, int it, int ls)
{
  int i;
  
  printf("Iteration %d:\n", it);
  printf("  fx = %f\n", fx);
  printf("new point   :    ");
  for (i = 0; i < n; i++)
    printf("%f ", x[i]);
  
  printf("\nnew gradient:    ");
  for (i = 0; i < n; i++)
    printf("%f ", g[i]);
  
  printf("\nxnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
  printf("\n");
  return 0;
}


/* Function:  min_LBFGS()
 * Incept:    ER, Mon Apr 22 09:22:46 EDT 2013 [janelia]
 *
 * Purpose:   n-dimensional minimization by L-BFGS.
 *           
 *            An initial point is provided by <x>, a vector of <n>
 *            components. The caller also provides a function <*func()> that 
 *            compute the objective function f(x) when called as 
 *            <(*func)(x, n, prm)>, and a function <*dfunc()> that can
 *            compute the gradient <dx> at <x> when called as 
 *            <(*dfunc)(x, n, prm, dx)>, given an allocated vector <dx>
 *            to put the derivative in. Any additional data or fixed
 *            parameters that these functions require are passed by
 *            the void pointer <prm>.
 *            
 *            The caller also provides an allocated workspace sufficient to
 *            hold four allocated n-vectors. (4 * sizeof(double) * n).
 *
 *            Iterations continue until the objective function has changed
 *            by less than a fraction <tol>. This should not be set to less than
 *            sqrt(<DBL_EPSILON>). 
 *
 *            Upon return, <x> is the minimum, and <ret_fx> is f(x),
 *            the function value at <x>.
 *            
 * Args:      x        - an initial guess n-vector; RETURN: x at the minimum
 *            u        - "units": maximum initial step size along gradient when bracketing.
 *            n        - dimensionality of all vectors
 *            *func()  - function for computing objective function f(x)
 *            *dfunc() - function for computing a gradient at x
 *            prm      - void ptr to any data/params func,dfunc need 
 *            tol      - convergence criterion applied to f(x)
 *            wrk      - allocated 4xn-vector for workspace
 *            ret_fx   - optRETURN: f(x) at the minimum
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOHALT> if it fails to converge in MAXITERATIONS.
 *            <eslERANGE> if the minimum is not finite, which may
 *            indicate a problem in the implementation or choice of <*func()>.
 *
 * Xref:      STL9/101.
 */
int
min_LBFGS(double *x, double *u, int n, 
	  double (*func)(double *, int, void *),
	  void (*dfunc)(double *, int, void *, double *),
	  void *prm, double tol, double *wrk, double *ret_fx)
{
  lbfgs_parameter_t      param;
  lbfgsfloatval_t        fx;
  struct lbfgs_data     *instance = NULL;
  struct optimize_data  *data = (struct optimize_data *)prm;
  double                 oldfx;
  double                *dx;
  int                    status;
  
  dx = wrk;
  oldfx = (*func)(x, n, prm);	/* init the objective function */
  
   /* Bail out if the function is +/-inf: this can happen if the caller
   * has screwed something up, or has chosen a bad start point.
   */
  if (oldfx == eslINFINITY || oldfx == -eslINFINITY)
	  ESL_EXCEPTION(eslERANGE, "minimum not finite");

  if (dfunc != NULL) 
    {
      (*dfunc)(x, n, prm, dx);	/* find the current negative gradient, - df(x)/dxi  */
      esl_vec_DScale(dx, n, -1.0);
    } 
  else numeric_derivative(x, u, n, func, prm, 1e-4, dx); /* resort to brute force */

#  /* data */
  ESL_ALLOC(instance, sizeof(struct lbfgs_data)); if (instance == NULL) return eslFAIL;
  instance->u     = u;
  instance->data  = data;
  instance->func  = func;
  instance->dfunc = dfunc;
  
   /* Initialize the parameters for the L-BFGS optimization
   */
  lbfgs_parameter_init(&param);
  /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
  param.max_iterations = MAXITERATIONS;
  param.epsilon = tol;
  param.max_step = 1e+50;
  param.min_step = 1e-50;
  
  /* Start the L-BFGS optimization; this will invoke the callback functions
   * evaluate() and progress() when necessary.
  */
  status = lbfgs(n, x, &fx, evaluate, progress, (void *)instance, &param);
  if (status < 0) { printf("L-BFGS optimization failed with status %d | LBFGS_ALREADY_MINIMIZED %d \n", status, LBFGSERR_ROUNDING_ERROR); status = eslFAIL; goto ERROR; }
   
  if (ret_fx != NULL) *ret_fx = (double)fx;
  
  if (instance) free(instance);
  return eslOK;

 ERROR:
  if (instance) free(instance);
  return status;
}




/*****************************************************************
 * Example main()
 *****************************************************************/
#ifdef LBFGS_MINIMIZER_EXAMPLE
/*::cexcerpt::LBFGS_minimizer_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DLBFGS_MINIMIZER_EXAMPLE lbfgs_minimizer.c esl_vectorops.c easel.c -lm
 * run:     ./example 
 */
#include <stdio.h>
#include <lbfgs.h>

#include "easel.h"
#include "esl_vectorops.h"
#include "esl_minimizer.h"

/* a simple multidimensional quadratic w/ a minimum at 0:
 *    $f(x) = a_1 x_1^2 + ... a_n x_n^2$
 */ 
static double
example_func(double *x, int n, void *prm)
{
  double *a;
  double  fx;
  int     i;

  a = (double *) prm;	/* cast the data vector */
  for (fx = 0., i = 0; i < n; i++)
    fx += a[i] * x[i] * x[i];
  return fx;
}
/* gradient of the f(x): d/dx_i = 2 a_i x_i
 */
static void
example_dfunc(double *x, int n, void *prm, double *dx)
{
  double *a;
  int     i;

  a = (double *) prm;	/* cast the data vector */
  for (i = 0; i < n; i++)
    dx[i] = 2.0 * a[i] * x[i];
}
int
main(int argc, char **argv)
{
  int    n = 6;
  double a[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  double x[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  double u[6] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  double fx;
  int    i;

  min_LBFGS(x, u, n, 
	    &example_func, &example_dfunc, (void *) a, 
	    0.0001, &fx);
  
  printf("At minimum: f(x) = %g\n", fx);
  printf("vector x = ");
  for (i = 0; i < 6; i++) printf("%g  ", x[i]);
  printf("\n");

  return 0;
}
/*::cexcerpt::minimizer_example::end::*/
#endif /*LBFGS_MINIMIZER_EXAMPLE*/






/*****************************************************************  
 * @LICENSE@
 * 
 * SVN $$
 * SVN $URL:$
 *****************************************************************/
