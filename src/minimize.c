/*  Minimize
 *
 * ER, Thu Mar  6 16:40:59 EST 2014 [Janelia] 
 * SVN $Id:$
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_vectorops.h"

#include "minimize.h"

static void   numeric_derivative(double *x, double *u, int n, double (*func)(double *, int, void*),
				 void *prm, double relstep, double *gx);
static double cubic_interpolation(double xa, double fa, double ga, double xb, double fb, double gb, double xmin, double xmax);
static int    Armijo(double *ori, double fori, double *gori, double *dori, int n, double firststep, double c1,
		     double (*bothfunc)(double *, int, void *, double *), void *prm,
		     double *x, double *g, double *ret_f, double *ret_dg, double *ret_step, double tol);
static int    Wolfe(double *ori, double fori, double *gori, double *dori, int n, double firststep, double c1, double c2,
		    double (*bothfunc)(double *, int, void *, double *), void *prm,
		    double *x, double *ret_step, double *ret_f, double *g, double tol);

/* Return the gradient at a point, determined numerically.
 */
static void
numeric_derivative(double *x, double *u, int n, 
		   double (*func)(double *, int, void*),
		   void *prm, double relstep,
		   double *gx)
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

      gx[i] = (-0.5 * (f1-f2)) / delta;

      ESL_DASSERT1((! isnan(gx[i])));
    }
}

/* cubic_interpolation():
 *
 * ER, Thu Oct 26 10:04:14 EDT 2017 [Cambridge]
 *
 * Purpose: Given two points with there fun values and derivatives,
 *          calculate the middle point by cubic interpolation.
 *
 *
 *
 * Args:    xa -  point
 *          fa -  function at xa
 *          ga -  gradient at xa
 *          xb -  point xb > xa (or swap)
 *          fb -  function at xb
 *          gb -  gradient at xb
 *        xmin -  [xim,xmax] interval
 *        xmax -  [xim,xmax] interval
 *
 * Returns:   xm - middle point 
 *
 */
static double cubic_interpolation(double xa, double fa, double ga, double xb, double fb, double gb, double xmin, double xmax)
{
  double da, db;
  double xc;
  double xm;       // the mid point
  double swapper;
  
  // we assume xb > xa
  if (xa == xb) ESL_EXCEPTION(eslENORESULT, "cubic interpolation(): xa has to be different from xb");
  if (xa >  xb)
    {
      swapper = xa; xa = xb; xb = swapper;
      swapper = fa; fa = fb; fb = swapper;
      swapper = ga; ga = gb; gb = swapper;
    }

  da = ga + gb - 3.*(fa-fb)/(xa-xb);
  db = da*da - ga*gb;
  
  if (db >= 0.) {
    db = sqrt(db);
    xc = xb - (xb-xa) * (gb + db - da) / (gb - ga + 2.*db);
    xm = ESL_MIN(ESL_MAX(xc,xmin),xmax);
  }
  else
    xm = 0.5 * (xmin + xmax);
  
  return xm;
}


/* Armijo():
 * ER, Wed Oct 25 23:42:41 EDT 2017 [Cambridge]
 *
 * Purpose:   Backtracking linesearch to satisfy Armijo condition.
 *            Derived from ArmijoBacktracking.m by Mark Schmidt
 *
 *
 *
 * Args:      ori         - original n-vector 
 *            fori        - f(ori)
 *            gori        - the gradien of f(x) at ori
 *            dori        - direction vector we're following from ori
 *            n           - dimensionality of x, g, and d
 *            firststep   - step (t) is initialized to this (positive) value
 *            c1          - coefficient for the Armijo condition
 *            c2          - coefficient for the Wolrf condition
 *            (*bothfunc) - ptr to caller's objective function
 *            prm         - ptr to any additional data (*func)() needs
 *            x           - RETURN: minimum, as an n-vector (caller allocated)
 *            g           - RETURN: gradient at x (caller allocated)
 *            ret_f       - optRETURN: the function at x
 *            ret_t       - optRETURN: the step size
 *
 * Returns:   <eslOK> on success.
 *
 * Reference: 
 */
static int Armijo(double *ori, double fori, double *gori, double *dori, int n,
		  double firststep, double c1,
		  double (*bothfunc)(double *, int, void *, double *), void *prm,
		  double *x, double *g, double *ret_f, double *ret_dg, double *ret_step, double tol)
{
  double dgori = esl_vec_DDot(dori, gori, n);  // initial d'*g
  double f;
  double dg;
  double t, t_prv;
  double min_step = 1e-8;
  double max_step = 0.6;
  int    nit = 0;
  int    status;
  
  // Check inputs 
  if (firststep <= 0.) ESL_EXCEPTION(eslENORESULT, "Step size is negative");
  if (dgori      > 0.) ESL_EXCEPTION(eslENORESULT, "Not a descent direction");
  
  // Calculate the first new point
  t_prv = 0.; 
  t     = firststep; 
  esl_vec_DCopy(ori, n, x);
  esl_vec_DAddScaled(x, dori, t, n);
  f = (*bothfunc)(x, n, prm, g);
  
  // Do until the Armijo condition is satisfied
  while (f > fori + c1 * t * dgori) {
    nit ++;
    
    // new dg
    dg = esl_vec_DDot(dori, g, n);
    
    // calculate a new step by cubic interpolation
    t = cubic_interpolation(0, fori, dgori, t, f, dg, 0, t);
    if (t < t_prv*min_step) t = t_prv*min_step;
    if (t > t_prv*max_step) t = t_prv*max_step;
    
    // calculate the new point
    esl_vec_DCopy(ori, n, x);
    esl_vec_DAddScaled(x, dori, t, n);
    f = (*bothfunc)(x, n, prm, g);
    
    //if (nit > MAXITER) printf("Armijo() reached is the max number of iterations\n");
    
    t_prv = t;
  }

  if (ret_f)    *ret_f    = f;
  if (ret_dg)   *ret_dg   = dg;
  if (ret_step) *ret_step = t;
  
  return eslOK;
}


/* Wolfe():
 * ER, Wed Oct 25 23:42:41 EDT 2017 [Cambridge]
 *
 * Purpose:   Backtracking linesearch to satisfy Armijo condition.
 *            Derived from WolfeLineSearch.m by Mark Schmidt
 *
 *
 *
 * Args:      ori         - original n-vector 
 *            fori        - f(ori)
 *            g           - the gradien of f(x) 
 *            dori        - direction vector we're following from ori
 *            n           - dimensionality of x, g, and d
 *            firststep   - step (t) is initialized to this (positive) value
 *            c1          - coefficient for the Armijo condition
 *            c2          - coefficient for the Wolrf condition
 *            (*bothfunc) - ptr to caller's objective function
 *            prm         - ptr to any additional data (*func)() needs
 *            x           - RETURN: minimum, as an n-vector (caller allocated)
 *            g           - RETURN: gradient at x (caller allocated)
 *            ret_f       - optRETURN: the function at x
 *            ret_t       - optRETURN: the step size
 *
 * Returns:   <eslOK> on success.
 *
 * Reference: 
 */
static int Wolfe(double *ori, double fori, double *gori, double *dori, int n,
		 double firststep, double c1, double c2,
		 double (*bothfunc)(double *, int, void *, double *), void *prm,
		 double *x, double *ret_step, double *ret_fx, double *g, double tol)
{
  double dgori = esl_vec_DDot(dori, gori, n);  // initial d'*g
  double f, f_prv;
  double dg, dg_prv;
  double t, t_prv, t_new;
  double min_step;
  double max_step;
  double ta, tb;
  double fa, fb;
  double dga, dgb;
  double tmax, tmin;
  int    nit = 0;
  int    found = FALSE;
  int    status;

  // Check inputs 
  if (firststep <= 0.) ESL_EXCEPTION(eslENORESULT, "Step size is negative");
  if (dgori      > 0.) ESL_EXCEPTION(eslENORESULT, "Not a descent direction");

  // init (t_prv,f_prv,g_prv)
  t_prv  = 0.;
  f_prv  = fori;
  dg_prv = dgori;
  
  // The first new point (t,f,g)
  // using x = xori + firststep*dori
  t  = firststep; 
  esl_vec_DCopy(ori, n, x);
  esl_vec_DAddScaled(x, dori, t, n);
  f  = (*bothfunc)(x, n, prm, g);
  dg = esl_vec_DDot(dori, g, n);

  while (nit < MAXITER) {

    if (f > fori + c1*t*dgori || (nit > 0 && f >= f_prv)) // Armijo not satisfied 
      break;
    else if  (fabs(dg) <= -c2*dgori) {                    // Armijo + strong_Wolfe satisfied, you are done
      found = TRUE;
      break;
    }
    else if (dg > 0.)
      break;

    if (t-t_prv < 1e-6) break; // not enough progress
    
    // we are still here (have not bailed out with either a solution or a bracket, then
    //
    // calculate a new step (t_new) by cubic interpolation between (t_prv,f_prv,g_prv) and (t,f,g)
    min_step = t + 0.01 * (t-t_prv);
    max_step = t * 10;
    t_new = cubic_interpolation(t_prv, f_prv, dg_prv, t, f, dg, min_step, max_step);

    // test we are making enough progress
    if (ESL_MIN(tmax-t,t-tmin) / (tmax-tmin) < 0.1) {
      if (fabs(tmax-t) < fabs(t-tmin)) t = tmax - 0.1*(tmax-tmin);
      else                             t = tmin + 0.1*(tmax-tmin);
    }

    // (t,f,g) becomes (t_prv,f_prv,g_prv)
    t_prv  = t;
    f_prv  = f;
    dg_prv = dg;
    
    // calculate the new point (t=t_new,f,g)
    // at x = xori + t_new * dori
    t  = t_new;
    esl_vec_DCopy(ori, n, x);
    esl_vec_DAddScaled(x, dori, t, n);
    f  = (*bothfunc)(x, n, prm, g);
    dg = esl_vec_DDot(dori, g, n);
    
    nit ++;
  }
 
  // now we either have a solution (found = TRUE)
  // or a bracket (ta,fa,ga) (tb,fb,gb)
  // such that fa < fb
  if (f < f_prv) {
    ta  = t;  tb  = t_prv;
    fa  = f;  fb  = f_prv;
    dga = dg; dgb = dg_prv;
  }
  else {
    tb  = t;  ta  = t_prv;
    fb  = f;  fa  = f_prv;
    dgb = dg; dga = dg_prv;
  }

  // refine the bracket
  while (found == FALSE && nit < MAXITER) {
    
    // calculate a new step (t) by cubic interpolation
    tmax = ESL_MAX(ta,tb);
    tmin = ESL_MIN(ta,tb);

    t = cubic_interpolation(ta, fa, dga, tb, fb, dgb, tmin, tmax);

    // test we are making enough progress
    if (ESL_MIN(tmax-t,t-tmin) / (tmax-tmin) < 0.1) {
      if (fabs(tmax-t) < fabs(t-tmin)) t = tmax - 0.1*(tmax-tmin);
      else                             t = tmin + 0.1*(tmax-tmin);
    }
  
    // calculate the new point (t,f,g)
    // at x = xori + t*dori
    esl_vec_DCopy(ori, n, x);
    esl_vec_DAddScaled(x, dori, t, n);
    f  = (*bothfunc)(x, n, prm, g);
    dg = esl_vec_DDot(dori, g, n);
    
    if (fabs(ta-tb)*dg < 1e-6) break; // no more progress is possible
	 
    if (f > fori + c1*t*dgori ||   // Armijo not satisfied 
	(nit > 0 && f >= fa)    )  // or new f is not lowest
      {
	tb  = t;
	fb  = f;
	dgb = dg;
      }
    else {
      if (fabs(dg) <= - c2*dgori) {  // strong Wolfe satisfied - we are done
	found = TRUE;
      }
      else if (dg*fabs(tb-ta) >= 0) {
	// old fb becomes new fa
	tb  = ta;
	fb  = fa;
	dgb = dga;
      }
         
      // new point becomes fa
      ta  = t;
      fa  = f;
      dga = dg;
    }

    nit ++;
  }
  //if (nit == MAXITER) printf("Wolfe() reached the max number of iterations\n");

  if (ret_fx)    *ret_fx  = fa;
  if (ret_step) *ret_step = ta;
  
  return eslOK;
}


/* Function:  min_ConjugateGradientDescent()
 * Incept:    ER, Fri Nov 10 10:18:00 EST 2017 [Cambridge]
 *
 * Purpose:   n-dimensional minimization by conjugate gradient descent.
 *           
 *            An initial point is provided by <x>, a vector of <n>
 *            components. The caller also provides a function <*func()> that 
 *            compute the objective function f(x) when called as 
 *            <(*func)(x, n, prm)>, and a function <*dfunc()> that can
 *            compute the gradient <gx> at <x> when called as 
 *            <(*dfunc)(x, n, prm, gx)>, given an allocated vector <gx>
 *            to put the derivative in. Any additional data or fixed
 *            parameters that these functions require are passed by
 *            the void pointer <prm>.
 *            
 *            The first step of each iteration is to try to bracket
 *            the minimum along the current direction. The initial step
 *            size is controlled by <u[]>; the first step will not exceed 
 *            <u[i]> for any dimension <i>. (You can think of <u> as
 *            being the natural "units" to use along a graph axis, if
 *            you were plotting the objective function.)
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
 *            stol     - convergence criterion applied to the line search
 *            wrk      - allocated 4xn-vector for workspace
 *            ret_fx   - optRETURN: f(x) at the minimum
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOHALT> if it fails to converge in MAXITER.
 *            <eslERANGE> if the minimum is not finite, which may
 *            indicate a problem in the implementation or choice of <*func()>.
 *
 * Xref:      STL9/101.
 */
int
min_ConjugateGradientDescent(double *x, double *u, int n, 
			     double (*func)(double *, int, void *),
			     double (*bothfunc)(double *, int, void *, double *),
			     void *prm, double tol, double stol, double *wrk, double *ret_fx)
{
  double oldfx;
  double coeff;
  double num, den;
  int    nit;
  int    i;
  double *gx, *cg, *w1, *w2;
  double cvg;
  double fa,fb,fc;
  double ax,bx,cx;
  double c1, c2;
  double fx;
  double gtd;
  double firststep;
  double t;        // step after Worfe algorithm
  double sum;

  gx = wrk;
  cg = wrk + n;
  w1 = wrk + 2*n;
  w2 = wrk + 3*n;

  /* init the objective function */
  if (bothfunc == NULL) 
    oldfx = (*func)(x, n, prm);	
  else
    oldfx = (*bothfunc)(x, n, prm, gx);	
  
  /* Bail out if the function is +/-inf or nan: this can happen if the caller
   * has screwed something up, or has chosen a bad start point.
   */
  if (! isfinite(oldfx)) ESL_EXCEPTION(eslERANGE, "minimum not finite");

  if (bothfunc == NULL) 
    numeric_derivative(x, u, n, func, prm, 1e-4, gx); /* resort to brute force */

  // First x, fx, gx, cg
  esl_vec_DCopy(gx, n, cg);  /* and make that the first conjugate direction, cg = -gx  */
  esl_vec_DScale(cg, n, -1.0);

  /* (failsafe) convergence test: a zero direction can happen, 
   * and it either means we're stuck or we're finished (most likely stuck)
   */
  for (i = 0; i < n; i++) 
    if (cg[i] != 0.) break;
  if  (i == n) {
    if (ret_fx != NULL) *ret_fx = oldfx;
    return eslOK;
  }
  
  for (nit = 0; nit < MAXITER; nit++)
    {
      
      // Initial step size
      if (nit == 0) {
	sum = 0.;
	for (i = 0; i < n; i ++)
	  sum += fabs(gx[i]);
	firststep = ESL_MIN(1.0, ((sum > 0.)? 1./sum:1.0) );
      }
      else {
	gtd = esl_vec_DDot(gx, cg, n);
	if (gtd > -tol) break;  // Check this is a good direction

	firststep = ESL_MIN(1.0, 2.*(fx-oldfx)/gtd);
	oldfx = fx;
      }

      // Strong Wolfe condition
      //
      // the new param are temporarily in w2
      // the negative gradient is temporarily in w1
      //
      c1 = 1e-4;  // parameter values in minFunc.m by Mark Schmidt
      c2 = 0.2;   
      Wolfe(x, oldfx, gx, cg, n, firststep, c1, c2, bothfunc, prm, w2, &t, &fx, w1, tol);
      esl_vec_DCopy(w2, n, x); //new parameters
      
      /* Main convergence test. 1e-10 factor is fudging the case where our
       * minimum is at exactly f()=0.
       */
      cvg = 2.0 * fabs((oldfx-fx)) / (1e-10 + fabs(oldfx) + fabs(fx));
      if (cvg <= tol) break;
      
      //fprintf(stdout, "(%d): Old f() = %.9f    New f() = %.9f    Convergence = %.9f\n", nit+1, oldfx, fx, cvg);

      if (nit == MAXITER-1) continue;

      /* Calculate the Hestenes-Stiefel coefficient */
      for (num = 0., i = 0; i < n; i++)
	num += (w1[i] - gx[i]) * w1[i];
      for (den = 0., i = 0; i < n; i++)
	den += (w1[i] - gx[i]) * cg[i];
      coeff = (den != 0)? num/den:0. ;
     
      /* Calculate the next conjugate gradient direction in w2 = -gx + coeff*cg */
      esl_vec_DCopy(w1, n, w2);
      esl_vec_DScale(w2, n, -1.0);
      esl_vec_DAddScaled(w2, cg, coeff, n);
      
      /* Finishing set up for next iteration: */
      esl_vec_DCopy(w1, n, gx);
      esl_vec_DCopy(w2, n, cg);
      
      /* Now: x is the current point; 
       *      fx is the function value at that point;
       *      gx is the current gradient at x;
       *      cg is the current conjugate gradient direction. 
       */

#if eslDEBUGLEVEL >= 2
      printf("\nesl_min_ConjugateGradientDescent():\n");
      printf("new point:     ");
      for (i = 0; i < n; i++)
	printf("%g ", x[i]);
      
      printf("\nnew gradient:    ");
      for (i = 0; i < n; i++)
	printf("%g ", gx[i]);
      
      numeric_derivative(x, u, n, func, prm, 1e-4, w1);
      printf("\n(numeric grad):  ");
      for (i = 0; i < n; i++)
	printf("%g ", w1[i]);
      
      printf("\nnew direction: ");
      for (i = 0; i < n; i++)
	printf("%g ", cg[i]);
      
      printf("\nOld f() = %g    New f() = %g    Convergence = %g\n\n", oldfx, fx, cvg);
#endif
      
      /* Second (failsafe) convergence test: a zero direction can happen, 
       * and it either means we're stuck or we're finished (most likely stuck)
       */
      for (i = 0; i < n; i++) 
	if (cg[i] != 0.) break;
      if  (i == n) break;

    }

  if (ret_fx != NULL) *ret_fx = fx;
  
  //if (nit == MAXITER) printf("min_ConjugateGradientDescent() reached the max number of iterations %d\n", nit);
  
  return eslOK;
}



