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

#include <lbfgs.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_vectorops.h"

#include "minimize.h"

/* bracket():
 * SRE, Wed Jul 27 11:43:32 2005 [St. Louis]
 *
 * Purpose:   Bracket a minimum. 
 *
 *            The minimization is quasi-one-dimensional, 
 *            starting from an initial <n>-dimension vector <ori>
 *            in the <n>-dimensional direction <d>.
 *            
 *            Caller passes a ptr to the objective function <*func()>,
 *            and a void pointer to any necessary conditional 
 *            parameters <prm>. The objective function will
 *            be evaluated at a point <x> by calling
 *            <(*func)(x, n, prm)>. The caller's function
 *            is responsible to casting <prm> to whatever it's
 *            supposed to be, which might be a ptr to a structure,
 *            for example; typically, for a parameter optimization
 *            problem, this holds the observed data.
 *            
 *            The routine works in scalar multipliers relative
 *            to origin <ori> and direction <d>; that is, a new <n>-dimensional
 *            point <b> is defined as <ori> + <bx><d>, for a scalar <bx>.
 *            
 *            The routine identifies a triplet <ax>, <bx>, <cx> such
 *            that $a < b < c$ and such that a minimum is known to
 *            exist in the $(a,b)$ interval because $f(b) < f(a),
 *            f(c)$. Also, the <a..b> and <b...c> intervals are in
 *            a golden ratio; the <b..c> interval is 1.618 times larger
 *            than <a..b>.
 *
 *            Since <d> is usually in the direction of the gradient,
 *            the points <ax>,<bx>,<cx> might be expected to be $\geq 0$;
 *            however, when <ori> is already close to the minimum, 
 *            it is often faster to bracket the minimum using
 *            a negative <ax>. The caller might then try to be "clever"
 *            and assume that the minimum is in the <bx..cx> interval
 *            when <ax> is negative, rather than the full <ax..cx>
 *            interval. That cleverness can fail, though, if <ori>
 *            is already in fact the minimum, because the line minimizer
 *            in brent() assumes a non-inclusive interval. Use
 *            <ax..cx> as the bracket.
 *            
 * Args:      ori       - n-dimensional starting vector
 *            d         - n-dimensional direction to minimize along
 *            n         - # of dimensions
 *            firststep - bx is initialized to this scalar multiplier
 *            *func()   - objective function to minimize
 *            prm       - void * to any constant data that *func() needs
 *            wrk       - workspace: 1 allocated n-dimensional vector
 *            ret_ax    - RETURN:  ax < bx < cx scalar bracketing triplet
 *            ret_bx    - RETURN:    ...ax may be negative
 *            ret_cx    - RETURN:    
 *            ret_fa    - RETURN:  function evaluated at a,b,c
 *            ret_fb    - RETURN:    ... f(b) < f(a),f(c)
 *            ret_fc    - RETURN:
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOHALT> if it fails to converge.
 *
 * Xref:      STL9/130.
 */
static int
bracket(double *ori, double *d, long n, double firststep,
	double (*func)(double *, long, void *), void *prm, 
	double *wrk, 
	double *ret_ax, double *ret_bx, double *ret_cx,
	double *ret_fa, double *ret_fb, double *ret_fc)
{
  double ax,bx,cx;		/* scalar multipliers */
  double fa,fb,fc;		/* f() evaluations at those points */
  double swapper;
  int    niter;
  
  /* Set and evaluate our first two points f(a) and f(b), which
   * are initially at 0.0 and <firststep>.
   */
  ax = 0.;  /* always start w/ ax at the origin, ax=0 */
  fa = (*func)(ori, n, prm);

  bx = firststep;
  esl_vec_DCopy(ori, n, wrk);
  esl_vec_DAddScaled(wrk, d, bx, n);
  fb = (*func)(wrk, n, prm);

  /* In principle, we usually know that the minimum m lies to the
   * right of a, m>=a, because d is likely to be a gradient.  You
   * might think we want 0 = a < b < c.  In practice, there's problems
   * with that. It's far easier to identify bad points (f(x) > f(a))
   * than to identify good points (f(x) < f(a)), because letting f(x)
   * blow up to infinity is fine as far as bracketing is concerned.
   * It can be almost as hard to identify a point b that f(b) < f(a)
   * as it is to find the minimum in the first place!
   * Counterintuitively, in cases where f(b)>f(a), it's better
   * to just swap the a,b labels and look for c on the wrong side
   * of a! This often works immediately, if f(a) was reasonably
   * close to the minimum and f(b) and f(c) are both terrible.
   */
  if (fb > fa)
    {
      swapper = ax; ax = bx; bx = swapper;
      swapper = fa; fa = fb; fb = swapper;
    }

  /* Make our first guess at c.
   * Remember, we don't know that b>a any more, and c might go negative.
   * We'll either have:      a..b...c with a=0;
   *                or:  c...b..a     with b=0.
   * In many cases, we'll immediately be done.
   */
  cx = bx + (bx-ax)*1.618;
  esl_vec_DCopy(ori, n, wrk);
  esl_vec_DAddScaled(wrk, d, cx, n);
  fc = (*func)(wrk, n, prm);
  
  /* We're not satisfied until fb < fa, fc; 
   * throughout the routine, we guarantee that fb < fa;
   * so we just check fc.
   */
  niter = 0;
  while (fc <= fb)
    {
      /* Slide over, discarding the a point; choose 
       * new c point even further away.
       */
      ax = bx; bx = cx;
      fa = fb; fb = fc;
      cx = bx+(bx-ax)*1.618;
      esl_vec_DCopy(ori, n, wrk);
      esl_vec_DAddScaled(wrk, d, cx, n);
      fc = (*func)(wrk, n, prm);

      /* This is a rare instance. We've reach the minimum
       * by trying to bracket it. Also check that not all
       * three points are the same.
       */
      if (ax != bx && bx != cx && fa == fb && fb == fc) break;

      niter++;
      if (niter > 100)
    	  ESL_EXCEPTION(eslENORESULT, "Failed to bracket a minimum.");
    }

  /* We're about to return. Assure the caller that the points
   * are in order a < b < c, not the other way.
   */
  if (ax > cx)
    {
      swapper = ax; ax = cx; cx = swapper;
      swapper = fa; fa = fc; fc = swapper;
    }

  /* Return.
   */
  ESL_DPRINTF2(("\nbracket(): %d iterations\n", niter));
  ESL_DPRINTF2(("bracket(): triplet is %g  %g  %g along current direction\n", 
		ax, bx, cx));
  ESL_DPRINTF2(("bracket(): f()'s there are: %g  %g  %g\n\n", 
		fa, fb, fc));

  *ret_ax = ax;  *ret_bx = bx;  *ret_cx = cx;
  *ret_fa = fa;  *ret_fb = fb;  *ret_fc = fc;
  return eslOK;
}


/* brent():
 * SRE, Sun Jul 10 19:07:05 2005 [St. Louis]
 *
 * Purpose:   Quasi-one-dimensional minimization of a function <*func()>
 *            in <n>-dimensions, along vector <dir> starting from a
 *            point <ori>. Identifies a scalar $x$ that approximates
 *            the position of the minimum along this direction, in a
 *            given bracketing interval (<a,b>).  The minimum must
 *            have been bracketed by the caller in the <(a,b)>
 *            interval.  <a> is often 0, because we often start at the
 *            <ori>.
 *
 *            A quasi-1D scalar coordinate $x$ (such as <a> or <b>) is
 *            transformed to a point $\mathbf{p}$ in n-space as:
 *            $\mathbf{p} = \mathbf{\mbox{ori}} + x
 *            \mathbf{\mbox{dir}}$.
 *
 *            Any extra (fixed) data needed to calculate <func> can be
 *            passed through the void <prm> pointer.
 *
 *            <eps> and <t> define the relative convergence tolerance,
 *            $\mbox{tol} = \mbox{eps} |x| + t$. <eps> should not be
 *            less than the square root of the machine precision.  The
 *            <DBL_EPSILON> is 2.2e-16 on many machines with 64-bit
 *            doubles, so <eps> is on the order of 1e-8 or more. <t>
 *            is a yet smaller number, used to avoid nonconvergence in
 *            the pathological case $x=0$.
 *
 *            Upon convergence (which is guaranteed), returns <xvec>,
 *            the n-dimensional minimum. Optionally, will also return
 *            <ret_x>, the scalar <x> that resulted in that
 *            n-dimensional minimum, and <ret_fx>, the objective
 *            function <*func(x)> at the minimum.
 *
 *            This is an implementation of the R.P. Brent (1973)
 *            algorithm for one-dimensional minimization without
 *            derivatives (modified from Brent's ALGOL60 code). Uses a
 *            combination of bisection search and parabolic
 *            interpolation; should exhibit superlinear convergence in
 *            most functions.
 *
 *
 * Args:      ori     - n-vector at origin
 *            dir     - direction vector (gradient) we're following from ori
 *            n       - dimensionality of ori, dir, and xvec
 *            (*func) - ptr to caller's objective function
 *            prm     - ptr to any additional data (*func)() needs
 *            a,b     - minimum is bracketed on interval [a,b]
 *            eps     - tol = eps |x| + t; eps >= 2 * relative machine precision
 *            t       - additional factor for tol to avoid x=0 case.
 *            xvec    - RETURN: minimum, as an n-vector (caller allocated)
 *            ret_x   - optRETURN: scalar multiplier that gave xvec
 *            ret_fx  - optRETURN: f(x)
 *
 * Returns:   (void)
 *
 * Reference: See [Brent73], Chapter 5. My version is derived directly
 *            from Brent's description and his ALGOL60 code. I've
 *            preserved his variable names as much as possible, to
 *            make the routine follow his published description
 *            closely. The Brent algorithm is also discussed in
 *            Numerical Recipes [Press88].
 */
static void
brent(double *ori, double *dir, long n,
      double (*func)(double *, long, void *), void *prm,
      double a, double b, double eps, double t,
      double *xvec, double *ret_x, double *ret_fx)
{
  double w,x,v,u;               /* with [a,b]: Brent's six points     */
  double m;                     /* midpoint of current [a,b] interval */
  double tol;                   /* tolerance = eps|x| + t */
  double fu,fv,fw,fx;           /* function evaluations */
  double p,q;                   /* numerator, denominator of parabolic interpolation */
  double r;
  double d,e;                   /* last, next-to-last values of p/q  */
  double c = 1. - (1./eslCONST_GOLD); /* Brent's c; 0.381966; golden ratio */
  int    niter;			/* number of iterations */

  x=v=w= a + c*(b-a);           /* initial guess of x by golden section */
  esl_vec_DCopy(ori, n, xvec);  /* build xvec from ori, dir, x */
  esl_vec_DAddScaled(xvec, dir, x, n);
  fx=fv=fw = (*func)(xvec, n, prm);   /* initial function evaluation */

  e     = 0.;
  niter = 0;
  while (1) /* algorithm is guaranteed to converge. */
    {
      m   = 0.5 * (a+b);
      tol = eps*fabs(x) + t;
      if (fabs(x-m) <= 2*tol - 0.5*(b-a)) break; /* convergence test. */
      niter++;

      p = q = r = 0.;
      if (fabs(e) > tol)
        { /* Compute parabolic interpolation, u = x + p/q */
          r = (x-w)*(fx-fv);
          q = (x-v)*(fx-fw);
          p = (x-v)*q - (x-w)*r;
          q = 2*(q-r);
          if (q > 0) { p = -p; } else {q = -q;}
          r = e;
          e=d;                  /* e is now the next-to-last p/q  */
        }

      if (fabs(p) < fabs(0.5*q*r) || p < q*(a-x) || p < q*(b-x))
        { /* Seems well-behaved? Use parabolic interpolation to compute new point u */
          d = p/q;              /* d remembers last p/q */
          u = x+d;              /* trial point, for now... */

          if (2.0*(u-a) < tol || 2.0*(b-u) < tol) /* don't evaluate func too close to a,b */
            d = (x < m)? tol : -tol;
        }
      else /* Badly behaved? Use golden section search to compute u. */
        {
          e = (x<m)? b-x : a-x;  /* e = largest interval */
          d = c*e;
        }

      /* Evaluate f(), but not too close to x.  */
      if      (fabs(d) >= tol) u = x+d;
      else if (d > 0)          u = x+tol;
      else                     u = x-tol;
      esl_vec_DCopy(ori, n, xvec);  /* build xvec from ori, dir, u */
      esl_vec_DAddScaled(xvec, dir, u, n);
      fu = (*func)(xvec, n, prm);   /* f(u) */

      /* Bookkeeping.  */
     if (fu <= fx)
        {
          if (u < x) b = x; else a = x;
          v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        }
      else
        {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x)
            { v = w; fv = fw; w = u; fw = fu; }
          else if (fu <= fv || v==x || v ==w)
            { v = u; fv = fu; }
        }
    }

  /* Return.
   */
  esl_vec_DCopy(ori, n, xvec);  /* build final xvec from ori, dir, x */
  esl_vec_DAddScaled(xvec, dir, x, n);
  if (ret_x  != NULL) *ret_x  = x;
  if (ret_fx != NULL) *ret_fx = fx;
  ESL_DPRINTF2(("\nbrent(): %d iterations\n", niter));
  ESL_DPRINTF2(("xx=%10.8f fx=%10.1f\n", x, fx));
}

/* Function:  min_Braket()
 * Incept:    ER, Wed Jun 22 08:49:42 2005 [janelia]
 * function from esl_minimize.c wrap as esl_min_GradientDescent()
 *
 * SRE, Wed Jul 27 11:43:32 2005 [St. Louis]
 *
 * Purpose:   Bracket a minimum. 
 *
 *            The minimization is quasi-one-dimensional, 
 *            starting from an initial <n>-dimension vector <ori>
 *            in the <n>-dimensional direction <d>.
 *            
 *            Caller passes a ptr to the objective function <*func()>,
 *            and a void pointer to any necessary conditional 
 *            parameters <prm>. The objective function will
 *            be evaluated at a point <x> by calling
 *            <(*func)(x, n, prm)>. The caller's function
 *            is responsible to casting <prm> to whatever it's
 *            supposed to be, which might be a ptr to a structure,
 *            for example; typically, for a parameter optimization
 *            problem, this holds the observed data.
 *            
 *            The routine works in scalar multipliers relative
 *            to origin <ori> and direction <d>; that is, a new <n>-dimensional
 *            point <b> is defined as <ori> + <bx><d>, for a scalar <bx>.
 *            
 *            The routine identifies a triplet <ax>, <bx>, <cx> such
 *            that $a < b < c$ and such that a minimum is known to
 *            exist in the $(a,b)$ interval because $f(b) < f(a),
 *            f(c)$. Also, the <a..b> and <b...c> intervals are in
 *            a golden ratio; the <b..c> interval is 1.618 times larger
 *            than <a..b>.
 *
 *            Since <d> is usually in the direction of the gradient,
 *            the points <ax>,<bx>,<cx> might be expected to be $\geq 0$;
 *            however, when <ori> is already close to the minimum, 
 *            it is often faster to bracket the minimum using
 *            a negative <ax>. The caller might then try to be "clever"
 *            and assume that the minimum is in the <bx..cx> interval
 *            when <ax> is negative, rather than the full <ax..cx>
 *            interval. That cleverness can fail, though, if <ori>
 *            is already in fact the minimum, because the line minimizer
 *            in brent() assumes a non-inclusive interval. Use
 *            <ax..cx> as the bracket.
 *            
 * Args:      ori       - n-dimensional starting vector
 *            d         - n-dimensional direction to minimize along
 *            n         - # of dimensions
 *            firststep - bx is initialized to this scalar multiplier
 *            *func()   - objective function to minimize
 *            prm       - void * to any constant data that *func() needs
 *            wrk       - workspace: 1 allocated n-dimensional vector
 *            ret_ax    - RETURN:  ax < bx < cx scalar bracketing triplet
 *            ret_bx    - RETURN:    ...ax may be negative
 *            ret_cx    - RETURN:    
 *            ret_fa    - RETURN:  function evaluated at a,b,c
 *            ret_fb    - RETURN:    ... f(b) < f(a),f(c)
 *            ret_fc    - RETURN:
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOHALT> if it fails to converge.
 *
 * Xref:      STL9/130.
 */
int
min_Bracket(double *x, double *dir, long n, double firststep,
	   double (*func)(double *, long, void *),
	   void *prm, double tol, double *wrk, double *ret_fx)
{
  double  oldfx;
  double *w1, *w2;
  double  fx;
  double  fa,fb,fc;
  double  ax,bx,cx;
  int     i;

  w1 = wrk + 2*n;
  w2 = wrk + 3*n;

  oldfx = (*func)(x, n, prm);	/* init the objective function */
  
  /* Bail out if the function is +/-inf or nan: this can happen if the caller
   * has screwed something up, or has chosen a bad start point.
   */
  if (! isfinite(oldfx)) ESL_EXCEPTION(eslERANGE, "minimum not finite");
  
  /* (failsafe) convergence test: a zero direction can happen, 
   * and it either means we're stuck or we're finished (most likely stuck)
   */
  for (i = 0; i < n; i++) 
    if (dir[i] != 0.) break;
  if  (i == n) {
    if (ret_fx != NULL) *ret_fx = oldfx;
    return eslOK;
  }
  
  /* Bracket the minimum.
   */
  bracket(x, dir, (int)n, firststep, func, prm, w1,
	  &ax, &bx, &cx, 
	  &fa, &fb, &fc);

#if 0
  /* guess of x by golden section */
  px = ax + cx*(bx-ax);      
  esl_vec_DAddScaled(x, dir, px, (int)n);
  fx = (*func)(x, n, prm);   /* function evaluation */
#else
  /* Minimize along the line given by the dir <dir> */
  brent(x, dir, (int)n, func, prm, ax, cx, tol, 1e-8, w2, NULL, &fx);
  esl_vec_DCopy(w2, n, x);
#endif

  /* Bail out if the function is now +/-inf: this can happen if the caller
   * has screwed something up.
   */
  if (fx == eslINFINITY || fx == -eslINFINITY)
    ESL_EXCEPTION(eslERANGE, "minimum not finite");
  
  if (ret_fx != NULL) *ret_fx = fx;
  
   return eslOK;
}



static int
NelderMead_initsimplex(double *x, double *u, long n, ESL_DMATRIX *simplex)
{
  int i;
  int j;

  for (i = 0; i <= n; i ++) 
    for (j = 0; j < n; j ++) 
      simplex->mx[i][j] = (j==i)? x[j] + u[j] : x[j];	
  
  return eslOK;
}

static int
NelderMead_indices(double *fx, ESL_DMATRIX *simplex, double (*func)(double *, long, void *), void *prm, 
		   int *ret_ih, int *ret_is, int *ret_il, double *ret_fh, double *ret_fs, double *ret_fl)
{
  double fl = +eslINFINITY;
  double fh = -eslINFINITY;
  double fs = -eslINFINITY;
  int    il;
  int    ih;
  int    is;
  int    i1;
  
  for (i1 = 0; i1 < simplex->n; i1 ++)  {
    fx[i1] = (*func)(simplex->mx[i1], simplex->m, prm); 
    if (fx[i1] < fl) { il = i1; fl = fx[i1]; }
    if (fx[i1] > fh) { ih = i1; fh = fx[i1]; }
  }    
  for (i1 = 0; i1 < simplex->n; i1 ++)  { // second highest
    if (fx[i1] > fs && i1 != ih) { is = i1; fs = fx[i1]; }
  }    
 
  *ret_ih = ih;
  *ret_is = is;
  *ret_il = il;
  *ret_fh = fh; // worst        (highest)
  *ret_fs = fs; // second worst (second highest)
  *ret_fl = fl; // best         (lowest)

  return eslOK;
}

static int
NelderMead_transform(double *w1, double *w2, double *fx, ESL_DMATRIX *simplex, double *c, double (*func)(double *, long, void *), void *prm, 
		     double alpha, double beta, double gamma, double delta, 
		     int *ret_ih, int *ret_is, int *ret_il, double *ret_fh, double *ret_fs, double *ret_fl, enum NMtransf_e  *ret_transf)
{
  enum NMtransf_e transf = TNONE;
  double          fr, fe, fc;
  double          fh = *ret_fh;
  double          fl = *ret_fl;
  double          fs = *ret_fs;
  int             ih = *ret_ih;
  int             il = *ret_il;
  int             is = *ret_is;
  int             n = simplex->m;
  int             i;
  int             i1;

  // try REFLECT
  for (i = 0; i < n; i++) 
    w1[i] = c[i] + alpha *(c[i] - simplex->mx[ih][i]);
  fr = (*func)(w1, n, prm); 

  if (fr < fs && fr >= fl) { // accept REFLECT
    transf = REFLECT;
    esl_vec_DCopy(w1, n, simplex->mx[ih]);
    fh = fr;
  }
  else if (fr < fl) { // try EXPAND
    for (i = 0; i < n; i++) 
      w2[i] = c[i] - gamma *(c[i] - w1[i]);
    fe = (*func)(w2, n, prm); 

    if (fe < fr) { // accept EXPAND
      transf = EXPAND;
      esl_vec_DCopy(w2, n, simplex->mx[ih]);
      fh = fe;
    }
    else {  // accept REFLECT
      transf = REFLECT;
      esl_vec_DCopy(w1, n, simplex->mx[ih]);
      fh = fr;
    }
  }
  else if (fr >= fs) { // try CONTRACT
    if (fr <  fh) { // outside CONTRACT
      for (i = 0; i < n; i++) 
	w2[i] = c[i] - beta *(c[i] - w1[i]);
      fc = (*func)(w2, n, prm); 
    }
    else { // inside CONTRACT
      for (i = 0; i < n; i++) 
	w2[i] = c[i] - beta *(c[i] - simplex->mx[ih][i]);
      fc = (*func)(w2, n, prm); 
    }

    if (fc < ESL_MIN(fr,fh)) { // accept CONTRACT
      transf = (fr < fh)? OCONTRACT : ICONTRACT;
      esl_vec_DCopy(w2, n, simplex->mx[ih]);
      fh = fc;
    }
    else { // SHRINK
      transf = SHRINK;
      for (i1 = 0; i1 < simplex->n; i1 ++) {
	if (i1 != il) {
	  for (i = 0; i < n; i++) 
	    simplex->mx[i1][i] = simplex->mx[il][i] + delta *(simplex->mx[i1][i] - simplex->mx[il][i]);
	  fx[i1] = (*func)(simplex->mx[i1], n, prm); 
	}
      }
      NelderMead_indices(fx, simplex, func, prm, &ih, &is, &il, &fh, &fs, &fl);
    }
  }

  if (transf == TNONE) { printf("transformation did not succeed\n"); return eslFAIL; }

  *ret_il     = il;
  *ret_is     = is;
  *ret_ih     = ih;
  *ret_fl     = fl;
  *ret_fs     = fs;
  *ret_fh     = fh;
  *ret_transf = transf;
  return eslOK;
}

static double
NelderMead_volume_update(double oldv, long n, double alpha, double beta, double gamma, double delta, enum NMtransf_e  transf)
{
  double v;

  switch(transf) {
  case REFLECT:
    v = oldv * alpha;
    break;
  case EXPAND:
    v = oldv * alpha * gamma;
    break;
  case OCONTRACT:
    v = oldv * alpha * beta;
   break;
  case ICONTRACT:
    v = oldv * beta;
    break;
  case SHRINK:
    v = exp(n*log(delta));
    break;
  default:
    printf("not a valid transformation %d\n", transf); exit(1);
  }
  return v;
}

static int
NelderMead_centroid(double *c, ESL_DMATRIX *simplex, int ih)
{
  int n1 = simplex->n;
  int n  = simplex->m;
  int i1, i;

  esl_vec_DSet(c, n, 0.0); // initialize to zero
  for (i = 0; i < n; i++) 
    for (i1 = 0; i1 < n1; i1++) 
      if (i1 != ih) c[i] += simplex->mx[i1][i];
  esl_vec_DScale(c, n, 1.0/(double)n);

  return eslOK;
}

static int
NelderMead_centroid_update(double *c, ESL_DMATRIX *simplex, int ih, int il, enum NMtransf_e transf)
{

  if      (transf == TNONE) 
    NelderMead_centroid(c, simplex, ih);
  else if (transf == SHRINK) { 
    NelderMead_centroid(c, simplex, ih);
  }
  else {
    NelderMead_centroid(c, simplex, ih);
  }
  
  return eslOK;
}

/* Function:  min_NelderMead()
 * Incept:    ER, Tue Apr  8 10:24:11 EDT 2014 [janelia]
 *
 *            Following Singer&Singer (2004)
 *
 * Purpose:   Multidimensional unconstrain optimization without
 *            derivatives by the Nelder-Mead algorithm.
 *
 * Args:      x        - an initial guess n-vector; RETURN: x at the minimum
 *            u        - stepsize to construct the initial simplex
 *            n        - dimensionality of all vectors
 *            *func()  - function for computing objective function f(x)
 *            prm      - void ptr to any data/params func,dfunc need 
 *            tol      - convergence criterion applied to f(x)
 *            wrk      - allocated 3x(n)+(n+1)-vector for workspace
 *            wrks     - allocated (n+1)x(n) dmatrix for workspace
 *            ret_fx   - optRETURN: f(x) at the minimum
 */
int
min_NelderMead(double *x, double *u, long n, 
	       double (*func)(double *, long, void *),
	       void *prm, double tol, double *wrk, ESL_DMATRIX *wrks, double *ret_fx)
{
  ESL_DMATRIX        *simplex;
  double             *w1, *w2, *fx, *cn;
  enum NMtransf_e     transf;
  enum NMtransf_e     oldtransf;
  double              alpha = 1.0;
  double              beta  = 0.5;
  double              gamma = 2.0;
  double              delta = 0.5;
  double              fh, fs, fl;
  double              oldfh, oldfs, oldfl;
  double              oldv, v;
  double              cvgf, cvgx;
  int                 n1 = n + 1;
  int                 ih, is, il;
  int                 oldih, oldis, oldil;
  int                 i;
  int                 it;

  oldfl = (*func)(x, n, prm);	/* init the objective function */
  
  /* Bail out if the function is +/-inf or nan: this can happen if the caller
   * has screwed something up, or has chosen a bad start point.
   */
  if (! isfinite(oldfl)) ESL_EXCEPTION(eslERANGE, "minimum not finite");
  
  /* assign the working space */
  fx      = wrk;
  w1      = wrk + n1;
  w2      = wrk + n1 + n;
  cn      = wrk + n1 + 2*n;
  simplex = wrks;
 
  /* (failsafe) convergence test: a zero direction can happen, 
   * and it either means we're stuck or we're finished (most likely stuck)
   */
  for (i = 0; i < n1; i++) 
    if (u[i] != 0.) break;
  if  (i == n1) {
    if (ret_fx != NULL) *ret_fx = oldfl;
    return eslOK;
  }
  
  NelderMead_initsimplex(x, u, n, simplex);
  oldih     = -1;
  oldis     = -1;
  oldil     = -1;
  oldv      = 1.0;
  oldtransf = TNONE;
  
  for (it = 0; it < MAXITERATIONS; it++)
    {
      /* determine indices h,s,l for worst, second-worst, and best points respectively */
      NelderMead_indices(fx, simplex, func, prm, &ih, &is, &il, &fh, &fs, &fl);
       
      /* calculate the centroid of the best side */
      NelderMead_centroid_update(cn, simplex, ih, il, oldtransf);
 
      /* compute the next working simplex */
      NelderMead_transform(w1, w2, fx, simplex, cn, func, prm, alpha, beta, gamma, delta, &ih, &is, &il, &fh, &fs, &fl, &transf);
      
      /* Main convergence test. 1e-10 factor is fudging the case where our
       * minimum is at exactly f()=0.
       */
      cvgf = 2.0 * fabs((fh-fl)) / (1e-10 + fabs(fh) + fabs(fl));
      if (cvgf <= tol) break;
     
      /* Second (failsafe) domain convergence test
       */
      v = NelderMead_volume_update(oldv, n, alpha, beta, gamma, delta, transf);
      cvgx = exp(1./(float)n * log(v));
      if (cvgx <= tol) break;
     
      oldih     = ih;  
      oldis     = is;  
      oldil     = il;  
      oldfh     = fh;  
      oldfs     = fs;  
      oldfl     = fl;  
      oldv      = v;
      oldtransf = transf;
    }
  
  /* Bail out if the function is now +/-inf: this can happen if the caller
   * has screwed something up.
   */
  if (fl == eslINFINITY || fl == -eslINFINITY)
    ESL_EXCEPTION(eslERANGE, "minimum not finite");
  
  esl_vec_DCopy(simplex->mx[il], n, x);
  if (ret_fx != NULL) *ret_fx = fl;
  
  if (it == MAXITERATIONS)
    ESL_FAIL(eslENOHALT, NULL, " ");
  
  return eslOK;
}


