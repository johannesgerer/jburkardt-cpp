# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "laguerre_test_int.hpp"

//****************************************************************************80

void laguerre_compute ( int order, double xtab[], double weight[], 
  double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_COMPUTE computes a Gauss-Laguerre quadrature rule.
//
//  Discussion:
//
//    In the simplest case, ALPHA is 0, and we are approximating the
//    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
//    it is easy to modify the rule to approximate the integral from
//    A to +oo as well.
//
//    If ALPHA is nonzero, then there is no simple way to extend the
//    rule to approximate the integral from A to +oo.  The simplest
//    procedures would be to approximate the integral from 0 to A.
//
//    The integration interval is [ A, +oo ) or [ 0, +oo ).
//
//    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.
//
//
//    If the integral to approximate is:
//
//        Integral ( A <= X < +oo ) EXP ( - X ) * F(X) dX
//      or
//        Integral ( 0 <= X < +oo ) EXP ( - X ) * X^ALPHA * F(X) dX
//
//    then the quadrature rule is:
//
//      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
//    or
//      sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//    If the integral to approximate is:
//
//        Integral ( A <= X < +oo ) F(X) dX
//      or
//        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX
//
//    then the quadrature rule is:
//
//      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) 
//        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
//    or
//      sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the quadrature rule to be computed.
//    ORDER must be at least 1.
//
//    Output, double XTAB[ORDER], the Gauss-Laguerre abscissas.
//
//    Output, double WEIGHT[ORDER], the Gauss-Laguerre weights.
//
//    Input, double ALPHA, the exponent of the X factor.
//    Set ALPHA = 0.0 for the simplest rule.
//    ALPHA must be nonnegative.
//
{
  double *b;
  double *c;
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double ratio;
  double x;

  b = new double[order];
  c = new double[order];
//
//  Set the recursion coefficients.
//
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( alpha + ( double ) ( 2 * i + 1 ) );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i ) * ( alpha + ( double ) ( i ) );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = r8_gamma ( alpha + 1.0 ) * prod;

  for ( i = 0; i < order; i++ )
  {
//
//  Compute an estimate for the root.
//
    if ( i == 0 )
    {
      x = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) / 
        ( 1.0 + 2.4 * ( double ) ( order ) + 1.8 * alpha );
    }
    else if ( i == 1 )
    {
      x = x + ( 15.0 + 6.25 * alpha ) / 
        ( 1.0 + 0.9 * alpha + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) ) 
        / ( 1.9 * ( double ) ( i - 1 ) );

      r2 = 1.26 * ( double ) ( i - 1 ) * alpha / 
        ( 1.0 + 3.5 * ( double ) ( i - 1 ) );

      ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );

      x = x + ratio * ( x - xtab[i-2] );
    }
//
//  Use iteration to find the root.
//
    laguerre_root ( &x, order, alpha, &dp2, &p1, b, c );
//
//  Set the abscissa and weight.
//
    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;
  }

  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_RECUR finds the value and derivative of a Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of L(ORDER)(X).
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponent of the X factor in the
//    integrand.
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - alpha - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
//****************************************************************************80

void laguerre_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_ROOT improves an approximate root of a Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponent of the X factor.
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    laguerre_recur ( &p2, dp2, p1, *x, order, alpha, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void legendre_compute ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_COMPUTE computes a Gauss-Legendre quadrature rule.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x) = 1.0.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 April 2006
//
//  Author:
//
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be greater than 0.
//
//    Output, double XTAB[ORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[ORDER], the weights of the rule.
//    The weights are positive, symmetric, and should sum to 2.
//
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pk;
  double pkm1;
  double pkp1;
  const double r8_pi = 3.1415926535897932385;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    cout << "\n";
    cout << "LEGENDRE_COMPUTE - Fatal error!\n";
    cout << "  Illegal value of ORDER = " << order << "\n";
    exit ( 1 );
  }

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * r8_pi / ( double ) ( 4 * order + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) ) 
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
//
//  Initial approximation H:
//
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
//
//  Refine H using one step of Newton's method:
//
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    weight[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
  }

  if ( ( order % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }
//
//  Shift the data up.
//
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
    weight[iback-1] = weight[iback-ncopy-1];
  }
//
//  Reflect values for the negative abscissas.
//
  for ( i = 1; i <= order - nmove; i++ )
  {
    xtab[i-1] = - xtab[order-i];
    weight[i-1] = weight[order-i];
  }

  return;
}
//****************************************************************************80

double p00_alpha ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_ALPHA returns the value of ALPHA for any problem.
//
//  Discussion:
//
//    ALPHA is the lower, finite limit of integration in the integral.
//
//    The typical or default value is 0.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Output, double P00_ALPHA, the value of ALPHA.
//
{
  double alpha;

  if ( problem == 1 )
  {
    alpha = p01_alpha ( );
  }
  else if ( problem == 2 )
  {
    alpha = p02_alpha ( );
  }
  else if ( problem == 3 )
  {
    alpha = p03_alpha ( );
  }
  else if ( problem == 4 )
  {
    alpha = p04_alpha ( );
  }
  else if ( problem == 5 )
  {
    alpha = p05_alpha ( );
  }
  else if ( problem == 6 )
  {
    alpha = p06_alpha ( );
  }
  else if ( problem == 7 )
  {
    alpha = p07_alpha ( );
  }
  else if ( problem == 8 )
  {
    alpha = p08_alpha ( );
  }
  else if ( problem == 9 )
  {
    alpha = p09_alpha ( );
  }
  else if ( problem == 10 )
  {
    alpha = p10_alpha ( );
  }
  else if ( problem == 11 )
  {
    alpha = p11_alpha ( );
  }
  else if ( problem == 12 )
  {
    alpha = p12_alpha ( );
  }
  else if ( problem == 13 )
  {
    alpha = p13_alpha ( );
  }
  else if ( problem == 14 )
  {
    alpha = p14_alpha ( );
  }
  else if ( problem == 15 )
  {
    alpha = p15_alpha ( );
  }
  else if ( problem == 16 )
  {
    alpha = p16_alpha ( );
  }
  else if ( problem == 17 )
  {
    alpha = p17_alpha ( );
  }
  else if ( problem == 18 )
  {
    alpha = p18_alpha ( );
  }
  else if ( problem == 19 )
  {
    alpha = p19_alpha ( );
  }
  else if ( problem == 20 )
  {
    alpha = p20_alpha ( );
  }
  else
  {
    cout << "\n";
    cout << "P00_ALPHA - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return alpha;
}
//****************************************************************************80

double p00_exact ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_EXACT returns the exact integral for any problem.
//
//  Discussion:
//
//    This routine provides a "generic" interface to the exact integral
//    routines for the various problems, and allows a problem to be called
//    by index (PROBLEM) rather than by name.
//
//    In most cases, the "exact" value of the integral is not given;
//    instead a "respectable" approximation is available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Output, double P00_EXACT, the exact value of the integral.
//
{
  double exact;

  if ( problem == 1 )
  {
    exact = p01_exact ( );
  }
  else if ( problem == 2 )
  {
    exact = p02_exact ( );
  }
  else if ( problem == 3 )
  {
    exact = p03_exact ( );
  }
  else if ( problem == 4 )
  {
    exact = p04_exact ( );
  }
  else if ( problem == 5 )
  {
    exact = p05_exact ( );
  }
  else if ( problem == 6 )
  {
    exact = p06_exact ( );
  }
  else if ( problem == 7 )
  {
    exact = p07_exact ( );
  }
  else if ( problem == 8 )
  {
    exact = p08_exact ( );
  }
  else if ( problem == 9 )
  {
    exact = p09_exact ( );
  }
  else if ( problem == 10 )
  {
    exact = p10_exact ( );
  }
  else if ( problem == 11 )
  {
    exact = p11_exact ( );
  }
  else if ( problem == 12 )
  {
    exact = p12_exact ( );
  }
  else if ( problem == 13 )
  {
    exact = p13_exact ( );
  }
  else if ( problem == 14 )
  {
    exact = p14_exact ( );
  }
  else if ( problem == 15 )
  {
    exact = p15_exact ( );
  }
  else if ( problem == 16 )
  {
    exact = p16_exact ( );
  }
  else if ( problem == 17 )
  {
    exact = p17_exact ( );
  }
  else if ( problem == 18 )
  {
    exact = p18_exact ( );
  }
  else if ( problem == 19 )
  {
    exact = p19_exact ( );
  }
  else if ( problem == 20 )
  {
    exact = p20_exact ( );
  }
  else
  {
    cout << "\n";
    cout << "P00_EXACT - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return exact;
}
//****************************************************************************80

double p00_exp_transform ( int problem, int order )

//****************************************************************************80
//
//  Purpose:
//
//    P00_EXP_TRANSFORM applies an exponential transform and Gauss-Legendre rule.
//
//  Discussion:
//
//    To approximate:
//
//      Integral ( alpha <= x < +oo ) f(x) dx
//
//    Transform:
//
//      u = exp ( -x )
//      du = - exp ( -x ) dx
//
//      x = - log ( u )
//      dx = - du / u
//
//      x = alpha    => u = exp ( -alpha )
//      x = Infinity => u = 0
//
//    Transformed integral:
//
//      Integral ( 0 <= u <= exp ( -alpha ) ) f ( -log(u) ) du / u
//
//    We apply a Gauss-Legendre rule here, but we could easily use any rule
//    that avoids evaluation at U = 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Input, int ORDER, the order of the Gauss-Legendre rule 
//    to apply.
//
//    Output, double P00_EXP_TRANSFORM, the approximate integral.
//
{
  double alpha;
  double *fu;
  int i;
  double result;
  double *u;
  double *u_log;
  double *weight;

  u = new double[order];
  u_log = new double[order];
  weight = new double[order];

  alpha = p00_alpha ( problem );
//
//  Get the abscissas and weights for Gauss-Legendre quadrature.
//
  legendre_compute ( order, u, weight );
//
//  Modify the weights from [-1,1] to [0,exp(-alpha)].
//
  for ( i = 0; i < order; i++ )
  {
    weight[i] = exp ( -alpha ) * weight[i] / 2.0;
  }
//
//  Linear transform of abscissas from [-1,1] to [0,exp(-alpha)].
//
  for ( i = 0; i < order; i++ )
  {
    u[i] = ( ( 1.0 + u[i] ) * exp ( - alpha )
           + ( 1.0 - u[i] ) * 0.0 )
           / ( 2.0        );
  }
//
//  Define U_LOG = - log ( U )
//
  for ( i = 0; i < order; i++ )
  {
    u_log[i] = - log ( u[i] );
  }
//
//  Evaluate F ( -LOG(U) ).
//
  fu = p00_fun ( problem, order, u_log );
//
//  The integrand is F ( -LOG(U) ) / U
//
  for ( i = 0; i < order; i++ )
  {
    fu[i] = fu[i] / u[i];
  }
//
//  Sum.
//
  result = r8vec_dot_product ( order, weight, fu );

  delete [] fu;
  delete [] u;
  delete [] u_log;
  delete [] weight;

  return result;
}
//****************************************************************************80

double *p00_fun ( int problem, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_FUN evaluates the integrand for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P00_FUN[N], the function values.
//
{
  double *f;

  if ( problem == 1 )
  {
    f = p01_fun ( n, x );
  }
  else if ( problem == 2 )
  {
    f = p02_fun ( n, x );
  }
  else if ( problem == 3 )
  {
    f = p03_fun ( n, x );
  }
  else if ( problem == 4 )
  {
    f = p04_fun ( n, x );
  }
  else if ( problem == 5 )
  {
    f = p05_fun ( n, x );
  }
  else if ( problem == 6 )
  {
    f = p06_fun ( n, x );
  }
  else if ( problem == 7 )
  {
    f = p07_fun ( n, x );
  }
  else if ( problem == 8 )
  {
    f = p08_fun ( n, x );
  }
  else if ( problem == 9 )
  {
    f = p09_fun ( n, x );
  }
  else if ( problem == 10 )
  {
    f = p10_fun ( n, x );
  }
  else if ( problem == 11 )
  {
    f = p11_fun ( n, x );
  }
  else if ( problem == 12 )
  {
    f = p12_fun ( n, x );
  }
  else if ( problem == 13 )
  {
    f = p13_fun ( n, x );
  }
  else if ( problem == 14 )
  {
    f = p14_fun ( n, x );
  }
  else if ( problem == 15 )
  {
    f = p15_fun ( n, x );
  }
  else if ( problem == 16 )
  {
    f = p16_fun ( n, x );
  }
  else if ( problem == 17 )
  {
    f = p17_fun ( n, x );
  }
  else if ( problem == 18 )
  {
    f = p18_fun ( n, x );
  }
  else if ( problem == 19 )
  {
    f = p19_fun ( n, x );
  }
  else if ( problem == 20 )
  {
    f = p20_fun ( n, x );
  }
  else
  {
    cout << "\n";
    cout << "P00_FUN - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return f;
}
//****************************************************************************80

double p00_gauss_laguerre ( int problem, int order )

//****************************************************************************80
//
//  Purpose:
//
//    P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Input, int ORDER, the order of the Gauss-Laguerre rule 
//    to apply.
//
//    Output, double P00_GAUSS_LAGUERRE, the approximate integral.
//
{
  double alpha;
  double alpha2;
  double *fx;
  int i;
  double result;
  double *weight;
  double *xtab;

  weight = new double[order];
  xtab = new double[order];

  alpha = p00_alpha ( problem );

  alpha2 = 0.0;
  laguerre_compute ( order, xtab, weight, alpha2 );

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = xtab[i] + alpha;
  }

   fx = p00_fun ( problem, order, xtab );
//
//  The Gauss-Laguerre rule assumes a weight of EXP(-X).
//
//  We need to multiply each F(X) by EXP(X) to implicitly 
//  adjust for this weight.
//
  for ( i = 0; i < order; i++ )
  {
    fx[i] = fx[i] * exp ( xtab[i] );
  }

  result = exp ( -alpha ) * r8vec_dot_product ( order, weight, fx );

  delete [] fx;
  delete [] weight;
  delete [] xtab;

  return result;
}
//****************************************************************************80

int p00_problem_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROBLEM_NUM returns the number of test integration problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P00_PROBLEM_NUM, the number of test problems.
//
{
  int problem_num;

  problem_num = 20;

  return problem_num;
}
//****************************************************************************80

double p00_rat_transform ( int problem, int order )

//****************************************************************************80
//
//  Purpose:
//
//    P00_RAT_TRANSFORM applies a rational transform and Gauss-Legendre rule.
//
//  Discussion:
//
//    To approximate:
//
//      Integral ( alpha <= x < +oo ) f(x) dx
//
//    Transform:
//
//      u = 1 / ( 1 + x )
//      du = - dx / ( 1 + x )^2
//
//      x = ( 1 - u ) / u
//      dx = - du / u^2
//
//      x = alpha    => u = 1 / ( 1 + alpha )
//      x = Infinity => u = 0
//
//    Transformed integral:
//
//      Integral ( 0 < u <= 1 / ( 1 + alpha ) ) f ( ( 1 - u ) / u ) du / u^2
//
//    We apply a Gauss-Legendre rule here, but we could easily use any rule
//    that avoids evaluation at U = 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Input, int ORDER, the order of the Gauss-Legendre rule 
//    to apply.
//
//    Output, double P00_RAT_TRANSFORM, the approximate integral.
//
{
  double alpha;
  double *fu;
  int i;
  double result;
  double *u;
  double *u_rat;
  double *weight;

  u = new double[order];
  u_rat = new double[order];
  weight = new double[order];

  alpha = p00_alpha ( problem );
//
//  Get the abscissas and weights for Gauss-Legendre quadrature.
//
  legendre_compute ( order, u, weight );
//
//  Modify the weights from [-1,1] to [0,1/(1+alpha)].
//
  for ( i = 0; i < order; i++ )
  {
    weight[i] = weight[i] / 2.0 / ( 1.0 + alpha );
  }
//
//  Linear transform of abscissas from [-1,1] to [0,exp(-alpha)].
//
  for ( i = 0; i < order; i++ )
  {
    u[i] = ( ( 1.0 + u[i] ) / ( 1.0 + alpha )
           + ( 1.0 - u[i] ) * 0.0 )
           / ( 2.0        );
  }
//
//  Define U_RAT = ( 1 - U ) / U.
//
  for ( i = 0; i < order; i++ )
  {
    u_rat[i] = ( 1.0 - u[i] ) / u[i];
  }
//
//  Evaluate F ( ( 1 - U ) / U ).
//
  fu = p00_fun ( problem, order, u_rat );
//
//  The integrand is F ( ( 1 - U ) / U ) / U^2
//
  for ( i = 0; i < order; i++ )
  {
    fu[i] = fu[i] / u[i] / u[i];
  }
//
//  Sum.
//
  result = r8vec_dot_product ( order, weight, fu );

  delete [] fu;
  delete [] u;
  delete [] u_rat;
  delete [] weight;

  return result;
}
//****************************************************************************80

string p00_title ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns the title for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Output, string P00_TITLE, the title of the problem.
//
{
  string title;

  if ( problem == 1 )
  {
    title = p01_title ( );
  }
  else if ( problem == 2 )
  {
    title = p02_title ( );
  }
  else if ( problem == 3 )
  {
    title = p03_title ( );
  }
  else if ( problem == 4 )
  {
    title = p04_title ( );
  }
  else if ( problem == 5 )
  {
    title = p05_title ( );
  }
  else if ( problem == 6 )
  {
    title = p06_title ( );
  }
  else if ( problem == 7 )
  {
    title = p07_title ( );
  }
  else if ( problem == 8 )
  {
    title = p08_title ( );
  }
  else if ( problem == 9 )
  {
    title = p09_title ( );
  }
  else if ( problem == 10 )
  {
    title = p10_title ( );
  }
  else if ( problem == 11 )
  {
    title = p11_title ( );
  }
  else if ( problem == 12 )
  {
    title = p12_title ( );
  }
  else if ( problem == 13 )
  {
    title = p13_title ( );
  }
  else if ( problem == 14 )
  {
    title = p14_title ( );
  }
  else if ( problem == 15 )
  {
    title = p15_title ( );
  }
  else if ( problem == 16 )
  {
    title = p16_title ( );
  }
  else if ( problem == 17 )
  {
    title = p17_title ( );
  }
  else if ( problem == 18 )
  {
    title = p18_title ( );
  }
  else if ( problem == 19 )
  {
    title = p19_title ( );
  }
  else if ( problem == 20 )
  {
    title = p20_title ( );
  }
  else
  {
    cout << "\n";
    cout << "P00_TITLE - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return title;
}
//****************************************************************************80

double p01_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_ALPHA returns ALPHA for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P01_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 2.0;

  return alpha;
}
//****************************************************************************80

double p01_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_EXACT returns the exact integral for problem 1.
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P01_EXACT, the value of the integral.
//
{
  double exact;

  exact = 0.19524754198276439152;

  return exact;
}
//****************************************************************************80

double *p01_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_FUN evaluates the integrand for problem 1.
//
//  Discussion:
//
//    D&R gives "exact" value as 0.19524753.
//    Mathematica returns        0.19524754198276439152...
//    D&R gives Laguerre(16) as  0.16623627...
//
//  Integral:
//
//    exp ( -2 ) Integral ( 2 <= x < +oo ) / ( x * log(x)^2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P01_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = exp ( -2.0 ) / ( x[i] * pow ( log ( x[i] ), 2 ) );
  }

  return f;
}
//****************************************************************************80

string p01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns the title for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P01_TITLE, the title of the problem.
//
{
  string title;

  title = "1 / ( x * log ( x )^2 )";

  return title;
}
//****************************************************************************80

double p02_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_ALPHA returns ALPHA for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P02_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 2.0;

  return alpha;
}
//****************************************************************************80

double p02_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_EXACT returns the exact integral for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P02_EXACT, the value of the integral.
//
{
  double exact;

  exact = 0.32510848278991335198;

  return exact;
}
//****************************************************************************80

double *p02_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_FUN evaluates the integrand for problem 2.
//
//  Discussion:
//
//    D&R gives "exact" value as 0.32510855.
//    Mathematica returns        0.32510848278991335198...
//    D&R gives Laguerre(16) as  0.19142399...
//
//  Integral:
//
//    exp ( -2 ) Integral ( 2 <= x < +oo ) / ( x * log(x)^(3/2) ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P02_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = exp ( -2.0 ) / ( x[i] * pow ( log ( x[i] ), 1.5 ) );
  }

  return f;
}
//****************************************************************************80

string p02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_TITLE returns the title for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P02_TITLE, the title of the problem.
//
{
  string title;

  title = "1 / ( x * log ( x )^(3/2) )";

  return title;
}
//****************************************************************************80

double p03_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_ALPHA returns ALPHA for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P03_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 2.0;

  return alpha;
}
//****************************************************************************80

double p03_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_EXACT returns the exact integral for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P03_EXACT, the value of the integral.
//
{
  double exact;

  exact = 13.628;

  return exact;
}
//****************************************************************************80

double *p03_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_FUN evaluates the integrand for problem 3.
//
//  Discussion:
//
//    D&R gives "exact" value as 13.628...
//    Mathematica returns        13.440045415012575106...
//    D&R gives Laguerre(16) as   0.44996932...
//
//  Integral:
//
//    exp ( -2 ) Integral ( 2 <= x < +oo ) / ( x^1.01 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P03_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = exp ( -2.0 ) * 1.0 / pow ( x[i], 1.01 );
  }

  return f;
}
//****************************************************************************80

string p03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_TITLE returns the title for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P03_TITLE, the title of the problem.
//
{
  string title;

  title = "1 / ( x^1.01 )";

  return title;
}
//****************************************************************************80

double p04_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_ALPHA returns ALPHA for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P04_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 2.0;

  return alpha;
}
//****************************************************************************80

double p04_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_EXACT returns the estimated integral for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P04_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = -0.0046848541335080643181;

  return exact;
}
//****************************************************************************80

double *p04_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_FUN evaluates the integrand for problem 4.
//
//  Discussion:
//
//    D&R gives "exact" value as -0.0046984...
//    Mathematica returns        -0.0046848541335080643181...
//    D&R gives Laguerre(16) as  -0.039258696...
//
//  Integral:
//
//    exp ( -2 ) Integral ( 2 <= x < +oo ) sin ( x ) / x dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P04_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      f[i] = exp ( -2.0 );
    }
    else
    {
      f[i] = exp ( -2.0 ) * sin ( x[i] ) / x[i];
    }
  }
  return f;
}
//****************************************************************************80

string p04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_TITLE returns the title for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P04_TITLE, the title of the problem.
//
{
  string title;

  title = "Sine integral";

  return title;
}
//****************************************************************************80

double p05_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_ALPHA returns ALPHA for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P05_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 2.0;

  return alpha;
}
//****************************************************************************80

double p05_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_EXACT returns the estimated integral for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P05_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.0015897286158592328774;

  return exact;
}
//****************************************************************************80

double *p05_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_FUN evaluates the integrand for problem 5.
//
//  Discussion:
//
//    D&R gives "exact" value as  0.00158973...
//    Mathematica returns         0.0015897286158592328774...
//    D&R gives Laguerre(16) as  -0.067859545...
//
//  Integral:
//
//    exp ( -2 ) Integral ( 2 <= x < +oo ) cos ( pi * x^2 / 2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P05_FUN[N], the function values.
//
{
  double *f;
  int i;
  const double r8_pi = 3.1415926535897932385;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = exp ( -2.0 ) * cos ( 0.5 * r8_pi * x[i] * x[i] );
  }

  return f;
}
//****************************************************************************80

string p05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_TITLE returns the title for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P05_TITLE, the title of the problem.
//
{
  string title;

  title = "Fresnel integral";

  return title;
}
//****************************************************************************80

double p06_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_ALPHA returns ALPHA for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P06_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 2.0;

  return alpha;
}
//****************************************************************************80

double p06_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_EXACT returns the exact integral for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P06_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.00056103711148387120640;

  return exact;
}
//****************************************************************************80

double *p06_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_FUN evaluates the integrand for problem 6.
//
//  Discussion:
//
//    D&R gives "exact" value as 0.0005610371...
//    Mathematica returns        0.00056103711148387120640...
//    D&R gives Laguerre(16) as  0.00056100775...
//
//  Integral:
//
//    exp ( -2 ) Integral ( 2 <= x < +oo ) exp ( -x^2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P06_FUN[N], the function values.
//
{
  double exponent_min = -80.0;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( - x[i] * x[i] < exponent_min )
    {
      f[i] = 0.0;
    }
    else
    {
      f[i] = exp ( -2.0 ) * exp ( - x[i] * x[i] );
    }
  }

  return f;
}
//****************************************************************************80

string p06_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_TITLE returns the title for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P06_TITLE, the title of the problem.
//
{
  string title;

  title = "Complementary error function";

  return title;
}
//****************************************************************************80

double p07_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_ALPHA returns ALPHA for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P07_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 2.0;

  return alpha;
}
//****************************************************************************80

double p07_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_EXACT returns the exact integral for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P07_EXACT, the value of the integral.
//
{
  double exact;

  exact = 0.16266891;

  return exact;
}
//****************************************************************************80

double *p07_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_FUN evaluates the integrand for problem 7.
//
//  Discussion:
//
//    D&R gives "exact" value as 0.16266891...
//    Mathematica does not return a value.
//    D&R gives Laguerre(16) as  0.097083064...
//
//  Integral:
//
//    exp ( -2 ) Integral ( 2 <= x < +oo ) sin ( x - 1 ) 
//      / sqrt ( x * ( x - 2 ) ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P07_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 2.0 )
    {
      f[i] = 0.0;
    }
    else
    {
      f[i] = exp ( -2.0 ) 
        * sin ( x[i] - 1.0 ) / sqrt ( x[i] * ( x[i] - 2.0 ) );
    }
  }

  return f;
}
//****************************************************************************80

string p07_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_TITLE returns the title for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P07_TITLE, the title of the problem.
//
{
  string title;

  title = "Bessel function";

  return title;
}
//****************************************************************************80

double p08_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_ALPHA returns ALPHA for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P08_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p08_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_EXACT returns the estimated integral for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P08_EXACT, the estimated value of the integral.
//
{
  double exact;
  const double r8_pi = 3.1415926535897932385;

  exact = r8_pi * r8_pi / 6.0;

  return exact;
} 
//****************************************************************************80

double *p08_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_FUN evaluates the integrand for problem 8.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) x / ( exp ( x ) - 1 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P08_FUN[N], the function values.
//
{
  double exponent_max = 80.0;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      f[i] = 1.0 / exp ( x[i] );
    }
    else if ( x[i] < exponent_max )
    {
      f[i] = x[i] / ( exp ( x[i] ) - 1.0 );
    }
    else
    {
      f[i] = 0.0;
    }
  }

  return f;
}
//****************************************************************************80

string p08_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_TITLE returns the title for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P08_TITLE, the title of the problem.
//
{
  string title;

  title = "Debye function";

  return title;
}
//****************************************************************************80

double p09_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_ALPHA returns ALPHA for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P09_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p09_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_EXACT returns the estimated integral for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P09_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 24.0;

  return exact;
}
//****************************************************************************80

double *p09_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P09_FUN evaluates the integrand for problem 9.
//
//  Discussion:
//
//    The integral is the definition of the Gamma function for
//    Z = 5, with exact value (Z-1)! = 24.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) x^4 exp ( -x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P09_FUN[N], the function values.
//
{
  double exponent_min = -80.0;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( -x[i] < exponent_min )
    {
      f[i] = 0.0;
    }
    else
    {
      f[i] = pow ( x[i], 4 ) * exp ( -x[i] );
    }
  }

  return f;
}
//****************************************************************************80

string p09_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_TITLE returns the title for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P09_TITLE, the title of the problem.
//
{
  string title;

  title = "Gamma(Z=5) function";

  return title;
}
//****************************************************************************80

double p10_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_ALPHA returns ALPHA for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P10_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p10_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_EXACT returns the estimated integral for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXACT, the estimated value of the integral.
//
{
  double exact;
  const double r8_pi = 3.1415926535897932385;

  exact = r8_pi / 2.0;

  return exact;
}
//****************************************************************************80

double *p10_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P10_FUN evaluates the integrand for problem 10.
//
//  Discussion:
//
//    S&S gives exact value as pi/2 = 1.5707963267948966192...
//    S&S gives Laguerre(16) as       1.5537377347...
//    S&S gives EXP_TRANSFORM(16) as  1.4293043007...
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) 1/(1+x*x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P10_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0 / ( 1.0 + x[i] * x[i] );
  }

  return f;
}
//****************************************************************************80

string p10_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_TITLE returns the title for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P10_TITLE, the title of the problem.
//
{
  string title;

  title = "1 / ( 1 + x*x )";

  return title;
}
//****************************************************************************80

double p11_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_ALPHA returns ALPHA for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P11_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p11_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_EXACT returns the estimated integral for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXACT, the estimated value of the integral.
//
{
  double exact;
  const double r8_pi = 3.1415926535897932385;

  exact = r8_pi;

  return exact;
}
//****************************************************************************80

double *p11_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P11_FUN evaluates the integrand for problem 11.
//
//  Discussion:
//
//    S&S gives exact value as pi =  3.1415926535897932385...
//    S&S gives Laguerre(16) as      2.6652685196...
//    S&S gives EXP_TRANSFORM(16) as 2.3629036166... 
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) 1/((1+x)*sqrt(x)) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P11_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      f[i] = 0.0;
    }
    else
    {
      f[i] = 1.0 / ( ( 1.0 + x[i] ) * sqrt ( x[i] ) );
    }
  }

  return f;
}
//****************************************************************************80

string p11_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_TITLE returns the title for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P11_TITLE, the title of the problem.
//
{
  string title;

  title = "1 / ( (1+x) * sqrt(x) )";

  return title;
}
//****************************************************************************80

double p12_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P12_ALPHA returns ALPHA for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P12_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p12_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P12_EXACT returns the estimated integral for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.5;

  return exact;
}
//****************************************************************************80

double *p12_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P12_FUN evaluates the integrand for problem 12.
//
//  Discussion:
//
//    S&S gives exact value as pi =  0.5
//    S&S gives Laguerre(16) as      0.5000000000...
//    S&S gives EXP_TRANSFORM(16) as 0.5019065783... 
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) exp ( -x ) * cos ( x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P12_FUN[N], the function values.
//
{
  double exponent_min = -80.0;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( - x[i] < exponent_min )
    {
      f[i] = 0.0;
    }
    else
    {
      f[i] = exp ( -x[i] ) * cos ( x[i] );
    }
  }

  return f;
}
//****************************************************************************80

string p12_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P12_TITLE returns the title for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P12_TITLE, the title of the problem.
//
{
  string title;

  title = "exp ( - x ) * cos ( x )";

  return title;
}
//****************************************************************************80

double p13_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P13_ALPHA returns ALPHA for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P13_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p13_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P13_EXACT returns the estimated integral for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXACT, the estimated value of the integral.
//
{
  double exact;
  const double r8_pi = 3.1415926535897932385;

  exact = r8_pi / 2.0;

  return exact;
}
//****************************************************************************80

double *p13_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P13_FUN evaluates the integrand for problem 13.
//
//  Discussion:
//
//    S&S gives exact value as pi/2 = 1.5707963267948966192...
//    S&S gives Laguerre(16) as       1.4399523793...
//    S&S gives EXP_TRANSFORM(16) as  1.3045186595...
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) sin ( x ) / x dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P13_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      f[i] = 1.0;
    }
    else
    {
      f[i] = sin ( x[i] ) / x[i];
    }
  }

  return f;
}
//****************************************************************************80

string p13_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P13_TITLE returns the title for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P13_TITLE, the title of the problem.
//
{
  string title;

  title = "sin(x) / x";

  return title;
}
//****************************************************************************80

double p14_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P14_ALPHA returns ALPHA for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P14_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p14_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P14_EXACT returns the estimated integral for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 1.0634618101722400407;

  return exact;
}
//****************************************************************************80

double *p14_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P14_FUN evaluates the integrand for problem 14.
//
//  Discussion:
//
//    S&S gives "exact" value as     1.0634618101...
//    Mathematica returns            1.0634618101722400407...
//    S&S gives Laguerre(16) as      1.0634713425...
//    S&S gives EXP_TRANSFORM(16) as 1.0634618101...
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) sin ( exp ( - x ) + exp ( - 4 x ) ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P14_FUN[N], the function values.
//
{
  double exponent_min = -80.0;
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( - x[i] < exponent_min )
    {
      f[i] = 0.0;
    }
    else if ( -4.0 * x[i] < exponent_min )
    {
      f[i] = sin ( exp ( -x[i] ) );
    }
    else
    {
      f[i] = sin ( exp ( -x[i] ) + exp ( -4.0 * x[i] ) );
    }
  }

  return f;
}
//****************************************************************************80

string p14_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P14_TITLE returns the title for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P14_TITLE, the title of the problem.
//
{
  string title;

  title = "sin ( exp(-x) + exp(-4x) )";

  return title;
}
//****************************************************************************80

double p15_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P15_ALPHA returns ALPHA for problem 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P15_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p15_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P15_EXACT returns the estimated integral for problem 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXACT, the estimated value of the integral.
//
{
  double exact;
  const double r8_pi = 3.1415926535897932385;

  exact = - r8_pi * log ( 10.0 ) / 20.0;

  return exact;
}
//****************************************************************************80

double *p15_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P15_FUN evaluates the integrand for problem 15.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) log(x) / (1+100*x*x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise deDoncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983,
//    ISBN: 3540125531,
//    LC: QA299.3.Q36.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P15_FUN[N], the function values.
//
{
  double *f;
  int i;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      f[i] = - r8_huge ( );
    }
    else
    {
      f[i] = log ( x[i] ) / ( 1.0 + 100.0 * x[i] * x[i] );
    }
  }
  return f;
}
//****************************************************************************80

string p15_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P15_TITLE returns the title for problem 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P15_TITLE, the title of the problem.
//
{
  string title;

  title = "log(x) / ( 1 + 100 x^2 )";

  return title;
}
//****************************************************************************80

double p16_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P16_ALPHA returns ALPHA for problem 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P16_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p16_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P16_EXACT returns the estimated integral for problem 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 1.0;

  return exact;
}
//****************************************************************************80

double *p16_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P16_FUN evaluates the integrand for problem 16.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) cos ( pi * x / 2 ) / sqrt ( x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise deDoncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983,
//    ISBN: 3540125531,
//    LC: QA299.3.Q36.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P16_FUN[N], the function values.
//
{
  double *f;
  int i;
  const double r8_pi = 3.1415926535897932385;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      f[i] = r8_huge ( );
    }
    else
    {
      f[i] = cos ( r8_pi * x[i] / 2.0 ) / sqrt ( x[i] );
    }
  }
  return f;
}
//****************************************************************************80

string p16_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P16_TITLE returns the title for problem 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P16_TITLE, the title of the problem.
//
{
  string title;

  title = "cos ( pi x / 2 ) / sqrt ( x )";

  return title;
}
//****************************************************************************80

double p17_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P17_ALPHA returns ALPHA for problem 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P17_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p17_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P17_EXACT returns the exact integral for problem 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P17_EXACT, the value of the integral.
//
{
  const double beta = 2.0;
  double exact;
  const double r8_pi = 3.1415926535897932385;

  exact = sqrt ( r8_pi ) * cos ( 0.5 * atan ( pow ( 2.0, beta ) ) ) 
    / sqrt ( sqrt ( 1.0 + pow ( 0.25, beta ) ) );

  return exact;
}
//****************************************************************************80

double *p17_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P17_FUN evaluates the integrand for problem 17.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo) exp ( - x / 2^beta ) * cos ( x ) / sqrt ( x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 84.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P17_FUN[N], the integrand values.
//
{
  static double beta = 2.0;
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = exp ( - x[i] / pow ( 2.0, beta ) ) * cos ( x[i] ) 
        / sqrt ( x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

string p17_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P17_TITLE returns the title for problem 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P17_TITLE, the title of the problem.
//
{
  string title;

  title = "exp ( - x / 2^beta ) * cos ( x ) / sqrt ( x )";

  return title;
}
//****************************************************************************80

double p18_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P18_ALPHA returns ALPHA for problem 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P18_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p18_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P18_EXACT returns the exact integral for problem 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P18_EXACT, the value of the integral.
//
{
  static double beta = 1.0;
  double exact;
 
  exact = pow ( 2.0, 3.0 * beta + 1.0 );

  return exact;
}
//****************************************************************************80

double *p18_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P18_FUN evaluates the integrand for problem 18.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) x^2 * exp ( - x / 2^beta ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 84.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P18_FUN[N], the integrand values.
//
{
  static double beta = 1.0;
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i] * x[i] * exp ( - x[i] / pow ( 2, beta ) );
  }
  return fx;
}
//****************************************************************************80

string p18_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P18_TITLE returns the title for problem 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P18_TITLE, the title of the problem.
//
{
  string title;

  title = "x^2 * exp ( - x / 2^beta )";

  return title;
}
//****************************************************************************80

double p19_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P19_ALPHA returns ALPHA for problem 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P19_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p19_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P19_EXACT returns the exact integral for problem 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P19_EXACT, the value of the integral.
//
{
  const double beta = 0.5;
  double exact;
  const double r8_pi = 3.1415926535897932385;

  if ( beta == 1.0 )
  {
    exact = 1.0 / 10.0;
  }
  else
  {
    exact = ( 1.0 - beta ) * r8_pi 
      / ( pow ( 10.0, beta ) * sin ( r8_pi * beta ) );
  }
  return exact;
}
//****************************************************************************80

double *p19_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P19_FUN evaluates the integrand for problem 19.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) x^(alpha-1) / ( 1 + 10 x )^2 dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 84.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P61_FUN[N], the integrand values.
//
{
  static double beta = 0.5;
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( beta == 1.0 )
    {
      fx[i] = 1.0 / pow ( 1.0 + 10.0 * x[i], 2 );
    }
    else if ( beta < 1.0 && x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = pow ( x[i], beta - 1.0 ) / pow ( 1.0 + 10.0 * x[i], 2 );
    }
  }
  return fx;
}
//****************************************************************************80

string p19_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P19_TITLE returns the title for problem 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P19_TITLE, the title of the problem.
//
{
  string title;

  title = "x^(beta-1) / ( 1 + 10 x )^2";

  return title;
}
//****************************************************************************80

double p20_alpha ( )

//****************************************************************************80
//
//  Purpose:
//
//    P20_ALPHA returns ALPHA for problem 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P20_ALPHA, the value of ALPHA.
//
{
  double alpha;

  alpha = 0.0;

  return alpha;
}
//****************************************************************************80

double p20_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P20_EXACT returns the exact integral for problem 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P20_EXACT, the value of the integral.
//
{
  static double beta = 1.0;
  double exact;

  exact = 
    ( 
      log ( 1.5 ) / pow ( 2.0, beta ) 
      - 1.0 / pow ( 2.0, beta + 1.0 ) * 
      log ( ( 16.0 + pow ( 0.25, beta ) ) / ( 1.0 + pow ( 0.25, beta ) ) ) 
      - atan ( pow ( 2.0, beta + 2.0 ) ) - atan ( pow ( 2.0, beta ) ) 
    ) / ( 1.0 + pow ( 0.25, beta ) );

  return exact;
}
//****************************************************************************80

double *p20_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P20_FUN evaluates the integrand for problem 20.
//
//  Integral:
//
//    Integral ( 0 <= x < +oo ) 
//      1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 84.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P20_FUN[N], the integrand values.
//
{
  static double beta = 1.0;
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( pow ( x[i] - 1.0, 2 ) + pow ( 0.25, beta ) == 0.0 || x[i] == 2.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = 1.0 / 
        ( pow ( 2.0, beta ) 
        * ( pow ( x[i] - 1.0, 2 ) + pow ( 0.25, beta ) ) 
        * ( x[i] - 2.0 ) );
    }
  }
  return fx;
}
//****************************************************************************80

string p20_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P20_TITLE returns the title for problem 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P20_TITLE, the title of the problem.
//
{
  string title;

  title = "1 / ( 2^beta * ( ( x - 1 )^2 + (1/4)^beta ) * ( x - 2 ) )";

  return title;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  const double r8_pi = 3.1415926535897932385;
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - r8_pi / sin ( r8_pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
