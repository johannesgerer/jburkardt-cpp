# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hermite_test_int.hpp"

//****************************************************************************80

void hermite_compute ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
//
//  Discussion:
//
//    The abscissas are the zeros of the N-th order Hermite polynomial.
//
//    The integration interval is ( -oo, +oo ).
//
//    The weight function is w(x) = exp ( - x*x ).
//
//    The integral to approximate:
//
//      Integral ( -oo < X < +oo ) exp ( - X*X ) * F(X) dX
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
//    30 April 2006
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
//    Input, int ORDER, the order of the formula to be computed.
//
//    Output, double XTAB[ORDER], the Gauss-Hermite abscissas.
//
//    Output, double WEIGHT[ORDER], the Gauss-Hermite weights.
//
{
  double cc;
  double dp2;
  int i;
  double p1;
  double s;
  double temp;
  double x;

  cc = 1.7724538509 * r8_gamma ( ( double ) ( order ) ) 
    / pow ( 2.0, order - 1 );

  s = pow ( 2.0 * ( double ) ( order ) + 1.0, 1.0 / 6.0 );

  for ( i = 0; i < ( order + 1 ) / 2; i++ )
  {
    if ( i == 0 )
    {
      x = s * s * s - 1.85575 / s;
    }
    else if ( i == 1 )
    {
      x = x - 1.14 * pow ( ( double ) ( order ), 0.426 ) / x;
    }
    else if ( i == 2 )
    {
      x = 1.86 * x - 0.86 * xtab[0];
    }
    else if ( i == 3 )
    {
      x = 1.91 * x - 0.91 * xtab[1];
    }
    else
    {
      x = 2.0 * x - xtab[i-2];
    }

    hermite_root ( &x, order, &dp2, &p1 );

    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;

    xtab[order-i-1] = -x;
    weight[order-i-1] = weight[i];
  }
//
//  Reverse the order of the XTAB values.
//
  r8vec_reverse ( order, xtab );

  if ( false )
  {
  for ( i = 0; i < order/2; i++ )
  {
    temp            = xtab[i];
    xtab[i]         = xtab[order-1-i];
    xtab[order-1-i] = temp;
  }
  }
  return;
}
//****************************************************************************80

double hermite_integral ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_INTEGRAL returns the value of a Hermite polynomial integral.
//
//  Discussion:
//
//    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x*x) dx
//
//    H(n) is 0 for n odd.
//
//    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the integral.  
//    0 <= N.
//
//    Output, double VALUE, the value of the integral.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  if ( n < 0 )
  {
    value = - r8_huge ( );
  }
  else if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = ( double ) ( i4_factorial2 ( n - 1 ) ) * sqrt ( r8_pi ) 
      / pow ( 2.0, n / 2 );
  }

  return value;
}
//****************************************************************************80

void hermite_recur ( double *p2, double *dp2, double *p1, double x, int order )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
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
//    Output, double *P2, the value of H(ORDER)(X).
//
//    Output, double *DP2, the value of H'(ORDER)(X).
//
//    Output, double *P1, the value of H(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
{
  int i;
  double dq0;
  double dq1;
  double dq2;
  double q0;
  double q1;
  double q2;

  q1 = 1.0;
  dq1 = 0.0;

  q2 = x;
  dq2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    q0 = q1;
    dq0 = dq1;

    q1 = q2;
    dq1 = dq2;

    q2  = x * q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * q0;
    dq2 = x * dq1 + q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * dq0;
  }

  *p2 = q2;
  *dp2 = dq2;
  *p1 = q1;

  return;
}
//****************************************************************************80

void hermite_root ( double *x, int order, double *dp2, double *p1 )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_ROOT improves an approximate root of a Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2006
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
//    Input, int ORDER, the order of the Hermite polynomial.
//
//    Output, double *DP2, the value of H'(ORDER)(X).
//
//    Output, double *P1, the value of H(ORDER-1)(X).
//
{
  double d;
  double eps = 1.0E-12;
  double p2;
  int step;
  int step_max = 10;

  for ( step = 1; step <= step_max; step++ )
  {
    hermite_recur ( &p2, dp2, p1, *x, order );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }

  return;
}
//****************************************************************************80

int i4_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL2 computes the double factorial function.
//
//  Formula:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    N!!
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial function.
//    If N is less than 1, I4_FACTORIAL2 is returned as 1.
//
//    Output, int I4_FACTORIAL2, the value of the double factorial of N.
//
{
  int value;

  if ( n < 1 )
  {
    return 1;
  }

  value = 1;

  while ( 1 < n )
  {
    value = value * n;
    n = n - 2;
  }

  return value;
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
//    31 Julyl 2010
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

void p00_fun ( int problem, int option, int n, double x[], double f[] )

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
//    31 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  if ( problem == 1 )
  {
    p01_fun ( option, n, x, f );
  }
  else if ( problem == 2 )
  {
    p02_fun ( option, n, x, f );
  }
  else if ( problem == 3 )
  {
    p03_fun ( option, n, x, f );
  }
  else if ( problem == 4 )
  {
    p04_fun ( option, n, x, f );
  }
  else if ( problem == 5 )
  {
    p05_fun ( option, n, x, f );
  }
  else if ( problem == 6 )
  {
    p06_fun ( option, n, x, f );
  }
  else if ( problem == 7 )
  {
    p07_fun ( option, n, x, f );
  }
  else if ( problem == 8 )
  {
    p08_fun ( option, n, x, f );
  }
  else
  {
    cout << "\n";
    cout << "P00_FUN - Fatal error!\n";
    cout << "  Illegal problem number = " << problem << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p00_gauss_hermite ( int problem, int order )

//****************************************************************************80
//
//  Purpose:
//
//    P00_GAUSS_HERMITE applies a Gauss-Hermite quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Input, int ORDER, the order of the rule to apply.
//
//    Output, double P00_GAUSS_HERMITE, the approximate integral.
//
{
  double *f_vec;
  int i;
  int option;
  double result;
  double *weight;
  double *xtab;

  f_vec = new double[order];
  weight = new double[order];
  xtab = new double[order];

  hermite_compute ( order, xtab, weight );

  option = 1;
  p00_fun ( problem, option, order, xtab, f_vec );

  result = r8vec_dot_product ( order, weight, f_vec );

  delete [] f_vec;
  delete [] weight;
  delete [] xtab;

  return result;
}
//****************************************************************************80

double p00_monte_carlo ( int problem, int order )

//****************************************************************************80
//
//  Purpose:
//
//    P00_MONTE_CARLO applies a Monte Carlo procedure to Hermite integrals.
//
//  Discussion:
//
//    We wish to estimate the integral:
//
//      I(f) = integral ( -oo < x < +oo ) f(x) exp ( - x * x ) dx
//
//    We do this by a Monte Carlo sampling procedure, in which 
//    we select N points X(1:N) from a standard normal distribution,
//    and estimate
//
//      Q(f) = sum ( 1 <= I <= N ) f(x(i)) / sqrt ( pi )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2010
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
//    Output, double P00_MONTE_CARLO, the approximate integral.
//
{
  double *f_vec;
  int option;
  const double r8_pi = 3.141592653589793;
  double result;
  int seed;
  double weight;
  double *x_vec;

  seed = 123456789;
  x_vec = r8vec_normal_01_new ( order, &seed );

  option = 2;
  f_vec = new double[order];

  p00_fun ( problem, option, order, x_vec, f_vec );

  weight = ( double ) ( order ) / sqrt ( r8_pi ) / sqrt ( 2.0 );

  result = r8vec_sum ( order, f_vec ) / weight;

  delete [] f_vec;
  delete [] x_vec;

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
//    31 July 2010
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

  problem_num = 8;

  return problem_num;
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
//    02 February 2010
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

double p00_turing ( int problem, double h, double tol, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TURING applies the Turing quadrature rule.
//
//  Discussion
//
//    We consider the approximation:
//
//      Integral ( -oo < x < +oo ) f(x) dx
//
//      = h * Sum ( -oo < i < +oo ) f(nh) + error term
//
//    Given H and a tolerance TOL, we start summing at I = 0, and
//    adding one more term in the positive and negative I directions,
//    until the absolute value of the next two terms being added 
//    is less than TOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Turing,
//    A Method for the Calculation of the Zeta Function,
//    Proceedings of the London Mathematical Society,
//    Volume 48, 1943, pages 180-197.
//
//  Parameters:
//
//    Input, int PROBLEM, the index of the problem.
//
//    Input, double H, the spacing to use.
//
//    Input, double TOL, the tolerance.  
//
//    Output, int N, the number of pairs of steps taken.
//    The actual number of function evaluations is 2*N+1.
//
//    Output, double P00_TURING, the approximate integral.
//
{
  double f_vec[2];
  int i;
  int n_too_many = 100000;
  int option;
  int order;
  double result;
  double xtab[2];

  option = 0;
  *n = 0;

  result = 0.0;
  order = 1;
  xtab[0] = 0.0;
  p00_fun ( problem, option, order, xtab, f_vec );
  result = result + h * f_vec[0];

  for ( ; ; )
  {
    *n = *n + 1;

    xtab[0] =   ( double ) ( *n ) * h;
    xtab[1] = - ( double ) ( *n ) * h;

    order = 2;
    p00_fun ( problem, option, order, xtab, f_vec );

    result = result + h * ( f_vec[0] + f_vec[1] );
//
//  Just do a simple-minded absolute error tolerance check to start with.
//
    if ( r8_abs ( f_vec[0] ) + r8_abs ( f_vec[1] ) <= tol )
    {
      break;
    }
//
//  Just in case things go crazy.
//
    if ( n_too_many <= *n )
    {
      cout << "\n";
      cout << "P00_TURING - Warning!\n";
      cout << "  Number of steps exceeded N_TOO_MANY = " << n_too_many << "\n";
      break;
    }

  }

  return result;
}
//****************************************************************************80

double p01_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_EXACT returns the exact integral for problem 1.
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
//    Output, double P01_EXACT, the value of the integral.
//
{
  double exact;
  double omega = 1.0;
  const double r8_pi = 3.141592653589793;

  exact = sqrt ( r8_pi ) * exp ( - omega * omega );

  return exact;
}
//****************************************************************************80

void p01_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_FUN evaluates the integrand for problem 1.
//
//  Discussion:
//
//    Squire gives exact value as sqrt(pi) * exp(-w*w).
//
//    Integral ( -oo < x < +oo ) exp(-x*x) cos(2*w*x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Squire,
//    Comparison of Gauss-Hermite and Midpoint Quadrature with Application
//    to the Voigt Function,
//    in Numerical Integration: 
//    Recent Developments, Software and Applications,
//    edited by Patrick Keast, Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  double omega = 1.0;

  for ( i = 0; i < n; i++ )
  {
    f[i] = cos ( 2.0 * omega * x[i] );
  }

  if ( option == 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( - x[i] * x[i] );
    }
  }
  else if ( option == 1 )
  {
  }
  else if ( option == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( - x[i] * x[i] / 2.0 );
    }
  }

  return;
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
//    26 May 2009
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

  title = "exp(-x*x) * cos(2*omega*x)";

  return title;
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
//    31 July 2007
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
  const double r8_pi = 3.141592653589793;

  exact = sqrt ( r8_pi );

  return exact;
}
//****************************************************************************80

void p02_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_FUN evaluates the integrand for problem 2.
//
//  Discussion:
//
//    The exact value is sqrt(pi).
//
//    Integral ( -oo < x < +oo ) exp(-x*x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0;
  }

  if ( option == 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( - x[i] * x[i] );
    }
  }
  else if ( option == 1 )
  {
  }
  else if ( option == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( - x[i] * x[i] / 2.0 );
    }
  }

  return;
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
//    26 May 2009
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

  title = "exp(-x*x)";

  return title;
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
  double p = 1.0;
  double q = 3.0;
  const double r8_pi = 3.141592653589793;

  exact = r8_pi / ( q * sin ( r8_pi * p / q ) );

  return exact;
}
//****************************************************************************80

void p03_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_FUN evaluates the integrand for problem 3.
//
//  Discussion:
//
//    The exact value is pi / (q*sin(pi*p/q) ), assuming 0 < p < q.
//
//    Integral ( -oo < x < +oo ) exp(-px) / ( 1 + exp ( -qx) ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  double p = 1.0;
  double q = 3.0;

  for ( i = 0; i < n; i++ )
  {
    f[i] = exp ( - p * x[i] ) / ( 1.0 + exp ( -q * x[i] ) );
  }

  if ( option == 0 )
  {
  }
  else if ( option == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( + x[i] * x[i] );
    }
  }
  else if ( option == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( + x[i] * x[i] / 2.0 );
    }
  }

  return;
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
//    26 May 2009
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

  title = "exp(-px) / ( 1 + exp(-qx) )";

  return title;
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
  const double r8_pi = 3.141592653589793;

  exact = sqrt ( r8_pi / 2.0 );

  return exact;
}
//****************************************************************************80

void p04_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_FUN evaluates the integrand for problem 4.
//
//  Discussion:
//
//    The exact value is sqrt ( pi / 2 )
//
//    Integral ( -oo < x < +oo ) sin ( x^2 ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = sin ( x[i] * x[i] );
  }

  if ( option == 0 )
  {
  }
  else if ( option == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( + x[i] * x[i] );
    }
  }
  else if ( option == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( + x[i] * x[i] / 2.0 );
    }
  }
  return;
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
//    26 May 2009
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

  title = "sin(x^2)";

  return title;
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
  const double r8_pi = 3.141592653589793;

  exact = r8_pi / 3.0;

  return exact;
}
//****************************************************************************80

void p05_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_FUN evaluates the integrand for problem 5.
//
//  Discussion:
//
//    The exact value is pi / 3.
//
//    Integral ( -oo < x < +oo ) dx / ( (1+x^2) sqrt(4+3x^2) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0 / ( ( 1.0 + x[i] * x[i] ) * sqrt ( 4.0 + 3.0 * x[i] * x[i] ) );
  }

  if ( option == 0 )
  {
  }
  else if ( option == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( + x[i] * x[i] );
    }
  }
  else if ( option == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( + x[i] * x[i] / 2.0 );
    }
  }

  return;
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
//    26 May 2009
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

  title = "1/( (1+x^2) sqrt(4+3x^2) )";

  return title;
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
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P06_EXACT, the value of the integral.
//
{
  double exact;
  int m;
  const double r8_pi = 3.141592653589793;

  p06_param ( 'G', 'M', &m );

  if ( m <= -1 )
  {
    exact = - r8_huge ( );
  }
  else if ( ( m % 2 ) == 1 )
  {
    exact = 0.0;
  }
  else
  {
    exact = ( double ) ( i4_factorial2 ( m - 1 ) ) * sqrt ( r8_pi ) 
      / pow ( 2.0, m / 2 );
  }

  return exact;
}
//****************************************************************************80

void p06_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_FUN evaluates the integrand for problem 6.
//
//  Discussion:
//
//    The exact value is (m-1)!! * sqrt ( pi ) / sqrt ( 2**m ).
//
//    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
//
//    The parameter M is set by calling P06_PARAM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;
  int m;

  p06_param ( 'G', 'M', &m );

  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], m );
  }

  if ( option == 0 ) 
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( - x[i] * x[i] );
    }
  }
  else if ( option == 1 )
  {
  }
  else if ( option == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( - x[i] * x[i] / 2.0 );
    }
  }

  return;
}
//****************************************************************************80

void p06_param ( char action, char name, int *value )

//****************************************************************************80
//
//  Purpose:
//
//    P06_PARAM gets or sets parameters for problem 6.
//
//  Discussion:
//
//    The parameter is named "M", and it represents the value of the exponent
//    in the integrand function:
//
//    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
//
//    M must be greater than -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char ACTION, the action.
//    'S' to set the value,
//    'G' to get the value.
//
//    Input, char NAME, the parameter name.
//    'M', the exponent.
//
//    Input/output, int *VALUE, the parameter value.
//    If ACTION = 'S', then VALUE is an input quantity, and M is set to VALUE.
//    If ACTION = 'G', then VALUE is an output quantity, and VALUE is set to M.
//
{
  static int m = 0;

  if ( action == 'S' || action == 's' )
  {
    if ( *value <= -1 )
    {
      cerr << "\n";
      cerr << "P06_PARAM - Fatal error!\n";
      cerr << "  Parameter M must be greater than -1.\n";
      exit ( 1 );
    }
    m = *value;
  }
  else if ( action == 'G' || action == 'g' )
  {
    *value = m;
  }
  else
  {
    cerr << "\n";
    cerr << "P06_PARAM - Fatal error!\n";
    cerr << "  Unrecognized value of ACTION = \"" << action << "\".\n";
    exit ( 1 );
  }
  return;
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
//    26 May 2009
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

  title = "x^m exp(-x*x)";

  return title;
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
//    02 February 2010
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
  double e_sqrt_sqrt = 1.2840254166877414841;
  double exact;
  const double r8_pi = 3.141592653589793;

  exact = 0.25 * sqrt ( r8_pi ) / e_sqrt_sqrt;

  return exact;
}
//****************************************************************************80

void p07_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_FUN evaluates the integrand for problem 7.
//
//  Discussion:
//
//    The exact value is (1/4) sqrt(pi) / sqrt(sqrt(e)).
//
//    Integral ( -oo < x < +oo ) x^2 cos(x) e^(-x^2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  if ( option == 0 ) 
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = pow ( x[i], 2 ) * cos ( x[i] ) * exp ( - x[i] * x[i] );
    }
  }
  else if ( option == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = pow ( x[i], 2 ) * cos ( x[i] );
    }
  }
  else if ( option == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = pow ( x[i], 2 ) * cos ( x[i] ) * exp ( - x[i] * x[i] / 2.0 );
    }
  }

  return;
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
//    02 February 2010
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

  title = "x^2 cos ( x ) exp(-x*x)";

  return title;
}
//****************************************************************************80

double p08_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_EXACT returns the exact integral for problem 8.
//
//  Discussion:
//
//    The 20 digit value of the answer was computed by Mathematica.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P08_EXACT, the value of the integral.
//
{
  double exact;

  exact = 3.0088235661136433510;

  return exact;
}
//****************************************************************************80

void p08_fun ( int option, int n, double x[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_FUN evaluates the integrand for problem 8.
//
//  Discussion:
//
//    The exact value is sqrt ( 2 pi ) * HypergeometricU ( -1/2, 0, 1 ).
//
//    Integral ( -oo < x < +oo ) sqrt(1+x*x/2) * exp(-x*x/2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OPTION:
//    0, integrand is f(x).
//    1, integrand is exp(-x*x) * f(x);
//    2, integrand is exp(-x*x/2) * f(x);
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double F[N], the function values.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = sqrt ( 1.0 + 0.5 * x[i] * x[i] );
  }

  if ( option == 0 ) 
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( - 0.5 * x[i] * x[i] );
    }
  }
  else if ( option == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      f[i] = f[i] * exp ( + 0.5 * x[i] * x[i] );
    }
  }
  else if ( option == 2 )
  {

  }

  return;
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
//    31 July 2010
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

  title = "sqrt(1+x*x/2) * exp(-x*x/2)";

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
  const double r8_pi = 3.141592653589793;
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
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
  else if ( y < 12.0 )
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

  value = 0.1797693134862E+309;

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
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
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

double *r8vec_normal_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, real R(N+1), is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, real Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
  int i;
  int m;
  static int made = 0;
  double *r;
  const double r8_pi = 3.141592653589793;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;

  x = new double[n];
//
//  I'd like to allow the user to reset the internal data.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * r8_pi * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
}
//****************************************************************************80

void r8vec_reverse ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_REVERSE reverses the elements of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of double precision values.
//
//    Input:
//
//      N = 5,
//      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 )
//
//    Output:
//
//      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, double A[N], the array to be reversed.
//
{
  int i;
  double temp;

  for ( i = 0; i < n / 2; i++ )
  {
    temp     = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = temp;
  }
  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
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
