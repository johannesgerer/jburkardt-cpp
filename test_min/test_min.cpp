# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "test_min.hpp"

//****************************************************************************80

double p00_f ( int problem, double x )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F evaluates the function for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P00_F, the value of the objective function.
//
{
  double f;

  if ( problem == 1 )
  {
    f = p01_f ( x );
  }
  else if ( problem == 2 )
  {
    f = p02_f ( x );
  }
  else if ( problem == 3 )
  {
    f = p03_f ( x );
  }
  else if ( problem == 4 )
  {
    f = p04_f ( x );
  }
  else if ( problem == 5 )
  {
    f = p05_f ( x );
  }
  else if ( problem == 6 )
  {
    f = p06_f ( x );
  }
  else if ( problem == 7 )
  {
    f = p07_f ( x );
  }
  else if ( problem == 8 )
  {
    f = p08_f ( x );
  }
  else if ( problem == 9 )
  {
    f = p09_f ( x );
  }
  else if ( problem == 10 )
  {
    f = p10_f ( x );
  }
  else if ( problem == 11 )
  {
    f = p11_f ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_F - Fatal error!\n";
    cerr << "  Illegal problem number PROBLEM = " << problem << "\n";
    exit ( 1 );
  }

  return f;
}
//****************************************************************************80

double p00_f1 ( int problem, double x )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F1 evaluates the first derivative for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Input, double X, the value of the variable.
//
//    Output, double F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  if ( problem == 1 )
  {
    f1 = p01_f1 ( x );
  }
  else if ( problem == 2 )
  {
    f1 = p02_f1 ( x );
  }
  else if ( problem == 3 )
  {
    f1 = p03_f1 ( x );
  }
  else if ( problem == 4 )
  {
    f1 = p04_f1 ( x );
  }
  else if ( problem == 5 )
  {
    f1 = p05_f1 ( x );
  }
  else if ( problem == 6 )
  {
    f1 = p06_f1 ( x );
  }
  else if ( problem == 7 )
  {
    f1 = p07_f1 ( x );
  }
  else if ( problem == 8 )
  {
    f1 = p08_f1 ( x );
  }
  else if ( problem == 9 )
  {
    f1 = p09_f1 ( x );
  }
  else if ( problem == 10 )
  {
    f1 = p10_f1 ( x );
  }
  else if ( problem == 11 )
  {
    f1 = p11_f1 ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_F1 - Fatal error!\n";
    cerr << "  Illegal problem number PROBLEM = " << problem << "\n";
    exit ( 1 );
  }

  return f1;
}
//****************************************************************************80

double p00_f1_dif ( int problem, double x )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F1_DIF approximates the first derivative via finite differences.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Input, double X, the point where the gradient is to 
//    be approximated.
//
//    Output, double F1_DIF, the approximated gradient vector.
//
{
  double dx;
  double eps;
  double f1_dif;
  double fminus;
  double fplus;
  double xi;

  eps = pow ( r8_epsilon ( ), 0.33 );

  if ( 0.0 <= x )
  {
    dx = eps * ( x + 1.0 );
  }
  else
  {
    dx = eps * ( x - 1.0 );
  }

  xi = x;
  x = xi + dx;
  fplus = p00_f ( problem, x );

  x = xi - dx;
  fminus = p00_f ( problem, x );

  f1_dif = ( fplus - fminus ) / ( 2.0 * dx );

  x = xi;

  return f1_dif;
}
//****************************************************************************80

double p00_f2 ( int problem, double x )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F2 evaluates the second derivative for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double F2, the second derivative.
//
{
  double f2;

  if ( problem == 1 )
  {
    f2 = p01_f2 ( x );
  }
  else if ( problem == 2 )
  {
    f2 = p02_f2 ( x );
  }
  else if ( problem == 3 )
  {
    f2 = p03_f2 ( x );
  }
  else if ( problem == 4 )
  {
    f2 = p04_f2 ( x );
  }
  else if ( problem == 5 )
  {
    f2 = p05_f2 ( x );
  }
  else if ( problem == 6 )
  {
    f2 = p06_f2 ( x );
  }
  else if ( problem == 7 )
  {
    f2 = p07_f2 ( x );
  }
  else if ( problem == 8 )
  {
    f2 = p08_f2 ( x );
  }
  else if ( problem == 9 )
  {
    f2 = p09_f2 ( x );
  }
  else if ( problem == 10 )
  {
    f2 = p10_f2 ( x );
  }
  else if ( problem == 11 )
  {
    f2 = p11_f2 ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_F2 - Fatal error!\n";
    cerr << "  Illegal problem number PROBLEM = " << problem << "\n";
    exit ( 1 );
  }

  return f2;
}
//****************************************************************************80

double p00_f2_dif ( int problem, double x )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F2_DIF approximates the second derivative via finite differences.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Input, double X, the value of the variable.
//
//    Output, double F2_DIF, the approximate second derivative.
//
{
  double eps;
  double f00;
  double f2_dif;
  double fmm;
  double fpp;
  double s;
  double xi;
//
//  Choose the stepsize.
//
  eps = pow ( r8_epsilon ( ), 0.33 );

  s = eps * ( r8_abs ( x ) + 1.0 );

  xi = x;

  f00 = p00_f ( problem, x );

  x = xi + s;
  fpp = p00_f ( problem, x );

  x = xi - s;
  fmm = p00_f ( problem, x );

  f2_dif = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s / s;

  x = xi;

  return f2_dif;
}
//****************************************************************************80

double p00_fmin ( double *a, double *b, int problem, double tol )

//****************************************************************************80
//
//  Purpose:
//
//    P00_FMIN seeks a minimizer of a scalar function of a scalar variable.
//
//  Discussion:
//
//    FMIN seeks an approximation to the point where F attains a minimum on
//    the interval (A,B).
//
//    The method used is a combination of golden section search and
//    successive parabolic interpolation.  Convergence is never much
//    slower than that for a Fibonacci search.  If F has a continuous
//    second derivative which is positive at the minimum (which is not
//    at A or B), then convergence is superlinear, and usually of the
//    order of about 1.324....
//
//    The function F is never evaluated at two points closer together
//    than EPS * ABS ( FMIN ) + (TOL/3), where EPS is approximately the
//    square root of the relative machine precision.  If F is a unimodal
//    function and the computed values of F are always unimodal when
//    separated by at least EPS * ABS ( XSTAR ) + (TOL/3), then FMIN
//    approximates the abcissa of the global minimum of F on the
//    interval [A, B] with an error less than 3 * EPS * ABS ( FMIN ) + TOL.
//    If F is not unimodal, then FMIN may approximate a local, but
//    perhaps non-global, minimum to the same accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization without Derivatives,
//    Prentice Hall, 1973.
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1988.
//
//  Parameters
//
//    Input/output, double A, B.  On input, the left and right
//    endpoints of the initial interval.  On output, the lower and upper 
//    bounds for the minimizer.
//
//    Input, int PROBLEM, the index of a problem.
//
//    Input, double TOL, the desired length of the interval of
//    uncertainty of the final result.  TOL must not be negative.
//
//    Output, double P00_FMIN, the abcissa approximating the 
//    minimizer of f.
//
{
  double c;
  double d;
  double e;
  double eps;
  double fu;
  double fv;
  double fw;
  double fx;
  double midpoint;
  double p;
  double q;
  double r;
  double tol1;
  double tol2;
  double u;
  double v;
  double w;
  double x;

  c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );
//
//  C is the squared inverse of the golden ratio.
//
//  EPS is the square root of the relative machine precision.
//
  eps = sqrt ( r8_epsilon ( ) );
//
//  Initialization.
//
  v = *a + c * ( *b - *a );
  w = v;
  x = v;
  e = 0.0;
  fx = p00_f ( problem, x );
  fv = fx;
  fw = fx;
//
//  The main loop starts here.
//
  for ( ; ; )
  {
    midpoint = 0.5 * ( *a + *b );
    tol1 = eps * r8_abs ( x ) + tol / 3.0;
    tol2 = 2.0 * tol1;
//
//  Check the stopping criterion.
//
    if ( r8_abs ( x - midpoint ) <= ( tol2 - 0.5 * ( *b - *a ) ) )
    {
      break;
    }
//
//  Is golden-section necessary?
//
    if ( r8_abs ( e ) <= tol1 )
    {
      if ( midpoint <= x )
      {
        e = *a - x;
      }
      else  
      {
        e = *b - x;
      }

      d = c * e;
    }
//
//  Consider fitting a parabola.
//
    else
    {
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = -p;
      }
      q = r8_abs ( q );
      r = e;
      e = d;
//
//  Choose a golden-section step if the parabola is not advised.
//
      if ( 
        ( r8_abs ( 0.5 * q * r ) <= r8_abs ( p ) ) ||
        ( p <= q * ( *a - x ) ) ||
        ( q * ( *b - x ) <= p ) )
      {

        if ( midpoint <= x )
        {
          e = *a - x;
        }
        else
        {
          e = *b - x;
        }
        d = c * e;
      }
//
//  Choose a parabolic interpolation step.
//
      else
      {
        d = p / q;
        u = x + d;

        if ( ( u - *a ) < tol2 )
        {
          d = r8_abs ( tol1 ) * r8_sign ( midpoint - x );
        }

        if ( ( *b - u ) < tol2 )
        {
          d = r8_abs ( tol1 ) * r8_sign ( midpoint - x );
        }
     }
   }
//
//  F must not be evaluated too close to X.
//
    if ( tol1 <= r8_abs ( d ) )
    {
      u = x + d;
    }

    if ( r8_abs ( d ) < tol1 )
    {
      u = x + r8_abs ( tol1 ) * r8_sign ( d );
    }

    fu = p00_f ( problem, u );
//
//  Update the data.
//
    if ( fu <= fx )
    {
      if ( x <= u )
      {
        *a = x;
      }
      else
      {
        *b = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
      continue;
    }

    if ( u < x )
    {
      *a = u;
    }
    else
    {
      *b = u;
    }

    if ( fu <= fw || w == x )
    {
      v = w;
      fv = fw;
      w = u;
      fw = fu;
    }
    else if ( fu <= fv || v == x || v == w )
    {
      v = u;
      fv = fu;
    }
  }

  return x;
}
//****************************************************************************80

void p00_interval ( int problem, double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P00_INTERVAL returns a bracketing interval for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Output, double *A, *B, two points, between which a local
//    minimizer should be sought.
//
{
  if ( problem == 1 )
  {
    p01_interval ( a, b );
  }
  else if ( problem == 2 )
  {
    p02_interval ( a, b );
  }
  else if ( problem == 3 )
  {
    p03_interval ( a, b );
  }
  else if ( problem == 4 )
  {
    p04_interval ( a, b );
  }
  else if ( problem == 5 )
  {
    p05_interval ( a, b );
  }
  else if ( problem == 6 )
  {
    p06_interval ( a, b );
  }
  else if ( problem == 7 )
  {
    p07_interval ( a, b );
  }
  else if ( problem == 8 )
  {
    p08_interval ( a, b );
  }
  else if ( problem == 9 )
  {
    p09_interval ( a, b );
  }
  else if ( problem == 10 )
  {
    p10_interval ( a, b );
  }
  else if ( problem == 11 )
  {
    p11_interval ( a, b );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_INTERVAL - Fatal error!\n";
    cerr << "  Illegal problem number PROBLEM = " << problem << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

int p00_problem_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROBLEM_NUM returns the number of problems available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//   Output, int P00_PROBLEM_NUM, the number of problems.
//
{
  int problem_num;

  problem_num = 11;

  return problem_num;
}
//****************************************************************************80

void p00_sol ( int problem, int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P00_SOL returns the solution for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  if ( problem == 1 )
  {
    p01_sol ( know, x );
  }
  else if ( problem == 2 )
  {
    p02_sol ( know, x );
  }
  else if ( problem == 3 )
  {
    p03_sol ( know, x );
  }
  else if ( problem == 4 )
  {
    p04_sol ( know, x );
  }
  else if ( problem == 5 )
  {
    p05_sol ( know, x );
  }
  else if ( problem == 6 )
  {
    p06_sol ( know, x );
  }
  else if ( problem == 7 )
  {
    p07_sol ( know, x );
  }
  else if ( problem == 8 )
  {
    p08_sol ( know, x );
  }
  else if ( problem == 9 )
  {
    p09_sol ( know, x );
  }
  else if ( problem == 10 )
  {
    p10_sol ( know, x );
  }
  else if ( problem == 11 )
  {
    p11_sol ( know, x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_SOL - Fatal error!\n";
    cerr << "  Illegal problem number PROBLEM = " << problem << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p00_start ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_START returns a starting point for optimization for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Output, double P00_START, a starting point for the optimization.
//
{
  double x;

  if ( problem == 1 )
  {
    x = p01_start ( );
  }
  else if ( problem == 2 )
  {
    x = p02_start ( );
  }
  else if ( problem == 3 )
  {
    x = p03_start ( );
  }
  else if ( problem == 4 )
  {
    x = p04_start ( );
  }
  else if ( problem == 5 )
  {
    x = p05_start ( );
  }
  else if ( problem == 6 )
  {
    x = p06_start ( );
  }
  else if ( problem == 7 )
  {
    x = p07_start ( );
  }
  else if ( problem == 8 )
  {
    x = p08_start ( );
  }
  else if ( problem == 9 )
  {
    x = p09_start ( );
  }
  else if ( problem == 10 )
  {
    x = p10_start ( );
  }
  else if ( problem == 11 )
  {
    x = p11_start ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_START - Fatal error!\n";
    cerr << "  Illegal problem number PROBLEM = " << problem << "\n";
    exit ( 1 );
  }

  return x;
}
//****************************************************************************80

string p00_title ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns a title for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Output, string P00_TITLE, a title for the problem.
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
  else
  {
    cerr << "\n";
    cerr << " 'P00_TITLE - Fatal error!\n";
    cerr << "  Illegal problem number PROBLEM = " << problem << "\n";
    exit ( 1 );
  }

  return title;
}
//****************************************************************************80

double p01_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P01_F evaluates the objective function for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P01_F, the value of the objective function.
//
{
  double f;
 
  f = ( x - 2.0 ) * ( x - 2.0 ) + 1.0;

  return f;
}
//****************************************************************************80

double p01_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P01_F1 evaluates the first derivative for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P01_F1, the first derivative of the
//    objective function.
//
{
  double f1;

  f1 = 2.0 * ( x - 2.0 );

  return f1;
}
//****************************************************************************80

double p01_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P01_F2 evaluates the second derivative for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P01_F2, the second derivative.
//
{
  double f2;

  f2 = 2.0;

  return f2;
}
//****************************************************************************80

void p01_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P01_INTERVAL returns a starting interval for optimization for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a = 0.0;
  *b = 3.141592653589793;

  return;
}
//****************************************************************************80

void p01_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P01_SOL returns the solution for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = 2.0;

  return;
}
//****************************************************************************80

double p01_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_START returns a starting point for optimization for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P01_START, a starting point for the optimization.
//
{
  double x;

  x = 3.141592653589793;

  return x;
}
//****************************************************************************80

string p01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns a title for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P01_TITLE, a title for the problem.
//
{
  string title;

  title = "Simple quadratic, (x-2)^2+1.";

  return title;
}
//****************************************************************************80

double p02_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P02_F evaluates the objective function for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P02_F, the value of the objective function.
//
{
  double f;

  f = x * x + exp ( - x );

  return f;
}
//****************************************************************************80

double p02_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P02_F1 evaluates the first derivative for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P02_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = 2.0 * x - exp ( -x );

  return f1;
}
//****************************************************************************80

double p02_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P02_F2 evaluates the second derivative for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P02_2, the second derivative.
//
{
  double f2;

  f2 = 2.0 + exp ( -x );

  return f2;
}
//****************************************************************************80

void p02_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P02_INTERVAL returns a starting interval for optimization for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  0.0;
  *b =  1.0;

  return;
}
//****************************************************************************80

void p02_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P02_SOL returns the solution for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = 0.351734;

  return;
}
//****************************************************************************80

double p02_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_START returns a starting point for optimization for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P02_START, a starting point for the optimization.
//
{
  double x;

  x = 0.8;

  return x;
}
//****************************************************************************80

string p02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_TITLE returns a title for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P02_TITLE, a title for the problem.
//
{
  string title;

  title = "Quadratic plus exponential, x^2 + e^(-x).";

  return title;
}
//****************************************************************************80

double p03_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P03_F evaluates the objective function for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P03_F, the value of the objective function.
//
{
  double f;

  f = ( ( x * x + 2.0 ) * x + 1.0 ) * x + 3.0;

  return f;
}
//****************************************************************************80

double p03_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P03_F1 evaluates the first derivative for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P03_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = ( 4.0 * x * x + 4.0 ) * x + 1.0;

  return f1;
}
//****************************************************************************80

double p03_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P03_F2 evaluates the second derivative for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P03_F2, the second derivative.
//
{
  double f2;

  f2 = 12.0 * x * x + 4.0;

  return f2;
}
//****************************************************************************80

void p03_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P03_INTERVAL returns a starting interval for optimization for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  -2.0;
  *b =  +2.0;

  return;
}
//****************************************************************************80

void p03_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P03_SOL returns the solution for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = -0.236733;

  return;
}
//****************************************************************************80

double p03_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_START returns a starting point for optimization for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P03_START, a starting point for the optimization.
//
{
  double x;

  x = 1.5;

  return x;
}
//****************************************************************************80

string p03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_TITLE returns a title for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P03_TITLE, a title for the problem.
//
{
  string title;

  title = "Quartic, x^4 + 2x^2 + x + 3.";

  return title;
}
//****************************************************************************80

double p04_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P04_F evaluates the objective function for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P04_F, the value of the objective function.
//
{
  double f;

  f = exp ( x ) + 0.01 / x;

  return f;
}
//****************************************************************************80

double p04_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P04_F1 evaluates the first derivative for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P04_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = exp ( x ) - 0.01 / x / x;

  return f1;
}
//****************************************************************************80

double p04_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P04_F2 evaluates the second derivative for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P04_F2, the second derivative.
//
{
  double f2;

  f2 = exp ( x ) + 0.02 / x / x / x;

  return f2;
}
//****************************************************************************80

void p04_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P04_INTERVAL returns a starting interval for optimization for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  0.0001;
  *b =  1.0;

  return;
}
//****************************************************************************80

void p04_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P04_SOL returns the solution for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double X, the solution, if known.
//
{
  *know = 1;
  *x = 0.0953446;

  return;
}
//****************************************************************************80

double p04_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_START returns a starting point for optimization for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P04_START, a starting point for the optimization.
//
{
  double x;

  x = 0.95;

  return x;
}
//****************************************************************************80

string p04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_TITLE returns a title for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P04_TITLE, a title for the problem.
//
{
  string title;

  title = "Steep valley, e^x + 1/(100x).";

  return title;
}
//****************************************************************************80

double p05_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P05_F evaluates the objective function for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P05_F, the value of the objective function.
//
{
  double f;

  f = exp ( x ) - 2.0 * x + 0.01 / x - 0.000001 / x / x;

  return f;
}
//****************************************************************************80

double p05_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P05_F1 evaluates the first derivative for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P05_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = exp ( x ) - 2.0 - 0.01 / x / x + 0.000002 / x / x / x;

  return f1;
}
//****************************************************************************80

double p05_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P05_F2 evaluates the second derivative for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P05_F2, the second derivative.
//
{
  double f2;

  f2 = exp ( x ) + 0.02 / x / x / x - 0.000006 / x / x / x / x;

  return f2;
}
//****************************************************************************80

void p05_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P05_INTERVAL returns a starting interval for optimization for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  0.0002;
  *b =  2.0;

  return;
}
//****************************************************************************80

void p05_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P05_SOL returns the solution for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = 0.703206;

  return;
}
//****************************************************************************80

double p05_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_START returns a starting point for optimization for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P05_START, a starting point for the optimization.
//
{
  double x;

  x = 1.5;

  return x;
}
//****************************************************************************80

string p05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_TITLE returns a title for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P05_TITLE, a title for the problem.
//
{
  string title;

  title = "Steep valley, e^x - 2x + 1/(100x) - 1/(1000000x^2).";

  return title;
} 
//****************************************************************************80

double p06_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P06_F evaluates the objective function for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Prentice Hall 1973,
//    Reprinted Dover, 2002
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P06_F, the value of the objective function.
//
{
  double f;

  f = 2.0 - x;

  return f;
}
//****************************************************************************80

double p06_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P06_F1 evaluates the first derivative for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P06_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = -1.0;

  return f1;
}
//****************************************************************************80

double p06_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P06_F2 evaluates the second derivative for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    LE Scales,
//    Introduction to Non-Linear Optimization,
//    Springer, 1985.
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P06_F2, the second derivative.
//
{
  double f2;

  f2 = 0.0;

  return f2;
}
//****************************************************************************80

void p06_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P06_INTERVAL returns a starting interval for optimization for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  7.0;
  *b =  9.0;

  return;
}
//****************************************************************************80

void p06_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P06_SOL returns the solution for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = 9.0;

  return;
}
//****************************************************************************80

double p06_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_START returns a starting point for optimization for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P06_START, a starting point for the optimization.
//
{
  double x;

  x = 7.2;

  return x;
}
//****************************************************************************80

string p06_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_TITLE returns a title for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P06_TITLE, a title for the problem.
//
{
  string title;

  title = "line, 2 - x.";

  return title;
}
//****************************************************************************80

double p07_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P07_F evaluates the objective function for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Prentice Hall 1973,
//    Reprinted Dover, 2002
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P07_F, the value of the objective function.
//
{
  double f;

  f = ( x + sin ( x ) ) * exp ( - x * x );

  return f;
}
//****************************************************************************80

double p07_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P07_F1 evaluates the first derivative for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P07_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = ( 1.0 - 2.0 * x * x + cos ( x ) 
         - 2.0 * x * sin ( x ) ) * exp ( - x * x );

  return f1;
}
//****************************************************************************80

double p07_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P07_F2 evaluates the second derivative for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P07_F2, the second derivative.
//
{
  double f2;

  f2 = ( - 4.0 - 2.0 * x - 4.0 * x * x * x 
    - 3.0 * sin ( x ) - 4.0 * x * cos ( x ) 
    + 4.0 * x * x * sin ( x ) ) * exp ( - x * x );

  return f2;
}
//****************************************************************************80

void p07_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P07_INTERVAL returns a starting interval for optimization for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  -10.0;
  *b =  +10.0;

  return;
}
//****************************************************************************80

void p07_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P07_SOL returns the solution for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = -0.6795786599525;

  return;
}
//****************************************************************************80

double p07_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_START returns a starting point for optimization for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P07_START, a starting point for the optimization.
//
{
  double x;

  x = -5.0;

  return x;
}
//****************************************************************************80

string p07_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_TITLE returns a title for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P07_TITLE, a title for the problem.
//
{
  string title;

  title = "The dying snake, ( x + sin(x) ) * e^(-x^2).";

  return title;
}
//****************************************************************************80

double p08_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P08_F evaluates the objective function for problem 8.
//
//  Discussion:
//
//    This function looks positive, but has a pole at x = pi,
//    near which f -> negative infinity, and has two zeroes nearby.  
//    None of this will show up computationally.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arnold Krommer, Chriexit ( 1 );h Ueberhuber,
//    Numerical Integration on Advanced Systems,
//    Springer, 1994, pages 185-186.
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P08_F, the value of the objective function.
//
{
  double f;
  double pi = 3.141592653589793;

  if ( x == pi )
  {
    f = - 10000.0;
  }
  else
  {
    f = 3.0 * x * x + 1.0 + ( log ( ( x - pi ) * ( x - pi ) ) ) / pow ( pi, 4 );
  }

  return f;
}
//****************************************************************************80

double p08_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P08_F1 evaluates the first derivative for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P08_F1, the first derivative of the 
//    objective function.
//
{
  double f1;
  double pi = 3.141592653589793;

  if ( x == pi )
  {
    f1 = 0.0;
  }
  else
  {
    f1 = 6.0 * x + ( 2.0 / ( x - pi ) ) / pow ( pi, 4 );
  }

  return f1;
}
//****************************************************************************80

double p08_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P08_F2 evaluates the second derivative for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P08_F2, the second derivative.
//
{
  double f2;
  double pi = 3.141592653589793;

  if ( x == pi )
  {
    f2 = 1.0;
  }
  else
  {
    f2 = 6.0 + ( - 2.0 / ( x - pi ) / ( x - pi ) ) / pow ( pi, 4 );
  }

  return f2;
}
//****************************************************************************80

void p08_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P08_INTERVAL returns a starting interval for optimization for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  2.0;
  *b =  4.0;

  return;
}
//****************************************************************************80

void p08_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P08_SOL returns the solution for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  double pi = 3.141592653589793;

  *know = 1;
  *x = pi;

  return;
}
//****************************************************************************80

double p08_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_START returns a starting point for optimization for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P08_START, a starting point for the optimization.
//
{
  double x;

  x = 3.1;

  return x;
}
//****************************************************************************80

string p08_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_TITLE returns a title for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P08_TITLE, a title for the problem.
//
{
  string title;

  title = "The \"Thin Pole\", x^2+1+log((pi-x)^2)/pi^4";

  return title;
}
//****************************************************************************80

double p09_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P09_F evaluates the objective function for problem 9.
//
//  Discussion:
//
//    This function is oscillatory, with many local minima.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P09_F, the value of the objective function.
//
{
  double f;

  f = x * x - 10.0 * sin ( x * x - 3.0 * x + 2.0 );

  return f;
}
//****************************************************************************80

double p09_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P09_F1 evaluates the first derivative for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P09_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = 2.0 * x 
    - 10.0 * cos ( x * x - 3.0 * x + 2.0 ) 
    * ( 2.0 * x - 3.0 );

  return f1;
}
//****************************************************************************80

double p09_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P09_F2 evaluates the second derivative for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P09_F2, the second derivative.
//
{
  double f2;

  f2 = 2.0  
    + 10.0 * sin ( x * x - 3.0 * x + 2.0 ) 
    * ( 2.0 * x - 3.0 ) * ( 2.0 * x - 3.0 ) 
    - 20.0 * cos ( x * x - 3.0 * x + 2.0 );

  return f2;
}
//****************************************************************************80

void p09_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P09_INTERVAL returns a starting interval for optimization for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  -5.0;
  *b =  +5.0;

  return;
}
//****************************************************************************80

void p09_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P09_SOL returns the solution for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = 0.146621498932095;

  return;
}
//****************************************************************************80

double p09_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_START returns a starting point for optimization for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P09_START, a starting point for the optimization.
//
{
  double x;

  x = -2.0;

  return x;
}
//****************************************************************************80

string p09_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_TITLE returns a title for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P09_TITLE, a title for the problem.
//
{
  string title;

  title = "The oscillatory parabola";

  return title;
}
//****************************************************************************80

double p10_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P10_F evaluates the objective function for problem 10.
//
//  Discussion:
//
//    This function is oscillatory.
//
//    The function has a local minimum at 1.7922 whose function value is
//    very close to the minimum value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Isabel Beichl, Dianne O'Leary, Francis Sullivan,
//    Monte Carlo Minimization and Counting: One, Two, Too Many,
//    Computing in Science and Engineering,
//    Volume 9, Number 1, January/February 2007.
//
//    Dianne O'Leary,
//    Scientific Computing with Case Studies,
//    SIAM, 2008,
//    ISBN13: 978-0-898716-66-5,
//    LC: QA401.O44.
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P10_F, the value of the objective function.
//
{
  double f;

  f =       cos (       x ) 
    + 5.0 * cos ( 1.6 * x ) 
    - 2.0 * cos ( 2.0 * x ) 
    + 5.0 * cos ( 4.5 * x ) 
    + 7.0 * cos ( 9.0 * x );

  return f;
}
//****************************************************************************80

double p10_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P10_F1 evaluates the first derivative for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P10_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  f1 = -             sin (       x ) 
       - 5.0 * 1.6 * sin ( 1.6 * x ) 
       + 2.0 * 2.0 * sin ( 2.0 * x ) 
       - 5.0 * 4.5 * sin ( 4.5 * x ) 
       - 7.0 * 9.0 * sin ( 9.0 * x );

  return f1;
}
//****************************************************************************80

double p10_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P10_F2 evaluates the second derivative for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P10_F2, the second derivative.
//
{
  double f2;

  f2 = -                   cos (       x ) 
       - 5.0 * 1.6 * 1.6 * cos ( 1.6 * x ) 
       + 2.0 * 2.0 * 2.0 * cos ( 2.0 * x ) 
       - 5.0 * 4.5 * 4.5 * cos ( 4.5 * x ) 
       - 7.0 * 9.0 * 9.0 * cos ( 9.0 * x );

  return f2;
}
//****************************************************************************80

void p10_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P10_INTERVAL returns a starting interval for optimization for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  0.0;
  *b =  7.0;

  return;
}
//****************************************************************************80

void p10_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P10_SOL returns the solution for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = 5.975691087433868;

  return;
}
//****************************************************************************80

double p10_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_START returns a starting point for optimization for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P10_START, a starting point for the optimization.
//
{
  double x;

  x = 0.5;

  return x;
}
//****************************************************************************80

string p10_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_TITLE returns a title for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P10_TITLE, a title for the problem.
//
{
  string title;

  title = "The cosine combo";

  return title;
}
//****************************************************************************80

double p11_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P11_F evaluates the objective function for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the objective function.
//
//    Output, double P11_F, the value of the objective function.
//
{
  double f;

  f = 1.0 + r8_abs ( 3.0 * x - 1.0 );

  return f;
}
//****************************************************************************80

double p11_f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P11_F1 evaluates the first derivative for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the variable.
//
//    Output, double P11_F1, the first derivative of the 
//    objective function.
//
{
  double f1;

  if ( 3.0 * x - 1.0 < 0.0 )
  {
    f1 = - 3.0;
  }
  else
  {
    f1 = + 3.0;
  }

  return f1;
}
//****************************************************************************80

double p11_f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    P11_F2 evaluates the second derivative for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the values of the variables.
//
//    Output, double P11_F2, the second derivative.
//
{
  double f2;

  f2 = 0.0;

  return f2;
}
//****************************************************************************80

void p11_interval ( double *a, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    P11_INTERVAL returns a starting interval for optimization for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double *A, *B, two points defining an interval in which
//    the local minimizer should be sought.
//
{
  *a =  0.0;
  *b =  1.0;

  return;
}
//****************************************************************************80

void p11_sol ( int *know, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    P11_SOL returns the solution for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *KNOW.
//    If KNOW is 0, then the solution is not known.
//    If KNOW is positive, then the solution is known, and is returned in X.
//
//    Output, double *X, the solution, if known.
//
{
  *know = 1;
  *x = 1.0 / 3.0;

  return;
}
//****************************************************************************80

double p11_start ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_START returns a starting point for optimization for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P11_START, a starting point for the optimization.
//
{
  double x;

  x = 0.75;

  return x;
}
//****************************************************************************80

string p11_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_TITLE returns a title for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P11_TITLE, a title for the problem.
//
{
  string title;

  title = "1 + |3x-1|";

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
    value = - x;
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

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
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
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  } 
  else
  {
    value = 1.0;
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

