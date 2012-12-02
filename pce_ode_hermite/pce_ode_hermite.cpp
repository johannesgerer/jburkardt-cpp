# include <cstdlib>
# include <iostream>
# include <cmath>
# include <iomanip>
# include <ctime>

using namespace std;

# include "pce_ode_hermite.hpp"

//****************************************************************************80

double he_double_product_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
//
//  Discussion:
//
//    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
//    Princeton, 2010,
//    ISBN13: 978-0-691-14212-8,
//    LC: QA274.23.X58.
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//
//    Output, double HE_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  double value;

  if ( i == j )
  {
    value = r8_factorial ( i );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double he_triple_product_integral ( int i, int j, int k )

//****************************************************************************80
//
//  Purpose:
//
//    HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
//
//  Discussion:
//
//    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
//    Princeton, 2010,
//    ISBN13: 978-0-691-14212-8,
//    LC: QA274.23.X58.
//
//  Parameters:
//
//    Input, int I, J, K, the polynomial indices.
//
//    Output, double HE_TRIPLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  int s;
  double value;

  s = ( i + j + k ) / 2;

  if ( s < i || s < j || s < k )
  {
    value = 0.0;
  }
  else if ( ( ( i + j + k ) % 2 ) != 0 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial ( i ) / r8_factorial ( s - i ) 
          * r8_factorial ( j ) / r8_factorial ( s - j ) 
          * r8_factorial ( k ) / r8_factorial ( s - k );
  }

  return value;
}
//****************************************************************************80

void pce_ode_hermite ( double ti, double tf, int nt, double ui, int np, 
  double alpha_mu, double alpha_sigma, double t[], double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    PCE_ODE_HERMITE applies the polynomial chaos expansion to a scalar ODE.
//
//  Discussion:
//
//    The deterministic equation is
//
//      du/dt = - alpha * u,
//      u(0) = u0
//
//    In the stochastic version, it is assumed that the decay coefficient
//    ALPHA is a Gaussian random variable with mean value ALPHA_MU and variance
//    ALPHA_SIGMA^2.
//
//    The exact expected value of the stochastic equation will be
//
//      u(t) = u0 * exp ( t^2/2)
//
//    This should be matched by the first component of the polynomial chaos
//    expansion.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TI, TF, the initial and final times.
//
//    Input, int  NT, the number of output points.
//
//    Input, double UI, the initial condition.
//
//    Input, int NP, the degree of the expansion.  Polynomials 
//    of degree 0 through NP will be used.
//
//    Input, double ALPHA_MU, ALPHA_SIGMA, the mean and standard 
//    deviation of the decay coefficient.
//
//    Output, double T[NT+1], U[(NT+1)*(NP+1)], the times and the PCE 
//    coefficients at the successive time steps.
//
{
  double dp;
  double dt;
  int i;
  int it;
  int j;
  int k;
  double t1;
  double t2;
  double term;
  double tp;
  double *u1;
  double *u2;

  u1 = new double[np+1];
  u2 = new double[np+1];

  dt = ( tf - ti ) / ( double ) ( nt );
//
//  Set the PCE coefficients for the initial time.
//
  t1 = ti;

  u1[0] = ui;
  for ( j = 1; j <= np; j++ )
  {
    u1[j] = 0.0;
  }
//
//  Copy into the output arrays.
//
  t[0] = t1;
  for ( j = 0; j <= np; j++ )
  {
    u[0+j*(nt+1)] = u1[j];
  }
//
//  Time integration.
//
  for ( it = 1; it <= nt; it++ )
  {
    t2 = ( ( double ) ( nt - it ) * ti   
         + ( double ) (      it ) * tf ) 
         / ( double ) ( nt      );

    for ( k = 0; k <= np; k++ )
    {
      dp = he_double_product_integral ( k, k );

      term = - alpha_mu * u1[k];

      i = 1;
      for ( j = 0; j <= np; j++)
      {
        tp = he_triple_product_integral ( i, j, k );
        term = term - alpha_sigma * u1[j] * tp / dp;
      }
      u2[k] = u1[k] + dt * term;
    }
//
//  Prepare for next step.
//
    t1 = t2;
    for ( j = 0; j <= np; j++ )
    {
      u1[j] = u2[j];
    }
//
//  Copy into the output arrays.
//
    t[it] = t1;
    for ( j = 0; j <= np; j++ )
    {
      u[it+j*(nt+1)] = u1[j];
    }
  }

  delete [] u1;
  delete [] u2;

  return;
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
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
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
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
