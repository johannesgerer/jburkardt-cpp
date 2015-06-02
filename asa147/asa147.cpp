# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa147.hpp"

//****************************************************************************80

void gamma_inc_values ( int *n_data, double *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function P(A,X) is defined as:
//
//      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= PN(A,X) <= 1
//
//    and
//
//      PN(A,INFINITY) = 1.0
//
//    In Mathematica, the function can be evaluated by:
//
//      1 - GammaRegularized[A,X]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 20

  double a_vec[N_MAX] = { 
     0.10E+00,  
     0.10E+00,  
     0.10E+00,  
     0.50E+00,  
     0.50E+00,  
     0.50E+00,  
     0.10E+01,  
     0.10E+01,  
     0.10E+01,  
     0.11E+01,  
     0.11E+01,  
     0.11E+01,  
     0.20E+01,  
     0.20E+01,  
     0.20E+01,  
     0.60E+01,  
     0.60E+01,  
     0.11E+02,  
     0.26E+02,  
     0.41E+02  };

  double fx_vec[N_MAX] = { 
     0.7382350532339351E+00,  
     0.9083579897300343E+00,  
     0.9886559833621947E+00,  
     0.3014646416966613E+00,  
     0.7793286380801532E+00,  
     0.9918490284064973E+00,  
     0.9516258196404043E-01,  
     0.6321205588285577E+00,  
     0.9932620530009145E+00,  
     0.7205974576054322E-01,  
     0.5891809618706485E+00,  
     0.9915368159845525E+00,  
     0.1018582711118352E-01,
     0.4421745996289254E+00,
     0.9927049442755639E+00,
     0.4202103819530612E-01,  
     0.9796589705830716E+00,  
     0.9226039842296429E+00,  
     0.4470785799755852E+00,  
     0.7444549220718699E+00 };

  double x_vec[N_MAX] = { 
     0.30E-01,  
     0.30E+00,  
     0.15E+01,  
     0.75E-01,  
     0.75E+00,  
     0.35E+01,  
     0.10E+00,  
     0.10E+01,  
     0.50E+01,  
     0.10E+00,   
     0.10E+01,  
     0.50E+01,  
     0.15E+00,  
     0.15E+01,  
     0.70E+01,  
     0.25E+01,  
     0.12E+02,  
     0.16E+02,  
     0.25E+02,  
     0.45E+02 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double gammds ( double x, double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMDS computes the incomplete Gamma integral.
//
//  Discussion:
//
//    The parameters must be positive.  An infinite series is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by Chi Leung Lau.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Chi Leung Lau,
//    Algorithm AS 147:
//    A Simple Series for the Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 29, Number 1, 1980, pages 113-114.
//
//  Parameters:
//
//    Input, double X, P, the arguments of the incomplete
//    Gamma integral.  X and P must be greater than 0.
//
//    Output, int *IFAULT, error flag.
//    0, no errors.
//    1, X <= 0 or P <= 0.
//    2, underflow during the computation.
//
//    Output, double GAMMDS, the value of the incomplete
//    Gamma integral.
//
{
  double a;
  double arg;
  double c;
  double e = 1.0E-09;
  double f;
  int ifault2;
  double uflo = 1.0E-37;
  double value;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  if ( p <= 0.0 ) 
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }
//
//  LGAMMA is the natural logarithm of the gamma function.
//
  arg = p * log ( x ) - lgamma ( p + 1.0 ) - x;

  if ( arg < log ( uflo ) )
  {
    value = 0.0;
    *ifault = 2;
    return value;
  }

  f = exp ( arg );

  if ( f == 0.0 )
  {
    value = 0.0;
    *ifault = 2;
    return value;
  }

  *ifault = 0;
//
//  Series begins.
//
  c = 1.0;
  value = 1.0;
  a = p;

  for ( ; ; )
  {
    a = a + 1.0;
    c = c * x / a;
    value = value + c;

    if ( c <= e * value )
    {
      break;
    }
  }

  value = value * f;

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
