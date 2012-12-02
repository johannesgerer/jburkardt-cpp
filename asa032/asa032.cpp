# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa032.hpp"

//****************************************************************************80

double alngam ( double xvalue, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    ALNGAM computes the logarithm of the gamma function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by Allan Macleod.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Allan Macleod,
//    Algorithm AS 245,
//    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
//    Applied Statistics,
//    Volume 38, Number 2, 1989, pages 397-402.
//
//  Parameters:
//
//    Input, double XVALUE, the argument of the Gamma function.
//
//    Output, int IFAULT, error flag.
//    0, no error occurred.
//    1, XVALUE is less than or equal to 0.
//    2, XVALUE is too big.
//
//    Output, double ALNGAM, the logarithm of the gamma function of X.
//
{
  double alr2pi = 0.918938533204673;
  double r1[9] = {
    -2.66685511495, 
    -24.4387534237, 
    -21.9698958928, 
     11.1667541262, 
     3.13060547623, 
     0.607771387771, 
     11.9400905721, 
     31.4690115749, 
     15.2346874070 };
  double r2[9] = {
    -78.3359299449, 
    -142.046296688, 
     137.519416416, 
     78.6994924154, 
     4.16438922228, 
     47.0668766060, 
     313.399215894, 
     263.505074721, 
     43.3400022514 };
  double r3[9] = {
    -2.12159572323E+05, 
     2.30661510616E+05, 
     2.74647644705E+04, 
    -4.02621119975E+04, 
    -2.29660729780E+03, 
    -1.16328495004E+05, 
    -1.46025937511E+05, 
    -2.42357409629E+04, 
    -5.70691009324E+02 };
  double r4[5] = {
     0.279195317918525, 
     0.4917317610505968, 
     0.0692910599291889, 
     3.350343815022304, 
     6.012459259764103 };
  double value;
  double x;
  double x1;
  double x2;
  double xlge = 510000.0;
  double xlgst = 1.0E+30;
  double y;

  x = xvalue;
  value = 0.0;
//
//  Check the input.
//
  if ( xlgst <= x )
  {
    *ifault = 2;
    return value;
  }

  if ( x <= 0.0 )
  {
    *ifault = 1;
    return value;
  }

  *ifault = 0;
//
//  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
//
  if ( x < 1.5 )
  {
    if ( x < 0.5 )
    {
      value = - log ( x );
      y = x + 1.0;
//
//  Test whether X < machine epsilon.
//
      if ( y == 1.0 )
      {
        return value;
      }
    }
    else
    {
      value = 0.0;
      y = x;
      x = ( x - 0.5 ) - 0.5;
    }

    value = value + x * (((( 
        r1[4]   * y 
      + r1[3] ) * y 
      + r1[2] ) * y 
      + r1[1] ) * y 
      + r1[0] ) / (((( 
                  y 
      + r1[8] ) * y 
      + r1[7] ) * y 
      + r1[6] ) * y 
      + r1[5] );

    return value;
  }
//
//  Calculation for 1.5 <= X < 4.0.
//
  if ( x < 4.0 )
  {
    y = ( x - 1.0 ) - 1.0;

    value = y * (((( 
        r2[4]   * x 
      + r2[3] ) * x 
      + r2[2] ) * x 
      + r2[1] ) * x 
      + r2[0] ) / (((( 
                  x 
      + r2[8] ) * x 
      + r2[7] ) * x 
      + r2[6] ) * x 
      + r2[5] );
  }
//
//  Calculation for 4.0 <= X < 12.0.
//
  else if ( x < 12.0 ) 
  {
    value = (((( 
        r3[4]   * x 
      + r3[3] ) * x 
      + r3[2] ) * x 
      + r3[1] ) * x 
      + r3[0] ) / (((( 
                  x 
      + r3[8] ) * x 
      + r3[7] ) * x 
      + r3[6] ) * x 
      + r3[5] );
  }
//
//  Calculation for 12.0 <= X.
//
  else
  {
    y = log ( x );
    value = x * ( y - 1.0 ) - 0.5 * y + alr2pi;

    if ( x <= xlge )
    {
      x1 = 1.0 / x;
      x2 = x1 * x1;

      value = value + x1 * ( ( 
             r4[2]   * 
        x2 + r4[1] ) * 
        x2 + r4[0] ) / ( ( 
        x2 + r4[4] ) * 
        x2 + r4[3] );
    }
  }

  return value;
}
//****************************************************************************80

double gamain ( double x, double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    GAMAIN computes the incomplete gamma ratio.
//
//  Discussion:
//
//    A series expansion is used if P > X or X <= 1.  Otherwise, a
//    continued fraction approximation is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by G Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    G Bhattacharjee,
//    Algorithm AS 32:
//    The Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 19, Number 3, 1970, pages 285-287.
//
//  Parameters:
//
//    Input, double X, P, the parameters of the incomplete 
//    gamma ratio.  0 <= X, and 0 < P.
//
//    Output, int *IFAULT, error flag.
//    0, no errors.
//    1, P <= 0.
//    2, X < 0.
//    3, underflow.
//    4, error return from the Log Gamma routine.
//
//    Output, double GAMAIN, the value of the incomplete gamma ratio.
//
{
  double a;
  double acu = 1.0E-08;
  double an;
  double arg;
  double b;
  double dif;
  double factor;
  double g;
  double gin;
  int i;
  double oflo = 1.0E+37;
  double pn[6];
  double rn;
  double term;
  double uflo = 1.0E-37;
  double value;
//
//  Check the input.
//
  if ( p <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  if ( x < 0.0 )
  {
    *ifault = 2;
    value = 0.0;
    return value;
  }

  if ( x == 0.0 )
  {
    *ifault = 0;
    value = 0.0;
    return value;
  }

  g = alngam ( p, ifault );

  if ( *ifault != 0 )
  {
    *ifault = 4;
    value = 0.0;
    return value;
  }

  arg = p * log ( x ) - x - g;

  if ( arg < log ( uflo ) )
  {
    *ifault = 3;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  factor = exp ( arg );
//
//  Calculation by series expansion.
//
  if ( x <= 1.0 || x < p )
  {
    gin = 1.0;
    term = 1.0;
    rn = p;

    for ( ; ; )
    {
      rn = rn + 1.0;
      term = term * x / rn;
      gin = gin + term;

      if ( term <= acu )
      {
        break;
      }
    }

    value = gin * factor / p;
    return value;
  }
//
//  Calculation by continued fraction.
//
  a = 1.0 - p;
  b = a + x + 1.0;
  term = 0.0;

  pn[0] = 1.0;
  pn[1] = x;
  pn[2] = x + 1.0;
  pn[3] = x * b;

  gin = pn[2] / pn[3];

  for ( ; ; )
  {
    a = a + 1.0;
    b = b + 2.0;
    term = term + 1.0;
    an = a * term;
    for ( i = 0; i <= 1; i++ )
    {
      pn[i+4] = b * pn[i+2] - an * pn[i];
    }

    if ( pn[5] != 0.0 )
    {
      rn = pn[4] / pn[5];
      dif = r8_abs ( gin - rn );
//
//  Absolute error tolerance satisfied?
//
      if ( dif <= acu )
      {
//
//  Relative error tolerance satisfied?
//
        if ( dif <= acu * rn )
        {
          value = 1.0 - factor * gin;
          break;
        }
      }
      gin = rn;
    }

    for ( i = 0; i < 4; i++ )
    {
      pn[i] = pn[i+2];
    }

    if ( oflo <= r8_abs ( pn[4] ) )
    {
      for ( i = 0; i < 4; i++ )
      {
        pn[i] = pn[i] / oflo;
      }
    }
  }

  return value;
}
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
//    17 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ABS, the absolute value of the argument.
//
{
  if ( 0.0 <= x )
  {
    return x;
  }
  else
  {
    return ( - x );
  }
}
//****************************************************************************80

void timestamp ( void )

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
