# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "toms291.hpp"

//****************************************************************************80

double alogam ( double x, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    ALOGAM computes the logarithm of the Gamma function.
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by Malcolm Pike, David Hill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Malcolm Pike, David Hill,
//    Algorithm 291:
//    Logarithm of Gamma Function,
//    Communications of the ACM,
//    Volume 9, Number 9, September 1966, page 684.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.
//    X should be greater than 0.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double ALOGAM, the logarithm of the Gamma
//    function of X.
//
{
  double f;
  double value;
  double y;
  double z;

  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  y = x;

  if ( x < 7.0 )
  {
    f = 1.0;
    z = y;

    while ( z < 7.0 )
    {
      f = f * z;
      z = z + 1.0;
    }
    y = z;
    f = - log ( f );
  }
  else
  {
    f = 0.0;
  }

  z = 1.0 / y / y;

  value = f + ( y - 0.5 ) * log ( y ) - y 
    + 0.918938533204673 + 
    ((( 
    - 0.000595238095238   * z 
    + 0.000793650793651 ) * z 
    - 0.002777777777778 ) * z 
    + 0.083333333333333 ) / y;

  return value;
}
//****************************************************************************80

void gamma_log_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG_VALUES returns some values of the Log Gamma function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Log[Gamma[x]]
//
//  Modified:
//
//    14 August 2004
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
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 20

  double fx_vec[N_MAX] = { 
      0.1524063822430784E+01,  
      0.7966778177017837E+00,  
      0.3982338580692348E+00,  
      0.1520596783998375E+00,  
      0.0000000000000000E+00,  
     -0.4987244125983972E-01,  
     -0.8537409000331584E-01,  
     -0.1081748095078604E+00,  
     -0.1196129141723712E+00,  
     -0.1207822376352452E+00,  
     -0.1125917656967557E+00,  
     -0.9580769740706586E-01,  
     -0.7108387291437216E-01,  
     -0.3898427592308333E-01,  
     0.00000000000000000E+00,  
     0.69314718055994530E+00,  
     0.17917594692280550E+01,  
     0.12801827480081469E+02,  
     0.39339884187199494E+02,  
     0.71257038967168009E+02 }; 

  double x_vec[N_MAX] = { 
      0.20E+00,  
      0.40E+00,  
      0.60E+00,  
      0.80E+00,  
      1.00E+00,  
      1.10E+00,  
      1.20E+00,  
      1.30E+00,  
      1.40E+00,  
      1.50E+00,  
      1.60E+00,  
      1.70E+00,  
      1.80E+00,  
      1.90E+00,  
      2.00E+00,  
      3.00E+00,  
      4.00E+00,  
     10.00E+00,  
     20.00E+00,  
     30.00E+00 }; 

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
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

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
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
