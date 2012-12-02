# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa111.hpp"

//****************************************************************************80

void normal_01_cdf_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NormalDistribution [ 0, 1 ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2004
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
# define N_MAX 17

  double fx_vec[N_MAX] = { 
     0.5000000000000000E+00,  
     0.5398278372770290E+00,  
     0.5792597094391030E+00,  
     0.6179114221889526E+00,  
     0.6554217416103242E+00,  
     0.6914624612740131E+00,  
     0.7257468822499270E+00,  
     0.7580363477769270E+00,  
     0.7881446014166033E+00,  
     0.8159398746532405E+00,  
     0.8413447460685429E+00,  
     0.9331927987311419E+00,  
     0.9772498680518208E+00,  
     0.9937903346742239E+00,  
     0.9986501019683699E+00,  
     0.9997673709209645E+00,  
     0.9999683287581669E+00 };

  double x_vec[N_MAX] = { 
     0.0000000000000000E+00,    
     0.1000000000000000E+00,  
     0.2000000000000000E+00,  
     0.3000000000000000E+00,  
     0.4000000000000000E+00,  
     0.5000000000000000E+00,  
     0.6000000000000000E+00,  
     0.7000000000000000E+00,  
     0.8000000000000000E+00,  
     0.9000000000000000E+00,  
     0.1000000000000000E+01,  
     0.1500000000000000E+01,  
     0.2000000000000000E+01,  
     0.2500000000000000E+01,  
     0.3000000000000000E+01,  
     0.3500000000000000E+01,  
     0.4000000000000000E+01 };

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

double ppnd ( double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    PPND produces the normal deviate value corresponding to lower tail area = P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by J Beasley, S Springer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    J Beasley, S Springer,
//    Algorithm AS 111:
//    The Percentage Points of the Normal Distribution,
//    Applied Statistics,
//    Volume 26, Number 1, 1977, pages 118-121.
//
//  Parameters:
//
//    Input, double P, the value of the cumulative probability
//    densitity function.  0 < P < 1.
//
//    Output, integer *IFAULT, error flag.
//    0, no error.
//    1, P <= 0 or P >= 1.  PPND is returned as 0.
//
//    Output, double PPND, the normal deviate value with the property that
//    the probability of a standard normal deviate being less than or
//    equal to PPND is P.
//
{
  double a0 = 2.50662823884;
  double a1 = -18.61500062529;
  double a2 = 41.39119773534;
  double a3 = -25.44106049637;
  double b1 = -8.47351093090;
  double b2 = 23.08336743743;
  double b3 = -21.06224101826;
  double b4 = 3.13082909833;
  double c0 = -2.78718931138;
  double c1 = -2.29796479134;
  double c2 = 4.85014127135;
  double c3 = 2.32121276858;
  double d1 = 3.54388924762;
  double d2 = 1.63706781897;
  double r;
  double split = 0.42;
  double value;

  *ifault = 0;
//
//  0.08 < P < 0.92
//
  if ( r8_abs ( p - 0.5 ) <= split )
  {
    r = ( p - 0.5 ) * ( p - 0.5 );

    value = ( p - 0.5 ) * ( ( ( 
        a3   * r 
      + a2 ) * r 
      + a1 ) * r 
      + a0 ) / ( ( ( ( 
        b4   * r 
      + b3 ) * r 
      + b2 ) * r 
      + b1 ) * r 
      + 1.0 );
  }
//
//  P < 0.08 or P > 0.92,
//  R = min ( P, 1-P )
//
  else if ( 0.0 < p && p < 1.0 )
  {
    if ( 0.5 < p )
    {
      r = sqrt ( - log ( 1.0 - p ) );
    }
    else
    {
      r = sqrt ( - log ( p ) );
    }

    value = ( ( ( 
        c3   * r 
      + c2 ) * r 
      + c1 ) * r 
      + c0 ) / ( ( 
        d2   * r 
      + d1 ) * r 
      + 1.0 );

    if ( p < 0.5 )
    {
      value = - value;
    }
  }
//
//  P <= 0.0 or 1.0 <= P
//
  else
  {
    *ifault = 1;
    value = 0.0;
  }

  return value;
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

void timestamp ( )

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
