# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "toms443.hpp"

//****************************************************************************80

void lambert_w_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAMBERT_W_VALUES returns some values of the Lambert W function.
//
//  Discussion:
//
//    The function W(X) is defined implicitly by:
//
//      W(X) * e^W(X) = X
//
//    The function is also known as the "Omega" function.
//
//    In Mathematica, the function can be evaluated by:
//
//      W = ProductLog [ X ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Brian Hayes,
//    "Why W?",
//    The American Scientist,
//    Volume 93, March-April 2005, pages 104-108.
//
//    Eric Weisstein,
//    "Lambert's W-Function",
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 1998.
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
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
    0.0000000000000000E+00,
    0.3517337112491958E+00,
    0.5671432904097839E+00,
    0.7258613577662263E+00,
    0.8526055020137255E+00,
    0.9585863567287029E+00,
    0.1000000000000000E+01,
    0.1049908894964040E+01,
    0.1130289326974136E+01,
    0.1202167873197043E+01,
    0.1267237814307435E+01,
    0.1326724665242200E+01,
    0.1381545379445041E+01,
    0.1432404775898300E+01,
    0.1479856830173851E+01,
    0.1524345204984144E+01,
    0.1566230953782388E+01,
    0.1605811996320178E+01,
    0.1745528002740699E+01,
    0.3385630140290050E+01,
    0.5249602852401596E+01,
    0.1138335808614005E+02 };
  static double x_vec[N_MAX] = {
    0.0000000000000000E+00,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.1500000000000000E+01,
    0.2000000000000000E+01,
    0.2500000000000000E+01,
    0.2718281828459045E+01,
    0.3000000000000000E+01,
    0.3500000000000000E+01,
    0.4000000000000000E+01,
    0.4500000000000000E+01,
    0.5000000000000000E+01,
    0.5500000000000000E+01,
    0.6000000000000000E+01,
    0.6500000000000000E+01,
    0.7000000000000000E+01,
    0.7500000000000000E+01,
    0.8000000000000000E+01,
    0.1000000000000000E+02,
    0.1000000000000000E+03,
    0.1000000000000000E+04,
    0.1000000000000000E+07 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
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
//****************************************************************************80

double wew_a ( double x, double &en )

//****************************************************************************80
//
//  Purpose:
//
//    WEW_A estimates Lambert's W function.
//
//  Discussion:
//
//    For a given X, this routine estimates the solution W of Lambert's 
//    equation:
//
//      X = W * EXP ( W )
//
//    This routine has higher accuracy than WEW_B.
//
//  Modified:
//
//    11 June 2014
//
//  Reference:
//
//    Fred Fritsch, R Shafer, W Crowley,
//    Algorithm 443: Solution of the transcendental equation w e^w = x,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 2, pages 123-124.
//
//  Parameters:
//
//    Input, double X, the argument of W(X)
//
//    Output, double &EN, the last relative correction to W(X).
//
//    Output, double WEW_A, the estimated value of W(X).
//
{
  const double c1 = 4.0 / 3.0;
  const double c2 = 7.0 / 3.0;
  const double c3 = 5.0 / 6.0;
  const double c4 = 2.0 / 3.0;
  double f;
  double temp;
  double temp2;
  double wn;
  double y;
  double zn;
//
//  Initial guess.
//
  f = log ( x );

  if ( x <= 6.46 )
  {
    wn = x * ( 1.0 + c1 * x ) / ( 1.0 + x * ( c2 + c3 * x ) );
    zn = f - wn - log ( wn );
  }
  else
  {
    wn = f;
    zn = - log ( wn );
  }
//
//  Iteration 1.
//
  temp = 1.0 + wn;
  y = 2.0 * temp * ( temp + c4 * zn ) - zn;
  wn = wn * ( 1.0 + zn * y / ( temp * ( y - zn ) ) );
//
//  Iteration 2.
//
  zn = f - wn - log ( wn );
  temp = 1.0 + wn;
  temp2 = temp + c4 * zn;
  en = zn * temp2 / ( temp * temp2 - 0.5 * zn );
  wn = wn * ( 1.0 + en );

  return wn;
}
//****************************************************************************80

double wew_b ( double x, double &en )

//****************************************************************************80
//
//  Purpose:
//
//    WEW_B estimates Lambert's W function.
//
//  Discussion:
//
//    For a given X, this routine estimates the solution W of Lambert's 
//    equation:
//
//      X = W * EXP ( W )
//
//    This routine has lower accuracy than WEW_A.
//
//  Modified:
//
//    11 June 2014
//
//  Reference:
//
//    Fred Fritsch, R Shafer, W Crowley,
//    Algorithm 443: Solution of the transcendental equation w e^w = x,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 2, pages 123-124.
//
//  Parameters:
//
//    Input, double X, the argument of W(X)
//
//    Output, double &EN, the last relative correction to W(X).
//
//    Output, double WEW_B, the estimated value of W(X).
//
{
  const double c1 = 4.0 / 3.0;
  const double c2 = 7.0 / 3.0;
  const double c3 = 5.0 / 6.0;
  const double c4 = 2.0 / 3.0;
  double f;
  double temp;
  double wn;
  double y;
  double zn;
//
//  Initial guess.
//
  f = log ( x );

  if ( x <= 0.7385 )
  {
    wn = x * ( 1.0 + c1 * x ) / ( 1.0 + x * ( c2 + c3 * x ) );
  }
  else
  {
    wn = f - 24.0 * ( ( f + 2.0 ) * f - 3.0 ) 
      / ( ( 0.7 * f + 58.0 ) * f + 127.0 );
  }
//
//  Iteration 1.
//
  zn = f - wn - log ( wn );
  temp = 1.0 + wn;
  y = 2.0 * temp * ( temp + c4 * zn ) - zn;
  en = zn * y / ( temp * ( y - zn ) );
  wn = wn * ( 1.0 + en );

  return wn;
}
