# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa109.hpp"

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

void beta_inc_values ( int *n_data, double *a, double *b, double *x, 
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 to X) T^(A-1) * (1-T)^(B-1) dT
//                      / Integral (0 to 1) T^(A-1) * (1-T)^(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0;
//      BETA_INC(A,B,1.0) = 1.0
//
//    The incomplete Beta function is also sometimes called the
//    "modified" Beta function, or the "normalized" Beta function
//    or the Beta CDF (cumulative density function.
//
//    In Mathematica, the function can be evaluated by:
//
//      BETA[X,A,B] / BETA[A,B]
//
//    The function can also be evaluated by using the Statistics package:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = BetaDistribution [ a, b ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2005
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
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
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
//    Output, double *A, B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 42

  double a_vec[N_MAX] = { 
      0.5E+00,  
      0.5E+00,  
      0.5E+00,  
      1.0E+00,  
      1.0E+00,  
      1.0E+00,  
      1.0E+00,  
      1.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      5.5E+00,  
     10.0E+00,  
     10.0E+00,  
     10.0E+00,  
     10.0E+00,  
     20.0E+00,  
     20.0E+00,  
     20.0E+00,  
     20.0E+00,  
     20.0E+00,  
     30.0E+00,  
     30.0E+00,  
     40.0E+00,
      0.1E+01, 
      0.1E+01, 
      0.1E+01, 
      0.1E+01, 
      0.1E+01, 
      0.1E+01, 
      0.1E+01, 
      0.1E+01, 
      0.2E+01, 
      0.3E+01, 
      0.4E+01, 
      0.5E+01 };

  double b_vec[N_MAX] = { 
      0.5E+00,  
      0.5E+00,  
      0.5E+00,  
      0.5E+00,  
      0.5E+00,  
      0.5E+00,  
      0.5E+00,  
      1.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      2.0E+00,  
      5.0E+00,  
      0.5E+00,  
      5.0E+00,  
      5.0E+00,  
     10.0E+00,  
      5.0E+00,  
     10.0E+00,  
     10.0E+00,  
     20.0E+00,  
     20.0E+00,  
     10.0E+00,  
     10.0E+00,  
     20.0E+00,
      0.5E+00, 
      0.5E+00, 
      0.5E+00, 
      0.5E+00, 
      0.2E+01, 
      0.3E+01, 
      0.4E+01, 
      0.5E+01, 
      0.2E+01, 
      0.2E+01, 
      0.2E+01, 
      0.2E+01 };

  double fx_vec[N_MAX] = { 
     0.6376856085851985E-01,  
     0.2048327646991335E+00,  
     0.1000000000000000E+01,  
     0.0000000000000000E+00,  
     0.5012562893380045E-02,  
     0.5131670194948620E-01,  
     0.2928932188134525E+00,  
     0.5000000000000000E+00,  
     0.2800000000000000E-01,  
     0.1040000000000000E+00,  
     0.2160000000000000E+00,  
     0.3520000000000000E+00,  
     0.5000000000000000E+00,  
     0.6480000000000000E+00,  
     0.7840000000000000E+00,  
     0.8960000000000000E+00,  
     0.9720000000000000E+00,  
     0.4361908850559777E+00,  
     0.1516409096347099E+00,  
     0.8978271484375000E-01,  
     0.1000000000000000E+01,  
     0.5000000000000000E+00,  
     0.4598773297575791E+00,  
     0.2146816102371739E+00,  
     0.9507364826957875E+00,  
     0.5000000000000000E+00,  
     0.8979413687105918E+00,  
     0.2241297491808366E+00,  
     0.7586405487192086E+00,  
     0.7001783247477069E+00,
     0.5131670194948620E-01, 
     0.1055728090000841E+00, 
     0.1633399734659245E+00, 
     0.2254033307585166E+00, 
     0.3600000000000000E+00, 
     0.4880000000000000E+00, 
     0.5904000000000000E+00, 
     0.6723200000000000E+00, 
     0.2160000000000000E+00, 
     0.8370000000000000E-01, 
     0.3078000000000000E-01, 
     0.1093500000000000E-01 };

  double x_vec[N_MAX] = { 
     0.01E+00,  
     0.10E+00,  
     1.00E+00,  
     0.00E+00,  
     0.01E+00,  
     0.10E+00,  
     0.50E+00,  
     0.50E+00,  
     0.10E+00,  
     0.20E+00,  
     0.30E+00,  
     0.40E+00,  
     0.50E+00,  
     0.60E+00,  
     0.70E+00,  
     0.80E+00,  
     0.90E+00,  
     0.50E+00,  
     0.90E+00,  
     0.50E+00,  
     1.00E+00,  
     0.50E+00,  
     0.80E+00,  
     0.60E+00,  
     0.80E+00,  
     0.50E+00,  
     0.60E+00,  
     0.70E+00,  
     0.80E+00,  
     0.70E+00,
     0.10E+00, 
     0.20E+00, 
     0.30E+00, 
     0.40E+00, 
     0.20E+00, 
     0.20E+00, 
     0.20E+00, 
     0.20E+00, 
     0.30E+00, 
     0.30E+00, 
     0.30E+00, 
     0.30E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0;
    *b = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double betain ( double x, double p, double q, double beta, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    BETAIN computes the incomplete Beta function ratio.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    KL Majumder, GP Bhattacharjee,
//    Algorithm AS 63:
//    The incomplete Beta Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 409-411.
//
//  Parameters:
//
//    Input, double X, the argument, between 0 and 1.
//
//    Input, double P, Q, the parameters, which
//    must be positive.
//
//    Input, double BETA, the logarithm of the complete
//    beta function.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    nonzero, an error occurred.
//
//    Output, double BETAIN, the value of the incomplete
//    Beta function ratio.
//
{
  double acu = 0.1E-14;
  double ai;
  double betain;
  double cx;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double value;
  double xx;

  value = x;
  *ifault = 0;
//
//  Check the input arguments.
//
  if ( p <= 0.0 || q <= 0.0 )
  {
    *ifault = 1;
    return value;
  }

  if ( x < 0.0 || 1.0 < x )
  {
    *ifault = 2;
    return value;
  }
//
//  Special cases.
//
  if ( x == 0.0 || x == 1.0 )
  {
    return value;
  }
//
//  Change tail if necessary and determine S.
//
  psq = p + q;
  cx = 1.0 - x;

  if ( p < psq * x )
  {
    xx = cx;
    cx = x;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    xx = x;
    pp = p;
    qq = q;
    indx = false;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = ( int ) ( qq + cx * psq );
//
//  Use the Soper reduction formula.
//
  rx = xx / cx;
  temp = qq - ai;
  if ( ns == 0 )
  {
    rx = xx;
  }

  for ( ; ; )
  {
    term = term * temp * rx / ( pp + ai );
    value = value + term;;
    temp = r8_abs ( term );

    if ( temp <= acu && temp <= acu * value )
    {
      value = value * exp ( pp * log ( xx ) 
      + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

      if ( indx )
      {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ai;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
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

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
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
//****************************************************************************80

double xinbta ( double p, double q, double beta, double alpha, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    XINBTA computes inverse of the incomplete Beta function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    GW Cran, KJ Martin, GE Thomas,
//    Remark AS R19 and Algorithm AS 109:
//    A Remark on Algorithms AS 63: The Incomplete Beta Integral
//    and AS 64: Inverse of the Incomplete Beta Integeral,
//    Applied Statistics,
//    Volume 26, Number 1, 1977, pages 111-114.
//
//  Parameters:
//
//    Input, double P, Q, the parameters of the incomplete
//    Beta function.
//
//    Input, double BETA, the logarithm of the value of
//    the complete Beta function.
//
//    Input, double ALPHA, the value of the incomplete Beta
//    function.  0 <= ALPHA <= 1.
//
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double XINBTA, the argument of the incomplete
//    Beta function which produces the value ALPHA.
//
//  Local Parameters:
//
//    Local, double SAE, the most negative decimal exponent
//    which does not cause an underflow.
//
{
  double a;
  double acu;
  double adj;
  double fpu;
  double g;
  double h;
  int iex;
  bool indx;
  double pp;
  double prev;
  double qq;
  double r;
  double s;
  double sae = -37.0;
  double sq;
  double t;
  double tx;
  double value;
  double w;
  double xin;
  double y;
  double yprev;

  fpu = pow ( 10.0, sae );

  *ifault = 0;
  value = alpha;
//
//  Test for admissibility of parameters.
//
  if ( p <= 0.0 || q <= 0.0 )
  {
    *ifault = 1;
    return value;
  }

  if ( alpha < 0.0 || 1.0 < alpha )
  {
    *ifault = 2;
    return value;
  }

  if ( alpha == 0.0 || alpha == 1.0 )
  {
    return value;
  }
//
//  Change tail if necessary.
//
  if ( 0.5 < alpha )
  {
    a = 1.0 - alpha;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    a = alpha;
    pp = p;
    qq = q;
    indx = false;
  }
//
//  Calculate the initial approximation.
//
  r = sqrt ( - log ( a * a ) );

  y = r - ( 2.30753 + 0.27061 * r ) 
    / ( 1.0 + ( 0.99229 + 0.04481 * r ) * r );

  if ( 1.0 < pp && 1.0 < qq )
  {
    r = ( y * y - 3.0 ) / 6.0;
    s = 1.0 / ( pp + pp - 1.0 );
    t = 1.0 / ( qq + qq - 1.0 );
    h = 2.0 / ( s + t );
    w = y * sqrt ( h + r ) / h - ( t - s ) 
      * ( r + 5.0 / 6.0 - 2.0 / ( 3.0 * h ) );
    value = pp / ( pp + qq * exp ( w + w ) );
  }
  else
  {
    r = qq + qq;
    t = 1.0 / ( 9.0 * qq );
    t = r * pow ( 1.0 - t + y * sqrt ( t ), 3 );

    if ( t <= 0.0 )
    {
      value = 1.0 - exp ( ( log ( ( 1.0 - a ) * qq ) + beta ) / qq );
    }
    else
    {
      t = ( 4.0 * pp + r - 2.0 ) / t;

      if ( t <= 1.0 )
      {
        value = exp ( ( log ( a * pp ) + beta ) / pp );
      }
      else
      {
        value = 1.0 - 2.0 / ( t + 1.0 );
      }
    }
  }
//
//  Solve for X by a modified Newton-Raphson method,
//  using the function BETAIN.
//
  r = 1.0 - pp;
  t = 1.0 - qq;
  yprev = 0.0;
  sq = 1.0;
  prev = 1.0;

  if ( value < 0.0001 )
  {
    value = 0.0001;
  }

  if ( 0.9999 < value )
  {
    value = 0.9999;
  }

  iex = r8_max ( - 5.0 / pp / pp - 1.0 / pow ( a, 0.2 ) - 13.0, sae );

  acu = pow ( 10.0, iex );

  for ( ; ; )
  {
    y = betain ( value, pp, qq, beta, ifault );

    if ( *ifault != 0 )
    {
      *ifault = 3;
      return value;
    }

    xin = value;
    y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.0 - xin ) );

    if ( y * yprev <= 0.0 )
    {
      prev = r8_max ( sq, fpu );
    }

    g = 1.0;

    for ( ; ; )
    {
      for ( ; ; )
      {
        adj = g * y;
        sq = adj * adj;

        if ( sq < prev )
        {
          tx = value - adj;

          if ( 0.0 <= tx && tx <= 1.0 )
          {
            break;
          }
        }
        g = g / 3.0;
      }

      if ( prev <= acu )
      {
        if ( indx )
        {
          value = 1.0 - value;
        }
        return value;
      }

      if ( y * y <= acu )
      {
        if ( indx )
        {
          value = 1.0 - value;
        }
        return value;
      }

      if ( tx != 0.0 && tx != 1.0 )
      {
        break;
      }

      g = g / 3.0;
    }

    if ( tx == value )
    {
      break;
    }

    value = tx;
    yprev = y;
  }

  if ( indx )
  {
    value = 1.0 - value;
  }

  return value;
}
