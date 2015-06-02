# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa109.hpp"

//****************************************************************************80

void beta_inc_values ( int &n_data, double &a, double &b, double &x,
  double &fx )

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
//      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
//                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
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
//    25 September 2014
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
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 45

  static double a_vec[N_MAX] = {
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
      0.5E+01,
      1.30625,
      1.30625,
      1.30625 };

  static double b_vec[N_MAX] = {
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
      0.2E+01,
     11.7562, 
     11.7562, 
     11.7562 };

  static double fx_vec[N_MAX] = {
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
     0.1093500000000000E-01,
     0.918884684620518,
     0.21052977489419,
     0.1824130512500673 };

  static double x_vec[N_MAX] = {
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
     0.30E+00,
     0.225609,
     0.0335568,
     0.0295222 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double betain ( double x, double p, double q, double beta, int &ifault )

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
//    25 September 2014
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
//    Output, int &IFAULT, error flag.
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
  ifault = 0;
//
//  Check the input arguments.
//
  if ( p <= 0.0 || q <= 0.0 )
  {
    cerr << "\n";
    cerr << "BETAIN - Fatal error!\n";
    cerr << "  P <= 0.0 or Q <= 0.0\n";
    ifault = 1;
    exit ( 1 );
  }

  if ( x < 0.0 || 1.0 < x )
  {
    cerr << "\n";
    cerr << "BETAIN - Fatal error!\n";
    cerr << "  X < 0.0 or 1 < X\n";
    ifault = 2;
    exit ( 1 );
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
    temp = fabs ( term );

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
//****************************************************************************80

double xinbta ( double p, double q, double beta, double alpha, int &ifault )

//****************************************************************************80
//
//  Purpose:
//
//    XINBTA computes inverse of the incomplete Beta function.
//
//  Discussion:
//
//    The accuracy exponent SAE was loosened from -37 to -30, because
//    the code would not otherwise accept the results of an iteration
//    with p = 0.3, q = 3.0, alpha = 0.2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2014
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
//    Output, int &IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double XINBTA, the argument of the incomplete
//    Beta function which produces the value ALPHA.
//
//  Local Parameters:
//
//    Local, double SAE, requests an accuracy of about 10^SAE.
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
  double sae = -30.0;
  double sq;
  double t;
  double tx;
  double value;
  double w;
  double xin;
  double y;
  double yprev;

  fpu = pow ( 10.0, sae );

  ifault = 0;
  value = alpha;
//
//  Test for admissibility of parameters.
//
  if ( p <= 0.0 )
  {
    cerr << "\n";
    cerr << "XINBTA - Fatal error!\n";
    cerr << "  P <= 0.0.\n";
    ifault = 1;
    exit ( 1 );
  }

  if ( q <= 0.0 )
  {
    cerr << "\n";
    cerr << "XINBTA - Fatal error!\n";
    cerr << "  Q <= 0.0.\n";
    ifault = 1;
    exit ( 1 );
  }

  if ( alpha < 0.0 || 1.0 < alpha )
  {
    cerr << "\n";
    cerr << "XINBTA - Fatal error!\n";
    cerr << "  ALPHA not between 0 and 1.\n";
    ifault = 2;
    exit ( 1 );
  }
//
//  If the answer is easy to determine, return immediately.
//
  if ( alpha == 0.0 )
  {
    value = 0.0;
    return value;
  }

  if ( alpha == 1.0 )
  {
    value = 1.0;
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
//
//  Iteration loop.
//
  for ( ; ; )
  {
    y = betain ( value, pp, qq, beta, ifault );

    if ( ifault != 0 )
    {
      cerr << "\n";
      cerr << "XINBTA - Fatal error!\n";
      cerr << "  BETAIN returned IFAULT = " << ifault << "\n";
      ifault = 1;
      exit ( 1 );
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
//
//  Check whether the current estimate is acceptable.
//  The change "VALUE = TX" was suggested by Ivan Ukhov.
//
      if ( prev <= acu || y * y <= acu )
      {
        value = tx;
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
