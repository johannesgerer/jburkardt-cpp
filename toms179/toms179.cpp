# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "toms179.hpp"

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
//    Malcolm Pike, David Hill
//    FORTRAN90 version by John Burkardt
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

void beta_cdf_values ( int *n_data, double *a, double *b, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CDF_VALUES returns some values of the Beta CDF.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
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

double mdbeta ( double x, double p, double q, int *ier )

//****************************************************************************80
//
//  Purpose:
//
//    MDBETA evaluates the incomplete beta function.
//
//  Modified:
//
//    30 January 2008
//
//  Author:
//
//    Oliver Ludwig
//    Modifications by John Burkardt
//
//  Reference:
//
//    Oliver Ludwig,
//    Algorithm 179:
//    Incomplete Beta Ratio,
//    Communications of the ACM,
//    Volume 6, Number 6, June 1963, page 314.
//
//  Parameters:
//
//    Input, double X, the value to which function is to be
//    integrated.  X must be in the range [0,1] inclusive.
//
//    Input, double P, the first parameter.  P must be greater
//    than 0.0.
//
//    Input, double Q, the second parameter.  Q must be greater
//    than 0.0.
//
//    Output, int *IER, error parameter.
//    0, normal exit.
//    1, X is not in the range [0,1] inclusive.
//    2, P or Q is less than or equal to 0.
//
//    Output, double MDBETA.  The probability that a random variable
//    from a Beta distribution having parameters P and Q will be less than
//    or equal to X.
//
//  Local parameters:
//
//    Local, double ALEPS, the logarithm of EPS1.
//
//    Local, double EPS, the machine precision.
//
//    Local, double EPS1, the smallest representable number.
//
{
  double aleps = - 179.6016;
  double c;
  double cnt;
  double d4;
  double dp;
  double dq;
  double eps = 2.2E-16;
  double eps1 = 1.0E-78;
  double finsum;
  int ib;
  int ifault;
  double infsum;
  int interval;
  double p1;
  double pq;
  double prob;
  double ps;
  double px;
  double temp;
  double wh;
  double xb;
  double y;
//
//  Check ranges of the arguments.
//
  prob = 0.0;
  y = x;

  if ( x < 0.0 || 1.0 < x )
  {
    *ier = 1;
    return prob;
  }

  if ( p <= 0.0 || q <= 0.0 )
  {
    *ier = 2;
    return prob;
  }

  *ier = 0;

  if ( x <= 0.5 )
  {
    interval = 0;
  }
  else
  {
    interval = 1;
    temp = p;
    p = q;
    q = temp;
    y = 1.0 - y;
  }

  if ( x == 0.0 || x == 1.0 )
  {
    prob = 0.0;

    if ( interval != 0 )
    {
      prob = 1.0 - prob;
      temp = p;
      p = q;
      q = temp;
    }
    return prob;
  }

  ib = q;
  temp = ib;
  ps = q - ( double ) ( ib );

  if ( q == temp )
  {
    ps = 1.0;
  }

  dp = p;
  dq = q;
  px = dp * log ( y );
  pq = alogam ( dp + dq, &ifault );
  p1 = alogam ( dp, &ifault );
  c = alogam ( dq, &ifault );
  d4 = log ( dp );
  xb = px + alogam ( ps + dp, &ifault ) - alogam ( ps, &ifault ) - d4 - p1;
//
//  Scaling
//
  ib = ( int ) ( xb / aleps );
  infsum = 0.0;
//
//  First term of a decreasing series will underflow.
//
  if ( ib == 0 )
  {
    infsum = exp ( xb );
    cnt = infsum * dp;
//
//  CNT will equal dexp ( temp ) * ( 1.d0 - ps ) * i * p * y**i / factorial ( i ).
//
    wh = 0.0;

    for ( ; ; )
    {
      wh = wh + 1.0;
      cnt = cnt * ( wh - ps ) * y / wh;
      xb = cnt / ( dp + wh );
      infsum = infsum + xb;

      if ( xb / eps < infsum )
      {
        break;
      }
    }
  }

  finsum = 0.0;

  if ( dq <= 1.0 )
  {
    prob = finsum + infsum;

    if ( interval != 0 )
    {
      prob = 1.0 - prob;
      temp = p;
      p = q;
      q = temp;
    }

    return prob;
  }

  xb = px + dq * log ( 1.0 - y ) + pq - p1 - log ( dq ) - c;
//
//  Scaling.
//
  ib = ( int ) ( xb / aleps );

  if ( ib < 0 )
  {
    ib = 0;
  }

  c = 1.0 / ( 1.0 - y );
  cnt = exp ( xb - ( double ) ( ib ) * aleps );
  ps = dq;
  wh = dq;

  for ( ; ; )
  {
    wh = wh - 1.0;

    if ( wh <= 0.0 )
    {
      prob = finsum + infsum;

      if ( interval != 0 )
      {
        prob = 1.0 - prob;
        temp = p;
        p = q;
        q = temp;
      }
      break;
    }

    px = ( ps * c ) / ( dp + wh );

    if ( px <= 1.0 )
    {
      if ( cnt / eps <= finsum || cnt <= eps1 / px )
      {
        prob = finsum + infsum;

        if ( interval != 0 )
        {
          prob = 1.0 - prob;
          temp = p;
          p = q;
          q = temp;
        }
        break;
      }
    }
    cnt = cnt * px;
//
//  Rescale.
//
    if ( 1.0 < cnt )
    {
      ib = ib - 1;
      cnt = cnt * eps1;
    }

    ps = wh;

    if ( ib == 0 )
    {
      finsum = finsum + cnt;
    }
  }

  return prob;
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
