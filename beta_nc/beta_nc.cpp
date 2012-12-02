# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "beta_nc.H"

//****************************************************************************80

double alogam ( double x, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    ALOGAM computes the logarithm of the Gamma function.
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

double beta_noncentral_cdf ( double a, double b, double lambda, double x,
  double error_max )

//****************************************************************************80
//
//  Purpose:
//
//  BETA_NONCENTRAL_CDF evaluates the noncentral Beta CDF.
//
//  Discussion:
//
//    The reference mistakenly phrases the opposite of the correct
//    stopping criterion, that is, it says:
//
//      "stop when PSUM < 1 - ERROR"
//
//    but this must be:
//
//      "stop when 1 - ERROR < PSUM."
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harry Posten,
//    An Effective Algorithm for the Noncentral Beta Distribution Function,
//    The American Statistician,
//    Volume 47, Number 2, May 1993, pages 129-131.
//
//  Parameters:
//
//    Input, double A, B, the shape parameters.
//
//    Input, double LAMBDA, the noncentrality parameter.
//
//    Input, double X, the argument of the function.
//
//    Input, double ERROR_MAX, the error control.
//
//    Output, double BETA_NONCENTRAL_CDF, the value of the noncentral Beta CDF.
//
{
  double beta_log;
  double bi;
  double bj;
  int i;
  int ifault;
  double p_sum;
  double pb_sum;
  double pi;
  double pj;
  double si;
  double sj;
  double value;

  i = 0;
  pi = exp ( - lambda / 2.0 );

  beta_log = alogam ( a, &ifault )
           + alogam ( b, &ifault )
           - alogam ( a + b, &ifault );

  bi = betain ( x, a, b, beta_log, &ifault );

  si = exp (
      a * log ( x )
    + b * log ( 1.0 - x )
    - beta_log
    - log ( a ) );

  p_sum = pi;
  pb_sum = pi * bi;

  while ( p_sum < 1.0 - error_max )
  {
    pj = pi;
    bj = bi;
    sj = si;

    i = i + 1;
    pi = 0.5 * lambda * pj / ( double ) ( i );
    bi = bj - sj;
    si = x * ( a + b + i - 1 ) * sj / ( a + i );

    p_sum = p_sum + pi;
    pb_sum = pb_sum + pi * bi;
  }

  value = pb_sum;

  return value;
}
//****************************************************************************80

void beta_noncentral_cdf_values ( int *n_data, double *a, double *b,
  double *lambda, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
//
//  Discussion:
//
//    The values presented here are taken from the reference, where they
//    were given to a limited number of decimal places.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R Chattamvelli, R Shanmugam,
//    Algorithm AS 310:
//    Computing the Non-central Beta Distribution Function,
//    Applied Statistics,
//    Volume 46, Number 1, 1997, pages 146-156.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, *B, the shape parameters.
//
//    Output, double *LAMBDA, the noncentrality parameter.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 25

  double a_vec[N_MAX] = {
        5.0,
        5.0,
        5.0,
       10.0,
       10.0,
       10.0,
       20.0,
       20.0,
       20.0,
       10.0,
       10.0,
       15.0,
       20.0,
       20.0,
       20.0,
       30.0,
       30.0,
       10.0,
       10.0,
       10.0,
       15.0,
       10.0,
       12.0,
       30.0,
       35.0 };
  double b_vec[N_MAX] = {
        5.0,
        5.0,
        5.0,
       10.0,
       10.0,
       10.0,
       20.0,
       20.0,
       20.0,
       20.0,
       10.0,
        5.0,
       10.0,
       30.0,
       50.0,
       20.0,
       40.0,
        5.0,
       10.0,
       30.0,
       20.0,
        5.0,
       17.0,
       30.0,
       30.0 };
  double fx_vec[N_MAX] = {
       0.4563021,
       0.1041337,
       0.6022353,
       0.9187770,
       0.6008106,
       0.0902850,
       0.9998655,
       0.9925997,
       0.9641112,
       0.9376626573,
       0.7306817858,
       0.1604256918,
       0.1867485313,
       0.6559386874,
       0.9796881486,
       0.1162386423,
       0.9930430054,
       0.0506899273,
       0.1030959706,
       0.9978417832,
       0.2555552369,
       0.0668307064,
       0.0113601067,
       0.7813366615,
       0.8867126477 };
  double lambda_vec[N_MAX] = {
        54.0,
       140.0,
       170.0,
        54.0,
       140.0,
       250.0,
        54.0,
       140.0,
       250.0,
       150.0,
       120.0,
        80.0,
       110.0,
        65.0,
       130.0,
        80.0,
       130.0,
        20.0,
        54.0,
        80.0,
       120.0,
        55.0,
        64.0,
       140.0,
        20.0 };
  double x_vec[N_MAX] = {
       0.8640,
       0.9000,
       0.9560,
       0.8686,
       0.9000,
       0.9000,
       0.8787,
       0.9000,
       0.9220,
       0.868,
       0.900,
       0.880,
       0.850,
       0.660,
       0.720,
       0.720,
       0.800,
       0.644,
       0.700,
       0.780,
       0.760,
       0.795,
       0.560,
       0.800,
       0.670 };

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
    *lambda = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
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
