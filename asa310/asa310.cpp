# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa310.hpp"

//****************************************************************************80

double alnorm ( double x, bool upper )

//****************************************************************************80
//
//  Purpose:
//
//    ALNORM computes the cumulative density of the standard normal distribution.
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
//    Original FORTRAN77 version by David Hill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Hill,
//    Algorithm AS 66:
//    The Normal Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 424-427.
//
//  Parameters:
//
//    Input, double X, is one endpoint of the semi-infinite interval
//    over which the integration takes place.
//
//    Input, bool UPPER, determines whether the upper or lower
//    interval is to be integrated:
//    .TRUE.  => integrate from X to + Infinity;
//    .FALSE. => integrate from - Infinity to X.
//
//    Output, double ALNORM, the integral of the standard normal
//    distribution over the desired interval.
//
{
  double a1 = 5.75885480458;
  double a2 = 2.62433121679;
  double a3 = 5.92885724438;
  double b1 = -29.8213557807;
  double b2 = 48.6959930692;
  double c1 = -0.000000038052;
  double c2 = 0.000398064794;
  double c3 = -0.151679116635;
  double c4 = 4.8385912808;
  double c5 = 0.742380924027;
  double c6 = 3.99019417011;
  double con = 1.28;
  double d1 = 1.00000615302;
  double d2 = 1.98615381364;
  double d3 = 5.29330324926;
  double d4 = -15.1508972451;
  double d5 = 30.789933034;
  double ltone = 7.0;
  double p = 0.398942280444;
  double q = 0.39990348504;
  double r = 0.398942280385;
  bool up;
  double utzero = 18.66;
  double value;
  double y;
  double z;

  up = upper;
  z = x;

  if ( z < 0.0 )
  {
    up = !up;
    z = - z;
  }

  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  y = 0.5 * z * z;

  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y
      / ( y + a1 + b1
      / ( y + a2 + b2
      / ( y + a3 ))));
  }
  else
  {
    value = r * exp ( - y )
      / ( z + c1 + d1
      / ( z + c2 + d2
      / ( z + c3 + d3
      / ( z + c4 + d4
      / ( z + c5 + d5
      / ( z + c6 ))))));
  }

  if ( !up )
  {
    value = 1.0 - value;
  }

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

double betanc ( double x, double a, double b, double lambda, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    BETANC computes the tail of the noncentral Beta distribution.
//
//  Discussion:
//
//    This routine returns the cumulative probability of X for the non-central
//    Beta distribution with parameters A, B and non-centrality LAMBDA.
//
//    Note that if LAMBDA = 0, the standard Beta distribution is defined.
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
//    Original FORTRAN77 version by Russell Lenth.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Russell Lenth,
//    Algorithm AS 226:
//    Computing Noncentral Beta Probabilities,
//    Applied Statistics,
//    Volume 36, Number 2, 1987, pages 241-244.
//
//    H Frick,
//    Algorithm AS R84:
//    A Remark on Algorithm AS 226:
//    Computing Noncentral Beta Probabilities,
//    Applied Statistics,
//    Volume 39, Number 2, 1990, pages 311-312.
//
//  Parameters:
//
//    Input, double X, the value defining the cumulative
//    probability lower tail.  Normally, 0 <= X <= 1, but any value
//    is allowed.
//
//    Input, double A, B, the parameters of the distribution.
//    0 < A, 0 < B.
//
//    Input, double LAMBDA, the noncentrality parameter
//    of the distribution.  0 <= LAMBDA.  The program can produce reasonably
//    accurate results for values of LAMBDA up to about 100.
//
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double BETANC, the cumulative probability
//    of X.
//
{
  double a0;
  double ax;
  double beta;
  double c;
  double errbd;
  double errmax = 1.0E-07;
  double gx;
  int itrmax = 150;
  double q;
  double sumq;
  double temp;
  double ualpha = 5.0;
  double value;
  double x0;
  double xj;

  *ifault = 0;

  if ( lambda < 0.0 ||
       a <= 0.0 ||
       b <= 0.0 )
  {
    *ifault = 2;
    value = -1.0;
    return value;
  }

  if ( x <= 0.0 )
  {
    value = 0.0;
    return value;
  }

  if ( 1.0 <= x )
  {
    value = 1.0;
    return value;
  }

  c = 0.5 * lambda;
//
//  Initialize the series.
//
  beta = lgamma ( a )
       + lgamma ( b )
       - lgamma ( a + b );

  temp = betain ( x, a, b, beta, ifault );

  gx = exp ( a * log ( x ) + b * log ( 1.0 - x ) - beta - log ( a ) );

  q = exp ( - c );

  xj = 0.0;
  ax = q * temp;
  sumq = 1.0 - q;
  value = ax;
//
//  Recur over subsequent terms until convergence is achieved.
//
  *ifault = 1;

  for ( ; ; )
  {
    xj = xj + 1.0;
    temp = temp - gx;
    gx = x * ( a + b + xj - 1.0 ) * gx / ( a + xj );
    q = q * c / xj;
    sumq = sumq - q;
    ax = temp * q;
    value = value + ax;
//
//  Check for convergence and act accordingly.
//
    errbd = fabs ( ( temp - gx ) * sumq );

    if ( errbd <= errmax )
    {
      *ifault = 0;
      break;
    }

    if (  itrmax < ( int ) xj )
    {
      break;
    }
  }
  return value;
}
//****************************************************************************80

double gammad ( double x, double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMAD computes the Incomplete Gamma Integral
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by B Shea.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    B Shea,
//    Algorithm AS 239:
//    Chi-squared and Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 37, Number 3, 1988, pages 466-473.
//
//  Parameters:
//
//    Input, double X, P, the parameters of the incomplete
//    gamma ratio.  0 <= X, and 0 < P.
//
//    Output, int IFAULT, error flag.
//    0, no error.
//    1, X < 0 or P <= 0.
//
//    Output, double GAMMAD, the value of the incomplete
//    Gamma integral.
//
{
  double a;
  double an;
  double arg;
  double b;
  double c;
  double elimit = - 88.0;
  double oflo = 1.0E+37;
  double plimit = 1000.0;
  double pn1;
  double pn2;
  double pn3;
  double pn4;
  double pn5;
  double pn6;
  double rn;
  double tol = 1.0E-14;
  bool upper;
  double value;
  double xbig = 1.0E+08;

  value = 0.0;
//
//  Check the input.
//
  if ( x < 0.0 )
  {
    *ifault = 1;
    return value;
  }

  if ( p <= 0.0 )
  {
    *ifault = 1;
    return value;
  }

  *ifault = 0;

  if ( x == 0.0 )
  {
    value = 0.0;
    return value;
  }
//
//  If P is large, use a normal approximation.
//
  if ( plimit < p )
  {
    pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 )
    + 1.0 / ( 9.0 * p ) - 1.0 );

    upper = false;
    value = alnorm ( pn1, upper );
    return value;
  }
//
//  If X is large set value = 1.
//
  if ( xbig < x )
  {
    value = 1.0;
    return value;
  }
//
//  Use Pearson's series expansion.
//
  if ( x <= 1.0 || x < p )
  {
    arg = p * log ( x ) - x - lgamma ( p + 1.0 );
    c = 1.0;
    value = 1.0;
    a = p;

    for ( ; ; )
    {
      a = a + 1.0;
      c = c * x / a;
      value = value + c;

      if ( c <= tol )
      {
        break;
      }
    }

    arg = arg + log ( value );

    if ( elimit <= arg )
    {
      value = exp ( arg );
    }
    else
    {
      value = 0.0;
    }
  }
//
//  Use a continued fraction expansion.
//
  else
  {
    arg = p * log ( x ) - x - lgamma ( p );
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    for ( ; ; )
    {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      an = a * c;
      pn5 = b * pn3 - an * pn1;
      pn6 = b * pn4 - an * pn2;

      if ( pn6 != 0.0 )
      {
        rn = pn5 / pn6;

        if ( fabs ( value - rn ) <= r8_min ( tol, tol * rn ) )
        {
          break;
        }
        value = rn;
      }

      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
//
//  Re-scale terms in continued fraction if terms are large.
//
      if ( oflo <= abs ( pn5 ) )
      {
        pn1 = pn1 / oflo;
        pn2 = pn2 / oflo;
        pn3 = pn3 / oflo;
        pn4 = pn4 / oflo;
      }
    }

    arg = arg + log ( value );

    if ( elimit <= arg )
    {
      value = 1.0 - exp ( arg );
    }
    else
    {
      value = 1.0;
    }
  }

  return value;
}
//****************************************************************************80

double ncbeta ( double a, double b, double lambda, double x, double errmax,
  int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NCBETA computes the noncentral Beta CDF.
//
//  Discussion:
//
//    Three corrections needed to be made to the text of this routine.
//    They are noted in the comments below.
//
//    Two of these corrections were errors in transcription made when
//    producing the online copy distributed by APSTAT.
//
//    One error, an error of omission, occurred in the original printed
//    copy of the routine, and was carried over into the online copy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R Chattamvelli, R Shanmugam.
//    C++ version by John Burkardt.
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
//    Input, double A, B, the shape parameters.
//    0 <= A, 0 <= B.
//
//    Input, double LAMBDA, the noncentrality parameter.
//    0 <= LAMBDA.
//
//    Input, double X, the value at which the CDF is desired.
//
//    Input, double ERRMAX, the precision tolerance.
//
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    1, X is 0 or 1.
//    2, X < 0 or 1 < X.
//    3, A, B or LAMBDA is less than 0.
//
//    Output, double NCBETA, the value of the noncentral Beta CDF.
//
{
  double beta;
  double c;
  double ebd;
  double errbd;
  double ftemp;
  double fx;
  double gx;
  int i;
  int iter1;
  int iter2;
  int iterhi;
  int iterlo;
  int j;
  int m;
  double mr;
  double psum;
  double q;
  double r;
  double s;
  double s0;
  double s1;
  double sum;
  double t;
  double t0;
  double t1;
  double temp;
  double value;
  int xj;

  *ifault = 0;
  value = x;
//
//  Check parameters.
//
  if ( lambda <= 0.0 )
  {
    *ifault = 3;
    return value;
  }

  if ( a <= 0.0 )
  {
    *ifault = 3;
    return value;
  }

  if ( b <= 0.0 )
  {
    *ifault = 3;
    return value;
  }

  if ( x <= 0.0 )
  {
    value = 0.0;
    return value;
  }

  if ( 1.0 <= x )
  {
    value = 1.0;
    return value;
  }

  c = 0.5 * lambda;
  xj = 0.0;
//
//  AS 226 as it stands is sufficient in this situation.
//
  if ( lambda < 54.0 )
  {
    value = betanc ( x, a, b, lambda, ifault );
    return value;
  }
  else
  {
    m = ( int ) ( c + 0.5 );
    mr = ( double ) ( m );
    iterlo = m - ( int ) ( 5.0 * sqrt ( mr ) );
    iterhi = m + ( int ) ( 5.0 * sqrt ( mr ) );
    t = - c + mr * log ( c ) - lgamma ( mr + 1.0 );
    q = exp ( t );
    r = q;
    psum = q;

    beta = lgamma ( a + mr )
         + lgamma ( b )
         - lgamma ( a + mr + b );

    s1 = ( a + mr ) * log ( x )
       + b * log ( 1.0 - x ) - log ( a + mr ) - beta;
    gx = exp ( s1 );
    fx = gx;
    temp = betain ( x, a + mr, b, beta, ifault );
    ftemp = temp;
    xj = xj + 1.0;
//
//  The online copy of AS 310 has "SUM = Q - TEMP" which is incorrect.
//
    sum = q * temp;
    iter1 = m;
//
//  The first set of iterations starts from M and goes downwards
//
    for ( ; ; )
    {
      if ( iter1 < iterlo )
      {
        break;
      }

      if ( q < errmax )
      {
        break;
      }
//
//  The online copy of AS 310 has "Q = Q - ITER1 / C" which is incorrect.
//
      q = q * iter1 / c;
      xj = xj + 1.0;
      gx = ( a + iter1 ) / ( x * ( a + b + iter1 - 1.0 ) ) * gx;
      iter1 = iter1 - 1;
      temp = temp + gx;
      psum = psum + q;
      sum = sum + q * temp;
    }

    t0 = lgamma ( a + b )
       - lgamma ( a + 1.0 )
       - lgamma ( b );

    s0 = a * log ( x ) + b * log ( 1.0 - x );
//
//  Both the online copy of AS 310 and the text printed in the reference
//  did not initialize the variable S to zero, which is incorrect.
//  JVB, 12 January 2008.
//
    s = 0.0;
    for ( i = 1; i <= iter1; i++ )
    {
      j = i - 1;
      s = s + exp ( t0 + s0 + j * log ( x ) );
      t1 = log ( a + b + j ) - log ( a + 1.0 + j ) + t0;
      t0 = t1;
    }
//
//  Compute the first part of error bound.
//
    errbd = ( 1.0 - gammad ( c, ( double ) ( iter1 ), ifault ) ) * ( temp + s );

    q = r;
    temp = ftemp;
    gx = fx;
    iter2 = m;

    for ( ; ; )
    {
      ebd = errbd + ( 1.0 - psum ) * temp;

      if ( ebd < errmax || iterhi <= iter2 )
      {
        break;
      }
      iter2 = iter2 + 1;
      xj = xj + 1.0;
      q = q * c / iter2;
      psum = psum + q;
      temp = temp - gx;
      gx = x * ( a + b + iter2 - 1.0 ) / ( a + iter2 ) * gx;
      sum = sum + q * temp;
    }
  }

  value = sum;

  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
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
