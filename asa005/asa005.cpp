# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa005.hpp"

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

double prncst ( double st, int idf, double d, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    PRNCST computes the lower tail of noncentral T distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by BE Cooper.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    BE Cooper,
//    Algorithm AS 5:
//    The Integral of the Non-Central T-Distribution,
//    Applied Statistics,
//    Volume 17, Number 2, 1968, page 193.
//
//  Parameters:
//
//    Input, double ST, the argument.
//
//    Input, int IDF, the number of degrees of freedom.
//
//    Input, double D, the noncentrality parameter.
//
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double PRNCST, the value of the lower tail of
//    the noncentral T distribution.
//
//  Local Parameters:
//
//    Local, double G1, 1.0 / sqrt(2.0 * pi)
//
//    Local, double G2, 1.0 / (2.0 * pi)
//
//    Local, double G3, sqrt(2.0 * pi)
//
{
  double a;
  double ak;
  double b;
  double da;
  double drb;
  double emin = 12.5;
  double f;
  double fk;
  double fkm1;
  double fmkm1;
  double fmkm2;
  double g1 = 0.3989422804;
  double g2 = 0.1591549431;
  double g3 = 2.5066282746;
  int ioe;
  int k;
  double rb;
  double sum;
  double value;

  f = ( double ) ( idf );
//
//  For very large IDF, use the normal approximation.
//
  if ( 100 < idf )
  {
    *ifault = 1;

    a = sqrt ( 0.5 * f ) * exp ( lgamma ( 0.5 * ( f - 1.0 ) ) 
      - lgamma ( 0.5 * f ) ) * d;

    value = alnorm ( ( st - a ) / sqrt ( f * ( 1.0 + d * d ) 
      / ( f - 2.0 ) - a * a ), false );
    return value;
  }

  *ifault = 0;
  ioe = ( idf % 2 );
  a = st / sqrt ( f );
  b = f / ( f + st * st );
  rb = sqrt ( b );
  da = d * a;
  drb = d * rb;

  if ( idf == 1 )
  {
    value = alnorm ( drb, true ) + 2.0 * tfn ( drb, a );
    return value;
  }

  sum = 0.0;

  if ( fabs ( drb ) < emin )
  {
    fmkm2 = a * rb * exp ( - 0.5 * drb * drb ) 
    * alnorm ( a * drb, false ) * g1;
  }
  else
  {
    fmkm2 = 0.0;
  }

  fmkm1 = b * da * fmkm2;

  if ( fabs ( d ) < emin )
  {
    fmkm1 = fmkm1 + b * a * g2 * exp ( - 0.5 * d * d );
  }

  if ( ioe == 0 )
  {
    sum = fmkm2;
  }
  else
  {
    sum = fmkm1;
  }

  ak = 1.0;
  fk = 2.0;

  for ( k = 2; k <= idf - 2; k = k + 2 )
  {
    fkm1 = fk - 1.0;
    fmkm2 = b * ( da * ak * fmkm1 + fmkm2 ) * fkm1 / fk;
    ak = 1.0 / ( ak * fkm1 );
    fmkm1 = b * ( da * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0 );

    if ( ioe == 0 )
    {
      sum = sum + fmkm2;
    }
    else
    {
      sum = sum + fmkm1;
    }
    ak = 1.0 / ( ak * fk );
    fk = fk + 2.0;
  }

  if ( ioe == 0 )
  {
    value = alnorm ( d, true ) + sum * g3;
  }
  else
  {
    value = alnorm ( drb, true ) + 2.0 * ( sum + tfn ( drb, a ) );
  }

  return value;
}
//****************************************************************************80

void student_noncentral_cdf_values ( int *n_data, int *df, double *lambda, 
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NoncentralStudentTDistribution [ df, lambda ]
//      CDF [ dist, x ]
//
//    Mathematica seems to have some difficulty computing this function
//    to the desired number of digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2004
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
//    Output, int *DF, double *LAMBDA, the parameters of the
//    function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 30

  int df_vec[N_MAX] = { 
     1,  2,  3, 
     1,  2,  3, 
     1,  2,  3, 
     1,  2,  3, 
     1,  2,  3, 
    15, 20, 25, 
     1,  2,  3, 
    10, 10, 10, 
    10, 10, 10, 
    10, 10, 10 };

  double fx_vec[N_MAX] = { 
     0.8975836176504333E+00,  
     0.9522670169E+00,  
     0.9711655571887813E+00,  
     0.8231218864E+00,  
     0.9049021510E+00,  
     0.9363471834E+00,  
     0.7301025986E+00,  
     0.8335594263E+00,  
     0.8774010255E+00,  
     0.5248571617E+00,  
     0.6293856597E+00,  
     0.6800271741E+00,  
     0.20590131975E+00,  
     0.2112148916E+00,  
     0.2074730718E+00,  
     0.9981130072E+00,  
     0.9994873850E+00,  
     0.9998391562E+00,  
     0.168610566972E+00,  
     0.16967950985E+00,  
     0.1701041003E+00,  
     0.9247683363E+00,  
     0.7483139269E+00,  
     0.4659802096E+00,  
     0.9761872541E+00,  
     0.8979689357E+00,  
     0.7181904627E+00,  
     0.9923658945E+00,  
     0.9610341649E+00,  
     0.8688007350E+00 };

  double lambda_vec[N_MAX] = { 
     0.0E+00,  
     0.0E+00,  
     0.0E+00,  
     0.5E+00,  
     0.5E+00,  
     0.5E+00,  
     1.0E+00,  
     1.0E+00,  
     1.0E+00,  
     2.0E+00,  
     2.0E+00,  
     2.0E+00,  
     4.0E+00,  
     4.0E+00,  
     4.0E+00,  
     7.0E+00,  
     7.0E+00,  
     7.0E+00,  
     1.0E+00,  
     1.0E+00,  
     1.0E+00,  
     2.0E+00,  
     3.0E+00,  
     4.0E+00,  
     2.0E+00,  
     3.0E+00,  
     4.0E+00,  
     2.0E+00,  
     3.0E+00,  
     4.0E+00 };

  double x_vec[N_MAX] = { 
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
      3.00E+00,  
     15.00E+00,  
     15.00E+00,  
     15.00E+00,  
      0.05E+00,  
      0.05E+00,  
      0.05E+00,  
      4.00E+00,  
      4.00E+00,  
      4.00E+00,  
      5.00E+00,  
      5.00E+00,  
      5.00E+00,  
      6.00E+00,  
      6.00E+00,  
      6.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *df = 0;
    *lambda = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *df = df_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double tfn ( double x, double fx )

//****************************************************************************80
//
//  Purpose:
//
//    TFN calculates the T-function of Owen.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by JC Young, Christoph Minder.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    MA Porter, DJ Winstanley,
//    Remark AS R30:
//    A Remark on Algorithm AS76:
//    An Integral Useful in Calculating Noncentral T and Bivariate
//    Normal Probabilities,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, page 113.
//
//    JC Young, Christoph Minder,
//    Algorithm AS 76: 
//    An Algorithm Useful in Calculating Non-Central T and 
//    Bivariate Normal Distributions,
//    Applied Statistics,
//    Volume 23, Number 3, 1974, pages 455-457.
//
//  Parameters:
//
//    Input, double X, FX, the parameters of the function.
//
//    Output, double TFN, the value of the T-function.
//
{
# define NG 5

  double fxs;
  int i;
  double r[NG] = {
    0.1477621, 
    0.1346334, 
    0.1095432, 
    0.0747257, 
    0.0333357 };
  double r1;
  double r2;
  double rt;
  double tp = 0.159155;
  double tv1 = 1.0E-35;
  double tv2 = 15.0;
  double tv3 = 15.0;
  double tv4 = 1.0E-05;
  double u[NG] = {
    0.0744372, 
    0.2166977, 
    0.3397048,
    0.4325317, 
    0.4869533 };
  double value;
  double x1;
  double x2;
  double xs;
//
//  Test for X near zero.
//
  if ( fabs ( x ) < tv1 )
  {
    value = tp * atan ( fx );
    return value;
  }
//
//  Test for large values of abs(X).
//
  if ( tv2 < fabs ( x ) )
  {
    value = 0.0;
    return value;
  }
//
//  Test for FX near zero.
//
  if ( fabs ( fx ) < tv1 )
  {
    value = 0.0;
    return value;
  }
//
//  Test whether abs ( FX ) is so large that it must be truncated.
//
  xs = - 0.5 * x * x;
  x2 = fx;
  fxs = fx * fx;
//
//  Computation of truncation point by Newton iteration.
//
  if ( tv3 <= log ( 1.0 + fxs ) - xs * fxs )
  {
    x1 = 0.5 * fx;
    fxs = 0.25 * fxs;

    for ( ; ; )
    {
      rt = fxs + 1.0;

      x2 = x1 + ( xs * fxs + tv3 - log ( rt ) ) 
      / ( 2.0 * x1 * ( 1.0 / rt - xs ) );

      fxs = x2 * x2;

      if ( fabs ( x2 - x1 ) < tv4 )
      {
        break;
      }
      x1 = x2;
    }
  }
//
//  Gaussian quadrature.
//
  rt = 0.0;
  for ( i = 0; i < NG; i++ )
  {
    r1 = 1.0 + fxs * pow ( 0.5 + u[i], 2 );
    r2 = 1.0 + fxs * pow ( 0.5 - u[i], 2 );

    rt = rt + r[i] * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 );
  }

  value = rt * x2 * tp;

  return value;
# undef NG
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
