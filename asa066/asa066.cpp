# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa066.hpp"

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

void normp ( double z, double *p, double *q, double *pdf )

//****************************************************************************80
//
//  Purpose:
//
//    NORMP computes the cumulative density of the standard normal distribution.
//
//  Discussion:
//
//    This is algorithm 5666 from Hart, et al.
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
//    Original FORTRAN77 version by Alan Miller.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
//    Charles Mesztenyi, John Rice, Henry Thacher, 
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double Z, divides the real line into two 
//    semi-infinite intervals, over each of which the standard normal 
//    distribution is to be integrated.
//
//    Output, double *P, *Q, the integrals of the standard normal
//    distribution over the intervals ( - Infinity, Z] and 
//    [Z, + Infinity ), respectively.
//
//    Output, double *PDF, the value of the standard normal distribution
//    at Z.
//
{
  double cutoff = 7.071;
  double expntl;
  double p0 = 220.2068679123761;
  double p1 = 221.2135961699311;
  double p2 = 112.0792914978709;
  double p3 = 33.91286607838300;
  double p4 = 6.373962203531650;
  double p5 = 0.7003830644436881;
  double p6 = 0.03526249659989109;
  double q0 = 440.4137358247522;
  double q1 = 793.8265125199484;
  double q2 = 637.3336333788311;
  double q3 = 296.5642487796737;
  double q4 = 86.78073220294608;
  double q5 = 16.06417757920695;
  double q6 = 1.755667163182642;
  double q7 = 0.08838834764831844;
  double root2pi = 2.506628274631001;
  double zabs;

  zabs = fabs ( z );
//
//  37 < |Z|.
//
  if ( 37.0 < zabs )
  {
    *pdf = 0.0;
    *p = 0.0;
  }
//
//  |Z| <= 37.
//
  else
  {
    expntl = exp ( - 0.5 * zabs * zabs );
    *pdf = expntl / root2pi;
//
//  |Z| < CUTOFF = 10 / sqrt(2).
//
    if ( zabs < cutoff )
    {
      *p = expntl * (((((( 
          p6   * zabs 
        + p5 ) * zabs 
        + p4 ) * zabs 
        + p3 ) * zabs 
        + p2 ) * zabs 
        + p1 ) * zabs 
        + p0 ) / ((((((( 
          q7   * zabs 
        + q6 ) * zabs 
        + q5 ) * zabs 
        + q4 ) * zabs 
        + q3 ) * zabs 
        + q2 ) * zabs 
        + q1 ) * zabs 
      + q0 );
    }
//
//  CUTOFF <= |Z|.
//
    else
    {
      *p = *pdf / ( 
        zabs + 1.0 / ( 
        zabs + 2.0 / ( 
        zabs + 3.0 / ( 
        zabs + 4.0 / ( 
        zabs + 0.65 )))));
    }
  }

  if ( z < 0.0 )
  {
    *q = 1.0 - *p;
  }
  else
  {
    *q = *p;
    *p = 1.0 - *q;
  }

  return;
}
//****************************************************************************80

void nprob ( double z, double *p, double *q, double *pdf )

//****************************************************************************80
//
//  Purpose:
//
//    NPROB computes the cumulative density of the standard normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by AG Adams.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    AG Adams,
//    Algorithm 39:
//    Areas Under the Normal Curve,
//    Computer Journal,
//    Volume 12, Number 2, May 1969, pages 197-198.
//
//  Parameters:
//
//    Input, double Z, divides the real line into 
//    two semi-infinite intervals, over each of which the standard normal 
//    distribution is to be integrated.
//
//    Output, double *P, *Q, the integrals of the standard normal
//    distribution over the intervals ( - Infinity, Z] and 
//    [Z, + Infinity ), respectively.
//
//    Output, double *PDF, the value of the standard normal
//    distribution at Z.
//
{
  double a0 = 0.5;
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 0.000000038052;
  double b2 = 1.00000615302;
  double b3 = 0.000398064794;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double y;
  double zabs;

  zabs = fabs ( z );
//
//  |Z| between 0 and 1.28
//
  if ( zabs <= 1.28 )
  {
    y = a0 * z * z;
    *pdf = exp ( - y ) * b0;

    *q = a0 - zabs * ( a1 - a2 * y 
      / ( y + a3 - a4 
      / ( y + a5 + a6 
      / ( y + a7 ))));
  }
//
//  |Z| between 1.28 and 12.7
//
  else if ( zabs <= 12.7 )
  {
    y = a0 * z * z;
    *pdf = exp ( - y ) * b0;

    *q = *pdf 
      / ( zabs - b1 + b2 
      / ( zabs + b3 + b4 
      / ( zabs - b5 + b6 
      / ( zabs + b7 - b8 
      / ( zabs + b9 + b10 
      / ( zabs + b11 ))))));
  }
//
//  Z far out in tail.
//
  else
  {
    *q = 0.0;
    *pdf = 0.0;
  }

  if ( z < 0.0 )
  {
    *p = *q;
    *q = 1.0 - *p;
  }
  else
  {
   *p = 1.0 - *q;
  }

  return;
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
