# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa152.hpp"

//****************************************************************************80

double alnfac ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    ALNFAC computes the logarithm of the factorial of N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial.
//
//    Output, double ALNFAC, the logarithm of the factorial of N.
//
{
  int ier;
  double value;

  value = alngam ( ( double ) ( n + 1 ), &ier );

  return value;
}
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

double chyper ( bool point, int kk, int ll, int mm, int nn, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    CHYPER computes point or cumulative hypergeometric probabilities.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by Richard Lund.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    PR Freeman,
//    Algorithm AS 59:
//    Hypergeometric Probabilities,
//    Applied Statistics,
//    Volume 22, Number 1, 1973, pages 130-133.
//
//    Richard Lund,
//    Algorithm AS 152:
//    Cumulative hypergeometric probabilities,
//    Applied Statistics,
//    Volume 29, Number 2, 1980, pages 221-223.
//
//    BL Shea,
//    Remark AS R77:
//    A Remark on Algorithm AS 152: Cumulative hypergeometric probabilities,
//    Applied Statistics,
//    Volume 38, Number 1, 1989, pages 199-204.
//
//  Parameters:
//
//    Input, bool POINT, is TRUE if the point probability is desired,
//    and FALSE if the cumulative probability is desired.
//
//    Input, int KK, the sample size.
//    0 <= KK <= MM.
//
//    Input, int LL, the number of successes in the sample.
//    0 <= LL <= KK.
//
//    Input, int MM, the population size that was sampled.
//    0 <= MM.
//
//    Input, int NN, the number of "successes" in the population.
//    0 <= NN <= MM.
//
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double CHYPER, the PDF (point probability) of
//    exactly LL successes out of KK samples, or the CDF (cumulative
//    probability) of up to LL successes out of KK samples.
//
{
  double arg;
  bool dir;
  double elimit = - 88.0;
  int i;
  int j;
  int k;
  int kl;
  int l;
  int m;
  int mbig = 600;
  double mean;
  int mnkl;
  int mvbig = 1000;
  int n;
  int nl;
  double p;
  double pt;
  double rootpi = 2.506628274631001;
  double scale = 1.0E+35;
  double sig;
  double value;

  *ifault = 0;

  k = kk + 1;
  l = ll + 1;
  m = mm + 1;
  n = nn + 1;

  dir = true;
//
//  Check arguments are within permitted limits.
//
  value = 0.0;

  if ( n < 1 || m < n || k < 1 || m < k )
  {
    *ifault = 1;
    return value;
  }

  if ( l < 1 || m - n < k - l )
  {
    *ifault = 2;
    return value;
  }

  if ( !point )
  {
    value = 1.0;
  }

  if ( n < l || k < l )
  {
    *ifault = 2;
    return value;
  }

  *ifault = 0;
  value = 1.0;

  if ( k == 1 || k == m || n == 1 || n == m )
  {
    return value;
  }

  if ( !point && ll == i4_min ( kk, nn ) )
  {
    return value;
  }

  p = ( double ) ( nn ) / ( double ) ( mm - nn );

  if ( 16.0 * r8_max ( p, 1.0 / p ) 
    < ( double ) ( i4_min ( kk, mm - kk ) ) &&
    mvbig < mm && - 100.0 < elimit )
  {
//
//  Use a normal approximation.
//
    mean = ( double ) ( kk * nn ) / ( double ) ( mm );

    sig = sqrt ( mean * ( ( double ) ( mm - nn ) / ( double ) ( mm ) ) 
    * ( ( double ) ( mm - kk ) / ( ( double ) ( mm - 1 ) ) ) );

    if ( point )
    {
      arg = - 0.5 * ( pow ( ( ( double ) ( ll ) - mean ) / sig, 2 ) );
      if ( elimit <= arg )
      {
        value = exp ( arg ) / ( sig * rootpi );
      }
      else
      {
        value = 0.0;
      }
    }
    else
    {
      value = alnorm ( ( ( double ) ( ll ) + 0.5 - mean ) / sig, false );
    }
  }
  else
  {
//
//  Calculate exact hypergeometric probabilities.
//  Interchange K and N if this saves calculations.
//
    if ( i4_min ( n - 1, m - n ) < i4_min ( k - 1, m - k ) )
    {
      i = k;
      k = n;
      n = i;
    }

    if ( m - k < k - 1 )
    {
      dir = !dir;
      l = n - l + 1;
      k = m - k + 1;
    }

    if ( mbig < mm )
    {
//
//  Take logarithms of factorials.
//
      p = alnfac ( nn ) 
        - alnfac ( mm ) 
        + alnfac ( mm - kk ) 
        + alnfac ( kk ) 
        + alnfac ( mm - nn ) 
        - alnfac ( ll ) 
        - alnfac ( nn - ll ) 
        - alnfac ( kk - ll ) 
        - alnfac ( mm - nn - kk + ll );

      if ( elimit <= p )
      {
        value = exp ( p );
      }
      else
      {
        value = 0.0;
      }
    }
    else
    {
//
//  Use Freeman/Lund algorithm.
//
      for ( i = 1; i <= l - 1; i++ )
      {
        value = value * ( double ) ( ( k - i ) * ( n - i ) ) 
        / ( double ) ( ( l - i ) * ( m - i ) );
      }

      if ( l != k )
      {
        j = m - n + l;
        for ( i = l; i <= k - 1; i++ )
        {
          value = value * ( double ) ( j - i ) / ( double ) ( m - i );
        }
      }
    }

    if ( point )
    {
      return value;
    }

    if ( value == 0.0 )
    {
//
//  We must recompute the point probability since it has underflowed.
//
      if ( mm <= mbig )
      {
        p = alnfac ( nn ) 
          - alnfac ( mm ) 
          + alnfac ( kk ) 
          + alnfac ( mm - nn ) 
          - alnfac ( ll ) 
          - alnfac ( nn - ll ) 
          - alnfac ( kk - ll ) 
          - alnfac ( mm - nn - kk + ll ) 
          + alnfac ( mm - kk );
      }

      p = p + log ( scale );

      if ( p < elimit )
      {
        *ifault = 3;
        if ( ( double ) ( nn * kk + nn + kk + 1 ) 
          / ( double ) ( mm + 2 ) < ( double ) ( ll ) )
        {
          value = 1.0;
        }
        return value;
      }
      else
      {
        p = exp ( p );
      }
    }
    else
//
//  Scale up at this point.
//
    {
      p = value * scale;
    }

    pt = 0.0;
    nl = n - l;
    kl = k - l;
    mnkl = m - n - kl + 1;

    if ( l <= kl )
    {
      for ( i = 1; i <= l - 1; i++ )
      {
        p = p * ( double ) ( ( l - i ) * ( mnkl - i ) ) / 
        ( double ) ( ( nl + i ) * ( kl + i ) );
        pt = pt + p;
      }
    }
    else
    {
      dir = !dir;
      for ( j = 0; j <= kl - 1; j++ )
      {
        p = p * ( double ) ( ( nl - j ) * ( kl - j ) ) 
        / ( double ) ( ( l + j ) * ( mnkl + j ) );
        pt = pt + p;
      }
    }

    if ( p == 0.0 )
    {
      *ifault = 3;
    }

    if ( dir )
    {
      value = value + ( pt / scale );
    }
    else
    {
      value = 1.0 - ( pt / scale );
    }
  }

  return value;
}
//****************************************************************************80

void hypergeometric_cdf_values ( int *n_data, int *sam, int *suc, int *pop, 
  int *n, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = HypergeometricDistribution [ sam, suc, pop ]
//      CDF [ dist, n ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2004
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
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *SAM, int *SUC, int *POP, the sample size, 
//    success size, and population parameters of the function.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 16

  double fx_vec[N_MAX] = { 
     0.6001858177500578E-01,  
     0.2615284665839845E+00,  
     0.6695237889132748E+00,  
     0.1000000000000000E+01,  
     0.1000000000000000E+01,  
     0.5332595856827856E+00,  
     0.1819495964117640E+00,  
     0.4448047017527730E-01,  
     0.9999991751316731E+00,  
     0.9926860896560750E+00,  
     0.8410799901444538E+00,  
     0.3459800113391901E+00,  
     0.0000000000000000E+00,  
     0.2088888139634505E-02,  
     0.3876752992448843E+00,  
     0.9135215248834896E+00 };

  int n_vec[N_MAX] = { 
     7,  8,  9, 10, 
     6,  6,  6,  6, 
     6,  6,  6,  6, 
     0,  0,  0,  0 };

  int pop_vec[N_MAX] = { 
    100, 100, 100, 100, 
    100, 100, 100, 100, 
    100, 100, 100, 100, 
    90,  200, 1000, 10000 };

  int sam_vec[N_MAX] = { 
    10, 10, 10, 10, 
     6,  7,  8,  9, 
    10, 10, 10, 10, 
    10, 10, 10, 10 };

  int suc_vec[N_MAX] = { 
    90, 90, 90, 90, 
    90, 90, 90, 90, 
    10, 30, 50, 70, 
    90, 90, 90, 90 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }
 
  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *sam = 0;
    *suc = 0;
    *pop = 0;
    *n = 0;
    *fx = 0.0;
  }
  else
  {
    *sam = sam_vec[*n_data-1];
    *suc = suc_vec[*n_data-1];
    *pop = pop_vec[*n_data-1];
    *n = n_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hypergeometric_pdf_values ( int *n_data, int *sam, int *suc, int *pop, 
  int *n, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      dist = HypergeometricDistribution [ sam, suc, pop ]
//      PDF [ dist, n ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2008
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
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *SAM, int *SUC, int *POP, the sample size, 
//    success size, and population parameters of the function.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 16

  double fx_vec[N_MAX] = { 
    0.05179370533242827E+00, 
    0.2015098848089788E+00, 
    0.4079953223292903E+00, 
    0.3304762110867252E+00, 
    0.5223047493549780E+00, 
    0.3889503452643453E+00, 
    0.1505614239732950E+00, 
    0.03927689321042477E+00, 
    0.00003099828465518108E+00, 
    0.03145116093938197E+00, 
    0.2114132170316862E+00, 
    0.2075776621999210E+00, 
    0.0000000000000000E+00, 
    0.002088888139634505E+00, 
    0.3876752992448843E+00, 
    0.9135215248834896E+00 };

  int n_vec[N_MAX] = { 
     7,  8,  9, 10, 
     6,  6,  6,  6, 
     6,  6,  6,  6, 
     0,  0,  0,  0 };

  int pop_vec[N_MAX] = { 
    100, 100, 100, 100, 
    100, 100, 100, 100, 
    100, 100, 100, 100, 
    90,  200, 1000, 10000 };

  int sam_vec[N_MAX] = { 
    10, 10, 10, 10, 
     6,  7,  8,  9, 
    10, 10, 10, 10, 
    10, 10, 10, 10 };

  int suc_vec[N_MAX] = { 
    90, 90, 90, 90, 
    90, 90, 90, 90, 
    10, 30, 50, 70, 
    90, 90, 90, 90 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }
 
  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *sam = 0;
    *suc = 0;
    *pop = 0;
    *n = 0;
    *fx = 0.0;
  }
  else
  {
    *sam = sam_vec[*n_data-1];
    *suc = suc_vec[*n_data-1];
    *pop = pop_vec[*n_data-1];
    *n = n_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
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
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( x < y )
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
