# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "asa266.hpp"

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

double digamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    DIGAMMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by Jose Bernardo.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jose Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, 1976, pages 315-317.
//
//  Parameters:
//
//    Input, double X, the argument of the digamma function.
//    0 < X.
//
//    Output, double DIGAMMA, the value of the digamma function at X.
//
{
  double euler_mascheroni = 0.57721566490153286060;
  double r;
  double value;
  double x2;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "DIGAMMA - Fatal error!\n";
    cerr << "  X <= 0.0.\n";
    exit ( 1 );
  }
//
//  Initialize.
//
  x2 = x;
  value = 0.0;
//
//  Use approximation for small argument.
//
  if ( x2 <= 0.00001 )
  {
    value = - euler_mascheroni - 1.0 / x2;
    return value;
  }
//
//  Reduce to DIGAMMA(X + N).
//
  while ( x2 < 8.5 )
  {
    value = value - 1.0 / x2;
    x2 = x2 + 1.0;
  }
//
//  Use Stirling's (actually de Moivre's) expansion.
//
  r = 1.0 / x2;
  value = value + log ( x2 ) - 0.5 * r;
  r = r * r;
  value = value 
    - r * ( 1.0 / 12.0
    - r * ( 1.0 / 120.0 
    - r *   1.0 / 252.0 ) );

  return value;
}
//****************************************************************************80

void dirichlet_check ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_CHECK checks the parameters of the Dirichlet PDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A[N], the probabilities for each component.
//    Each A(I) should be nonnegative, and at least one should be positive.
//
{
  int i;
  int positive;

  positive = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 ) 
    {
      cerr << "\n";
      cerr << "DIRICHLET_CHECK - Fatal error!\n";
      cerr << "  A(I) < 0.\n";
      cerr << "  For I = " << i << "\n";
      cerr << "  A[I] = " << a[i] << "\n";
      exit ( 1 );
    }
    else if ( 0.0 < a[i] )
    {
      positive = 1;
    }
  }

  if ( ! positive )
  {
    cerr << "\n";
    cerr << "DIRICHLET_CHECK - Fatal error!\n";
    cerr << "  All parameters are zero!\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void dirichlet_estimate ( int k, int n, double x[], int ix, int init, 
  double alpha[], double &rlogl, double v[], double g[], int &niter, 
  double &s, double &eps, int &ifault )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_ESTIMATE estimates the parameters of a Dirichlet distribution.
//
//  Discussion:
//
//    This routine requires several auxilliary routines:
//
//      ALOGAM (CACM algorithm 291 or AS 245),
//      DIGAMA (AS 103),
//      GAMMAD (AS 239),
//      PPCHI2 (AS 91),
//      TRIGAM (AS 121).
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
//    Original FORTRAN77 version by A Naryanan.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    A. Naryanan,
//    Algorithm AS 266:
//    Maximum Likelihood Estimation of the Parameters of the
//    Dirichlet Distribution,
//    Applied Statistics,
//    Volume 40, Number 2, 1991, pages 365-374.
//
//  Parameters:
//
//    Input, int K, the number of parameters.
//    2 <= K.
//
//    Input, int N, the number of observations.
//    K < N.
//
//    Input, double X[IX*K], contains the N by K array of samples
//    from the distribution.  X(I,J) is the J-th component of
//    the I-th sample.
//
//    Input, int IX, the leading dimension of the array X.
//    N <= IX.
//
//    Input, int INIT, specifies how the parameter estimates
//    are to be initialized:
//    1, use the method of moments;
//    2, initialize each ALPHA to the minimum of X;
//    otherwise, the input values of ALPHA already contain estimates.
//
//    Input/output, double ALPHA[K].
//    On input, if INIT is not 1 or 2, then ALPHA must contain
//    initial estimates for the parameters.
//    On output, with IFAULT = 0, ALPHA contains the computed
//    estimates for the parameters.
//
//    Output, double &RLOGL, the value of the log-likelihood function
//    at the solution point.
//
//    Output, double V[K*K]; the covariance between ALPHA(I) and ALPHA(J).
//
//    Output, double G[K], contains an estimate of the derivative of
//    the log-likelihood with respect to each component of ALPHA.
//
//    Output, int &NITER, contains the number of Newton-Raphson
//    iterations performed.
//
//    Output, double &S, the value of the chi-squared statistic.
//
//    Output, double &EPS, contains the probability that the 
//    chi-squared statistic is less than S.
//
//    Output, int &IFAULT, error indicator.
//    0, no error, the results were computed successfully;
//    1, K < 2;
//    2, N <= K;
//    3, IX < N;
//    4, if X(I,J) <= 0 for any I or J, or if
//       ABS ( Sum ( 1 <= J <= K ) X(I,J) - 1 ) >= GAMMA = 0.001;
//    5, if IFAULT is returned nonzero from the chi-square
//       routine PPCHI2;
//    6, if ALPHA(J) <= 0 for any J during any step of the iteration;
//    7, if MAXIT iterations were carried out but convergence
//       was not achieved.
//
{
  double alpha_min = 0.00001;
  double an;
  double beta;
  double chi2;
  double gamma = 0.0001;
  double gg;
  int i;
  int i2;
  int ifault2;
  int it_max = 100;
  int it_num;
  int j;
  int kk;
  double rk;
  double sum1;
  double sum2;
  double temp;
  double varp1;
  double *work;
  double *work2;
  double x_min;
  double x11;
  double x12;

  ifault = 0;
//
//  Check the input arguments.
//
  if ( k < 2 )
  {
    ifault = 1;
    cerr << "\n";
    cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
    cerr << "  K < 2.\n";
    exit ( 1 );
  }

  if ( n <= k )
  {
    ifault = 2;
    cerr << "\n";
    cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
    cerr << "  N <= K.\n";
    exit ( 1 );
  }

  if ( ix < n )
  {
    ifault = 3;
    cerr << "\n";
    cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
    cerr << "  IX < N.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < k; j++ )
    {
      if ( x[i+j*ix] <= 0.0 )
      {
        niter = i;
        ifault = 4;
        cerr << "\n";
        cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
        cerr << "  X(I,J) <= 0.\n";
        exit ( 1 );
      }
    }

    sum2 = 0.0;
    for ( j = 0; j < k; j++ )
    {
      sum2 = sum2 + x[i+j*ix];
    }

    if ( gamma <= fabs ( sum2 - 1.0 ) )
    {
      ifault = 4;
      niter = i;
      cerr << "\n";
      cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
      cerr << "  Row I does not sum to 1.\n";
      exit ( 1 );
    }
  }

  ifault = 0;

  an = ( double ) ( n );
  rk = ( double ) ( k );
  niter = 0;
//
//  Calculate initial estimates using the method of moments.
//
  if ( init == 1 )
  {
    sum2 = 0.0;
    for ( j = 0; j < k - 1; j++ )
    {
      sum1 = 0.0;
      for ( i = 0; i < n; i++ )
      {
        sum1 = sum1 + x[i+j*ix];
      }
      alpha[j] = sum1 / an;
      sum2 = sum2 + alpha[j];
    }

    alpha[k-1] = 1.0 - sum2;

    x12 = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x12 = x12 + pow ( x[i+0*ix], 2 );
    }

    x12 = x12 / an;
    varp1 = x12 - pow ( alpha[0], 2 );

    x11 = ( alpha[0] - x12 ) / varp1;
    for ( j = 0; j < k; j++ )
    {
      alpha[j] = x11 * alpha[j];
    }
  }
//
//  Calculate initial estimates using Ronning's suggestion.
//
  else if ( init == 2 )
  {
    x_min = x[0+0*ix];
    for ( j = 0; j < k; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        x_min = r8_min ( x_min, x[i+j*ix] );
      }
    }

    x_min = r8_max ( x_min, alpha_min );

    for ( j = 0; j < k; j++ )
    {
      alpha[j] = x_min;
    }
  }
//
//  Check whether any ALPHA's are negative or zero.
//
  for ( j = 0; j < k; j++ )
  {
    if ( alpha[j] <= 0.0 )
    {
      ifault = 6;
      cerr << "\n";
      cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
      cerr << "  For J = " << j << "\n";
      cerr << "  ALPHA(J) = " << alpha[j] << "\n";
      cerr << "  but ALPHA(J) must be positive.\n";
      exit ( 1 );
    }
  }
//
//  Calculate N * log ( G(J) ) for J = 1,2,...,K and store in WORK array.
//
  work = new double[k];

  for ( j = 0; j < k; j++ )
  {
    work[j] = 0.0;
    for ( i = 0; i < n; i++ )
    {
      work[j] = work[j] + log ( x[i+j*ix] );
    }
  }
//
//  Call Algorithm AS 91 to compute CHI2, the chi-squared value.
//
  gg = lgamma ( rk / 2.0 );

  chi2 = ppchi2 ( gamma, rk, gg, &ifault2 );

  if ( ifault2 != 0 )
  {
    ifault = 5;
    cerr << "\n";
    cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
    cerr << "  PPCHI2 returns error code.\n";
    exit ( 1 );
  }
//
//  Carry out the Newton iteration.
//
  work2 = new double[k];

  for ( it_num = 1; it_num <= it_max; it_num++ )
  {
    sum2 = r8vec_sum ( k, alpha );

    sum1 = 0.0;
    for ( j = 0; j < k; j++ )
    {
      work2[j] = trigamma ( alpha[j], &ifault2 );
      sum1 = sum1 + 1.0 / work2[j];
    }

    beta = trigamma ( sum2, &ifault2 );
    beta = an * beta / ( 1.0 - beta * sum1 );

    temp = digamma ( sum2 );

    for ( j = 0; j < k; j++ )
    {
      g[j] = an * ( temp - digamma ( alpha[j] ) ) + work[j];
    }
//
//  Calculate the lower triangle of the Variance-Covariance matrix V.
//
    sum2 = beta / an / an;
    for ( i = 0; i < k; i++ )
    {
      for ( j = 0; j < k; j++ )
      {
        v[i+j*k] = sum2 / ( work2[i] * work2[j] );
        if ( i == j )
        {
          v[i+j*k] = v[i+j*k] + 1.0 / ( an * work2[j] );
        }
      }
    }
//
//  Postmultiply the Variance-Covariance matrix V by G and store in WORK2.
//
    delete [] work2;
    work2 = r8mat_mv_new ( k, k, v, g );
//
//  Update the ALPHA's.
//
    niter = it_num;

    for ( j = 0; j < k; j++ )
    {
      alpha[j] = alpha[j] + work2[j];
      alpha[j] = r8_max ( alpha[j], alpha_min );
    }

    for ( j = 0; j < k; j++ )
    {
      if ( alpha[j] <= 0.0 )
      {
        ifault = 6;
        cerr << "\n";
        cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
        cerr << "  Newton iteration " << it_num << "\n";
        cerr << "  Computed ALPHA[J] <= 0.\n";
        cerr << "  J = " << j << "\n";
        cerr << "  ALPHA[J] = " << alpha[j] << "\n";
        exit ( 1 );
      }
    }
//
//  Test for convergence.
//
    s = r8vec_dot_product ( k, g, work2 );

    if ( s < chi2 )
    {
      eps = gammad ( s / 2.0, rk / 2.0, &ifault2 );

      sum2 = r8vec_sum ( k, alpha );

      rlogl = 0.0;
      for ( j = 0; j < k; j++ )
      {
        rlogl = rlogl + ( alpha[j] - 1.0 ) * work[j] - an * 
          lgamma ( alpha[j] );
      }

      rlogl = rlogl + an * lgamma ( sum2 );

      delete [] work;
      delete [] work2;

      return;
    }
  }

  ifault = 7;

  cerr << "\n";
  cerr << "DIRICHLET_ESTIMATE - Fatal error!\n";
  cerr << "  No convergence.\n";

  exit ( 1 );
}
//****************************************************************************80

double *dirichlet_mean ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MEAN returns the means of the Dirichlet PDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A[N], the probabilities for each component.
//    Each A[I] should be nonnegative, and at least one should be positive.
//
//    Output, double DIRICHLET_MEAN[N], the means of the PDF.
//
{
  double a_sum;
  int i;
  double *mean;

  dirichlet_check ( n, a );

  mean = new double[n];

  a_sum = r8vec_sum ( n, a );

  for ( i = 0; i < n; i++ )
  {
    mean[i] = a[i] / a_sum;
  }

  return mean;
}
//****************************************************************************80

double *dirichlet_mix_mean ( int comp_max, int comp_num, int elem_num, 
  double a[], double comp_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_MEAN returns the means of a Dirichlet mixture PDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, int COMP_MAX, the leading dimension of A, which
//    must be at least COMP_NUM.
//
//    Input, int COMP_NUM, the number of components in the
//    Dirichlet mixture density, that is, the number of distinct Dirichlet PDF's
//    that are mixed together.
//
//    Input, int ELEM_NUM, the number of elements of an
//    observation.
//
//    Input, double A[COMP_MAX*ELEM_NUM], the probabilities for
//    element ELEM_NUM in component COMP_NUM.
//    Each A(I,J) should be greater than or equal to 0.0.
//
//    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of
//    the densities.
//    These do not need to be normalized.  The weight of a given component is
//    the relative probability that that component will be used to generate
//    the sample.
//
//    Output, double DIRICHLET_MIX_MEAN[ELEM_NUM], the means for each element.
//
{
  double *a_sum;
  int comp_i;
  double comp_weight_sum;
  int elem_i;
  double *mean;
//
//  Check.
//
  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      if ( a[comp_i+elem_i*comp_max] < 0.0 )
      {
        cerr << "\n";
        cerr << "DIRICHLET_MIX_MEAN - Fatal error!\n";
        cerr << "  A(COMP,ELEM) < 0.\n";
        cerr << "  COMP = " << comp_i << "\n";
        cerr << "  ELEM = " << elem_i << "\n";
        cerr << "  A[COMP,ELEM] = " << a[comp_i+elem_i*comp_max] << "\n";
        exit ( 1 );
      }
    }
  }

  comp_weight_sum = r8vec_sum ( comp_num, comp_weight );

  a_sum = new double[comp_num];

  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    a_sum[comp_i] = 0.0;
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      a_sum[comp_i] = a_sum[comp_i] + a[comp_i+elem_i*comp_max];
    }
  }

  mean = new double[elem_num];

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    mean[elem_i] = 0.0;
    for ( comp_i = 0; comp_i < comp_num; comp_i++ )
    {
      mean[elem_i] = mean[elem_i] 
        + comp_weight[comp_i] * a[comp_i+elem_i*comp_max] / a_sum[comp_i];
    }
    mean[elem_i] = mean[elem_i] / comp_weight_sum;
  }

  delete [] a_sum;

  return mean;
}
//****************************************************************************80

double *dirichlet_mix_sample ( int comp_max, int comp_num, int elem_num, 
  double a[], double comp_weight[], int &seed, int &comp )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_SAMPLE samples a Dirichlet mixture PDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, int COMP_MAX, the leading dimension of A, which 
//    must be at least COMP_NUM.
//
//    Input, int COMP_NUM, the number of components in the
//    Dirichlet mixture density, that is, the number of distinct Dirichlet PDF's
//    that are mixed together.
//
//    Input, int ELEM_NUM, the number of elements of an
//    observation.
//
//    Input, double A[COMP_MAX*ELEM_NUM], the probabilities for
//    element ELEM_NUM in component COMP_NUM.
//    Each A(I,J) should be greater than or equal to 0.0.
//
//    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of
//    the densities.
//    These do not need to be normalized.  The weight of a given component is
//    the relative probability that that component will be used to generate
//    the sample.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, int &COMP, the index of the component of the
//    Dirichlet mixture that was chosen to generate the sample.
//
//    Output, double DIRICHLET_MIX_SAMPLE[ELEM_NUM], a sample of the PDF.
//
{
  int comp_i;
  double comp_weight_sum;
  int elem_i;
  double r;
  double sum2;
  double *x;
  double x_sum;
//
//  Check.
//
  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      if ( a[comp_i+elem_i*comp_max] < 0.0 )
      {
        cerr << "\n";
        cerr << "DIRICHLET_MIX_SAMPLE - Fatal error!\n";
        cerr << "  A(COMP,ELEM) < 0.'\n";
        cerr << "  COMP = " << comp_i << "\n";
        cerr << "  ELEM = " << elem_i << "\n";
        cerr << "  A(COMP,ELEM) = " << a[comp_i+elem_i*comp_max] << "\n";
        exit ( 1 );
      }
    }
  }
//
//  Choose a particular density MIX.
//
  comp_weight_sum = r8vec_sum ( comp_num, comp_weight );

  r = r8_uniform_ab ( 0.0, comp_weight_sum, seed );

  comp = 0;
  sum2 = 0.0;

  while ( comp < comp_num )
  {
    sum2 = sum2 + comp_weight[comp];

    if ( r <= sum2 )
    {
      break;
    }
    comp = comp + 1;
  }
//
//  Sample density COMP.
//
  x = new double[elem_num];

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    x[elem_i] = gamma_sample ( a[comp+elem_i*comp_max], 1.0, seed );
  }
//
//  Normalize the result.
//
  x_sum = r8vec_sum ( elem_num, x );

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    x[elem_i] = x[elem_i] / x_sum;
  }

  return x;
}
//****************************************************************************80

double *dirichlet_sample ( int n, double a[], int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_SAMPLE samples the Dirichlet PDF.
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
//    John Burkardt
//
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 169.
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A[N], the probabilities for each component.
//    Each A(I) should be nonnegative, and at least one should be
//    positive.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double DIRICHLET_SAMPLE[N], a sample of the PDF.  The entries 
//    of X should sum to 1.
//
{
  int i;
  double *x;
  double x_sum;

  dirichlet_check ( n, a );

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = gamma_sample ( a[i], 1.0, seed );
  }
//
//  Normalize the result.
//
  x_sum = r8vec_sum ( n, x );

  for ( i = 0; i < n; i++ )
  {
    x[i] = x[i] / x_sum;
  }
  return x;
}
//****************************************************************************80

double *dirichlet_variance ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_VARIANCE returns the variances of the Dirichlet PDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A[N], the probabilities for each component.
//    Each A[I] should be nonnegative, and at least one should be positive.
//
//    Output, double DIRICHLET_VARIANCE[N], the variances of the PDF.
//
{
  double a_sum;
  int i;
  double *variance;

  dirichlet_check ( n, a );

  a_sum = r8vec_sum ( n, a );

  variance = new double[n];

  for ( i = 0; i < n; i++ )
  {
    variance[i] = a[i] * ( a_sum - a[i] ) / ( a_sum * a_sum * ( a_sum + 1.0 ) );
  }

  return variance;
}
//****************************************************************************80

double exponential_01_sample ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_SAMPLE samples the Exponential PDF with parameters 0, 1.
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
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double EXPONENTIAL_01_SAMPLE, a sample of the PDF.
//
{
  double a;
  double b;
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  a = 0.0;
  b = 1.0;
  x = exponential_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double exponential_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF_INV inverts the Exponential CDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXPONENTIAL_CDF_INV, the corresponding argument.
//
{
  double x;
//
//  Check.
//
  if ( b <= 0.0 )
  {
    cerr << "\n";
    cerr << "EXPONENTIAL_CDF_INV - Fatal error!\n";
    cerr << "  B <= 0.0\n";
    exit ( 1 );
  }

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cerr << "\n";
    cerr << "EXPONENTIAL_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a - b * log ( 1.0 - cdf );

  return x;
}
//****************************************************************************80

double gamain ( double x, double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    GAMAIN computes the incomplete gamma ratio.
//
//  Discussion:
//
//    A series expansion is used if P > X or X <= 1.  Otherwise, a
//    continued fraction approximation is used.
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
//    Original FORTRAN77 version by G Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    G Bhattacharjee,
//    Algorithm AS 32:
//    The Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 19, Number 3, 1970, pages 285-287.
//
//  Parameters:
//
//    Input, double X, P, the parameters of the incomplete 
//    gamma ratio.  0 <= X, and 0 < P.
//
//    Output, int *IFAULT, error flag.
//    0, no errors.
//    1, P <= 0.
//    2, X < 0.
//    3, underflow.
//    4, error return from the Log Gamma routine.
//
//    Output, double GAMAIN, the value of the incomplete gamma ratio.
//
{
  double a;
  double acu = 1.0E-08;
  double an;
  double arg;
  double b;
  double dif;
  double factor;
  double g;
  double gin;
  int i;
  double oflo = 1.0E+37;
  double pn[6];
  double rn;
  double term;
  double uflo = 1.0E-37;
  double value;

  *ifault = 0;
//
//  Check the input.
//
  if ( p <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  if ( x < 0.0 )
  {
    *ifault = 2;
    value = 0.0;
    return value;
  }

  if ( x == 0.0 )
  {
    *ifault = 0;
    value = 0.0;
    return value;
  }

  g = lgamma ( p );

  arg = p * log ( x ) - x - g;

  if ( arg < log ( uflo ) )
  {
    *ifault = 3;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  factor = exp ( arg );
//
//  Calculation by series expansion.
//
  if ( x <= 1.0 || x < p )
  {
    gin = 1.0;
    term = 1.0;
    rn = p;

    for ( ; ; )
    {
      rn = rn + 1.0;
      term = term * x / rn;
      gin = gin + term;

      if ( term <= acu )
      {
        break;
      }
    }

    value = gin * factor / p;
    return value;
  }
//
//  Calculation by continued fraction.
//
  a = 1.0 - p;
  b = a + x + 1.0;
  term = 0.0;

  pn[0] = 1.0;
  pn[1] = x;
  pn[2] = x + 1.0;
  pn[3] = x * b;

  gin = pn[2] / pn[3];

  for ( ; ; )
  {
    a = a + 1.0;
    b = b + 2.0;
    term = term + 1.0;
    an = a * term;
    for ( i = 0; i <= 1; i++ )
    {
      pn[i+4] = b * pn[i+2] - an * pn[i];
    }

    if ( pn[5] != 0.0 )
    {
      rn = pn[4] / pn[5];
      dif = fabs ( gin - rn );
//
//  Absolute error tolerance satisfied?
//
      if ( dif <= acu )
      {
//
//  Relative error tolerance satisfied?
//
        if ( dif <= acu * rn )
        {
          value = 1.0 - factor * gin;
          break;
        }
      }
      gin = rn;
    }

    for ( i = 0; i < 4; i++ )
    {
      pn[i] = pn[i+2];
    }

    if ( oflo <= fabs ( pn[4] ) )
    {
      for ( i = 0; i < 4; i++ )
      {
        pn[i] = pn[i] / oflo;
      }
    }
  }

  return value;
}
//****************************************************************************80

double gamma_sample ( double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_SAMPLE samples the Gamma PDF.
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
//    Original FORTRAN77 version by Joachim Ahrens, Ulrich Dieter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Computer Methods for Sampling from Gamma, Beta, Poisson and
//    Binomial Distributions,
//    Computing,
//    Volume 12, Number 3, September 1974, pages 223-246.
//
//    Joachim Ahrens, Ulrich Dieter,
//    Generating Gamma Variates by a Modified Rejection Technique,
//    Communications of the ACM,
//    Volume 25, Number 1, January 1982, pages 47-54.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double GAMMA_SAMPLE, a sample of the PDF.
//
{
  const double a1 =  0.3333333;
  const double  a2 = -0.2500030;
  const double  a3 =  0.2000062;
  const double  a4 = -0.1662921;
  const double  a5 =  0.1423657;
  const double  a6 = -0.1367177;
  const double  a7 =  0.1233795;
  const double  e1 = 1.0;
  const double  e2 = 0.4999897;
  const double  e3 = 0.1668290;
  const double  e4 = 0.0407753;
  const double  e5 = 0.0102930;
  const double  q1 =  0.04166669;
  const double  q2 =  0.02083148;
  const double  q3 =  0.00801191;
  const double  q4 =  0.00144121;
  const double  q5 = -0.00007388;
  const double  q6 =  0.00024511;
  const double  q7 =  0.00024240;

  double bcoef;
  double c;
  double d;
  double e;
  double p;
  double q;
  double q0;
  double r;
  double s;
  double s2;
  double si;
  double t;
  double u;
  double v;
  double w;
  double x;
//
//  Allow A = 0.
//
  if ( a == 0.0 )
  {
    x = 0.0;
    return x;
  }
//
//  A < 1.
//
  if ( a < 1.0 )
  {
    for ( ; ; )
    {
      p = r8_uniform_01 ( seed );
      p = ( 1.0 + 0.3678794 * a ) * p;

      e = exponential_01_sample ( seed );

      if ( 1.0 <= p )
      {
        x = - log ( ( 1.0 + 0.3678794 * a - p ) / a );

        if ( ( 1.0 - a ) * log ( x ) <= e )
        {
          x = x / b;
          return x;
        }
      }
      else
      {
        x = exp ( log ( p ) / a );

        if ( x <= e )
        {
          x = x / b;
          return x;
        }
      }
    }
  }
//
//  1 <= A.
//
  else
  {
    s2 = a - 0.5;
    s = sqrt ( a - 0.5 );
    d = sqrt ( 32.0 ) - 12.0 * sqrt ( a - 0.5 );

    t = r8_normal_01 ( seed );
    x = pow ( sqrt ( a - 0.5 ) + 0.5 * t, 2 );

    if ( 0.0 <= t )
    {
      x = x / b;
      return x;
    }

    u = r8_uniform_01 ( seed );

    if ( d * u <= t * t * t )
    {
      x = x / b;
      return x;
    }

    r = 1.0 / a;
    q0 = ( ( ( ( ( ( 
           q7   * r 
         + q6 ) * r 
         + q5 ) * r 
         + q4 ) * r 
         + q3 ) * r 
         + q2 ) * r 
         + q1 ) * r;

    if ( a <= 3.686 )
    {
      bcoef = 0.463 + s - 0.178 * s2;
      si = 1.235;
      c = 0.195 / s - 0.079 + 0.016 * s;
    }
    else if ( a <= 13.022 )
    {
      bcoef = 1.654 + 0.0076 * s2;
      si = 1.68 / s + 0.275;
      c = 0.062 / s + 0.024;
    }
    else
    {
      bcoef = 1.77;
      si = 0.75;
      c = 0.1515 / s;
    }

    if ( 0.0 < sqrt ( a - 0.5 ) + 0.5 * t )
    {
      v = 0.5 * t / s;

      if ( 0.25 < fabs ( v ) )
      {
        q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
      }
      else
      {
        q = q0 + 0.5 * t * t * ( ( ( ( ( ( 
               a7   * v 
             + a6 ) * v 
             + a5 ) * v 
             + a4 ) * v 
             + a3 ) * v 
             + a2 ) * v 
             + a1 ) * v;
      }

      if ( log ( 1.0 - u ) <= q )
      {
        x = x / b;
        return x;
      }
    }

    for ( ; ; )
    {
      e = exponential_01_sample ( seed );

      u = r8_uniform_01 ( seed );
      if ( u < 0.5 )
      {
        t = bcoef - fabs ( si * e );
      }
      else
      {
        t = bcoef + fabs ( si * e );
      }

      if ( - 0.7187449 <= t )
      {
        v = 0.5 * t / s;

        if ( 0.25 < fabs ( v ) )
        {
          q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
        }
        else
        {
          q = q0 + 0.5 * t * t * ( ( ( ( ( ( 
                 a7   * v 
               + a6 ) * v 
               + a5 ) * v 
               + a4 ) * v 
               + a3 ) * v 
               + a2 ) * v 
               + a1 ) * v;
        }

        if ( 0.0 < q )
        {
          if ( 0.5 < q )
          {
            w = exp ( q ) - 1.0;
          }
          else
          {
            w = ( ( ( ( 
                   e5   * q 
                 + e4 ) * q 
                 + e3 ) * q 
                 + e2 ) * q 
                 + e1 ) * q;
          }

          if ( c * fabs ( u ) <= w * exp ( e - 0.5 * t * t ) )
          {
            x = pow ( s + 0.5 * t, 2 ) / b;
            return x;
          }
        }
      }
    }
  }
  return x;
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
      if ( oflo <= fabs ( pn5 ) )
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

double gammds ( double x, double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMDS computes the incomplete Gamma integral.
//
//  Discussion:
//
//    The parameters must be positive.  An infinite series is used.
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
//    Original FORTRAN77 version by Chi Leung Lau.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Chi Leung Lau,
//    Algorithm AS 147:
//    A Simple Series for the Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 29, Number 1, 1980, pages 113-114.
//
//  Parameters:
//
//    Input, double X, P, the arguments of the incomplete
//    Gamma integral.  X and P must be greater than 0.
//
//    Output, int *IFAULT, error flag.
//    0, no errors.
//    1, X <= 0 or P <= 0.
//    2, underflow during the computation.
//
//    Output, double GAMMDS, the value of the incomplete
//    Gamma integral.
//
{
  double a;
  double arg;
  double c;
  double e = 1.0E-09;
  double f;
  int ifault2;
  double uflo = 1.0E-37;
  double value;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  if ( p <= 0.0 ) 
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }
//
//  LGAMMA is the natural logarithm of the gamma function.
//
  arg = p * log ( x ) - lgamma ( p + 1.0 ) - x;

  if ( arg < log ( uflo ) )
  {
    value = 0.0;
    *ifault = 2;
    return value;
  }

  f = exp ( arg );

  if ( f == 0.0 )
  {
    value = 0.0;
    *ifault = 2;
    return value;
  }

  *ifault = 0;
//
//  Series begins.
//
  c = 1.0;
  value = 1.0;
  a = p;

  for ( ; ; )
  {
    a = a + 1.0;
    c = c * x / a;
    value = value + c;

    if ( c <= e * value )
    {
      break;
    }
  }

  value = value * f;

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
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
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
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

double lngamma ( double z, int *ier )

//****************************************************************************80
//
//  Purpose:
//
//    LNGAMMA computes Log(Gamma(X)) using a Lanczos approximation.
//
//  Discussion:
//
//    This algorithm is not part of the Applied Statistics algorithms.
//    It is slower but gives 14 or more significant decimal digits
//    accuracy, except around X = 1 and X = 2.   The Lanczos series from
//    which this algorithm is derived is interesting in that it is a
//    convergent series approximation for the gamma function, whereas
//    the familiar series due to De Moivre (and usually wrongly called
//    the Stirling approximation) is only an asymptotic approximation, as
//    is the true and preferable approximation due to Stirling.
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
//    Original FORTRAN77 version by Alan Miller.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Cornelius Lanczos,
//    A precision approximation of the gamma function,
//    SIAM Journal on Numerical Analysis, B,
//    Volume 1, 1964, pages 86-96.
//
//  Parameters:
//
//    Input, double Z, the argument of the Gamma function.
//
//    Output, int *IER, error flag.
//    0, no error occurred.
//    1, Z is less than or equal to 0.
//
//    Output, double LNGAMMA, the logarithm of the gamma function of Z.
//
{
  double a[9] = {
         0.9999999999995183, 
       676.5203681218835, 
    - 1259.139216722289, 
       771.3234287757674, 
     - 176.6150291498386, 
        12.50734324009056, 
       - 0.1385710331296526, 
         0.9934937113930748E-05, 
         0.1659470187408462E-06 };
  int j;
  double lnsqrt2pi = 0.9189385332046727;
  double tmp;
  double value;

  if ( z <= 0.0 )
  {
    *ier = 1;
    value = 0.0;
    return value;
  }

  *ier = 0;

  value = 0.0;
  tmp = z + 7.0;
  for ( j = 8; 1 <= j; j-- )
  {
    value = value + a[j] / tmp;
    tmp = tmp - 1.0;
  }

  value = value + a[0];
  value = log ( value ) + lnsqrt2pi - ( z + 6.5 ) 
    + ( z - 0.5 ) * log ( z + 6.5 );

  return value;
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
//    16 January 2008
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

double ppchi2 ( double p, double v, double g, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    PPCHI2 evaluates the percentage points of the Chi-squared PDF.
//
//  Discussion
//
//    Incorporates the suggested changes in AS R85 (vol.40(1),
//    pages 233-5, 1991) which should eliminate the need for the limited
//    range for P, though these limits have not been removed
//    from the routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Donald Best, DE Roberts.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Donald Best, DE Roberts,
//    Algorithm AS 91:
//    The Percentage Points of the Chi-Squared Distribution,
//    Applied Statistics,
//    Volume 24, Number 3, 1975, pages 385-390.
//
//  Parameters:
//
//    Input, double P,  value of the chi-squared cumulative
//    probability density function.
//    0.000002 <= P <= 0.999998.
//
//    Input, double V, the parameter of the chi-squared probability
//    density function.
//    0 < V.
//
//    Input, double G, the value of log ( Gamma ( V / 2 ) ).
//
//    Output, int *IFAULT, is nonzero if an error occurred.
//    0, no error.
//    1, P is outside the legal range.
//    2, V is not positive.
//    3, an error occurred in GAMMAD.
//    4, the result is probably as accurate as the machine will allow.
//
//    Output, double PPCHI2, the value of the chi-squared random
//    deviate with the property that the probability that a chi-squared random
//    deviate with parameter V is less than or equal to PPCHI2 is P.
//
{
  double a;
  double aa = 0.6931471806;
  double b;
  double c;
  double c1 = 0.01;
  double c2 = 0.222222;
  double c3 = 0.32;
  double c4 = 0.4;
  double c5 = 1.24;
  double c6 = 2.2;
  double c7 = 4.67;
  double c8 = 6.66;
  double c9 = 6.73;
  double c10 = 13.32;
  double c11 = 60.0;
  double c12 = 70.0;
  double c13 = 84.0;
  double c14 = 105.0;
  double c15 = 120.0;
  double c16 = 127.0;
  double c17 = 140.0;
  double c18 = 175.0;
  double c19 = 210.0;
  double c20 = 252.0;
  double c21 = 264.0;
  double c22 = 294.0;
  double c23 = 346.0;
  double c24 = 420.0;
  double c25 = 462.0;
  double c26 = 606.0;
  double c27 = 672.0;
  double c28 = 707.0;
  double c29 = 735.0;
  double c30 = 889.0;
  double c31 = 932.0;
  double c32 = 966.0;
  double c33 = 1141.0;
  double c34 = 1182.0;
  double c35 = 1278.0;
  double c36 = 1740.0;
  double c37 = 2520.0;
  double c38 = 5040.0;
  double ch;
  double e = 0.5E-06;
  int i;
  int if1;
  int maxit = 20;
  double pmax = 0.999998;
  double pmin = 0.000002;
  double p1;
  double p2;
  double q;
  double s1;
  double s2;
  double s3;
  double s4;
  double s5;
  double s6;
  double t;
  double value;
  double x;
  double xx;
//
//  Test arguments and initialize.
//
  value = - 1.0;

  if ( p < pmin || pmax < p )
  {
    *ifault = 1;
    return value;
  }

  if ( v <= 0.0 )
  {
    *ifault = 2;
    return value;
  }

  *ifault = 0;
  xx = 0.5 * v;
  c = xx - 1.0;
//
//  Starting approximation for small chi-squared
//
  if ( v < - c5 * log ( p ) )
  {
    ch = pow ( p * xx * exp ( g + xx * aa ), 1.0 / xx );

    if ( ch < e )
    {
      value = ch;
      return value;
    }
  }
//
//  Starting approximation for V less than or equal to 0.32
//
  else if ( v <= c3 )
  {
    ch = c4;
    a = log ( 1.0 - p );

    for ( ; ; )
    {
      q = ch;
      p1 = 1.0 + ch * ( c7 + ch );
      p2 = ch * (c9 + ch * ( c8 + ch ) );

      t = - 0.5 + (c7 + 2.0 * ch ) / p1 - ( c9 + ch * ( c10 + 
      3.0 * ch ) ) / p2;

      ch = ch - ( 1.0 - exp ( a + g + 0.5 * ch + c * aa ) * p2 / p1) / t;

      if ( fabs ( q / ch - 1.0 ) <= c1 )
      {
        break;
      }
    }
  }
  else
  {
//
//  Call to algorithm AS 111 - note that P has been tested above.
//  AS 241 could be used as an alternative.
//
    x = ppnd ( p, ifault );
//
//  Starting approximation using Wilson and Hilferty estimate
//
    p1 = c2 / v;
    ch = v * pow ( x * sqrt ( p1 ) + 1.0 - p1, 3 );
//
//  Starting approximation for P tending to 1.
//
    if ( c6 * v + 6.0 < ch )
    {
      ch = - 2.0 * ( log ( 1.0 - p ) - c * log ( 0.5 * ch ) + g );
    }
  }
//
//  Call to algorithm AS 239 and calculation of seven term
//  Taylor series
//
  for ( i = 1; i <= maxit; i++ )
  {
    q = ch;
    p1 = 0.5 * ch;
    p2 = p - gammad ( p1, xx, &if1 );

    if ( if1 != 0 )
    {
      *ifault = 3;
      return value;
    }

    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) );
    b = t / ch;
    a = 0.5 * t - b * c;
    s1 = ( c19 + a * ( c17 + a * ( c14 + a * ( c13 + a * ( c12 + 
    c11 * a ))))) / c24;
    s2 = ( c24 + a * ( c29 + a * ( c32 + a * ( c33 + c35 * a )))) / c37;
    s3 = ( c19 + a * ( c25 + a * ( c28 + c31 * a ))) / c37;
    s4 = ( c20 + a * ( c27 + c34 * a) + c * ( c22 + a * ( c30 + c36 * a ))) / c38;
    s5 = ( c13 + c21 * a + c * ( c18 + c26 * a )) / c37;
    s6 = ( c15 + c * ( c23 + c16 * c )) / c38;
    ch = ch + t * ( 1.0 + 0.5 * t * s1 - b * c * ( s1 - b * 
    ( s2 - b * ( s3 - b * ( s4 - b * ( s5 - b * s6 ))))));

    if ( e < fabs ( q / ch - 1.0 ) )
    {
       value = ch;
       return value;
    }
  }

 *ifault = 4;
 value = ch;

  return value;
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
  if ( fabs ( p - 0.5 ) <= split )
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

double ppnd16 ( double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    PPND16 inverts the standard normal CDF.
//
//  Discussion:
//
//    The result is accurate to about 1 part in 10**16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 March 2010
//
//  Author:
//
//    Original FORTRAN77 version by Michael Wichura.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Wichura,
//    The Percentage Points of the Normal Distribution,
//    Algorithm AS 241,
//    Applied Statistics,
//    Volume 37, Number 3, pages 477-484, 1988.
//
//  Parameters:
//
//    Input, double P, the value of the cumulative probability 
//    densitity function.  0 < P < 1.  If P is outside this range, an "infinite"
//    value is returned.
//
//    Output, integer *IFAULT, error flag.
//    0, no error.
//    1, P <= 0 or P >= 1.  PPND is returned as 0.
//
//    Output, double PPND16, the normal deviate value 
//    with the property that the probability of a standard normal deviate being 
//    less than or equal to this value is P.
//
{
  static double a[8] = { 
    3.3871328727963666080,     1.3314166789178437745e+2,
    1.9715909503065514427e+3,  1.3731693765509461125e+4,
    4.5921953931549871457e+4,  6.7265770927008700853e+4,
    3.3430575583588128105e+4,  2.5090809287301226727e+3 };
  static double b[8] = {
    1.0,                       4.2313330701600911252e+1,
    6.8718700749205790830e+2,  5.3941960214247511077e+3,
    2.1213794301586595867e+4,  3.9307895800092710610e+4,
    2.8729085735721942674e+4,  5.2264952788528545610e+3 };
  static double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770e-1,
    2.27238449892691845833e-2,  7.74545014278341407640e-4 };
  static double const1 = 0.180625;
  static double const2 = 1.6;
  static double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550e-1,
    1.48103976427480074590e-1,  1.51986665636164571966e-2,
    5.47593808499534494600e-4,  1.05075007164441684324e-9 };
  static double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230e-1,
    2.65321895265761230930e-2,  1.24266094738807843860e-3,
    2.71155556874348757815e-5,  2.01033439929228813265e-7 };
  static double f[8] = {
    1.0,                        5.99832206555887937690e-1,
    1.36929880922735805310e-1,  1.48753612908506148525e-2,
    7.86869131145613259100e-4,  1.84631831751005468180e-5,
    1.42151175831644588870e-7,  2.04426310338993978564e-15 };
  double q;
  double r;
  static double split1 = 0.425;
  static double split2 = 5.0;
  double value;

  *ifault = 0;

  if ( p <= 0.0 )
  {
    *ifault = 1;
    value = -r8_huge ( );
    return value;
  }

  if ( 1.0 <= p )
  {
    *ifault = 1;
    value = r8_huge ( );
    return value;
  }

  q = p - 0.5;

  if ( fabs ( q ) <= split1 )
  {
    r = const1 - q * q;
    value = q * r8poly_value ( 8, a, r ) / r8poly_value ( 8, b, r );
  }
  else
  {
    if ( q < 0.0 )
    {
      r = p;
    }
    else
    {
      r = 1.0 - p;
    }

    if ( r <= 0.0 )
    {
      value = -1.0;
      exit ( 1 );
    }

    r = sqrt ( -log ( r ) );

    if ( r <= split2 )
    {
      r = r - const2;
      value = r8poly_value ( 8, c, r ) / r8poly_value ( 8, d, r ); 
     }
     else
     {
       r = r - split2;
       value = r8poly_value ( 8, e, r ) / r8poly_value ( 8, f, r );
    }

    if ( q < 0.0 )
    {
      value = -value;
    }

  }

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG evaluates the logarithm of the gamma function.
//
//  Discussion:
//
//    This routine calculates the LOG(GAMMA) function for a positive real
//    argument X.  Computation is based on an algorithm outlined in
//    references 1 and 2.  The program uses rational functions that
//    theoretically approximate LOG(GAMMA) to at least 18 significant
//    decimal digits.  The approximation for X > 12 is from reference
//    3, while approximations for X < 12.0 are similar to those in
//    reference 1, but are unpublished.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2013
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody, Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the
//    Gamma Function,
//    Mathematics of Computation,
//    Volume 21, Number 98, April 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA_LOG, the value of the function.
//
{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04,
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double corr;
  static double d1 = -5.772156649015328605195174E-01;
  static double d2 = 4.227843350984671393993777E-01;
  static double d4 = 1.791759469228055000094023;
  static double frtbig = 2.25E+76;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = { 
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+10, 
    1.702665737765398868392998E+11, 
    4.926125793377430887588120E+11, 
    5.606251856223951465078242E+11 };
  double q1[8] = { 
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = { 
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = { 
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+10, 
    1.016803586272438228077304E+11, 
    3.417476345507377132798597E+11, 
    4.463158187419713286462081E+11 };
  double res;
  static double sqrtpi = 0.9189385332046727417803297;
  static double xbig = 2.55E+305;
  double xden;
  static double xinf = 1.79E+308;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double y;
  double ysq;

  y = x;

  if ( 0.0 < y && y <= xbig )
  {
    if ( y <= r8_epsilon ( ) )
    {
      res = - log ( y );
    }
//
//  EPS < X <= 1.5.
//
    else if ( y <= 1.5 )
    {
      if ( y < 0.6796875 )
      {
        corr = -log ( y );
        xm1 = y;
      }
      else
      {
        corr = 0.0;
        xm1 = ( y - 0.5 ) - 0.5;
      }

      if ( y <= 0.5 || 0.6796875 <= y )
      {
        xden = 1.0;
        xnum = 0.0;
        for ( i = 0; i < 8; i++ )
        {
          xnum = xnum * xm1 + p1[i];
          xden = xden * xm1 + q1[i];
        }
        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
      }
      else
      {
        xm2 = ( y - 0.5 ) - 0.5;
        xden = 1.0;
        xnum = 0.0;
        for ( i = 0; i < 8; i++ )
        {
          xnum = xnum * xm2 + p2[i];
          xden = xden * xm2 + q2[i];
        }
        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );
      }
    }
//
//  1.5 < X <= 4.0.
//
    else if ( y <= 4.0 )
    {
      xm2 = y - 2.0;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }
      res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
    }
//
//  4.0 < X <= 12.0.
//
    else if ( y <= 12.0 )
    {
      xm4 = y - 4.0;
      xden = -1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm4 + p4[i];
        xden = xden * xm4 + q4[i];
      }
      res = d4 + xm4 * ( xnum / xden );
    }
//
//  Evaluate for 12 <= argument.
//
    else
    {
      res = 0.0;

      if ( y <= frtbig )
      {
        res = c[6];
        ysq = y * y;
        for ( i = 0; i < 6; i++ )
        {
          res = res / ysq + c[i];
        }
      }
      res = res / y;
      corr = log ( y );
      res = res + sqrtpi - 0.5 * corr;
      res = res + y * ( corr - 1.0 );
    }
  }
//
//  Return for bad arguments.
//
  else
  {
    res = xinf;
  }
//
//  Final adjustments and return.
//
  return res;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

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

double r8_normal_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
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
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8_NORMAL_01, a normally distributed random value.
//
{
  double pi = 3.141592653589793;
  double r1;
  double r2;
  static int seed2 = 0;
  static int used = 0;
  double x;
  static double y = 0.0;
//
//  If we've used an even number of values so far, generate two more, 
//  return one, and save one.
//
  if ( ( used % 2 )== 0 )
  {
    r1 = r8_uniform_01 ( seed );
    seed2 = seed;
    r2 = r8_uniform_01 ( seed2 );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );
  }
  else
  {
    seed = seed2;
    x = y;
  }

  used = used + 1;

  return x;
}
//****************************************************************************80

double r8_psi ( double xx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI evaluates the psi function.
//
//  Discussion:
//
//    This routine evaluates the logarithmic derivative of the
//    Gamma function,
//
//      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
//             = d/dX LN ( GAMMA(X) )
//
//    for real X, where either
//
//      - XMAX1 < X < - XMIN, and X is not a negative integer,
//
//    or
//
//      XMIN < X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Anthony Strecok, Henry Thacher.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody, Anthony Strecok, Henry Thacher,
//    Chebyshev Approximations for the Psi Function,
//    Mathematics of Computation,
//    Volume 27, Number 121, January 1973, pages 123-127.
//
//  Parameters:
//
//    Input, double XX, the argument of the psi function.
//
//    Output, double R8_PSI, the value of the psi function.  
//    PSI is assigned the value 0 when the psi function is undefined.
//
{
  double aug;
  double den;
  double fourth = 0.25;
  int i;
  int n;
  int nq;
  double p1[9] = {
   4.5104681245762934160E-03, 
   5.4932855833000385356E+00, 
   3.7646693175929276856E+02, 
   7.9525490849151998065E+03, 
   7.1451595818951933210E+04, 
   3.0655976301987365674E+05, 
   6.3606997788964458797E+05, 
   5.8041312783537569993E+05, 
   1.6585695029761022321E+05 };
  double p2[7] = {
  -2.7103228277757834192E+00, 
  -1.5166271776896121383E+01, 
  -1.9784554148719218667E+01, 
  -8.8100958828312219821E+00, 
  -1.4479614616899842986E+00, 
  -7.3689600332394549911E-02, 
  -6.5135387732718171306E-21 };
  double piov4 = 0.78539816339744830962;
  double q1[8] = {
   9.6141654774222358525E+01, 
   2.6287715790581193330E+03, 
   2.9862497022250277920E+04, 
   1.6206566091533671639E+05, 
   4.3487880712768329037E+05, 
   5.4256384537269993733E+05, 
   2.4242185002017985252E+05, 
   6.4155223783576225996E-08 };
  double q2[6] = {
   4.4992760373789365846E+01, 
   2.0240955312679931159E+02, 
   2.4736979003315290057E+02, 
   1.0742543875702278326E+02, 
   1.7463965060678569906E+01, 
   8.8427520398873480342E-01 };
  double sgn;
  double three = 3.0;
  double upper;
  double value;
  double w;
  double x;
  double x01 = 187.0E+00;
  double x01d = 128.0E+00;
  double x02 = 6.9464496836234126266E-04;
  double xinf = 1.70E+38;
  double xlarge = 2.04E+15;
  double xmax1 = 3.60E+16;
  double xmin1 = 5.89E-39;
  double xsmall = 2.05E-09;
  double z;

  x = xx;
  w = fabs ( x );
  aug = 0.0;
//
//  Check for valid arguments, then branch to appropriate algorithm.
//
  if ( xmax1 <= - x || w < xmin1 )
  {
    if ( 0.0 < x )
    {
      value = - r8_huge ( );
    }
    else
    {
      value = r8_huge ( );
    }
    return value;
  }
//
//  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
//  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
//
  if ( x < 0.5 )
  {
    if ( w <= xsmall )
    {
      aug = - 1.0 / x;
    }
//
//  Argument reduction for cotangent.
//
    else
    {
      if ( x < 0.0 )
      {
        sgn = piov4;
      }
      else
      {
        sgn = - piov4;
      }

      w = w - ( double ) ( ( int ) ( w ) );
      nq = ( int ) ( w * 4.0 );
      w = 4.0 * ( w - ( double ) ( nq ) * fourth );
//
//  W is now related to the fractional part of 4.0 * X.
//  Adjust argument to correspond to values in the first
//  quadrant and determine the sign.
//
      n = nq / 2;

      if ( n + n != nq )
      {
        w = 1.0 - w;
      }

      z = piov4 * w;

      if ( ( n % 2 ) != 0 )
      {
        sgn = - sgn;
      }
//
//  Determine the final value for  -pi * cotan(pi*x).
//
      n = ( nq + 1 ) / 2;
      if ( ( n % 2 ) == 0 )
      {
//
//  Check for singularity.
//
        if ( z == 0.0 )
        {
          if ( 0.0 < x )
          {
            value = - r8_huge ( );
          }
          else
          {
            value = r8_huge ( );
          }
          return value;
        }
        aug = sgn * ( 4.0 / tan ( z ) );
      }
      else
      {
        aug = sgn * ( 4.0 * tan ( z ) );
      }
    }
    x = 1.0 - x;
  }
//
//  0.5 <= X <= 3.0.
//
  if ( x <= three )
  {
    den = x;
    upper = p1[0] * x;
    for ( i = 0; i < 7; i++ )
    {
      den = ( den + q1[i] ) * x;
      upper = ( upper + p1[i+1] ) * x;
    }
    den = ( upper + p1[8] ) / ( den + q1[7] );
    x = ( x - x01 / x01d ) - x02;
    value = den * x + aug;
    return value;
  }
//
//  3.0 < X.
//
  if ( x < xlarge )
  {
    w = 1.0 / ( x * x );
    den = w;
    upper = p2[0] * w;
    for ( i = 0; i < 5; i++ )
    {
      den = ( den + q2[i] ) * w;
      upper = ( upper + p2[i+1] ) * w;
    }
    aug = ( upper + p2[6] ) / ( den + q2[5] ) - 0.5 / x + aug;
  }

  value = aug + log ( x );

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

double *r8col_mean ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_MEAN returns the column means of an R8COL.
//
//  Discussion:
//
//    An R8COL is an M by N array of R8's, regarded as an array of N columns,
//    each of length M.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    R8COL_MEAN =
//      1.5  4.0  5.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], the array to be examined.
//
//    Output, double R8COL_MEAN[N], the means, or averages, of the columns.
//
{
  int i;
  int j;
  double *mean;

  mean = new double[n];

  for ( j = 0; j < n; j++ )
  {
    mean[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      mean[j] = mean[j] + a[i+j*m];
    }
    mean[j] = mean[j] / ( double ) ( m );
  }

  return mean;
}
//****************************************************************************80

double *r8col_variance ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_VARIANCE returns the variances of an R8COL.
//
//  Discussion:
//
//    An R8COL is an M by N array of R8's, regarded as an array of N columns,
//    each of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double A[M*N], the array whose variances are desired.
//
//    Output, double R8COL_VARIANCE[N], the variances of the rows.
//
{
  int i;
  int j;
  double mean;
  double *variance;

  variance = new double[n];

  for ( j = 0; j < n; j++ )
  {
    mean = 0.0;
    for ( i = 0; i < m; i++ )
    {
      mean = mean + a[i+j*m];
    }
    mean = mean / ( double ) ( m );

    variance[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      variance[j] = variance[j] + pow ( a[i+j*m] - mean, 2 );
    }

    if ( 1 < m )
    {
      variance[j] = variance[j] / ( double ) ( m - 1 );
    }
    else
    {
      variance[j] = 0.0;
    }
  }

  return variance;
}
//****************************************************************************80

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8MAT_MV_NEW[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double r8poly_value ( int m, double c[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VALUE evaluates a polynomial.
//
//  Discussion:
//
//    The polynomial 
//
//      p(x) = c1 + c2 * x + c3 * x^2 + ... + cm * x^(m-1)
//
//    is to be evaluated at the vector of values X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the degree.
//
//    Input, double C[M+1], the polynomial coefficients.  
//    C[0] is the constant term.
//
//    Input, double X, the evaluation points.
//
//    Output, double R8POLY_VALUE, the value of the polynomial at the 
//    evaluation points.
//
{
  int i;
  int j;
  double p;

  p = c[m];

  for ( i = m - 1; 0 <= i; i-- )
  {
    p = p * x + c[i];
  }
  return p;
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
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
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

double trigamma ( double x, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    TRIGAMMA calculates trigamma(x) = d**2 log(gamma(x)) / dx**2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by BE Schneider.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    BE Schneider,
//    Algorithm AS 121:
//    Trigamma Function,
//    Applied Statistics,
//    Volume 27, Number 1, pages 97-99, 1978.
//
//  Parameters:
//
//    Input, double X, the argument of the trigamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double TRIGAMMA, the value of the trigamma function at X.
//
{
  double a = 0.0001;
  double b = 5.0;
  double b2 =  0.1666666667;
  double b4 = -0.03333333333;
  double b6 =  0.02380952381;
  double b8 = -0.03333333333;
  double value;
  double y;
  double z;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  z = x;
//
//  Use small value approximation if X <= A.
//
  if ( x <= a )
  {
    value = 1.0 / x / x;
    return value;
  }
//
//  Increase argument to ( X + I ) >= B.
//
  value = 0.0;

  while ( z < b )
  {
    value = value + 1.0 / z / z;
    z = z + 1.0;
  }
//
//  Apply asymptotic formula if argument is B or greater.
//
  y = 1.0 / z / z;

  value = value + 0.5 *
      y + ( 1.0
    + y * ( b2
    + y * ( b4
    + y * ( b6
    + y *   b8 )))) / z;

  return value;
}
