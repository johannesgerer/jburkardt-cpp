# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "stroud.hpp"

//****************************************************************************80

double arc_sine ( double s )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_SINE computes the arc sine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ASIN routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//
//    In particular, you may get the value NaN returned.
//
//    This routine truncates arguments outside the range, avoiding the problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the argument.
//
//    Output, double ARC_SINE, an angle whose sine is S.
//
{
  double value;

  s = r8_max ( s, -1.0 );
  s = r8_min ( s, +1.0 );

  value = asin ( s );

  return value;
}
//****************************************************************************80

double ball_f1_nd ( double func ( int n, double x[] ), int n, double center[], 
  double r )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_F1_ND approximates an integral inside a ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
//
//  Discussion:
//
//    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the user supplied
//    function which evaluates F at the N-vector X.
//
//    Input, int N, the dimension of the space.
//
//    Input, double CENTER[N], the center of the ball.
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_F1_ND, the approximate integral of the function.
//
{
  int i;
  int ihi;
  int itemp;
  int j;
  int k;
  int khi;
  int ktemp;
  double quad;
  double result;
  double t;
  double temp;
  double u;
  double u2;
  double v;
  double volume;
  double w;
  double *x;
  double y;

  if ( r == 0.0 )
  {
    result = 0.0;
    return result;
  }

  x = new double[n];

  u2 = ( 1.0 - 2.0 * sqrt ( 1.0 / ( double ) ( n + 4 ) ) ) 
    / ( double ) ( n + 2 );
  u = sqrt ( u2 );
  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i] - r * u;
  }
  w = 1.0 / ( double ) ( ( n + 1 ) * i4_power ( 2, n ) );

  quad = 0.0;
  ihi = i4_power ( 2, n );

  for ( i = 0; i < ihi; i++ )
  {
    itemp = i;
    for ( j = 0; j < n; j++ )
    {
      u = ( center[j] - x[j] ) / r;

      if ( ( itemp % 2 ) == 1 )
      {
        x[j] = center[j] - r8_abs ( x[j] - center[j] );
      }
      else
      {
        x[j] = center[j] + r8_abs ( x[j] - center[j] );
      }
      itemp = itemp / 2;
    }
    quad = quad + w * func ( n, x );
  }

  temp = sqrt ( ( double ) ( n + 4 ) );

  t = sqrt ( 2.0 * ( double ) ( n + 1 ) / ( double ) ( n + 2 ) ) 
    / ( ( double ) ( n ) * temp );;

  y = ( 1.0 + 2.0 / ( ( double ) ( n ) * temp ) ) / ( double ) ( n + 2 );
  v = sqrt ( y - t );
  u = sqrt ( y + ( double ) ( n - 1 ) * t );

  khi = i4_power ( 2, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = center[j] - r * v;
    }
    x[i] = center[i] - r * u;

    for ( k = 0; k < khi; k++ )
    {
      ktemp = k;
      for ( j = 0; j < n; j++ )
      {
        if ( ( ktemp % 2 ) == 1 )
        {
          x[j] = center[j] - r8_abs ( x[j] - center[j] );
        }
        else
        {
          x[j] = center[j] + r8_abs ( x[j] - center[j] );
        }
        ktemp = ktemp / 2;
      }
      quad = quad + w * func ( n, x );
    }
    x[i] = center[i] - r * v;
  }

  volume = ball_volume_nd ( n, r );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double ball_f3_nd ( double func ( int n, double x[] ), int n, double center[], 
  double r )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_F3_ND approximates an integral inside a ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N) - CENTER(1:N) )^2 <= R * R.
//
//  Discussion:
//
//    A 2**(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the user supplied
//    function which evaluates F at the N-vector X.
//
//    Input, int N, the dimension of the space.
//
//    Input, double CENTER[N], the center of the ball.
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_F3_ND, the approximate integral of the function.
//
{
  int i;
  int j;
  int jhi;
  int jtemp;
  int k;
  double quad;
  double result;
  double ri;
  double s;
  double volume;
  double weight;
  double *x;

  if ( r == 0.0 )
  {
    result = 0.0;
    return result;
  }

  x = new double[n];

  quad = 0.0;
//
//  The first point is the center of the ball.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i];
  }
  weight = 4.0 / ( double ) ( i4_power ( n + 2, 2 ) );
  quad = quad + weight * func ( n, x );

  s = 1.0 / sqrt ( ( double ) ( n + 4 ) );

  for ( i = 0; i < n; i++ )
  {
    ri = sqrt ( ( double ) ( i + 3 ) / ( double ) ( n + 4 ) );
//
//  Set up the first point, with (I) zeroes, RI, and then N-I-1 S's.
//
    for ( j = 0; j < n; j++ )
    {
      if ( j < i )
      {
        x[j] = center[j];
      }
      else if ( j == i )
      {
        x[j] = center[j] + r * ri;
      }
      else
      {
        x[j] = center[j] + r * s;
      }
    }

    weight = pow ( 2.0, i + 1 - n ) * ( double ) ( n + 4 ) 
      / ( double ) ( ( i + 2 ) * ( i + 3 ) * ( n + 2 ) );
//
//  Now go through all sign permutations of the basic point.
//
    jhi = i4_power ( 2, n - i );

    for ( j = 0; j < jhi; j++ )
    {
      jtemp = j;

      for ( k = i; k < n; k++ )
      {
        if ( ( jtemp % 2 ) == 1 )
        {
          x[k] = center[k] - r8_abs ( x[k] - center[k] );
        }
        else
        {
          x[k] = center[k] + r8_abs ( x[k] - center[k] );
        }
        jtemp = jtemp / 2;
      }
      quad = quad + weight * func ( n, x );
    }
  }

  volume = ball_volume_nd ( n, r );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double ball_monomial_nd ( int n, int p[], double r )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_MONOMIAL_ND integrates a monomial on a ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gerald Folland,
//    How to Integrate a Polynomial Over a Sphere,
//    American Mathematical Monthly,
//    Volume 108, Number 5, May 2001, pages 446-448.
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, int P[N], the exponents of X(1) through X(N) in the monomial.
//    The exponents must be nonnegative.
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_MONOMIAL_ND, the integral of
//    X1**P(1)*X2**P(2)*...*XN**P(N) over the ball.
//
{
  int i;
  double power;
  double value;

  power = ( double ) ( n );
  for ( i = 0; i < n; i++ )
  {
    power = power + ( double ) p[i];
  }

  value = sphere_unit_monomial_nd ( n, p ) * pow ( r, power ) / power;

  return value;
}
//****************************************************************************80

double ball_unit_07_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_UNIT_07_3D approximates an integral inside the unit ball in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z <= 1.
//
//  Discussion:
//
//    A 64 point 7-th degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of 
//    the user supplied function which evaluates F(X,Y,Z).
//
//    Output, double BALL_UNIT_07_3D, the approximate integral of the function.
//
{
  int order = 4;

  double angle;
  int i;
  int j;
  int k;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double volume;
  double w;
  double weight1[4] = {
    0.19455533421780251826, 
    0.13877799911553081506, 
    0.13877799911553081506, 
    0.19455533421780251826 };
  double weight2[4];
  double weight3[4];
  double x;
  double xtab1[4] = {
    -0.906179845938663992797626878299, 
    -0.538469310105683091036314420700, 
     0.538469310105683091036314420700, 
     0.906179845938663992797626878299 };
  double xtab2[4];
  double xtab3[4];
  double y;
  double z;
//
//  Set XTAB2 and WEIGHT2.
//
  for ( j = 0; j < order; j++ )
  {
    angle = pi * ( double ) ( 2 * j - 1 ) / ( double ) ( 2 * order );
    xtab2[j] = cos ( angle );
  }

  for ( j = 0; j < order; j++ )
  {
    weight2[j] = 1.0;
  }
//
//  Set XTAB3 and WEIGHT3 for the interval [-1,1].
//
  legendre_set ( order, xtab3, weight3 );

  w = 3.0 / 16.0;

  quad = 0.0;

  for ( i = 0; i < order; i++ )
  {
    for ( j = 0; j < order; j++ )
    {
      for ( k = 0; k < order; k++ )
      {
        x = xtab1[i] * sqrt ( 1.0 - xtab2[j] * xtab2[j] ) 
                     * sqrt ( 1.0 - xtab3[k] * xtab3[k] );
        y = xtab1[i] * xtab2[j] * sqrt ( 1.0 - xtab3[k] * xtab3[k] );
        z = xtab1[i] * xtab3[k];

        quad = quad + w * weight1[i] * weight2[j] * weight3[k] 
          * func ( x, y, z );
      }
    }
  }

  volume = ball_unit_volume_3d ( );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double ball_unit_14_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_UNIT_14_3D approximates an integral inside the unit ball in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z <= 1.
//
//  Discussion:
//
//    A 288 point 14-th degree formula is used, Stroud number S3:14-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of 
//    the user supplied function which evaluates F(X,Y,Z).
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  int i;
  int j;
  int k;
  int l;
  int m;
  int n;
  double quad;
  double r[4] = {
    0.968160240, 0.836031107, 0.613371433, 0.324253423 };
  double result;
  double temp;
  double volume;
  double w1;
  double w2;
  double weight[4] = { 
    0.076181268, 0.126263673, 0.098048133, 0.032840260 };
  double x;
  double xtab[5] = { 
    -0.151108275, 0.315838353, 0.346307112, 
    -0.101808787, -0.409228403 };
  double y;
  double ytab[5] = { 
    0.155240600, 0.257049387, 0.666277790, 
    0.817386065, 0.501547712 };
  double z;
  double ztab[5] = { 
    0.976251323, 0.913330032, 0.660412970, 
    0.567022920, 0.762221757 };

  quad = 0.0;

  for ( m = 0; m < 4; m++ )
  {
    w1 = 125.0 * weight[m] / 3360.0;
    x = 0.525731112 * r[m];
    y = 0.850650808 * r[m];
    z = 0.0;


    for ( j = 0; j < 2; j++ )
    {
      x = -x;
      for ( k = 0; k < 2; k++ )
      {
        y = -y;
        for ( l = 0; l < 3; l++ )
        {
          temp = z;
          z = y;
          y = x;
          x = temp;
          quad = quad + w1 * func ( x, y, z );
        }
      }
    }

    w2 = 143.0 * weight[m] / 3360.0;

    for ( n = 0; n < 5; n++ )
    {
      x = xtab[n] * r[m];
      y = ytab[n] * r[m];
      z = ztab[n] * r[m];

      for ( i = 0; i < 3; i++ )
      {
        temp = x;
        x = z;
        z = -y;
        y = -temp;

        for ( j = 0; j < 3; j++ )
        {
          temp = z;
          z = y;
          y = x;
          x = temp;

          quad = quad + w2 * func ( x, y, z );
        }
        y = -y;
        z = -z;
        quad = quad + w2 * func ( x, y, z );
      }
    }
  }

  volume = ball_unit_volume_3d ( );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double ball_unit_15_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_UNIT_15_3D approximates an integral inside the unit ball in 3D.
//
//  Integration region:
//
//    X * X + Y * Y + Z * Z <= 1.
//
//  Discussion:
//
//    A 512 point 15-th degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 October 2000
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which evaluates F(X,Y,Z).
//
//    Output, double BALL_UNIT_15_3D, the approximate integral of the function.
//
{
  double cj;
  double ck;
  int i;
  int j;
  int k;
  int order1 = 4;
  int order2 = 8;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double sj;
  double sk;
  double volume;
  double w;
  double weight1[4] = {
    0.0328402599, 0.0980481327, 0.1262636728, 0.0761812678 };
  double *weight2;
  double x;
  double xtab1[4] = {
    0.3242534234, 0.6133714327, 0.8360311073, 0.9681602395 };
  double *xtab2;
  double y;
  double z;

  xtab2 = new double[order2];
  weight2 = new double[order2];

  legendre_set ( order2, xtab2, weight2 );

  w = 3.0 / 32.0;

  quad = 0.0;

  for ( i = 0; i < order1; i++ )
  {
    for ( j = 0; j < order2; j++ )
    {
      sj = xtab2[j];
      cj = sqrt ( 1.0 - sj * sj );

      for ( k = 1; k <= 16; k++ )
      {
        sk = sin ( ( double ) ( k ) * pi / 8.0 );
        ck = cos ( ( double ) ( k ) * pi / 8.0 );
        x = xtab1[i] * cj * ck;
        y = xtab1[i] * cj * sk;
        z = xtab1[i] * sj;
        quad = quad + w * weight1[i] * weight2[j] * func ( x, y, z );
      }
    }
  }

  volume = ball_unit_volume_3d ( );
  result = quad * volume;

  delete [] xtab2;
  delete [] weight2;

  return result;
}
//****************************************************************************80

double ball_unit_f1_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_UNIT_F1_ND approximates an integral inside the unit ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= 1.
//
//  Discussion:
//
//    An (N+1)*2^N point 5-th degree formula is used, Stroud number SN:5-6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function which evaluates F at the N-vector X.
//
//    Input, int N, the dimension of the space.
//
//    Output, double BALL_UNIT_F1_ND, the approximate integral of the function.
//
{
  int i;
  int ihi;
  int itemp;
  int j;
  int k;
  int khi;
  int ktemp;
  double quad;
  double result;
  double t;
  double temp;
  double u;
  double u2;
  double v;
  double volume;
  double w;
  double *x;
  double y;

  x = new double[n];

  u2 = ( 1.0 - 2.0 * sqrt ( 1.0 / ( double ) ( n + 4 ) ) ) 
    / ( double ) ( n + 2 );
  u = sqrt ( u2 );
  for ( i = 0; i < n; i++ )
  {
    x[i] = - u;
  }
  w = 1.0 / ( double ) ( ( n + 1 ) * i4_power ( 2, n ) );

  quad = 0.0;
  ihi = i4_power ( 2, n );

  for ( i = 0; i < ihi; i++ )
  {
    itemp = i;

    for ( j = 0; j < n; j++ )
    {
      if ( ( itemp % 2 ) == 1 )
      {
        x[j] = - r8_abs ( x[j] );
      }
      else
      {
        x[j] = r8_abs ( x[j] );
      }
      itemp = itemp / 2;
    }
    quad = quad + w * func ( n, x );
  }

  temp = sqrt ( ( double ) ( n + 4 ) );

  t = sqrt ( 2.0 * ( double ) ( n + 1 ) / ( double ) ( n + 2 ) ) 
    / ( ( double ) ( n ) * temp );

  y = ( 1.0 + 2.0 / ( ( double ) ( n ) * temp ) ) 
    / ( double ) ( n + 2 );
  v = sqrt ( y - t );
  u = sqrt ( y + ( double ) ( n - 1 ) * t );

  khi = i4_power ( 2, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = -v;
    }
    x[i] = - u;

    for ( k = 0; k < khi; k++ )
    {
      ktemp = k;

      for ( j = 0; j < n; j++ )
      {
        if ( ( ktemp % 2 ) == 1 )
        {
          x[j] = - r8_abs ( x[j] );
        }
        else
        {
          x[j] = r8_abs ( x[j] );
        }
        ktemp = ktemp / 2;
      }
      quad = quad + w * func ( n, x );
    }
    x[i] = - v;
  }

  volume = ball_unit_volume_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double ball_unit_f3_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_UNIT_F3_ND approximates an integral inside the unit ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= 1.
//
//  Discussion:
//
//    A 2^(N+1)-1 point 5-th degree formula is used, Stroud number SN:5-4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function which evaluates F at the N-vector X.
//
//    Input, int N, the dimension of the space.
//
//    Output, double BALL_UNIT_F3_ND, the approximate integral of the function.
//
{
  int i;
  int j;
  int jtemp;
  int k;
  double quad;
  double result;
  double ri;
  double s;
  double volume;
  double weight;
  double *x;

  quad = 0.0;
//
//  The first point is the center of the ball.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }
  weight = 4.0 / ( double ) ( ( n + 2 ) * ( n + 2 ) );
  quad = quad + weight * func ( n, x );

  s = 1.0 / sqrt ( ( double ) ( n + 4 ) );

  for ( i = 0; i < n; i++ )
  {
    ri = sqrt ( ( double ) ( i + 3 ) / ( double ) ( n + 4 ) );
//
//  Set up the first point, with (I-1) zeroes, RI, and then N-I S's.
//
    for ( j = 0; j < n; j++ )
    {
      if ( j < i )
      {
        x[j] = 0.0;
      }
      else if ( j == i )
      {
        x[j] = ri;
      }
      else
      {
        x[j] = s;
      }
    }

    weight = pow ( 2.0, i + 1 - n ) * ( double ) ( n + 4 ) 
      / ( double ) ( ( i + 2 ) * ( i + 3 ) * ( n + 2 ) );
//
//  Now go through all sign permutations of the basic point.
//
    for ( j = 0; j < i4_power ( 2, n - i ); j++ )
    {
      jtemp = j;

      for ( k = i; k < n; k++ )
      {
        if ( ( jtemp % 2 ) == 1 )
        {
          x[k] = - r8_abs ( x[k] );
        }
        else
        {
          x[k] = r8_abs ( x[k] );
        }
        jtemp = jtemp / 2;
      }
      quad = quad + weight * func ( n, x );
    }
  }

  volume = ball_unit_volume_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double ball_unit_volume_3d ( )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_UNIT_VOLUME_3D computes the volume of the unit ball in 3D.
//
//  Integration region:
//
//    X * X + Y * Y + Z * Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double BALL_UNIT_VOLUME_3D, the volume of the ball.
//
{
  double pi = 3.141592653589793;
  double value;

  value = ( 4.0 / 3.0 ) * pi;

  return value;
}
//****************************************************************************80

double ball_unit_volume_nd ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_UNIT_VOLUME_ND computes the volume of the unit ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= 1.
//
//  Discussion:
//
//    N  Volume
//
//    2             PI
//    3  (4/3)    * PI
//    4  (1/2)    * PI^2
//    5  (8/15)   * PI^2
//    6  (1/6)    * PI^3
//    7  (16/105) * PI^3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Output, double BALL_UNIT_VOLUME_ND, the volume of the ball.
//
{
  int i;
  int m;
  double pi = 3.141592653589793;
  double volume;

  if ( ( n % 2 ) == 0 )
  {
    m = n / 2;
    volume = pow ( pi, m );
    for ( i = 1; i <= m; i++ )
    {
      volume = volume / ( double ) ( i );
    }
  }
  else
  {
    m = ( n - 1 ) / 2;
    volume = pow ( pi, m ) * i4_power ( 2, n );
    for ( i = m + 1; i <= 2 * m + 1; i++ )
    {
      volume = volume / ( double ) ( i );
    }
  }

  return volume;
}
//****************************************************************************80

double ball_volume_3d ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_VOLUME_3D computes the volume of a ball in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z <= R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_VOLUME_3D, the volume of the ball.
//
{
  double pi = 3.141592653589793;
  double volume;

  volume = ( 4.0 / 3.0 ) * pi * r * r * r;

  return volume;
}
//****************************************************************************80

double ball_volume_nd ( int n, double r )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_VOLUME_ND computes the volume of a ball in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) <= R * R
//
//  Discussion:
//
//    N  Volume
//
//    2             PI   * R^2
//    3  (4/3)    * PI   * R^3
//    4  (1/2)    * PI^2 * R^4
//    5  (8/15)   * PI^2 * R^5
//    6  (1/6)    * PI^3 * R^6
//    7  (16/105) * PI^3 * R^7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double R, the radius of the ball.
//
//    Output, double BALL_VOLUME_ND, the volume of the ball.
//
{
  double volume;

  volume = ball_unit_volume_nd ( n ) * pow ( r, n );

  return volume;
}
//****************************************************************************80

double c1_geg_monomial_integral ( double alpha, int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on C1.
//
//  Discussion:
//
//    C1_GEG is the interval [-1,+1] with the Gegenbauer weight function
//
//      w(alpha;x) = (1-x^2)^alpha
//
//    with -1.0 < alpha.
//
//    value = integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the exponent of (1-X^2).
//    - 1.0 < ALPHA.
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double C1_GEG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double value;
  double value1;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "C1_GEG_MONOMIAL_INTEGRAL - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  c = ( double ) ( expon );

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = 2.0 * r8_gamma ( 1.0 + c ) * r8_gamma ( 1.0 + alpha ) 
    * value1 / r8_gamma ( 2.0 + alpha + c );

  return value;
}
//****************************************************************************80

double c1_jac_monomial_integral ( double alpha, double beta, int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over C1.
//
//  Discussion:
//
//    value = integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the exponent of (1-X) in the weight factor.
//
//    Input, double BETA, the exponent of (1+X) in the weight factor.
//
//    Input, int EXPON, the exponent.
//
//    Output, double C1_JAC_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;
  double value2;

  c = ( double ) ( expon );

  if ( ( expon % 2 ) == 0 )
  {
    s = +1.0;
  }
  else
  {
    s = -1.0;
  }

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + beta + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  arg1 = - beta;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value2 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = r8_gamma ( 1.0 + c ) * ( 
      s * r8_gamma ( 1.0 + beta  ) * value1 
    / r8_gamma ( 2.0 + beta  + c ) 
    +     r8_gamma ( 1.0 + alpha ) * value2 
    / r8_gamma ( 2.0 + alpha + c ) );

  return value;
}
//****************************************************************************80

double c1_leg_monomial_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    C1_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on C1.
//
//  Discussion:
//
//    C1_LEG is the interval [-1,+1] with the Legendre weight function
//
//      w(x) = 1.
//
//    value = integral ( -1 <= x <= +1 ) x^expon dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double C1_LEG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double value;

  if ( expon < 0 )
  {
    cerr << "\n";
    cerr << "C1_LEG_MONOMIAL_INTEGRAL - Fatal error!\n";
    cerr << "  EXPON < 0.\n";
    exit ( 1 );
  }

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  value = 2.0 / ( double) ( expon + 1 );

  return value;
}
//****************************************************************************80

double circle_annulus ( double func ( double x, double y ), double center[2], 
  double radius1, double radius2, int nr )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_ANNULUS approximates an integral in an annulus.
//
//  Integration region:
//
//    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Peirce,
//    Numerical Integration Over the Planar Annulus,
//    Journal of the Society for Industrial and Applied Mathematics,
//    Volume 5, Number 2, June 1957, pages 66-73.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the user supplied 
//    function of two variables which is to be integrated.
//
//    Input, double CENTER[2], the center of the circle.
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Input, int NR, the order of the rule.  This quantity specifies
//    the number of distinct radii to use.  The number of angles used will
//    be 4*NR, for a total of 4*NR*NR points.
//
//    Output, double CIRCLE_ANNULUS, the approximation to the integral.
//
{
  double a;
  double area;
  double b;
  double c;
  double d;
  int i;
  int j;
  int nt;
  double pi = 3.141592653589793;
  double quad;
  double *ra;
  double result;
  double *rw;
  double t;
  double tw;
  double x;
  double y;
//
//  Choose radial abscissas and weights.
//
  ra = new double[nr];
  rw = new double[nr];

  legendre_set ( nr, ra, rw );
  a = -1.0;
  b = +1.0;
  c = radius1 * radius1;
  d = radius2 * radius2;

  rule_adjust ( a, b, c, d, nr, ra, rw );

  for ( i = 0; i < nr; i++ )
  {
    ra[i] = sqrt ( ra[i] );
  }
  for ( i = 0; i < nr; i++ )
  {
    rw[i] = rw[i] / ( radius2 - radius1 ) / ( radius2 + radius1 );
  }
//
//  Set angular abscissas and weights.
//
  nt = 4 * nr;

  tw = 1.0 / ( double ) ( nt );
//
//  Approximate the integral.
//
  quad = 0.0;
  for ( i = 0; i < nt; i++ )
  {
    t = 2.0 * pi * ( double ) ( i - 1 ) / ( double ) ( nt );
    for ( j = 0; j < nr; j++ )
    {
      x = center[0] + ra[j] * cos ( t );
      y = center[1] + ra[j] * sin ( t );
      quad = quad + tw * rw[j] * func ( x, y );
    }
  }

  area = circle_annulus_area_2d ( radius1, radius2 );
  result = quad * area;

  delete [] ra;
  delete [] rw;

  return result;
}
//****************************************************************************80

double circle_annulus_area_2d ( double radius1, double radius2 )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_ANNULUS_AREA_2D returns the area of a circular annulus in 2D.
//
//  Integration region:
//
//    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Output, double CIRCLE_ANNULUS_AREA_2D, the area of the annulus.
//
{
  double pi = 3.141592653589793;
  double value;

  value = pi * ( radius1 + radius2 ) * ( radius2 - radius1 );

  return value;
}
//****************************************************************************80

double circle_annulus_sector ( double func ( double x, double y ), 
  double center[2], double radius1, double radius2, double theta1, 
  double theta2, int nr )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_ANNULUS_SECTOR approximates an integral in a circular annulus sector.
//
//  Discussion:
//
//    A circular annulus sector comprises the area between two concentric
//    circles and two concentric rays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Peirce,
//    Numerical Integration Over the Planar Annulus,
//    Journal of the Society for Industrial and Applied Mathematics,
//    Volume 5, Number 2, June 1957, pages 66-73.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the center of the circle.
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Input, double THETA1, THETA2, the angles defining the sector.
//    The sector is measured from THETA1 to THETA2.
//
//    Input, int NR, the order of the rule.  This quantity specifies
//    the number of distinct radii to use.  The number of angles used will
//    be 4*NR, for a total of 4*NR*NR points.
//
//    Output, double CIRCLE_ANNULUS_SECTOR, the approximation to the integral.
//
{
  double a;
  double area;
  double b;
  double c;
  double d;
  int i;
  int j;
  int nt;
  double quad;
  double *ra;
  double result;
  double *rw;
  double *ta;
  double *tw;
  double x;
  double y;
//
//  Set the radial abscissas and weights.
//
  ra = new double[nr];
  rw = new double[nr];

  legendre_set ( nr, ra, rw );

  a = -1.0;
  b = +1.0;
  c = radius1 * radius1;
  d = radius2 * radius2;

  rule_adjust ( a, b, c, d, nr, ra, rw );

  for ( i = 0; i < nr; i++ )
  {
    ra[i] = sqrt ( ra[i] );
  }
  for ( i = 0; i < nr; i++ )
  {
    rw[i] = rw[i] / ( radius2 - radius1 ) / ( radius2 + radius1 );
  }
//
//  Pick angles evenly spaced between THETA1 and THETA2, but do not
//  include the endpoints, and use a half interval for the first and last.
//
  nt = 4 * nr;

  ta = tvec_even_bracket3 ( nt, theta1, theta2 );
  tw = new double[nt];
  for ( i = 0; i < nt; i++ )
  {
    tw[i] = 1.0 / ( double ) ( nt );
  }
//
//  Approximate the integral.
//
  quad = 0.0;
  for ( i = 0; i < nt; i++ )
  {
    for ( j = 0; j < nr; j++ )
    {
      x = center[0] + ra[j] * cos ( ta[i] );
      y = center[1] + ra[j] * sin ( ta[i] );
      quad = quad + tw[i] * rw[j] * func ( x, y );
    }
  }

  area = circle_annulus_sector_area_2d ( radius1, radius2, theta1, theta2 );

  result = quad * area;

  delete [] ra;
  delete [] rw;
  delete [] ta;
  delete [] tw;

  return result;
}
//****************************************************************************80

double circle_annulus_sector_area_2d ( double radius1, double radius2, 
  double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_ANNULUS_SECTOR_AREA_2D returns the area of a circular annulus sector in 2D.
//
//  Discussion:
//
//    A circular annulus sector comprises the area between two concentric
//    circles and two concentric rays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RADIUS1, RADIUS2, the radii of the circles.
//
//    Input, double THETA1, THETA2, the angles of the rays.
//    Ordinarily, (THETA2-THETA1) is between 0 and 2*PI.
//
//    Output, double CIRCLE_ANNULUS_SECTOR_AREA_2D, the area of the
//    circular annulus sector.
//
{
  double area;

  area = 0.5 * ( radius1 + radius2 ) * ( radius2 - radius1 ) 
    * ( theta2 - theta1 );

  return area;
}
//****************************************************************************80

double circle_area_2d ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_AREA_2D returns the area of a circle in 2D.
//
//  Integration region:
//
//    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Output, double CIRCLE_AREA_2D, the area of the circle.
//
{
  double area;
  double pi = 3.141592653589793;

  area = pi * r * r;

  return area;
}
//****************************************************************************80

double circle_cap_area_2d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_CAP_AREA_2D computes the area of a circle cap in 2D.
//
//  Discussion:
//
//    Draw any radius R of the circle and denote as P the point where the
//    radius intersects the circle.  Now consider the point Q which lies
//    on the radius and which is H units from P.  The line which is
//    perpendicular to the radius R and passes through Q divides the
//    circle into two pieces.  The piece including the point P is the
//    circular cap of height (or thickness) H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Input, double H, the "height" of the circle cap.  
//
//    Output, double CIRCLE_CAP_AREA_2D, the area of the circle cap.
//
{
  double area;
  double pi = 3.141592653589793;
  double theta;

  if ( h <= 0.0 )
  {
    area = 0.0;
  }
  else if ( h <= r )
  {
    theta = 2.0 * arc_sine ( sqrt ( h * ( 2.0 * r - h ) ) / r );
    area = r * r * ( theta - sin ( theta ) ) / 2.0;
  }
  else if ( h <= 2.0 * r )
  {
    theta = 2.0 * arc_sine ( sqrt ( h * ( 2.0 * r - h ) ) / r );
    area = r * r * ( pi - ( theta - sin ( theta ) ) / 2.0 );
  }
  else if ( 2.0 * r <= h )
  {
    area = pi * r * r;
  }

  return area;
}
//****************************************************************************80

double circle_cum ( double func ( double x, double y ), double center[2], 
  double radius, int order )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_CUM approximates an integral on the circumference of a circle in 2D.
//
//  Integration region:
//
//    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
//
//  Discussion:
//
//    An ORDER point, (ORDER-1)-th degree formula is used, 
//    Stroud number U2:M-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the coordinates of the center of 
//    the circle.
//
//    Input, double RADIUS, the radius of the circle.
//
//    Input, int ORDER, the number of points to use.
//
//    Output, double CIRCLE_CUM, the approximate integral of the function.
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double volume;
  double x;
  double y;

  quad = 0.0;

  for ( i = 0; i < order; i++ )
  {
    angle = ( double ) ( 2 * i ) * pi / ( double ) ( order );
    x = center[0] + radius * cos ( angle );
    y = center[1] + radius * sin ( angle );
    quad = quad + func ( x, y );
  }

  quad = quad / ( double ) ( order );

  volume = pi * radius * radius;
  result = quad * volume;

  return result;
}
//****************************************************************************80

void circle_rt_set ( int rule, int nr, int nt, int nc, double ra[], 
  double rw[], double ta[], double tw[], double *cw )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_RT_SET sets an R, THETA product quadrature rule in the unit circle.
//
//  Discussion:
//
//    For a given value of RULE, here are the number of points used at the
//    center (NC), the number of points along the radial direction (NR) and
//    the number of points along the theta direction (NT).  The total number
//    of points in the rule will be 
//
//      Total = NC + NR * NT.
//
//    The user, when choosing RULE, must allocate enough space in the arrays
//    RA, RW, TA and TW for the resulting values of NR and NT.
//
//    RULE  NC  NR  NT  Total
//    ----  --  --  --  -----
//       1   1   0   0      1
//       2   0   1   4      4
//       3   1   1   4      5
//       4   1   1   6      7
//       5   1   2   4      9
//       6   0   3   4     12
//       7   1   2  10     21
//       8   0   4  16     64
//       9   0   5  20    120
//
//    The integral of F(X,Y) over the unit circle is approximated by
//
//      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
//      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
//      = approximately
//        CW * F(0,0) 
//        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
//        RW(I) * TW(J) * F ( R(I) * cos ( TA(J) ), R(I) * sin ( TA(J) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
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
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule desired.
//
//    Input, int NR, the number of R abscissas.
//
//    Input, int NT, the number of Theta abscissas.
//
//    Input, int NC, the number of center abscissas (0 or 1 ).
//
//    Output, double RA[NR], RW[NR], the R abscissas and weights.
//
//    Output, double TA[NT], TW[NT], the THETA abscissas and weights.
//
//    Output, double *ZW, the weight to use for the center.
//
{
  double a;
  double b;
  double c;
  double d;
  int i;
  double pi = 3.141592653589793;
  double u;
  double v;
  double w;

  if ( rule == 1 )
  {
    *cw = 1.0;
  }
  else if ( rule == 2 )
  {
    ra[0] = 0.5;
    rw[0] = 1.0;

    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 1.0 / ( double ) ( nt );
    }
    *cw = 0.0;
  }
  else if ( rule == 3 )
  {
    ra[0] = 1.0;
    rw[0] = 1.0;

    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 0.125;
    }
    *cw = 0.5;
  }
  else if ( rule == 4 )
  {
    ra[0] = sqrt ( 2.0 / 3.0 );
    rw[0] = 1.0;

    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 0.125;
    }
    *cw = 0.25;
  }
  else if ( rule == 5 )
  {
    a = 1.0;
    b = sqrt ( 2.0 ) / 2.0;
    u = 1.0 / 6.0;
    v = 4.0 / 6.0;

    ra[0] = a;
    ra[1] = b;
    rw[0] = u;
    rw[1] = v;

    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 1.0 / ( double ) ( nt );
    }
    *cw = 4.0 / 24.0;
  }
  else if ( rule == 6 )
  {
    a = sqrt ( 3.0 ) / 2.0;
    b = sqrt ( ( 27.0 - 3.0 * sqrt ( 29.0 ) ) / 52.0 );
    c = sqrt ( ( 27.0 + 3.0 * sqrt ( 29.0 ) ) / 52.0 );

    u = 8.0 / 27.0;
    v = ( 551.0 + 41.0 * sqrt ( 29.0 ) ) / 1566.0;
    w = ( 551.0 - 41.0 * sqrt ( 29.0 ) ) / 1566.0;

    ra[0] = a;
    ra[1] = b;
    ra[2] = c;
    rw[0] = u;
    rw[1] = v;
    rw[2] = w;

    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 1.0 / ( double ) ( nt );
    }
    *cw = 0.0;
  }
  else if ( rule == 7 )
  {
    a = sqrt ( ( 6.0 - sqrt ( 6.0 ) ) / 10.0 );
    b = sqrt ( ( 6.0 + sqrt ( 6.0 ) ) / 10.0 );
    u = ( 16.0 + sqrt ( 6.0 ) ) / 36.0;
    v = ( 16.0 - sqrt ( 6.0 ) ) / 36.0;

    ra[0] = a;
    ra[1] = b;
    rw[0] = u;
    rw[1] = v;

    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 1.0 / ( double ) ( nt );
    }
    *cw = 1.0 / 9.0;
  }
  else if ( rule == 8 )
  {
    legendre_set ( nr, ra, rw );
    a = -1.0;
    b = +1.0;
    c =  0.0;
    d = +1.0;
    rule_adjust ( a, b, c, d, nr, ra, rw );

    for ( i = 0; i < nr; i++ )
    {
      ra[i] = sqrt ( ra[i] );
    }
    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 1.0 / ( double ) ( nt );
    }
    *cw = 0.0;
  }
  else if ( rule == 9 )
  {
    legendre_set ( nr, ra, rw );
    a = -1.0;
    b = +1.0;
    c =  0.0;
    d = +1.0;
    rule_adjust ( a, b, c, d, nr, ra, rw );

    for ( i = 0; i < nr; i++ )
    {
      ra[i] = sqrt ( ra[i] );
    }
    for ( i = 0; i < nt; i++ )
    {
      ta[i] = ( double ) ( 2 * i ) * pi / ( double ) ( nt );
    }
    for ( i = 0; i < nt; i++ )
    {
      tw[i] = 1.0 / ( double ) ( nt );
    }
    *cw = 0.0;
  }
  else
  {
    cerr << "\n";
    cerr << "CIRCLE_RT_SET - Fatal error!\n";
    cerr << "  There is no rule of index " << rule << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void circle_rt_size ( int rule, int *nr, int *nt, int *nc )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_RT_SIZE sizes an R, THETA product quadrature rule in the unit circle.
//
//  Discussion:
//
//    For a given value of RULE, here are the number of points used at the
//    center (NC), the number of points along the radial direction (NR) and
//    the number of points along the theta direction (NT).  The total number
//    of points in the rule will be 
//
//      Total = NC + NR * NT.
//
//    The user, when choosing RULE, must allocate enough space in the arrays
//    RA, RW, TA and TW for the resulting values of NR and NT.
//
//    RULE  NC  NR  NT  Total
//    ----  --  --  --  -----
//       1   1   0   0      1
//       2   0   1   4      4
//       3   1   1   4      5
//       4   1   1   6      7
//       5   1   2   4      9
//       6   0   3   4     12
//       7   1   2  10     21
//       8   0   4  16     64
//       9   0   5  20    120
//
//    The integral of F(X,Y) over the unit circle is approximated by
//
//      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
//      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
//      = approximately
//        ZW * F(0,0) 
//        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
//        RW(I) * TW(J) * F ( R(I) * cos ( TA(J) ), R(I) * sin ( TA(J) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2008
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
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule desired.
//
//    Output, int *NR, the number of R abscissas.
//    
//    Output, int *NT, the number of Theta abscissas.
//
//    Output, int *NC, the number of center abscissas (0 or 1).
//
{
  if ( rule == 1 )
  {
    *nr = 0;
    *nt = 0;
    *nc = 1;
  }
  else if ( rule == 2 )
  {
    *nr = 1;
    *nt = 4;
    *nc = 0;
  }
  else if ( rule == 3 )
  {
    *nr = 1;
    *nt = 4;
    *nc = 1;
  }
  else if ( rule == 4 )
  {
    *nr = 1;
    *nt = 6;
    *nc = 1;
  }
  else if ( rule == 5 )
  {
    *nr = 2;
    *nt = 4;
    *nc = 1;
  }
  else if ( rule == 6 )
  {
    *nr = 3;
    *nt = 4;
    *nc = 0;
  }
  else if ( rule == 7 )
  {
    *nr = 2;
    *nt = 10;
    *nc = 1;
  }
  else if ( rule == 8 )
  {
    *nr = 4;
    *nt = 16;
    *nc = 0;
  }
  else if ( rule == 9 )
  {
    *nr = 5;
    *nt = 20;
    *nc = 0;
  }
  else
  {
    cerr << "\n";
    cerr << "CIRCLE_RT_SIZE - Fatal error!\n";
    cerr << "  There is no rule of index " << rule << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double circle_rt_sum ( double func ( double x, double y ), double center[2], 
  double radius, int nr, double ra[], double rw[], int nt, double ta[], 
  double tw[], double zw )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_RT_SUM applies an R, THETA product quadrature rule inside a circle.
//
//  Integration region:
//
//    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= RADIUS^2.
//
//  Discussion:
//
//    The product rule is assumed to be have the form:
//
//      Integral_Approx = ZW * F(CENTER(1),CENTER(2)) +
//        sum ( 1 <= IR <= NR ) Sum ( 1 <= IT <= NT )
//        RW(IR) * TW(IT) * F ( CENTER(1) + R(IR) * RADIUS * Cos ( TA(IT) ),
//                              CENTER(2) + R(IR) * RADIUS * Sin ( TA(IT) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the center of the circle.
//
//    Input, double RADIUS, the radius of the circle.
//
//    Input, int NR, the number of R abscissas.
//
//    Input, double RA[NR], RW[NR], the R abscissas and weights.
//
//    Input, int NT, the number of Theta abscissas.
//
//    Input, double TA[NT], TW[NT], the THETA abscissas and weights.
//
//    Input, double ZW, the weight to use for the center.
//
//    Output, double CIRCLE_RT_SUM, the approximate integral of the function.
//
{
  int ir;
  int it;
  double quad;
  double rct;
  double result;
  double rst;
  double volume;
  double x;
  double y;

  quad = 0.0;

  if ( zw != 0.0 )
  {
    x = center[0];
    y = center[1];
    quad = quad + zw * func ( x, y );
  }

  for ( it = 0; it < nt; it++ )
  {
    rct = radius * cos ( ta[it] );
    rst = radius * sin ( ta[it] );
    for ( ir = 0; ir < nr; ir++ )
    {
      x = center[0] + ra[ir] * rct;
      y = center[1] + ra[ir] * rst;
      quad = quad + tw[it] * rw[ir] * func ( x, y );
    }
  }
  volume = circle_area_2d ( radius );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double circle_sector ( double func ( double x, double y ), double center[2], 
  double radius, double theta1, double theta2, int nr )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SECTOR approximates an integral in a circular sector.
//
//  Discussion:
//
//    A sector is contained within a circular arc and the lines joining each
//    endpoint of the arc to the center of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the center of the circle.
//
//    Input, double RADIUS, the radius of the circle.
//
//    Input, double THETA1, THETA2, the angles defining the sector.
//    The sector is measured from THETA1 to THETA2.
//
//    Input, int NR, the number of radial values used in the approximation
//    of the integral.  NR must be at least 1.  Higher values improve the
//    accuracy of the integration, at the cost of more function evaluations.
//
//    Output, double CIRCLE_SECTOR, the approximation to the integral.
//
{
  double a;
  double area;
  double b;
  double c;
  double d;
  int i;
  int j;
  int nt;
  double quad;
  double *ra;
  double result;
  double *rw;
  double *ta;
  double *tw;
  double x;
  double y;
//
//  Set the radial abscissas and weights.
//
  ra = new double[nr];
  rw = new double[nr];

  legendre_set ( nr, ra, rw );

  a = -1.0;
  b = +1.0;
  c =  0.0;
  d =  radius * radius;

  rule_adjust ( a, b, c, d, nr, ra, rw );

  for ( i = 0; i < nr; i++ )
  {
    ra[i] = sqrt ( ra[i] );
  }
  for ( i = 0; i <  nr; i++ )
  {
    rw[i] = rw[i] / radius / radius;
  }
//
//  Pick angles evenly spaced between THETA1 and THETA2, but do not
//  include the endpoints, and use a half interval for the first and last.
//
  nt = 4 * nr;

  ta = tvec_even_bracket3 ( nt, theta1, theta2 );

  tw = new double[nt];
  for ( i = 0; i < nt; i++ )
  {
    tw[i] = 1.0 / ( double ) ( nt );
  }
//
//  Approximate the integral.
//
  quad = 0.0;
  for ( i = 0; i < nr; i++ )
  {
    for ( j = 0; j < nt; j++ )
    {
      x = center[0] + ra[i] * cos ( ta[j] );
      y = center[1] + ra[i] * sin ( ta[j] );
      quad = quad + rw[i] * tw[j] * func ( x, y );
    }
  }

  area = circle_sector_area_2d ( radius, theta1, theta2 );
  result = quad * area;

  delete [] ra;
  delete [] rw;
  delete [] ta;
  delete [] tw;

  return result;
}
//****************************************************************************80

double circle_sector_area_2d ( double r, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SECTOR_AREA_2D returns the area of a circular sector in 2D.
//
//  Discussion:
//
//    A sector is contained within a circular arc and the lines joining each
//    endpoint of the arc to the center of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Input, double THETA1, THETA2, the angles of the rays
//    that delimit the sector.
//
//    Output, double CIRCLE_SECTOR_AREA_2D, the area of the sector.
//
{
  double value;

  value = 0.5 * r * r * ( theta2 - theta1 );

  return value;
}
//****************************************************************************80

double circle_triangle_area_2d ( double r, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
//
//  Discussion:
//
//    A circle triangle is formed by drawing a circular arc, and considering
//    the triangle formed by the endpoints of the arc plus the center of
//    the circle.
//
//    The normal situation is that 0 < ( THETA2 - THETA1 ) < PI.  Outside
//    this range, the triangle can actually have NEGATIVE area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Input, double THETA1, THETA2, the angles of the rays that
//    delimit the arc.
//
//    Output, double CIRCLE_TRIANGLE_AREA_2D, the (signed) area
//    of the triangle.
//
{
  double value;

  value = 0.5 * r * r * sin ( theta2 - theta1 );

  return value;
}
//****************************************************************************80

void circle_xy_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_XY_SET sets an XY quadrature rule inside the unit circle in 2D.
//
//  Integration region:
//
//    X*X + Y*Y <= 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
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
//    Frank Lether,
//    A Generalized Product Rule for the Circle,
//    SIAM Journal on Numerical Analysis,
//    Volume 8, Number 2, June 1971, pages 249-253.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule desired.
//      1, 1 point 1-st degree;
//      2, 4 point 3-rd degree, Stroud S2:3-1;
//      3, 4 point 3-rd degree, Lether #1;
//      4, 4 point 3-rd degree, Stroud S2:3-2;
//      5, 5 point 3-rd degree;
//      6, 7 point 5-th degree;
//      7, 9 point 5-th degree;
//      8, 9 point 5-th degree, Lether #2;
//      9, 12 point 7-th degree;
//     10, 16 point 7-th degree, Lether #3;
//     11, 21 point 9-th degree, Stroud S2:9-3;
//     12, 25 point 9-th degree, Lether #4 (after correcting error);
//     13, 64 point 15-th degree Gauss product rule.
//
//    Input, int ORDER, the order of the desired rule.
//
//    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas of 
//    the rule.
//
//    Output, double WEIGHT[ORDER], the ORDER weights of the rule.
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  double g;
  double h;
  int i;
  int j;
  int k;
  int nr;
  double pi = 3.141592653589793;
  double r;
  double *ra;
  double *rw;
  double s;
  double u;
  double v;
  double w;
  double w1;
  double w2;
  double w3;
  double w4;
  double w5;
  double w6;
  double w7;
  double w8;
  double w9;
  double z;

  if ( rule == 1 )
  {
    xtab[0] = 0.0;
    ytab[0] = 0.0;
    weight[0] = 1.0;
  }
  else if ( rule == 2 )
  {
    a = 0.5;
    b = 0.25;
    z = 0.0;

    xtab[0] =   a;
    xtab[1] = - a;
    xtab[2] =   z;
    xtab[3] =   z;

    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   a;
    ytab[3] = - a;

    weight[0] = b;
    weight[1] = b;
    weight[2] = b;
    weight[3] = b;
  }
  else if ( rule == 3 )
  {
    a = 0.5;
    b = 0.25;

    xtab[0] =   a;
    xtab[1] = - a;
    xtab[2] = - a;
    xtab[3] =   a;

    ytab[0] =   a;
    ytab[1] =   a;
    ytab[2] = - a;
    ytab[3] = - a;

    weight[0] = b;
    weight[1] = b;
    weight[2] = b;
    weight[3] = b;
  }
  else if ( rule == 4 )
  {
    a = sqrt ( 2.0 ) / 2.0;
    b = 0.25;

    xtab[0] =   a;
    xtab[1] = - a;
    xtab[2] = - a;
    xtab[3] =   a;

    ytab[0] =   a;
    ytab[1] =   a;
    ytab[2] = - a;
    ytab[3] = - a;

    weight[0] = b;
    weight[1] = b;
    weight[2] = b;
    weight[3] = b;
  }
  else if ( rule == 5 )
  {
    a = 1.0;
    b = 0.5;
    c = 0.125;
    z = 0.0;

    xtab[0] =   z;
    xtab[1] =   a;
    xtab[2] =   z;
    xtab[3] = - a;
    xtab[4] =   z;

    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   a;
    ytab[3] =   z;
    ytab[4] = - a;

    weight[0] = b;
    weight[1] = c;
    weight[2] = c;
    weight[3] = c;
    weight[4] = c;
  }
  else if ( rule == 6 )
  {
    a = sqrt ( 2.0 / 3.0 );
    b = sqrt ( 1.0 / 6.0 );
    c = sqrt ( 2.0 ) / 2.0;
    d = 0.125;
    e = 0.25;
    z = 0.0;

    xtab[0] =   z;
    xtab[1] =   a;
    xtab[2] = - a;
    xtab[3] =   b;
    xtab[4] = - b;
    xtab[5] =   b;
    xtab[6] = - b;

    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   z;
    ytab[3] =   c;
    ytab[4] =   c;
    ytab[5] = - c;
    ytab[6] = - c;

    weight[0] = e;
    weight[1] = d;
    weight[2] = d;
    weight[3] = d;
    weight[4] = d;
    weight[5] = d;
    weight[6] = d;
  }
  else if ( rule == 7 )
  {
    a = 0.5;
    b = 1.0;
    c = 4.0 / 24.0;
    d = 1.0 / 24.0;
    z = 0.0;

    xtab[0] =   z;
    xtab[1] =   b;
    xtab[2] = - b;
    xtab[3] =   z;
    xtab[4] =   z;
    xtab[5] =   a;
    xtab[6] = - a;
    xtab[7] = - a;
    xtab[8] =   a;

    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   z;
    ytab[3] =   b;
    ytab[4] = - b;
    ytab[5] =   a;
    ytab[6] =   a;
    ytab[7] = - a;
    ytab[8] = - a;

    weight[0] = c;
    weight[1] = d;
    weight[2] = d;
    weight[3] = d;
    weight[4] = d;
    weight[5] = c;
    weight[6] = c;
    weight[7] = c;
    weight[8] = c;
  }
  else if ( rule == 8 )
  {
    a = sqrt ( 2.0 ) / 2.0;
    b = sqrt ( 3.0 / 5.0 );
    c = sqrt ( 3.0 / 10.0 );

    w1 = 16.0 / 72.0;
    w2 =  8.0 / 72.0;
    w3 = 10.0 / 72.0;
    w4 =  5.0 / 72.0;

    z = 0.0;

    xtab[0] =   z;
    xtab[1] =   a;
    xtab[2] = - a;
    xtab[3] =   z;
    xtab[4] =   z;
    xtab[5] =   a;
    xtab[6] =   a;
    xtab[7] = - a;
    xtab[8] = - a;

    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   z;
    ytab[3] =   b;
    ytab[4] = - b;
    ytab[5] =   c;
    ytab[6] = - c;
    ytab[7] =   c;
    ytab[8] = - c;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w3;
    weight[4] = w3;
    weight[5] = w4;
    weight[6] = w4;
    weight[7] = w4;
    weight[8] = w4;
  }
  else if ( rule == 9 )
  {
    a = sqrt ( 3.0 ) / 2.0;
    b = sqrt ( ( 27.0 - 3.0 * sqrt ( 29.0 ) ) / 104.0 );
    c = sqrt ( ( 27.0 + 3.0 * sqrt ( 29.0 ) ) / 104.0 );
    u = 2.0 / 27.0;
    v = ( 551.0 + 41.0 * sqrt ( 29.0 ) ) / 6264.0;
    w = ( 551.0 - 41.0 * sqrt ( 29.0 ) ) / 6264.0;
    z = 0.0;

    xtab[0]  =   a;
    xtab[1]  = - a;
    xtab[2]  =   z;
    xtab[3]  =   z;
    xtab[4]  =   b;
    xtab[5]  = - b;
    xtab[6]  =   b;
    xtab[7]  = - b;
    xtab[8]  =   c;
    xtab[9]  =   c;
    xtab[10] = - c;
    xtab[11] = - c;

    ytab[0]  =   z;
    ytab[1]  =   z;
    ytab[2]  =   a;
    ytab[3]  = - a;
    ytab[4]  =   b;
    ytab[5]  =   b;
    ytab[6]  = - b;
    ytab[7]  = - b;
    ytab[8]  =   c;
    ytab[9]  = - c;
    ytab[10] =   c;
    ytab[11] = - c;

    weight[0]  = u;
    weight[1]  = u;
    weight[2]  = u;
    weight[3]  = u;
    weight[4]  = v;
    weight[5]  = v;
    weight[6]  = v;
    weight[7]  = v;
    weight[8]  = w;
    weight[9]  = w;
    weight[10] = w;
    weight[11] = w;
  }
  else if ( rule == 10 )
  {
    a = sqrt ( ( 3.0 - sqrt ( 5.0 ) ) / 8.0 );
    b = sqrt ( ( 15.0 + 3.0 * sqrt ( 5.0 ) 
      - 2.0 * sqrt ( 30.0 ) - 2.0 * sqrt ( 6.0 ) ) / 56.0 );
    c = sqrt ( ( 15.0 + 3.0 * sqrt ( 5.0 ) 
      + 2.0 * sqrt ( 30.0 ) + 2.0 * sqrt ( 6.0 ) ) / 56.0 );
    d = sqrt ( ( 3.0 + sqrt ( 5.0 ) ) / 8.0 );
    e = sqrt ( ( 15.0 - 3.0 * sqrt ( 5.0 ) 
      - 2.0 * sqrt ( 30.0 ) + 2.0 * sqrt ( 6.0 ) ) / 56.0 );
    f = sqrt ( ( 15.0 - 3.0 * sqrt ( 5.0 ) 
      + 2.0 * sqrt ( 30.0 ) - 2.0 * sqrt ( 6.0 ) ) / 56.0 );
    w1 = ( 90.0 + 5.0 * sqrt ( 30.0 ) + 18.0 * sqrt ( 5.0 ) 
       + 5.0 * sqrt ( 6.0 ) ) / 1440.0;
    w2 = ( 90.0 - 5.0 * sqrt ( 30.0 ) + 18.0 * sqrt ( 5.0 ) 
       - 5.0 * sqrt ( 6.0 ) ) / 1440.0;
    w3 = ( 90.0 + 5.0 * sqrt ( 30.0 ) - 18.0 * sqrt ( 5.0 ) 
       - 5.0 * sqrt ( 6.0 ) ) / 1440.0;
    w4 = ( 90.0 - 5.0 * sqrt ( 30.0 ) - 18.0 * sqrt ( 5.0 ) 
       + 5.0 * sqrt ( 6.0 ) ) / 1440.0;

    xtab[0]  =   a;
    xtab[1]  =   a;
    xtab[2]  = - a;
    xtab[3]  = - a;
    xtab[4]  =   a;
    xtab[5]  =   a;
    xtab[6]  = - a;
    xtab[7]  = - a;
    xtab[8]  =   d;
    xtab[9]  =   d;
    xtab[10] = - d;
    xtab[11] = - d;
    xtab[12] =   d;
    xtab[13] =   d;
    xtab[14] = - d;
    xtab[15] = - d;

    ytab[0]  =   b;
    ytab[1]  = - b;
    ytab[2]  =   b;
    ytab[3]  = - b;
    ytab[4]  =   c;
    ytab[5]  = - c;
    ytab[6]  =   c;
    ytab[7]  = - c;
    ytab[8]  =   e;
    ytab[9]  = - e;
    ytab[10] =   e;
    ytab[11] = - e;
    ytab[12] =   f;
    ytab[13] = - f;
    ytab[14] =   f;
    ytab[15] = - f;

    weight[0]  = w1;
    weight[1]  = w1;
    weight[2]  = w1;
    weight[3]  = w1;
    weight[4]  = w2;
    weight[5]  = w2;
    weight[6]  = w2;
    weight[7]  = w2;
    weight[8]  = w3;
    weight[9]  = w3;
    weight[10] = w3;
    weight[11] = w3;
    weight[12] = w4;
    weight[13] = w4;
    weight[14] = w4;
    weight[15] = w4;
  }
  else if ( rule == 11 )
  {
    xtab[0] = 0.0;
    ytab[0] = 0.0;
    weight[0] = 1.0 / 9.0;
    
    for ( i = 1; i < 11; i++ )
    {
      weight[i] = ( 16.0 + sqrt ( 6.0 ) ) / 360.0;
    }
    for ( i = 11; i < 21; i++ )
    {
      weight[i] = ( 16.0 - sqrt ( 6.0 ) ) / 360.0;
    }

    r = sqrt ( ( 6.0 - sqrt ( 6.0 ) ) / 10.0 );

    for ( i = 0; i < 10; i++ )
    {
      a = 2.0 * pi * ( double ) ( i ) / 10.0;
      xtab[i+1] = r * cos ( a );
      ytab[i+1] = r * sin ( a );
    }

    r = sqrt ( ( 6.0 + sqrt ( 6.0 ) ) / 10.0 );

    for ( i = 0; i < 10; i++ )
    {
      a = 2.0 * pi * ( double ) ( i ) / 10.0;
      xtab[i+11] = r * cos ( a );
      ytab[i+11] = r * sin ( a );
    }
  }
//
//  There was apparently a misprint in the Lether paper.  The quantity
//  which here reads "322" was printed there as "332".
//
  else if ( rule == 12 )
  {
    a = 0.5;
    b = sqrt ( 3.0 ) / 2.0;
    c = sqrt ( ( 35.0 + 2.0 * sqrt ( 70.0 ) ) / 252.0 );
    d = sqrt ( ( 35.0 - 2.0 * sqrt ( 70.0 ) ) / 252.0 );
    e = sqrt ( ( 35.0 + 2.0 * sqrt ( 70.0 ) ) / 84.0 );
    f = sqrt ( ( 35.0 - 2.0 * sqrt ( 70.0 ) ) / 84.0 );
    g = sqrt ( ( 35.0 + 2.0 * sqrt ( 70.0 ) ) / 63.0 );
    h = sqrt ( ( 35.0 - 2.0 * sqrt ( 70.0 ) ) / 63.0 );

    w1 = 64.0 / 675.0;
    w2 = 16.0 / 225.0;
    w3 = 16.0 / 675.0;
    w4 = ( 322.0 - 13.0 * sqrt ( 70.0 ) ) / 21600.0;
    w5 = ( 322.0 + 13.0 * sqrt ( 70.0 ) ) / 21600.0;
    w6 = ( 322.0 - 13.0 * sqrt ( 70.0 ) ) / 7200.0;
    w7 = ( 322.0 + 13.0 * sqrt ( 70.0 ) ) / 7200.0;
    w8 = ( 322.0 - 13.0 * sqrt ( 70.0 ) ) / 5400.0;
    w9 = ( 322.0 + 13.0 * sqrt ( 70.0 ) ) / 5400.0;
    z = 0.0;

    xtab[0]  =   z;
    xtab[1]  =   a;
    xtab[2]  = - a;
    xtab[3]  =   b;
    xtab[4]  = - b;
    xtab[5]  =   b;
    xtab[6]  =   b;
    xtab[7]  = - b;
    xtab[8]  = - b;
    xtab[9]  =   b;
    xtab[10] =   b;
    xtab[11] = - b;
    xtab[12] = - b;
    xtab[13] =   a; 
    xtab[14] =   a;
    xtab[15] = - a;
    xtab[16] = - a;
    xtab[17] =   a;
    xtab[18] =   a;
    xtab[19] = - a;
    xtab[20] = - a;
    xtab[21] =   z;
    xtab[22] =   z;
    xtab[23] =   z;
    xtab[24] =   z;

    ytab[0]  =   z;
    ytab[1]  =   z;
    ytab[2]  =   z;
    ytab[3]  =   z;
    ytab[4]  =   z;
    ytab[5]  =   c;
    ytab[6]  = - c;
    ytab[7]  =   c;
    ytab[8]  = - c;
    ytab[9]  =   d;
    ytab[10] = - d;
    ytab[11] =   d;
    ytab[12] = - d;
    ytab[13] =   e; 
    ytab[14] = - e;
    ytab[15] =   e;
    ytab[16] = - e;
    ytab[17] =   f;
    ytab[18] = - f;
    ytab[19] =   f;
    ytab[20] = - f;
    ytab[21] =   g;
    ytab[22] = - g;
    ytab[23] =   h;
    ytab[24] = - h;

    weight[0]  = w1;
    weight[1]  = w2;
    weight[2]  = w2;
    weight[3]  = w3;
    weight[4]  = w3;
    weight[5]  = w4;
    weight[6]  = w4;
    weight[7]  = w4;
    weight[8]  = w4;
    weight[9]  = w5;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w6;
    weight[17] = w7;
    weight[18] = w7;
    weight[19] = w7;
    weight[20] = w7;
    weight[21] = w8;
    weight[22] = w8;
    weight[23] = w9;
    weight[24] = w9;
  }
  else if ( rule == 13 )
  {
    nr = 4;
    ra = new double[nr];
    rw = new double[nr];

    legendre_set ( nr, ra, rw );

    a = -1.0;
    b = +1.0;
    c =  0.0;
    d = +1.0;
    rule_adjust ( a, b, c, d, nr, ra, rw );

    for ( i = 0; i < nr; i++ )
    {
      ra[i] = sqrt ( ra[i] );
    }

    i = 0;
    for ( j = 0; j < 16; j++ )
    {
      c = cos ( pi * ( double ) ( j ) / 8.0 );
      s = sin ( pi * ( double ) ( j ) / 8.0 );

      for ( k = 0; k < nr; k++ )
      {
        xtab[i] = c * ra[k];
        ytab[i] = s * ra[k];
        weight[i] = rw[k] / 16.0;
        i = i + 1;
      }
    }
    delete [] ra;
    delete [] rw;
  }
  else
  {
    cerr << "\n";
    cerr << "CIRCLE_XY_SET - Fatal error!\n";
    cerr << "  There is no rule of index " << rule << "\n";
    exit ( 1 );
  }

  return;
}

//****************************************************************************80

int circle_xy_size ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_XY_SIZE sizes an XY quadrature rule inside the unit circle in 2D.
//
//  Integration region:
//
//    X*X + Y*Y <= 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
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
//    Frank Lether,
//    A Generalized Product Rule for the Circle,
//    SIAM Journal on Numerical Analysis,
//    Volume 8, Number 2, June 1971, pages 249-253.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule desired.
//      1, 1 point 1-st degree;
//      2, 4 point 3-rd degree, Stroud S2:3-1;
//      3, 4 point 3-rd degree, Lether #1;
//      4, 4 point 3-rd degree, Stroud S2:3-2;
//      5, 5 point 3-rd degree;
//      6, 7 point 5-th degree;
//      7, 9 point 5-th degree;
//      8, 9 point 5-th degree, Lether #2;
//      9, 12 point 7-th degree;
//     10, 16 point 7-th degree, Lether #3;
//     11, 21 point 9-th degree, Stroud S2:9-3;
//     12, 25 point 9-th degree, Lether #4 (after correcting error);
//     13, 64 point 15-th degree Gauss product rule.
//
//    Output, int CIRCLE_XY_SIZE, the order of the desired rule.
//
{
  int order;

  if ( rule == 1 )
  {
    order = 1;
  }
  else if ( rule == 2 )
  {
    order = 4;
  }
  else if ( rule == 3 )
  {
    order = 4;
  }
  else if ( rule == 4 )
  {
    order = 4;
  }
  else if ( rule == 5 )
  {
    order = 5;
  }
  else if ( rule == 6 )
  {
    order = 7;
  }
  else if ( rule == 7 )
  {
    order = 9;
  }
  else if ( rule == 8 )
  {
    order = 9;
  }
  else if ( rule == 9 )
  {
    order = 12;
  }
  else if ( rule == 10 )
  {
    order = 16;
  }
  else if ( rule == 11 )
  {
    order = 21;
  }
  else if ( rule == 12 )
  {
    order = 25;
  }
  else if ( rule == 13 )
  {
    order = 64;
  }
  else
  {
    order = -1;
    cerr << "\n";
    cerr << "CIRCLE_XY_SIZE - Fatal error!\n";
    cerr << "  There is no rule of index " << rule << "\n";
    exit ( 1 );
  }

  return order;
}
//****************************************************************************80

double circle_xy_sum ( double func ( double x, double y ), double center[2], 
  double r, int order, double xtab[], double ytab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_XY_SUM applies an XY quadrature rule inside a circle in 2D.
//
//  Integration region:
//
//    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= R * R.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the coordinates of the center of 
//    the circle.
//
//    Input, double R, the radius of the circle.
//
//    Input, int ORDER, the order of the rule.  The rule is
//    assumed to be defined on the unit circle.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the XY
//    coordinates of the abscissas of the quadrature rule for the unit circle.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;
  double x;
  double y;

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    x = center[0] + r * xtab[i];
    y = center[1] + r * ytab[i];
    quad = quad + weight[i] * func ( x, y );
  }

  volume = circle_area_2d ( r );
  result = quad * volume;

  return result;
}
//****************************************************************************80

void cn_geg_00_1 ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_00_1 implements the midpoint rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int k;
  double volume;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_00_1 - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  expon = 0;
  volume = c1_geg_monomial_integral ( alpha, expon );
  volume = pow ( volume, n );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_geg_00_1_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_00_1_SIZE sizes the midpoint rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, int CN_GEG_00_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_00_1_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void cn_geg_01_1 ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_01_1 implements a precision 1 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[[O], the weights.
//
{
  int expon;
  int i;
  int j;
  int k;
  double value1;
  double value2;
  double volume;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_01_1 - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  expon = 0;
  value1 = c1_geg_monomial_integral ( alpha, expon );
  volume = pow ( value1, n );

  expon = 1;
  value2 = c1_geg_monomial_integral ( alpha, expon );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / value1;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_geg_01_1_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_01_1_SIZE sizes a precision 1 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, int CN_GEG_01_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_01_1_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void cn_geg_02_xiu ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_02_XIU implements the Xiu precision 2 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_02_XIU - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 1.0 / ( 2.0 * alpha + 3.0 );

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = c1_geg_monomial_integral ( alpha, expon );
  volume = pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_geg_02_xiu_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_02_XIU_SIZE sizes the Xiu rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, int CN_GEG_02_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_02_XIU_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  o = n + 1;

  return o;
}
//****************************************************************************80

void cn_geg_03_xiu ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_03_XIU implements the Xiu precision 3 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int expon;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_03_XIU - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  expon = 0;
  volume = c1_geg_monomial_integral ( alpha, expon );
  volume = pow ( volume, n );

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * j ) * pi / ( double ) ( n );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg ) / sqrt ( 2.0 * alpha + 3.0 );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg ) / sqrt ( 2.0 * alpha + 3.0 );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = sqrt ( 2.0 ) * r8_mop ( j ) / sqrt ( 2.0 * alpha + 3.0 );
      if ( n == 1 )
      {
        x[i+j*n] = x[i+j*n] / sqrt ( 2.0 );
      }
      i = i + 1;
    }
  }

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int cn_geg_03_xiu_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_GEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the parameter.
//    -1.0 < ALPHA.
//
//    Output, int CN_GEG_03_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_03_XIU_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  o = 2 * n;

  return o;
}
//****************************************************************************80

double cn_geg_monomial_integral ( int n, double alpha, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_GEG_MONOMIAL_INTEGRAL: integral of monomial with Gegenbauer weight on CN.
//
//  Discussion:
//
//    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
//
//      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
//
//    with -1.0 < alpha.
//
//    value = integral ( CN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i)^2)^alpha dx(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of (1-X).
//    -1.0 < ALPHA.
//
//    Input, int EXPON[N], the exponents.
//
//    Output, double CN_GEG_MONOMIAL_INTEGRA, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_GEG_MONOMIAL_INTEGRAL - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = c1_geg_monomial_integral ( alpha, expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void cn_jac_00_1 ( int n, double alpha, double beta, int o, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_00_1 implements the midpoint rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int k;
  double pi = 3.141592653589793;
  double volume;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_00_1 - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_00_1 - Fatal error!\n";
    cerr << "  BETA <= -1.0\n";
    exit ( 1 );
  }

  expon = 0;
  volume = c1_jac_monomial_integral ( alpha, beta, expon );
  volume = pow ( volume, n );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_jac_00_1_size ( int n, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_00_1_SIZE sizes the midpoint rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Output, int CN_JAC_00_1_SIZE the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_00_1_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_00_1_SIZE - Fatal error!\n";
    cerr << "  BETA <= -1.0\n";
    exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void cn_jac_01_1 ( int n, double alpha, double beta, int o, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_01_1 implements a precision 1 rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double value1;
  double value2;
  double volume;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_01_1 - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_01_1 - Fatal error!\n";
    cerr << "  BETA <= -1.0\n";
    exit ( 1 );
  }

  expon = 0;
  value1 = c1_jac_monomial_integral ( alpha, beta, expon );
  volume = pow ( value1, n );

  expon = 1;
  value2 = c1_jac_monomial_integral ( alpha, beta, expon );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / value1;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_jac_01_1_size ( int n, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_01_1_SIZE sizes a precision 1 rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha. 
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Output, int CN_JAC_01_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_01_1_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_01_1_SIZE - Fatal error!\n";
    cerr << "  BETA <= -1.0\n";
    exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void cn_jac_02_xiu ( int n, double alpha, double beta, int o, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_02_XIU implements the Xiu precision 2 rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_02_XIU - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_02_XIU - Fatal error!\n";
    cerr << "  BETA <= -1.0\n";
    exit ( 1 );
  }

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = ( alpha + beta + 2.0 ) / 2.0;
  delta0 = ( alpha - beta ) / 2.0;
  c1 = 2.0 * ( alpha + 1.0 ) * ( beta + 1.0 ) / ( alpha + beta + 3.0 )
    / ( alpha + beta + 2.0 );

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = c1_jac_monomial_integral ( alpha, beta, expon );
  volume = pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_jac_02_xiu_size ( int n, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_02_XIU_SIZE sizes the Xiu rule for region CN_JAC.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN is the cube [-1,+1]^N with the Jacobi (beta) weight function
//
//      w(alpha,beta;x) = product ( 1 <= i <= n ) (1-x(i))^beta (1+x(i))^alpha.
//
//    with -1 < alpha, -1 < beta.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, BETA, the parameters.
//    -1.0 < ALPHA, -1.0 < BETA.
//
//    Output, int CN_JAC_02_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_02_XIU_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "CN_JAC_02_XIU_SIZE - Fatal error!\n";
    cerr << "  BETA <= -1.0\n";
    exit ( 1 );
  }

  o = n + 1;

  return o;
}
//****************************************************************************80

double cn_jac_monomial_integral ( int n, double alpha, double beta, 
  int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_JAC_MONOMIAL_INTEGRAL: integral of a monomial with Jacobi weight over CN.
//
//  Discussion:
//
//    value = integral ( CN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) (1-x(i))^alpha (1+x(i))^beta dx(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of (1-X) in the weight factor.
//
//    Input, double BETA, the exponent of (1+X) in the weight factor.
//
//    Input, int EXPON[N], the exponents.
//
//    Output, double CN_JAC_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = c1_jac_monomial_integral ( alpha, beta, expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void cn_leg_01_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_01_1 implements the midpoint rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int k;
  double volume;

  expon = 0;
  volume = c1_leg_monomial_integral ( expon );
  volume = pow ( volume, n );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = volume;

  return;
}
//****************************************************************************80

int cn_leg_01_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_01_1_SIZE sizes the midpoint rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_01_1_SIZE, the order.
//
{
  int o;

  o = 1;

  return o;
}
//****************************************************************************80

void cn_leg_02_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_02_XIU implements the Xiu precision 2 rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = 1.0;
  delta0 = 0.0;
  c1 = 1.0 / 3.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = c1_leg_monomial_integral ( expon );
  volume = pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int cn_leg_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_02_XIU_SIZE sizes the Xiu rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_02_XIU_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

void cn_leg_03_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_1 implements the Stroud rule CN:3-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//    The necessary treatment of the final coordinate of points when
//    N is odd seems to vary from what Stroud declares! 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int expon;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;

  expon = 0;
  volume = c1_leg_monomial_integral ( expon );
  volume = pow ( volume, n );

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * ( j + 1 ) ) * pi / ( double ) ( n );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg ) / sqrt ( 3.0 );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg ) / sqrt ( 3.0 );
      i = i + 1;
    }
//
//  The following code does not correspond to what Stroud declares.
//
    if ( i < n )
    {
      if ( n == 1 )
      {
        x[i+j*n] =                r8_mop ( j + 1 ) / sqrt ( 3.0 );
      }
      else
      {
        x[i+j*n] = sqrt ( 2.0 ) * r8_mop ( j + 1 ) / sqrt ( 3.0 );
      }
      i = i + 1;
    }
  }

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_leg_03_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_1_SIZE sizes the Stroud rule CN:3-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_03_1_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void cn_leg_03_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_XIU implements the Xiu precision 3 rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int expon;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;

  expon = 0;
  volume = c1_leg_monomial_integral ( expon );
  volume = pow ( volume, n );

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * ( j + 1 ) ) * pi / ( double ) ( n );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg ) / sqrt ( 3.0 );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg ) / sqrt ( 3.0 );
      i = i + 1;
    }

    if ( i < n )
    {
      if ( n == 1 )
      {
        x[i+j*n] =                r8_mop ( j + 1 ) / sqrt ( 3.0 );
      }
      else
      {
        x[i+j*n] = sqrt ( 2.0 ) * r8_mop ( j + 1 ) / sqrt ( 3.0 );
      }
      i = i + 1;
    }
  }

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int cn_leg_03_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_03_XIU_SIZE sizes the Xiu precision 3 rule for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_03_XIU_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void cn_leg_05_1 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_1 implements the Stroud rule CN:5-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N^2 + N + 2.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    N must be 4, 5, or 6.
//
//    Input, int OPTION, is only used in case N = 4 or 5.
//    In that case, OPTION should be 1 or 2 to select the
//    two available variants of the rule.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double arg;
  double b;
  double c;
  double eta;
  int expon;
  double gamma;
  int i;
  int i1;
  int i2;
  int k;
  double lambda;
  double mu;
  double volume;
  double xsi;

  if ( n < 4 || 6 < n )
  {
    cerr << "\n";
    cerr << "CN_LEG_05_1 - Fatal error!\n";
    cerr << "  The value of N must be 4, 5, or 6.\n";
    exit ( 1 );
  }

  if ( n == 4 || n == 5 )
  {
    if ( option < 1 || 2 < option )
    {
      cerr << "\n";
      cerr << "CN_LEG_05_1 - Fatal error!\n";
      cerr << "  When N = 4 or 5, the value of OPTION must be 1 or 2.\n";
      exit ( 1 );
    }
  }

  expon = 0;
  volume = c1_leg_monomial_integral ( expon );
  volume = pow ( volume, n );

  if ( n == 4 && option == 1 )
  {
    eta    =   0.778984505799815E+00;
    lambda =   1.284565137874656E+00;
    xsi =    - 0.713647298819253E+00;
    mu =     - 0.715669761974162E+00;
    gamma =    0.217089151000943E+00;
    a =        0.206186096875899E-01 * volume;
    b =        0.975705820221664E-02 * volume;
    c =        0.733921929172573E-01 * volume;
  }
  else if ( n == 4 && option == 2 )
  {
    eta    =   0.546190755827425E+00;
    lambda =   0.745069130115661E+00;
    xsi =    - 0.413927294508700E+00;
    mu =     - 0.343989637454535E+00;
    gamma =    1.134017894600344E+00;
    a =        0.853094758323323E-01 * volume;
    b =        0.862099000096395E-01 * volume;
    c =        0.116418206881849E-01 * volume;
  }
  else if ( n == 5 && option == 1 )
  {
    eta    =   0.522478547481276E+00;
    lambda =   0.936135175985774E+00;
    xsi =    - 0.246351362101519E+00;
    mu =     - 0.496308106093758E+00;
    gamma =    0.827180176822930E+00;
    a =        0.631976901960153E-01 * volume;
    b =        0.511464127430166E-01 * volume;
    c =        0.181070246088902E-01 * volume;
  }
  else if ( n == 5 && option == 2 )
  {
    eta    =   0.798317301388741E+00;
    lambda =   0.637344273885728E+00;
    xsi =    - 0.455245909918377E+00;
    mu =     - 1.063446229997311E+00;
    gamma =    0.354482076665770E+00;
    a =        0.116952384292206E-01 * volume;
    b =        0.701731258612708E-01 * volume;
    c =        0.137439132264426E-01 * volume;
  }
  else if ( n == 6 )
  {
    eta    =   0.660225291773525E+00;
    lambda =   1.064581294844754E+00;
    xsi =      0.000000000000000E+00;
    mu =     - 0.660225291773525E+00;
    gamma =    0.660225291773525E+00;
    a =        0.182742214532872E-01 * volume;
    b =        0.346020761245675E-01 * volume;
    c =        0.182742214532872E-01 * volume;
  }
  k = -1;

  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = eta;
  }
  w[k] = a;

  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = - eta;
  }
  w[k] = a;

  for ( i1 = 0; i1 < n; i1++ )
  {
    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = xsi;
    }
    x[i1+k*n] = lambda;
    w[k] = b;
  }

  for ( i1 = 0; i1 < n; i1++ )
  {
    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = - xsi;
    }
    x[i1+k*n] = - lambda;
    w[k] = b;
  }

  for ( i1 = 0; i1 < n - 1; i1++ )
  {
    for ( i2 = i1 + 1; i2 < n; i2++ )
    {
      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = gamma;
      }
      x[i1+k*n] = mu;
      x[i2+k*n] = mu;
      w[k] = c;
    }
  }

  for ( i1 = 0; i1 < n - 1; i1++ )
  {
    for ( i2 = i1 + 1; i2 < n; i2++ )
    {
      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = - gamma;
      }
      x[i1+k*n] = - mu;
      x[i2+k*n] = - mu;
      w[k] = c;
    }
  }

  return;
}
//****************************************************************************80

int cn_leg_05_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_1_SIZE sizes the Stroud rule CN:5-1 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N^2 + N + 2.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_05_1_SIZE, the order.
//
{
  int o;

  o = n * n + n + 2;

  return o;
}
//****************************************************************************80

void cn_leg_05_2 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_2 implements the Stroud rule CN:5-2 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 N^2 + 1.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    N must be at least 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double b0;
  double b1;
  double b2;
  int expon;
  int i;
  int i1;
  int i2;
  int k;
  double r;
  double volume;

  if ( n < 2 )
  {
    cerr << "\n";
    cerr << "CN_LEG_05_2 - Fatal error!\n";
    cerr << "  N must be at least 2.\n";
    exit ( 1 );
  }

  expon = 0;
  volume = c1_leg_monomial_integral ( expon );
  volume = pow ( volume, n );

  b0 = ( double ) ( 25 * n * n - 115 * n + 162 ) * volume / 162.0;
  b1 = ( double ) ( 70 - 25 * n ) * volume / 162.0;
  b2 = 25.0 * volume / 324.0;

  r = sqrt ( 3.0 / 5.0 );

  k = - 1;

  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = 0.0;
  }
  w[k] = b0;

  for ( i1 = 0; i1 < n; i1++ )
  {
    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = 0.0;
    }
    x[i1+k*n] = + r;
    w[k] = b1;

    k = k + 1;
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = 0.0;
    }
    x[i1+k*n] = - r;
    w[k] = b1;
  }

  for ( i1 = 0; i1 < n - 1; i1++ )
  {
    for ( i2 = i1 + 1; i2 < n; i2++ )
    {
      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = + r;
      x[i2+k*n] = + r;
      w[k] = b2;

      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = + r;
      x[i2+k*n] = - r;
      w[k] = b2;

      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = - r;
      x[i2+k*n] = + r;
      w[k] = b2;

      k = k + 1;
      for ( i = 0; i < n; i++ )
      {
        x[i+k*n] = 0.0;
      }
      x[i1+k*n] = - r;
      x[i2+k*n] = - r;
      w[k] = b2;
    }
  }
  return;
}
//****************************************************************************80

int cn_leg_05_2_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_05_2_SIZE sizes the Stroud rule CN:5-2 for region CN_LEG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 N^2 + 1.
//
//    The rule has precision P = 5.
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int CN_LEG_05_2_SIZE, the order.
//
{
  int o;

  o = 2 * n * n + 1;

  return o;
}
//****************************************************************************80

double cn_leg_monomial_integral ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    CN_LEG_MONOMIAL_INTEGRAL: integral of monomial with Legendre weight on CN.
//
//  Discussion:
//
//    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
//
//      w(x) = 1.
//
//    value = integral ( CN ) product ( 1 <= i <= n ) x(I)^expon(i) dx(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int EXPON(N), the exponents.
//
//    Output, double CN_LEG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = c1_leg_monomial_integral ( expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

double cone_unit_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    CONE_UNIT_3D approximates an integral inside the unit cone in 3D.
//
//  Integration Region:
//
//      X^2 + Y^2 <= 1 - Z  
//
//    and
//
//      0 <= Z <= 1.
//
//  Discussion:
//
//    An 48 point degree 7 formula, Stroud CN:S2:7-1, is used.
//
//    (There is a typographical error in the S2:7-1 formula for B3.)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which evaluates F(X,Y,Z).
//
//    Output, double CONE_UNIT_3D, the approximate integral of the function.
//
{
  double a;
  double b;
  double c;
  double h;
  int i;
  double quad;
  double r;
  double result;
  double u[4] = {
    0.04850054945, 0.2386007376, 
    0.5170472951,  0.7958514179 };
  double volume;
  double w1[4] = {
    0.1108884156,  0.1434587878, 
    0.06863388717, 0.01035224075 };
  double w2[3];
  double x;
  double y;
  double z;

  a = sqrt ( 3.0 ) / 2.0;
  b = sqrt ( ( 27.0 - 3.0 * sqrt ( 29.0 ) ) / 104.0 );
  c = sqrt ( ( 27.0 + 3.0 * sqrt ( 29.0 ) ) / 104.0 );
  w2[0] = 2.0 / 9.0;
  w2[1] = 3.0 * ( 551.0 + 4.0 * sqrt ( 29.0 ) ) / 6264.0;
  w2[2] = 3.0 * ( 551.0 - 4.0 * sqrt ( 29.0 ) ) / 6264.0;

  quad = 0.0;

  for ( i = 0; i < 4; i++ )
  {
    x = a * ( 1.0 - u[i] );
    y = 0.0;
    z = u[i];
    quad = quad + w1[i] * w2[0] * func ( x, y, z );

    x = -a * ( 1.0 - u[i] );
    y = 0.0;
    z = u[i];
    quad = quad + w1[i] * w2[0] * func ( x, y, z );

    x = 0.0;
    y = a * ( 1.0 - u[i] );
    z = u[i];
    quad = quad + w1[i] * w2[0] * func ( x, y, z );

    x = 0.0;
    y = -a * ( 1.0 - u[i] );
    z = u[i];
    quad = quad + w1[i] * w2[0] * func ( x, y, z );
  }

  for ( i = 0; i < 4; i++ )
  {
    x =  b * ( 1.0 - u[i] );
    y =  b * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[1] * func ( x, y, z );

    x = -b * ( 1.0 - u[i] );
    y =  b * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[1] * func ( x, y, z );

    x = -b * ( 1.0 - u[i] );
    y = -b * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[1] * func ( x, y, z );

    x =  b * ( 1.0 - u[i] );
    y = -b * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[1] * func ( x, y, z );

    x =  c * ( 1.0 - u[i] );
    y =  c * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[2] * func ( x, y, z );

    x = -c * ( 1.0 - u[i] );
    y =  c * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[2] * func ( x, y, z );

    x = -c * ( 1.0 - u[i] );
    y = -c * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[2] * func ( x, y, z );

    x =  c * ( 1.0 - u[i] );
    y = -c * ( 1.0 - u[i] );
    z =  u[i];
    quad = quad + w1[i] * w2[2] * func ( x, y, z );

  }
  r = 1.0;
  h = 1.0;

  volume = cone_volume_3d ( r, h );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double cone_volume_3d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    CONE_VOLUME_3D returns the volume of a cone in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the base of the cone.
//
//    Input, double H, the height of the cone.
//
//    Output, double CONE_VOLUME_3D, the volume of the cone.
//
{
  double pi = 3.141592653589793;
  double value;

  value = ( pi / 3.0 ) * h * r * r;

  return value;
}
//****************************************************************************80

double cube_shell_nd ( double func ( int n, double x[] ), int n, double r1, 
  double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_SHELL_ND approximates an integral inside a cubic shell in N dimensions.
//
//  Integration region:
//
//    R1 <= abs ( X(1:N) ) <= R2
//
//  Discussion:
//
//    An N*2^N point third degree formula is used, Stroud number CNSHELL:3-4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//  user supplied function.
//
//    Input, int N, the dimension of the space.
//
//    Input, double R1, R2, the inner and outer radii of the cubical
//    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
//    2*R1.
//
//    Output, double CUBE_SHELL_ND, the approximate integral of the function.
//
{
  bool done;
  int i;
  int j;
  double quad;
  double rmax;
  double rmin;
  double result;
  double u;
  double v;
  double volume;
  double *x;

  if ( r1 == r2 )
  {
    result = 0.0;
    return result;
  }

  rmax = r8_max ( r1, r2 );
  rmin = r8_min ( r1, r2 );  

  u = sqrt ( ( double ) ( n ) * ( pow ( rmax, n + 2 ) - pow ( rmin, n + 2 ) ) 
    / ( ( double ) ( n + 2 ) * ( pow ( rmax, n ) - pow ( rmin, n ) ) ) );

  v = u / sqrt ( 3.0 );

  x = new double[n];

  quad = 0.0;
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = v;
    }
    x[i] = u;

    for ( ; ; )
    {
      quad = quad + func ( n, x );

      done = r8vec_mirror_next ( n, x );

      if ( done )
      {
        break;
      }
    }
  }

  quad = quad / ( double ) ( n * i4_power ( 2, n ) );

  volume = cube_shell_volume_nd ( n, r1, r2 );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double cube_shell_volume_nd ( int n, double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_SHELL_VOLUME_ND computes the volume of a cubic shell in ND.
//
//  Integration region:
//
//    R1 <= abs ( X(1:N) ) <= R2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double R1, R2, the inner and outer radii of the cubic
//    shell.  The outer cube is of side 2*R2, the inner, missing cube of side
//    2*R1.
//
//    Output, double CUBE_SHELL_VOLUME_ND, the volume of the cubic
//    shell.
//
{
  double value;


  value = ( pow ( r2, n ) - pow ( r1, n ) ) * i4_power ( 2, n );

  return value;
}
//****************************************************************************80

double cube_unit_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_UNIT_3D approximates an integral inside the unit cube in 3D.
//
//  Integration region:
//
//      -1 <= X <= 1,
//    and
//      -1 <= Y <= 1,
//    and
//      -1 <= Z <= 1.
//
//  Discussion:
//
//    An 8 point third degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the
//    user supplied routine to evaluate F(X,Y,Z).
//
//    Output, double CUBE_UNIT_3D, the approximate integral of the function.
//
{
  double quad;
  double result;
  double s;
  double volume;
  double w;
  double x;
  double y;
  double z;

  s = 1.0 / sqrt ( 3.0 );
  w = 1.0 / 8.0;

  x = s;
  y = s;
  z = s;

  quad = w * ( 
      func (  x,  y,  z ) + func (  x,  y, -z ) 
    + func (  x, -y,  z ) + func (  x, -y, -z ) 
    + func ( -x,  y,  z ) + func ( -x,  y, -z ) 
    + func ( -x, -y,  z ) + func ( -x, -y, -z ) );

  volume = cube_unit_volume_nd ( 3 );
  result = quad * volume;

  return result;
}
//****************************************************************************80

void cube_unit_nd ( double func ( int n, double x[] ), double qa[], 
  double qb[], int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_UNIT_ND approximates an integral inside the unit cube in ND.
//
//  Integration region:
//
//    -1 <= X(1:N) <= 1
//
//  Discussion:
//
//    A K**N point product formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Lyness, BJJ McHugh,
//    Integration Over Multidimensional Hypercubes, 
//    A Progressive Procedure,
//    The Computer Journal,
//    Volume 6, 1963, pages 264-270.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Output, double QA[K], QB[K], two sets of estimates for
//    the integral.  The QB entries are obtained from the
//    QA entries by Richardson extrapolation, and QB(K) is
//    the best estimate for the integral.
//
//    Input, int N, the dimension of the cube.
//
//    Input, int K, the highest order of integration, and the order
//    of Richardson extrapolation.  K can be no greater than 10.
//
{
  double g[10*10];
  int i;
  int j;
  int kmax = 10;

  g[0+0*10] =  1.0E+00;
  g[1+0*10] = -0.3333333333333E+00;
  g[1+1*10] =  0.1333333333333E+01;
  g[2+0*10] =  0.4166666666667E-01;
  g[2+1*10] = -0.1066666666667E+01;
  g[2+2*10] =  0.2025000000000E+01;
  g[3+0*10] = -0.2777777777778E-02;
  g[3+1*10] =  0.3555555555556E+00;
  g[3+2*10] = -0.2603571428571E+01;
  g[3+3*10] =  0.3250793650794E+01;
  g[4+0*10] =  0.1157407407407E-03;
  g[4+1*10] = -0.6772486772487E-01;
  g[4+2*10] =  0.1464508928571E+01;
  g[4+3*10] = -0.5779188712522E+01;
  g[4+4*10] =  0.5382288910935E+01;
  g[5+0*10] = -0.3306878306878E-05;
  g[5+1*10] =  0.8465608465608E-02;
  g[5+2*10] = -0.4881696428571E+00;
  g[5+3*10] =  0.4623350970018E+01;
  g[5+4*10] = -0.1223247479758E+02;
  g[5+5*10] =  0.9088831168831E+01;
  g[6+0*10] =  0.6889329805996E-07;
  g[6+1*10] = -0.7524985302763E-03;
  g[6+2*10] =  0.1098381696429E+00;
  g[6+3*10] = -0.2241624712736E+01;
  g[6+4*10] =  0.1274216124748E+02;
  g[6+5*10] = -0.2516907092907E+02;
  g[6+6*10] =  0.1555944865432E+02;
  g[7+0*10] = -0.1093544413650E-08;
  g[7+1*10] =  0.5016656868509E-04;
  g[7+2*10] = -0.1797351866883E-01;
  g[7+3*10] =  0.7472082375786E+00;
  g[7+4*10] = -0.8168052081717E+01;
  g[7+5*10] =  0.3236023405166E+02;
  g[7+6*10] = -0.5082753227079E+02;
  g[7+7*10] =  0.2690606541646E+02;
  g[8+0*10] =  0.1366930517063E-10;
  g[8+1*10] = -0.2606055516108E-05;
  g[8+2*10] =  0.2246689833604E-02;
  g[8+3*10] = -0.1839281815578E+00;
  g[8+4*10] =  0.3646451822195E+01;
  g[8+5*10] = -0.2588818724133E+02;
  g[8+6*10] =  0.7782965878964E+02;
  g[8+7*10] = -0.1012934227443E+03;
  g[8+8*10] =  0.4688718347156E+02;
  g[9+0*10] = -0.1380737896023E-12;
  g[9+1*10] =  0.1085856465045E-06;
  g[9+2*10] = -0.2222000934334E-03;
  g[9+3*10] =  0.3503393934435E-01;
  g[9+4*10] = -0.1215483940732E+01;
  g[9+5*10] =  0.1456210532325E+02;
  g[9+6*10] = -0.7477751530769E+02;
  g[9+7*10] =  0.1800771959898E+03;
  g[9+8*10] = -0.1998874663788E+03;
  g[9+9*10] =  0.8220635246624E+02;

  if ( kmax < k ) 
  {
    cerr << "\n";
    cerr << "CUBE_UNIT_ND - Fatal error!\n";
    cerr << "  K must be no greater than KMAX = " << kmax << "\n";
    cerr << "  but the input K is " << k << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < k; i++ )
  {
    qa[i] = qmdpt ( func, n, i+1 );
  }

  qb[0] = qa[0];

  for ( i = 1; i < k; i++ )
  {
     qb[i] = 0.0;
     for ( j = 0; j <= i; j++ )
     {
       qb[i] = qb[i] + g[i+j*10] * qa[j];
     }
  }
  return;
}
//****************************************************************************80

double cube_unit_volume_nd ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_UNIT_VOLUME_ND returns the volume of the unit cube in ND.
//
//  Integration region:
//
//    -1 <= X(1:N) <= 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Output, double CUBE_UNIT_VOLUME_ND, the volume of the unit
//    cube in ND.
//
{
  double value;

  value = pow ( 2.0, n );

  return value;
}
//****************************************************************************80

double ellipse_area_2d ( double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
//
//  Integration region:
//
//    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R1, R2, the major and minor semi-axes.
//
//    Output, double ELLIPSE_AREA_2D, the area of the ellipse.
//
{
  double pi = 3.141592653589793;
  double value;

  value = pi * r1 * r2;

  return value;
}
//****************************************************************************80

double ellipse_circumference_2d ( double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_CIRCUMFERENCE_2D returns the circumference of an ellipse in 2D.
//
//  Discussion:
//
//    There is no closed formula for the circumference of an ellipse.
//
//    Defining the eccentricity by
//
//      E = sqrt ( 1 - ( r2 / r1 )^2 )
//
//    where R1 and R2 are the major and minor axes, then
//
//      circumference
//        = 4 * R1 * E(K,2*PI)
//        = R1 * Integral ( 0 <= T <= 2*PI ) sqrt ( 1 - E^2 * sin^2 ( T ) ) dT
//
//    This integral can be approximated by the Gauss-Kummer formula.
//
//  Integration region:
//
//    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Harris, Horst Stocker,
//    Handbook of Mathematics and Computational Science,
//    Springer, 1998,
//    ISBN: 0-387-94746-9,
//    LC: QA40.S76.
//
//  Parameters:
//
//    Input, double R1, R2, the major and minor semi-axes.
//
//    Output, double ELLIPSE_CIRCUMFERENCE_2D, the
//    circumference of the ellipse.
//
{
  double e;
  int i;
  double pi = 3.141592653589793;
  double term;
  double value;

  if ( r1 == r2 )
  {
    value = 2.0 * pi * r1;
    return value;
  }
//
//  Compute the eccentricity of the ellipse.
//
  e = sqrt ( 1.0 - pow ( r8_min ( r1, r2 ) / r8_max ( r1, r2 ), 2 ) );

  value = 1.0;
  term = value;
  i = 0;

  for ( ; ; )
  {
    i = i + 1;
    term = term * ( 2 * i - 3 ) * ( 2 * i - 1 ) * e * e 
      / ( double ) ( 2 * 2 * i * i );

    if ( r8_abs ( term ) <= r8_epsilon ( ) * ( r8_abs ( value ) + 1.0 ) )
    {
      break;
    }
    value = value + term;
  }

  value = 2.0 * pi * r8_max ( r1, r2 ) * value;

  return value;
}
//****************************************************************************80

double ellipse_eccentricity_2d ( double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_ECCENTRICITY_2D returns the eccentricity of an ellipse in 2D.
//
//  Integration region:
//
//    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R1, R2, the major and minor semi-axes.
//
//    Output, double ELLIPSE_ECCENTRICITY_2D, the eccentricity 
//    of the ellipse.
//
{
  double major;
  double minor;
  double value;

  minor = r8_min ( r8_abs ( r1 ), r8_abs ( r2 ) );
  major = r8_max ( r8_abs ( r1 ), r8_abs ( r2 ) );

  if ( major == 0.0 )
  {
    value = - r8_huge ( );
    return value;
  }

  value = sqrt ( 1.0 - pow ( minor / major, 2 ) );

  return value;
}
//****************************************************************************80

double ellipsoid_volume_3d ( double r1, double r2, double r3 )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSOID_VOLUME_3D returns the volume of an ellipsoid in 3d.
//
//  Discussion:
//
//    This is not a general ellipsoid, but one for which each of the 
//    axes lies along a coordinate axis.
//
//  Integration region:
//
//      ( ( X - CENTER(1) ) / R1 )^2 
//    + ( ( Y - CENTER(2) ) / R2 )^2
//    + ( ( Z - CENTER(3) ) / R3 )^2 <= 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R1, R2, R3, the semi-axes of the ellipsoid.
//
//    Output, double ELLIPSOID_VOLUME_3D, the volume of the ellipsoid.
//
{
  double pi = 3.141592653589793;
  double value;

  value = ( 4.0 / 3.0 ) * pi * r1 * r2 * r3;

  return value;
}
//****************************************************************************80

void en_r2_01_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_01_1 implements the Stroud rule 1.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int k;
  double pi = 3.141592653589793E+00;
  double volume;

  volume = sqrt ( pow ( pi, n ) );

  r8vec_zero ( n * o, x );

  k = 0;
//
//  1 point.
//
  w[k] = volume;

  return;
}
//****************************************************************************80

int en_r2_01_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_01_1_SIZE sizes the Stroud rule 1.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_01_1_SIZE, the order.
//
{
  int o;

  o = 1;

  return o;
}
//****************************************************************************80

void en_r2_02_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_02_XIU implements the Xiu rule for region EN_R2.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double delta0;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = 2.0;
  delta0 = 0.0;
  c1 = 1.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  volume_1d = sqrt ( pi );
  volume = pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int en_r2_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_02_XIU_SIZE sizes the Xiu for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = N + 1;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_01_1_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

void en_r2_03_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_03_1 implements the Stroud rule 3.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2 * N.
//
//    The rule has precision P = 3.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  int i;
  int k;
  double pi = 3.141592653589793;
  double r;
  double volume;

  volume = sqrt ( pow ( pi, n ) );

  a = volume / ( double ) ( o );
  r = sqrt ( ( double ) ( n ) / 2.0 );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - r;
    w[k] = a;
    k = k + 1;
    x[i+k*n] = + r;
    w[k] = a;
  }

  return;
}
//****************************************************************************80

int en_r2_03_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_03_1_SIZE sizes the Stroud rule 3.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2 * N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_03_1_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void en_r2_03_2 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_03_2 implements the Stroud rule 3.2 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^N.
//
//    The rule has precision P = 3.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  int i;
  int i1;
  int k;
  bool more;
  double pi = 3.141592653589793E+00;
  double r;
  double volume;

  volume = sqrt ( pow ( pi, n ) );

  a = volume / ( double ) ( o );
  r = sqrt ( 0.5 );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  2^N points.
//
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - r;
  }
  w[k] = a;
  more = true;

  while ( more )
  {
    more = false;
    for ( i = n - 1; 0 <= i; i-- )
    {
      if ( x[i+k*n] < 0.0 )
      {
        k = k + 1;
        for ( i1 = 0; i1 < i; i1++ )
        {
          x[i1+k*n] = x[i1+(k-1)*n];
        }
        x[i+k*n]     = + r;
        for ( i1 = i + 1; i1 < n; i1++ )
        {
          x[i1+k*n] = - r;
        }
        w[k] = a;
        more = true;
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

int en_r2_03_2_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_03_2_SIZE sizes the Stroud rule 3.2 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_03_2_SIZE, the order.
//
{
  int o;

  o = i4_power ( 2, n );

  return o;
}
//****************************************************************************80

void en_r2_03_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_03_XIU implements the Xiu precision 3 rule for region EN_R2.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;

  volume = sqrt ( pow ( pi, n ) );

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( ( 2 * r - 1 ) * j ) * pi / ( double ) ( n );
      x[i+j*n] = cos ( arg );
      i = i + 1;
      x[i+j*n] = sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      if ( n == 1 )
      {
        x[i+j*n] = x[i+j*n] / sqrt ( 2.0 );
      }
      i = i + 1;
    }
  }
  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int en_r2_03_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_R2.
//
//  Discussion:
//
//    The rule has order 
//
//      O = 2 * N.
//
//    The rule has precision P = 3.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_XIU_SIZE, the order.
//
{
  int o;

  o = 2 * n;

  return o;
}
//****************************************************************************80

void en_r2_05_1 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_1 implements the Stroud rule 5.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = N^2 + N + 2.
//
//    The rule has precision P = 5.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
//    the OPTION variable to 1 or 2.
//
//    Versions of this rule are only available for N = 2 through 7.
//
//    There is a typographical error in the reference.
//    For the second version of the rule for N = 2, the line
//      gamma =    0.313300683022281E+00
//    should read
//      gamma =    0.312200683022281E+00
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    2 <= N <= 7.
//
//    Input, int OPTION, selects option 1 or 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double b;
  double c;
  double eta;
  double gamma;
  int i;
  int i1;
  int j;
  int k;
  double lambda;
  double mu;
  double pi = 3.141592653589793E+00;
  double volume;
  double xsi;

  if ( n < 2 || 7 < n )
  {
    cerr << "\n";
    cerr << "EN_R2_05_1 - Fatal error!\n";
    cerr << "  2 <= N <= 7 required.\n";
    exit ( 1 );
  }

  if ( option < 1 || 2 < option )
  {
    cerr << "\n";
    cerr << "EN_R2_05_1 - Fatal error!\n";
    cerr << "  1 <= OPTION <= 2 required.\n";
    exit ( 1 );
  }

  if ( option == 2 )
  {
    if ( n != 3 && n != 5 && n != 6 )
    {
      cerr << "\n";
      cerr << "EN_R2_05_1 - Fatal error!\n";
      cerr << "  OPTION = 2 requires N = 3, 5 or 6.\n";
      exit ( 1 );
    }
  }

  volume = sqrt ( pow ( pi, n ) );

  if ( n == 2 )
  {
    eta =      0.446103183094540E+00;
    lambda =   0.136602540378444E+01;
    xsi =    - 0.366025403784439E+00;
    mu =       0.198167882945871E+01;
    gamma =    0.000000000000000E+00;
    a =        0.328774019778636E+00 * volume;
    b =        0.833333333333333E-01 * volume;
    c =        0.455931355469736E-02 * volume;
  }
  else if ( n == 3 && option == 1 )
  {
    eta =      0.476731294622796E+00;
    lambda =   0.935429018879534E+00;
    xsi =    - 0.731237647787132E+00;
    mu =       0.433155309477649E+00;
    gamma =    0.266922328697744E+01;
    a =        0.242000000000000E+00 * volume;
    b =        0.810000000000000E-01 * volume;
    c =        0.500000000000000E-02 * volume;
  }
//
//  The value of gamma that follows corrects an error in the reference.
//
  else if ( n == 3 && option == 2 )
  {
    eta =      0.476731294622796E+00;
    lambda =   0.128679320334269E+01;
    xsi =    - 0.379873463323979E+00;
    mu =     - 0.192386729447751E+01;
    gamma =    0.312200683022281E+00;
    a =        0.242000000000000E+00 * volume;
    b =        0.810000000000000E-01 * volume;
    c =        0.500000000000000E-02 * volume;
  }
  else if ( n == 4 ) 
  {
    eta =      0.523945658287507E+00;
    lambda =   0.119433782552719E+01;
    xsi =    - 0.398112608509063E+00;
    mu =     - 0.318569372920112E+00;
    gamma =    0.185675837424096E+01;
    a =        0.155502116982037E+00 * volume;
    b =        0.777510584910183E-01 * volume;
    c =        0.558227484231506E-02 * volume;
  }
  else if ( n == 5 && option == 1 )
  {
    eta =      0.214972564378798E+01;
    lambda =   0.464252986016289E+01;
    xsi =    - 0.623201054093728E+00;
    mu =     - 0.447108700673434E+00;
    gamma =    0.812171426076311E+00;
    a =        0.487749259189752E-03 * volume;
    b =        0.487749259189752E-03 * volume;
    c =        0.497073504444862E-01 * volume;
  }
  else if ( n == 5 && option == 2 )
  {
    eta =      0.615369528365158E+00;
    lambda =   0.132894698387445E+01;
    xsi =    - 0.178394363877324E+00;
    mu =     - 0.745963266507289E+00;
    gamma =    0.135503972310817E+01;
    a =        0.726415024414905E-01 * volume;
    b =        0.726415024414905E-01 * volume;
    c =        0.641509853510569E-02 * volume;
  }
  else if ( n == 6 && option == 1 )
  {
    eta =      0.100000000000000E+01;
    lambda =   0.141421356237309E+01;
    xsi =      0.000000000000000E+00;
    mu =     - 0.100000000000000E+01;
    gamma =    0.100000000000000E+01;
    a =        0.781250000000000E-02 * volume;
    b =        0.625000000000000E-01 * volume;
    c =        0.781250000000000E-02 * volume;
  }
  else if ( n == 6 && option == 2 )
  {
    eta =      0.100000000000000E+01;
    lambda =   0.942809041582063E+00;
    xsi =    - 0.471404520791032E+00;
    mu =     - 0.166666666666667E+01;
    gamma =    0.333333333333333E+00;
    a =        0.781250000000000E-02 * volume;
    b =        0.625000000000000E-01 * volume;
    c =        0.781250000000000E-02 * volume;
  }
  else if ( n == 7 )
  {
    eta =      0.000000000000000E+00;
    lambda =   0.959724318748357E+00;
    xsi =    - 0.772326488820521E+00;
    mu =     - 0.141214270131942E+01;
    gamma =    0.319908106249452E+00;
    a =        0.111111111111111E+00 * volume;
    b =        0.138888888888889E-01 * volume;
    c =        0.138888888888889E-01 * volume;
  }

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  2 points.
//
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - eta;
  }
  w[k] = a;
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = + eta;
  }
  w[k] = a;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    for ( i1 = 0; i1 < n; i1++ )
    {
      x[i1+k*n] = - xsi;
    }
    x[i+k*n] = - lambda;
    w[k] = b;
    k = k + 1;
    for ( i1 = 0; i1 < n; i1++ )
    {
      x[i1+k*n] = + xsi;
    }
    x[i+k*n] = + lambda;
    w[k] = b;
  }
//
//  2 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      for ( i1 = 0; i1 < n; i1++ )
      {
        x[i1+k*n] = - gamma;
      }
      x[i+k*n] = - mu;
      x[j+k*n] = - mu;
      w[k] = c;
      k = k + 1;
      for ( i1 = 0; i1 < n; i1++ )
      {
        x[i1+k*n] = + gamma;
      }
      x[i+k*n] = + mu;
      x[j+k*n] = + mu;
      w[k] = c;
    }
  }
  return;
}
//****************************************************************************80

int en_r2_05_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_1_SIZE sizes the Stroud rule 5.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = N^2 + N + 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_05_1_SIZE, the order.
//
{
  int o;

  o = n * n + n + 2;

  return o;
}
//****************************************************************************80

void en_r2_05_2 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_2 implements the Stroud rule 5.2 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2 * N^2 + 1.
//
//    The rule has precision P = 5.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double b;
  double c;
  int i;
  int j;
  int k;
  double pi = 3.141592653589793E+00;
  double r;
  double s;
  double volume;

  volume = sqrt ( pow ( pi, n ) );

  a = 2.0E+00 * volume / ( double ) ( n + 2 );
  b = ( double ) ( 4 - n ) * volume / 2.0E+00 
    / ( double ) ( ( n + 2 ) * ( n + 2 ) );
  c = volume / ( double ) ( ( n + 2 ) * ( n + 2 ) );

  r = sqrt ( ( double ) ( n + 2 ) / 2.0E+00 );
  s = sqrt ( ( double ) ( n + 2 ) / 4.0E+00 );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = a;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - r;
    w[k] = b;
    k = k + 1;
    x[i+k*n] = + r;
    w[k] = b;
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - s;
      x[j+k*n] = - s;
      w[k] = c;
      k = k + 1;
      x[i+k*n] = - s;
      x[j+k*n] = + s;
      w[k] = c;
      k = k + 1;
      x[i+k*n] = + s;
      x[j+k*n] = - s;
      w[k] = c;
      k = k + 1;
      x[i+k*n] = + s;
      x[j+k*n] = + s;
      w[k] = c;
    }
  }
  return;
}
//****************************************************************************80

int en_r2_05_2_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_2_SIZE sizes the Stroud rule 5.2 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2 * N^2 + 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_01_1_SIZE, the order.
//
{
  int o;

  o = 2 * n * n + 1;

  return o;
}
//****************************************************************************80

void en_r2_05_3 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_3 implements the Stroud rule 5.3 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^N + 2 * N.
//
//    The rule has precision P = 5.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    The rule requires 3 <= N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    3 <= N.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double b;
  int i;
  int i1;
  int k;
  bool more;
  double pi = 3.141592653589793E+00;
  double r;
  double s;
  double volume;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "EN_R2_05_3 - Fatal error!\n";
    cerr << "  3 <= N is required.\n";
    exit ( 1 );
  }

  volume = sqrt ( pow ( pi, n ) );

  a = 4.0E+00 * volume / ( double ) ( ( n + 2 ) * ( n + 2 ) );
  b = ( double ) ( ( n - 2 ) * ( n - 2 ) ) * volume / ( double ) ( i4_power ( 2, n ) ) 
    / ( double ) ( ( n + 2 ) * ( n + 2 ) );
  r = sqrt ( ( double ) ( n + 2 ) / 4.0E+00 );
  s = sqrt ( ( double ) ( n + 2 ) / 2.0E+00 / ( double ) ( n - 2 ) );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - r;
    w[k] = a;
    k = k + 1;
    x[i+k*n] = + r;
    w[k] = a;
  }
//
//  2^N points.
//
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - s;
  }
  w[k] = b;
  more = true;
  while ( more )
  {
    more = false;
    for ( i = n - 1; 0 <= i; i-- )
    {
      if ( x[i+k*n] < 0.0E+00 )
      {
        k = k + 1;
        for ( i1 = 0; i1 < n; i1++ )
        {
          x[i1+k*n] = x[i1+(k-1)*n];
        }
        x[i+k*n] = + s;
        for ( i1 = i + 1; i1 < n; i1++ )
        {
          x[i1+k*n] = - s;
        }
        w[k] = b;
        more = true;
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

int en_r2_05_3_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_3_SIZE sizes the Stroud rule 5.3 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^N + 2 * N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_05_3_SIZE, the order.
//
{
  int o;

  o = i4_power ( 2, n ) + 2 * n;

  return o;
}
//****************************************************************************80

void en_r2_05_4 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_4 implements the Stroud rule 5.4 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^(N+1) - 1.
//
//    The rule has precision P = 5.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double b;
  int i;
  int i1;
  int j;
  int k;
  bool more;
  double pi = 3.141592653589793E+00;
  double r;
  double s;
  double volume;

  volume = sqrt ( pow ( pi, n ) );

  s = sqrt ( 0.5E+00 );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  2^N + 2^(N-1) + 2^(N-2) + ... + 1 = 2^(N+1)-1 points.
//  but do the last point separately.
//
  for ( i = 0; i < n; i++ )
  {
    r = sqrt ( ( double ) ( i + 3 ) / 2.0E+00 );
    b = pow ( 2.0E+00, i + 1 - n ) * volume / ( double ) ( i + 2 ) 
      / ( double ) ( i + 3 );

    k = k + 1;
    x[i+k*n] = - r;
    for ( i1 = i + 1; i1 < n; i1++ )
    {
      x[i1+k*n] = - s;
    }
    w[k] = b;
    more = true;
    while ( more )
    {
      more = false;
      for ( j = n - 1; 0 <= j; j-- )
      {
        if ( x[j+k*n] < 0.0E+00 )
        {
          k = k + 1;
          for ( i1 = 0; i1 < n; i1++ )
          {
            x[i1+k*n] = x[i1+(k-1)*n];
          }
          x[j+k*n] = r8_abs ( x[j+k*n] );
          for ( i1 = j + 1; i1 < n; i1++ )
          {
            x[i1+k*n] = - r8_abs ( x[i1+k*n] );
          }
          w[k] = b;
          more = true;
          break;
        }
      }
    }
  }
//
//  Last point.
//
  k = k + 1;
  w[k] = 2.0E+00 * volume / ( double ) ( n + 2 );

  return;
}
//****************************************************************************80

int en_r2_05_4_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_4_SIZE sizes the Stroud rule 5.4 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^(N+1) - 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_05_4_SIZE, the order.
//
{
  int o;

  o = i4_power ( 2, n + 1 ) - 1;

  return o;
}
//****************************************************************************80

void en_r2_05_5 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_5 implements the Stroud rule 5.5 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = N * 2^N + 1.
//
//    The rule has precision P = 5.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    There is a second version of this rule however it results in
//    complex abscissas, and so it has been disabled.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double b;
  int i;
  int i1;
  int j;
  int k;
  bool more;
  double n_r8;
  int option;
  double pi = 3.141592653589793E+00;
  double r;
  double s;
  double volume;

  volume = sqrt ( pow ( pi, n ) );

  n_r8 = ( double ) ( n );

  a = 2.0E+00 * volume / ( n_r8 + 2.0E+00 );
  b =           volume / ( n_r8 + 2.0E+00 ) / pow ( 2.0, n );

  option = 1;

  if ( option == 1 )
  {
    r = sqrt ( ( n_r8 + 2.0E+00 
      + ( n_r8 - 1.0E+00 ) * sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) 
      / 2.0E+00 / n_r8 );
    s = sqrt ( ( n_r8 + 2.0E+00 
      -                      sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) 
      / 2.0E+00 / n_r8 );
  }
  else if ( option == 2 )
  {
    r = sqrt ( ( n_r8 + 2.0E+00 
      - ( n_r8 - 1.0E+00 ) * sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) 
      / 2.0E+00 / n_r8 );
    s = sqrt ( ( n_r8 + 2.0E+00 
      +                      sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) 
      / 2.0E+00 / n_r8 );
  }

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = a;
//
//  N * 2^N points:
//  N choices for location of R, 2^N choices of sign pattern.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    for ( i1 = 0; i1 < n; i1++ )
    {
      x[i1+k*n] = - s;
    }
    x[i+k*n] = - r;
    w[k] = b;

    more = true;

    while ( more )
    {
      more = false;
      for ( j = n - 1; 0 <= j; j-- )
      {
        if ( x[j+k*n] < 0.0E+00 )
        {
          k = k + 1;
          for ( i1 = 0; i1 < n; i1++ )
          {
            x[i1+k*n] = x[i1+(k-1)*n];
          }
          x[j+k*n]     =   r8_abs ( x[j+k*n] );
          for ( i1 = j + 1; i1 < n; i1++ )
          {
            x[i1+k*n] = - r8_abs ( x[i1+k*n] );
          }
          w[k] = b;
          more = true;
          break;
        }
      }
    }
  }
  return;
}
//****************************************************************************80

int en_r2_05_5_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_5_SIZE sizes the Stroud rule 5.5 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = N * 2^N + 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_05_5_SIZE, the order.
//
{
  int o;

  o = n * i4_power ( 2, n ) + 1;

  return o;
}
//****************************************************************************80

void en_r2_05_6 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_6 implements the Stroud rule 5.6 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = ( N + 1 ) * 2^N.
//
//    The rule has precision P = 5.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    The rule requires 5 <= N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    5 <= N.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  int i;
  int i1;
  int j;
  int k;
  bool more;
  double n_r8;
  double pi = 3.141592653589793E+00;
  double r;
  double s;
  double t;
  double volume;

  if ( n < 5 )
  {
    cerr << "\n";
    cerr << "EN_R2_05_6 - Fatal error!\n";
    cerr << "  5 <= N is required.\n";
    exit ( 1 );
  }

  volume = sqrt ( pow ( pi, n ) );

  n_r8 = ( double ) ( n );

  a = volume / pow ( 2.0, n ) / ( n_r8 + 1.0E+00 );

  r = sqrt ( ( n_r8 - sqrt ( 2.0E+00 ) 
    + ( n_r8 - 1.0E+00 ) * sqrt ( 2.0E+00 * ( n_r8 + 1.0E+00 ) ) ) 
    / 2.0E+00 / n_r8 );
  s = sqrt ( ( n_r8 - sqrt ( 2.0E+00 ) 
    -                      sqrt ( 2.0E+00 * ( n_r8 + 1.0E+00 ) ) ) 
    / 2.0E+00 / n_r8 );
  t = sqrt ( ( 1.0E+00 + sqrt ( 2.0E+00 ) ) / 2.0E+00 );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  N * 2^N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    for ( i1 = 0; i1 < n; i1++ )
    {
      x[i1+k*n] = - s;
    }
    x[i+k*n] = - r;
    w[k] = a;

    more = true;

    while ( more )
    {
      more = false;
      for ( j = n - 1; 0 <= j; j-- )
      {
        if ( x[j+k*n] < 0.0E+00 )
        {
          k = k + 1;
          for ( i1 = 0; i1 < n; i1++ )
          {
            x[i1+k*n] = x[i1+(k-1)*n];
          }
          x[j+k*n] = r8_abs ( x[j+k*n] );
          for ( i1 = j + 1; i1 < n; i1++ )
          {
            x[i1+k*n] = - r8_abs ( x[i1+k*n] );
          }
          w[k] = a;
          more = true;
          break;
        }
      }
    }
  }
//
//  2^N points.
//
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - t;
  }
  w[k] = a;
  more = true;
  while ( more )
  {
    more = false;
    for ( j = n - 1; 0 <= j; j-- )
    {
      if ( x[j+k*n] < 0.0E+00 )
      {
        k = k + 1;
        for ( i1 = 0; i1 < n; i1++ )
        {
          x[i1+k*n] = x[i1+(k-1)*n];
        }
        x[j+k*n] = r8_abs ( x[j+k*n] );
        for ( i1 = j + 1; i1 < n; i1++ )
        {
          x[i1+k*n] = - r8_abs ( x[i1+k*n] );
        }
        w[k] = a;
        more = true;
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

int en_r2_05_6_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_05_6_SIZE sizes the Stroud rule 5.6 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = ( N + 1 ) * 2^N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_05_6_SIZE, the order.
//
{
  int o;

  o = ( n + 1 ) * i4_power ( 2, n );

  return o;
}
//****************************************************************************80

void en_r2_07_1 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_07_1 implements the Stroud rule 7.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^N + 2 * N^2 + 1.
//
//    The rule has precision P = 7.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    There are two versions of the rule, chosen by setting the
//    OPTION variable to 1 or 2.  
//
//    Option 1 is only valid for N = 3, 4, 6 or 7.
//    Option 2 is only valid for N = 3 or 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    N = 3, 4, 6 or 7.
//
//    Input, int OPTION, chooses rule option 1 or 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a;
  double b;
  double c;
  double d;
  int i;
  int i1;
  int j;
  int k;
  bool more;
  double n_r8;
  double pi = 3.141592653589793E+00;
  double r;
  double s;
  double t;
  double volume;

  if ( option < 1 || 2 < option )
  {
    cerr << "\n";
    cerr << "EN_R2_07_1 - Fatal error!\n";
    cerr << "  1 <= OPTION <= 2 required.\n";
    exit ( 1 );
  }

  if ( option == 1 )
  {
    if ( n != 3 && n != 4 && n != 6 && n != 7 )
    {
      cerr << "\n";
      cerr << "EN_R2_07_1 - Fatal error!\n";
      cerr << "  OPTION 1 requires N =  3, 4, 6 or 7.\n";
      exit ( 1 );
    }
  }

  if ( option == 2 )
  {
    if ( n != 3 && n != 4 )
    {
      cerr << "\n";
      cerr << "EN_R2_07_1 - Fatal error!\n";
      cerr << "  OPTION 2 requires N =  3 or 4.\n";
      exit ( 1 );
    }
  }

  volume = sqrt ( pow ( pi, n ) );

  n_r8 = ( double ) ( n );

  if ( option == 1 )
  {
    r = sqrt ( ( 3.0E+00 * ( 8.0E+00 - n_r8 ) - ( n_r8 - 2.0E+00 ) 
      * sqrt ( 3.0E+00 * ( 8.0E+00 - n_r8 ) ) ) / 2.0E+00 / ( 5.0E+00 - n_r8 ) );
    s = sqrt ( ( 3.0E+00 *             n_r8   -          2.0E+00   
      * sqrt ( 3.0E+00 * ( 8.0E+00 - n_r8 ) ) ) / 2.0E+00 
      / ( 3.0E+00 * n_r8 - 8.0E+00 ) );
    t = sqrt ( ( 6.0E+00 + sqrt ( 3.0E+00 * ( 8.0E+00 - n_r8 ) ) ) / 2.0E+00 );
  }
  else if ( option == 2 )
  {
    r = sqrt ( ( 3.0E+00 * ( 8.0E+00 - n_r8 ) + ( n_r8 - 2.0E+00 ) 
      * sqrt ( 3.0E+00 * ( 8.0E+00 - n_r8 ) ) ) / 2.0E+00 / ( 5.0E+00 - n_r8 ) );
    s = sqrt ( ( 3.0E+00 *             n_r8   +          2.0E+00   
      * sqrt ( 3.0E+00 * ( 8.0E+00 - n_r8 ) ) ) / 2.0E+00 
      / ( 3.0E+00 * n_r8 - 8.0E+00 ) );
    t = sqrt ( ( 6.0E+00 - sqrt ( 3.0E+00 * ( 8.0E+00 - n_r8 ) ) ) / 2.0E+00 );
  }

  b = ( 8.0E+00 - n_r8 ) * volume / 8.0E+00 / pow ( r, 6 );
  c = volume / pow ( 2.0E+00, n + 3 ) / pow ( s, 6 );
  d = volume / 16.0E+00 / pow ( t, 6 );
  a = volume - 2.0E+00 * n_r8 * b - pow ( 2.0E+00, n ) * c - 2.0E+00 * n_r8 
    * ( n_r8 - 1.0E+00 ) * d;

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = a;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - r;
    w[k] = b;
    k = k + 1;
    x[i+k*n] = + r;
    w[k] = b;
  }
//
//  2^N points.
//
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - s;
  }
  w[k] = c;
  more = true;
  while ( more )
  {
    more = false;
    for ( i = n - 1; 0 <= i; i-- )
    {
      if ( x[i+k*n] < 0.0E+00 )
      {
        k = k + 1;
        for ( i1 = 0; i1 < n; i1++ )
        {
          x[i1+k*n] = x[i1+(k-1)*n];
        }
        x[i+k*n] = r8_abs ( x[i+k*n] );
        for ( i1 = i + 1; i1 < n; i1++ )
        {
          x[i1+k*n] = - r8_abs ( x[i1+k*n] );
        }
        w[k] = c;
        more = true;
        break;
      }
    }
  }
//
//  2 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - t;
      x[j+k*n] = - t;
      w[k] = d;
      k = k + 1;
      x[i+k*n] = - t;
      x[j+k*n] = + t;
      w[k] = d;
      k = k + 1;
      x[i+k*n] = + t;
      x[j+k*n] = - t;
      w[k] = d;
      k = k + 1;
      x[i+k*n] = + t;
      x[j+k*n] = + t;
      w[k] = d;
    }
  }
  return;
}
//****************************************************************************80

int en_r2_07_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_07_1_SIZE sizes the Stroud rule 7.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^N + 2 * N^2 + 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_07_1_SIZE, the order.
//
{
  int o;

  o = i4_power ( 2, n ) + 2 * n * n + 1;

  return o;
}
//****************************************************************************80

void en_r2_07_2 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_07_2 implements the Stroud rule 7.2 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^(N+1) + 4 * N^2.
//
//    The rule has precision P = 7.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    The rule requires 3 <= N.
//
//    The reference has a typographical error in the description of this rule.
//    The formula:
//
//      (t,t,t,...,t)FS
//
//    should read
//
//      (t,t,0,...,0)FS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    3 <= N.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double a1;
  double a2;
  double b;
  double c;
  double d;
  int i;
  int i1;
  int j;
  int k;
  bool more;
  double n_r8;
  double pi = 3.141592653589793E+00;
  double r;
  double rho1;
  double rho2;
  double s;
  double t;
  double volume;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "EN_R2_07_2 - Fatal error!\n";
    cerr << "  3 <= N is required.\n";
    exit ( 1 );
  }

  volume = sqrt ( pow ( pi, n ) );

  n_r8 = ( double ) ( n );

  rho1 = sqrt ( ( n_r8 + 2.0E+00 - sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) 
    / 2.0E+00 );
  rho2 = sqrt ( ( n_r8 + 2.0E+00 + sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) 
    / 2.0E+00 );
  a1 = ( n_r8 + 2.0E+00 + sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) / 2.0E+00 
    / ( n_r8 + 2.0E+00 );
  a2 = ( n_r8 + 2.0E+00 - sqrt ( 2.0E+00 * ( n_r8 + 2.0E+00 ) ) ) / 2.0E+00 
    / ( n_r8 + 2.0E+00 );

  r = 1.0E+00;
  s = sqrt ( 1.0E+00 / n_r8 );
  t = sqrt ( 0.5E+00 );
  b = ( 8.0E+00 - n_r8 ) * volume / n_r8 / ( n_r8 + 2.0E+00 ) / ( n_r8 + 4.0E+00 );
  c = pow ( n_r8, 3 ) * volume / pow ( 2.0E+00, n ) / n_r8 / ( n_r8 + 2.0E+00 ) 
    / ( n_r8 + 4.0E+00 );
  d = 4.0E+00 * volume / n_r8 / ( n_r8 + 2.0E+00 ) / ( n_r8 + 4.0E+00 );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  2 * 2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - rho1 * r;
    w[k] = a1 * b;
    k = k + 1;
    x[i+k*n] = - rho2 * r;
    w[k] = a2 * b;
    k = k + 1;
    x[i+k*n] = + rho1 * r;
    w[k] = a1 * b;
    k = k + 1;
    x[i+k*n] = + rho2 * r;
    w[k] = a2 * b;
  }
//
//  2 * 2^N points.
//
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - rho1 * s;
  }
  w[k] = a1 * c;
  k = k + 1;
  for ( i1 = 0; i1 < n; i1++ )
  {
    x[i1+k*n] = - rho2 * s;
  }
  w[k] = a2 * c;
  more = true;
  while ( more )
  {
    more = false;
    for ( i = n - 1; 0 <= i; i-- )
    {
      if ( x[i+k*n] < 0.0E+00 )
      {
        k = k + 1;
        for ( i1 = 0; i1 < n; i1++ )
        {
          x[i1+k*n] = x[i1+(k-2)*n];;
        }
        x[i+k*n] = r8_abs ( x[i+k*n] );
        for ( i1 = i + 1; i1 < n; i1++ )
        {
          x[i1+k*n] = - r8_abs ( x[i1+k*n] );
        }
        w[k] = a1 * c;
        k = k + 1;
        for ( i1 = 0; i1 < n; i1++ )
        {
          x[i1+k*n] = x[i1+(k-2)*n];;
        }
        x[i+k*n] = r8_abs ( x[i+k*n] );
        for ( i1 = i + 1; i1 < n; i1++ )
        {
          x[i1+k*n] = - r8_abs ( x[i1+k*n] );
        }
        w[k] = a2 * c;
        more = true;
        break;
      }
    }
  }
//
//  2 * 4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - rho1 * t;
      x[j+k*n] = - rho1 * t;
      w[k] = a1 * d;
      k = k + 1;
      x[i+k*n] = - rho1 * t;
      x[j+k*n] = + rho1 * t;
      w[k] = a1 * d;
      k = k + 1;
      x[i+k*n] = + rho1 * t;
      x[j+k*n] = - rho1 * t;
      w[k] = a1 * d;
      k = k + 1;
      x[i+k*n] = + rho1 * t;
      x[j+k*n] = + rho1 * t;
      w[k] = a1 * d;
      k = k + 1;
      x[i+k*n] = - rho2 * t;
      x[j+k*n] = - rho2 * t;
      w[k] = a2 * d;
      k = k + 1;
      x[i+k*n] = - rho2 * t;
      x[j+k*n] = + rho2 * t;
      w[k] = a2 * d;
      k = k + 1;
      x[i+k*n] = + rho2 * t;
      x[j+k*n] = - rho2 * t;
      w[k] = a2 * d;
      k = k + 1;
      x[i+k*n] = + rho2 * t;
      x[j+k*n] = + rho2 * t;
      w[k] = a2 * d;
    }
  }
  return;
}
//****************************************************************************80

int en_r2_07_2_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_07_2_SIZE sizes the Stroud rule 7.2 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = 2^(N+1) + 4 * N^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_07_2_SIZE, the order.
//
{
  int o;

  o = i4_power ( 2, n + 1 ) + 4 * n * n;

  return o;
}
//****************************************************************************80

void en_r2_07_3 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_07_3 implements the Stroud rule 7.3 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
//
//    The rule has precision P = 7.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    There are two versions of each rule, chosen by setting the
//    OPTION variable to 1 or 2.
//
//    The rule as tabulated by Stenger is available for N = 2 through 20.
//    This function accepts N = 3 through 6.
//
//     N    O
//    __  ___
//     3   45
//     4   97
//     5  181
//     6  305
//
//    The reference has a typographical error for N = 5, OPTION 1, B4:
//
//      -(1)0.736330882774831
//
//    should read
//
//      (-1)0.736330882774831
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    3 <= N <= 6.
//
//    Input, int OPTION, chooses rule option 1 or 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double b0;
  double b1;
  double b2;
  double b3;
  double b4;
  double b5;
  int i;
  int j;
  int k;
  int l;
  bool more;
  double pi = 3.141592653589793E+00;
  double u;
  double v;
  double volume;

  if ( n < 3 || 6 < n )
  {
    cerr << "\n";
    cerr << "EN_R2_07_3 - Fatal error!\n";
    cerr << "  3 <= N <= 6 required.\n";
    exit ( 1 );
  }

  if ( option < 1 || 2 < option )
  {
    cerr << "\n";
    cerr << "EN_R2_07_3 - Fatal error!\n";
    cerr << "  1 <= OPTION <= 2 required.\n";
    exit ( 1 );
  }

  volume = sqrt ( pow ( pi, n ) );

  if ( n == 3 && option == 1 )
  {
    u =    0.524647623275290E+00;
    v =    0.165068012388578E+01;
    b0 = - 0.166705761599566E+02;
    b1 =   0.100296981655678E+02;
    b2 =   0.161699246687754E+00;
    b3 = - 0.604719151221535E+01;
    b4 =   0.234381399489666E-01;
    b5 =   0.417194501880647E+01;
  }
  else if ( n == 3 && option == 2 )
  {
    u =    0.165068012388578E+01;
    v =    0.524647623275290E+00;
    b0 =   0.166705761599566E+02;
    b1 =   0.178903161957074E+00;
    b2 = - 0.665808190965810E+01;
    b3 =   0.148361823143070E-01;
    b4 =   0.229669852539758E+01;
    b5 =   0.430097881732984E-02;
  }
  else if ( n == 4 && option == 1 )
  {
    u  =   0.524647623275290E+00;
    v  =   0.165068012388578E+01;
    b0 = - 0.167539329651562E+03;
    b1 =   0.687922329603575E+02;
    b2 =   0.203518409659014E+00;
    b3 = - 0.255075279116885E+02;
    b4 =   0.415430214106084E-01;
    b5 =   0.739458001434961E+01;
  }
  else if ( n == 4 && option == 2 )
  {
    u =    0.165068012388578E+01;
    v =    0.524647623275290E+00;
    b0 =   0.688432856406677E+02;
    b1 =   0.294997847268286E+00;
    b2 = - 0.199427272118378E+02;
    b3 =   0.110498755408511E-01;
    b4 =   0.407079214570997E+01;
    b5 =   0.762328646743931E-02;
  }
  else if ( n == 5 && option == 1 )
  {
    u  =   0.524647623275290E+00;
    v  =   0.165068012388578E+01;
    b0 = - 0.826940846964452E+03;
    b1 =   0.264779097660331E+03;
    b2 =   0.213460812375320E+00;
    b3 = - 0.714240197186780E+02;
    b4 =   0.736330882774831E-01;
    b5 =   0.131065518222629E+02;
  }
  else if ( n == 5 && option == 2 )
  {
    u =    0.165068012388578E+01;
    v =    0.524647623275290E+00;
    b0 =   0.220502344940121E+03;
    b1 =   0.537746975313769E+00;
    b2 = - 0.497781460739792E+02;
    b3 = - 0.743845245712926E-02;
    b4 =   0.721529121489956E+01;
    b5 =   0.135119234557687E-01;
  }
  else if ( n == 6 && option == 1 )
  {
    u  =   0.524647623275290E+00;
    v  =   0.165068012388578E+01;
    b0 = - 0.309679578630802E+04;
    b1 =   0.815423321880237E+03;
    b2 =   0.117326937169073E+00;
    b3 = - 0.173057295296448E+03;
    b4 =   0.130511250871491E+00;
    b5 =   0.232307582494626E+02;
  }
  else if ( n == 6 && option == 2 )
  {
    u =    0.165068012388578E+01;
    v =    0.524647623275290E+00;
    b0 =   0.616293651884027E+03;
    b1 =   0.107529736766179E+01;
    b2 = - 0.113807008098269E+03;
    b3 = - 0.610828352270520E-01;
    b4 =   0.127887706992535E+02;
    b5 =   0.239492607623178E-01;
  }

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = b0;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - u;
    w[k] = b1;
    k = k + 1;
    x[i+k*n] = + u;
    w[k] = b1;
  }
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - v;
    w[k] = b2;
    k = k + 1;
    x[i+k*n] = + v;
    w[k] = b2;
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = - u;
      w[k] = b3;
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = + u;
      w[k] = b3;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = - u;
      w[k] = b3;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = + u;
      w[k] = b3;
    }
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = - v;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = + v;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = - v;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = + v;
      w[k] = b4;
    }
  }
//
//  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
//
  for ( i = 0; i < n - 2; i++ )
  {
    for ( j = i + 1; j < n - 1; j++ )
    {
      for ( l = j + 1; l < n; l++ )
      {
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b5;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b5;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b5;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b5;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b5;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b5;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b5;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b5;
      }
    }
  }
  return;
}
//****************************************************************************80

int en_r2_07_3_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_07_3_SIZE sizes the Stroud rule 7.3 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_07_3_SIZE, the order.
//
{
  int o;

  o = ( 4 * n * n * n + 8 * n + 3 ) / 3;

  return o;
}
//****************************************************************************80

void en_r2_09_1 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_09_1 implements the Stroud rule 9.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
//
//    The rule has precision P = 9.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    There are two versions of each rule, chosen by setting the 
//    OPTION variable to 1 or 2.
//
//    The rule as tabulated by Stenger is available for N = 2 through 20.
//    This function accepts N = 3 through 6.
//
//     N    O
//    __  ___
//     3   77
//     4  193
//     5  421
//     6  825
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    3 <= N <= 6.
//
//    Input, int OPTION, chooses rule option 1 or 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double b0;
  double b1;
  double b2;
  double b3;
  double b4;
  double b5;
  double b6;
  double b7;
  double b8;
  int i;
  int j;
  int k;
  int l;
  int m;
  bool more;
  double pi = 3.141592653589793E+00;
  double u;
  double v;
  double volume;

  if ( n < 3 || 6 < n )
  {
    cerr << "\n";
    cerr << "EN_R2_09_1 - Fatal error!\n";
    cerr << "  3 <= N <= 6 required.\n";
    exit ( 1 );
  }

  if ( option < 1 || 2 < option )
  {
    cerr << "\n";
    cerr << "EN_R2_09_1 - Fatal error!\n";
    cerr << "  1 <= OPTION <= 2 required.\n";
    exit ( 1 );
  }

  volume = sqrt ( pow ( pi, n ) );

  if ( n == 3 )
  {
    u =    0.202018287045609E+01;
    v =    0.958572464613819E+00;
    b0 =   0.676448734429924E+00;
    b1 =   0.511989106291551E-02;
    b2 =   0.448595723493744E+00;
    b3 =   0.235223454595606E-03;
    b4 =   0.915390713080005E-01;
    b5 =   0.139208199920793E-01;
    b6 =   0.235223454595606E-03;
    b7 =   0.915390713080008E-01;
    b8 =   0.000000000000000E+00;
  }
  else if ( n == 4 && option == 1 )
  {
    u =    0.202018287045609E+01;
    v =    0.958572464613819E+00;
    b0 = - 0.860452945007048E+00;
    b1 = - 0.405511998533795E-01;
    b2 =   0.107026475449715E+01;
    b3 =   0.138974239307092E-03;
    b4 = - 0.162248779448181E+00;
    b5 =   0.246740110027234E-01;
    b6 =   0.138974239307094E-03;
    b7 =   0.162248779448181E+00;
    b8 =   0.138974239307094E-03;
  }
  else if ( n == 4 && option == 2 )
  {
    u =    0.958572464613819E+00;
    v =    0.202018287045609E+01;
    b0 =   0.265029088766810E-02;
    b1 =   0.637601342635332E+00;
    b2 = - 0.394394059389228E-01;
    b3 =   0.540829264827264E-01;
    b4 = - 0.416922717921281E-03;
    b5 =   0.246740110027234E-01;
    b6 =   0.540829264827270E-01;
    b7 =   0.416922717921281E-03;
    b8 =   0.540829264827269E-01;
  }
  else if ( n == 5 && option == 1 )
  {
    u =    0.202018287045609E+01;
    v =    0.958572464613819E+00;
    b0 = - 0.827347006200826E+01;
    b1 = - 0.160820174530905E+00;
    b2 =   0.353499863758467E+01;
    b3 =   0.738976276909564E-03;
    b4 = - 0.862735421812943E+00;
    b5 =   0.437335458190621E-01;
    b6 = - 0.246325425636523E-03;
    b7 =   0.287578473937648E+00;
    b8 =   0.246325425636523E-03;
  }
  else if ( n == 5 && option == 2 )
  {
    u =    0.958572464613819E+00;
    v =    0.202018287045609E+01;
    b0 = - 0.624416791055272E+00;
    b1 =   0.467494915583104E+00;
    b2 = - 0.152937760910536E+00;
    b3 =   0.287578473937646E+00;
    b4 = - 0.221692883072871E-02;
    b5 =   0.437335458190621E-01;
    b6 = - 0.958594913125490E-01;
    b7 =   0.738976276909568E-03;
    b8 =   0.958594913125492E-01;
  }
  else if ( n == 6 && option == 1 )
  {
    u =    0.202018287045609E+01;
    v =    0.958572464613819E+00;
    b0 = - 0.361840434143098E+02;
    b1 = - 0.447936529138517E+00;
    b2 =   0.112077863004144E+02;
    b3 =   0.392940404320855E-02;
    b4 = - 0.254859786784158E+01;
    b5 =   0.775156917007496E-01;
    b6 = - 0.130980134773619E-02;
    b7 =   0.509719573568315E+00;
    b8 =   0.436600449245395E-03;
  }
  else if ( n == 6 && option == 2 )
  {
    u =    0.958572464613819E+00;
    v =    0.202018287045609E+01;
    b0 =   0.448873836333650E+01;
    b1 = - 0.238473566140736E+01;
    b2 = - 0.413008493198885E+00;
    b3 =   0.152915872070494E+01;
    b4 = - 0.654900673868093E-02;
    b5 =   0.775156917007496E-01;
    b6 = - 0.509719573568314E+00;
    b7 =   0.130980134773618E-02;
    b8 =   0.169906524522772E+00;
  }

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = b0;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - u;
    w[k] = b1;
    k = k + 1;
    x[i+k*n] = + u;
    w[k] = b1;
  }
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - v;
    w[k] = b2;
    k = k + 1;
    x[i+k*n] = + v;
    w[k] = b2;
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = - u;
      w[k] = b3;
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = + u;
      w[k] = b3;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = - u;
      w[k] = b3;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = + u;
      w[k] = b3;
    }
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = - v;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = + v;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = - v;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = + v;
      w[k] = b4;
    }
  }
//
//  4 * ( N * ( N - 1 ) ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = - v;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = + v;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = - v;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = + v;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = - u;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = + u;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = - u;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = + u;
      w[k] = b5;
    }
  }
//
//  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
//
  for ( i = 0; i < n - 2; i++ )
  {
    for ( j = i + 1; j < n - 1; j++ )
    {
      for ( l = j + 1; l < n; l++ )
      {
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b6;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b6;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b6;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b6;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b6;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b6;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b6;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b6;
      }
    }
  }
//
//  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
//
  for ( i = 0; i < n - 2; i++ )
  {
    for ( j = i + 1; j < n - 1; j++ )
    {
      for ( l = j + 1; l < n; l++ )
      {
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = - v;
        x[l+k*n] = - v;
        w[k] = b7;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = - v;
        x[l+k*n] = + v;
        w[k] = b7;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = + v;
        x[l+k*n] = - v;
        w[k] = b7;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = + v;
        x[l+k*n] = + v;
        w[k] = b7;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = - v;
        x[l+k*n] = - v;
        w[k] = b7;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = - v;
        x[l+k*n] = + v;
        w[k] = b7;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = + v;
        x[l+k*n] = - v;
        w[k] = b7;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = + v;
        x[l+k*n] = + v;
        w[k] = b7;
      }
    }
  }
//
//  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
//
  for ( i = 0; i < n - 3; i++ )
  {
    for ( j = i + 1; j < n - 2; j++ )
    {
      for ( l = j + 1; l < n - 1; l++ )
      {
        for ( m = l + 1; m < n; m++ )
        {
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b8;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b8;
        }
      }
    }
  }
  return;
}
//****************************************************************************80

int en_r2_09_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_09_1_SIZE sizes the Stroud rule 9.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_09_1_SIZE, the order.
//
{
  int o;

  o = (  2 * i4_power ( n, 4 )
      -  4 * i4_power ( n, 3 )
      + 22 * i4_power ( n, 2 )
      -  8 *            n 
      +  3 ) / 3;

  return o;
}
//****************************************************************************80

void en_r2_11_1 ( int n, int option, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_11_1 implements the Stroud rule 11.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order 
//
//      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
//
//    The rule has precision P = 11.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    There are two versions of each rule, chosen by setting the
//    OPTION variable to 1 or 2.
//
//    The rule as tabulated by Stenger is available for N = 2 through 20.
//    This function accepts N = 3 through 5.
//
//     N    O
//    __  ___
//     3  151
//     4  417
//     5  983
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    3 <= N <= 5.
//
//    Input, int OPTION, chooses rule option 1 or 2.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double b0;
  double b1;
  double b2;
  double b3;
  double b4;
  double b5;
  double b6;
  double b7;
  double b8;
  double b9;
  double b10;
  double b11;
  double b12;
  double b13;
  double b14;
  double b15;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int j;
  int k;
  int l;
  int m;
  bool more;
  double pi = 3.141592653589793E+00;
  double u;
  double v;
  double volume;
  double w2;

  if ( n < 3 || 5 < n )
  {
    cerr << "\n";
    cerr << "EN_R2_11_1 - Fatal error!\n";
    cerr << "  3 <= N <= 5 required.\n";
    exit ( 1 );
  }

  if ( option < 1 || 2 < option )
  {
    cerr << "\n";
    cerr << "EN_R2_11_1 - Fatal error!\n";
    cerr << "  1 <= OPTION <= 2 required.\n";
    exit ( 1 );
  }

  volume = sqrt ( pow ( pi, n ) );

  if ( n == 3 && option == 1 )
  {
    u =     0.235060497367449E+01;
    v =     0.436077411927617E+00;
    w2 =    0.133584907401370E+01;
    b0 =  - 0.881591029957858E+01;
    b1 =  - 0.751996143360650E-01;
    b2 =    0.621743189471515E+01;
    b3 =    0.241426451456494E+00;
    b4 =  - 0.120709739276065E-02;
    b5 =  - 0.427751221210138E+01;
    b6 =    0.550169924840163E-01;
    b7 =    0.237084999634707E-01;
    b8 =  - 0.169791992887741E-02;
    b9 =  - 0.252266276123350E-04;
    b10 =   0.326777873717691E+01;
    b11 =   0.968469949206802E-02;
    b12 =   0.789754514877422E-03;
    b13 =   0.000000000000000E+00;
    b14 =   0.000000000000000E+00;
    b15 =   0.000000000000000E+00;
  }
  else if ( n == 3 && option == 2 )
  {
    u =     0.235060497367449E+01;
    v =     0.133584907401370E+01;
    w2 =    0.436077411927617E+00;
    b0 =  - 0.141214037032900E+02;
    b1 =  - 0.803730274707282E-01;
    b2 =    0.235546545595906E+00;
    b3 =    0.888123191556611E+01;
    b4 =    0.142467131155533E-03;
    b5 =    0.582993124006494E-01;
    b6 =  - 0.561099173155661E+01;
    b7 =  - 0.204028691521686E-02;
    b8 =    0.252880089932256E-01;
    b9 =  - 0.814378678627283E-04;
    b10 =   0.804353953375146E-02;
    b11 =   0.393451849690453E+01;
    b12 =   0.171183493169724E-03;
    b13 =   0.000000000000000E+00;
    b14 =   0.000000000000000E+00;
    b15 =   0.000000000000000E+00;
  }
  else if ( n == 4 && option == 1 )
  {
    u =     0.235060497367449E+01;
    v =     0.436077411927617E+00;
    w2 =    0.133584907401370E+01;
    b0 =    0.241502736147339E+03;
    b1 =  - 0.196095938531478E+00;
    b2 =  - 0.128675737999280E+03;
    b3 =    0.307568784278696E+00;
    b4 =  - 0.480908422319460E-02;
    b5 =    0.698087019367085E+02;
    b6 =    0.631837143743771E-01;
    b7 =    0.392226151971179E-01;
    b8 =  - 0.300948471646799E-02;
    b9 =  - 0.650235306755170E-04;
    b10 = - 0.386951974646715E+02;
    b11 =   0.171656829095787E-01;
    b12 =   0.139980343116450E-02;
    b13 =   0.101552487093372E-04;
    b14 =   0.222435922356439E+02;
    b15 =   0.000000000000000E+00;
  }
  else if ( n == 4 && option == 2 )
  {
    u =     0.235060497367449E+01;
    v =     0.133584907401370E+01;
    w2 =    0.436077411927617E+00;
    b0 =  - 0.151944464736584E+03;
    b1 =  - 0.223498438689039E+00;
    b2 =    0.243574919068010E+00;
    b3 =    0.634373877008693E+02;
    b4 =  - 0.782065187814018E-04;
    b5 =    0.911833754536616E-01;
    b6 =  - 0.238927288245914E+02;
    b7 =  - 0.422314408318853E-02;
    b8 =    0.448218289217760E-01;
    b9 =  - 0.138053374667391E-03;
    b10 =   0.607473265800655E-02;
    b11 =   0.697375246129742E+01;
    b12 =   0.303414841680135E-03;
    b13 = - 0.314574391771792E-05;
    b14 =   0.409103498175100E-02;
    b15 =   0.000000000000000E+00;
  }
  else if ( n == 5 && option == 1 )
  {
    u =     0.235060497367449E+01;
    v =     0.436077411927617E+00;
    w2 =    0.133584907401370E+01;
    b0 =    0.255885269311763E+04;
    b1 =  - 0.439598677491526E+00;
    b2 =  - 0.106541406144610E+04;
    b3 =    0.453540909054264E+00;
    b4 =  - 0.132100905623778E-01;
    b5 =    0.418606568954203E+03;
    b6 =    0.511394563043680E-01;
    b7 =    0.645581013845604E-01;
    b8 =  - 0.533417277494500E-02;
    b9 =  - 0.137981626254496E-03;
    b10 = - 0.147436933189884E+03;
    b11 =   0.304253807765057E-01;
    b12 =   0.248108698207828E-02;
    b13 =   0.113652094546015E-04;
    b14 =   0.394257407160391E+02;
    b15 =   0.331725011358320E-05;
  }
  else if ( n == 5 && option == 2 )
  {
    u =     0.235060497367449E+01;
    v =     0.133584907401370E+01;
    w2 =    0.436077411927617E+00;
    b0 =  - 0.761305347548192E+03;
    b1 =  - 0.536360805019297E+00;
    b2 =    0.110669832078736E+00;
    b3 =    0.246421088923968E+03;
    b4 =  - 0.773649327968607E-03;
    b5 =    0.169088641205970E+00;
    b6 =  - 0.670700680243651E+02;
    b7 =  - 0.856090560229205E-02;
    b8 =    0.794446232770302E-01;
    b9 =  - 0.220272863263544E-03;
    b10 = - 0.373515812228225E-02;
    b11 =   0.123606544052884E+02;
    b12 =   0.537788804557843E-03;
    b13 = - 0.122101861480881E-04;
    b14 =   0.725117070759373E-02;
    b15 =   0.331725011358320E-05;
  }

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  w[k] = b0;
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - u;
    w[k] = b1;
    k = k + 1;
    x[i+k*n] = + u;
    w[k] = b1;
  }
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - v;
    w[k] = b2;
    k = k + 1;
    x[i+k*n] = + v;
    w[k] = b2;
  }
//
//  2 * N points.
//
  for ( i = 0; i < n; i++ )
  {
    k = k + 1;
    x[i+k*n] = - w2;
    w[k] = b3;
    k = k + 1;
    x[i+k*n] = + w2;
    w[k] = b3;
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = - u;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = + u;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = - u;
      w[k] = b4;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = + u;
      w[k] = b4;
    }
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = - v;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = + v;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = - v;
      w[k] = b5;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = + v;
      w[k] = b5;
    }
  }
//
//  4 * ( N * ( N - 1 ) / 2 ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - w2;
      x[j+k*n] = - w2;
      w[k] = b6;
      k = k + 1;
      x[i+k*n] = - w2;
      x[j+k*n] = + w2;
      w[k] = b6;
      k = k + 1;
      x[i+k*n] = + w2;
      x[j+k*n] = - w2;
      w[k] = b6;
      k = k + 1;
      x[i+k*n] = + w2;
      x[j+k*n] = + w2;
      w[k] = b6;
    }
  }
//
//  4 * ( N * ( N - 1 ) ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = - v;
      w[k] = b7;
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = + v;
      w[k] = b7;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = - v;
      w[k] = b7;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = + v;
      w[k] = b7;
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = - u;
      w[k] = b7;
      k = k + 1;
      x[i+k*n] = - v;
      x[j+k*n] = + u;
      w[k] = b7;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = - u;
      w[k] = b7;
      k = k + 1;
      x[i+k*n] = + v;
      x[j+k*n] = + u;
      w[k] = b7;
    }
  }
//
//  4 * ( N * ( N - 1 ) ) points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = - w2;
      w[k] = b8;
      k = k + 1;
      x[i+k*n] = - u;
      x[j+k*n] = + w2;
      w[k] = b8;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = - w2;
      w[k] = b8;
      k = k + 1;
      x[i+k*n] = + u;
      x[j+k*n] = + w2;
      w[k] = b8;
      k = k + 1;
      x[i+k*n] = - w2;
      x[j+k*n] = - u;
      w[k] = b8;
      k = k + 1;
      x[i+k*n] = - w2;
      x[j+k*n] = + u;
      w[k] = b8;
      k = k + 1;
      x[i+k*n] = + w2;
      x[j+k*n] = - u;
      w[k] = b8;
      k = k + 1;
      x[i+k*n] = + w2;
      x[j+k*n] = + u;
      w[k] = b8;
    }
  }
//
//  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
//
  for ( i = 0; i < n - 2; i++ )
  {
    for ( j = i + 1; j < n - 1; j++ )
    {
      for ( l = j + 1; l < n; l++ )
      {
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b9;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b9;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b9;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b9;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b9;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b9;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b9;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b9;
      }
    }
  }
//
//  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
//
  for ( i = 0; i < n - 2; i++ )
  {
    for ( j = i + 1; j < n - 1; j++ )
    {
      for ( l = j + 1; l < n; l++ )
      {
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = - v;
        x[l+k*n] = - v;
        w[k] = b10;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = - v;
        x[l+k*n] = + v;
        w[k] = b10;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = + v;
        x[l+k*n] = - v;
        w[k] = b10;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = + v;
        x[l+k*n] = + v;
        w[k] = b10;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = - v;
        x[l+k*n] = - v;
        w[k] = b10;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = - v;
        x[l+k*n] = + v;
        w[k] = b10;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = + v;
        x[l+k*n] = - v;
        w[k] = b10;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = + v;
        x[l+k*n] = + v;
        w[k] = b10;
      }
    }
  }
//
//  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
//
  for ( i = 0; i < n - 2; i++ )
  {
    for ( j = i + 1; j < n - 1; j++ )
    {
      for ( l = j + 1; l < n; l++ )
      {
        k = k + 1;
        x[i+k*n] = - w2;
        x[j+k*n] = - w2;
        x[l+k*n] = - w2;
        w[k] = b11;
        k = k + 1;
        x[i+k*n] = - w2;
        x[j+k*n] = - w2;
        x[l+k*n] = + w2;
        w[k] = b11;
        k = k + 1;
        x[i+k*n] = - w2;
        x[j+k*n] = + w2;
        x[l+k*n] = - w2;
        w[k] = b11;
        k = k + 1;
        x[i+k*n] = - w2;
        x[j+k*n] = + w2;
        x[l+k*n] = + w2;
        w[k] = b11;
        k = k + 1;
        x[i+k*n] = + w2;
        x[j+k*n] = - w2;
        x[l+k*n] = - w2;
        w[k] = b11;
        k = k + 1;
        x[i+k*n] = + w2;
        x[j+k*n] = - w2;
        x[l+k*n] = + w2;
        w[k] = b11;
        k = k + 1;
        x[i+k*n] = + w2;
        x[j+k*n] = + w2;
        x[l+k*n] = - w2;
        w[k] = b11;
        k = k + 1;
        x[i+k*n] = + w2;
        x[j+k*n] = + w2;
        x[l+k*n] = + w2;
        w[k] = b11;
      }
    }
  }
//
//  8 * ( N * ( N - 1 ) * ( N - 2 ) / 2 ) points.
//
  for ( i = 0; i < n - 2; i++ )
  {
    for ( j = i + 1; j < n - 1; j++ )
    {
      for ( l = j + 1; l < n; l++ )
      {
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = - v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - u;
        x[l+k*n] = + v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = - v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + u;
        x[l+k*n] = + v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = - v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - u;
        x[l+k*n] = + v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = - v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + u;
        x[l+k*n] = + v;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - v;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = - v;
        x[l+k*n] = + u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + v;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - u;
        x[j+k*n] = + v;
        x[l+k*n] = + u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - v;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = - v;
        x[l+k*n] = + u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + v;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + u;
        x[j+k*n] = + v;
        x[l+k*n] = + u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = - v;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = - u;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = - u;
        x[l+k*n] = + u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = + u;
        x[l+k*n] = - u;
        w[k] = b12;
        k = k + 1;
        x[i+k*n] = + v;
        x[j+k*n] = + u;
        x[l+k*n] = + u;
        w[k] = b12;
      }
    }
  }
//
//  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
//
  for ( i = 0; i < n - 3; i++ )
  {
    for ( j = i + 1; j < n - 2; j++ )
    {
      for ( l = j + 1; l < n - 1; l++ )
      {
        for ( m = l + 1; m < n; m++ )
        {
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = - u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = - u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = - u;
          x[m+k*n] = + u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = - u;
          w[k] = b13;
          k = k + 1;
          x[i+k*n] = + u;
          x[j+k*n] = + u;
          x[l+k*n] = + u;
          x[m+k*n] = + u;
          w[k] = b13;
        }
      }
    }
  }
//
//  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
//
  for ( i = 0; i < n - 3; i++ )
  {
    for ( j = i + 1; j < n - 2; j++ )
    {
      for ( l = j + 1; l < n - 1; l++ )
      {
        for ( m = l + 1; m < n; m++ )
        {
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = - v;
          x[l+k*n] = - v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = - v;
          x[l+k*n] = - v;
          x[m+k*n] = + v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = - v;
          x[l+k*n] = + v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = - v;
          x[l+k*n] = + v;
          x[m+k*n] = + v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = + v;
          x[l+k*n] = - v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = + v;
          x[l+k*n] = - v;
          x[m+k*n] = + v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = + v;
          x[l+k*n] = + v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = - v;
          x[j+k*n] = + v;
          x[l+k*n] = + v;
          x[m+k*n] = + v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = - v;
          x[l+k*n] = - v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = - v;
          x[l+k*n] = - v;
          x[m+k*n] = + v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = - v;
          x[l+k*n] = + v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = - v;
          x[l+k*n] = + v;
          x[m+k*n] = + v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = + v;
          x[l+k*n] = - v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = + v;
          x[l+k*n] = - v;
          x[m+k*n] = + v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = + v;
          x[l+k*n] = + v;
          x[m+k*n] = - v;
          w[k] = b14;
          k = k + 1;
          x[i+k*n] = + v;
          x[j+k*n] = + v;
          x[l+k*n] = + v;
          x[m+k*n] = + v;
          w[k] = b14;
        }
      }
    }
  }
//
//  All quintuples UUUUU with 32 sign combinations.
//
  for ( i1 = 0; i1 < n - 4; i1++ )
  {
    for ( i2 = i1 + 1; i2 < n - 3; i2++ )
    {
      for ( i3 = i2 + 1; i3 < n - 2; i3++ )
      {
        for ( i4 = i3 + 1; i4 < n - 1; i4++ )
        {
          for ( i5 = i4 + 1; i5 < n; i5++ )
          {
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = - u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = - u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = - u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = - u;
            x[i5+k*n] = + u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = - u;
            w[k] = b15;
            k = k + 1;
            x[i1+k*n] = + u;
            x[i2+k*n] = + u;
            x[i3+k*n] = + u;
            x[i4+k*n] = + u;
            x[i5+k*n] = + u;
            w[k] = b15;
          }
        }
      }
    }
  }
  return;
}
//****************************************************************************80

int en_r2_11_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_11_1_SIZE sizes the Stroud rule 11.1 for region EN_R2.
//
//  Discussion:
//
//    The rule has order 
//    O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EN_R2_11_1_SIZE, the order.
//
{
  int o;

  o = (   4 * i4_power ( n, 5 )
      -  20 * i4_power ( n, 4 )
      + 140 * i4_power ( n, 3 )
      - 130 * i4_power ( n, 2 ) 
      +  96 *            n 
      +  15 ) / 15;

  return o;
}
//****************************************************************************80

double en_r2_monomial_integral ( int n, int alpha[] )

//****************************************************************************80
//
//  Purpose:
//
//    EN_R2_MONOMIAL_INTEGRAL evaluates monomial integrals in EN_R2.
//
//  Discussion:
//
//    ALPHA is the set of polynomial exponents.
//
//    EN_R2 is the entire N-dimensional space with weight function
//
//      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
//
//    The integral to be evaluated is
//
//      value = integral ( EN ) x(1)^alpha(1) * x(2)^alpha(2) * ... 
//        * x(n)^alpha(n) * w(x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int ALPHA[N], the polynomial exponents.
//    0 <= ALPHA[*].
//
//    Output, double EN_R2_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg;
  int i;
  double value;

  for ( i = 0; i < n; i++ )
  {
    if ( alpha[i] < 0 )
    {
      cerr << "\n";
      cerr << "EN_R2_MONOMIAL_INTEGRAL - Fatal error!\n";
      cerr << "  ALPHA[" << i << "] < 0.\n";
      exit ( 1 );
    }
  }

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    if ( ( alpha[i] % 2 == 1 ) )
    {
      value = 0.0;
      break;
    }
    else
    {
      arg = ( ( double ) ( alpha[i] + 1 ) ) / 2.0;
      value = value * r8_gamma ( arg );
    }
  }

  return value;
}
//****************************************************************************80

double ep1_glg_monomial_integral ( int expon, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EP1_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EP1.
//
//  Discussion:
//
//    EP1_GLG is the interval [0,+oo) with generalized Laguerre weight function:
//
//      w(alpha;x) = x^alpha exp ( - x )
//
//    value = integral ( 0 <= x < +oo ) x^expon x^alpha exp ( - x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, double EP1_GLG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double arg;
  double exact;

  arg = alpha + ( double ) ( expon + 1 );

  exact = r8_gamma ( arg );

  return exact;
}
//****************************************************************************80

double ep1_lag_monomial_integral ( int expon )

//****************************************************************************80
//
//  Purpose:
//
//    EP1_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EP1.
//
//  Discussion:
//
//    EP1 is the interval [0,+oo) with exponential or Laguerre weight function:
//
//      w(x) = exp ( - x )
//
//    value = integral ( 0 <= x < oo ) x^expon exp ( - x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//    0 <= EXPON.
//
//    Output, double EP1_LAG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  double value;

  value = r8_factorial ( expon );

  return value;
}
//****************************************************************************80

void epn_glg_00_1 ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_00_1 implements the "midpoint rule" for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double volume;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "EPN_GLG_00_1 - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  expon = 0;
  volume = ep1_glg_monomial_integral ( expon, alpha );
  volume = pow ( volume, n );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = 1.0;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int epn_glg_00_1_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_00_1_SIZE sizes the midpoint rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, int EPN_GLG_00_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "EPN_GLG_00_1_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void epn_glg_01_1 ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_01_1 implements a precision 1 rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double value1;
  double value2;
  double volume;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "EPN_GLG_01_1 - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  expon = 0;
  value1 = ep1_glg_monomial_integral ( expon, alpha );
  volume = pow ( value1, n );

  expon = 1;
  value2 = ep1_glg_monomial_integral ( expon, alpha );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / value1;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int epn_glg_01_1_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_01_1_SIZE sizes a precision 1 rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, int EPN_GLG_01_1_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "EPN_GLG_01_1_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  o = 1;

  return o;
}
//****************************************************************************80

void epn_glg_02_xiu ( int n, double alpha, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_02_XIU implements the Xiu precision 2 rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "EPN_GLG_02_XIU - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = - 1.0;
  delta0 = alpha + 1.0;
  c1 = - alpha - 1.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = ep1_glg_monomial_integral ( expon, alpha );
  volume = pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int epn_glg_02_xiu_size ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_02_XIU_SIZE sizes the Xiu rule for region EPN_GLG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, int EPN_GLG_02_XIU_SIZE, the order.
//
{
  int o;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "EPN_GLG_02_XIUI_SIZE - Fatal error!\n";
    cerr << "  ALPHA <= -1.0\n";
    exit ( 1 );
  }

  o = n + 1;

  return o;
}
//****************************************************************************80

double epn_glg_monomial_integral ( int n, int expon[], double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_GLG_MONOMIAL_INTEGRAL: integral of monomial with GLG weight on EPN.
//
//  Discussion:
//
//    EPN_GLG is the N-dimensional positive space [0,+oo)^N with generalized
//    Laguerre weight function:
//
//      w(alpha;x) = product ( 1 <= i <= n ) x(i)^alpha exp ( - x(i) )
//
//    value = integral ( EPN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) x(i)^alpha exp ( - x(i) ) dx(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int EXPON[N], the exponents.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//    -1.0 < ALPHA.
//
//    Output, double EPN_GLG_MONOMIAL_INTEGRAL, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = ep1_glg_monomial_integral ( expon[i], alpha );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void epn_lag_00_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_00_1 implements the "midpoint rule" for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double volume;

  expon = 0;
  volume = ep1_lag_monomial_integral ( expon );
  volume = pow ( volume, n );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = 1.0;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int epn_lag_00_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_00_1_SIZE sizes the midpoint rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 0.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EPN_LAG_00_1_SIZE, the order.
//
{
  int o;

  o = 1;

  return o;
}
//****************************************************************************80

void epn_lag_01_1 ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_01_1 implements a precision 1 rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  int expon;
  int i;
  int k;
  double value1;
  double value2;
  double volume;

  expon = 0;
  value1 = ep1_lag_monomial_integral ( expon );
  volume = pow ( value1, n );

  expon = 1;
  value2 = ep1_lag_monomial_integral ( expon );

  r8vec_zero ( n * o, x );

  k = - 1;
//
//  1 point.
//
  k = k + 1;
  for ( i = 0; i < n; i++ )
  {
    x[i+k*n] = value2 / value1;
  }
  w[k] = volume;

  return;
}
//****************************************************************************80

int epn_lag_01_1_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_01_1_SIZE sizes a precision 1 rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order O = 1.
//
//    The rule has precision P = 1.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EPN_LOG_01_1_SIZE, the order.
//
{
  int o;

  o = 1;

  return o;
}
//****************************************************************************80

void epn_lag_02_xiu ( int n, int o, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_02_XIU implements the Xiu precision 2 rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  double c1;
  double coef;
  double delta0;
  int expon;
  double gamma0;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;
  double volume;
  double volume_1d;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );

      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      i = i + 1;
    }
  }

  gamma0 = - 1.0;
  delta0 = 1.0;
  c1 = - 1.0;

  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }

  expon = 0;
  volume_1d = ep1_lag_monomial_integral ( expon );
  volume = pow ( volume_1d, n );

  for ( j = 0; j < o; j++ )
  {
    w[j] = volume / ( double ) ( o );
  }
  return;
}
//****************************************************************************80

int epn_lag_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_02_XIU_SIZE sizes the Xiu rule for region EPN_LAG.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int EPN_LAG_02_XIU_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

double epn_lag_monomial_integral ( int n, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    EPN_LAG_MONOMIAL_INTEGRAL: integral of monomial with Laguerre weight on EPN.
//
//  Discussion:
//
//    EPN is the N-dimensional positive space [0,+oo)^N with exponential 
//    or Laguerre weight function:
//
//      w(x(1:n)) = exp ( - sum ( x(1:n) ) )
//
//    value = integral ( EPN ) 
//      product ( 1 <= i <= n ) x(I)^expon(i) exp ( -x(i) ) dx(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int EXPON(N), the exponents.
//
//    Output, double EPN_LAG_MONOMIAL_VALUE, the value of the integral.
//
{
  int i;
  double value;
  double value2;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value2 = ep1_lag_monomial_integral ( expon[i] );
    value = value * value2;
  }

  return value;
}
//****************************************************************************80

void gw_02_xiu ( int n, int o, double gamma0, double delta0, double c1, 
  double volume_1d, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GW_02_XIU implements the Golub-Welsch version of the Xiu rule.
//
//  Discussion:
//
//    The rule has order 
//
//      O = N + 1.
//
//    The rule has precision P = 2.
//
//    It is assumed that the integral is over an N-dimensional region,
//    and has the form
//
//      Integral f(x) w(x) dx
//
//    where w(x) is separable into identical and independent components:
//
//      w(x) = v(x1) * v(x2) * ... * v(xn)
//
//    Associated with the weight function v(x), we assume there is a
//    family of orthogonal polynomials satisfying a three-term recurrence
//    of the form:
//
//      x P(n,x) = An * P(n+1,x) + Bn * P(n,x) + Cn * P(n-1,x)
//
//    with P(0,x) = 1, and P(-1,x) = 0.
//
//    This routine can construct the desired quadrature rule by knowing
//    the values of C1, used in the definition of P2, the values
//    GAMMA0 = 1/A0 and DELTA0 = - B0/A0, for which it is the case that
//    P(1,X) = GAMMA0 * X + DELTA0, and the value of VOLUME_1D, that is,
//    the 1D integral of v(x) over the region.
//
//    Note the values for the following standard polynomial families:
//
//    Chebyshev Type 1
//      V(X) =      1 / sqrt ( 1 - X^2 )
//      Interval =  [-1,+1]
//      GAMMA0 =    1.0
//      DELTA0 =    0.0
//      C1 =        1/2
//      VOLUME_1D = PI
//
//    Chebyshev Type 2
//      V(X) =      sqrt ( 1 - X^2 )
//      Interval =  [-1,+1]
//      GAMMA0 =    2.0
//      DELTA0 =    0.0
//      C1 =        1/2
//      VOLUME_1D = PI / 2
//
//    Gegenbauer
//      V(X) =      ( 1 - X^2 )^A
//      Interval =  [-1,+1]
//      GAMMA0 =    2 * A + 1
//      DELTA0 =    0.0
//      C1 =        ( 2 * A + 1 ) / ( 2 A + 3 )
//      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
//
//    Gegenbauer* (Removes singularity at ALPHA = -0.5):
//      V(X) =      ( 1 - X^2 )^A
//      Interval =  [-1,+1]
//      GAMMA0 =    1
//      DELTA0 =    0.0
//      C1 =        1 / ( 2 A + 3 )
//      VOLUME_1D = sqrt ( PI ) * Gamma(A+1) / Gamma(A+3/2)
//
//    Generalized Hermite
//      V(X) = |x|^A exp ( - x^2 )
//      Interval = (-oo,+oo)
//      GAMMA0 =    2
//      DELTA0 =    0
//      C1 =        2+2A
//      VOLUME_1D = Gamma((A+1)/2)
//
//    Generalized Laguerre
//      V(X) =       x^A exp ( - x )
//      Interval =  [0,+oo)
//      GAMMA0 =    -1.0
//      DELTA0 =     A+1.0
//      C1 =        -A-1.0
//      VOLUME_1D =  Gamma(A+1)
//
//    Hermite (physicist)
//      V(X) =      exp ( - x^2 )
//      Interval =  (-oo,+oo)
//      GAMMA0 =    2.0
//      DELTA0 =    0.0
//      C1 =        1.0
//      VOLUME_1D = sqrt ( PI )
//
//    Hermite (probabilist)
//      V(X) =      exp ( - x^2 / 2 )
//      Interval =  (-oo,+oo)
//      GAMMA0 =    1.0
//      DELTA0 =    0.0
//      C1 =        1.0
//      VOLUME_1D = sqrt ( 2 PI )
//
//    Jacobi
//      V(X) =      (1-x)^A (1+x)^B
//      Interval =  [-1,+1]
//      GAMMA0 =    (A+B+2)/2  
//      DELTA0 =    (A-B)/2
//      C1 =        2(A+1)(B+1)/(A+B+3)/(A+B+2)
//      VOLUME_1D = 2^(A+B+1) * Gamma(A+1) * Gamma(B+1) / ( A+B+1) / Gamma(A+B+1)
//
//    Laguerre
//      V(X) =       exp ( - x )
//      Interval =  [0,+oo)
//      GAMMA0 =    -1.0
//      DELTA0 =     1.0
//      C1 =        -1.0
//      VOLUME_1D =  1.0
//
//    Legendre
//      V(X) =      1.0
//      Interval =  [-1,+1]
//      GAMMA0 =    1.0
//      DELTA0 =    0.0
//      C1 =        1/3
//      VOLUME_1D = 2.0
//                                  
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, int O, the order.
//
//    Input, double GAMMA0, the ratio 1 / A0.
//
//    Input, double DELTA0, the ratio B0 / A0.
//
//    Input, double C1, the coefficient of P(0,X) in the definition of P(2,X).
//
//    Input, double VOLUME_1D, the 1D integral of V(X).
//
//    Output, double X[N*O], the abscissas.
//
//    Output, double W[O], the weights.
//
{
  double arg;
  int i;
  int j;
  double pi = 3.141592653589793;
  int r;

  for ( j = 0; j < o; j++ )
  {
    i = 0;
    for ( r = 1; r <= ( n / 2 ); r++ )
    {
      arg = ( double ) ( 2 * r * j ) * pi / ( double ) ( n + 1 );
      x[i+j*n] = sqrt ( 2.0 ) * cos ( arg );
      i = i + 1;
      x[i+j*n] = sqrt ( 2.0 ) * sin ( arg );
      i = i + 1;
    }

    if ( i < n )
    {
      x[i+j*n] = r8_mop ( j );
      i = i + 1;
    }
  }
//
//  Adjust for the GW rule.
//
  for ( j = 0; j < o; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = ( sqrt ( gamma0 * c1 ) * x[i+j*n] - delta0 ) / gamma0;
    }
  }
//
//  The weights are equal.
//
  for ( j = 0; j < o; j++ )
  {
    w[j] = pow ( volume_1d, n ) / ( double ) ( o );
  }

  return;
}
//****************************************************************************80

int gw_02_xiu_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    GW_02_XIU_SIZE sizes the Golub Welsch version of the Xiu rule.
//
//  Discussion:
//
//    The rule has order O = N + 1;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical integration formulas of degree two,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 1515-1520.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, int GW_02_XIU_SIZE, the order.
//
{
  int o;

  o = n + 1;

  return o;
}
//****************************************************************************80

double hexagon_area_2d ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
//
//  Discussion:
//
//    The formula for the area only requires the radius, and does
//    not depend on the location of the center, or the orientation
//    of the hexagon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the hexagon.
//
//    Output, double HEXAGON_AREA_2D, the area of the hexagon.
//
{
  double value;

  value = r * r * hexagon_unit_area_2d ( );

  return value;
}
//****************************************************************************80

double hexagon_sum ( double func ( double x, double y ), double center[2], 
  double r, int order, double xtab[], double ytab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEXAGON_SUM applies a quadrature rule inside a hexagon in 2D.
//
//  Discussion:
//
//    The input quadrature rule is assumed to be defined for a unit hexagon.
//
//    The input quadrature rule may be defined by calling HEXAGON_UNIT_SET.
//
//  Integration region:
//
//    The definition is given in terms of THETA, the angle in degrees of the
//    vector (X-CENTER(1),Y-CENTER(2)).  The following six conditions apply,
//    respectively, between the bracketing values of THETA of 0, 60, 120, 
//    180, 240, 300, and 360.
//
//      0 <= Y-CENTER(2) <= -SQRT(3) * (X-CENTER(1)) + R * SQRT(3)
//      0 <= Y-CENTER(2) <=                     R * SQRT(3)/2
//      0 <= Y-CENTER(2) <=  SQRT(3) * (X-CENTER(1)) + R * SQRT(3) 
//      -SQRT(3) * (X-CENTER(1)) - R * SQRT(3)	<= Y-CENTER(2) <= 0
//                        - R * SQRT(3)/2 <= Y-CENTER(2) <= 0
//       SQRT(3) * (X-CENTER(1)) - R * SQRT(3)   <= Y-CENTER(2) <= 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function of two variables which is to be integrated.
//
//    Input, double CENTER[2], the center of the hexagon.
//
//    Input, double R, the radius of the hexagon.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;
  double x;
  double y;

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    x = center[0] + r * xtab[i];
    y = center[1] + r * ytab[i];
    quad = quad + weight[i] * func ( x, y );
  }

  volume = hexagon_area_2d ( r );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double hexagon_unit_area_2d ( )

//****************************************************************************80
//
//  Purpose:
//
//    HEXAGON_UNIT_AREA_2D returns the area of the unit regular hexagon in 2D.
//
//  Integration region:
//
//    The definition is given in terms of THETA, the angle in degrees of the
//    vector (X,Y).  The following six conditions apply, respectively,
//    between the bracketing values of THETA of 0, 60, 120, 180, 240,
//    300, and 360.
//
//                              0 <= Y <= -SQRT(3) * X + SQRT(3)
//                              0 <= Y <=                 SQRT(3)/2
//                              0 <= Y <=  SQRT(3) * X + SQRT(3)
//      - SQRT(3) * X - SQRT(3)   <= Y <= 0
//                    - SQRT(3)/2 <= Y <= 0
//        SQRT(3) * X - SQRT(3)   <= Y <= 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double HEXAGON_UNIT_AREA_2D, the area of the hexagon.
//
{
  double value;

  value = 3.0 * sqrt ( 3.0 ) / 2.0;

  return value;
}
//****************************************************************************80

void hexagon_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEXAGON_UNIT_SET sets a quadrature rule inside the unit hexagon in 2D.
//
//  Integration region:
//
//    The definition is given in terms of THETA, the angle in degrees of the
//    vector (X,Y).  The following six conditions apply, respectively,
//    between the bracketing values of THETA of 0, 60, 120, 180, 240,
//    300, and 360.
//
//                              0 <= Y <= -SQRT(3) * X + SQRT(3)
//                              0 <= Y <=                 SQRT(3)/2
//                              0 <= Y <=  SQRT(3) * X + SQRT(3)
//       -SQRT(3) * X - SQRT(3)   <= Y <= 0
//                    - SQRT(3)/2 <= Y <= 0
//        SQRT(3) * X - SQRT(3)   <= Y <= 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
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
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule desired.
//      1, 1 point,  degree 1;
//      2, 4 points, degree 3;
//      3, 7 points, degree 3;
//      4, 7 points, degree 5;
//
//    Input, int ORDER, the order of the desired rule.
//
//    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[ORDER], the weights of the rule.
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double z;

  if ( rule == 1 )
  {
    xtab[0] = 0.0;
    ytab[0] = 0.0;
    weight[0] = 1.0;
  }
//
//  Stroud rule H2:3-1.
//
  else if ( rule == 2 )
  {
    a = sqrt ( 5.0 / 12.0 );
    b = 1.0 / 4.0;
    z = 0.0;

    xtab[0] =   a;
    xtab[1] = - a;
    xtab[2] =   z;
    xtab[3] =   z;

    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   a;
    ytab[3] = - a;

    weight[0] = b;
    weight[1] = b;
    weight[2] = b;
    weight[3] = b;
  }
//
//  Stroud rule H2:3-2.
//
  else if ( rule == 3 )
  {
    a = sqrt ( 3.0 ) / 2.0;
    b =  0.5;
    c =  1.0;
    d =  5.0 / 72.0;
    e = 42.0 / 72.0;
    z =  0.0;

    xtab[0] =   z;
    xtab[1] =   c;
    xtab[2] = - c;
    xtab[3] =   b;
    xtab[4] = - b;
    xtab[5] =   b;
    xtab[6] = - b;
 
    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   z;
    ytab[3] =   a;
    ytab[4] =   a;
    ytab[5] = - a;
    ytab[6] = - a;

    weight[0] =   e;
    weight[1] =   d;
    weight[2] =   d;
    weight[3] =   d;
    weight[4] =   d;
    weight[5] =   d;
    weight[6] =   d;
  }
//
//  Stroud rule H2:5-1.
//
  else if ( rule == 4 )
  {
    a = sqrt ( 14.0 ) / 5.0;
    b = sqrt ( 14.0 ) / 10.0;
    c = sqrt ( 42.0 ) / 10.0;
    d = 125.0 / 1008.0;
    e = 258.0 / 1008.0;
    z = 0.0;

    xtab[0] =   z;
    xtab[1] =   a;
    xtab[2] = - a;
    xtab[3] =   b;
    xtab[4] = - b;
    xtab[5] =   b;
    xtab[6] = - b;

    ytab[0] =   z;
    ytab[1] =   z;
    ytab[2] =   z;
    ytab[3] =   c;
    ytab[4] =   c;
    ytab[5] = - c;
    ytab[6] = - c;

    weight[0] =   e;
    weight[1] =   d;
    weight[2] =   d;
    weight[3] =   d;
    weight[4] =   d;
    weight[5] =   d;
    weight[6] =   d;
  }
  else
  {
    cerr << "\n";
    cerr << "HEXAGON_UNIT_SET - Fatal error!\n";
    cerr << "  Illegal input value of RULE = " << rule << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

int hexagon_unit_size ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    HEXAGON_UNIT_SIZE sizes a quadrature rule inside the unit hexagon in 2D.
//
//  Integration region:
//
//    The definition is given in terms of THETA, the angle in degrees of the
//    vector (X,Y).  The following six conditions apply, respectively,
//    between the bracketing values of THETA of 0, 60, 120, 180, 240,
//    300, and 360.
//
//                              0 <= Y <= -SQRT(3) * X + SQRT(3)
//                              0 <= Y <=                 SQRT(3)/2
//                              0 <= Y <=  SQRT(3) * X + SQRT(3)
//       -SQRT(3) * X - SQRT(3)   <= Y <= 0
//                    - SQRT(3)/2 <= Y <= 0
//        SQRT(3) * X - SQRT(3)   <= Y <= 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
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
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule desired.
//      1, 1 point,  degree 1;
//      2, 4 points, degree 3;
//      3, 7 points, degree 3;
//      4, 7 points, degree 5;
//
//    Output, int HEXAGON_UNIT_SIZE, the order of the desired rule.
//    If RULE is not legal, then ORDER is returned as -1.
//
{
  int order;

  if ( rule == 1 )
  {
    order = 1;
  }
//
//  Stroud rule H2:3-1.
//
  else if ( rule == 2 )
  {
    order = 4;
  }
//
//  Stroud rule H2:3-2.
//
  else if ( rule == 3 )
  {
    order = 7;
  }
//
//  Stroud rule H2:5-1.
//
  else if ( rule == 4 )
  {
    order = 7;
  }
  else
  {
    order = -1;
  }

  return order;
}
//****************************************************************************80

int i4_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL returns N!.
//
//  Discussion:
//
//    N! = Product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    0 <= N.
//
//    Output, int I4_FACTORIAL, the factorial of N.
//
{
  int fact;
  int i;
//
//  Check.
//
  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "I4_FACTORIAL - Fatal error!\n";
    cerr << "  N < 0.\n";
    return 0;
  }

  fact = 1;

  for ( i = 2; i <= n; i++ )
  {
    fact = fact * i;
  }

  return fact;
}
//****************************************************************************80

int i4_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL2 computes the double factorial function.
//
//  Discussion:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial function.
//    If N is less than 1, I4_FACTORIAL2 is returned as 1.
//
//    Output, int I4_FACTORIAL2, the value of the double factorial function.
//
{
  int n_copy;
  int value;

  if ( n < 1 )
  {
    value = 1;
    return value;
  }

  n_copy = n;
  value = 1;

  while ( 1 < n_copy )
  {
    value = value * n_copy;
    n_copy = n_copy - 2;
  }

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

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

void i4vec_zero ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
//****************************************************************************80

void ksub_next2 ( int n, int k, int a[], int *in, int *iout )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_NEXT2 generates the subsets of size K from a set of size N.
//
//  Discussion:
//
//    This routine uses the revolving door method.  It has no "memory".
//    It simply calculates the successor of the input set,
//    and will start from the beginning after the last set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//    N must be positive.
//
//    Input, int K, the size of the desired subset.  K must be
//    between 0 and N.
//
//    Input/output, int A[K].  On input, the user must
//    supply a subset of size K in A.  That is, A must
//    contain K unique numbers, in order, between 1 and N.  On
//    output, A(I) is the I-th element of the output subset.
//    The output array is also in sorted order.
//
//    Output, int *IN, the element of the output subset which
//    was not in the input set.  Each new subset differs from the
//    last one by adding one element and deleting another.
//
//    Output, int *IOUT, the element of the input subset which
//    is not in the output subset.
//
{
  int j;
  int m;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT2 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  but 0 < N is required!\n";
    exit ( 1 );
  }

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_NEXT2 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  j = 0;

  for ( ; ; )
  {
    if ( 0 < j || ( k % 2 ) == 0 )
    {
      j = j + 1;

      if ( k < j )
      {
        a[k-1] = k;
        *in = k;
        *iout = n;
        return;
      }

      if ( a[j-1] != j )
      {
        *iout = a[j-1];
        *in = *iout - 1;
        a[j-1] = *in;

        if ( j != 1 )
        {
          *in = j - 1;
          a[j-2] = *in;
        }

        return;

      }

    }

    j = j + 1;
    m = n;

    if ( j < k )
    {
      m = a[j] - 1;
    }

    if ( m != a[j-1] )
    {
      break;
    }

  }

  *in = a[j-1] + 1;
  a[j-1] = *in;
  *iout = *in - 1;

  if ( j != 1 )
  {
    a[j-2] = *iout;
    *iout = j - 1;
  }

  return;
}
//****************************************************************************80

void legendre_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function w(x-1] = 1.0;
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    Quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//    The quadrature rule will integrate exactly all polynomials up to
//    X**(2*N-1).
//
//    The abscissas of the rule are the zeroes of the Legendre polynomial
//    P(N)(X).
//
//    The integral produced by a Gauss-Legendre rule is equal to the
//    integral of the unique polynomial of degree N-1 which
//    agrees with the function at the N abscissas of the rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2009
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
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//    N must be between 1 and 33, 63, 64, 65, 127 or 255.
//
//    Output, double X[N], the abscissas of the rule.
//
//    Output, double W[N], the weights of the rule.
//    The weights are positive, symmetric and should sum to 2.
//
{
  if ( n == 1 )
  {
    x[0] =   0.0;

    w[0] = 2.0;
  }
  else if ( n == 2 )
  {
    x[0] = - 0.577350269189625764509148780502;
    x[1] =   0.577350269189625764509148780502;

    w[0] = 1.0;
    w[1] = 1.0;
  }
  else if ( n == 3 )
  {
    x[0] = - 0.774596669241483377035853079956;
    x[1] =   0.0;
    x[2] =   0.774596669241483377035853079956;

    w[0] = 5.0 / 9.0;
    w[1] = 8.0 / 9.0;
    w[2] = 5.0 / 9.0;
  }
  else if ( n == 4 )
  {
    x[0] = - 0.861136311594052575223946488893;
    x[1] = - 0.339981043584856264802665759103;
    x[2] =   0.339981043584856264802665759103;
    x[3] =   0.861136311594052575223946488893;

    w[0] = 0.347854845137453857373063949222;
    w[1] = 0.652145154862546142626936050778;
    w[2] = 0.652145154862546142626936050778;
    w[3] = 0.347854845137453857373063949222;
  }
  else if ( n == 5 )
  {
    x[0] = - 0.906179845938663992797626878299;
    x[1] = - 0.538469310105683091036314420700;
    x[2] =   0.0;
    x[3] =   0.538469310105683091036314420700;
    x[4] =   0.906179845938663992797626878299;

    w[0] = 0.236926885056189087514264040720;
    w[1] = 0.478628670499366468041291514836;
    w[2] = 0.568888888888888888888888888889;
    w[3] = 0.478628670499366468041291514836;
    w[4] = 0.236926885056189087514264040720;
  }
  else if ( n == 6 )
  {
    x[0] = - 0.932469514203152027812301554494;
    x[1] = - 0.661209386466264513661399595020;
    x[2] = - 0.238619186083196908630501721681;
    x[3] =   0.238619186083196908630501721681;
    x[4] =   0.661209386466264513661399595020;
    x[5] =   0.932469514203152027812301554494;

    w[0] = 0.171324492379170345040296142173;
    w[1] = 0.360761573048138607569833513838;
    w[2] = 0.467913934572691047389870343990;
    w[3] = 0.467913934572691047389870343990;
    w[4] = 0.360761573048138607569833513838;
    w[5] = 0.171324492379170345040296142173;
  }
  else if ( n == 7 )
  {
    x[0] = - 0.949107912342758524526189684048;
    x[1] = - 0.741531185599394439863864773281;
    x[2] = - 0.405845151377397166906606412077;
    x[3] =   0.0;
    x[4] =   0.405845151377397166906606412077;
    x[5] =   0.741531185599394439863864773281;
    x[6] =   0.949107912342758524526189684048;

    w[0] = 0.129484966168869693270611432679;
    w[1] = 0.279705391489276667901467771424;
    w[2] = 0.381830050505118944950369775489;
    w[3] = 0.417959183673469387755102040816;
    w[4] = 0.381830050505118944950369775489;
    w[5] = 0.279705391489276667901467771424;
    w[6] = 0.129484966168869693270611432679;
  }
  else if ( n == 8 )
  {
    x[0] = - 0.960289856497536231683560868569;
    x[1] = - 0.796666477413626739591553936476;
    x[2] = - 0.525532409916328985817739049189;
    x[3] = - 0.183434642495649804939476142360;
    x[4] =   0.183434642495649804939476142360;
    x[5] =   0.525532409916328985817739049189;
    x[6] =   0.796666477413626739591553936476;
    x[7] =   0.960289856497536231683560868569;

    w[0] = 0.101228536290376259152531354310;
    w[1] = 0.222381034453374470544355994426;
    w[2] = 0.313706645877887287337962201987;
    w[3] = 0.362683783378361982965150449277;
    w[4] = 0.362683783378361982965150449277;
    w[5] = 0.313706645877887287337962201987;
    w[6] = 0.222381034453374470544355994426;
    w[7] = 0.101228536290376259152531354310;
  }
  else if ( n == 9 )
  {
    x[0] = - 0.968160239507626089835576202904;
    x[1] = - 0.836031107326635794299429788070;
    x[2] = - 0.613371432700590397308702039341;
    x[3] = - 0.324253423403808929038538014643;
    x[4] =   0.0;
    x[5] =   0.324253423403808929038538014643;
    x[6] =   0.613371432700590397308702039341;
    x[7] =   0.836031107326635794299429788070;
    x[8] =   0.968160239507626089835576202904;

    w[0] = 0.812743883615744119718921581105E-01;
    w[1] = 0.180648160694857404058472031243;
    w[2] = 0.260610696402935462318742869419;
    w[3] = 0.312347077040002840068630406584;
    w[4] = 0.330239355001259763164525069287;
    w[5] = 0.312347077040002840068630406584;
    w[6] = 0.260610696402935462318742869419;
    w[7] = 0.180648160694857404058472031243;
    w[8] = 0.812743883615744119718921581105E-01;
  }
  else if ( n == 10 )
  {
    x[0] =  - 0.973906528517171720077964012084;
    x[1] =  - 0.865063366688984510732096688423;
    x[2] =  - 0.679409568299024406234327365115;
    x[3] =  - 0.433395394129247190799265943166;
    x[4] =  - 0.148874338981631210884826001130;
    x[5] =    0.148874338981631210884826001130;
    x[6] =    0.433395394129247190799265943166;
    x[7] =    0.679409568299024406234327365115;
    x[8] =    0.865063366688984510732096688423;
    x[9] =   0.973906528517171720077964012084;

    w[0] =  0.666713443086881375935688098933E-01;
    w[1] =  0.149451349150580593145776339658;
    w[2] =  0.219086362515982043995534934228;
    w[3] =  0.269266719309996355091226921569;
    w[4] =  0.295524224714752870173892994651;
    w[5] =  0.295524224714752870173892994651;
    w[6] =  0.269266719309996355091226921569;
    w[7] =  0.219086362515982043995534934228;
    w[8] =  0.149451349150580593145776339658;
    w[9] = 0.666713443086881375935688098933E-01;
  }
  else if ( n == 11 )
  {
    x[0] =  - 0.978228658146056992803938001123;
    x[1] =  - 0.887062599768095299075157769304;
    x[2] =  - 0.730152005574049324093416252031;
    x[3] =  - 0.519096129206811815925725669459;
    x[4] =  - 0.269543155952344972331531985401;
    x[5] =    0.0;
    x[6] =    0.269543155952344972331531985401;
    x[7] =    0.519096129206811815925725669459;
    x[8] =    0.730152005574049324093416252031;
    x[9] =   0.887062599768095299075157769304;
    x[10] =   0.978228658146056992803938001123;

    w[0] =  0.556685671161736664827537204425E-01;
    w[1] =  0.125580369464904624634694299224;
    w[2] =  0.186290210927734251426097641432;
    w[3] =  0.233193764591990479918523704843;
    w[4] =  0.262804544510246662180688869891;
    w[5] =  0.272925086777900630714483528336;
    w[6] =  0.262804544510246662180688869891;
    w[7] =  0.233193764591990479918523704843;
    w[8] =  0.186290210927734251426097641432;
    w[9] = 0.125580369464904624634694299224;
    w[10] = 0.556685671161736664827537204425E-01;
  }
  else if ( n == 12 )
  {
    x[0] =  - 0.981560634246719250690549090149;
    x[1] =  - 0.904117256370474856678465866119;
    x[2] =  - 0.769902674194304687036893833213;
    x[3] =  - 0.587317954286617447296702418941;
    x[4] =  - 0.367831498998180193752691536644;
    x[5] =  - 0.125233408511468915472441369464;
    x[6] =    0.125233408511468915472441369464;
    x[7] =    0.367831498998180193752691536644;
    x[8] =    0.587317954286617447296702418941;
    x[9] =   0.769902674194304687036893833213;
    x[10] =   0.904117256370474856678465866119;
    x[11] =   0.981560634246719250690549090149;

    w[0] =  0.471753363865118271946159614850E-01;
    w[1] =  0.106939325995318430960254718194;
    w[2] =  0.160078328543346226334652529543;
    w[3] =  0.203167426723065921749064455810;
    w[4] =  0.233492536538354808760849898925;
    w[5] =  0.249147045813402785000562436043;
    w[6] =  0.249147045813402785000562436043;
    w[7] =  0.233492536538354808760849898925;
    w[8] =  0.203167426723065921749064455810;
    w[9] = 0.160078328543346226334652529543;
    w[10] = 0.106939325995318430960254718194;
    w[11] = 0.471753363865118271946159614850E-01;
  }
  else if ( n == 13 )
  {
    x[0] =  - 0.984183054718588149472829448807;
    x[1] =  - 0.917598399222977965206547836501;
    x[2] =  - 0.801578090733309912794206489583;
    x[3] =  - 0.642349339440340220643984606996;
    x[4] =  - 0.448492751036446852877912852128;
    x[5] =  - 0.230458315955134794065528121098;
    x[6] =    0.0;
    x[7] =    0.230458315955134794065528121098;
    x[8] =    0.448492751036446852877912852128;
    x[9] =   0.642349339440340220643984606996;
    x[10] =   0.801578090733309912794206489583;
    x[11] =   0.917598399222977965206547836501;
    x[12] =   0.984183054718588149472829448807;

    w[0] =  0.404840047653158795200215922010E-01;
    w[1] =  0.921214998377284479144217759538E-01;
    w[2] =  0.138873510219787238463601776869;
    w[3] =  0.178145980761945738280046691996;
    w[4] =  0.207816047536888502312523219306;
    w[5] =  0.226283180262897238412090186040;
    w[6] =  0.232551553230873910194589515269;
    w[7] =  0.226283180262897238412090186040;
    w[8] =  0.207816047536888502312523219306;
    w[9] = 0.178145980761945738280046691996;
    w[10] = 0.138873510219787238463601776869;
    w[11] = 0.921214998377284479144217759538E-01;
    w[12] = 0.404840047653158795200215922010E-01;
  }
  else if ( n == 14 )
  {
    x[0] =  - 0.986283808696812338841597266704;
    x[1] =  - 0.928434883663573517336391139378;
    x[2] =  - 0.827201315069764993189794742650;
    x[3] =  - 0.687292904811685470148019803019;
    x[4] =  - 0.515248636358154091965290718551;
    x[5] =  - 0.319112368927889760435671824168;
    x[6] =  - 0.108054948707343662066244650220;
    x[7] =    0.108054948707343662066244650220;
    x[8] =    0.319112368927889760435671824168;
    x[9] =   0.515248636358154091965290718551;
    x[10] =   0.687292904811685470148019803019;
    x[11] =   0.827201315069764993189794742650;
    x[12] =   0.928434883663573517336391139378;
    x[13] =   0.986283808696812338841597266704;

    w[0] =  0.351194603317518630318328761382E-01;
    w[1] =  0.801580871597602098056332770629E-01;
    w[2] =  0.121518570687903184689414809072;
    w[3] =  0.157203167158193534569601938624;
    w[4] =  0.185538397477937813741716590125;
    w[5] =  0.205198463721295603965924065661;
    w[6] =  0.215263853463157790195876443316;
    w[7] =  0.215263853463157790195876443316;
    w[8] =  0.205198463721295603965924065661;
    w[9] = 0.185538397477937813741716590125;
    w[10] = 0.157203167158193534569601938624;
    w[11] = 0.121518570687903184689414809072;
    w[12] = 0.801580871597602098056332770629E-01;
    w[13] = 0.351194603317518630318328761382E-01;
  }
  else if ( n == 15 )
  {
    x[0] =  - 0.987992518020485428489565718587;
    x[1] =  - 0.937273392400705904307758947710;
    x[2] =  - 0.848206583410427216200648320774;
    x[3] =  - 0.724417731360170047416186054614;
    x[4] =  - 0.570972172608538847537226737254;
    x[5] =  - 0.394151347077563369897207370981;
    x[6] =  - 0.201194093997434522300628303395;
    x[7] =    0.0;
    x[8] =    0.201194093997434522300628303395;
    x[9] =   0.394151347077563369897207370981;
    x[10] =   0.570972172608538847537226737254;
    x[11] =   0.724417731360170047416186054614;
    x[12] =   0.848206583410427216200648320774;
    x[13] =   0.937273392400705904307758947710;
    x[14] =   0.987992518020485428489565718587;

    w[0] =  0.307532419961172683546283935772E-01;
    w[1] =  0.703660474881081247092674164507E-01;
    w[2] =  0.107159220467171935011869546686;
    w[3] =  0.139570677926154314447804794511;
    w[4] =  0.166269205816993933553200860481;
    w[5] =  0.186161000015562211026800561866;
    w[6] =  0.198431485327111576456118326444;
    w[7] =  0.202578241925561272880620199968;
    w[8] =  0.198431485327111576456118326444;
    w[9] = 0.186161000015562211026800561866;
    w[10] = 0.166269205816993933553200860481;
    w[11] = 0.139570677926154314447804794511;
    w[12] = 0.107159220467171935011869546686;
    w[13] = 0.703660474881081247092674164507E-01;
    w[14] = 0.307532419961172683546283935772E-01;
  }
  else if ( n == 16 )
  {
    x[0] =  - 0.989400934991649932596154173450;
    x[1] =  - 0.944575023073232576077988415535;
    x[2] =  - 0.865631202387831743880467897712;
    x[3] =  - 0.755404408355003033895101194847;
    x[4] =  - 0.617876244402643748446671764049;
    x[5] =  - 0.458016777657227386342419442984;
    x[6] =  - 0.281603550779258913230460501460;
    x[7] =  - 0.950125098376374401853193354250E-01;
    x[8] =    0.950125098376374401853193354250E-01;
    x[9] =   0.281603550779258913230460501460;
    x[10] =   0.458016777657227386342419442984;
    x[11] =   0.617876244402643748446671764049;
    x[12] =   0.755404408355003033895101194847;
    x[13] =   0.865631202387831743880467897712;
    x[14] =   0.944575023073232576077988415535;
    x[15] =   0.989400934991649932596154173450;

    w[0] =  0.271524594117540948517805724560E-01;
    w[1] =  0.622535239386478928628438369944E-01;
    w[2] =  0.951585116824927848099251076022E-01;
    w[3] =  0.124628971255533872052476282192;
    w[4] =  0.149595988816576732081501730547;
    w[5] =  0.169156519395002538189312079030;
    w[6] =  0.182603415044923588866763667969;
    w[7] =  0.189450610455068496285396723208;
    w[8] =  0.189450610455068496285396723208;
    w[9] = 0.182603415044923588866763667969;
    w[10] = 0.169156519395002538189312079030;
    w[11] = 0.149595988816576732081501730547;
    w[12] = 0.124628971255533872052476282192;
    w[13] = 0.951585116824927848099251076022E-01;
    w[14] = 0.622535239386478928628438369944E-01;
    w[15] = 0.271524594117540948517805724560E-01;
  }
  else if ( n == 17 )
  {
    x[0] =  - 0.990575475314417335675434019941;
    x[1] =  - 0.950675521768767761222716957896;
    x[2] =  - 0.880239153726985902122955694488;
    x[3] =  - 0.781514003896801406925230055520;
    x[4] =  - 0.657671159216690765850302216643;
    x[5] =  - 0.512690537086476967886246568630;
    x[6] =  - 0.351231763453876315297185517095;
    x[7] =  - 0.178484181495847855850677493654;
    x[8] =    0.0;
    x[9] =   0.178484181495847855850677493654;
    x[10] =   0.351231763453876315297185517095;
    x[11] =   0.512690537086476967886246568630;
    x[12] =   0.657671159216690765850302216643;
    x[13] =   0.781514003896801406925230055520;
    x[14] =   0.880239153726985902122955694488;
    x[15] =   0.950675521768767761222716957896;
    x[16] =   0.990575475314417335675434019941;

    w[0] =  0.241483028685479319601100262876E-01;
    w[1] =  0.554595293739872011294401653582E-01;
    w[2] =  0.850361483171791808835353701911E-01;
    w[3] =  0.111883847193403971094788385626;
    w[4] =  0.135136368468525473286319981702;
    w[5] =  0.154045761076810288081431594802;
    w[6] =  0.168004102156450044509970663788;
    w[7] =  0.176562705366992646325270990113;
    w[8] =  0.179446470356206525458265644262;
    w[9] = 0.176562705366992646325270990113;
    w[10] = 0.168004102156450044509970663788;
    w[11] = 0.154045761076810288081431594802;
    w[12] = 0.135136368468525473286319981702;
    w[13] = 0.111883847193403971094788385626;
    w[14] = 0.850361483171791808835353701911E-01;
    w[15] = 0.554595293739872011294401653582E-01;
    w[16] = 0.241483028685479319601100262876E-01;
  }
  else if ( n == 18 )
  {
    x[0] =  - 0.991565168420930946730016004706;
    x[1] =  - 0.955823949571397755181195892930;
    x[2] =  - 0.892602466497555739206060591127;
    x[3] =  - 0.803704958972523115682417455015;
    x[4] =  - 0.691687043060353207874891081289;
    x[5] =  - 0.559770831073947534607871548525;
    x[6] =  - 0.411751161462842646035931793833;
    x[7] =  - 0.251886225691505509588972854878;
    x[8] =  - 0.847750130417353012422618529358E-01;
    x[9] =   0.847750130417353012422618529358E-01;
    x[10] =   0.251886225691505509588972854878;
    x[11] =   0.411751161462842646035931793833;
    x[12] =   0.559770831073947534607871548525;
    x[13] =   0.691687043060353207874891081289;
    x[14] =   0.803704958972523115682417455015;
    x[15] =   0.892602466497555739206060591127;
    x[16] =   0.955823949571397755181195892930;
    x[17] =   0.991565168420930946730016004706;

    w[0] =  0.216160135264833103133427102665E-01;
    w[1] =  0.497145488949697964533349462026E-01;
    w[2] =  0.764257302548890565291296776166E-01;
    w[3] =  0.100942044106287165562813984925;
    w[4] =  0.122555206711478460184519126800;
    w[5] =  0.140642914670650651204731303752;
    w[6] =  0.154684675126265244925418003836;
    w[7] =  0.164276483745832722986053776466;
    w[8] =  0.169142382963143591840656470135;
    w[9] = 0.169142382963143591840656470135;
    w[10] = 0.164276483745832722986053776466;
    w[11] = 0.154684675126265244925418003836;
    w[12] = 0.140642914670650651204731303752;
    w[13] = 0.122555206711478460184519126800;
    w[14] = 0.100942044106287165562813984925;
    w[15] = 0.764257302548890565291296776166E-01;
    w[16] = 0.497145488949697964533349462026E-01;
    w[17] = 0.216160135264833103133427102665E-01;
  }
  else if ( n == 19 )
  {
    x[0] =  - 0.992406843843584403189017670253;
    x[1] =  - 0.960208152134830030852778840688;
    x[2] =  - 0.903155903614817901642660928532;
    x[3] =  - 0.822714656537142824978922486713;
    x[4] =  - 0.720966177335229378617095860824;
    x[5] =  - 0.600545304661681023469638164946;
    x[6] =  - 0.464570741375960945717267148104;
    x[7] =  - 0.316564099963629831990117328850;
    x[8] =  - 0.160358645640225375868096115741;
    x[9] =   0.0;
    x[10] =   0.160358645640225375868096115741;
    x[11] =   0.316564099963629831990117328850;
    x[12] =   0.464570741375960945717267148104;
    x[13] =   0.600545304661681023469638164946;
    x[14] =   0.720966177335229378617095860824;
    x[15] =   0.822714656537142824978922486713;
    x[16] =   0.903155903614817901642660928532;
    x[17] =   0.960208152134830030852778840688;
    x[18] =   0.992406843843584403189017670253;

    w[0] =  0.194617882297264770363120414644E-01;
    w[1] =  0.448142267656996003328381574020E-01;
    w[2] =  0.690445427376412265807082580060E-01;
    w[3] =  0.914900216224499994644620941238E-01;
    w[4] =  0.111566645547333994716023901682;
    w[5] =  0.128753962539336227675515784857;
    w[6] =  0.142606702173606611775746109442;
    w[7] =  0.152766042065859666778855400898;
    w[8] =  0.158968843393954347649956439465;
    w[9] = 0.161054449848783695979163625321;
    w[10] = 0.158968843393954347649956439465;
    w[11] = 0.152766042065859666778855400898;
    w[12] = 0.142606702173606611775746109442;
    w[13] = 0.128753962539336227675515784857;
    w[14] = 0.111566645547333994716023901682;
    w[15] = 0.914900216224499994644620941238E-01;
    w[16] = 0.690445427376412265807082580060E-01;
    w[17] = 0.448142267656996003328381574020E-01;
    w[18] = 0.194617882297264770363120414644E-01;
  }
  else if ( n == 20 )
  {
    x[0] =  - 0.993128599185094924786122388471;
    x[1] =  - 0.963971927277913791267666131197;
    x[2] =  - 0.912234428251325905867752441203;
    x[3] =  - 0.839116971822218823394529061702;
    x[4] =  - 0.746331906460150792614305070356;
    x[5] =  - 0.636053680726515025452836696226;
    x[6] =  - 0.510867001950827098004364050955;
    x[7] =  - 0.373706088715419560672548177025;
    x[8] =  - 0.227785851141645078080496195369;
    x[9] = - 0.765265211334973337546404093988E-01;
    x[10] =   0.765265211334973337546404093988E-01;
    x[11] =   0.227785851141645078080496195369;
    x[12] =   0.373706088715419560672548177025;
    x[13] =   0.510867001950827098004364050955;
    x[14] =   0.636053680726515025452836696226;
    x[15] =   0.746331906460150792614305070356;
    x[16] =   0.839116971822218823394529061702;
    x[17] =   0.912234428251325905867752441203;
    x[18] =   0.963971927277913791267666131197;
    x[19] =   0.993128599185094924786122388471;

    w[0] =  0.176140071391521183118619623519E-01;
    w[1] =  0.406014298003869413310399522749E-01;
    w[2] =  0.626720483341090635695065351870E-01;
    w[3] =  0.832767415767047487247581432220E-01;
    w[4] =  0.101930119817240435036750135480;
    w[5] =  0.118194531961518417312377377711;
    w[6] =  0.131688638449176626898494499748;
    w[7] =  0.142096109318382051329298325067;
    w[8] =  0.149172986472603746787828737002;
    w[9] = 0.152753387130725850698084331955;
    w[10] = 0.152753387130725850698084331955;
    w[11] = 0.149172986472603746787828737002;
    w[12] = 0.142096109318382051329298325067;
    w[13] = 0.131688638449176626898494499748;
    w[14] = 0.118194531961518417312377377711;
    w[15] = 0.101930119817240435036750135480;
    w[16] = 0.832767415767047487247581432220E-01;
    w[17] = 0.626720483341090635695065351870E-01;
    w[18] = 0.406014298003869413310399522749E-01;
    w[19] = 0.176140071391521183118619623519E-01;
  }
  else if ( n == 21 )
  {
    x[ 0] =  -0.9937521706203896E+00;
    x[ 1] =  -0.9672268385663063E+00;
    x[ 2] =  -0.9200993341504008E+00;
    x[ 3] =  -0.8533633645833173E+00;
    x[ 4] =  -0.7684399634756779E+00;
    x[ 5] =  -0.6671388041974123E+00;
    x[ 6] =  -0.5516188358872198E+00;
    x[ 7] =  -0.4243421202074388E+00;
    x[ 8] =  -0.2880213168024011E+00;
    x[ 9] =  -0.1455618541608951E+00;
    x[10] =   0.0000000000000000E+00;
    x[11] =   0.1455618541608951E+00;
    x[12] =   0.2880213168024011E+00;
    x[13] =   0.4243421202074388E+00;
    x[14] =   0.5516188358872198E+00;
    x[15] =   0.6671388041974123E+00;
    x[16] =   0.7684399634756779E+00;
    x[17] =   0.8533633645833173E+00;
    x[18] =   0.9200993341504008E+00;
    x[19] =   0.9672268385663063E+00;
    x[20] =   0.9937521706203896E+00;

    w[ 0] =   0.1601722825777420E-01;
    w[ 1] =   0.3695378977085242E-01;
    w[ 2] =   0.5713442542685715E-01;
    w[ 3] =   0.7610011362837928E-01;
    w[ 4] =   0.9344442345603393E-01;
    w[ 5] =   0.1087972991671484E+00;
    w[ 6] =   0.1218314160537285E+00;
    w[ 7] =   0.1322689386333373E+00;
    w[ 8] =   0.1398873947910731E+00;
    w[ 9] =   0.1445244039899700E+00;
    w[10] =   0.1460811336496904E+00;
    w[11] =   0.1445244039899700E+00;
    w[12] =   0.1398873947910731E+00;
    w[13] =   0.1322689386333373E+00;
    w[14] =   0.1218314160537285E+00;
    w[15] =   0.1087972991671484E+00;
    w[16] =   0.9344442345603393E-01;
    w[17] =   0.7610011362837928E-01;
    w[18] =   0.5713442542685715E-01;
    w[19] =   0.3695378977085242E-01;
    w[20] =   0.1601722825777420E-01;
  }
  else if ( n == 22 )
  {
    x[ 0] =  -0.9942945854823994E+00;
    x[ 1] =  -0.9700604978354287E+00;
    x[ 2] =  -0.9269567721871740E+00;
    x[ 3] =  -0.8658125777203002E+00;
    x[ 4] =  -0.7878168059792081E+00;
    x[ 5] =  -0.6944872631866827E+00;
    x[ 6] =  -0.5876404035069116E+00;
    x[ 7] =  -0.4693558379867570E+00;
    x[ 8] =  -0.3419358208920842E+00;
    x [9] =  -0.2078604266882213E+00;
    x[10] =  -0.6973927331972223E-01;
    x[11] =   0.6973927331972223E-01;
    x[12] =   0.2078604266882213E+00;
    x[13] =   0.3419358208920842E+00;
    x[14] =   0.4693558379867570E+00;
    x[15] =   0.5876404035069116E+00;
    x[16] =   0.6944872631866827E+00;
    x[17] =   0.7878168059792081E+00;
    x[18] =   0.8658125777203002E+00;
    x[19] =   0.9269567721871740E+00;
    x[20] =   0.9700604978354287E+00;
    x[21] =   0.9942945854823994E+00;

    w[ 0] =   0.1462799529827203E-01;
    w[ 1] =   0.3377490158481413E-01;
    w[ 2] =   0.5229333515268327E-01;
    w[ 3] =   0.6979646842452038E-01;
    w[ 4] =   0.8594160621706777E-01;
    w[ 5] =   0.1004141444428809E+00;
    w[ 6] =   0.1129322960805392E+00;
    w[ 7] =   0.1232523768105124E+00;
    w[ 8] =   0.1311735047870623E+00;
    w[ 9] =   0.1365414983460152E+00;
    w[10] =   0.1392518728556321E+00;
    w[11] =   0.1392518728556321E+00;
    w[12] =   0.1365414983460152E+00;
    w[13] =   0.1311735047870623E+00;
    w[14] =   0.1232523768105124E+00;
    w[15] =   0.1129322960805392E+00;
    w[16] =   0.1004141444428809E+00;
    w[17] =   0.8594160621706777E-01;
    w[18] =   0.6979646842452038E-01;
    w[19] =   0.5229333515268327E-01;
    w[20] =   0.3377490158481413E-01;
    w[21] =   0.1462799529827203E-01;
  }
  else if ( n == 23 )
  {
    x[ 0] =  -0.9947693349975522E+00;
    x[ 1] =  -0.9725424712181152E+00;
    x[ 2] =  -0.9329710868260161E+00;
    x[ 3] =  -0.8767523582704416E+00;
    x[ 4] =  -0.8048884016188399E+00;
    x[ 5] =  -0.7186613631319502E+00;
    x[ 6] =  -0.6196098757636461E+00;
    x[ 7] =  -0.5095014778460075E+00;
    x[ 8] =  -0.3903010380302908E+00;
    x[ 9] =  -0.2641356809703449E+00;
    x[10] =  -0.1332568242984661E+00;
    x[11] =   0.0000000000000000E+00;
    x[12] =   0.1332568242984661E+00;
    x[13] =   0.2641356809703449E+00;
    x[14] =   0.3903010380302908E+00;
    x[15] =   0.5095014778460075E+00;
    x[16] =   0.6196098757636461E+00;
    x[17] =   0.7186613631319502E+00;
    x[18] =   0.8048884016188399E+00;
    x[19] =   0.8767523582704416E+00;
    x[20] =   0.9329710868260161E+00;
    x[21] =   0.9725424712181152E+00;
    x[22] =   0.9947693349975522E+00;

    w[ 0] =   0.1341185948714167E-01;
    w[ 1] =   0.3098800585697944E-01;
    w[ 2] =   0.4803767173108464E-01;
    w[ 3] =   0.6423242140852586E-01;
    w[ 4] =   0.7928141177671895E-01;
    w[ 5] =   0.9291576606003514E-01;
    w[ 6] =   0.1048920914645414E+00;
    w[ 7] =   0.1149966402224114E+00;
    w[ 8] =   0.1230490843067295E+00;
    w[ 9] =   0.1289057221880822E+00;
    w[10] =   0.1324620394046967E+00;
    w[11] =   0.1336545721861062E+00;
    w[12] =   0.1324620394046967E+00;
    w[13] =   0.1289057221880822E+00;
    w[14] =   0.1230490843067295E+00;
    w[15] =   0.1149966402224114E+00;
    w[16] =   0.1048920914645414E+00;
    w[17] =   0.9291576606003514E-01;
    w[18] =   0.7928141177671895E-01;
    w[19] =   0.6423242140852586E-01;
    w[20] =   0.4803767173108464E-01;
    w[21] =   0.3098800585697944E-01;
    w[22] =   0.1341185948714167E-01;
  }
  else if ( n == 24 )
  {
    x[ 0] =  -0.9951872199970213E+00;
    x[ 1] =  -0.9747285559713095E+00;
    x[ 2] =  -0.9382745520027327E+00;
    x[ 3] =  -0.8864155270044011E+00;
    x[ 4] =  -0.8200019859739029E+00;
    x[ 5] =  -0.7401241915785544E+00;
    x[ 6] =  -0.6480936519369755E+00;
    x[ 7] =  -0.5454214713888396E+00;
    x[ 8] =  -0.4337935076260451E+00;
    x[ 9] =  -0.3150426796961634E+00;
    x[10] =  -0.1911188674736163E+00;
    x[11] =  -0.6405689286260562E-01;
    x[12] =   0.6405689286260562E-01;
    x[13] =   0.1911188674736163E+00;
    x[14] =   0.3150426796961634E+00;
    x[15] =   0.4337935076260451E+00;
    x[16] =   0.5454214713888396E+00;
    x[17] =   0.6480936519369755E+00;
    x[18] =   0.7401241915785544E+00;
    x[19] =   0.8200019859739029E+00;
    x[20] =   0.8864155270044011E+00;
    x[21] =   0.9382745520027327E+00;
    x[22] =   0.9747285559713095E+00;
    x[23] =   0.9951872199970213E+00;

    w[ 0] =   0.1234122979998730E-01;
    w[ 1] =   0.2853138862893375E-01;
    w[ 2] =   0.4427743881741982E-01;
    w[ 3] =   0.5929858491543672E-01;
    w[ 4] =   0.7334648141108031E-01;
    w[ 5] =   0.8619016153195320E-01;
    w[ 6] =   0.9761865210411380E-01;
    w[ 7] =   0.1074442701159656E+00;
    w[ 8] =   0.1155056680537256E+00;
    w[ 9] =   0.1216704729278035E+00;
    w[10] =   0.1258374563468283E+00;
    w[11] =   0.1279381953467521E+00;
    w[12] =   0.1279381953467521E+00;
    w[13] =   0.1258374563468283E+00;
    w[14] =   0.1216704729278035E+00;
    w[15] =   0.1155056680537256E+00;
    w[16] =   0.1074442701159656E+00;
    w[17] =   0.9761865210411380E-01;
    w[18] =   0.8619016153195320E-01;
    w[19] =   0.7334648141108031E-01;
    w[20] =   0.5929858491543672E-01;
    w[21] =   0.4427743881741982E-01;
    w[22] =   0.2853138862893375E-01;
    w[23] =   0.1234122979998730E-01;
  }
  else if ( n == 25 )
  {
    x[ 0] =  -0.9955569697904981E+00;
    x[ 1] =  -0.9766639214595175E+00;
    x[ 2] =  -0.9429745712289743E+00;
    x[ 3] =  -0.8949919978782754E+00;
    x[ 4] =  -0.8334426287608340E+00;
    x[ 5] =  -0.7592592630373577E+00;
    x[ 6] =  -0.6735663684734684E+00;
    x[ 7] =  -0.5776629302412229E+00;
    x[ 8] =  -0.4730027314457150E+00;
    x[ 9] =  -0.3611723058093879E+00;
    x[10] =  -0.2438668837209884E+00;
    x[11] =  -0.1228646926107104E+00;
    x[12] =   0.0000000000000000E+00;
    x[13] =   0.1228646926107104E+00;
    x[14] =   0.2438668837209884E+00;
    x[15] =   0.3611723058093879E+00;
    x[16] =   0.4730027314457150E+00;
    x[17] =   0.5776629302412229E+00;
    x[18] =   0.6735663684734684E+00;
    x[19] =   0.7592592630373577E+00;
    x[20] =   0.8334426287608340E+00;
    x[21] =   0.8949919978782754E+00;
    x[22] =   0.9429745712289743E+00;
    x[23] =   0.9766639214595175E+00;
    x[24] =   0.9955569697904981E+00;

    w[ 0] =   0.1139379850102617E-01;
    w[ 1] =   0.2635498661503214E-01;
    w[ 2] =   0.4093915670130639E-01;
    w[ 3] =   0.5490469597583517E-01;
    w[ 4] =   0.6803833381235694E-01;
    w[ 5] =   0.8014070033500101E-01;
    w[ 6] =   0.9102826198296370E-01;
    w[ 7] =   0.1005359490670506E+00;
    w[ 8] =   0.1085196244742637E+00;
    w[ 9] =   0.1148582591457116E+00;
    w[10] =   0.1194557635357847E+00;
    w[11] =   0.1222424429903101E+00;
    w[12] =   0.1231760537267154E+00;
    w[13] =   0.1222424429903101E+00;
    w[14] =   0.1194557635357847E+00;
    w[15] =   0.1148582591457116E+00;
    w[16] =   0.1085196244742637E+00;
    w[17] =   0.1005359490670506E+00;
    w[18] =   0.9102826198296370E-01;
    w[19] =   0.8014070033500101E-01;
    w[20] =   0.6803833381235694E-01;
    w[21] =   0.5490469597583517E-01;
    w[22] =   0.4093915670130639E-01;
    w[23] =   0.2635498661503214E-01;
    w[24] =   0.1139379850102617E-01;
  }
  else if ( n == 26 )
  {
    x[ 0] =  -0.9958857011456169E+00;
    x[ 1] =  -0.9783854459564710E+00;
    x[ 2] =  -0.9471590666617142E+00;
    x[ 3] =  -0.9026378619843071E+00;
    x[ 4] =  -0.8454459427884981E+00;
    x[ 5] =  -0.7763859488206789E+00;
    x[ 6] =  -0.6964272604199573E+00;
    x[ 7] =  -0.6066922930176181E+00;
    x[ 8] =  -0.5084407148245057E+00;
    x[ 9] =  -0.4030517551234863E+00;
    x[10] =  -0.2920048394859569E+00;
    x[11] =  -0.1768588203568902E+00;
    x[12] =  -0.5923009342931320E-01;
    x[13] =   0.5923009342931320E-01;
    x[14] =   0.1768588203568902E+00;
    x[15] =   0.2920048394859569E+00;
    x[16] =   0.4030517551234863E+00;
    x[17] =   0.5084407148245057E+00;
    x[18] =   0.6066922930176181E+00;
    x[19] =   0.6964272604199573E+00;
    x[20] =   0.7763859488206789E+00;
    x[21] =   0.8454459427884981E+00;
    x[22] =   0.9026378619843071E+00;
    x[23] =   0.9471590666617142E+00;
    x[24] =   0.9783854459564710E+00;
    x[25] =   0.9958857011456169E+00;

    w[ 0] =   0.1055137261734304E-01;
    w[ 1] =   0.2441785109263173E-01;
    w[ 2] =   0.3796238329436282E-01;
    w[ 3] =   0.5097582529714782E-01;
    w[ 4] =   0.6327404632957484E-01;
    w[ 5] =   0.7468414976565967E-01;
    w[ 6] =   0.8504589431348521E-01;
    w[ 7] =   0.9421380035591416E-01;
    w[ 8] =   0.1020591610944255E+00;
    w[ 9] =   0.1084718405285765E+00;
    w[10] =   0.1133618165463197E+00;
    w[11] =   0.1166604434852967E+00;
    w[12] =   0.1183214152792622E+00;
    w[13] =   0.1183214152792622E+00;
    w[14] =   0.1166604434852967E+00;
    w[15] =   0.1133618165463197E+00;
    w[16] =   0.1084718405285765E+00;
    w[17] =   0.1020591610944255E+00;
    w[18] =   0.9421380035591416E-01;
    w[19] =   0.8504589431348521E-01;
    w[20] =   0.7468414976565967E-01;
    w[21] =   0.6327404632957484E-01;
    w[22] =   0.5097582529714782E-01;
    w[23] =   0.3796238329436282E-01;
    w[24] =   0.2441785109263173E-01;
    w[25] =   0.1055137261734304E-01;
  }
  else if ( n == 27 )
  {
    x[ 0] =  -0.9961792628889886E+00;
    x[ 1] =  -0.9799234759615012E+00;
    x[ 2] =  -0.9509005578147051E+00;
    x[ 3] =  -0.9094823206774911E+00;
    x[ 4] =  -0.8562079080182945E+00;
    x[ 5] =  -0.7917716390705082E+00;
    x[ 6] =  -0.7170134737394237E+00;
    x[ 7] =  -0.6329079719464952E+00;
    x[ 8] =  -0.5405515645794569E+00;
    x[ 9] =  -0.4411482517500269E+00;
    x[10] =  -0.3359939036385089E+00;
    x[11] =  -0.2264593654395369E+00;
    x[12] =  -0.1139725856095300E+00;
    x[13] =   0.0000000000000000E+00;
    x[14] =   0.1139725856095300E+00;
    x[15] =   0.2264593654395369E+00;
    x[16] =   0.3359939036385089E+00;
    x[17] =   0.4411482517500269E+00;
    x[18] =   0.5405515645794569E+00;
    x[19] =   0.6329079719464952E+00;
    x[20] =   0.7170134737394237E+00;
    x[21] =   0.7917716390705082E+00;
    x[22] =   0.8562079080182945E+00;
    x[23] =   0.9094823206774911E+00;
    x[24] =   0.9509005578147051E+00;
    x[25] =   0.9799234759615012E+00;
    x[26] =   0.9961792628889886E+00;

    w[ 0] =   0.9798996051294232E-02;
    w[ 1] =   0.2268623159618062E-01;
    w[ 2] =   0.3529705375741969E-01;
    w[ 3] =   0.4744941252061504E-01;
    w[ 4] =   0.5898353685983366E-01;
    w[ 5] =   0.6974882376624561E-01;
    w[ 6] =   0.7960486777305781E-01;
    w[ 7] =   0.8842315854375689E-01;
    w[ 8] =   0.9608872737002842E-01;
    w[ 9] =   0.1025016378177459E+00;
    w[10] =   0.1075782857885332E+00;
    w[11] =   0.1112524883568452E+00;
    w[12] =   0.1134763461089651E+00;
    w[13] =   0.1142208673789570E+00;
    w[14] =   0.1134763461089651E+00;
    w[15] =   0.1112524883568452E+00;
    w[16] =   0.1075782857885332E+00;
    w[17] =   0.1025016378177459E+00;
    w[18] =   0.9608872737002842E-01;
    w[19] =   0.8842315854375689E-01;
    w[20] =   0.7960486777305781E-01;
    w[21] =   0.6974882376624561E-01;
    w[22] =   0.5898353685983366E-01;
    w[23] =   0.4744941252061504E-01;
    w[24] =   0.3529705375741969E-01;
    w[25] =   0.2268623159618062E-01;
    w[26] =   0.9798996051294232E-02;
  }
  else if ( n == 28 )
  {
    x[ 0] =  -0.9964424975739544E+00;
    x[ 1] =  -0.9813031653708728E+00;
    x[ 2] =  -0.9542592806289382E+00;
    x[ 3] =  -0.9156330263921321E+00;
    x[ 4] =  -0.8658925225743951E+00;
    x[ 5] =  -0.8056413709171791E+00;
    x[ 6] =  -0.7356108780136318E+00;
    x[ 7] =  -0.6566510940388650E+00;
    x[ 8] =  -0.5697204718114017E+00;
    x[ 9] =  -0.4758742249551183E+00;
    x[10] =  -0.3762515160890787E+00;
    x[11] =  -0.2720616276351780E+00;
    x[12] =  -0.1645692821333808E+00;
    x[13] =  -0.5507928988403427E-01;
    x[14] =   0.5507928988403427E-01;
    x[15] =   0.1645692821333808E+00;
    x[16] =   0.2720616276351780E+00;
    x[17] =   0.3762515160890787E+00;
    x[18] =   0.4758742249551183E+00;
    x[19] =   0.5697204718114017E+00;
    x[20] =   0.6566510940388650E+00;
    x[21] =   0.7356108780136318E+00;
    x[22] =   0.8056413709171791E+00;
    x[23] =   0.8658925225743951E+00;
    x[24] =   0.9156330263921321E+00;
    x[25] =   0.9542592806289382E+00;
    x[26] =   0.9813031653708728E+00;
    x[27] =   0.9964424975739544E+00;

    w[ 0] =   0.9124282593094672E-02;
    w[ 1] =   0.2113211259277118E-01;
    w[ 2] =   0.3290142778230441E-01;
    w[ 3] =   0.4427293475900429E-01;
    w[ 4] =   0.5510734567571667E-01;
    w[ 5] =   0.6527292396699959E-01;
    w[ 6] =   0.7464621423456877E-01;
    w[ 7] =   0.8311341722890127E-01;
    w[ 8] =   0.9057174439303289E-01;
    w[ 9] =   0.9693065799792999E-01;
    w[10] =   0.1021129675780608E+00;
    w[11] =   0.1060557659228464E+00;
    w[12] =   0.1087111922582942E+00;
    w[13] =   0.1100470130164752E+00;
    w[14] =   0.1100470130164752E+00;
    w[15] =   0.1087111922582942E+00;
    w[16] =   0.1060557659228464E+00;
    w[17] =   0.1021129675780608E+00;
    w[18] =   0.9693065799792999E-01;
    w[19] =   0.9057174439303289E-01;
    w[20] =   0.8311341722890127E-01;
    w[21] =   0.7464621423456877E-01;
    w[22] =   0.6527292396699959E-01;
    w[23] =   0.5510734567571667E-01;
    w[24] =   0.4427293475900429E-01;
    w[25] =   0.3290142778230441E-01;
    w[26] =   0.2113211259277118E-01;
    w[27] =   0.9124282593094672E-02;
  }
  else if ( n == 29 )
  {
    x[ 0] =  -0.9966794422605966E+00;
    x[ 1] =  -0.9825455052614132E+00;
    x[ 2] =  -0.9572855957780877E+00;
    x[ 3] =  -0.9211802329530588E+00;
    x[ 4] =  -0.8746378049201028E+00;
    x[ 5] =  -0.8181854876152524E+00;
    x[ 6] =  -0.7524628517344771E+00;
    x[ 7] =  -0.6782145376026865E+00;
    x[ 8] =  -0.5962817971382278E+00;
    x[ 9] =  -0.5075929551242276E+00;
    x[10] =  -0.4131528881740087E+00;
    x[11] =  -0.3140316378676399E+00;
    x[12] =  -0.2113522861660011E+00;
    x[13] =  -0.1062782301326792E+00;
    x[14] =   0.0000000000000000E+00;
    x[15] =   0.1062782301326792E+00;
    x[16] =   0.2113522861660011E+00;
    x[17] =   0.3140316378676399E+00;
    x[18] =   0.4131528881740087E+00;
    x[19] =   0.5075929551242276E+00;
    x[20] =   0.5962817971382278E+00;
    x[21] =   0.6782145376026865E+00;
    x[22] =   0.7524628517344771E+00;
    x[23] =   0.8181854876152524E+00;
    x[24] =   0.8746378049201028E+00;
    x[25] =   0.9211802329530588E+00;
    x[26] =   0.9572855957780877E+00;
    x[27] =   0.9825455052614132E+00;
    x[28] =   0.9966794422605966E+00;

    w[ 0] =   0.8516903878746365E-02;
    w[ 1] =   0.1973208505612276E-01;
    w[ 2] =   0.3074049220209360E-01;
    w[ 3] =   0.4140206251868281E-01;
    w[ 4] =   0.5159482690249799E-01;
    w[ 5] =   0.6120309065707916E-01;
    w[ 6] =   0.7011793325505125E-01;
    w[ 7] =   0.7823832713576385E-01;
    w[ 8] =   0.8547225736617248E-01;
    w[ 9] =   0.9173775713925882E-01;
    w[10] =   0.9696383409440862E-01;
    w[11] =   0.1010912737599150E+00;
    w[12] =   0.1040733100777293E+00;
    w[13] =   0.1058761550973210E+00;
    w[14] =   0.1064793817183143E+00;
    w[15] =   0.1058761550973210E+00;
    w[16] =   0.1040733100777293E+00;
    w[17] =   0.1010912737599150E+00;
    w[18] =   0.9696383409440862E-01;
    w[19] =   0.9173775713925882E-01;
    w[20] =   0.8547225736617248E-01;
    w[21] =   0.7823832713576385E-01;
    w[22] =   0.7011793325505125E-01;
    w[23] =   0.6120309065707916E-01;
    w[24] =   0.5159482690249799E-01;
    w[25] =   0.4140206251868281E-01;
    w[26] =   0.3074049220209360E-01;
    w[27] =   0.1973208505612276E-01;
    w[28] =   0.8516903878746365E-02;
  }
  else if ( n == 30 )
  {
    x[ 0] =  -0.9968934840746495E+00;
    x[ 1] =  -0.9836681232797472E+00;
    x[ 2] =  -0.9600218649683075E+00;
    x[ 3] =  -0.9262000474292743E+00;
    x[ 4] =  -0.8825605357920526E+00;
    x[ 5] =  -0.8295657623827684E+00;
    x[ 6] =  -0.7677774321048262E+00;
    x[ 7] =  -0.6978504947933158E+00;
    x[ 8] =  -0.6205261829892429E+00;
    x[ 9] =  -0.5366241481420199E+00;
    x[10] =  -0.4470337695380892E+00;
    x[11] =  -0.3527047255308781E+00;
    x[12] =  -0.2546369261678899E+00;
    x[13] =  -0.1538699136085835E+00;
    x[14] =  -0.5147184255531770E-01;
    x[15] =   0.5147184255531770E-01;
    x[16] =   0.1538699136085835E+00;
    x[17] =   0.2546369261678899E+00;
    x[18] =   0.3527047255308781E+00;
    x[19] =   0.4470337695380892E+00;
    x[20] =   0.5366241481420199E+00;
    x[21] =   0.6205261829892429E+00;
    x[22] =   0.6978504947933158E+00;
    x[23] =   0.7677774321048262E+00;
    x[24] =   0.8295657623827684E+00;
    x[25] =   0.8825605357920526E+00;
    x[26] =   0.9262000474292743E+00;
    x[27] =   0.9600218649683075E+00;
    x[28] =   0.9836681232797472E+00;
    x[29] =   0.9968934840746495E+00;

    w[ 0] =   0.7968192496166648E-02;
    w[ 1] =   0.1846646831109099E-01;
    w[ 2] =   0.2878470788332330E-01;
    w[ 3] =   0.3879919256962704E-01;
    w[ 4] =   0.4840267283059405E-01;
    w[ 5] =   0.5749315621761905E-01;
    w[ 6] =   0.6597422988218052E-01;
    w[ 7] =   0.7375597473770516E-01;
    w[ 8] =   0.8075589522942023E-01;
    w[ 9] =   0.8689978720108314E-01;
    w[10] =   0.9212252223778619E-01;
    w[11] =   0.9636873717464424E-01;
    w[12] =   0.9959342058679524E-01;
    w[13] =   0.1017623897484056E+00;
    w[14] =   0.1028526528935587E+00;
    w[15] =   0.1028526528935587E+00;
    w[16] =   0.1017623897484056E+00;
    w[17] =   0.9959342058679524E-01;
    w[18] =   0.9636873717464424E-01;
    w[19] =   0.9212252223778619E-01;
    w[20] =   0.8689978720108314E-01;
    w[21] =   0.8075589522942023E-01;
    w[22] =   0.7375597473770516E-01;
    w[23] =   0.6597422988218052E-01;
    w[24] =   0.5749315621761905E-01;
    w[25] =   0.4840267283059405E-01;
    w[26] =   0.3879919256962704E-01;
    w[27] =   0.2878470788332330E-01;
    w[28] =   0.1846646831109099E-01;
    w[29] =   0.7968192496166648E-02;
  }
  else if ( n == 31 )
  {
    x[ 0] =  -0.99708748181947707454263838179654;    
    x[ 1] =  -0.98468590966515248400211329970113;    
    x[ 2] =  -0.96250392509294966178905249675943;    
    x[ 3] =  -0.93075699789664816495694576311725;    
    x[ 4] =  -0.88976002994827104337419200908023;    
    x[ 5] =  -0.83992032014626734008690453594388;    
    x[ 6] =  -0.78173314841662494040636002019484;    
    x[ 7] =  -0.71577678458685328390597086536649;    
    x[ 8] =  -0.64270672292426034618441820323250;    
    x[ 9] =  -0.56324916140714926272094492359516;    
    x[10] =  -0.47819378204490248044059403935649;    
    x[11] =  -0.38838590160823294306135146128752;    
    x[12] =  -0.29471806998170161661790389767170;    
    x[13] =  -0.19812119933557062877241299603283;    
    x[14] =  -0.99555312152341520325174790118941E-01;
    x[15] =   0.00000000000000000000000000000000;   
    x[16] =   0.99555312152341520325174790118941E-01;
    x[17] =   0.19812119933557062877241299603283;    
    x[18] =   0.29471806998170161661790389767170;    
    x[19] =   0.38838590160823294306135146128752;    
    x[20] =   0.47819378204490248044059403935649;    
    x[21] =   0.56324916140714926272094492359516;    
    x[22] =   0.64270672292426034618441820323250;    
    x[23] =   0.71577678458685328390597086536649;    
    x[24] =   0.78173314841662494040636002019484;    
    x[25] =   0.83992032014626734008690453594388;    
    x[26] =   0.88976002994827104337419200908023;    
    x[27] =   0.93075699789664816495694576311725;    
    x[28] =   0.96250392509294966178905249675943;    
    x[29] =   0.98468590966515248400211329970113;    
    x[30] =   0.99708748181947707454263838179654;
 
    w[ 0] =   0.74708315792487746093913218970494E-02;
    w[ 1] =   0.17318620790310582463552990782414E-01;
    w[ 2] =   0.27009019184979421800608642617676E-01;
    w[ 3] =   0.36432273912385464024392008749009E-01;
    w[ 4] =   0.45493707527201102902315857856518E-01;
    w[ 5] =   0.54103082424916853711666259085477E-01;
    w[ 6] =   0.62174786561028426910343543686657E-01;
    w[ 7] =   0.69628583235410366167756126255124E-01;
    w[ 8] =   0.76390386598776616426357674901331E-01;
    w[ 9] =   0.82392991761589263903823367431962E-01;
    w[10] =   0.87576740608477876126198069695333E-01;
    w[11] =   0.91890113893641478215362871607150E-01;
    w[12] =   0.95290242912319512807204197487597E-01;
    w[13] =   0.97743335386328725093474010978997E-01;
    w[14] =   0.99225011226672307874875514428615E-01;
    w[15] =   0.99720544793426451427533833734349E-01;
    w[16] =   0.99225011226672307874875514428615E-01;
    w[17] =   0.97743335386328725093474010978997E-01;
    w[18] =   0.95290242912319512807204197487597E-01;
    w[19] =   0.91890113893641478215362871607150E-01;
    w[20] =   0.87576740608477876126198069695333E-01;
    w[21] =   0.82392991761589263903823367431962E-01;
    w[22] =   0.76390386598776616426357674901331E-01;
    w[23] =   0.69628583235410366167756126255124E-01;
    w[24] =   0.62174786561028426910343543686657E-01;
    w[25] =   0.54103082424916853711666259085477E-01;
    w[26] =   0.45493707527201102902315857856518E-01;
    w[27] =   0.36432273912385464024392008749009E-01;
    w[28] =   0.27009019184979421800608642617676E-01;
    w[29] =   0.17318620790310582463552990782414E-01;
    w[30] =   0.74708315792487746093913218970494E-02;
  }
  else if ( n == 32 )
  {
    x[0] =  - 0.997263861849481563544981128665;
    x[1] =  - 0.985611511545268335400175044631;
    x[2] =  - 0.964762255587506430773811928118;
    x[3] =  - 0.934906075937739689170919134835;
    x[4] =  - 0.896321155766052123965307243719;
    x[5] =  - 0.849367613732569970133693004968;
    x[6] =  - 0.794483795967942406963097298970;
    x[7] =  - 0.732182118740289680387426665091;
    x[8] =  - 0.663044266930215200975115168663;
    x[9] =  - 0.587715757240762329040745476402;
    x[10] = - 0.506899908932229390023747474378;
    x[11] = - 0.421351276130635345364119436172;
    x[12] = - 0.331868602282127649779916805730;
    x[13] = - 0.239287362252137074544603209166;
    x[14] = - 0.144471961582796493485186373599;
    x[15] = - 0.483076656877383162348125704405E-01;
    x[16] =   0.483076656877383162348125704405E-01;
    x[17] =   0.144471961582796493485186373599;
    x[18] =   0.239287362252137074544603209166;
    x[19] =   0.331868602282127649779916805730;
    x[20] =   0.421351276130635345364119436172;
    x[21] =   0.506899908932229390023747474378;
    x[22] =   0.587715757240762329040745476402;
    x[23] =   0.663044266930215200975115168663;
    x[24] =   0.732182118740289680387426665091;
    x[25] =   0.794483795967942406963097298970;
    x[26] =   0.849367613732569970133693004968;
    x[27] =   0.896321155766052123965307243719;
    x[28] =   0.934906075937739689170919134835;
    x[29] =   0.964762255587506430773811928118;
    x[30] =   0.985611511545268335400175044631;
    x[31] =   0.997263861849481563544981128665;

    w[0] =  0.701861000947009660040706373885E-02;
    w[1] =  0.162743947309056706051705622064E-01;
    w[2] =  0.253920653092620594557525897892E-01;
    w[3] =  0.342738629130214331026877322524E-01;
    w[4] =  0.428358980222266806568786466061E-01;
    w[5] =  0.509980592623761761961632446895E-01;
    w[6] =  0.586840934785355471452836373002E-01;
    w[7] =  0.658222227763618468376500637069E-01;
    w[8] =  0.723457941088485062253993564785E-01;
    w[9] =  0.781938957870703064717409188283E-01;
    w[10] = 0.833119242269467552221990746043E-01;
    w[11] = 0.876520930044038111427714627518E-01;
    w[12] = 0.911738786957638847128685771116E-01;
    w[13] = 0.938443990808045656391802376681E-01;
    w[14] = 0.956387200792748594190820022041E-01;
    w[15] = 0.965400885147278005667648300636E-01;
    w[16] = 0.965400885147278005667648300636E-01;
    w[17] = 0.956387200792748594190820022041E-01;
    w[18] = 0.938443990808045656391802376681E-01;
    w[19] = 0.911738786957638847128685771116E-01;
    w[20] = 0.876520930044038111427714627518E-01;
    w[21] = 0.833119242269467552221990746043E-01;
    w[22] = 0.781938957870703064717409188283E-01;
    w[23] = 0.723457941088485062253993564785E-01;
    w[24] = 0.658222227763618468376500637069E-01;
    w[25] = 0.586840934785355471452836373002E-01;
    w[26] = 0.509980592623761761961632446895E-01;
    w[27] = 0.428358980222266806568786466061E-01;
    w[28] = 0.342738629130214331026877322524E-01;
    w[29] = 0.253920653092620594557525897892E-01;
    w[30] = 0.162743947309056706051705622064E-01;
    w[31] = 0.701861000947009660040706373885E-02;
  }
  else if ( n == 33 )
  {
    x[ 0] =  -0.9974246942464552;    
    x[ 1] =  -0.9864557262306425;
    x[ 2] =  -0.9668229096899927;
    x[ 3] =  -0.9386943726111684;    
    x[ 4] =  -0.9023167677434336;    
    x[ 5] =  -0.8580096526765041;    
    x[ 6] =  -0.8061623562741665;    
    x[ 7] =  -0.7472304964495622;    
    x[ 8] =  -0.6817319599697428;    
    x[ 9] =  -0.6102423458363790;    
    x[10] =  -0.5333899047863476;    
    x[11] =  -0.4518500172724507;    
    x[12] =  -0.3663392577480734;    
    x[13] =  -0.2776090971524970;    
    x[14] =  -0.1864392988279916;    
    x[15] =  -0.09363106585473338;
    x[16] =   0.000000000000000;
    x[17] =   0.09363106585473338;
    x[18] =   0.1864392988279916;    
    x[19] =   0.2776090971524970;    
    x[20] =   0.3663392577480734;    
    x[21] =   0.4518500172724507;    
    x[22] =   0.5333899047863476;    
    x[23] =   0.6102423458363790;    
    x[24] =   0.6817319599697428;    
    x[25] =   0.7472304964495622;    
    x[26] =   0.8061623562741665;    
    x[27] =   0.8580096526765041;    
    x[28] =   0.9023167677434336;    
    x[29] =   0.9386943726111684;    
    x[30] =   0.9668229096899927;    
    x[31] =   0.9864557262306425;    
    x[32] =   0.9974246942464552;    
 
    w[ 0] =   0.6606227847587558E-02;
    w[ 1] =   0.1532170151293465E-01;
    w[ 2] =   0.2391554810174960E-01;
    w[ 3] =   0.3230035863232891E-01;
    w[ 4] =   0.4040154133166965E-01;
    w[ 5] =   0.4814774281871162E-01;
    w[ 6] =   0.5547084663166357E-01;
    w[ 7] =   0.6230648253031755E-01;
    w[ 8] =   0.6859457281865676E-01;
    w[ 9] =   0.7427985484395420E-01;
    w[10] =   0.7931236479488685E-01;
    w[11] =   0.8364787606703869E-01;
    w[12] =   0.8724828761884425E-01;
    w[13] =   0.9008195866063859E-01;
    w[14] =   0.9212398664331678E-01;
    w[15] =   0.9335642606559612E-01;
    w[16] =   0.9376844616020999E-01;
    w[17] =   0.9335642606559612E-01;
    w[18] =   0.9212398664331678E-01;
    w[19] =   0.9008195866063859E-01;
    w[20] =   0.8724828761884425E-01;
    w[21] =   0.8364787606703869E-01;
    w[22] =   0.7931236479488685E-01;
    w[23] =   0.7427985484395420E-01;
    w[24] =   0.6859457281865676E-01;
    w[25] =   0.6230648253031755E-01;
    w[26] =   0.5547084663166357E-01;
    w[27] =   0.4814774281871162E-01;
    w[28] =   0.4040154133166965E-01;
    w[29] =   0.3230035863232891E-01;
    w[30] =   0.2391554810174960E-01;
    w[31] =   0.1532170151293465E-01;
    w[32] =   0.6606227847587558E-02;
  }
  else if ( n == 64 )
  {
    x[0] =  - 0.999305041735772139456905624346;
    x[1] =  - 0.996340116771955279346924500676;
    x[2] =  - 0.991013371476744320739382383443;
    x[3] =  - 0.983336253884625956931299302157;
    x[4] =  - 0.973326827789910963741853507352;
    x[5] =  - 0.961008799652053718918614121897;
    x[6] =  - 0.946411374858402816062481491347;
    x[7] =  - 0.929569172131939575821490154559;
    x[8] =  - 0.910522137078502805756380668008;
    x[9] =  - 0.889315445995114105853404038273;
    x[10] = - 0.865999398154092819760783385070;
    x[11] = - 0.840629296252580362751691544696;
    x[12] = - 0.813265315122797559741923338086;
    x[13] = - 0.783972358943341407610220525214;
    x[14] = - 0.752819907260531896611863774886;
    x[15] = - 0.719881850171610826848940217832;
    x[16] = - 0.685236313054233242563558371031;
    x[17] = - 0.648965471254657339857761231993;
    x[18] = - 0.611155355172393250248852971019;
    x[19] = - 0.571895646202634034283878116659;
    x[20] = - 0.531279464019894545658013903544;
    x[21] = - 0.489403145707052957478526307022;
    x[22] = - 0.446366017253464087984947714759;
    x[23] = - 0.402270157963991603695766771260;
    x[24] = - 0.357220158337668115950442615046;
    x[25] = - 0.311322871990210956157512698560;
    x[26] = - 0.264687162208767416373964172510;
    x[27] = - 0.217423643740007084149648748989;
    x[28] = - 0.169644420423992818037313629748;
    x[29] = - 0.121462819296120554470376463492;
    x[30] = - 0.729931217877990394495429419403E-01;
    x[31] = - 0.243502926634244325089558428537E-01;
    x[32] =   0.243502926634244325089558428537E-01;
    x[33] =   0.729931217877990394495429419403E-01;
    x[34] =   0.121462819296120554470376463492;
    x[35] =   0.169644420423992818037313629748;
    x[36] =   0.217423643740007084149648748989;
    x[37] =   0.264687162208767416373964172510;
    x[38] =   0.311322871990210956157512698560;
    x[39] =   0.357220158337668115950442615046;
    x[40] =   0.402270157963991603695766771260;
    x[41] =   0.446366017253464087984947714759;
    x[42] =   0.489403145707052957478526307022;
    x[43] =   0.531279464019894545658013903544;
    x[44] =   0.571895646202634034283878116659;
    x[45] =   0.611155355172393250248852971019;
    x[46] =   0.648965471254657339857761231993;
    x[47] =   0.685236313054233242563558371031;
    x[48] =   0.719881850171610826848940217832;
    x[49] =   0.752819907260531896611863774886;
    x[50] =   0.783972358943341407610220525214;
    x[51] =   0.813265315122797559741923338086;
    x[52] =   0.840629296252580362751691544696;
    x[53] =   0.865999398154092819760783385070;
    x[54] =   0.889315445995114105853404038273;
    x[55] =   0.910522137078502805756380668008;
    x[56] =   0.929569172131939575821490154559;
    x[57] =   0.946411374858402816062481491347;
    x[58] =   0.961008799652053718918614121897;
    x[59] =   0.973326827789910963741853507352;
    x[60] =   0.983336253884625956931299302157;
    x[61] =   0.991013371476744320739382383443;
    x[62] =   0.996340116771955279346924500676;
    x[63] =   0.999305041735772139456905624346;

    w[0] =  0.178328072169643294729607914497E-02;
    w[1] =  0.414703326056246763528753572855E-02;
    w[2] =  0.650445796897836285611736039998E-02;
    w[3] =  0.884675982636394772303091465973E-02;
    w[4] =  0.111681394601311288185904930192E-01;
    w[5] =  0.134630478967186425980607666860E-01;
    w[6] =  0.157260304760247193219659952975E-01;
    w[7] =  0.179517157756973430850453020011E-01;
    w[8] =  0.201348231535302093723403167285E-01;
    w[9] =  0.222701738083832541592983303842E-01;
    w[10] = 0.243527025687108733381775504091E-01;
    w[11] = 0.263774697150546586716917926252E-01;
    w[12] = 0.283396726142594832275113052002E-01;
    w[13] = 0.302346570724024788679740598195E-01;
    w[14] = 0.320579283548515535854675043479E-01;
    w[15] = 0.338051618371416093915654821107E-01;
    w[16] = 0.354722132568823838106931467152E-01;
    w[17] = 0.370551285402400460404151018096E-01;
    w[18] = 0.385501531786156291289624969468E-01;
    w[19] = 0.399537411327203413866569261283E-01;
    w[20] = 0.412625632426235286101562974736E-01;
    w[21] = 0.424735151236535890073397679088E-01;
    w[22] = 0.435837245293234533768278609737E-01;
    w[23] = 0.445905581637565630601347100309E-01;
    w[24] = 0.454916279274181444797709969713E-01;
    w[25] = 0.462847965813144172959532492323E-01;
    w[26] = 0.469681828162100173253262857546E-01;
    w[27] = 0.475401657148303086622822069442E-01;
    w[28] = 0.479993885964583077281261798713E-01;
    w[29] = 0.483447622348029571697695271580E-01;
    w[30] = 0.485754674415034269347990667840E-01;
    w[31] = 0.486909570091397203833653907347E-01;
    w[32] = 0.486909570091397203833653907347E-01;
    w[33] = 0.485754674415034269347990667840E-01;
    w[34] = 0.483447622348029571697695271580E-01;
    w[35] = 0.479993885964583077281261798713E-01;
    w[36] = 0.475401657148303086622822069442E-01;
    w[37] = 0.469681828162100173253262857546E-01;
    w[38] = 0.462847965813144172959532492323E-01;
    w[39] = 0.454916279274181444797709969713E-01;
    w[40] = 0.445905581637565630601347100309E-01;
    w[41] = 0.435837245293234533768278609737E-01;
    w[42] = 0.424735151236535890073397679088E-01;
    w[43] = 0.412625632426235286101562974736E-01;
    w[44] = 0.399537411327203413866569261283E-01;
    w[45] = 0.385501531786156291289624969468E-01;
    w[46] = 0.370551285402400460404151018096E-01;
    w[47] = 0.354722132568823838106931467152E-01;
    w[48] = 0.338051618371416093915654821107E-01;
    w[49] = 0.320579283548515535854675043479E-01;
    w[50] = 0.302346570724024788679740598195E-01;
    w[51] = 0.283396726142594832275113052002E-01;
    w[52] = 0.263774697150546586716917926252E-01;
    w[53] = 0.243527025687108733381775504091E-01;
    w[54] = 0.222701738083832541592983303842E-01;
    w[55] = 0.201348231535302093723403167285E-01;
    w[56] = 0.179517157756973430850453020011E-01;
    w[57] = 0.157260304760247193219659952975E-01;
    w[58] = 0.134630478967186425980607666860E-01;
    w[59] = 0.111681394601311288185904930192E-01;
    w[60] = 0.884675982636394772303091465973E-02;
    w[61] = 0.650445796897836285611736039998E-02;
    w[62] = 0.414703326056246763528753572855E-02;
    w[63] = 0.178328072169643294729607914497E-02;
  }
  else if ( n == 65 )
  {
    x[ 0] =  -0.9993260970754129;    
    x[ 1] =  -0.9964509480618492;    
    x[ 2] =  -0.9912852761768016;    
    x[ 3] =  -0.9838398121870350;    
    x[ 4] =  -0.9741315398335512;    
    x[ 5] =  -0.9621827547180553;    
    x[ 6] =  -0.9480209281684076;    
    x[ 7] =  -0.9316786282287494;    
    x[ 8] =  -0.9131934405428462;    
    x[ 9] =  -0.8926078805047389;    
    x[10] =  -0.8699692949264071;    
    x[11] =  -0.8453297528999303;    
    x[12] =  -0.8187459259226514;    
    x[13] =  -0.7902789574921218;    
    x[14] =  -0.7599943224419998;    
    x[15] =  -0.7279616763294247;    
    x[16] =  -0.6942546952139916;    
    x[17] =  -0.6589509061936252;    
    x[18] =  -0.6221315090854003;    
    x[19] =  -0.5838811896604873;    
    x[20] =  -0.5442879248622271;    
    x[21] =  -0.5034427804550069;    
    x[22] =  -0.4614397015691450;    
    x[23] =  -0.4183752966234090;    
    x[24] =  -0.3743486151220660;    
    x[25] =  -0.3294609198374864;    
    x[26] =  -0.2838154539022487;    
    x[27] =  -0.2375172033464168;    
    x[28] =  -0.1906726556261428;    
    x[29] =  -0.1433895546989752;    
    x[30] =  -0.9577665320919751E-01;
    x[31] =  -0.4794346235317186E-01;
    x[32] =    0.000000000000000;    
    x[33] =   0.4794346235317186E-01;
    x[34] =   0.9577665320919751E-01;
    x[35] =   0.1433895546989752;    
    x[36] =   0.1906726556261428;    
    x[37] =   0.2375172033464168;    
    x[38] =   0.2838154539022487;    
    x[39] =   0.3294609198374864;    
    x[40] =   0.3743486151220660;    
    x[41] =   0.4183752966234090;    
    x[42] =   0.4614397015691450;    
    x[43] =   0.5034427804550069;    
    x[44] =   0.5442879248622271;    
    x[45] =   0.5838811896604873;    
    x[46] =   0.6221315090854003;    
    x[47] =   0.6589509061936252;    
    x[48] =   0.6942546952139916;    
    x[49] =   0.7279616763294247;    
    x[50] =   0.7599943224419998;    
    x[51] =   0.7902789574921218;    
    x[52] =   0.8187459259226514;    
    x[53] =   0.8453297528999303;    
    x[54] =   0.8699692949264071;    
    x[55] =   0.8926078805047389;    
    x[56] =   0.9131934405428462;    
    x[57] =   0.9316786282287494;    
    x[58] =   0.9480209281684076;    
    x[59] =   0.9621827547180553;    
    x[60] =   0.9741315398335512;    
    x[61] =   0.9838398121870350;    
    x[62] =   0.9912852761768016;    
    x[63] =   0.9964509480618492;    
    x[64] =   0.9993260970754129;    
 
    w[ 0] =   0.1729258251300218E-02;
    w[ 1] =   0.4021524172003703E-02;
    w[ 2] =   0.6307942578971821E-02;
    w[ 3] =   0.8580148266881443E-02;
    w[ 4] =   0.1083267878959798E-01;
    w[ 5] =   0.1306031163999490E-01;
    w[ 6] =   0.1525791214644825E-01;
    w[ 7] =   0.1742042199767025E-01;
    w[ 8] =   0.1954286583675005E-01;
    w[ 9] =   0.2162036128493408E-01;
    w[10] =   0.2364812969128723E-01;
    w[11] =   0.2562150693803776E-01;
    w[12] =   0.2753595408845034E-01;
    w[13] =   0.2938706778931066E-01;
    w[14] =   0.3117059038018911E-01;
    w[15] =   0.3288241967636860E-01;
    w[16] =   0.3451861839854901E-01;
    w[17] =   0.3607542322556527E-01;
    w[18] =   0.3754925344825770E-01;
    w[19] =   0.3893671920405121E-01;
    w[20] =   0.4023462927300549E-01;
    w[21] =   0.4143999841724028E-01;
    w[22] =   0.4255005424675579E-01;
    w[23] =   0.4356224359580051E-01;
    w[24] =   0.4447423839508296E-01;
    w[25] =   0.4528394102630023E-01;
    w[26] =   0.4598948914665173E-01;
    w[27] =   0.4658925997223349E-01;
    w[28] =   0.4708187401045461E-01;
    w[29] =   0.4746619823288551E-01;
    w[30] =   0.4774134868124067E-01;
    w[31] =   0.4790669250049590E-01;
    w[32] =   0.4796184939446662E-01;
    w[33] =   0.4790669250049590E-01;
    w[34] =   0.4774134868124067E-01;
    w[35] =   0.4746619823288551E-01;
    w[36] =   0.4708187401045461E-01;
    w[37] =   0.4658925997223349E-01;
    w[38] =   0.4598948914665173E-01;
    w[39] =   0.4528394102630023E-01;
    w[40] =   0.4447423839508296E-01;
    w[41] =   0.4356224359580051E-01;
    w[42] =   0.4255005424675579E-01;
    w[43] =   0.4143999841724028E-01;
    w[44] =   0.4023462927300549E-01;
    w[45] =   0.3893671920405121E-01;
    w[46] =   0.3754925344825770E-01;
    w[47] =   0.3607542322556527E-01;
    w[48] =   0.3451861839854901E-01;
    w[49] =   0.3288241967636860E-01;
    w[50] =   0.3117059038018911E-01;
    w[51] =   0.2938706778931066E-01;
    w[52] =   0.2753595408845034E-01;
    w[53] =   0.2562150693803776E-01;
    w[54] =   0.2364812969128723E-01;
    w[55] =   0.2162036128493408E-01;
    w[56] =   0.1954286583675005E-01;
    w[57] =   0.1742042199767025E-01;
    w[58] =   0.1525791214644825E-01;
    w[59] =   0.1306031163999490E-01;
    w[60] =   0.1083267878959798E-01;
    w[61] =   0.8580148266881443E-02;
    w[62] =   0.6307942578971821E-02;
    w[63] =   0.4021524172003703E-02;
    w[64] =   0.1729258251300218E-02;
  }
  else if ( n == 127 ) 
  {
    x[  0] =  -0.99982213041530614629963254927125E+00;
    x[  1] =  -0.99906293435531189513828920479421E+00;    
    x[  2] =  -0.99769756618980462107441703193392E+00;    
    x[  3] =  -0.99572655135202722663543337085008E+00;    
    x[  4] =  -0.99315104925451714736113079489080E+00;  
    x[  5] =  -0.98997261459148415760778669967548E+00;   
    x[  6] =  -0.98619317401693166671043833175407E+00;    
    x[  7] =  -0.98181502080381411003346312451200E+00;    
    x[  8] =  -0.97684081234307032681744391886221E+00;    
    x[  9] =  -0.97127356816152919228894689830512E+00;    
    x[ 10] =  -0.96511666794529212109082507703391E+00;    
    x[ 11] =  -0.95837384942523877114910286998060E+00;    
    x[ 12] =  -0.95104920607788031054790764659636E+00;   
    x[ 13] =  -0.94314718462481482734544963026201E+00;    
    x[ 14] =  -0.93467258232473796857363487794906E+00;    
    x[ 15] =  -0.92563054405623384912746466814259E+00;    
    x[ 16] =  -0.91602655919146580931308861741716E+00;   
    x[ 17] =  -0.90586645826182138280246131760282E+00;    
    x[ 18] =  -0.89515640941708370896904382642451E+00;   
    x[ 19] =  -0.88390291468002656994525794802849E+00;    
    x[ 20] =  -0.87211280599856071141963753428864E+00;    
    x[ 21] =  -0.85979324109774080981203134414483E+00;   
    x[ 22] =  -0.84695169913409759845333931085437E+00;    
    x[ 23] =  -0.83359597615489951437955716480123E+00;    
    x[ 24] =  -0.81973418036507867415511910167470E+00;   
    x[ 25] =  -0.80537472720468021466656079404644E+00;   
    x[ 26] =  -0.79052633423981379994544995252740E+00;   
    x[ 27] =  -0.77519801587020238244496276354566E+00;  
    x[ 28] =  -0.75939907785653667155666366659810E+00;   
    x[ 29] =  -0.74313911167095451292056688997595E+00;   
    x[ 30] =  -0.72642798867407268553569290153270E+00;    
    x[ 31] =  -0.70927585412210456099944463906757E+00;   
    x[ 32] =  -0.69169312100770067015644143286666E+00; 
    x[ 33] =  -0.67369046373825048534668253831602E+00;
    x[ 34] =  -0.65527881165548263027676505156852E+00;
    x[ 35] =  -0.63646934240029724134760815684175E+00;
    x[ 36] =  -0.61727347512685828385763916340822E+00; 
    x[ 37] =  -0.59770286357006522938441201887478E+00; 
    x[ 38] =  -0.57776938897061258000325165713764E+00; 
    x[ 39] =  -0.55748515286193223292186190687872E+00; 
    x[ 40] =  -0.53686246972339756745816636353452E+00;
    x[ 41] =  -0.51591385950424935727727729906662E+00; 
    x[ 42] =  -0.49465204002278211739494017368636E+00;
    x[ 43] =  -0.47308991924540524164509989939699E+00;
    x[ 44] =  -0.45124058745026622733189858020729E+00;
    x[ 45] =  -0.42911730928019337626254405355418E+00;
    x[ 46] =  -0.40673351568978256340867288124339E+00;
    x[ 47] =  -0.38410279579151693577907781452239E+00;
    x[ 48] =  -0.36123888860586970607092484346723E+00;
    x[ 49] =  -0.33815567472039850137600027657095E+00;
    x[ 50] =  -0.31486716786289498148601475374890E+00; 
    x[ 51] =  -0.29138750639370562079451875284568E+00; 
    x[ 52] =  -0.26773094472238862088834352027938E+00;
    x[ 53] =  -0.24391184465391785797071324453138E+00;
    x[ 54] =  -0.21994466666968754245452337866940E+00;
    x[ 55] =  -0.19584396114861085150428162519610E+00;
    x[ 56] =  -0.17162435953364216500834492248954E+00; 
    x[ 57] =  -0.14730056544908566938932929319807E+00;
    x[ 58] =  -0.12288734577408297172603365288567E+00;
    x[ 59] =  -0.98399521677698970751091751509101E-01;
    x[ 60] =  -0.73851959621048545273440409360569E-01;
    x[ 61] =  -0.49259562331926630315379321821927E-01;
    x[ 62] =  -0.24637259757420944614897071846088E-01;
    x[ 63] =   0.00000000000000000000000000000000E+00;
    x[ 64] =   0.24637259757420944614897071846088E-01;
    x[ 65] =   0.49259562331926630315379321821927E-01;
    x[ 66] =   0.73851959621048545273440409360569E-01;
    x[ 67] =   0.98399521677698970751091751509101E-01;
    x[ 68] =   0.12288734577408297172603365288567E+00;
    x[ 69] =   0.14730056544908566938932929319807E+00;
    x[ 70] =   0.17162435953364216500834492248954E+00;
    x[ 71] =   0.19584396114861085150428162519610E+00;
    x[ 72] =   0.21994466666968754245452337866940E+00;    
    x[ 73] =   0.24391184465391785797071324453138E+00;   
    x[ 74] =   0.26773094472238862088834352027938E+00;   
    x[ 75] =   0.29138750639370562079451875284568E+00;   
    x[ 76] =   0.31486716786289498148601475374890E+00;    
    x[ 77] =   0.33815567472039850137600027657095E+00;   
    x[ 78] =   0.36123888860586970607092484346723E+00;    
    x[ 79] =   0.38410279579151693577907781452239E+00;    
    x[ 80] =   0.40673351568978256340867288124339E+00;  
    x[ 81] =   0.42911730928019337626254405355418E+00;    
    x[ 82] =   0.45124058745026622733189858020729E+00;   
    x[ 83] =   0.47308991924540524164509989939699E+00;   
    x[ 84] =   0.49465204002278211739494017368636E+00; 
    x[ 85] =   0.51591385950424935727727729906662E+00; 
    x[ 86] =   0.53686246972339756745816636353452E+00; 
    x[ 87] =   0.55748515286193223292186190687872E+00;   
    x[ 88] =   0.57776938897061258000325165713764E+00;  
    x[ 89] =   0.59770286357006522938441201887478E+00;  
    x[ 90] =   0.61727347512685828385763916340822E+00;  
    x[ 91] =   0.63646934240029724134760815684175E+00;    
    x[ 92] =   0.65527881165548263027676505156852E+00;  
    x[ 93] =   0.67369046373825048534668253831602E+00;   
    x[ 94] =   0.69169312100770067015644143286666E+00;   
    x[ 95] =   0.70927585412210456099944463906757E+00;   
    x[ 96] =   0.72642798867407268553569290153270E+00;   
    x[ 97] =   0.74313911167095451292056688997595E+00;    
    x[ 98] =   0.75939907785653667155666366659810E+00;   
    x[ 99] =   0.77519801587020238244496276354566E+00;    
    x[100] =   0.79052633423981379994544995252740E+00;   
    x[101] =   0.80537472720468021466656079404644E+00;   
    x[102] =   0.81973418036507867415511910167470E+00;  
    x[103] =   0.83359597615489951437955716480123E+00;   
    x[104] =   0.84695169913409759845333931085437E+00;   
    x[105] =   0.85979324109774080981203134414483E+00; 
    x[106] =   0.87211280599856071141963753428864E+00;  
    x[107] =   0.88390291468002656994525794802849E+00;   
    x[108] =   0.89515640941708370896904382642451E+00;    
    x[109] =   0.90586645826182138280246131760282E+00;   
    x[110] =   0.91602655919146580931308861741716E+00;  
    x[111] =   0.92563054405623384912746466814259E+00; 
    x[112] =   0.93467258232473796857363487794906E+00; 
    x[113] =   0.94314718462481482734544963026201E+00;  
    x[114] =   0.95104920607788031054790764659636E+00; 
    x[115] =   0.95837384942523877114910286998060E+00; 
    x[116] =   0.96511666794529212109082507703391E+00;
    x[117] =   0.97127356816152919228894689830512E+00; 
    x[118] =   0.97684081234307032681744391886221E+00; 
    x[119] =   0.98181502080381411003346312451200E+00;  
    x[120] =   0.98619317401693166671043833175407E+00;
    x[121] =   0.98997261459148415760778669967548E+00;
    x[122] =   0.99315104925451714736113079489080E+00; 
    x[123] =   0.99572655135202722663543337085008E+00; 
    x[124] =   0.99769756618980462107441703193392E+00; 
    x[125] =   0.99906293435531189513828920479421E+00;
    x[126] =   0.99982213041530614629963254927125E+00; 

    w[  0] =   0.45645726109586654495731936146574E-03;
    w[  1] =   0.10622766869538486959954760554099E-02;
    w[  2] =   0.16683488125171936761028811985672E-02;
    w[  3] =   0.22734860707492547802810838362671E-02;
    w[  4] =   0.28772587656289004082883197417581E-02;
    w[  5] =   0.34792893810051465908910894094105E-02;
    w[  6] =   0.40792095178254605327114733456293E-02;
    w[  7] =   0.46766539777779034772638165662478E-02;
    w[  8] =   0.52712596565634400891303815906251E-02;
    w[  9] =   0.58626653903523901033648343751367E-02;
    w[ 10] =   0.64505120486899171845442463868748E-02;
    w[ 11] =   0.70344427036681608755685893032552E-02;
    w[ 12] =   0.76141028256526859356393930849227E-02;
    w[ 13] =   0.81891404887415730817235884718726E-02;
    w[ 14] =   0.87592065795403145773316804234385E-02;
    w[ 15] =   0.93239550065309714787536985834029E-02;
    w[ 16] =   0.98830429087554914716648010899606E-02;
    w[ 17] =   0.10436130863141005225673171997668E-01;
    w[ 18] =   0.10982883090068975788799657376065E-01;
    w[ 19] =   0.11522967656921087154811609734510E-01;
    w[ 20] =   0.12056056679400848183529562144697E-01;
    w[ 21] =   0.12581826520465013101514365424172E-01;
    w[ 22] =   0.13099957986718627426172681912499E-01;
    w[ 23] =   0.13610136522139249906034237533759E-01;
    w[ 24] =   0.14112052399003395774044161633613E-01;
    w[ 25] =   0.14605400905893418351737288078952E-01;
    w[ 26] =   0.15089882532666922992635733981431E-01;
    w[ 27] =   0.15565203152273955098532590262975E-01;
    w[ 28] =   0.16031074199309941802254151842763E-01;
    w[ 29] =   0.16487212845194879399346060358146E-01;
    w[ 30] =   0.16933342169871654545878815295200E-01;
    w[ 31] =   0.17369191329918731922164721250350E-01;
    w[ 32] =   0.17794495722974774231027912900351E-01;
    w[ 33] =   0.18208997148375106468721469154479E-01;
    w[ 34] =   0.18612443963902310429440419898958E-01;
    w[ 35] =   0.19004591238555646611148901044533E-01;
    w[ 36] =   0.19385200901246454628112623489471E-01;
    w[ 37] =   0.19754041885329183081815217323169E-01;
    w[ 38] =   0.20110890268880247225644623956287E-01;
    w[ 39] =   0.20455529410639508279497065713301E-01;
    w[ 40] =   0.20787750081531811812652137291250E-01;
    w[ 41] =   0.21107350591688713643523847921658E-01;
    w[ 42] =   0.21414136912893259295449693233545E-01;
    w[ 43] =   0.21707922796373466052301324695331E-01;
    w[ 44] =   0.21988529885872983756478409758807E-01;
    w[ 45] =   0.22255787825930280235631416460158E-01;
    w[ 46] =   0.22509534365300608085694429903050E-01;
    w[ 47] =   0.22749615455457959852242553240982E-01;
    w[ 48] =   0.22975885344117206754377437838947E-01;
    w[ 49] =   0.23188206663719640249922582981729E-01;
    w[ 50] =   0.23386450514828194170722043496950E-01;
    w[ 51] =   0.23570496544381716050033676844306E-01;
    w[ 52] =   0.23740233018760777777714726703424E-01;
    w[ 53] =   0.23895556891620665983864481754172E-01;
    w[ 54] =   0.24036373866450369675132086026456E-01;
    w[ 55] =   0.24162598453819584716522917710986E-01;
    w[ 56] =   0.24274154023278979833195063936748E-01;
    w[ 57] =   0.24370972849882214952813561907241E-01;
    w[ 58] =   0.24452996155301467956140198471529E-01;
    w[ 59] =   0.24520174143511508275183033290175E-01;
    w[ 60] =   0.24572466031020653286354137335186E-01;
    w[ 61] =   0.24609840071630254092545634003360E-01;
    w[ 62] =   0.24632273575707679066033370218017E-01;
    w[ 63] =   0.24639752923961094419579417477503E-01;
    w[ 64] =   0.24632273575707679066033370218017E-01;
    w[ 65] =   0.24609840071630254092545634003360E-01;
    w[ 66] =   0.24572466031020653286354137335186E-01;
    w[ 67] =   0.24520174143511508275183033290175E-01;
    w[ 68] =   0.24452996155301467956140198471529E-01;
    w[ 69] =   0.24370972849882214952813561907241E-01;
    w[ 70] =   0.24274154023278979833195063936748E-01;
    w[ 71] =   0.24162598453819584716522917710986E-01;
    w[ 72] =   0.24036373866450369675132086026456E-01;
    w[ 73] =   0.23895556891620665983864481754172E-01;
    w[ 74] =   0.23740233018760777777714726703424E-01;
    w[ 75] =   0.23570496544381716050033676844306E-01;
    w[ 76] =   0.23386450514828194170722043496950E-01;
    w[ 77] =   0.23188206663719640249922582981729E-01;
    w[ 78] =   0.22975885344117206754377437838947E-01;
    w[ 79] =   0.22749615455457959852242553240982E-01;
    w[ 80] =   0.22509534365300608085694429903050E-01;
    w[ 81] =   0.22255787825930280235631416460158E-01;
    w[ 82] =   0.21988529885872983756478409758807E-01;
    w[ 83] =   0.21707922796373466052301324695331E-01;
    w[ 84] =   0.21414136912893259295449693233545E-01;
    w[ 85] =   0.21107350591688713643523847921658E-01;
    w[ 86] =   0.20787750081531811812652137291250E-01;
    w[ 87] =   0.20455529410639508279497065713301E-01;
    w[ 88] =   0.20110890268880247225644623956287E-01;
    w[ 89] =   0.19754041885329183081815217323169E-01;
    w[ 90] =   0.19385200901246454628112623489471E-01;
    w[ 91] =   0.19004591238555646611148901044533E-01;
    w[ 92] =   0.18612443963902310429440419898958E-01;
    w[ 93] =   0.18208997148375106468721469154479E-01;
    w[ 94] =   0.17794495722974774231027912900351E-01;
    w[ 95] =   0.17369191329918731922164721250350E-01;
    w[ 96] =   0.16933342169871654545878815295200E-01;
    w[ 97] =   0.16487212845194879399346060358146E-01;
    w[ 98] =   0.16031074199309941802254151842763E-01;
    w[ 99] =   0.15565203152273955098532590262975E-01;
    w[100] =   0.15089882532666922992635733981431E-01;
    w[101] =   0.14605400905893418351737288078952E-01;
    w[102] =   0.14112052399003395774044161633613E-01;
    w[103] =   0.13610136522139249906034237533759E-01;
    w[104] =   0.13099957986718627426172681912499E-01;
    w[105] =   0.12581826520465013101514365424172E-01;
    w[106] =   0.12056056679400848183529562144697E-01;
    w[107] =   0.11522967656921087154811609734510E-01;
    w[108] =   0.10982883090068975788799657376065E-01;
    w[109] =   0.10436130863141005225673171997668E-01;
    w[110] =   0.98830429087554914716648010899606E-02;
    w[111] =   0.93239550065309714787536985834029E-02;
    w[112] =   0.87592065795403145773316804234385E-02;
    w[113] =   0.81891404887415730817235884718726E-02;
    w[114] =   0.76141028256526859356393930849227E-02;
    w[115] =   0.70344427036681608755685893032552E-02;
    w[116] =   0.64505120486899171845442463868748E-02;
    w[117] =   0.58626653903523901033648343751367E-02;
    w[118] =   0.52712596565634400891303815906251E-02;
    w[119] =   0.46766539777779034772638165662478E-02;
    w[120] =   0.40792095178254605327114733456293E-02;
    w[121] =   0.34792893810051465908910894094105E-02;
    w[122] =   0.28772587656289004082883197417581E-02;
    w[123] =   0.22734860707492547802810838362671E-02;
    w[124] =   0.16683488125171936761028811985672E-02;
    w[125] =   0.10622766869538486959954760554099E-02;
    w[126] =   0.45645726109586654495731936146574E-03;
  }
  else if ( n == 255 )
  {
    x[ 0] =      -0.9999557053175637;
    x[ 1] =      -0.9997666213120006;
    x[ 2] =        -0.99942647468017;
    x[ 3] =      -0.9989352412846546;
    x[ 4] =      -0.9982929861369679;
    x[ 5] =      -0.9974998041266158;
    x[ 6] =      -0.9965558144351986;
    x[ 7] =      -0.9954611594800263;
    x[ 8] =      -0.9942160046166302;
    x[ 9] =      -0.9928205380219891;
    x[10] =      -0.9912749706303856;
    x[11] =      -0.9895795360859201;
    x[12] =      -0.9877344906997324;
    x[13] =      -0.9857401134074193;
    x[14] =      -0.9835967057247763;
    x[15] =      -0.9813045917010171;
    x[16] =      -0.9788641178690681;
    x[17] =       -0.976275653192736;
    x[18] =      -0.9735395890106436;
    x[19] =      -0.9706563389768804;
    x[20] =      -0.9676263389983388;
    x[21] =      -0.9644500471687263;
    x[22] =      -0.9611279436992478;
    x[23] =       -0.957660530845962;
    x[24] =      -0.9540483328338163;
    x[25] =      -0.9502918957773683;
    x[26] =      -0.9463917875982043;
    x[27] =      -0.9423485979390644;
    x[28] =      -0.9381629380746873;
    x[29] =      -0.9338354408193861;
    x[30] =      -0.9293667604313699;
    x[31] =      -0.9247575725138244;
    x[32] =      -0.9200085739127664;
    x[33] =       -0.915120482611687;
    x[34] =      -0.9100940376230008;
    x[35] =       -0.904929998876315;
    x[36] =      -0.8996291471035368;
    x[37] =      -0.8941922837208367;
    x[38] =      -0.8886202307074841;
    x[39] =      -0.8829138304815741;
    x[40] =      -0.8770739457726654;
    x[41] =      -0.8711014594913465;
    x[42] =      -0.8649972745957512;
    x[43] =       -0.858762313955043;
    x[44] =      -0.8523975202098902;
    x[45] =      -0.8459038556299511;
    x[46] =       -0.839282301968391;
    x[47] =      -0.8325338603134556;
    x[48] =      -0.8256595509371186;
    x[49] =      -0.8186604131408319;
    x[50] =      -0.8115375050983958;
    x[51] =      -0.8042919036959787;
    x[52] =      -0.7969247043693057;
    x[53] =      -0.7894370209380444;
    x[54] =      -0.7818299854374094;
    x[55] =      -0.7741047479470157;
    x[56] =      -0.7662624764170006;
    x[57] =      -0.7583043564914468;
    x[58] =      -0.7502315913291283;
    x[59] =      -0.7420454014216102;
    x[60] =      -0.7337470244087263;
    x[61] =      -0.7253377148914649;
    x[62] =      -0.7168187442422908;
    x[63] =      -0.7081914004129306;
    x[64] =      -0.6994569877396524;
    x[65] =      -0.6906168267460676;
    x[66] =      -0.6816722539434864;
    x[67] =      -0.6726246216288551;
    x[68] =       -0.663475297680307;
    x[69] =      -0.6542256653503588;
    x[70] =      -0.6448771230567811;
    x[71] =      -0.6354310841711771;
    x[72] =      -0.6258889768052999;
    x[73] =      -0.6162522435951415;
    x[74] =      -0.6065223414828266;
    x[75] =      -0.5967007414963417;
    x[76] =      -0.5867889285271373;
    x[77] =      -0.5767884011056313;
    x[78] =      -0.5667006711746527;
    x[79] =      -0.5565272638608558;
    x[80] =      -0.5462697172441424;
    x[81] =      -0.5359295821251249;
    x[82] =      -0.5255084217906666;
    x[83] =      -0.5150078117775342;
    x[84] =      -0.5044293396341982;
    x[85] =       -0.493774604680817;
    x[86] =       -0.483045217767442;
    x[87] =      -0.4722428010304787;
    x[88] =      -0.4613689876474424;
    x[89] =      -0.4504254215900437;
    x[90] =      -0.4394137573756426;
    x[91] =      -0.4283356598171081;
    x[92] =      -0.4171928037711214;
    x[93] =      -0.4059868738849605;
    x[94] =      -0.3947195643418044;
    x[95] =      -0.3833925786045958;
    x[96] =      -0.3720076291585012;
    x[97] =      -0.3605664372520062;
    x[98] =      -0.3490707326366864;
    x[99] =      -0.3375222533056927;
    x[100] =      -0.3259227452309905;
    x[101] =      -0.3142739620993925;
    x[102] =      -0.3025776650474256;
    x[103] =      -0.2908356223950708;
    x[104] =      -0.2790496093784178;
    x[105] =      -0.2672214078812731;
    x[106] =      -0.2553528061657641;
    x[107] =       -0.243445598601978;
    x[108] =      -0.2315015853966777;
    x[109] =      -0.2195225723211354;
    x[110] =      -0.2075103704381242;
    x[111] =      -0.1954667958281108;
    x[112] =      -0.1833936693146885;
    x[113] =      -0.1712928161892939;
    x[114] =      -0.1591660659352477;
    x[115] =       -0.147015251951162;
    x[116] =      -0.1348422112737553;
    x[117] =      -0.1226487843001178;
    x[118] =      -0.1104368145094688;
    x[119] =     -0.09820814818444755;
    x[120] =     -0.08596463413198061;
    x[121] =     -0.07370812340376778;
    x[122] =     -0.06144046901642827;
    x[123] =     -0.04916352567134998;
    x[124] =     -0.03687914947428402;
    x[125] =     -0.02458919765472701;
    x[126] =     -0.01229552828513332;
    x[127] =                        0;
    x[128] =      0.01229552828513332;
    x[129] =      0.02458919765472701;
    x[130] =      0.03687914947428402;
    x[131] =      0.04916352567134998;
    x[132] =      0.06144046901642827;
    x[133] =      0.07370812340376778;
    x[134] =      0.08596463413198061;
    x[135] =      0.09820814818444755;
    x[136] =       0.1104368145094688;
    x[137] =       0.1226487843001178;
    x[138] =       0.1348422112737553;
    x[139] =        0.147015251951162;
    x[140] =       0.1591660659352477;
    x[141] =       0.1712928161892939;
    x[142] =       0.1833936693146885;
    x[143] =       0.1954667958281108;
    x[144] =       0.2075103704381242;
    x[145] =       0.2195225723211354;
    x[146] =       0.2315015853966777;
    x[147] =        0.243445598601978;
    x[148] =       0.2553528061657641;
    x[149] =       0.2672214078812731;
    x[150] =       0.2790496093784178;
    x[151] =       0.2908356223950708;
    x[152] =       0.3025776650474256;
    x[153] =       0.3142739620993925;
    x[154] =       0.3259227452309905;
    x[155] =       0.3375222533056927;
    x[156] =       0.3490707326366864;
    x[157] =       0.3605664372520062;
    x[158] =       0.3720076291585012;
    x[159] =       0.3833925786045958;
    x[160] =       0.3947195643418044;
    x[161] =       0.4059868738849605;
    x[162] =       0.4171928037711214;
    x[163] =       0.4283356598171081;
    x[164] =       0.4394137573756426;
    x[165] =       0.4504254215900437;
    x[166] =       0.4613689876474424;
    x[167] =       0.4722428010304787;
    x[168] =        0.483045217767442;
    x[169] =        0.493774604680817;
    x[170] =       0.5044293396341982;
    x[171] =       0.5150078117775342;
    x[172] =       0.5255084217906666;
    x[173] =       0.5359295821251249;
    x[174] =       0.5462697172441424;
    x[175] =       0.5565272638608558;
    x[176] =       0.5667006711746527;
    x[177] =       0.5767884011056313;
    x[178] =       0.5867889285271373;
    x[179] =       0.5967007414963417;
    x[180] =       0.6065223414828266;
    x[181] =       0.6162522435951415;
    x[182] =       0.6258889768052999;
    x[183] =       0.6354310841711771;
    x[184] =       0.6448771230567811;
    x[185] =       0.6542256653503588;
    x[186] =        0.663475297680307;
    x[187] =       0.6726246216288551;
    x[188] =       0.6816722539434864;
    x[189] =       0.6906168267460676;
    x[190] =       0.6994569877396524;
    x[191] =       0.7081914004129306;
    x[192] =       0.7168187442422908;
    x[193] =       0.7253377148914649;
    x[194] =       0.7337470244087263;
    x[195] =       0.7420454014216102;
    x[196] =       0.7502315913291283;
    x[197] =       0.7583043564914468;
    x[198] =       0.7662624764170006;
    x[199] =       0.7741047479470157;
    x[200] =       0.7818299854374094;
    x[201] =       0.7894370209380444;
    x[202] =       0.7969247043693057;
    x[203] =       0.8042919036959787;
    x[204] =       0.8115375050983958;
    x[205] =       0.8186604131408319;
    x[206] =       0.8256595509371186;
    x[207] =       0.8325338603134556;
    x[208] =        0.839282301968391;
    x[209] =       0.8459038556299511;
    x[210] =       0.8523975202098902;
    x[211] =        0.858762313955043;
    x[212] =       0.8649972745957512;
    x[213] =       0.8711014594913465;
    x[214] =       0.8770739457726654;
    x[215] =       0.8829138304815741;
    x[216] =       0.8886202307074841;
    x[217] =       0.8941922837208367;
    x[218] =       0.8996291471035368;
    x[219] =        0.904929998876315;
    x[220] =       0.9100940376230008;
    x[221] =        0.915120482611687;
    x[222] =       0.9200085739127664;
    x[223] =       0.9247575725138244;
    x[224] =       0.9293667604313699;
    x[225] =       0.9338354408193861;
    x[226] =       0.9381629380746873;
    x[227] =       0.9423485979390644;
    x[228] =       0.9463917875982043;
    x[229] =       0.9502918957773683;
    x[230] =       0.9540483328338163;
    x[231] =        0.957660530845962;
    x[232] =       0.9611279436992478;
    x[233] =       0.9644500471687263;
    x[234] =       0.9676263389983388;
    x[235] =       0.9706563389768804;
    x[236] =       0.9735395890106436;
    x[237] =        0.976275653192736;
    x[238] =       0.9788641178690681;
    x[239] =       0.9813045917010171;
    x[240] =       0.9835967057247763;
    x[241] =       0.9857401134074193;
    x[242] =       0.9877344906997324;
    x[243] =       0.9895795360859201;
    x[244] =       0.9912749706303856;
    x[245] =       0.9928205380219891;
    x[246] =       0.9942160046166302;
    x[247] =       0.9954611594800263;
    x[248] =       0.9965558144351986;
    x[249] =       0.9974998041266158;
    x[250] =       0.9982929861369679;
    x[251] =       0.9989352412846546;
    x[252] =         0.99942647468017;
    x[253] =       0.9997666213120006;
    x[254] =       0.9999557053175637;

    w[ 0] =    0.0001136736199914808;
    w[ 1] =    0.0002645938711908564;
    w[ 2] =    0.0004156976252681932;
    w[ 3] =    0.0005667579456482639;
    w[ 4] =    0.0007177364780061286;
    w[ 5] =    0.0008686076661194581;
    w[ 6] =     0.001019347976427318;
    w[ 7] =       0.0011699343729388;
    w[ 8] =     0.001320343990022177;
    w[ 9] =     0.001470554042778403;
    w[10] =     0.001620541799041545;
    w[11] =     0.001770284570660304;
    w[12] =     0.001919759711713187;
    w[13] =     0.002068944619501569;
    w[14] =     0.002217816736754017;
    w[15] =     0.002366353554396287;
    w[16] =      0.00251453261459971;
    w[17] =     0.002662331513971696;
    w[18] =      0.00280972790682046;
    w[19] =     0.002956699508457498;
    w[20] =     0.003103224098519095;
    w[21] =     0.003249279524294296;
    w[22] =     0.003394843704053401;
    w[23] =     0.003539894630372244;
    w[24] =     0.003684410373449933;
    w[25] =     0.003828369084417135;
    w[26] =     0.003971748998634907;
    w[27] =     0.004114528438981242;
    w[28] =     0.004256685819126112;
    w[29] =     0.004398199646792759;
    w[30] =      0.00453904852700618;
    w[31] =     0.004679211165326077;
    w[32] =     0.004818666371065699;
    w[33] =      0.00495739306049505;
    w[34] =     0.005095370260027839;
    w[35] =     0.005232577109391968;
    w[36] =     0.005368992864783177;
    w[37] =     0.005504596902000804;
    w[38] =     0.005639368719565862;
    w[39] =     0.005773287941820301;
    w[40] =     0.005906334322007422;
    w[41] =     0.006038487745332765;
    w[42] =     0.006169728232005295;
    w[43] =     0.006300035940257733;
    w[44] =     0.006429391169346602;
    w[45] =     0.006557774362530328;
    w[46] =     0.006685166110026254;
    w[47] =     0.006811547151944815;
    w[48] =     0.006936898381201466;
    w[49] =     0.007061200846405536;
    w[50] =     0.007184435754724984;
    w[51] =     0.007306584474728122;
    w[52] =     0.007427628539199977;
    w[53] =     0.007547549647934514;
    w[54] =     0.007666329670501377;
    w[55] =     0.007783950648986801;
    w[56] =     0.007900394800708624;
    w[57] =     0.008015644520904983;
    w[58] =     0.008129682385395602;
    w[59] =     0.008242491153216323;
    w[60] =     0.008354053769225508;
    w[61] =     0.008464353366682819;
    w[62] =     0.008573373269798925;
    w[63] =     0.008681096996256795;
    w[64] =     0.008787508259703609;
    w[65] =     0.008892590972213036;
    w[66] =     0.008996329246717397;
    w[67] =     0.009098707399409718;
    w[68] =     0.009199709952114802;
    w[69] =     0.009299321634629343;
    w[70] =     0.009397527387030594;
    w[71] =     0.009494312361953241;
    w[72] =     0.009589661926834022;
    w[73] =     0.009683561666124043;
    w[74] =     0.009775997383468165;
    w[75] =     0.009866955103851452;
    w[76] =     0.009956421075711706;
    w[77] =      0.01004438177301882;
    w[78] =      0.01013082389731963;
    w[79] =      0.01021573437974821;
    w[80] =       0.0102991003830022;
    w[81] =      0.01038090930328312;
    w[82] =      0.01046114877220228;
    w[83] =      0.01053980665865038;
    w[84] =      0.01061687107063194;
    w[85] =      0.01069233035706287;
    w[86] =      0.01076617310953212;
    w[87] =      0.01083838816402652;
    w[88] =      0.01090896460261843;
    w[89] =      0.01097789175511656;
    w[90] =      0.01104515920067912;
    w[91] =      0.01111075676938929;
    w[92] =      0.01117467454379268;
    w[93] =      0.01123690286039691;
    w[94] =      0.01129743231113249;
    w[95] =      0.01135625374477508;
    w[96] =      0.01141335826832922;
    w[97] =      0.01146873724837283;
    w[98] =      0.01152238231236217;
    w[99] =      0.01157428534989815;
    w[100] =      0.01162443851395193;
    w[101] =      0.01167283422205182;
    w[102] =      0.01171946515742932;
    w[103] =      0.01176432427012535;
    w[104] =      0.01180740477805627;
    w[105] =      0.01184870016803913;
    w[106] =      0.01188820419677619;
    w[107] =      0.01192591089179929;
    w[108] =      0.01196181455237226;
    w[109] =      0.01199590975035326;
    w[110] =      0.01202819133101508;
    w[111] =      0.01205865441382472;
    w[112] =      0.01208729439318107;
    w[113] =      0.01211410693911137;
    w[114] =      0.01213908799792579;
    w[115] =      0.01216223379283022;
    w[116] =      0.01218354082449738;
    w[117] =      0.01220300587159574;
    w[118] =      0.01222062599127671;
    w[119] =      0.01223639851961942;
    w[120] =      0.01225032107203351;
    w[121] =      0.01226239154361966;
    w[122] =      0.01227260810948789;
    w[123] =      0.01228096922503318;
    w[124] =      0.01228747362616942;
    w[125] =      0.01229212032952021;
    w[126] =      0.01229490863256759;
    w[127] =      0.01229583811375833;
    w[128] =      0.01229490863256759;
    w[129] =      0.01229212032952021;
    w[130] =      0.01228747362616942;
    w[131] =      0.01228096922503318;
    w[132] =      0.01227260810948789;
    w[133] =      0.01226239154361966;
    w[134] =      0.01225032107203351;
    w[135] =      0.01223639851961942;
    w[136] =      0.01222062599127671;
    w[137] =      0.01220300587159574;
    w[138] =      0.01218354082449738;
    w[139] =      0.01216223379283022;
    w[140] =      0.01213908799792579;
    w[141] =      0.01211410693911137;
    w[142] =      0.01208729439318107;
    w[143] =      0.01205865441382472;
    w[144] =      0.01202819133101508;
    w[145] =      0.01199590975035326;
    w[146] =      0.01196181455237226;
    w[147] =      0.01192591089179929;
    w[148] =      0.01188820419677619;
    w[149] =      0.01184870016803913;
    w[150] =      0.01180740477805627;
    w[151] =      0.01176432427012535;
    w[152] =      0.01171946515742932;
    w[153] =      0.01167283422205182;
    w[154] =      0.01162443851395193;
    w[155] =      0.01157428534989815;
    w[156] =      0.01152238231236217;
    w[157] =      0.01146873724837283;
    w[158] =      0.01141335826832922;
    w[159] =      0.01135625374477508;
    w[160] =      0.01129743231113249;
    w[161] =      0.01123690286039691;
    w[162] =      0.01117467454379268;
    w[163] =      0.01111075676938929;
    w[164] =      0.01104515920067912;
    w[165] =      0.01097789175511656;
    w[166] =      0.01090896460261843;
    w[167] =      0.01083838816402652;
    w[168] =      0.01076617310953212;
    w[169] =      0.01069233035706287;
    w[170] =      0.01061687107063194;
    w[171] =      0.01053980665865038;
    w[172] =      0.01046114877220228;
    w[173] =      0.01038090930328312;
    w[174] =       0.0102991003830022;
    w[175] =      0.01021573437974821;
    w[176] =      0.01013082389731963;
    w[177] =      0.01004438177301882;
    w[178] =     0.009956421075711706;
    w[179] =     0.009866955103851452;
    w[180] =     0.009775997383468165;
    w[181] =     0.009683561666124043;
    w[182] =     0.009589661926834022;
    w[183] =     0.009494312361953241;
    w[184] =     0.009397527387030594;
    w[185] =     0.009299321634629343;
    w[186] =     0.009199709952114802;
    w[187] =     0.009098707399409718;
    w[188] =     0.008996329246717397;
    w[189] =     0.008892590972213036;
    w[190] =     0.008787508259703609;
    w[191] =     0.008681096996256795;
    w[192] =     0.008573373269798925;
    w[193] =     0.008464353366682819;
    w[194] =     0.008354053769225508;
    w[195] =     0.008242491153216323;
    w[196] =     0.008129682385395602;
    w[197] =     0.008015644520904983;
    w[198] =     0.007900394800708624;
    w[199] =     0.007783950648986801;
    w[200] =     0.007666329670501377;
    w[201] =     0.007547549647934514;
    w[202] =     0.007427628539199977;
    w[203] =     0.007306584474728122;
    w[204] =     0.007184435754724984;
    w[205] =     0.007061200846405536;
    w[206] =     0.006936898381201466;
    w[207] =     0.006811547151944815;
    w[208] =     0.006685166110026254;
    w[209] =     0.006557774362530328;
    w[210] =     0.006429391169346602;
    w[211] =     0.006300035940257733;
    w[212] =     0.006169728232005295;
    w[213] =     0.006038487745332765;
    w[214] =     0.005906334322007422;
    w[215] =     0.005773287941820301;
    w[216] =     0.005639368719565862;
    w[217] =     0.005504596902000804;
    w[218] =     0.005368992864783177;
    w[219] =     0.005232577109391968;
    w[220] =     0.005095370260027839;
    w[221] =      0.00495739306049505;
    w[222] =     0.004818666371065699;
    w[223] =     0.004679211165326077;
    w[224] =      0.00453904852700618;
    w[225] =     0.004398199646792759;
    w[226] =     0.004256685819126112;
    w[227] =     0.004114528438981242;
    w[228] =     0.003971748998634907;
    w[229] =     0.003828369084417135;
    w[230] =     0.003684410373449933;
    w[231] =     0.003539894630372244;
    w[232] =     0.003394843704053401;
    w[233] =     0.003249279524294296;
    w[234] =     0.003103224098519095;
    w[235] =     0.002956699508457498;
    w[236] =      0.00280972790682046;
    w[237] =     0.002662331513971696;
    w[238] =      0.00251453261459971;
    w[239] =     0.002366353554396287;
    w[240] =     0.002217816736754017;
    w[241] =     0.002068944619501569;
    w[242] =     0.001919759711713187;
    w[243] =     0.001770284570660304;
    w[244] =     0.001620541799041545;
    w[245] =     0.001470554042778403;
    w[246] =     0.001320343990022177;
    w[247] =       0.0011699343729388;
    w[248] =     0.001019347976427318;
    w[249] =    0.0008686076661194581;
    w[250] =    0.0007177364780061286;
    w[251] =    0.0005667579456482639;
    w[252] =    0.0004156976252681932;
    w[253] =    0.0002645938711908564;
    w[254] =    0.0001136736199914808;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1 through 33, 63, 64, 65, 127 or 255.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_x1 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_X1 sets a Gauss-Legendre rule for ( 1 + X ) * F(X) on [-1,1].
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x-1] = 1 + x.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) ( 1 + X ) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be between 1 and 9.
//
//    Output, double XTAB[ORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[ORDER], the weights of the rule.
//
{
  if ( order == 1 )
  {
    xtab[1-1] =  0.333333333333333333333333333333E+00;

    weight[1-1] = 2.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = -0.289897948556635619639456814941E+00;
    xtab[2-1] =  0.689897948556635619639456814941E+00;

    weight[1-1] =  0.727834473024091322422523991699E+00;
    weight[2-1] =  1.27216552697590867757747600830E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = -0.575318923521694112050483779752E+00;
    xtab[2-1] =  0.181066271118530578270147495862E+00;
    xtab[3-1] =  0.822824080974592105208907712461E+00;

    weight[1-1] =  0.279307919605816490135525088716E+00;
    weight[2-1] =  0.916964425438344986775682378225E+00;
    weight[3-1] =  0.803727654955838523088792533058E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = -0.720480271312438895695825837750E+00;
    xtab[2-1] = -0.167180864737833640113395337326E+00;
    xtab[3-1] =  0.446313972723752344639908004629E+00;
    xtab[4-1] =  0.885791607770964635613757614892E+00;

    weight[1-1] =  0.124723883800032328695500588386E+00;
    weight[2-1] =  0.519390190432929763305824811559E+00;
    weight[3-1] =  0.813858272041085443165617903743E+00;
    weight[4-1] =  0.542027653725952464833056696312E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = -0.802929828402347147753002204224E+00;
    xtab[2-1] = -0.390928546707272189029229647442E+00;
    xtab[3-1] =  0.124050379505227711989974959990E+00;
    xtab[4-1] =  0.603973164252783654928415726409E+00;
    xtab[5-1] =  0.920380285897062515318386619813E+00;

    weight[1-1] =  0.0629916580867691047411692662740E+00;
    weight[2-1] =  0.295635480290466681402532877367E+00;
    weight[3-1] =  0.585547948338679234792151477424E+00;
    weight[4-1] =  0.668698552377478261966702492391E+00;
    weight[5-1] =  0.387126360906606717097443886545E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = -0.853891342639482229703747931639E+00;
    xtab[2-1] = -0.538467724060109001833766720231E+00;
    xtab[3-1] = -0.117343037543100264162786683611E+00;
    xtab[4-1] =  0.326030619437691401805894055838E+00;
    xtab[5-1] =  0.703842800663031416300046295008E+00;
    xtab[6-1] =  0.941367145680430216055899446174E+00;

    weight[1-1] =  0.0349532072544381270240692132496E+00;
    weight[2-1] =  0.175820662202035902032706497222E+00;
    weight[3-1] =  0.394644603562621056482338042193E+00;
    weight[4-1] =  0.563170215152795712476307356284E+00;
    weight[5-1] =  0.542169988926074467362761586552E+00;
    weight[6-1] =  0.289241322902034734621817304499E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = -0.887474878926155707068695617935E+00;
    xtab[2-1] = -0.639518616526215270024840114382E+00;
    xtab[3-1] = -0.294750565773660725252184459658E+00;
    xtab[4-1] =  0.0943072526611107660028971153047E+00;
    xtab[5-1] =  0.468420354430821063046421216613E+00;
    xtab[6-1] =  0.770641893678191536180719525865E+00;
    xtab[7-1] =  0.955041227122575003782349000858E+00;

    weight[1-1] =  0.0208574488112296163587654972151E+00;
    weight[2-1] =  0.109633426887493901777324193433E+00;
    weight[3-1] =  0.265538785861965879934591955055E+00;
    weight[4-1] =  0.428500262783494679963649011999E+00;
    weight[5-1] =  0.509563589198353307674937943100E+00;
    weight[6-1] =  0.442037032763498409684482945478E+00;
    weight[7-1] =  0.223869453693964204606248453720E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = -0.910732089420060298533757956283E+00;
    xtab[2-1] = -0.711267485915708857029562959544E+00;
    xtab[3-1] = -0.426350485711138962102627520502E+00;
    xtab[4-1] = -0.0903733696068532980645444599064E+00;
    xtab[5-1] =  0.256135670833455395138292079035E+00;
    xtab[6-1] =  0.571383041208738483284917464837E+00;
    xtab[7-1] =  0.817352784200412087992517083851E+00;
    xtab[8-1] =  0.964440169705273096373589797925E+00;

    weight[1-1] =  0.0131807657689951954189692640444E+00;
    weight[2-1] =  0.0713716106239448335742111888042E+00;
    weight[3-1] =  0.181757278018795592332221684383E+00;
    weight[4-1] =  0.316798397969276640481632757440E+00;
    weight[5-1] =  0.424189437743720042818124385645E+00;
    weight[6-1] =  0.450023197883549464687088394417E+00;
    weight[7-1] =  0.364476094545494505382889847132E+00;
    weight[8-1] =  0.178203217446223725304862478136E+00;
  }
  else if ( order == 9 )
  {
    xtab[1-1] = -0.927484374233581078117671398464E+00;
    xtab[2-1] = -0.763842042420002599615429776011E+00;
    xtab[3-1] = -0.525646030370079229365386614293E+00;
    xtab[4-1] = -0.236234469390588049278459503207E+00;
    xtab[5-1] =  0.0760591978379781302337137826389E+00;
    xtab[6-1] =  0.380664840144724365880759065541E+00;
    xtab[7-1] =  0.647766687674009436273648507855E+00;
    xtab[8-1] =  0.851225220581607910728163628088E+00;
    xtab[9-1] =  0.971175180702246902734346518378E+00;

    weight[1-1] =  0.00872338834309252349019620448007E+00;
    weight[2-1] =  0.0482400171391415162069086091476E+00;
    weight[3-1] =  0.127219285964216005046760427743E+00;
    weight[4-1] =  0.233604781180660442262926091607E+00;
    weight[5-1] =  0.337433287379681397577000079834E+00;
    weight[6-1] =  0.401235236773473158616600898930E+00;
    weight[7-1] =  0.394134968689382820640692081477E+00;
    weight[8-1] =  0.304297020437232650320317215016E+00;
    weight[9-1] =  0.145112014093119485838598391765E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_X1 - Fatal error!\n";
    cerr << "  Illegal input value of ORDER = " << order << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

void legendre_set_x2 ( int order, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_X2 sets Gauss-Legendre rules for ( 1 + X )^2*F(X) on [-1,1].
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x-1] = ( 1 + x )^2.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) ( 1 + X )^2 * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be between 1 and 9.
//
//    Output, double XTAB[ORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[ORDER], the weights of the rule.
//
{
  if ( order == 1 )
  {
    xtab[1-1] =  0.5E+00;

    weight[1-1] =  2.66666666666666666666666666666E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = -0.0883036880224505775998524725910E+00;
    xtab[2-1] =  0.754970354689117244266519139258E+00;

    weight[1-1] =  0.806287056638603444666851075928E+00;
    weight[2-1] =  1.86037961002806322199981559074E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = -0.410004419776996766244796955168E+00;
    xtab[2-1] =  0.305992467923296230556472913192E+00;
    xtab[3-1] =  0.854011951853700535688324041976E+00;

    weight[1-1] =  0.239605624068645584091811926047E+00;
    weight[2-1] =  1.16997015407892817602809616291E+00;
    weight[3-1] =  1.25709088851909290654675857771E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = -0.591702835793545726606755921586E+00;
    xtab[2-1] = -0.0340945902087350046811467387661E+00;
    xtab[3-1] =  0.522798524896275389882037174551E+00;
    xtab[4-1] =  0.902998901106005341405865485802E+00;

    weight[1-1] =  0.0828179259993445222751812523731E+00;
    weight[2-1] =  0.549071097383384602539010760334E+00;
    weight[3-1] =  1.14767031839371367238662411421E+00;
    weight[4-1] =  0.887107324890223869465850539752E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = -0.702108425894032836232448374820E+00;
    xtab[2-1] = -0.268666945261773544694327777841E+00;
    xtab[3-1] =  0.220227225868961343518209179230E+00;
    xtab[4-1] =  0.653039358456608553790815164028E+00;
    xtab[5-1] =  0.930842120163569816951085142737E+00;

    weight[1-1] =  0.0329106016247920636689299329544E+00;
    weight[2-1] =  0.256444805783695354037991444453E+00;
    weight[3-1] =  0.713601289772720001490035944563E+00;
    weight[4-1] =  1.00959169519929190423066348132E+00;
    weight[5-1] =  0.654118274286167343239045863379E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = -0.773611232355123732602532012021E+00;
    xtab[2-1] = -0.431362254623427837535325249187E+00;
    xtab[3-1] = -0.0180728263295041680220798103354E+00;
    xtab[4-1] =  0.395126163954217534500188844163E+00;
    xtab[5-1] =  0.736872116684029732026178298518E+00;
    xtab[6-1] =  0.948190889812665614490712786006E+00;

    weight[1-1] =  0.0146486064549543818622276447204E+00;
    weight[2-1] =  0.125762377479560410622810097040E+00;
    weight[3-1] =  0.410316569036929681761034600615E+00;
    weight[4-1] =  0.756617493988329628546336413760E+00;
    weight[5-1] =  0.859011997894245060846045458784E+00;
    weight[6-1] =  0.500309621812647503028212451747E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = -0.822366333126005527278634734418E+00;
    xtab[2-1] = -0.547034493182875002223997992852E+00;
    xtab[3-1] = -0.200043026557985860387937545780E+00;
    xtab[4-1] =  0.171995710805880507163425502299E+00;
    xtab[5-1] =  0.518891747903884926692601716998E+00;
    xtab[6-1] =  0.793821941703901970495546427988E+00;
    xtab[7-1] =  0.959734452453198985538996625765E+00;

    weight[1-1] =  0.00714150426951365443207221475404E+00;
    weight[2-1] =  0.0653034050584375560578544725498E+00;
    weight[3-1] =  0.235377690316228918725962815880E+00;
    weight[4-1] =  0.505171029671130381676271523850E+00;
    weight[5-1] =  0.733870426238362032891332767175E+00;
    weight[6-1] =  0.725590596901489156295739839779E+00;
    weight[7-1] =  0.394212014211504966587433032679E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = -0.857017929919813794402037235698E+00;
    xtab[2-1] = -0.631543407166567521509503573952E+00;
    xtab[3-1] = -0.339104543648722903660229021109E+00;
    xtab[4-1] = -0.0111941563689783438801237300122E+00;
    xtab[5-1] =  0.316696017045595559454075475675E+00;
    xtab[6-1] =  0.609049663022520165351466780939E+00;
    xtab[7-1] =  0.834198765028697794599267293239E+00;
    xtab[8-1] =  0.967804480896157932935972899807E+00;

    weight[1-1] =  0.00374814227227757804631954025851E+00;
    weight[2-1] =  0.0357961737041152639660521680263E+00;
    weight[3-1] =  0.137974910241879862433949246199E+00;
    weight[4-1] =  0.326515411108352185491692769217E+00;
    weight[5-1] =  0.547577467373226177976217604887E+00;
    weight[6-1] =  0.682278153375510121675529810121E+00;
    weight[7-1] =  0.614544746137780998436053880546E+00;
    weight[8-1] =  0.318231662453524478640851647411E+00;
  }
  else if ( order == 9 )
  {
    xtab[1-1] = -0.882491728426548422828684254270E+00;
    xtab[2-1] = -0.694873684026474640346360850039E+00;
    xtab[3-1] = -0.446537143480670863635920316400E+00;
    xtab[4-1] = -0.159388112702326252531544826624E+00;
    xtab[5-1] =  0.141092709224374414981503995427E+00;
    xtab[6-1] =  0.428217823321559204544020866175E+00;
    xtab[7-1] =  0.676480966471850715860378175342E+00;
    xtab[8-1] =  0.863830940812464825046988286026E+00;
    xtab[9-1] =  0.973668228805771018909618924364E+00;

    weight[1-1] =  0.00209009877215570354392734918986E+00;
    weight[2-1] =  0.0205951891648697848186537272448E+00;
    weight[3-1] =  0.0832489326348178964194106978875E+00;
    weight[4-1] =  0.210746247220398685903797568021E+00;
    weight[5-1] =  0.388325022916052063676224499399E+00;
    weight[6-1] =  0.554275165518437673725822282791E+00;
    weight[7-1] =  0.621388553284444032628761363828E+00;
    weight[8-1] =  0.523916296267173054255512857631E+00;
    weight[9-1] =  0.262081160888317771694556320674E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET_X2 - Fatal error!\n";
    cerr << "  Illegal input value of ORDER = " << order << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double lens_half_2d ( double func ( double x, double y ), double center[2], 
  double r, double theta1, double theta2, int order )

//****************************************************************************80
//
//  Purpose:
//
//    LENS_HALF_2D approximates an integral in a circular half lens in 2D.
//
//  Discussion:
//
//    A circular half lens is formed by drawing a circular arc,
//    and joining its endpoints.
//
//    This rule for a circular half lens simply views the region as 
//    a product region, with a coordinate "S" that varies along the
//    radial direction, and a coordinate "T" that varies in the perpendicular
//    direction, and whose extent varies as a function of S.
//
//    A Gauss-Legendre rule is used to construct a product rule that is
//    applied to the region.  The accuracy of the Gauss-Legendre rule,
//    which is valid for a rectangular product rule region, does not
//    apply straightforwardly to this region, since the limits in the
//    "T" coordinate are being handled implicitly.
//
//    This is simply an application of the QMULT_2D algorithm of Stroud.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double CENTER[2], the center of the circle.
//
//    Input, double R, the radius of the circle.
//
//    Input, double THETA1, THETA2, the angles of the rays
//    that begin and end the arc.
//
//    Input, int ORDER, the order of the Gauss-Legendre rule
//    to be used.  Legal values include 1 through 20, 32 or 64.
//
//    Output, double LENS_HALF_2D, the approximate value
//    of the integral of the function over the half lens.
//
{
  double ax;
  double ay;
  double bx;
  double by;
  double cx;
  double cy;
  double dx;
  double dy;
  int i;
  int j;
  double quad;
  double s_length;
  double sx;
  double sy;
  double t_length;
  double tdirx;
  double tdiry;
  double thi;
  double tx;
  double ty;
  double w1;
  double w2;
  double *weight;
  double *xtab;
//
//  Determine the points A (on the secant) and B (on the circumference)
//  that will form the "S" direction.
//
  ax = center[0] + r * 0.5 * ( cos ( theta1 ) + cos ( theta2 ) );
  ay = center[1] + r * 0.5 * ( sin ( theta1 ) + sin ( theta2 ) );

  bx = center[0] + r * cos ( 0.5 * ( theta1 + theta2 ) );
  by = center[1] + r * sin ( 0.5 * ( theta1 + theta2 ) );
//
//  Find the length of the line between A and B.
//
  s_length = sqrt ( pow ( ax - bx, 2 ) + pow ( ay - by, 2 ) );

  if ( s_length == 0.0 )
  {
    quad = 0.0;
    return quad;
  }
//
//  Retrieve the Legendre rule of the given order.
//
  xtab = new double[order];
  weight = new double[order];

  legendre_set ( order, xtab, weight );
//
//  Determine the unit vector in the T direction.
//
  tdirx = ( ay - by ) / s_length;
  tdiry = ( bx - ax ) / s_length;

  quad = 0.0;

  for ( i = 0; i < order; i++ )
  {
    w1 = 0.5 * s_length * weight[i];
//
//  Map the quadrature point to an S coordinate.
//
    sx = ( ( 1.0 - xtab[i] ) * ax   
         + ( 1.0 + xtab[i] ) * bx ) 
         /   2.0;
    sy = ( ( 1.0 - xtab[i] ) * ay   
         + ( 1.0 + xtab[i] ) * by ) 
         /   2.0;
//
//  Determine the length of the line in the T direction, from the
//  S axis to the circle circumference.
//
    thi = sqrt ( ( r - 0.25 * ( 1.0 - xtab[i] ) * s_length ) 
                            * ( 1.0 - xtab[i] ) * s_length );
// 
//  Determine the maximum and minimum T coordinates by going
//  up and down in the T direction from the S axis.
//
    cx = sx + tdirx * thi;
    cy = sy + tdiry * thi;
    dx = sx - tdirx * thi;
    dy = sy - tdiry * thi;
//
//  Find the length of the T direction.
//
    t_length = sqrt ( pow ( cx - dx, 2 ) + pow ( cy - dy, 2 ) );

    for ( j = 0; j < order; j++ )
    {
      w2 = 0.5 * t_length * weight[j];
//
//  Map the quadrature point to a T coordinate.
//
      tx = ( ( 1.0 - xtab[j] ) * cx   
           + ( 1.0 + xtab[j] ) * dx ) 
           /   2.0;
      ty = ( ( 1.0 - xtab[j] ) * cy   
           + ( 1.0 + xtab[j] ) * dy ) 
           /   2.0;

      quad = quad + w1 * w2 * func ( tx, ty );
    }
  }

  delete [] xtab;
  delete [] weight;

  return quad;
}
//****************************************************************************80

double lens_half_area_2d ( double r, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    LENS_HALF_AREA_2D returns the area of a circular half lens in 2D.
//
//  Discussion:
//
//    A circular half lens is formed by drawing a circular arc, 
//    and joining its endpoints.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Input, double THETA1, THETA2, the angles of the rays
//    that begin and end the arc.
//
//    Output, double LENS_HALF_AREA_2D, the area of the half lens.
//
{
  double sector;
  double triangle;
  double value;

  sector = circle_sector_area_2d ( r, theta1, theta2 );
  triangle = circle_triangle_area_2d ( r, theta1, theta2 );
  value = sector - triangle;

  return value;
}
//****************************************************************************80

double lens_half_h_area_2d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    LENS_HALF_H_AREA_2D returns the area of a circular half lens in 2D.
//
//  Discussion:
//
//    A circular half lens is formed by drawing a circular arc, and joining 
//    its endpoints.
//
//    This particular half lens is described by the "height" of the region.  
//    In other words, the half lens is the region that would be submerged 
//    if a circle of radius R were standing in water of depth H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Input, double H, the height of the half lens region.
//
//    Output, double LENS_HALF_H_AREA_2D, the area of the half lens.
//
{
  double angle;
  double area;
  double half_width;
  double pi = 3.141592653589793;
  double sector;
  double triangle;

  if ( h <= 0.0 )
  {
    area = 0.0;
  }
  else if ( 2.0 * r <= h )
  {
    area = pi * r * r;
  }
  else
  {
    half_width = sqrt ( h * ( 2.0 * r - h ) );
    angle = 2.0 * atan2 ( half_width, r - h );
    sector = r * r * angle / 2.0;
    triangle = ( r - h ) * half_width;
    area = sector - triangle;
  }
  return area;
}
//****************************************************************************80

double lens_half_w_area_2d ( double r, double w )

//****************************************************************************80
//
//  Purpose:
//
//    LENS_HALF_W_AREA_2D returns the area of a circular half lens in 2D.
//
//  Discussion:
//
//    A half lens is formed by drawing a circular arc, and joining its endpoints.
//    This half lens is described by the "width" of the region.  In other words,
//    it is the portion of the circle under water if the width
//    of the water surface is W.  There are two possible values for this
//    area, A and (PI*R*R-A).  The routine returns the smaller of the 
//    two values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//
//    Input, double W, the width of the half lens region.
//
//    Output, double LENS_HALF_W_AREA_2D, the area of the half lens.
//
{
  double angle;
  double area;
  double h;
  double half_width;
  double pi = 3.141592653589793;
  double sector;
  double triangle;

  if ( w <= 0.0 )
  {
    area = 0.0;
  }
  else if ( 2.0 * r <= w ) 
  {
    area = 0.5 * pi * r * r;
  }
  else
  {
    half_width = 0.5 * w;
    h = r - sqrt ( r * r - half_width * half_width );
    angle = 2.0 * atan2 ( half_width, r - h );
    sector = r * r * angle / 2.0;
    triangle = ( r - h ) * half_width;
    area = sector - triangle;
  }
  return area;
}
//****************************************************************************80

double *monomial_value ( int dim_num, int point_num, double x[], int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, double X[DIM_NUM*POINT_NUM], the point coordinates.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Output, double MONOMIAL_VALUE[POINT_NUM], the value of the monomial.
//
{
  int dim;
  int point;
  double *value;

  value = new double[point_num];

  for ( point = 0; point < point_num; point++ )
  {
    value[point] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0 != expon[dim] )
    {
      for ( point = 0; point < point_num; point++ )
      {
        value[point] = value[point] * pow ( x[dim+point*dim_num], expon[dim] );
      }
    }
  }

  return value;
}
//****************************************************************************80

double octahedron_unit_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    OCTAHEDRON_UNIT_ND approximates integrals in the unit octahedron in ND.
//
//  Integration region:
//
//    sum ( abs ( X(1:N) ) ) <= 1.
//
//  Discussion:
//
//    A 2*N point 3rd degree formula is used, Stroud number GN:3-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the octahedron.
//
//    Output, double OCTAHEDRON_UNIT_ND, the approximate integral of the function.
//
{
  int i;
  int j;
  double quad;
  double r;
  double result;
  double volume;
  double w;
  double *x;

  x = new double[n];

  w = 1.0 / ( double ) ( 2 * n );

  r = sqrt ( ( double ) ( 2 * n ) 
    / ( double ) ( ( n + 1 ) * ( n + 2 ) ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  quad = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x[i] = r;
    for ( j = 0; j < 2; j++ )
    {
      quad = quad + w * func ( n, x );
      x[i] = - x[i];
    }
    x[i] = 0.0;
  }

  volume = octahedron_unit_volume_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}

//****************************************************************************80

double octahedron_unit_volume_nd ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    OCTAHEDRON_UNIT_VOLUME_ND returns the volume of the unit octahedron in ND.
//
//  Integration region:
//
//    sum ( abs ( X(1:N) ) ) <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Output, double OCTAHEDRON_UNIT_VOLUME_ND, the volume of
//    the unit octahedron.
//
{
  int i;
  double volume;

  volume = 1.0;
  for ( i = 1; i <= n; i++ )
  {
    volume = volume * 2.0 / ( double ) ( i );
  }

  return volume;
}
//****************************************************************************80

double parallelipiped_volume_3d ( double x[4], double y[4], double z[4] )

//****************************************************************************80
//
//  Purpose:
//
//    PARALLELIPIPED_VOLUME_3D returns the volume of a parallelipiped in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X[4], Y[4], Z[4], the coordinates of one corner
//    of the parallelipiped, and its 3 immediate neighbors.
//
//    Output, double PARALLELIPIPED_VOLUME_3D, the volume of
//    the parallelipiped.
//
{
  double volume;

  volume = r8_abs ( 
    ( z[1] - z[0] ) * ( y[3] * x[2] - y[2] * x[3] ) + 
    ( z[2] - z[0] ) * ( x[3] * y[1] - x[1] * y[3] ) + 
    ( z[3] - z[0] ) * ( x[1] * y[2] - x[2] * y[1] ) + 
    ( z[2] - z[1] ) * ( y[3] * x[0] - y[0] * x[3] ) + 
    ( z[3] - z[1] ) * ( x[2] * y[0] - x[0] * y[2] ) + 
    ( z[3] - z[2] ) * ( x[0] * y[1] - x[1] * y[0] ) );

  return volume;
}
//****************************************************************************80

double parallelipiped_volume_nd ( int n, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PARALLELIPIPED_VOLUME_ND returns the volume of a parallelipiped in ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double V[N*(N+1)], the N+1 columns of V contains the 
//    coordinates of the "corners" of the parallelipiped.
//
//    Output, double PARALLELIPIPED_VOLUME_ND, the volume of
//    the parallelipiped.
//
{
  double det;
  int i;
  int info;
  int j;
  int *pivot;
  double volume;
  double *w;
//
//  Compute the volume of the N-dimensional parallelipiped.
//
  w = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      w[i+j*n] = v[i+(j+1)*n] - v[i+0*n];
    }
  }

  pivot = new int[n];

  info = r8ge_fa ( n, w, pivot );

  if ( info != 0 )
  {
    volume = 0.0;
  }
  else
  {
    det = r8ge_det ( n, w, pivot );

    volume = r8_abs ( det );
  }

  delete [] pivot;
  delete [] w;

  return volume;
}
//****************************************************************************80

double polygon_1_2d ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_1_2D integrates the function 1 over a polygon in 2D.
//
//  Integration region:
//
//    The polygon bounded by the points (X(1:N), Y(1:N)).
//
//  Formula:
//
//    INTEGRAL = 0.5 * sum ( 1 <= I <= N )
//      ( X(I) + X(I-1) ) * ( Y(I) - Y(I-1) )
//
//    where X(0) and Y(0) should be replaced by X(N) and Y(N).
//
//  Discussion:
//
//    The integral of 1 over a polygon is the area of the polygon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    SF Bockman,
//    Generalizing the Formula for Areas of Polygons to Moments,
//    American Mathematical Society Monthly,
//    1989, pages 131-132.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//    N should be at least 3 for a nonzero result.
//
//    Input, double X[N], Y[N], the coordinates of the vertices
//    of the polygon.  These vertices should be given in counter-clockwise order.
//
//    Output, double POLYGON_1_2D, the value of the integral.
//
{
  int i;
  int im1;
  double result;

  result = 0.0;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_1_2D - Fatal error!\n";
    cerr << "  The number of vertices must be at least 3.\n";
    cerr << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }
  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      im1 = n - 1;
    }
    else
    {
      im1 = i - 1;
    }
    result = result + 0.5 * ( x[i] + x[im1] ) * ( y[i] - y[im1] );
  }

  return result;
}
//****************************************************************************80

double polygon_x_2d ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_X_2D integrates the function X over a polygon in 2D.
//
//  Integration region:
//
//    The polygon bounded by the points (X(1:N), Y(1:N)).
//
//  Formula:
//
//    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
//      ( X(I)^2 + X(I) * X(I-1) + X(I-1)^2 ) * ( Y(I) - Y(I-1) )
//
//    where X(0) and Y(0) should be replaced by X(N) and Y(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    SF Bockman,
//    Generalizing the Formula for Areas of Polygons to Moments,
//    American Mathematical Society Monthly,
//    1989, pages 131-132.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//    N should be at least 3 for a nonzero result.
//
//    Input, double X[N], Y[N], the coordinates of the vertices
//    of the polygon.  These vertices should be given in counter-clockwise order.
//
//    Output, double POLYGON_X_2D, the value of the integral.
//
{
  int i;
  int im1;
  double result;

  result = 0.0;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_X_2D - Fatal error!\n";
    cerr << "  The number of vertices must be at least 3.\n";
    cerr << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      im1 = n - 1;
    }
    else
    {
      im1 = i - 1;
    }
    result = result + ( x[i] * x[i] + x[i] * x[im1] + x[im1] * x[im1] ) 
      * ( y[i] - y[im1] );
  }

  result = result / 6.0;

  return result;
}
//****************************************************************************80

double polygon_xx_2d ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
//
//  Integration region:
//
//    The polygon bounded by the points (X(1:N), Y(1:N)).
//
//  Formula:
//
//    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
//      ( X(I)^3 + X(I)^2 * X(I-1) + X(I) * X(I-1)^2 + X(I-1)^3 )
//      * ( Y(I) - Y(I-1) )
//
//    where X(0) and Y(0) should be replaced by X(N) and Y(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    SF Bockman,
//    Generalizing the Formula for Areas of Polygons to Moments,
//    American Mathematical Society Monthly,
//    Volume 96, Number 2, February 1989, pages 131-132.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//    N should be at least 3 for a nonzero result.
//
//    Input, double X[N], Y[N], the coordinates of the vertices
//    of the polygon.  These vertices should be given in
//    counter-clockwise order.
//
//    Output, double RESULT, the value of the integral.
//
{
  int i;
  int im1;
  double result;

  result = 0.0;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_XX_2D - Fatal error!\n";
    cerr << "  The number of vertices must be at least 3.\n";
    cerr << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      im1 = n - 1;
    }
    else
    {
      im1 = i - 1;
    }

    result = result + ( x[i] * x[i] * x[i] + x[i] * x[i] * x[im1] 
      + x[i] * x[im1] * x[im1] + x[im1] * x[im1] * x[im1] ) 
      * ( y[i] - y[im1] );
  }

  result = result / 12.0;

  return result;
}
//****************************************************************************80

double polygon_xy_2d ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
//
//  Integration region:
//
//    The polygon bounded by the points (X(1:N), Y(1:N)).
//
//  Formula:
//
//    INTEGRAL = (1/24) * sum ( 1 <= I <= N )
//      ( Y(I)   * ( 3 * X(I)**2 + 2 * X(I) * X(I-1) +     X(I-1)**2 )
//      + Y(I-1) * (     X(I)**2 + 2 * X(I) * X(I-1) + 3 * X(I-1)**2 ) )
//      * ( Y(I) - Y(I-1) )
//
//    where X(0) and Y(0) should be replaced by X(N) and Y(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    SF Bockman,
//    Generalizing the Formula for Areas of Polygons to Moments,
//    American Mathematical Society Monthly,
//    Volume 96, Number 2, February 1989, pages 131-132.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//    N should be at least 3 for a nonzero result.
//
//    Input, double X[N], Y[N], the coordinates of the vertices
//    of the polygon.  These vertices should be given in
//    counter-clockwise order.
//
//    Output, double POLYGON_XY_2D, the value of the integral.
//
{
  int i;
  int im1;
  double result;

  result = 0.0;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_XY_2D - Fatal error!\n";
    cerr << "  The number of vertices must be at least 3.\n";
    cerr << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      im1 = n - 1;
    }
    else
    {
      im1 = i - 1;
    }

    result = result + ( 
      y[i] * ( 3.0 * x[i] * x[i] + 2.0 * x[i] * x[im1] + x[im1] * x[im1] ) 
      + y[im1] * ( x[i] * x[i] + 2.0 * x[i] * x[im1] + 3.0 * x[im1] * x[im1] ) 
      ) * ( y[i] - y[im1] );
  }

  result = result / 24.0;

  return result;
}
//****************************************************************************80

double polygon_y_2d ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_Y_2D integrates the function Y over a polygon in 2D.
//
//  Integration region:
//
//    The polygon bounded by the points (X(1:N), Y(1:N)).
//
//  Formula:
//
//    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
//      - ( Y(I)^2 + Y(I) * Y(I-1) + Y(I-1)^2 ) * ( X(I) - X(I-1) )
//
//    where X(0) and Y(0) should be replaced by X(N) and Y(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    SF Bockman,
//    Generalizing the Formula for Areas of Polygons to Moments,
//    American Mathematical Society Monthly,
//    Volume 96, Number 2, February 1989, pages 131-132.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//    N should be at least 3 for a nonzero result.
//
//    Input, double X[N], Y[N], the coordinates of the vertices
//    of the polygon.  These vertices should be given in
//    counter-clockwise order.
//
//    Output, double POLYGON_Y_2D, the value of the integral.
//
{
  int i;
  int im1;
  double result;

  result = 0.0;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_Y_2D - Fatal error!\n";
    cerr << "  The number of vertices must be at least 3.\n";
    cerr << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      im1 = n - 1;
    }
    else
    {
      im1 = i - 1;
    }

    result = result - ( y[i] * y[i] + y[i] * y[im1] + y[im1] * y[im1] ) 
      * ( x[i] - x[im1] );
  }

  result = result / 6.0;

  return result;
}
//****************************************************************************80

double polygon_yy_2d ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
//
//  Integration region:
//
//    The polygon bounded by the points (X(1:N), Y(1:N)).
//
//  Formula:
//
//    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
//      - ( Y(I)^3 + Y(I)^2 * Y(I-1) + Y(I) * Y(I-1)^2 + Y(I-1)^3 )
//      * ( X(I) - X(I-1) )
//
//    where X(0) and Y(0) should be replaced by X(N) and Y(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    SF Bockman,
//    Generalizing the Formula for Areas of Polygons to Moments,
//    American Mathematical Society Monthly,
//    Volume 96, Number 2, February 1989, pages 131-132.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//    N should be at least 3 for a nonzero result.
//
//    Input, double X[N], Y[N], the coordinates of the vertices
//    of the polygon.  These vertices should be given in
//    counter-clockwise order.
//
//    Output, double POLYGON_YY_2D, the value of the integral.
//
{
  int i;
  int im1;
  double result;

  result = 0.0;

  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_YY_2D - Fatal error!\n";
    cerr << "  The number of polygonal vertices must be\n";
    cerr << "  at least 3, but the input polygon has N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      im1 = n - 1;
    }
    else
    {
      im1 = i - 1;
    }

    result = result - ( 
        y[i] * y[i] * y[i] 
      + y[i] * y[i] * y[im1] 
      + y[i] * y[im1] * y[im1] 
      + y[im1] * y[im1] * y[im1]
    ) * ( x[i] - x[im1] );
  }

  result = result / 12.0;

  return result;
}
//****************************************************************************80

double pyramid_unit_o01_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O01_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    A 1 point degree 1 formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  double quad;
  double result;
  double volume;
  double w;
  double x;
  double y;
  double z;
//
//  Quadrature.
//
  quad = 0.0;

  x = 0.0;
  y = 0.0;
  z = 1.0 / 4.0;
  w = 1.0;

  quad = quad + w * func ( x, y, z );
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o05_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O05_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    A 5 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O05_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 5;
  double quad;
  double result;
  double volume;
  double w[5] = {
   0.21093750000000000000, 
   0.21093750000000000000, 
   0.21093750000000000000, 
   0.21093750000000000000, 
   0.15625000000000000000 };
  double x[5] = {
  -0.48686449556014765641, 
   0.48686449556014765641, 
   0.48686449556014765641, 
  -0.48686449556014765641, 
   0.00000000000000000000 };
  double y[5] = {
  -0.48686449556014765641, 
  -0.48686449556014765641, 
   0.48686449556014765641, 
   0.48686449556014765641, 
   0.00000000000000000000 };
  double z[5] = {
   0.16666666666666666667, 
   0.16666666666666666667, 
   0.16666666666666666667, 
   0.16666666666666666667, 
   0.70000000000000000000 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o06_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O06_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    A 6 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O06_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 6;
  double quad;
  double result;
  double volume;
  double w[6] = {
   0.21000000000000000000,
   0.21000000000000000000,
   0.21000000000000000000,
   0.21000000000000000000,
   0.06000000000000000000,
   0.10000000000000000000 };
  double x[6] = {
  -0.48795003647426658968,
   0.48795003647426658968,
   0.48795003647426658968,
  -0.48795003647426658968,
   0.00000000000000000000,
   0.00000000000000000000 };
  double y[6] = {
  -0.48795003647426658968,
  -0.48795003647426658968,
   0.48795003647426658968,
   0.48795003647426658968,
   0.00000000000000000000,
   0.00000000000000000000 };
  double z[6] = {
   0.16666666666666666667,
   0.16666666666666666667,
   0.16666666666666666667,
   0.16666666666666666667,
   0.58333333333333333333,
   0.75000000000000000000 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o08_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O08_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    An 8 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O08_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 8;
  double quad;
  double result;
  double volume;
  double w[8] = {
   0.075589411559869072938,
   0.075589411559869072938,
   0.075589411559869072938,
   0.075589411559869072938,
   0.17441058844013092706,
   0.17441058844013092706,
   0.17441058844013092706,
   0.17441058844013092706 };
  double x[8] = {
  -0.26318405556971359557,
   0.26318405556971359557,
   0.26318405556971359557,
  -0.26318405556971359557,
  -0.50661630334978742377,
   0.50661630334978742377,
   0.50661630334978742377,
  -0.50661630334978742377 };
  double y[8] = {
  -0.26318405556971359557,
  -0.26318405556971359557,
   0.26318405556971359557,
   0.26318405556971359557,
  -0.50661630334978742377,
  -0.50661630334978742377,
   0.50661630334978742377,
   0.50661630334978742377 };
  double z[8] = {
   0.54415184401122528880,
   0.54415184401122528880,
   0.54415184401122528880,
   0.54415184401122528880,
   0.12251482265544137787,
   0.12251482265544137787,
   0.12251482265544137787,
   0.12251482265544137787 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o08b_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O08B_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    An 8 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O08B_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 8;
  double quad;
  double result;
  double volume;
  double w[8] = {
   0.16438287736328777572,
   0.16438287736328777572,
   0.16438287736328777572,
   0.16438287736328777572,
   0.085617122636712224276,
   0.085617122636712224276,
   0.085617122636712224276,
   0.085617122636712224276 };
  double x[8] = {
  -0.51197009372656270107,
   0.51197009372656270107,
   0.51197009372656270107,
  -0.51197009372656270107,
  -0.28415447557052037456,
   0.28415447557052037456,
   0.28415447557052037456,
  -0.28415447557052037456 };
  double y[8] = {
  -0.51197009372656270107,
  -0.51197009372656270107,
   0.51197009372656270107,
   0.51197009372656270107,
  -0.28415447557052037456,
  -0.28415447557052037456,
   0.28415447557052037456,
   0.28415447557052037456 };
  double z[8] = {
   0.11024490204163285720,
   0.11024490204163285720,
   0.11024490204163285720,
   0.11024490204163285720,
   0.518326526529795714229,
   0.518326526529795714229,
   0.518326526529795714229,
   0.518326526529795714229 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o09_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O09_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    A 9 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O09_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 9;
  double quad;
  double result;
  double volume;
  double w[9] = {
   0.13073389672275944791,
   0.13073389672275944791,
   0.13073389672275944791,
   0.13073389672275944791,
   0.10989110327724055209,
   0.10989110327724055209,
   0.10989110327724055209,
   0.10989110327724055209,
   0.03750000000000000000 };
  double x[9] = {
  -0.52966422253852215131,
   0.52966422253852215131,
   0.52966422253852215131,
  -0.52966422253852215131,
  -0.34819753825720418039,
   0.34819753825720418039,
   0.34819753825720418039,
  -0.34819753825720418039,
   0.00000000000000000000 };
  double y[9] = {
  -0.52966422253852215131,
  -0.52966422253852215131,
   0.52966422253852215131,
   0.52966422253852215131,
  -0.34819753825720418039,
  -0.34819753825720418039,
   0.34819753825720418039,
   0.34819753825720418039,
   0.00000000000000000000 };
  double z[9] = {
   0.08176876558246862335,
   0.08176876558246862335,
   0.08176876558246862335,
   0.08176876558246862335,
   0.400374091560388519511,
   0.400374091560388519511,
   0.400374091560388519511,
   0.400374091560388519511,
   0.83333333333333333333 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o13_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O13_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    A 13 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O13_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 13;
  double quad;
  double result;
  double volume;
  double w[13] = {
   0.063061594202898550725,
   0.063061594202898550725,
   0.063061594202898550725,
   0.063061594202898550725,
   0.042101946815575556199,
   0.042101946815575556199,
   0.042101946815575556199,
   0.042101946815575556199,
   0.13172030707666776585,
   0.13172030707666776585,
   0.13172030707666776585,
   0.13172030707666776585,
   0.05246460761943250889 };
  double x[13] = {
  -0.38510399211870384331,
   0.38510399211870384331,
   0.38510399211870384331,
  -0.38510399211870384331,
  -0.40345831960728204766,
   0.40345831960728204766,
   0.00000000000000000000,
   0.00000000000000000000,
  -0.53157877436961973359,
   0.53157877436961973359,
   0.53157877436961973359,
  -0.53157877436961973359,
   0.00000000000000000000 };
  double y[13] = {
  -0.38510399211870384331,
  -0.38510399211870384331,
   0.38510399211870384331,
   0.38510399211870384331,
   0.00000000000000000000,
   0.00000000000000000000,
  -0.40345831960728204766,
   0.40345831960728204766,
  -0.53157877436961973359,
  -0.53157877436961973359,
   0.53157877436961973359,
   0.53157877436961973359,
   0.00000000000000000000 };
  double z[13] = {
  0.428571428571428571429,
  0.428571428571428571429,
  0.428571428571428571429,
  0.428571428571428571429,
  0.33928571428571428571,
  0.33928571428571428571,
  0.33928571428571428571,
  0.33928571428571428571,
  0.08496732026143790850,
  0.08496732026143790850,
  0.08496732026143790850,
  0.08496732026143790850,
  0.76219701803768503595 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o18_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O18_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    An 18 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O18_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 18;
  double quad;
  double result;
  double volume;
  double w[18] = {
   0.023330065296255886709,
   0.037328104474009418735,
   0.023330065296255886709,
   0.037328104474009418735,
   0.059724967158415069975,
   0.037328104474009418735,
   0.023330065296255886709,
   0.037328104474009418735,
   0.023330065296255886709,
   0.05383042853090460712,
   0.08612868564944737139,
   0.05383042853090460712,
   0.08612868564944737139,
   0.13780589703911579422,
   0.08612868564944737139,
   0.05383042853090460712,
   0.08612868564944737139,
   0.05383042853090460712 };
  double x[18] = {
  -0.35309846330877704481,
   0.00000000000000000000,
   0.35309846330877704481,
  -0.35309846330877704481,
   0.00000000000000000000,
   0.35309846330877704481,
  -0.35309846330877704481,
   0.00000000000000000000,
   0.35309846330877704481,
  -0.67969709567986745790,
   0.00000000000000000000,
   0.67969709567986745790,
  -0.67969709567986745790,
   0.00000000000000000000,
   0.67969709567986745790,
  -0.67969709567986745790,
   0.00000000000000000000,
   0.67969709567986745790 };
  double y[18] = {
  -0.35309846330877704481,
  -0.35309846330877704481,
  -0.35309846330877704481,
   0.00000000000000000000,
   0.00000000000000000000,
   0.00000000000000000000,
   0.35309846330877704481,
   0.35309846330877704481,
   0.35309846330877704481,
  -0.67969709567986745790,
  -0.67969709567986745790,
  -0.67969709567986745790,
   0.00000000000000000000,
   0.00000000000000000000,
   0.00000000000000000000,
   0.67969709567986745790,
   0.67969709567986745790,
   0.67969709567986745790 };
  double z[18] = {
  0.544151844011225288800,
  0.544151844011225288800,
  0.544151844011225288800,
  0.544151844011225288800,
  0.544151844011225288800,
  0.544151844011225288800,
  0.544151844011225288800,
  0.544151844011225288800,
  0.544151844011225288800,
  0.12251482265544137787,
  0.12251482265544137787,
  0.12251482265544137787,
  0.12251482265544137787,
  0.12251482265544137787,
  0.12251482265544137787,
  0.12251482265544137787,
  0.12251482265544137787,
  0.12251482265544137787 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o27_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O27_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    A 27 point formula is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O27_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 27;
  double quad;
  double result;
  double volume;
  double w[27] = {
   0.036374157653908938268,
   0.05819865224625430123,
   0.036374157653908938268,
   0.05819865224625430123,
   0.09311784359400688197,
   0.05819865224625430123,
   0.036374157653908938268,
   0.05819865224625430123,
   0.036374157653908938268,
   0.033853303069413431019,
   0.054165284911061489631,
   0.033853303069413431019,
   0.054165284911061489631,
   0.08666445585769838341,
   0.054165284911061489631,
   0.033853303069413431019,
   0.054165284911061489631,
   0.033853303069413431019,
   0.006933033103838124540,
   0.011092852966140999264,
   0.006933033103838124540,
   0.011092852966140999264,
   0.017748564745825598822,
   0.011092852966140999264,
   0.006933033103838124540,
   0.011092852966140999264,
   0.006933033103838124540 };
  double x[27] = {
  -0.7180557413198889387,
   0.00000000000000000000,
   0.7180557413198889387,
  -0.7180557413198889387,
   0.00000000000000000000,
   0.7180557413198889387,
  -0.7180557413198889387,
   0.00000000000000000000,
   0.7180557413198889387,
  -0.50580870785392503961,
   0.00000000000000000000,
   0.50580870785392503961,
  -0.50580870785392503961,
   0.00000000000000000000,
   0.50580870785392503961,
  -0.50580870785392503961,
   0.00000000000000000000,
   0.50580870785392503961,
  -0.22850430565396735360,
   0.00000000000000000000,
   0.22850430565396735360,
  -0.22850430565396735360,
   0.00000000000000000000,
   0.22850430565396735360,
  -0.22850430565396735360,
   0.00000000000000000000,
   0.22850430565396735360 };
  double y[27] = {
  -0.7180557413198889387,
  -0.7180557413198889387,
  -0.7180557413198889387,
   0.00000000000000000000,
   0.00000000000000000000,
   0.00000000000000000000,
   0.7180557413198889387,
   0.7180557413198889387,
   0.7180557413198889387,
  -0.50580870785392503961,
  -0.50580870785392503961,
  -0.50580870785392503961,
   0.00000000000000000000,
   0.00000000000000000000,
   0.00000000000000000000,
   0.50580870785392503961,
   0.50580870785392503961,
   0.50580870785392503961,
  -0.22850430565396735360,
  -0.22850430565396735360,
  -0.22850430565396735360,
   0.00000000000000000000,
   0.00000000000000000000,
   0.00000000000000000000,
   0.22850430565396735360,
   0.22850430565396735360,
   0.22850430565396735360 };
  double z[27] = {
  0.07299402407314973216,
  0.07299402407314973216,
  0.07299402407314973216,
  0.07299402407314973216,
  0.07299402407314973216,
  0.07299402407314973216,
  0.07299402407314973216, 
  0.07299402407314973216,
  0.07299402407314973216,
  0.34700376603835188472,
  0.34700376603835188472,
  0.34700376603835188472,
  0.34700376603835188472,
  0.34700376603835188472,
  0.34700376603835188472,
  0.34700376603835188472,
  0.34700376603835188472,
  0.34700376603835188472,
  0.70500220988849838312,
  0.70500220988849838312,
  0.70500220988849838312,
  0.70500220988849838312,
  0.70500220988849838312,
  0.70500220988849838312,
  0.70500220988849838312,
  0.70500220988849838312,
  0.70500220988849838312 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_o48_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_O48_3D approximates an integral inside the unit pyramid in 3D.
//
//  Discussion:
//
//    An 48 point degree 7 formula, Stroud CN:C2:7-1, is used.
//
//    The (X,Y,Z) integration region can be represented as:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function which is to be integrated.
//
//    Output, double PYRAMID_UNIT_O48_3D, the approximate integral
//    of the function.
//
{
  int i;
  int order = 48;
  double quad;
  double result;
  double volume;
  double w[48] = {
  2.01241939442682455E-002, 
  2.01241939442682455E-002, 
  2.01241939442682455E-002, 
  2.01241939442682455E-002, 
  2.60351137043010779E-002, 
  2.60351137043010779E-002, 
  2.60351137043010779E-002, 
  2.60351137043010779E-002, 
  1.24557795239745531E-002, 
  1.24557795239745531E-002, 
  1.24557795239745531E-002, 
  1.24557795239745531E-002, 
  1.87873998794808156E-003, 
  1.87873998794808156E-003, 
  1.87873998794808156E-003, 
  1.87873998794808156E-003, 
  4.32957927807745280E-002, 
  4.32957927807745280E-002, 
  4.32957927807745280E-002, 
  4.32957927807745280E-002, 
  1.97463249834127288E-002, 
  1.97463249834127288E-002, 
  1.97463249834127288E-002, 
  1.97463249834127288E-002, 
  5.60127223523590526E-002, 
  5.60127223523590526E-002, 
  5.60127223523590526E-002, 
  5.60127223523590526E-002, 
  2.55462562927473852E-002, 
  2.55462562927473852E-002, 
  2.55462562927473852E-002, 
  2.55462562927473852E-002, 
  2.67977366291788643E-002, 
  2.67977366291788643E-002, 
  2.67977366291788643E-002, 
  2.67977366291788643E-002, 
  1.22218992265373354E-002, 
  1.22218992265373354E-002, 
  1.22218992265373354E-002, 
  1.22218992265373354E-002, 
  4.04197740453215038E-003, 
  4.04197740453215038E-003, 
  4.04197740453215038E-003, 
  4.04197740453215038E-003, 
  1.84346316995826843E-003, 
  1.84346316995826843E-003, 
  1.84346316995826843E-003, 
  1.84346316995826843E-003 };
  double x[48] = {
  0.88091731624450909E+00,     
 -0.88091731624450909E+00,     
   0.0000000000000000E+00,     
   0.0000000000000000E+00,     
  0.70491874112648223E+00,     
 -0.70491874112648223E+00,     
   0.0000000000000000E+00,     
   0.0000000000000000E+00,     
  0.44712732143189760E+00,     
 -0.44712732143189760E+00,     
   0.0000000000000000E+00,     
   0.0000000000000000E+00,     
  0.18900486065123448E+00,     
 -0.18900486065123448E+00,     
   0.0000000000000000E+00,     
   0.0000000000000000E+00,     
  0.36209733410322176E+00,     
 -0.36209733410322176E+00,     
 -0.36209733410322176E+00,     
  0.36209733410322176E+00,     
  0.76688932060387538E+00,     
 -0.76688932060387538E+00,     
 -0.76688932060387538E+00,     
  0.76688932060387538E+00,     
  0.28975386476618070E+00,     
 -0.28975386476618070E+00,     
 -0.28975386476618070E+00,     
  0.28975386476618070E+00,     
  0.61367241226233160E+00,     
 -0.61367241226233160E+00,     
 -0.61367241226233160E+00,     
  0.61367241226233160E+00,     
  0.18378979287798017E+00,     
 -0.18378979287798017E+00,     
 -0.18378979287798017E+00,     
  0.18378979287798017E+00,     
  0.38925011625173161E+00,     
 -0.38925011625173161E+00,     
 -0.38925011625173161E+00,     
  0.38925011625173161E+00,     
  7.76896479525748113E-02, 
 -7.76896479525748113E-02, 
 -7.76896479525748113E-02, 
  7.76896479525748113E-02, 
  0.16453962988669860E+00,     
 -0.16453962988669860E+00,     
 -0.16453962988669860E+00,     
  0.16453962988669860E+00 };
  double y[48] = {
   0.0000000000000000E+00,      
   0.0000000000000000E+00,      
  0.88091731624450909E+00,      
 -0.88091731624450909E+00,      
   0.0000000000000000E+00,      
   0.0000000000000000E+00,      
  0.70491874112648223E+00,      
 -0.70491874112648223E+00,     
   0.0000000000000000E+00,      
   0.0000000000000000E+00,      
  0.44712732143189760E+00,      
 -0.44712732143189760E+00,      
   0.0000000000000000E+00,      
   0.0000000000000000E+00,      
  0.18900486065123448E+00,      
 -0.18900486065123448E+00,      
  0.36209733410322176E+00,      
  0.36209733410322176E+00,      
 -0.36209733410322176E+00,      
 -0.36209733410322176E+00,      
  0.76688932060387538E+00,      
  0.76688932060387538E+00,      
 -0.76688932060387538E+00,      
 -0.76688932060387538E+00,      
  0.28975386476618070E+00,      
  0.28975386476618070E+00,      
 -0.28975386476618070E+00,      
 -0.28975386476618070E+00,      
  0.61367241226233160E+00,      
  0.61367241226233160E+00,      
 -0.61367241226233160E+00,      
 -0.61367241226233160E+00,      
  0.18378979287798017E+00,      
  0.18378979287798017E+00,      
 -0.18378979287798017E+00,      
 -0.18378979287798017E+00,      
  0.38925011625173161E+00,      
  0.38925011625173161E+00,      
 -0.38925011625173161E+00,      
 -0.38925011625173161E+00,      
  7.76896479525748113E-02, 
  7.76896479525748113E-02, 
 -7.76896479525748113E-02, 
 -7.76896479525748113E-02, 
  0.16453962988669860E+00,      
  0.16453962988669860E+00,      
 -0.16453962988669860E+00, 
 -0.16453962988669860E+00 };
  double z[48] = {
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  4.85005494469969989E-02, 
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.23860073755186201E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.51704729510436798E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00,      
  0.79585141789677305E+00, 
  0.79585141789677305E+00 };
//
//  Quadrature.
//
  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + w[i] * func ( x[i], y[i], z[i] );
  }
//
//  Volume.
//
  volume = pyramid_unit_volume_3d ( );
//
//  Result.
//
  result = quad * volume;

  return result;
}
//****************************************************************************80

double pyramid_unit_monomial_3d ( int alpha, int beta, int gamma )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_MONOMIAL_3D: monomial integral in a unit pyramid in 3D.
//
//  Discussion:
//
//    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
//    over the unit pyramid.
//
//    The unit pyramid is defined as:
//
//    - ( 1 - Z ) <= X <= 1 - Z
//    - ( 1 - Z ) <= Y <= 1 - Z
//              0 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int ALPHA, BETA, GAMMA, the exponents of
//    X, Y and Z in the monomial.
//
//    Output, double PYRAMID_UNIT_MONOMIAL_3D, the volume of the pyramid.
//
{
  int i;
  int i_hi;
  double value;

  value = 0.0;

  if ( ( alpha % 2 ) == 0 && ( beta % 2 ) == 0 )
  {
    i_hi = 2 + alpha + beta;

    for ( i = 0; i <= i_hi; i++ )
    {
      value = value + r8_mop ( i ) * r8_choose ( i_hi, i ) 
      / ( double ) ( i + gamma + 1 );
    }

    value = value 
          * 2.0 / ( double ) ( alpha + 1 )
          * 2.0 / ( double ) ( beta + 1 );
  }

  return value;
}
//****************************************************************************80

double pyramid_unit_volume_3d ( )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_VOLUME_3D: volume of a unit pyramid with square base in 3D.
//
//  Integration region:
//
//    - ( 1 - Z ) <= X <= 1 - Z
//    - ( 1 - Z ) <= Y <= 1 - Z
//              0 <= Z <= 1.
//
//  Discussion:
//
//    The volume of this unit pyramid is 4/3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double PYRAMID_UNIT_VOLUME_3D, the volume of the pyramid.
//
{
  double volume;

  volume = 4.0 / 3.0;

  return volume;
}
//****************************************************************************80

double pyramid_volume_3d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_VOLUME_3D returns the volume of a pyramid with square base in 3D.
//
//  Integration region:
//
//    - ( H - Z ) * R <= X <= ( H - Z ) * R
//    - ( H - Z ) * R <= Y <= ( H - Z ) * R
//                  0 <= Z <= H.
//
//  Discussion:
//
//    A pyramid with square base can be regarded as the upper half of a
//    3D octahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the "radius" of the pyramid, that is, half the
//    length of one of the sides of the square base.
//
//    Input, double H, the height of the pyramid.
//
//    Output, double PYRAMID_VOLUME_3D, the volume of the pyramid.
//
{
  double value;

  value = ( 4.0 / 3.0 ) * h * r * r;

  return value;
}
//****************************************************************************80

double qmdpt ( double func ( int n, double x[] ), int n, int nsub )

//****************************************************************************80
//
//  Purpose:
//
//    QMDPT carries out product midpoint quadrature for the unit cube in ND.
//
//  Integration region:
//
//    -1 <= X(1:N) <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the user supplied
//    function to be integrated.
//
//    Input, int N, the dimension of the cube.
//
//    Input, int NSUB, the number of subdivisions (in each dimension).
//
//    Output, double QMDPT, the approximate integral of the function.
//
{
  int i;
  int ihi;
  int *ix;
  int j;
  bool more;
  double quad;
  double result;
  double volume;
  double w;
  double *x;

  w = 1.0 / ( double ) ( i4_power ( nsub, n ) );
  quad = 0.0;

  more = false;
  ihi = i4_power ( nsub, n );

  ix = new int[n];
  x = new double[n];

  for ( i = 0; i < ihi; i++ )
  {
    vec_lex_next ( n, nsub, ix, &more );

    for ( j = 0; j < n; j++ )
    {
      x[j] = ( double ) ( 2 * ix[j] + 1 - nsub ) / ( double ) ( nsub );
    }
    quad = quad + w * func ( n, x );
  }

  volume = pow ( 2.0, n );
  result = quad * volume;

  delete [] ix;
  delete [] x;

  return result;
}
//****************************************************************************80

double qmult_1d ( double func ( double x ), double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    QMULT_1D approximates an integral over an interval in 1D.
//
//  Integration region:
//
//    A <= X <= B.
//
//  Discussion:
//
//    A 16 point 31-st degree Gauss-Legendre formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x ), the name of the user supplied
//    function which evaluates F(X).
//
//    Input, double A, B, the lower and upper limits of integration.
//
//    Output, double QMULT_1D, the approximate integral of 
//    the function.
//
{
  int i;
  int order = 16;
  double quad;
  double result;
  double volume;
  double *weight;
  double x;
  double *xtab;

  xtab = new double[order];
  weight = new double[order];

  legendre_set ( order, xtab, weight );

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    x = 0.5 * ( b - a ) * xtab[i] + 0.5 * ( a + b );
    quad = quad + 0.5 * weight[i] * func ( x );
  }

  volume = b - a;
  result = quad * volume;

  delete [] xtab;
  delete [] weight;

  return result;
}
//****************************************************************************80

double qmult_2d ( double func ( double x, double y ), double a, double b, 
  double fup ( double x ), double flo ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    QMULT_2D approximates an integral with varying Y dimension in 2D.
//
//  Integration region:
//
//      A <= X <= B
//
//    and
//
//      FLO(X) <= Y <= FHI(X).
//
//  Discussion:
//
//    A 256 point product of two 16 point 31-st degree Gauss-Legendre
//    quadrature formulas is used.
//
//    This routine could easily be modified to use a different
//    order product rule by changing the value of ORDER.
//
//    Another easy change would allow the X and Y directions to
//    use quadrature rules of different orders.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function which evaluates F(X,Y).
//
//    Input, double A, B, the lower and upper limits of X integration.
//
//    Input, double FUP ( double x ), double FLO ( double x ), 
//    the names of the user supplied functions which evaluate the upper 
//    and lower limits of the Y integration.
//
{
  double c;
  double d;
  int i;
  int j;
  int order = 16;
  double quad;
  double w1;
  double w2;
  double *weight;
  double x;
  double *xtab;
  double y;

  xtab = new double[order];
  weight = new double[order];

  legendre_set ( order, xtab, weight );

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    w1 = 0.5 * ( b - a ) * weight[i];
    x = 0.5 * ( b - a ) * xtab[i] + 0.5 * ( b + a );
    c = flo ( x );
    d = fup ( x );

    for ( j = 0; j < order; j++ )
    {
      w2 = 0.5 * ( d - c ) * weight[j];
      y = 0.5 * ( d - c ) * xtab[j] + 0.5 * ( d + c );
      quad = quad + w1 * w2 * func ( x, y );
    }
  }

  delete [] xtab;
  delete [] weight;

  return quad;
}
//****************************************************************************80

double qmult_3d ( double func ( double x, double y, double z ), double a, 
  double b, double fup1 ( double x ), double flo1 ( double x ), 
  double fup2 ( double x, double y ), double flo2 ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    QMULT_3D approximates an integral with varying Y and Z dimension in 3D.
//
//  Integration region:
//
//      A         <= X <= B,
//    and
//      FLO(X)    <= Y <= FHI(X),
//    and
//      FLO2(X,Y) <= Z <= FHI2(X,Y).
//
//  Discussion:
//
//    A 4096 point product of three 16 point 31-st degree Gauss-Legendre
//    quadrature formulas is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied unction which evaluates F(X,Y,Z).
//
//    Input, double A, B, the lower and upper limits of X integration.
//
//    Input, double FUP1 ( double x ), double FLO1 ( double x ), the names 
//    of the user supplied functions which evaluate the upper and lower
//    limits of the Y integration.
//
//    Input, double FUP2 ( double x, double y ), 
//    double FLO2 ( double x, double y ), the names of the user
//    supplied functions which evaluate the upper and lower
//    limits of the Z integration.
//
//    Output, double QMULT_3D, the approximate integral of 
//    the function.
//
{
  double c;
  double d;
  double e;
  double f;
  int i;
  int j;
  int k;
  int order = 16;
  double quad;
  double result;
  double volume;
  double w1;
  double w2;
  double w3;
  double *weight;
  double x;
  double *xtab;
  double y;
  double z;

  xtab = new double[order];
  weight = new double[order];

  legendre_set ( order, xtab, weight );

  quad = 0.0;

  for ( i = 0; i < order; i++ )
  {
    x = 0.5 * ( b - a ) * xtab[i] + 0.5 * ( b + a );
    w1 = 0.5 * weight[i];
    c = flo1 ( x );
    d = fup1 ( x );

    for ( j = 0; j < order; j++ )
    {
      w2 = 0.5 * ( d - c ) * weight[j];
      y = 0.5 * ( d - c ) * xtab[j] + 0.5 * ( d + c );
      e = flo2 ( x, y );
      f = fup2 ( x, y );

      for ( k = 0; k < order; k++ )
      {
        w3 = 0.5 * ( f - e ) * weight[k];
        z = 0.5 * ( f - e ) * xtab[k] + 0.5 * ( f + e );
        quad = quad + w1 * w2 * w3 * func ( x, y, z );
      }
    }
  }

  volume = b - a;
  result = quad * volume;

  delete [] xtab;
  delete [] weight;

  return result;
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

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
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
//    01 July 2004
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
  double value;

  value = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + value )  )
  {
    value = value / 2.0;
  }

  value = 2.0 * value;

  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
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
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = false;
  fact = 1.0;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
//
//  Discussion:
//
//    Computation is based on an algorithm outlined in references 1 and 2.
//    The program uses rational functions that theoretically approximate
//    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
//    approximation for 12 < X is from reference 3, while approximations
//    for X < 12.0 are similar to those in reference 1, but are unpublished.
//
//    The accuracy achieved depends on the arithmetic system, the compiler,
//    intrinsic functions, and proper selection of the machine-dependent
//    constants.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody, Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
//    Mathematics of Computation,
//    Volume 21, Number 98, April 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maely, Charles Mesztenyi,
//    John Rice, Henry Thatcher, Christop Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.  X must be positive.
//
//    Output, double R8_GAMMA_LOG, the logarithm of the Gamma function of X.
//    If X <= 0.0, or if overflow would occur, the program returns the
//    value XINF, the largest representable double precision number.
//
//
//  Explanation of machine-dependent constants
//
//  BETA   - radix for the real number representation.
//
//  MAXEXP - the smallest positive power of BETA that overflows.
//
//  XBIG   - largest argument for which LN(GAMMA(X)) is representable
//           in the machine, i.e., the solution to the equation
//             LN(GAMMA(XBIG)) = BETA**MAXEXP.
//
//  FRTBIG - Rough estimate of the fourth root of XBIG
//
//
//  Approximate values for some important machines are:
//
//                            BETA      MAXEXP         XBIG
//
//  CRAY-1        (S.P.)        2        8191       9.62E+2461
//  Cyber 180/855
//    under NOS   (S.P.)        2        1070       1.72E+319
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)        2         128       4.08E+36
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)        2        1024       2.55D+305
//  IBM 3033      (D.P.)       16          63       4.29D+73
//  VAX D-Format  (D.P.)        2         127       2.05D+36
//  VAX G-Format  (D.P.)        2        1023       1.28D+305
//
//
//                           FRTBIG
//
//  CRAY-1        (S.P.)   3.13E+615
//  Cyber 180/855
//    under NOS   (S.P.)   6.44E+79
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)   1.42E+9
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)   2.25D+76
//  IBM 3033      (D.P.)   2.56D+18
//  VAX D-Format  (D.P.)   1.20D+9
//  VAX G-Format  (D.P.)   1.89D+76
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
  double d1 = - 5.772156649015328605195174E-01;
  double d2 =   4.227843350984671393993777E-01;
  double d4 =   1.791759469228055000094023E+00;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {
    4.945235359296727046734888E+00, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = {
    4.974607845568932035012064E+00, 
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
    2.940378956634553899906876E+010,
    1.702665737765398868392998E+011,
    4.926125793377430887588120E+011, 
    5.606251856223951465078242E+011 };
  double pnt68 = 0.6796875E+00;
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
    1.488613728678813811542398E+010, 
    1.016803586272438228077304E+011, 
    3.417476345507377132798597E+011, 
    4.463158187419713286462081E+011 };
  double res;
  double sqrtpi = 0.9189385332046727417803297E+00;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
//
//  Return immediately if the argument is out of range.
//
  if ( x <= 0.0 || xbig < x )
  {
    return r8_huge ( );
  }

  if ( x <= r8_epsilon ( ) )
  {
    res = -log ( x );
  }
  else if ( x <= 1.5 )
  {
    if ( x < pnt68 )
    {
      corr = - log ( x );
      xm1 = x;
    }
    else
    {
      corr = 0.0;
      xm1 = ( x - 0.5 ) - 0.5;
    }

    if ( x <= 0.5 || pnt68 <= x )
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
      xm2 = ( x - 0.5 ) - 0.5;
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
  else if ( x <= 4.0 )
  {
    xm2 = x - 2.0;
    xden = 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm2 + p2[i];
      xden = xden * xm2 + q2[i];
    }

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
  }
  else if ( x <= 12.0 )
  {
    xm4 = x - 4.0;
    xden = - 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm4 + p4[i];
      xden = xden * xm4 + q4[i];
    }

    res = d4 + xm4 * ( xnum / xden );
  }
  else
  {
    res = 0.0;

    if ( x <= frtbig )
    {

      res = c[6];
      xsq = x * x;

      for ( i = 0; i < 6; i++ )
      {
        res = res / xsq + c[i];
      }

    }

    res = res / x;
    corr = log ( x );
    res = res + sqrtpi - 0.5 * corr;
    res = res + x * ( corr - 1.0 );

  }

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

double r8_hyper_2f1 ( double a, double b, double c, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HYPER_2F1 evaluates the hypergeometric function 2F1(A,B,C,X).
//
//  Discussion:
//
//    A bug was corrected.  A line which read
//      c1 = - ( - 1.0, m ) * gc / ( gam * gbm * rm );
//    was corrected to read
//      c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );
//    JVB, 05 July 2009.
//
//    A minor bug was corrected.  The HW variable, used in several places as
//    the "old" value of a quantity being iteratively improved, was not
//    being initialized.  JVB, 11 February 2008.
//
//    The FORTRAN77 original version of this routine is copyrighted by
//    Shanjie Zhang and Jianming Jin.  However, they give permission to
//    incorporate this routine into a user program provided that the copyright
//    is acknowledged.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Shanjie Zhang, Jianming Jin,
//    Computation of Special Functions,
//    Wiley, 1996,
//    ISBN: 0-471-11963-6,
//    LC: QA351.C45
//
//  Parameters:
//
//    Input, double A, B, C, X, the arguments of the function.
//    C must not be equal to a nonpositive integer.
//    X < 1.
//
//    Output, double R8_HYPER_2F1, the value of the function.
//
{
  double a0;
  double aa;
  double bb;
  double c0;
  double c1;
  double el = 0.5772156649015329;
  double eps;
  double f0;
  double f1;
  double g0;
  double g1;
  double g2;
  double g3;
  double ga;
  double gabc;
  double gam;
  double gb;
  double gbm;
  double gc;
  double gca;
  double gcab;
  double gcb;
  double gm;
  double hf;
  double hw;
  int j;
  int k;
  bool l0;
  bool l1;
  bool l2;
  bool l3;
  bool l4;
  bool l5;
  int m;
  int nm;
  double pa;
  double pb;
  double pi = 3.141592653589793;
  double r;
  double r0;
  double r1;
  double rm;
  double rp;
  double sm;
  double sp;
  double sp0;
  double x1;

  l0 = ( c == ( int ) ( c ) ) && ( c < 0.0 );
  l1 = ( 1.0 - x < 1.0E-15 ) && ( c - a - b <= 0.0 );
  l2 = ( a == ( int ) ( a ) ) && ( a < 0.0 );
  l3 = ( b == ( int ) ( b ) ) && ( b < 0.0 );
  l4 = ( c - a == ( int ) ( c - a ) ) && ( c - a <= 0.0 );
  l5 = ( c - b == ( int ) ( c - b ) ) && ( c - b <= 0.0 );

  if ( l0 )
  {
    cerr << "\n";
    cerr << "R8_HYPER_2F1 - Fatal error!\n";
    cerr << "  The hypergeometric series is divergent.\n";
    cerr << "  C is integral and negative.\n";
    cerr << "  C = " << c << "\n";
    exit ( 1 );
  }

  if ( l1 )
  {
    cerr << "\n";
    cerr << "R8_HYPER_2F1 - Fatal error!\n";
    cerr << "  The hypergeometric series is divergent.\n";
    cerr << "  1 - X < 0, C - A - B <= 0\n";
    cerr << "  A = " << a << "\n";
    cerr << "  B = " << b << "\n";
    cerr << "  C = " << c << "\n";
    cerr << "  X = " << x << "\n";
    exit ( 1 );
  }

  if ( 0.95 < x )
  {
    eps = 1.0E-08;
  }
  else
  {
    eps = 1.0E-15;
  }

  if ( x == 0.0 || a == 0.0 || b == 0.0 )
  {
    hf = 1.0;
    return hf;
  }
  else if ( 1.0 - x == eps && 0.0 < c - a - b )
  {
    gc = r8_gamma ( c );
    gcab = r8_gamma ( c - a - b );
    gca = r8_gamma ( c - a );
    gcb = r8_gamma ( c - b );
    hf = gc * gcab / ( gca * gcb );
    return hf;
  }
  else if ( 1.0 + x <= eps && r8_abs ( c - a + b - 1.0 ) <= eps )
  {
    g0 = sqrt ( pi ) * pow ( 2.0, - a );
    g1 = r8_gamma ( c );
    g2 = r8_gamma ( 1.0 + a / 2.0 - b );
    g3 = r8_gamma ( 0.5 + 0.5 * a );
    hf = g0 * g1 / ( g2 * g3 );
    return hf;
  }
  else if ( l2 || l3 )
  {
    if ( l2 )
    {
      nm = ( int ) ( r8_abs ( a ) );
    }

    if ( l3 )
    {
      nm = ( int ) ( r8_abs ( b ) );
    }

    hf = 1.0;
    r = 1.0;

    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }

    return hf;
  }
  else if ( l4 || l5 )
  {
    if ( l4 )
    {
      nm = ( int ) ( r8_abs ( c - a ) );
    }

    if ( l5 )
    {
      nm = ( int ) ( r8_abs ( c - b ) );
    }

    hf = 1.0;
    r  = 1.0;
    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }
    hf = pow ( 1.0 - x, c - a - b ) * hf;
    return hf;
  }

  aa = a;
  bb = b;
  x1 = x;

  if ( x < 0.0 )
  {
    x = x / ( x - 1.0 );
    if ( a < c && b < a && 0.0 < b )
    {
      a = bb;
      b = aa;
    }
    b = c - b;
  }

  if ( 0.75 <= x )
  {
    gm = 0.0;

    if ( r8_abs ( c - a - b - ( int ) ( c - a - b ) ) < 1.0E-15 )
    {
      m = ( int ) ( c - a - b );
      ga = r8_gamma ( a );
      gb = r8_gamma ( b );
      gc = r8_gamma ( c );
      gam = r8_gamma ( a + m );
      gbm = r8_gamma ( b + m );

      pa = r8_psi ( a );
      pb = r8_psi ( b );

      if ( m != 0 )
      {
        gm = 1.0;
      }

      for ( j = 1; j <= abs ( m ) - 1; j++ )
      {
        gm = gm * j;
      }

      rm = 1.0;
      for ( j = 1; j <= abs ( m ); j++ )
      {
        rm = rm * j;
      }

      f0 = 1.0;
      r0 = 1.0;;
      r1 = 1.0;
      sp0 = 0.0;;
      sp = 0.0;

      if ( 0 <= m )
      {
        c0 = gm * gc / ( gam * gbm );
        c1 = - gc * pow ( x - 1.0, m ) / ( ga * gb * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( a + k - 1.0 ) + 1.0 / ( b + k - 1.0 ) 
          - 1.0 / ( double ) ( k );
        }

        f1 = pa + pb + sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + ( 1.0 - a ) 
              / ( ( j + k ) * ( a + j + k - 1.0 ) ) 
              + 1.0 / ( b + j + k - 1.0 );
          }

          rp = pa + pb + 2.0 * el + sp + sm + log ( 1.0 - x );

          r1 = r1 * ( a + m + k - 1.0 ) * ( b + m + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }
          hw = f1;
        }
        hf = f0 * c0 + f1 * c1;
      }
      else if ( m < 0 )
      {
        m = - m;
        c0 = gm * gc / ( ga * gb * pow ( 1.0 - x, m ) );
        c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a - m + k - 1.0 ) * ( b - m + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( double ) ( k );
        }

        f1 = pa + pb - sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) 
            / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + 1.0 / ( double ) ( j + k );
          }

          rp = pa + pb + 2.0 * el + sp - sm + log ( 1.0 - x );

          r1 = r1 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }

          hw = f1;
        }

        hf = f0 * c0 + f1 * c1;
      }
    }
    else
    {
      ga = r8_gamma ( a );
      gb = r8_gamma ( b );
      gc = r8_gamma ( c );
      gca = r8_gamma ( c - a );
      gcb = r8_gamma ( c - b );
      gcab = r8_gamma ( c - a - b );
      gabc = r8_gamma ( a + b - c );
      c0 = gc * gcab / ( gca * gcb );
      c1 = gc * gabc / ( ga * gb ) * pow ( 1.0 - x, c - a - b );
      hf = 0.0;
      hw = hf;
      r0 = c0;
      r1 = c1;

      for ( k = 1; k <= 250; k++ )
      {
        r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
          / ( k * ( a + b - c + k ) ) * ( 1.0 - x );

        r1 = r1 * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
          / ( k * ( c - a - b + k ) ) * ( 1.0 - x );

        hf = hf + r0 + r1;

        if ( r8_abs ( hf - hw ) < r8_abs ( hf ) * eps )
        {
          break;
        }
        hw = hf;
      }
      hf = hf + c0 + c1;
    }
  }
  else
  {
    a0 = 1.0;

    if ( a < c && c < 2.0 * a && b < c && c < 2.0 * b )
    {
      a0 = pow ( 1.0 - x, c - a - b );
      a = c - a;
      b = c - b;
    }

    hf = 1.0;
    hw = hf;
    r = 1.0;

    for ( k = 1; k <= 250; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;

      hf = hf + r;

      if ( r8_abs ( hf - hw ) <= r8_abs ( hf ) * eps )
      {
        break;
      }

      hw = hf;
    }
    hf = a0 * hf;
  }

  if ( x1 < 0.0 )
  {
    x = x1;
    c0 = 1.0 / pow ( 1.0 - x, aa );
    hf = c0 * hf;
  }

  a = aa;
  b = bb;

  if ( 120 < k )
  {
    cerr << "\n";
    cerr << "R8_HYPER_2F1 - Warning!\n";
    cerr << "  A large number of iterations were needed.\n";
    cerr << "  The accuracy of the results should be checked.\n";
  }

  return hf;
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

double r8_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8 value.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

double r8_psi ( double xx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI evaluates the function Psi(X).
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
//    09 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
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
//    Input, double XX, the argument of the function.
//
//    Output, double R8_PSI, the value of the function.
//
{
  double aug;
  double den;
  int i;
  int n;
  int nq;
  double one = 1.0;
  double p1[9] = { 
   4.5104681245762934160E-03, 
   5.4932855833000385356, 
   3.7646693175929276856E+02, 
   7.9525490849151998065E+03, 
   7.1451595818951933210E+04, 
   3.0655976301987365674E+05, 
   6.3606997788964458797E+05, 
   5.8041312783537569993E+05, 
   1.6585695029761022321E+05 };
  double p2[7] = { 
  -2.7103228277757834192, 
  -1.5166271776896121383E+01, 
  -1.9784554148719218667E+01, 
  -8.8100958828312219821, 
  -1.4479614616899842986, 
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
  double upper;
  double value;
  double w;
  double x;
  double x01 = 187.0;
  double x01d = 128.0;
  double x02 = 6.9464496836234126266E-04;
  double xinf = 1.70E+38;
  double xlarge = 2.04E+15;
  double xmax1 = 3.60E+16;
  double xmin1 = 5.89E-39;
  double xsmall = 2.05E-09;
  double z;
  double zero = 0.0;

  x = xx;
  w = r8_abs ( x );
  aug = zero;
//
//  Check for valid arguments, then branch to appropriate algorithm.
//
  if ( xmax1 <= - x || w < xmin1 )
  {
    if ( zero < x )
    {
      value = - xinf;
    }
    else
    {
      value = xinf;
    }
    return value;
  }

  if ( x < 0.5 )
  {
//
//  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
//  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
//
    if ( w <= xsmall )
    {
      aug = - one / x;
    }
//
//  Argument reduction for cotangent.
//
    else
    {
      if ( x < zero )
      {
        sgn = piov4;
      }
      else
      {
        sgn = - piov4;
      }

      w = w - ( double ) ( ( int ) ( w ) );
      nq = ( int ) ( w * 4.0 );
      w = 4.0 * ( w - ( double ) ( nq ) * 0.25 );
//
//  W is now related to the fractional part of 4.0 * X.
//  Adjust argument to correspond to values in the first
//  quadrant and determine the sign.
//
      n = nq / 2;

      if ( n + n != nq )
      {
        w = one - w;
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
        if ( z == zero )
        {
          if ( zero < x )
          {
            value = -xinf;
          }
          else
          {
            value = xinf;
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
    x = one - x;
  }
//
//  0.5 <= X <= 3.0.
//
  if ( x <= 3.0 )
  {
    den = x;
    upper = p1[0] * x;
    for ( i = 1; i <= 7; i++ )
    {
      den = ( den + q1[i-1] ) * x;
      upper = ( upper + p1[i]) * x;
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
    w = one / ( x * x );
    den = w;
    upper = p2[0] * w;
    for ( i = 1; i <= 5; i++ )
    {
      den = ( den + q2[i-1] ) * w;
      upper = ( upper + p2[i] ) * w;
    }
    aug = ( upper + p2[6] ) / ( den + q2[5] ) - 0.5 / x + aug;
  }

  value = aug + log ( x );

  return value;
}
//****************************************************************************80

void r8_swap ( double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP switches two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = *x;
  *x = *y;
  *y = z;
 
  return;
}
//****************************************************************************80

void r8_swap3 ( double *x, double *y, double *z )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP3 swaps three R8's.
//
//  Example:
//
//    Input:
//
//      X = 1, Y = 2, Z = 3
//
//    Output:
//
//      X = 2, Y = 3, Z = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double *X, *Y, *Z, three values to be swapped.
//
{
  double w;

   w = *x;
  *x = *y;
  *y = *z;
  *z =  w;

  return;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8ge_det ( int n, double a_lu[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_DET computes the determinant of a matrix factored by R8GE_FA or R8GE_TRF.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA or R8GE_TRF.
//
//    Input, int PIVOT[N], as computed by R8GE_FA or R8GE_TRF.
//
//    Output, double R8GE_DET, the determinant of the matrix.
//
{
  double det;
  int i;

  det = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    det = det * a_lu[i-1+(i-1)*n];
    if ( pivot[i-1] != i )
    {
      det = -det;
    }
  }

  return det;
}
//****************************************************************************80

int r8ge_fa ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FA performs a LINPACK-style PLU factorization of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
//
//    The two dimensional array is stored by columns in a one dimensional
//    array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N], the matrix to be factored.
//    On output, A contains an upper triangular matrix and the multipliers
//    which were used to obtain it.  The factorization can be written
//    A = L * U, where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Output, int PIVOT[N], a vector of pivot indices.
//
//    Output, int R8GE_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int i;
  int j;
  int k;
  int l;
  double t;
//
  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L, the index of the pivot row.
//
    l = k;

    for ( i = k+1; i <= n; i++ )
    {
      if ( r8_abs ( a[l-1+(k-1)*n] ) < r8_abs ( a[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
//
//  If the pivot index is zero, the algorithm has failed.
//
    if ( a[l-1+(k-1)*n] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8GE_FA - Fatal error!\n";
      cerr << "  Zero pivot on step " << k << "\n";
      return k;
    }
//
//  Interchange rows L and K if necessary.
//
    if ( l != k )
    {
      t              = a[l-1+(k-1)*n];
      a[l-1+(k-1)*n] = a[k-1+(k-1)*n];
      a[k-1+(k-1)*n] = t;
    }
//
//  Normalize the values that lie below the pivot entry A(K,K).
//
    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        t              = a[l-1+(j-1)*n];
        a[l-1+(j-1)*n] = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = t;
      }

      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }

    }

  }

  pivot[n-1] = n;

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    cerr << "\n";
    cerr << "R8GE_FA - Fatal error!\n";
    cerr << "  Zero pivot on step " << n << "\n";
    return n;
  }

  return 0;
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

double r8vec_even_select ( int n, double xlo, double xhi, int ival )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
//
//  Discussion:
//
//    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / ( N - 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, double XLO, XHI, the low and high values.
//
//    Input, int IVAL, the index of the desired point.
//    IVAL is normally between 1 and N, but may be any integer value.
//
//    Output, double R8VEC_EVEN_SELECT, the IVAL-th of N evenly spaced values
//    between XLO and XHI.
//    Unless N = 1, X(1) = XLO and X(N) = XHI.
//    If N = 1, then X(1) = 0.5*(XLO+XHI).
//
{
  double xval;

  if ( n == 1 )
  {
    xval = 0.5 * ( xlo + xhi );
  }
  else
  {
    xval = ( ( double ) ( n - ival     ) * xlo 
           + ( double ) (     ival - 1 ) * xhi ) 
           / ( double ) ( n        - 1 );
  }

  return xval;
}
//****************************************************************************80

bool r8vec_mirror_next ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
//
//  Discussion:
//
//    In normal use, the user would set every element of A to be positive.
//    The routine will take the input value of A, and output a copy in
//    which the signs of one or more entries have been changed.  Repeatedly
//    calling the routine with the output from the previous call will generate
//    every distinct "variation" of A; that is, all possible sign variations.
//
//    When the output variable DONE is TRUE (or equal to 1), then the
//    output value of A_NEW is the last in the series.
//
//    Note that A may have some zero values.  The routine will essentially
//    ignore such entries; more exactly, it will not stupidly assume that -0
//    is a proper "variation" of 0.
//
//    Also, it is possible to call this routine with the signs of A set
//    in any way you like.  The routine will operate properly, but it
//    will nonethess terminate when it reaches the value of A in which
//    every nonzero entry has negative sign.
//
//
//    More efficient algorithms using the Gray code seem to require internal
//    memory in the routine, which is not one of MATLAB's strong points,
//    or the passing back and forth of a "memory array", or the use of
//    global variables, or unnatural demands on the user.  This form of
//    the routine is about as clean as I can make it.
//
//  Example:
//
//      Input         Output
//    ---------    --------------
//    A            A         DONE
//    ---------    --------  ----
//     1  2  3     -1  2  3  false
//    -1  2  3      1 -2  3  false
//     1 -2  3     -1 -2  3  false
//    -1 -2  3      1  2 -3  false
//     1  2 -3     -1  2 -3  false
//    -1  2 -3      1 -2 -3  false
//     1 -2 -3     -1 -2 -3  false
//    -1 -2 -3      1  2  3  true
//
//     1  0  3     -1  0  3  false
//    -1  0  3      1  0 -3  false
//     1  0 -3     -1  0 -3  false
//    -1  0 -3      1  0  3  true
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, double A[N], a vector of real numbers.  On 
//    output, some signs have been changed.
//
//    Output, bool R8VEC_MIRROR_NEXT, is TRUE if the input vector A was 
//    the last element
//    in the series (every entry was nonpositive); the output vector is reset
//    so that all entries are nonnegative, but presumably the ride is over.
//
{
  bool done;
  int i;
  int positive;
//
//  Seek the first strictly positive entry of A.
//
  positive = -1;
  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      positive = i;
      break;
    }
  }
//
//  If there is no strictly positive entry of A, there is no successor.
//
  if ( positive == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = - a[i];
    }
    done = true;
    return done;
  }
//
//  Otherwise, negate A up to the positive entry.
//
  for ( i = 0; i <= positive; i++ )
  {
    a[i] = - a[i];
  }
  done = false;

  return done;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    Input, char *TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << setw(6)  << i + 1 << "  " 
         << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO zeroes an R8VEC.
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
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
//****************************************************************************80

double rectangle_3d ( double func ( double x, double y, double z ),
  double a[3], double b[3] )

//****************************************************************************80
//
//  Purpose:
//
//    RECTANGLE_3D approximates an integral inside a rectangular block in 3D.
//
//  Integration region:
//
//      A(1) <= X <= B(1),
//    and
//      A(2) <= Y <= B(2),
//    and
//      A(3) <= Z <= B(3).
//
//  Discussion:
//
//    An 8 point third degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function.
//
//    Input, double A[3], B[3], the lower and upper limits
//    for X, Y and Z.
//
//    Output, double RECTANGLE_3D, the approximate integral of the function.
//
{
  int i;
  int j;
  int k;
  double quad;
  double result;
  double sqr3;
  double volume;
  double w;
  double x;
  double y;
  double z;

  sqr3 = 1.0 / sqrt ( 3.0 );
  w = 1.0 / 8.0;

  quad = 0.0;

  for ( i = 1; i <= 2; i++ )
  {
    x = sqr3 * i4_power ( -1, i );
    x = 0.5 * ( ( 1.0 - x ) * b[0] + ( 1.0 + x ) * a[0] );

    for ( j = 1; j <= 2; j++ )
    {
      y = sqr3 * i4_power ( -1, j );
      y = 0.5 * ( ( 1.0 - y ) * b[1] + ( 1.0 + y ) * a[1] );

      for ( k = 1; k <= 2; k++ )
      {
        z = sqr3 * i4_power ( -1, k );
        z = 0.5 * ( ( 1.0 - z ) * b[2] + ( 1.0 + z ) * a[2] );

        quad = quad + w * func ( x, y, z );
      }
    }
  }
  volume = ( b[0] - a[0] ) * ( b[1] - a[1] ) * ( b[2] - a[2] );
  result = volume * quad;

  return result;
}
//****************************************************************************80

double rectangle_sub_2d ( double func ( double x, double y ), double xval[2], 
  double yval[2], int nsub[2], int order, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    RECTANGLE_SUB_2D carries out a composite quadrature over a rectangle in 2D.
//
//  Integration region:
//
//      XVAL(1) <= X <= XVAL(2),
//    and
//      YVAL(1) <= Y <= YVAL(2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the function 
//    to be integrated.
//
//    Input, double XVAL[2], the left and right X coordinates.
//
//    Input, double YVAL[2], the lower and upper Y coordinates.
//
//    Input, int NSUB[2], the number of subintervals to use in the X 
//    and Y directions.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double RECTANGLE_SUB_2D, the approximate integral of the function.
//
{
  double a[2];
  double b[2];
  int i;
  int j;
  int k;
  double quad_sub;
  double result;
  double result_sub;
  double volume;
  double volume_sub;
  double x;
  double xhi;
  double xlo;
  double y;
  double yhi;
  double ylo;

  a[0] = xval[0];
  a[1] = yval[0];
  b[0] = xval[1];
  b[1] = yval[1];

  for ( i = 0; i < 2; i++ )
  {
    if ( a[i] == b[i] )
    {
      result = 0.0;
      return result;
    }
  }

  for ( i = 0; i < 2; i++ )
  {
    if ( nsub[i] < 1 )
    {
      cerr << "\n";
      cerr << "RECTANGLE_SUB_2D - Fatal error!\n";
      cerr << "  Nonpositive value of NSUB[" << i 
           << "] = " << nsub[i] << "\n";
      exit ( 1 );
    }
  }
//
//  Break up the X interval into NSUB(1) subintervals.
//
  volume = 0.0;
  result = 0.0;

  for ( i = 1; i <= nsub[0]; i++ )
  {
     xlo = r8vec_even_select ( nsub[0]+1, a[0], b[0], i   );
     xhi = r8vec_even_select ( nsub[0]+1, a[0], b[0], i + 1 );
//
//  Break up the Y interval into NSUB(2) subintervals.
//
    for ( j = 1; j <= nsub[1]; j++ )
    {
      ylo = r8vec_even_select ( nsub[1]+1, a[1], b[1], j   );
      yhi = r8vec_even_select ( nsub[1]+1, a[1], b[1], j+1 );

      quad_sub = 0.0;
      for ( k = 0; k < order; k++ )
      {
        x = xlo + 0.5 * ( xtab[k] + 1.0 ) * ( xhi - xlo );
        y = ylo + 0.5 * ( ytab[k] + 1.0 ) * ( yhi - ylo );

        quad_sub = quad_sub + weight[k] * func ( x, y ) / 4.0;
      }

      volume_sub = ( xhi - xlo ) * ( yhi - ylo );
      result_sub = quad_sub * volume_sub;

      volume = volume + volume_sub;
      result = result + result_sub;
    }
  }

  return result;
}
//****************************************************************************80

void rule_adjust ( double a, double b, double c, double d, int order, 
  double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
//
//  Discussion:
//
//    Most quadrature rules are defined on a special interval, like
//    [-1,1] or [0,1].  To integrate over an interval, the abscissas
//    and weights must be adjusted.  This can be done on the fly,
//    or by calling this routine.
//
//    If the weight function W(X) is not 1, then the W vector will
//    require further adjustment by the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the definition interval.
//
//    Input, double C, D, the endpoints of the integration interval.
//
//    Input, int ORDER, the number of abscissas and weights.
//
//    Input/output, double X[ORDER], W[ORDER], the abscissas
//    and weights.
//
{
  int i;

  for ( i = 0; i < order; i++ )
  {
    x[i] = ( ( b - x[i]     ) * c   
           + (     x[i] - a ) * d ) 
           / ( b              - a );
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = ( ( d - c ) / ( b - a ) ) * w[i];
  }

  return;
}
//****************************************************************************80

int s_len_trim ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n ) 
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

double simplex_nd ( double func ( int n, double x[] ), int n, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_ND approximates an integral inside a simplex in ND.
//
//  Discussion:
//
//    An N+1 point second degree formula is used.
//
//    The integration region is the simplex bounded by the origin and a 
//    convex combination of N points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function which evaluates F(X).
//
//    Input, int N, the dimension of the space.
//
//    Input/output, double V[N*(N+1)].  On input, each of the
//    N+1 columns of V contains the N coordinates of one of the
//    "corners" of the simplex in entries 1 through N, with
//    the last column being left free.
//    On output, V has been overwritten in the process of
//    computing the volume of the simplex.
//
//    Output, double SIMPLEX_ND, the approximate integral of the function.
//
{
  double c;
  int i;
  int j;
  double quad;
  double result;
  double s;
  double volume;
  double w;
  double *x;

  x = new double[n];

  c = 1.0 / sqrt ( ( double ) ( n + 2 ) );
  w = 1.0 / ( double ) ( n + 1 );

  for ( j = 0; j < n; j++ )
  {
    s = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      s = s + v[i+j*n];
    }
    x[j] = w * ( 1.0 - c ) * s;
  }

  quad = 0.0;

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = x[i] + c * v[i+j*n];
    }
    quad = quad + w * func ( n, x );
    for ( i = 0; i < n; i++ )
    {
      x[i] = x[i] - c * v[i+j*n];
    }
  }

  volume = simplex_volume_nd ( n, v );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double simplex_unit_01_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_01_ND approximates an integral inside the unit simplex in ND.
//
//  Integration region:
//
//      0 <= X(1:N),
//    and
//      sum ( X(1:N) ) <= 1.
//
//  Discussion:
//
//    A 1 point formula of degree 1 is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Axel Grundmann, Michael Moeller,
//    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
//    SIAM Journal on Numerical Analysis,
//    Volume 15, Number 2, April 1978, pages 282-290.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.  
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  double coef = 1.0;
  int i;
  double quad;
  double result;
  double volume;
  double *x;

  x = new double[n];

  quad = 0.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0 / ( double ) ( n );
  }
  quad = quad + coef * func ( n, x );

  volume = simplex_unit_volume_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double simplex_unit_03_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_03_ND approximates an integral inside the unit simplex in ND.
//
//  Integration region:
//
//      0 <= X(1:N),
//    and
//      sum ( X(1:N) ) <= 1.
//
//  Discussion:
//
//    An N+2 point formula of degree 3 is used.  This is Stroud TN:3-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Axel Grundmann, Michael Moeller,
//    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
//    SIAM Journal on Numerical Analysis,
//    Volume 15, Number 2, April 1978, pages 282-290.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function which is to be integrated.
//
//    Input, int N, the dimension of the space.  
//
//    Output, double SIMPLEX_UNIT_03_ND, the approximate integral of the function.
//
{
  double a;
  double b;
  double coef;
  int i;
  double quad;
  double result;
  double volume;
  double *x;

  x = new double[n];

  quad = 0.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0 / ( double ) ( n + 1 );
  }
  coef = -0.25 * ( double ) ( ( n + 1 ) * ( n + 1 ) )  / ( double ) ( n + 2 );
  quad = quad + coef * func ( n, x );

  a = 1.0 / ( double ) ( n + 3 );
  b = 3.0 / ( double ) ( n + 3 );

  for ( i = 0; i < n; i++ )
  {
    x[i] = a;
  }
  coef = 0.25 * ( double ) ( ( n + 3 ) * ( n + 3 ) ) 
    / ( double ) ( ( n + 1 ) * ( n + 2 ) );
  quad = quad + coef * func ( n, x );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b;
    quad = quad + coef * func ( n, x );
    x[i] = a;
  }

  volume = simplex_unit_volume_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double simplex_unit_05_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_05_ND approximates an integral inside the unit simplex in ND.
//
//  Integration region:
//
//      0 <= X(1:N),
//    and
//      sum ( X(1:N) ) <= 1.
//
//  Discussion:
//
//    An N^2 + 3 N + 3 point formula of degree 5 is used.  This is
//    Stroud formula TN:5-1.
//
//    (For N = 2, the number of points is actually only 7, and
//     for N = 3, the number of points is actually only 15.)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    A Fifth Degree Integration Formula for the N-Simplex,
//    SIAM Journal on Numerical Analysis,
//    Volume 6, Number 1, March 1969.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function is to be integrated.
//
//    Input, int N, the dimension of the space.  For this routine,
//    it must be the case that 2 <= N <= 16.
//
//    Output, double SIMPLEX_UNIT_05_ND, the approximate integral of the function.
//
{
  double coef1[16] = {
    0.0,
    0.225, 
    0.118518518519, 
    0.0631521898883, 
    0.235714285714, 
    0.791575476992, 
    1.85798728021, 
    3.53666958042, 
    5.90844340844, 
    9.03765432098, 
    12.9758241758, 
    17.7645108738, 
    23.4375030259, 
    30.0224941950, 
    37.5423613501, 
    46.0161454949 };
  double coef21[16] = { 
    0.0,
    0.12593918054483, 
    0.0719370837790, 
    0.0470456145702, 
    0.0333009774677, 
    0.0248633014592, 
    0.0192679696358, 
    0.0153322153879, 
    0.0124316229901, 
    0.0102112988361, 
    0.00845730697460, 
    0.00703433430999, 
    0.00585330520067, 
    0.00485356735291, 
    0.00399261092720, 
    0.00323988713017 };
  double coef22[16] = {
    0.0,
    0.13239415278851, 
    0.0690682072263, 
    0.0371530185868, 
   -0.0719253160920, 
   -0.264323879461, 
   -0.537926779961, 
   -0.886895605701, 
   -1.30409181465, 
   -1.78227048964, 
   -2.31462336314, 
   -2.89499045158, 
   -3.51790849765, 
   -4.17858310668, 
   -4.87282884913, 
   -5.59699944261 };
  double coef31[16] = {
    0.0,
    0.0, 
    0.0529100529100, 
    0.0261368740713, 
    0.0499020181331, 
    0.0782233395867, 
    0.109041040862, 
    0.140874828568, 
    0.172735353396, 
    0.203992490408, 
    0.234263814181, 
    0.263332763315, 
    0.291091849264, 
    0.317504208212, 
    0.342577872069, 
    0.366348654344 };
  double coef32[16] = {
    0.0,
    0.0, 
    0.0, 
    0.0254485903613, 
    0.0165000982690, 
    0.0115218303668, 
    0.00850478779483, 
    0.00655297510968, 
    0.00522372456259, 
    0.00428017828134, 
    0.00358722367033, 
    0.00306362964360, 
    0.00265836687133, 
    0.00233816221525, 
    0.00208061510846, 
    0.00187022027571 };
  int i;
  int j;
  double quad;
  double r1;
  double r2;
  double result;
  double s1;
  double s2;
  double u1;
  double u2;
  double v1;
  double v2;
  double volume;
  double *x;

  result = 0.0;

  if ( n < 2 || 16 < n )
  {
    cerr << "\n";
    cerr << "SIMPLEX_UNIT_05_ND - Fatal error!\n";
    cerr << "  Input spatial dimension N out of range.\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }

  x = new double[n];

  quad = 0.0;
//
//  S1
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0 / ( double ) ( n + 1 );
  }
  quad = quad + coef1[n-1] * func ( n, x );
//
//  S21
//
  r1 = ( ( double ) ( n + 4 ) - sqrt ( 15.0 ) ) 
    / ( double ) ( n * n + 8 * n + 1 );
  s1 = 1.0 - ( double ) ( n ) * r1;

  for ( i = 0; i < n; i++ )
  {
    x[i] = r1;
  }

  for ( i = 0; i < n + 1; i++ )
  {
    quad = quad + coef21[n-1] * func ( n, x );

    if ( 0 < i ) 
    {
      x[i-1] = r1;
    }
    if ( i < n )
    {
      x[i] = s1;
    }
  }
//
//  S22
//
  r2 = ( ( double ) ( n + 4 ) + sqrt ( 15.0 ) ) 
    / ( double ) ( n * n + 8 * n + 1 );
  s2 = 1.0 - ( double ) ( n ) * r2;

  for ( i = 0; i < n; i++ )
  {
    x[i] = r2;
  }
  for ( i = 0; i < n + 1; i++ )
  {
    quad = quad + coef22[n-1] * func ( n, x );

    if ( 0 < i )
    {
      x[i-1] = r2;
    }
    if ( i < n )
    {
      x[i] = s2;
    }
  }
//
//  S31
//
  u1 = ( ( double ) ( n + 7 ) + 2.0 * sqrt ( 15.0 ) ) 
    / ( double ) ( n * n + 14 * n - 11 );
  v1 = ( ( double ) ( 4 * n - 2 ) 
    - ( double ) ( n - 1 ) * sqrt ( 15.0 ) ) 
    / ( double ) ( n * n + 14 * n - 11 );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = u1;
    }
    x[i] = v1;

    for ( j = i; j < n; j++ )
    {
      if ( i < j - 1 )
      {
        x[j-1] = u1;
      }

      x[j] = v1;

      quad = quad + coef31[n-1] * func ( n, x );
    }
  }
//
//  S32
//
  u2 = ( ( double ) ( n + 7 ) - 2.0 * sqrt ( 15.0 ) ) 
    / ( double ) ( n * n + 14 * n - 11 );
  v2 = ( ( double ) ( 4 * n - 2 ) 
    + ( double ) ( n - 1 ) * sqrt ( 15.0 ) ) 
    / ( double ) ( n * n + 14 * n - 11 );

  for ( i = 0; i < n; i++ )
  {

    for ( j = 0; j < n; j++ )
    {
      x[j] = u2;
    }
    x[i] = v2;

    for ( j = i; j < n; j++ )
    {
      if ( i < j - 1 )
      {
        x[j-1] = u2;
      }

      x[j] = v2;

      quad = quad + coef32[n-1] * func ( n, x );
    }
  }

  volume = simplex_unit_volume_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double simplex_unit_05_2_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_05_2_ND approximates an integral inside the unit simplex in ND.
//
//  Integration region:
//
//      0 <= X(1:N),
//    and
//      sum ( X(1:N) ) <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Axel Grundmann, Michael Moeller,
//    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
//    SIAM Journal on Numerical Analysis,
//    Volume 15, Number 2, April 1978, pages 282-290.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Output, double SIMPLEX_UNIT_05_2_ND, the approximate integral of the function.
//
{
  double a;
  double b;
  double coef;
  int i;
  int j;
  double quad;
  double result;
  double volume;
  double *x;

  x = new double[n];

  quad = 0.0;
//
//  Group 1
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0 / ( double ) ( n + 1 );
  }
  coef = ( double ) ( i4_power ( n + 1, 4 ) ) 
    / ( double ) ( 32 * ( n + 2 ) * ( n + 3 ) );
  quad = quad + coef * func ( n, x );
//
//  Group 2
//
  a = 1.0 / ( double ) ( n + 3 );
  b = 3.0 / ( double ) ( n + 3 );

  for ( i = 0; i < n; i++ )
  {
    x[i] = a;
  }
  coef = -( double ) ( i4_power ( n + 3, 4 ) ) 
    / ( double ) ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 4 ) );
  quad = quad + coef * func ( n, x );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b;
    quad = quad + coef * func ( n, x );
    x[i] = a;
  }
//
//  Group 3
//
  a = 1.0 / ( double ) ( n + 5 );
  b = 5.0 / ( double ) ( n + 5 );

  for ( i = 0; i < n; i++ )
  {
    x[i] = a;
  }
  coef = ( double ) ( i4_power ( n + 5, 4 ) ) 
    / ( double ) ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) );
  quad = quad + coef * func ( n, x );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b;
    quad = quad + coef * func ( n, x );
    x[i] = a;
  }
//
//  Group 4
//
  a = 1.0 / ( double ) ( n + 5 );
  b = 3.0 / ( double ) ( n + 5 );

  coef = ( double ) ( i4_power ( n + 5, 4 ) ) 
    / ( double ) ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = a;
    }
    x[i] = b;
    quad = quad + coef * func ( n, x );

    for ( j = i + 1; j < n; j++ )
    {
      x[j] = b;
      quad = quad + coef * func ( n, x );
      x[j] = a;
    }
  }

  volume = simplex_unit_volume_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double simplex_unit_volume_nd ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_UNIT_VOLUME_ND returns the volume of the unit simplex in ND.
//
//  Integration region:
//
//      0 <= X(1:N),
//    and
//      sum ( X(1:N) ) <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Output, double SIMPLEX_UNIT_VOLUME_ND, the volume of the
//    unit simplex.
//
{
  double value;

  value = 1.0 / ( double ) ( i4_factorial ( n ) );

  return value;
}
//****************************************************************************80

double simplex_volume_nd ( int n, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_VOLUME_ND returns the volume of a simplex in ND.
//
//  Integration region:
//
//    The simplex bounded by the origin and a convex combination of N points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double V[N*(N+1)], the coordinates of the vertices.
//
//    Output, double SIMPLEX_VOLUME_ND, the volume of 
//    the unit simplex.
//
{
  double det;
  int i;
  int info;
  int j;
  int *pivot;
  double volume;
  double *w;

  w = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      w[i+j*n] = v[i+(j+1)*n] - v[i+0*n];
    }
  }

  pivot = new int[n];

  info = r8ge_fa ( n, w, pivot );

  if ( info != 0 )
  {
    volume = 0.0;
  }
  else
  {
    det = r8ge_det ( n, w, pivot );
//
//  Multiply by the volume of the unit simplex, which serves as a
//  conversion factor between a parallelipiped and the simplex.
//
    volume = r8_abs ( det ) * simplex_unit_volume_nd ( n );
  }

  delete [] pivot;
  delete [] w;

  return volume;
}
//****************************************************************************80

double sin_power_int ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_POWER_INT evaluates the sine power integral.
//
//  Discussion:
//
//    The function is defined by
//
//      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
//
//    The algorithm uses the following fact:
//
//      Integral sin^n ( t ) = (1/n) * (
//        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, double A, B, the limits of integration.
//
//    Input, integer N, the power of the sine function.
//
//    Output, double SIN_POWER_INT, the value of the integral.
//
{
  double ca;
  double cb;
  int m;
  int mlo;
  double sa;
  double sb;
  double value;

  value = 0.0;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "SIN_POWER_INT - Fatal error!\n";
    cerr << "  Power N < 0.\n";
    exit ( 1 );
  }

  sa = sin ( a );
  sb = sin ( b );
  ca = cos ( a );
  cb = cos ( b );

  if ( ( n % 2 ) == 0 )
  {
    value = b - a;
    mlo = 2;
  }
  else
  {
    value = ca - cb;
    mlo = 3;
  }

  for ( m = mlo; m <= n; m = m + 2 )
  {
    value = ( ( double ) ( m - 1 ) * value 
      + pow ( sa, (m-1) ) * ca - pow ( sb, (m-1) ) * cb ) 
      / ( double ) ( m );
  }

  return value;
}
//****************************************************************************80

double sphere_05_nd ( double func ( int n, double x[] ), int n, double center[], 
  double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_05_ND approximates an integral on the surface of a sphere in ND.
//
//  Integration region:
//
//    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2
//
//  Discussion:
//
//    A 2*N+2^N points 5-th degree formula is used, Stroud number UN:5-2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Input, double CENTER[N], the center of the sphere.
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_05_ND, the approximate integral of the function.
//
{
  int i;
  int iadd;
  int ihi;
  int *ix;
  bool more;
  int ncard;
  double quad;
  double result;
  double volume;
  double w1;
  double w2;
  double *x;
  double x1;
  double x2;

  x1 = 1.0;
  x2 = 1.0 / sqrt ( ( double ) ( n ) );

  w1 = 1.0 / ( double ) ( n * ( n + 2 ) );
  w2 = ( double ) ( n ) / ( double ) ( ( n + 2 ) * i4_power ( 2, n ) );

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i];
  }
  quad = 0.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i] + r * x1;
    quad = quad + w1 * func ( n, x );
    x[i] = center[i] - r * x1;
    quad = quad + w1 * func ( n, x );
    x[i] = center[i];
  }

  more = false;
  ihi = i4_power ( 2, n );

  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i] - r * x2;
  }

  ix = new int[n];

  for ( i = 0; i < ihi; i++ )
  {
    subset_gray_next ( n, ix, &more, &ncard, &iadd );

    if ( iadd != -1 )
    {
      x[iadd-1] = center[iadd-1] - ( x[iadd-1] - center[iadd-1] );
    }
    quad = quad + w2 * func ( n, x );
  }

  volume = sphere_area_nd ( n, r );
  result = quad * volume;

  delete [] ix;
  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_07_1_nd ( double func ( int n, double x[] ), int n, 
  double center[], double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_07_1_ND approximates an integral on the surface of a sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N) - CENTER(1:N) )^2 = R * R.
//
//  Discussion:
//
//    A 2^N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Input, double CENTER[N], the center of the sphere.
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_07_1_ND, the approximate integral of the function.
//
{
  int i;
  int iadd;
  int *ix;
  int j;
  int jhi;
  bool more;
  int ncard;
  double quad;
  double result;
  double volume;
  double w1;
  double w2;
  double w3;
  double *x;
  double x1;
  double x2;
  double x3;

  ix = new int[n];
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i];
  }
  w1 = ( double ) ( 8 - n ) 
     / ( double ) ( n * ( n + 2 ) * ( n + 4 ) );

  w2 = ( double ) ( n * n * n ) 
     / ( double ) ( i4_power ( 2, n ) * n * ( n + 2 ) * ( n + 4 ) );

  w3 = 4.0 / ( double ) ( n * ( n + 2 ) * ( n + 4 ) );

  x1 = 1.0;
  x2 = 1.0 / sqrt ( ( double ) ( n ) );
  x3 = 1.0 / sqrt ( 2.0 );

  quad = 0.0;
//
//  First term.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i] + r * x1;
    quad = quad + w1 * func ( n, x );
    x[i] = center[i] - r * x1;
    quad = quad + w1 * func ( n, x );
    x[i] = center[i];
  }
//
//  Second term.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i] - r * x2;
  }
  more = false;
  jhi = i4_power ( 2, n );

  for ( j = 0; j < jhi; j++ )
  {
    subset_gray_next ( n, ix, &more, &ncard, &iadd );

    if ( iadd != -1 )
    {
      x[iadd-1] = center[iadd-1] - ( x[iadd-1] - center[iadd-1] );
    }
    quad = quad + w2 * func ( n, x );
  }
//
//  Third term.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i];
  }
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      x[i] = center[i] + r * x3;
      x[j] = center[j] + r * x3;
      quad = quad + w3 * func ( n, x );
      x[i] = center[i] - r * x3;
      x[j] = center[j] + r * x3;
      quad = quad + w3 * func ( n, x );
      x[i] = center[i] + r * x3;
      x[j] = center[j] - r * x3;
      quad = quad + w3 * func ( n, x );
      x[i] = center[i] - r * x3;
      x[j] = center[j] - r * x3;
      quad = quad + w3 * func ( n, x );
      x[i] = center[i];
      x[j] = center[j];
    }
  }

  volume = sphere_area_nd ( n, r );
  result = quad * volume;

  delete [] ix;
  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_area_3d ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_AREA_3D computes the area of a sphere in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z = R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_AREA_3D, the area of the sphere.
//
{
  double pi = 3.141592653589793;
  double value;

  value = 4.0 * pi * r * r;

  return value;
}
//****************************************************************************80

double sphere_area_nd ( int n, double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_AREA_ND computes the area of a sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = R * R
//
//  Discussion:
//
//    N   Area
//
//    2   2       * PI   * R
//    3   4       * PI   * R^2
//    4   2       * PI^2 * R^3
//    5   (8/3)   * PI^2 * R^4
//    6             PI^3 * R^5
//    7   (16/15) * PI^3 * R^6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_AREA_ND, the area of the sphere.
//
{
  double value;

  value = sphere_unit_area_nd ( n ) * pow ( r, n - 1 );

  return value;
}
//****************************************************************************80

double sphere_cap_area_2d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CAP_AREA_2D computes the surface area of a spherical cap in 2D.
//
//  Discussion:
//
//    Draw any radius of the sphere and note the point P where the radius
//    intersects the sphere.  Consider the point on the radius line which is
//    H units from P.  Draw the circle that lies in the plane perpendicular to
//    the radius, and which intersects the sphere.  The circle divides the sphere
//    into two pieces, and the corresponding disk divides the solid sphere into
//    two pieces.  The spherical cap is the part of the solid sphere that
//    includes the point P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Input, double H, the "height" of the spherical cap. 
//
//    Output, double SPHERE_CAP_AREA_2D, the area of the spherical cap.
//
{
  double area;
  double pi = 3.141592653589793;
  double theta;

  if ( h <= 0.0 )
  {
    area = 0.0;
  }
  else if ( 2.0 * r <= h )
  {
    area = 2.0 * pi * r;
  }
  else
  {
    theta = 2.0 * arc_sine ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );

    area = r * theta;

    if ( r <= h )
    {
      area = 2.0 * pi * r - area;
    }
  }

  return area;
}
//****************************************************************************80

double sphere_cap_area_3d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CAP_AREA_3D computes the surface area of a spherical cap in 3D.
//
//  Discussion:
//
//    Draw any radius of the sphere and note the point P where the radius
//    intersects the sphere.  Consider the point on the radius line which is
//    H units from P.  Draw the circle that lies in the plane perpendicular to
//    the radius, and which intersects the sphere.  The circle divides the sphere
//    into two pieces, and the corresponding disk divides the solid sphere into
//    two pieces.  The spherical cap is the part of the solid sphere that
//    includes the point P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Input, double H, the "height" of the spherical cap. 
//
//    Output, double SPHERE_CAP_AREA_3D, the area of the spherical cap.
//
{
  double area;
  double pi = 3.141592653589793;

  if ( h <= 0.0 )
  {
    area = 0.0;
  }
  else if ( 2.0 * r <= h )
  {
    area = 4.0 * pi * r * r;
  }
  else
  {
    area = 2.0 * pi * r * h;
  }

  return area;
}
//****************************************************************************80

double sphere_cap_area_nd ( int dim_num, double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CAP_AREA_ND computes the area of a spherical cap in ND.
//
//  Discussion:
//
//    The spherical cap is a portion of the surface of the sphere:
//
//      sum ( X(1:N)^2 ) = R * R
//
//    which is no more than H units from the uppermost point on the sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Thomas Ericson, Victor Zinoviev,
//    Codes on Euclidean Spheres,
//    Elsevier, 2001
//    QA166.7 E75
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, double R, the radius of the sphere.
//
//    Input, double H, the "thickness" of the spherical cap,
//    which is normally between 0 and 2 * R.
//
//    Output, double SPHERE_CAP_AREA_ND, the area of the spherical cap.
//
{
  double area;
  double haver_sine;
  int i;
  double theta;
  double ti;
  double tj;
  double tk;

  if ( h <= 0.0 )
  {
    area = 0.0;
    return area;
  }

  if ( 2.0 * r <= h )
  {
    area = sphere_area_nd ( dim_num, r );
    return area;
  }
//
//  For cases where R < H < 2 * R, work with the complementary region.
//
  haver_sine = sqrt ( ( 2.0 * r - h ) * h );

  theta = arc_sine ( haver_sine / r );

  if ( dim_num < 1 )
  {
    area = -1.0;
  }
  else if ( dim_num == 1 )
  {
    area = 0.0;
  }
  else if ( dim_num == 2 )
  {
    area = 2.0 * theta * r;
  }
  else
  {
    ti = theta;

    tj = ti;
    ti = 1.0 - cos ( theta );

    for ( i = 2; i <= dim_num - 2; i++ )
    {
      tk = tj;
      tj = ti;
      ti = ( ( double ) ( i - 1 ) * tk 
        - cos ( theta ) * pow ( sin ( theta ), i - 1 ) ) 
        / ( double ) ( i );
    }
    area = sphere_k ( dim_num-1 ) * ti * pow ( r, dim_num - 1 );
  }
//
//  Adjust for cases where R < H < 2R.
//
  if ( r < h )
  {
    area = sphere_area_nd ( dim_num, r ) - area;
  }
  return area;
}
//****************************************************************************80

double sphere_cap_volume_2d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CAP_VOLUME_2D computes the volume of a spherical cap in 2D.
//
//  Discussion:
//
//    Draw any radius R of the circle and denote as P the point where the
//    radius intersects the circle.  Now consider the point Q which lies
//    on the radius and which is H units from P.  The line which is
//    perpendicular to the radius R and passes through Q divides the
//    circle into two pieces.  The piece including the point P is the
//    spherical (circular) cap of height (or thickness) H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Input, double H, the "height" of the spherical cap.
//
//    Output, double VOLUME, the volume (area) of the spherical cap.
//
{
  double pi = 3.141592653589793;
  double theta;
  double volume;

  if ( h <= 0.0 )
  {
    volume = 0.0;
  }
  else if ( 2.0 * r <= h )
  {
    volume = pi * r * r;
  }
  else
  {
    theta = 2.0 * arc_sine ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
    volume = r * r * ( theta - sin ( theta ) ) / 2.0;
    if ( r < h )
    {
      volume = pi * r * r - volume;
    }
  }
  return volume;
}
//****************************************************************************80

double sphere_cap_volume_3d ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CAP_VOLUME_3D computes the volume of a spherical cap in 3D.
//
//  Discussion:
//
//    Draw any radius of the sphere and note the point P where the radius
//    intersects the sphere.  Consider the point on the radius line which is
//    H units from P.  Draw the circle that lies in the plane perpendicular to
//    the radius, and which intersects the sphere.  The circle divides the sphere
//    into two pieces, and the corresponding disk divides the solid sphere into
//    two pieces.  The part of the solid sphere that includes the point P
//    is the spherical cap of height (or thickness) H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Input, double H, the "height" of the spherical cap.
//
//    Output, double SPHERE_CAP_VOLUME_3D, the volume of the spherical cap.
//
{
  double pi = 3.141592653589793;
  double volume;

  if ( h <= 0.0 )
  {
    volume = 0.0;
  }
  else if ( 2.0 * r <= h )
  {
    volume = ( 4.0 / 3.0 ) * pi * r * r * r;
  }
  else
  {
    volume = ( 1.0 / 3.0 ) * pi * h * h * ( 3.0 * r - h );
  }
  return volume;
}
//****************************************************************************80

double sphere_cap_volume_nd ( int dim_num, double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CAP_VOLUME_ND computes the volume of a spherical cap in ND.
//
//  Discussion:
//
//    The spherical cap is a portion of the surface and interior of the sphere:
//
//      sum ( X(1:N)^2 ) <= R * R
//
//    which is no more than H units from some point P on the sphere.
//
//
//    The algorithm proceeds from the observation that the N-dimensional
//    sphere can be parameterized by a quantity RC that runs along the
//    radius from the center to the point P.  The value of RC at the
//    base of the spherical cap is (R-H) and at P it is R.  We intend to
//    use RC as our integration parameeter.
//
//    The volume of the spherical cap is then the integral, as RC goes
//    from (R-H) to R, of the N-1 dimensional volume of the sphere
//    of radius RS, where RC * RC + RS * RS = R * R.
//
//    The volume of the N-1 dimensional sphere of radius RS is simply 
//    some constants times RS**(N-1).
// 
//    After factoring out the constant terms, and writing RC = R * cos ( T ),
//    and RS = R * sin ( T ), and letting 
//      T_MAX = arc_sine ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) ),
//    the "interesting part" of our integral becomes
//
//      constants * R**N * Integral ( T = 0 to T_MAX ) sin**N ( T ) dT
//
//    The integral of sin**N ( T ) dT can be handled by recursion.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the space.
//
//    Input, double R, the radius of the sphere.
//
//    Input, double H, the "thickness" of the spherical cap,
//    which is normally between 0 and 2 * R.
//
//    Output, double SPHERE_CAP_VOLUME_ND, the volume of the spherical cap.
//
{
  double angle;
  double factor1;
  double factor2;
  double volume;
  double volume2;

  if ( h <= 0.0 )
  {
    volume = 0.0;
    return volume;
  }

  if ( 2.0 * r <= h )
  {
    volume = sphere_volume_nd ( dim_num, r );
    return volume;
  }

  if ( dim_num < 1 )
  {
    volume = -1.0;
  }
  else if ( dim_num == 1 )
  {
    volume = h;
  }
  else
  {
    factor1 = sphere_unit_volume_nd ( dim_num - 1 );

    angle = arc_sine ( sqrt ( ( 2.0 * r - h ) * h / r ) );

    factor2 = sin_power_int ( 0.0, angle, dim_num );

    volume = factor1 * factor2 * pow ( r, dim_num );

    if ( r < h )
    {
      volume2 = sphere_volume_nd ( dim_num, r );
      volume = volume2 - volume;
    }
  }
  return volume;
}
//****************************************************************************80

double sphere_k ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_K computes a factor useful for spherical computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Thomas Ericson, Victor Zinoviev,
//    Codes on Euclidean Spheres,
//    Elsevier, 2001
//    QA166.7 E75
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Output, double SPHERE_K, the factor.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( ( n % 2 ) == 0 )
  {
    value = pow ( 2.0 * pi, n / 2 );
  }
  else
  {
    value = 2.0 * pow ( 2.0 * pi, ( n - 1 ) / 2 );
  }
  value = value / ( double ) ( i4_factorial2 ( n - 2 ) );

  return value;
}
//****************************************************************************80

double sphere_monomial_int_nd ( int n, double r, int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_MONOMIAL_INT_ND integrates a monomial on the surface of a sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = R * R.
//
//  Discussion:
//
//    The sphere may have nonunit radius, but it must be centered at 0.
//
//    The monomial is F(X) = X(1)^E(1) * X(2)^E(2) * ... * X(N)^E(N).
//
//    This routine is useful for testing the accuracy of quadrature
//    rules on the sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double R, the radius of the sphere.
//
//    Input, int E[N], the exponents of X, Y and Z in the monomial.
//    Each exponent must be nonnegative.
//
//    Output, double SPHERE_MONOMIAL_INT_ND, the integral.
//
{
  bool all_zero;
  bool any_odd;
  int e_sum;
  int i;
  double integral;
  double pi = 3.141592653589793;

  integral = 0.0;

  for ( i = 0; i < n; i++ )
  {
    if ( e[i] < 0 )
    {
      cerr << "\n";
      cerr << "SPHERE_MONOMIAL_INT_ND - Fatal error!\n";
      cerr << "  All exponents must be nonnegative.\n";
      exit ( 1 );
    }
  }

  all_zero = true;
  for ( i = 0; i < n; i++ )
  {
    if ( e[i] != 0 )
    {
      all_zero = false;
      break;
    }
  }

  any_odd = false;
  for ( i = 0; i < n; i++ )
  {
    if ( ( e[i] % 2 ) == 1 )
    {
      any_odd = true;
      break;
    }
  }

  e_sum = 0;
  for ( i = 0; i < n; i++ )
  {
    e_sum = e_sum + e[i];
  }

  if ( all_zero )
  {
    integral = 2.0 * sqrt ( pow ( pi, n ) ) 
      / r8_gamma ( 0.5 * ( double ) ( n ) );
  }
  else if ( any_odd )
  {
    integral = 0.0;
  }
  else
  {
    integral = 2.0;

    for ( i = 0; i < n; i++ )
    {
      integral = integral * r8_gamma ( 0.5 * ( double ) ( e[i] + 1 ) );
    }

    integral = integral / r8_gamma ( 0.5 * ( ( double ) ( e_sum + n ) ) );
  }

  integral = integral * pow ( r, e_sum + 2 );

  return integral;
}
//****************************************************************************80

double sphere_shell_03_nd ( double func ( int n, double x[] ), int n, 
  double center[], double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_SHELL_03_ND approximates an integral inside a spherical shell in ND.
//
//  Integration region:
//
//    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
//
//  Discussion:
//
//    An 2*N point 3-rd degree formula is used, Stroud number SN-Shell:3-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Input, double CENTER[N], the center of the spheres.
//
//    Input, double R1, R2, the inner and outer radiuses that
//    define the spherical shell.
//
//    Output, double SPHERE_SHELL_03_ND, the approximate integral of the function.
//
{
  int i;
  double quad;
  double r;
  double result;
  double rho;
  double volume;
  double w;
  double *x;

  if ( r1 == r2 )
  {
    result = 0.0;
    return result;
  }

  rho = r1 / r2;

  r = ( double ) ( n ) * ( 1.0 - pow ( rho, n + 2 ) ) 
    / ( ( double ) ( n + 2 ) * ( 1.0 - pow ( rho, n ) ) );
  r = sqrt ( r );
  w = 1.0 / ( double ) ( 2 * n );

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i];
  }

  quad = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x[i] = center[i] + r * r2;
    quad = quad + w * func ( n, x );
    x[i] = center[i] - r * r2;
    quad = quad + w * func ( n, x );
    x[i] = center[i];
  }

  volume = sphere_shell_volume_nd ( n, r1, r2 );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_shell_volume_nd ( int n, double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_SHELL_VOLUME_ND computes the volume of a spherical shell in ND.
//
//  Integration region:
//
//    R1*R1 <= sum ( X(1:N) - CENTER(1:N) )^2 <= R2*R2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, double R1, R2, the radiuses of the inner and 
//    outer spheres.
//
//    Output, double SPHERE_SHELL_VOLUME_ND, the volume of the
//    spherical shell.
//
{
  double volume;

  volume = ball_volume_nd ( n, r2 ) - ball_volume_nd ( n, r1 );

  return volume;
}
//****************************************************************************80

double sphere_unit_03_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_03_ND approximates an integral on the surface of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = 1.
//
//  Discussion:
//
//    A 2*N point 3rd degree formula is used, Stroud number UN:3-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Output, double SPHERE_UNIT_03_ND, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;
  double w;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  w = 1.0 / ( double ) ( 2 * n );

  quad = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0;
    quad = quad + w * func ( n, x );
    x[i] = -1.0;
    quad = quad + w * func ( n, x );
    x[i] = 0.0;
  }
  volume = sphere_unit_area_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_unit_04_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_04_ND approximates an integral on the surface of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = 1.
//
//  Discussion:
//
//    A 2*N*N point 5th degree formula is used, Stroud number UN:5-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Output, double SPHERE_UNIT_04_ND, the approximate integral of the function.
//
{
  int i;
  int j;
  double quad;
  double result;
  double s;
  double volume;
  double w1;
  double w2;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  w1 = ( double ) ( 4 - n ) / ( double ) ( 2 * n * ( n + 2 ) );

  quad = 0.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0;
    quad = quad + w1 * func ( n, x );
    x[i] = -1.0;
    quad = quad + w1 * func ( n, x );
    x[i] = 0.0;
  }

  s = 1.0 / sqrt ( 2.0 );
  w2 = 1.0 / ( double ) ( n * ( n + 2 ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = s;

    for ( j = i + 1; j < n; j++ )
    {
      x[j] = s;
      quad = quad + w2 * func ( n, x );
      x[j] = -s;
      quad = quad + w2 * func ( n, x );
      x[j] = 0.0;
    }

    x[i] = -s;

    for ( j = i + 1; j < n; j++ )
    {
      x[j] = s;
      quad = quad + w2 * func ( n, x );
      x[j] = -s;
      quad = quad + w2 * func ( n, x );
      x[j] = 0.0;
    }
    x[i] = 0.0;
  }
  volume = sphere_unit_area_nd ( n );
  result = quad * volume;

  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_unit_05_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_05_ND approximates an integral on the surface of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = 1.
//
//  Discussion:
//
//    A 2*N+2**N points 5-th degree formula is used, Stroud number UN:5-2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Output, double SPHERE_UNIT_05_ND, the approximate integral of the function.
//
{
  int i;
  int iadd;
  int ihi;
  int *ix;
  bool more;
  int ncard;
  double quad;
  double result;
  double volume;
  double w1;
  double w2;
  double *x;
  double x1;
  double x2;

  x1 = 1.0;
  x2 = 1.0 / sqrt ( ( double ) ( n ) );

  w1 = 1.0 / ( double ) ( n * ( n + 2 ) );
  w2 = ( double ) ( n ) / ( double ) ( ( n + 2 ) * i4_power ( 2, n ) );

  ix = new int[n];
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  quad = 0.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = x1;
    quad = quad + w1 * func ( n, x );
    x[i] = -x1;
    quad = quad + w1 * func ( n, x );
    x[i] = 0.0;
  }

  more = false;
  ihi = i4_power ( 2, n );

  for ( i = 0; i < n; i++ )
  {
    x[i] = -x2;
  }

  for ( i = 0; i < ihi; i++ )
  {
    subset_gray_next ( n, ix, &more, &ncard, &iadd );

    if ( iadd != -1 )
    {
      x[iadd-1] = -x[iadd-1];
    }
    quad = quad + w2 * func ( n, x );
  }
  volume = sphere_unit_area_nd ( n );
  result = quad * volume;

  delete [] ix;
  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_unit_07_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_07_3D approximates an integral on the surface of the unit sphere in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z = 1.
//
//  Discussion:
//
//    A 32 point 7-th degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Output, double SPHERE_UNIT_07_3D, the approximate integral of the function.
//
{
  double angle;
  int i;
  int j;
  int k;
  int order1 = 2;
  int order2 = 4;
  int order3 = 4;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double volume;
  double weight1[2];
  double weight2[4];
  double weight3[4];
  double x;
  double xtab1[2];
  double xtab2[4];
  double xtab3[4];
  double y;
  double z;
//
//  Set XTAB1 and WATE1.
//
  xtab1[0] = -1.0;
  xtab1[1] =  1.0;
  weight1[0] = 1.0;
  weight1[1] = 1.0;
//
//  Set XTAB2 and WATE2.
//
  for ( j = 0; j < order2; j++ )
  {
    angle = pi * ( double ) ( 2 * j + 1 ) / ( double ) ( 2 * order2 );
    xtab2[j] = cos ( angle );
  }

  for ( j = 0; j < order2; j++ )
  {
    weight2[j] = 1.0 / ( double ) ( 4 * order2 );
  }
//
//  Set XTAB3 and WATE3.
//
  legendre_set ( order3, xtab3, weight3 );

  quad = 0.0;
  for ( i = 0; i < order1; i++ )
  {
    for ( j = 0; j < order2; j++ )
    {
      for ( k = 0; k < order3; k++ )
      {
        x = xtab1[i] * sqrt ( 1.0 - xtab2[j] * xtab2[j] ) 
                     * sqrt ( 1.0 - xtab3[k] * xtab3[k] );
        y = xtab1[i] * xtab2[j] * sqrt ( 1.0 - xtab3[k] * xtab3[k] );
        z = xtab1[i] * xtab3[k];

        quad = quad + weight1[i] * weight2[j] * weight3[k] * func ( x, y, z );
      }
    }
  }
  volume = sphere_unit_area_3d ( );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double sphere_unit_07_1_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_07_1_ND approximates an integral on the surface of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = 1.
//
//  Discussion:
//
//    A 2**N + 2*N*N point 7th degree formula is used, Stroud number UN:7-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Output, double SPHERE_UNIT_07_1_ND, the approximate integral of the function.
//
{
  int i;
  int iadd;
  int *ix;
  int j;
  int jhi;
  bool more;
  int ncard;
  double quad;
  double result;
  double volume;
  double w1;
  double w2;
  double w3;
  double *x;
  double x1;
  double x2;
  double x3;

  ix = new int[n];
  x = new double[n];

  w1 = ( double ) ( 8 - n ) / ( double ) ( n * ( n + 2 ) * ( n + 4 ) );
  w2 = ( double ) ( n * n * n ) 
    / ( double ) ( i4_power ( 2, n ) * n * ( n + 2 ) * ( n + 4 ) );
  w3 = 4.0 / ( double ) ( n * ( n + 2 ) * ( n + 4 ) );

  x1 = 1.0;
  x2 = 1.0 / sqrt ( ( double ) ( n ) );
  x3 = 1.0 / sqrt ( 2.0 );

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  quad = 0.0;
//
//  First term.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = x1;
    quad = quad + w1 * func ( n, x );
    x[i] = -x1;
    quad = quad + w1 * func ( n, x );
    x[i] = 0.0;
  }
//
//  Second term.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = -x2;
  }
  more = false;
  jhi = i4_power ( 2, n );

  for ( j = 0; j < jhi; j++ )
  {
    subset_gray_next ( n, ix, &more, &ncard, &iadd );

    if ( iadd != -1 )
    {
      x[iadd-1] = -x[iadd-1];
    }
    quad = quad + w2 * func ( n, x );
  }
//
//  Third term.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      x[i] = x3;
      x[j] = x3;
      quad = quad + w3 * func ( n, x );
      x[i] = -x3;
      x[j] = x3;
      quad = quad + w3 * func ( n, x );
      x[i] = x3;
      x[j] = -x3;
      quad = quad + w3 * func ( n, x );
      x[i] = -x3;
      x[j] = -x3;
      quad = quad + w3 * func ( n, x );
      x[i] = 0.0;
      x[j] = 0.0;
    }
  }
  volume = sphere_unit_area_nd ( n );
  result = quad * volume;

  delete [] ix;
  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_unit_07_2_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_07_2_ND approximates an integral on the surface of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = 1.
//
//  Discussion:
//
//    A 2^N * ( N + 1 ) point 7th degree formula is used, Stroud number UN:7-2.
//
//    Some of the weights in this quadrature formula are negative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.
//
//    Output, double SPHERE_UNIT_07_2_ND, the approximate integral of the function.
//
{
  int iadd;
  int i;
  int *ix;
  int j;
  int jhi;
  bool more;
  int ncard;
  double quad;
  double result;
  double volume;
  double w1;
  double w2;
  double *x;
  double x1;
  double x2;
  double x3;

  ix = new int[n];
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }
  w1 = - ( double ) ( n * n ) 
    / ( double ) ( i4_power ( 2, n + 3 ) * ( n + 2 ) );
  w2 = ( double ) ( ( n + 4 ) * ( n + 4 ) ) 
    / ( double ) ( i4_power ( 2, n + 3 ) * n * ( n + 2 ) );
  x1 = 1.0 / sqrt ( ( double ) ( n ) );
  x2 = sqrt ( 5.0 / ( double ) ( n + 4 ) );
  x3 = 1.0 / sqrt ( ( double ) ( n + 4 ) );

  quad = 0.0;

  for ( j = 0; j < n; j++ )
  {
    x[j] = - x1;
  }
  more = false;
  jhi = i4_power ( 2, n );

  for ( j = 0; j < jhi; j++ )
  {
    subset_gray_next ( n, ix, &more, &ncard, &iadd );

    if ( iadd != -1 )
    {
      x[iadd-1] = - x[iadd-1];
    }
    quad = quad + w1 * func ( n, x );
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = - x3;
    }
    x[i] = - x2;
    more = false;

    for ( j = 0; j < jhi; j++ )
    {
      subset_gray_next ( n, ix, &more, &ncard, &iadd );

      if ( iadd != -1 )
      {
        x[iadd-1] = - x[iadd-1];
      }
      quad = quad + w2 * func ( n, x );
    }
  }
  volume = sphere_unit_area_nd ( n );
  result = quad * volume;

  delete [] ix;
  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_unit_11_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_11_3D approximates an integral on the surface of the unit sphere in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z = 1.
//
//  Discussion:
//
//    A 50 point 11-th degree formula is used, Stroud number U3:11-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    AD McLaren,
//    Mathematics of Computation,
//    Volume 17, pages 361-383, 1963.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Output, double SPHERE_UNIT_11_3D, the approximate integral of the function.
//
{
  int i;
  int j;
  int k;
  int l;
  double quad;
  double result;
  double volume;
  double w1;
  double w2;
  double w3;
  double w4;
  double x;
  double y;
  double z;

  quad = 0.0;

  w1 = 9216.0 / 725760.0;
  x = 1.0;
  y = 0.0;
  z = 0.0;

  for ( i = 0; i < 2; i++ )
  {
    x = -x;
    for ( j = 0; j < 3; j++ )
    {
      r8_swap3 ( &x, &y, &z );
      quad = quad + w1 * func ( x, y, z );
    }
  }

  w2 = 16384.0 / 725760.0;
  x = sqrt ( 0.5 );
  y = sqrt ( 0.5 );
  z = 0.0;

  for ( i = 0; i < 2; i++ )
  {
    x = -x;
    for ( j = 0; j < 2; j++ )
    {
      y = -y;
      for ( k = 0; k < 3; k++ )
      {
        r8_swap3 ( &x, &y, &z );
        quad = quad + w2 * func ( x, y, z );
      }
    }
  }

  w3 = 15309.0 / 725760.0;
  x = sqrt ( 1.0 / 3.0 );
  y = sqrt ( 1.0 / 3.0 );
  z = sqrt ( 1.0 / 3.0 );

  for ( i = 0; i < 2; i++ )
  {
    x = -x;
    for ( j = 0; j < 2; j++ )
    {
      y = -y;
      for ( k = 0; k < 2; k++ )
      {
        z = -z;
        quad = quad + w3 * func ( x, y, z );
      }
    }
  }

  w4 = 14641.0 / 725760.0;
  x = sqrt ( 1.0 / 11.0 );
  y = sqrt ( 1.0 / 11.0 );
  z = 3.0 * sqrt ( 1.0 / 11.0 );

  for ( i = 0; i < 2; i++ )
  {
    x = -x;
    for ( j = 0; j < 2; j++ )
    {
      y = -y;
      for ( k = 0; k < 2; k++ )
      {
        z = -z;
        for ( l = 0; l < 3; l++ )
        {
          r8_swap3 ( &x, &y, &z );
          quad = quad + w4 * func ( x, y, z );
        }
      }
    }
  }

  volume = sphere_unit_area_3d ( );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double sphere_unit_11_nd ( double func ( int n, double x[] ), int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_11_ND approximates an integral on the surface of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) = 1
//
//  Discussion:
//
//    An 2^N * ( N^2 + N + 1 ) point formula of degree 5 is used.
//
//    (For N = 3, the number of points is actually only 56, and
//     for N = 4, the number of points is actually only 240.)
//
//    One element of COEF31 was changed from
//      0.0236339091329 to
//      0.0236639091329
//    by Stroud, when going from his paper to his later textbook.
//    This correction was pointed out by David Wright, 16 February 2010.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    A Fifth Degree Integration Formula for the N-Simplex,
//    SIAM Journal on Numerical Analysis,
//    Volume 6, Number 1, March 1969.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( int n, double x[] ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int N, the dimension of the space.  For this routine,
//    it must be the case that 3 <= N <= 16.
//
//    Output, double SPHERE_UNIT_11_ND, the approximate integral of the function.
//
{
  double area;
  double coef1[16] = { 
    0.0,
    0.0,
    0.128571428571, 
    0.0518518518518, 
    0.0211979378646, 
    0.281250000000, 
    1.11934731935, 
    2.82751322751, 
    5.68266145619, 
    9.93785824515, 
    15.8196616478, 
    23.5285714285, 
    33.2409299392, 
    45.1113811729, 
    59.2754264177, 
    75.8518518518 };
  double coef21[16] = {
    0.0,
    0.0,
    0.163795782462, 
    0.0967270533860, 
    0.0638253880175, 
    0.0452340041459, 
    0.0336329118818, 
    0.0261275095270, 
    0.0208331595340, 
    0.0169937111647, 
    0.0141147212492, 
    0.0118949128383, 
    0.0101424250926, 
    0.00873046796644, 
    0.00757257014768, 
    0.00660813369775 };
  double coef22[16] = {
    0.0,
    0.0,
    0.126680408014, 
    0.0514210947621, 
    0.0213579471658, 
   -0.108726067638, 
   -0.371589499738, 
   -0.786048144448, 
   -1.36034060198, 
   -2.09547695631, 
   -2.98784764467, 
   -4.03107480702, 
   -5.21726499521, 
   -6.53783099707, 
   -7.98401677102, 
   -9.54722261180 };
  double coef31[16] = {
    0.0,
    0.0,
    0.0, 
    0.0592592592592, 
    0.0236639091329, 
    0.0525940190875, 
    0.0925052768546, 
    0.141316953438, 
    0.196818580052, 
    0.257027634179, 
    0.320299222258, 
    0.385326226441, 
    0.451098131789, 
    0.516849445559, 
    0.582010515746, 
    0.646165210110 };
  double coef32[16] = {
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0316246294890, 
    0.0207194729760, 
    0.0144303800811, 
    0.0105348984135, 
    0.00798435122193, 
    0.00623845929545, 
    0.00499896882962, 
    0.00409176297655, 
    0.00341037426698, 
    0.00288710646943, 
    0.00247745182907, 
    0.00215128820597 };
  int i;
  int iadd;
  int *ix;
  int j;
  int k;
  bool more;
  int ncard;
  double quad;
  double r1;
  double r2;
  double result;
  double s1;
  double s2;
  double u1;
  double u2;
  double v1;
  double v2;
  double *x;

  result = 0.0;

  if ( n < 3 || 16 < n )
  {
    cerr << "\n";
    cerr << "SPHERE_UNIT_11_ND - Fatal error!\n";
    cerr << "  Input spatial dimension N out of range.\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }

  ix = new int[n];
  x = new double[n];

  quad = 0.0;
//
//  S1
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0 / sqrt ( ( double ) ( n ) );
  }

  more = false;

  for ( ; ; )
  {
    subset_gray_next ( n, ix, &more, &ncard, &iadd );

    if ( iadd != -1 )
    {
      x[iadd-1] = -x[iadd-1];
    }
    quad = quad + coef1[n-1] * func ( n, x );

    if ( !more )
    {
      break;
    }
  }
//
//  S21
//
  r1 = ( ( double ) ( n + 6 ) - 4.0 * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 12 * n - 12 );
  r1 = sqrt ( r1 );

  s1 = ( ( double ) ( 7 * n - 6 ) 
    + ( double ) ( 4 * ( n - 1 ) ) * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 12 * n - 12 );
  s1 = sqrt ( s1 );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = r1;
    }
    x[i] = s1;

    more = false;

    for ( ; ; )
    {
      subset_gray_next ( n, ix, &more, &ncard, &iadd );

      if ( iadd != -1 )
      {
        x[iadd-1] = -x[iadd-1];
      }

      quad = quad + coef21[n-1] * func ( n, x );

      if ( !more )
      {
        break;
      }
    }
  }
//
//  S22
//
  r2 = ( ( double ) ( n + 6 ) + 4.0 * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 12 * n - 12 );
  r2 = sqrt ( r2 );

  s2 = ( ( double ) ( 7 * n - 6 ) 
    - ( double ) ( 4 * ( n - 1 ) ) * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 12 * n - 12 );
  s2 = sqrt ( s2 );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = r2;
    }
    x[i] = s2;

    more = false;

    for ( ; ; )
    {
      subset_gray_next ( n, ix, &more, &ncard, &iadd );

      if ( iadd != -1 )
      {
        x[iadd-1] = -x[iadd-1];
      }

      quad = quad + coef22[n-1] * func ( n, x );

      if ( !more )
      {
        break;
      }
    }
  }
//
//  S31
//
  u1 = ( ( double ) ( n + 12 ) + 8.0 * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 24 * n - 48 );
  u1 = sqrt ( u1 );

  v1 = ( ( double ) ( 7 * n - 12 ) 
    - ( double ) ( 4 * n - 8 ) * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 24 * n - 48 );
  v1 = sqrt ( v1 );

  for ( i = 0; i < n; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        x[k] = u1;
      }
      x[i] = v1;
      x[j] = v1;

      more = false;

      for ( ; ; )
      {
        subset_gray_next ( n, ix, &more, &ncard, &iadd );

        if ( iadd != -1 )
        {
          x[iadd-1] = -x[iadd-1];
        }

        quad = quad + coef31[n-1] * func ( n, x );

        if ( !more )
        {
          break;
        }
      }
    }
  }
//
//  S32
//
  u2 = ( ( double ) ( n + 12 ) - 8.0 * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 24 * n - 48 );
  u2 = sqrt ( u2 );

  v2 = ( ( double ) ( 7 * n - 12 ) 
    + ( double ) ( 4 * n - 8 ) * sqrt ( 3.0 ) ) 
    / ( double ) ( n * n + 24 * n - 48 );
  v2 = sqrt ( v2 );

  for ( i = 0; i < n; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        x[k] = u2;
      }
      x[i] = v2;
      x[j] = v2;

      more = false;

      for ( ; ; )
      {
        subset_gray_next ( n, ix, &more, &ncard, &iadd );

        if ( iadd != -1 )
        {
          x[iadd-1] = -x[iadd-1];
        }

        quad = quad + coef32[n-1] * func ( n, x );

        if ( !more )
        {
          break;
        }
      }
    }
  }

  area = sphere_unit_area_nd ( n );
  result = quad * area / pow ( 2.0, n );

  delete [] ix;
  delete [] x;

  return result;
}
//****************************************************************************80

double sphere_unit_14_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_14_3D approximates an integral on the surface of the unit sphere in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z = 1.
//
//  Discussion:
//
//    A 72 point 14-th degree formula is used, Stroud number U3:14-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    AD McLaren,
//    Mathematics of Computation,
//    Volume 17, pages 361-383, 1963.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Output, double SPHERE_UNIT_14_3D, the approximate integral of the function.
//
{
  int i;
  int j;
  int k;
  double quad;
  double result;
  double temp;
  double volume;
  double w1;
  double w2;
  double x;
  double xtab[5] = {
    -0.151108275, 0.315838353, 0.346307112, -0.101808787, -0.409228403 };
  double y;
  double ytab[5] = {
    0.155240600, 0.257049387, 0.666277790,  0.817386065, 0.501547712 };
  double z;
  double ztab[5] = {
    0.976251323, 0.913330032, 0.660412970,  0.567022920, 0.762221757 };

  quad = 0.0;

  w1 = 125.0 / 10080.0;
  x = 0.525731112;
  y = 0.850650808;
  z = 0.0;

  for ( i = 0; i < 2; i++ )
  {
    x = -x;
    for ( j = 0; j < 2; j++ )
    {
      y = -y;
      for ( k = 0; k < 3; k++ )
      {
        r8_swap3 ( &x, &y, &z );
        quad = quad + w1 * func ( x, y, z );
      }
    }
  }

  w2 = 143.0 / 10080.0;

  for ( i = 0; i < 5; i++ )
  {
    x = xtab[i];
    y = ytab[i];
    z = ztab[i];

    for ( j = 0; j < 3; j++ )
    {
      temp = x;
      x = z;
      z = -y;
      y = -temp;

      for ( k = 0; k < 3; k++ )
      {
        r8_swap3 ( &x, &y, &z );
        quad = quad + w2 * func ( x, y, z );
      }

      y = -y;
      z = -z;
      quad = quad + w2 * func ( x, y, z );
    }
  }

  volume = sphere_unit_area_3d ( );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double sphere_unit_15_3d ( double func ( double x, double y, double z ) )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_15_3D approximates an integral on the surface of the unit sphere in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z = 1.
//
//  Discussion:
//
//    A 128 point 15-th degree spherical product Gauss formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Output, double SPHERE_UNIT_15_3D, the approximate integral of the function.
//
{
  double angle;
  int i;
  int j;
  int k;
  int order = 8;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double volume;
  double *weight;
  double x;
  double *xtab;
  double y;
  double z;

  weight = new double[order];
  xtab = new double[order];

  legendre_set ( order, xtab, weight );

  for ( i = 0; i < order; i++ )
  {
    weight[i] = weight[i] / 32.0;
  }

  quad = 0.0;

  for ( j = 0; j < order; j++ )
  {
    for ( k = 0; k < 16; k++ )
    {
      angle = ( double ) ( k ) * pi / 8.0;
      x = sqrt ( 1.0 - xtab[j] * xtab[j] ) * cos ( angle );
      y = sqrt ( 1.0 - xtab[j] * xtab[j] ) * sin ( angle );
      z = xtab[j];

      quad = quad + weight[j] * func ( x, y, z );
    }
  }
  volume = sphere_unit_area_3d ( );
  result = quad * volume;

  delete [] weight;
  delete [] xtab;

  return result;
}
//****************************************************************************80

double sphere_unit_area_3d ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_AREA_3D computes the surface area of the unit sphere in 3D.
//
//  Integration region:
//
//    X*X + Y*Y + Z*Z = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double SPHERE_UNIT_AREA_3D, the area of the sphere.
//
{
  double area;
  double pi = 3.141592653589793;

  area = 4.0 * pi;

  return area;
}
//****************************************************************************80

double sphere_unit_area_nd ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_AREA_ND computes the surface area of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( ( X(1:N) - CENTER(1:N) )^2 ) = R * R.
//
//  Discussion:
//
//    N   Area
//
//    2   2       * PI
//    3   4       * PI
//    4   2       * PI^2
//    5   (8/3)   * PI^2
//    6             PI^3
//    7   (16/15) * PI^3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Output, double SPHERE_UNIT_AREA_ND, the area of the sphere.
//
{
  double area;
  int i;
  int m;
  double pi = 3.141592653589793;

  if ( ( n % 2 ) == 0 )
  {
    m = n / 2;
    area = 2.0 * pow ( pi, m );
    for ( i = 1; i <= m - 1; i++ )
    {
      area = area / ( double ) ( i );
    }
  }
  else
  {
    m = ( n - 1 ) / 2;
    area = pow ( 2.0, n ) * pow ( pi, m );
    for ( i = m + 1; i <= 2 * m; i++ )
    {
      area = area / ( double ) ( i );
    }
  }

  return area;
}
//****************************************************************************80

void sphere_unit_area_values ( int *n_data, int *n, double *area )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
//
//  Discussion:
//
//    The formula for the surface area of the unit sphere in N dimensions is:
//
//      Sphere_Unit_Area ( N ) = 2 * PI**(N/2) / Gamma ( N / 2 )
//
//    Some values of the function include:
//
//       N   Area
//
//       2    2        * PI
//       3  ( 4 /    ) * PI
//       4  ( 2 /   1) * PI^2
//       5  ( 8 /   3) * PI^2
//       6  ( 1 /   1) * PI^3
//       7  (16 /  15) * PI^3
//       8  ( 1 /   3) * PI^4
//       9  (32 / 105) * PI^4
//      10  ( 1 /  12) * PI^5
//
//    For the unit sphere, Area(N) = N * Volume(N)
//
//    In Mathematica, the function can be evaluated by:
//
//      2 * Pi^(n/2) / Gamma[n/2]
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the spatial dimension.
//
//    Output, double *AREA, the area of the unit sphere 
//    in that dimension.
//
{
# define N_MAX 20

  double area_vec[N_MAX] = { 
     0.2000000000000000E+01,  
     0.6283185307179586E+01,  
     0.1256637061435917E+02,  
     0.1973920880217872E+02,  
     0.2631894506957162E+02,  
     0.3100627668029982E+02,  
     0.3307336179231981E+02,  
     0.3246969701133415E+02,  
     0.2968658012464836E+02,  
     0.2550164039877345E+02,  
     0.2072514267328890E+02,  
     0.1602315322625507E+02,  
     0.1183817381218268E+02,  
     0.8389703410491089E+01,  
     0.5721649212349567E+01,  
     0.3765290085742291E+01,  
     0.2396678817591364E+01,  
     0.1478625959000308E+01,  
     0.8858104195716824E+00,  
     0.5161378278002812E+00 };

  int n_vec[N_MAX] = { 
     1, 
     2, 
     3, 
     4, 
     5, 
     6, 
     7, 
     8, 
     9, 
    10, 
    11, 
    12, 
    13, 
    14, 
    15, 
    16, 
    17, 
    18, 
    19, 
    20 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *area = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *area = area_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double sphere_unit_monomial_nd ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_MONOMIAL_ND integrates a monomial on the surface of the unit sphere in ND.
//
//  Integration region:
//
//    sum ( X(1:N)^2 ) == 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gerald Folland,
//    How to Integrate a Polynomial Over a Sphere,
//    American Mathematical Monthly,
//    Volume 108, May 2001, pages 446-448.
//
//  Parameters:
//
//    Input, int N, the dimension of the space.
//
//    Input, int P[N], the exponents of X(1) through X(N) in the monomial.
//    The exponents P(N) must be nonnegative.
//
//    Output, double SPHERE_UNIT_MONOMIAL_ND, the integral of
//    X1**P(1)*X2**P(2)*...*XN**P(N) over the unit sphere.
//
{
  double arg1;
  double arg2;
  int i;
  double temp;
  double value;

  for ( i = 0; i < n; i++ )
  {
    if ( p[i] % 2 == 1 )
    {
      value = 0.0;
      return value;
    }
  }

  temp = 0.0;
  arg2 = 0.0;

  for ( i = 0; i < n; i++ )
  {
    arg1 = ( double ) ( p[i] + 1 ) / 2.0;
    temp = temp + r8_gamma_log ( arg1 );
    arg2 = arg2 + arg1;
  }
  temp = temp - r8_gamma_log ( arg2 );
  
  value = 2.0 * exp ( temp );

  return value;
}
//****************************************************************************80

double sphere_unit_volume_nd ( int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
//
//  Discussion:
//
//    The unit sphere in ND satisfies:
//
//      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
//
//    Results for the first few values of DIM_NUM are:
//
//     DIM_NUM  Volume
//
//     1    2
//     2    1        * PI
//     3  ( 4 /   3) * PI
//     4  ( 1 /   2) * PI^2
//     5  ( 8 /  15) * PI^2
//     6  ( 1 /   6) * PI^3
//     7  (16 / 105) * PI^3
//     8  ( 1 /  24) * PI^4
//     9  (32 / 945) * PI^4
//    10  ( 1 / 120) * PI^5
//
//    For the unit sphere, Volume(DIM_NUM) = 2 * PI * Volume(DIM_NUM-2)/ DIM_NUM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Output, double SPHERE_UNIT_VOLUME_ND, the volume of the sphere.
//
{
  int i;
  int m;
  double pi = 3.141592653589793;
  double volume;

  if ( ( dim_num % 2 ) == 0 )
  {
    m = dim_num / 2;
    volume = pow ( pi, m );
    for ( i = 1; i <= m; i++ )
    {
      volume = volume / ( double ) ( i );
    }
  }
  else
  {
    m = ( dim_num - 1 ) / 2;
    volume = pow ( pi, m ) * pow ( 2.0, dim_num );
    for ( i = m + 1; i <= 2 * m + 1; i++ )
    {
      volume = volume / ( double ) ( i );
    }
  }
  return volume;
}
//****************************************************************************80

void sphere_unit_volume_values ( int *n_data, int *n, double *volume )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
//
//  Discussion:
//
//    The formula for the volume of the unit sphere in N dimensions is
//
//      Volume(N) = 2 * PI**(N/2) / ( N * Gamma ( N / 2 ) )
//
//    This function satisfies the relationships:
//
//      Volume(N) = 2 * PI * Volume(N-2) / N
//      Volume(N) = Area(N) / N
//
//    Some values of the function include:
//
//       N  Volume
//
//       1    1
//       2    1        * PI
//       3  ( 4 /   3) * PI
//       4  ( 1 /   2) * PI^2
//       5  ( 8 /  15) * PI^2
//       6  ( 1 /   6) * PI^3
//       7  (16 / 105) * PI^3
//       8  ( 1 /  24) * PI^4
//       9  (32 / 945) * PI^4
//      10  ( 1 / 120) * PI^5
//
//    In Mathematica, the function can be evaluated by:
//
//      2 * Pi^(n/2) / ( n * Gamma[n/2] )
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the spatial dimension.
//
//    Output, double *VOLUME, the volume of the unit 
//    sphere in that dimension.
//
{
# define N_MAX 20

  int n_vec[N_MAX] = { 
     1,  2, 
     3,  4, 
     5,  6, 
     7,  8, 
     9, 10, 
    11, 12, 
    13, 14, 
    15, 16, 
    17, 18, 
    19, 20 };

  double volume_vec[N_MAX] = { 
     0.2000000000000000E+01,   
     0.3141592653589793E+01,  
     0.4188790204786391E+01,  
     0.4934802200544679E+01,  
     0.5263789013914325E+01,  
     0.5167712780049970E+01,  
     0.4724765970331401E+01,  
     0.4058712126416768E+01,  
     0.3298508902738707E+01,  
     0.2550164039877345E+01,  
     0.1884103879389900E+01,   
     0.1335262768854589E+01,  
     0.9106287547832831E+00,  
     0.5992645293207921E+00,  
     0.3814432808233045E+00,  
     0.2353306303588932E+00,  
     0.1409811069171390E+00,  
     0.8214588661112823E-01,  
     0.4662160103008855E-01,  
     0.2580689139001406E-01  };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *volume = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *volume = volume_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double sphere_volume_2d ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_VOLUME_2D computes the volume of an implicit sphere in 2D.
//
//  Discussion:
//
//    An implicit sphere in 2D satisfies the equation:
//
//      sum ( ( P(1:DIM_NUM) - CENTER(1:DIM_NUM) )^2 ) = R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_VOLUME_2D, the volume of the sphere.
//
{
  double pi = 3.141592653589793;
  double volume;

  volume = pi * r * r;

  return volume;
}
//****************************************************************************80

double sphere_volume_3d ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_VOLUME_3D computes the volume of an implicit sphere in 3D.
//
//  Discussion:
//
//    An implicit sphere in 3D satisfies the equation:
//
//      sum ( ( P(1:DIM_NUM) - CENTER(1:DIM_NUM) )^2 ) = R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_VOLUME_3D, the volume of the sphere.
//
{
  double pi = 3.141592653589793;
  double volume;

  volume = ( 4.0 / 3.0 ) * pi * r * r * r;

  return volume;
}
//****************************************************************************80

double sphere_volume_nd ( int dim_num, double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_VOLUME_ND computes the volume of an implicit sphere in ND.
//
//  Discussion:
//
//    An implicit sphere in ND satisfies the equation:
//
//      sum ( ( X(1:N) - CENTER(1:N) )^2 ) = R * R
//
//    where R is the radius and CENTER is the center.
//
//    Results for the first few values of N are:
//
//    DIM_NUM  Volume
//    -     -----------------------
//    2                PI   * R^2
//    3     (4/3)    * PI   * R^3
//    4     (1/2)    * PI^2 * R^4
//    5     (8/15)   * PI^2 * R^5
//    6     (1/6)    * PI^3 * R^6
//    7     (16/105) * PI^3 * R^7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_VOLUME_ND, the volume of the sphere.
//
{
  double volume;

  volume = pow ( r, dim_num ) * sphere_unit_volume_nd ( dim_num );

  return volume;
}
//****************************************************************************80

double square_sum ( double func ( double x, double y ), double center[2], 
  double r, int order, double xtab[], double ytab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_SUM carries out a quadrature rule over a square.
//
//  Integration region:
//
//      abs ( X - CENTER(1) ) <= R 
//    and
//      abs ( Y - CENTER(2) ) <= R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the function 
//    to be integrated.  
//
//    Input, double CENTER[2], the center of the square.
//
//    Input, double R, the radius of the square.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas of 
//    the rule.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double SQUARE_SUM, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;
  double x;
  double y;

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    x = center[0] + r * xtab[i];
    y = center[1] + r * ytab[i];
    quad = quad + 0.25 * weight[i] * func ( x, y );
  }

  volume = 4.0 * r * r;
  result = quad * volume;

  return result;
}
//****************************************************************************80

void square_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_UNIT_SET sets quadrature weights and abscissas in the unit square.
//
//  Discussion;
//
//    To get the value of ORDER associated with a given rule, 
//    call SQUARE_UNIT_SIZE first.
//
//  Integration region:
//
//      -1 <= X <= 1,
//    and
//      -1 <= Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gilbert Strang, George Fix,
//    An Analysis of the Finite Element Method,
//    Cambridge, 1973,
//    ISBN: 096140888X,
//    LC: TA335.S77.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule number.
//    1, order 1, degree 1 rule.
//    2, order 4, degree 3, rule.
//    3, order 9, degree 5 rule.
//    4, order 12 degree 7 rule, Stroud number C2:7-1.
//    5, order 13 degree 7 rule, Stroud number C2:7-3.
//    6, order 64 degree 15 product rule.
//
//    Input, int ORDER, the order of the rule.
//
//    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double a;
  double c;
  int i;
  int j;
  int k;
  int order2 = 8;
  double r;
  double s;
  double t;
  double w1;
  double w2;
  double w3;
  double *weight2;
  double *xtab2;
  double z;

  if ( rule == 1 )
  {
    weight[0] = 4.0;

    xtab[0] = 0.0;
    ytab[0] = 0.0;
  }
  else if ( rule == 2 )
  {
    a = 1.0;
    s = 1.0 / sqrt ( 3.0 );

    xtab[0] = - s;
    xtab[1] =   s;
    xtab[2] = - s;
    xtab[3] =   s;

    ytab[0] = - s;
    ytab[1] = - s;
    ytab[2] =   s;
    ytab[3] =   s;

    weight[0] = a;
    weight[1] = a;
    weight[2] = a;
    weight[3] = a;
  }
  else if ( rule == 3 )
  {
    s = sqrt ( 0.6 );
    z = 0.0;
    w1 = 64.0 / 81.0;
    w2 = 25.0 / 81.0;
    w3 = 40.0 / 81.0;

    xtab[0] =   z;
    xtab[1] = - s;
    xtab[2] =   s;
    xtab[3] = - s;
    xtab[4] =   s;
    xtab[5] =   z;
    xtab[6] = - s;
    xtab[7] =   s;
    xtab[8] =   z;

    ytab[0] =   z;
    ytab[1] = - s;
    ytab[2] = - s;
    ytab[3] =   s;
    ytab[4] =   s;
    ytab[5] = - s;
    ytab[6] =   z;
    ytab[7] =   z;
    ytab[8] =   s;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w2;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w3;
    weight[8] = w3;
  }
  else if ( rule == 4 )
  {
    r = sqrt ( 6.0 / 7.0 );
    c = 3.0 * sqrt ( 583.0 );
    s = sqrt ( ( 114.0 - c ) / 287.0 );
    t = sqrt ( ( 114.0 + c ) / 287.0 );
    w1 = 4.0 * 49.0 / 810.0;
    w2 = 4.0 * ( 178981.0 + 923.0 * c ) / 1888920.0;
    w3 = 4.0 * ( 178981.0 - 923.0 * c ) / 1888920.0;
    z = 0.0;

    xtab[0]  =   r;
    xtab[1]  =   z;
    xtab[2]  = - r;
    xtab[3]  =   z;
    xtab[4]  =   s;
    xtab[5]  = - s;
    xtab[6]  = - s;
    xtab[7]  =   s;
    xtab[8]  =   t;
    xtab[9]  = - t;
    xtab[10] = - t;
    xtab[11] =   t;

    ytab[0]  =   z;
    ytab[1]  =   r;
    ytab[2]  =   z;
    ytab[3]  = - r;
    ytab[4]  =   s;
    ytab[5]  =   s;
    ytab[6]  = - s;
    ytab[7]  = - s;
    ytab[8]  =   t;
    ytab[9]  =   t;
    ytab[10] = - t;
    ytab[11] = - t;

    weight[0]  = w1;
    weight[1]  = w1;
    weight[2]  = w1;
    weight[3]  = w1;
    weight[4]  = w2;
    weight[5]  = w2;
    weight[6]  = w2;
    weight[7]  = w2;
    weight[8]  = w3;
    weight[9]  = w3;
    weight[10] = w3;
    weight[11] = w3;
  }
  else if ( rule == 5 )
  {
    r = sqrt ( 12.0 / 35.0 );
    c = 3.0 * sqrt ( 186.0 );
    s = sqrt ( ( 93.0 + c ) / 155.0 );
    t = sqrt ( ( 93.0 - c ) / 155.0 );
    w1 =  8.0 / 162.0;
    w2 = 98.0 / 162.0;
    w3 = 31.0 / 162.0;
    z = 0.0;

    xtab[0]  =   z;
    xtab[1]  =   r;
    xtab[2]  = - r;
    xtab[3]  =   z;
    xtab[4]  =   z;
    xtab[5]  =   s;  
    xtab[6]  =   s;
    xtab[7]  = - s;
    xtab[8]  = - s;
    xtab[9]  =   t;
    xtab[10] =   t;
    xtab[11] = - t;
    xtab[12] = - t;

    ytab[0]  =   z;
    ytab[1]  =   z;
    ytab[2]  =   z;
    ytab[3]  =   r;
    ytab[4]  = - r;
    ytab[5]  =   t;
    ytab[6]  = - t;
    ytab[7]  =   t;
    ytab[8]  = - t;
    ytab[9]  =   s;
    ytab[10] = - s;
    ytab[11] =   s;
    ytab[12] = - s;

    weight[0]  = w1;
    weight[1]  = w2;
    weight[2]  = w2;
    weight[3]  = w2;
    weight[4]  = w2;
    weight[5]  = w3;
    weight[6]  = w3;
    weight[7]  = w3;
    weight[8]  = w3;
    weight[9]  = w3;
    weight[10] = w3;
    weight[11] = w3;
    weight[12] = w3;
  }
  else if ( rule == 6 )
  {
    xtab2 = new double[order2];
    weight2 = new double[order2];

    legendre_set ( order2, xtab2, weight2 );

    k = 0;

    for ( i = 0; i < order2; i++ )
    {
      for ( j = 0; j < order2; j++ )
      {
        xtab[k] = xtab2[i];
        ytab[k] = xtab2[j];
        weight[k] = weight2[i] * weight2[j];
        k = k + 1;
      }
    }
    delete [] xtab2;
    delete [] weight2;
  }
  else
  {
    cerr << "\n";
    cerr << "SQUARE_UNIT_SET - Fatal error!\n";
    cerr << "  Illegal value of RULE = " << rule << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int square_unit_size ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_UNIT_SIZE sizes a quadrature rule in the unit square.
//
//  Integration region:
//
//      -1 <= X <= 1,
//    and
//      -1 <= Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gilbert Strang, George Fix,
//    An Analysis of the Finite Element Method,
//    Cambridge, 1973,
//    ISBN: 096140888X,
//    LC: TA335.S77.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the rule number.
//    1, a 1 point 1st degree rule.
//    2, a 4 point 3rd degree rule.
//    3, a 9 point 5th degree rule.
//    4, a 12 point 7-th degree rule, Stroud number C2:7-1.
//    5, a 13 point 7-th degree rule, Stroud number C2:7-3.
//    6, a 64 point 15-th degree product rule.
//
//    Output, int SQUARE_UNIT_SIZE, the order of the rule.
//
{
  int order;

  if ( rule == 1 ) 
  {
    order = 1;
  }
  else if ( rule == 2 )
  {
    order = 4;
  }
  else if ( rule == 3 )
  {
    order = 9;
  }
  else if ( rule == 4 )
  {
    order = 12;
  }
  else if ( rule == 5 )
  {
    order = 13;
  }
  else if ( rule == 6 )
  {
    order = 64;
  }
  else
  {
    order = - 1;
  }

  return order;
}
//****************************************************************************80

double square_unit_sum ( double func ( double x, double y ), int order, 
  double xtab[], double ytab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_UNIT_SUM carries out a quadrature rule over the unit square.
//
//  Integration region:
//
//      -1 <= X <= 1, 
//    and
//      -1 <= Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the function
//    to be integrated. 
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights.
//
//    Output, double SQUARE_UNIT_SUM, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + weight[i] * func ( xtab[i], ytab[i] ) / 4.0;
  }

  volume = 1.0;
  result = quad * volume;

  return result;
}
//****************************************************************************80

void subset_gray_next ( int n, int a[], bool *more, int *ncard, int *iadd )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
//
//  Discussion:
//
//    It generates the subsets one at a time, by adding or subtracting
//    exactly one element on each step.
//
//    The user should set MORE = .FALSE. and the value of N before
//    the first call.  On return, the user may examine A which contains
//    the definition of the new subset, and must check .MORE., because
//    as soon as it is .FALSE. on return, all the subsets have been
//    generated and the user probably should cease calling.
//
//    The first set returned is the empty set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2003
//
//  Author:
//
//    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the total set from which
//    subsets will be drawn.
//
//    Input/output, int A[N].  On each return, the Gray code for the newly
//    generated subset.  A[I] = 0 if element I is in the subset, 1 otherwise.
//
//    Input/output, bool *MORE.  Set this variable FALSE before
//    the first call.  Normally, MORE will be returned TRUE but once
//    all the subsets have been generated, MORE will be
//    reset FALSE on return and you should stop calling the program.
//
//    Input/output, int *NCARD, the cardinality of the set returned,
//    which may be any value between 0 (the empty set) and N (the
//    whole set).
//
//    Output, int *IADD, the element which was added or removed to the
//    previous subset to generate the current one.  Exception:
//    the empty set is returned on the first call, and IADD is set to -1.
{
  int i;
//
//  First set returned is the empty set.
//
  if ( !(*more) )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0;
    }

    *iadd = 0;
    *ncard = 0;
    *more = true;
  }
  else
  {
    *iadd = 1;

    if ( ( *ncard % 2 ) != 0 )
    {
      for ( ; ; )
      {
        *iadd = *iadd + 1;
        if ( a[*iadd-2] != 0 )
        {
          break;
        }

      }

    }

    a[*iadd-1] = 1 - a[*iadd-1];
    *ncard = *ncard + 2 * a[*iadd-1] - 1;
//
//  Last set returned is the singleton A(N).
//
    if ( *ncard == a[n-1] )
    {
      *more = false;
    }

  }

  return;
}
//****************************************************************************80

double tetra_07 ( double func ( double x, double y, double z ), double x[], 
  double y[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_07 approximates an integral inside a tetrahedron in 3D.
//
//  Integration region:
//
//    Points inside a tetrahedron whose four corners are given.
//
//  Discussion:
//
//    A 64 point 7-th degree conical product Gauss formula is used,
//    Stroud number T3:7-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966, pages 42-43,
//    LC: QA299.4G3S7
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double X[4], Y[4], Z[4], the coordinates of 
//    the vertices.
//
//    Output, double TETRAQ_07, the approximate integral of the function.
//
{
  double a;
  double b;
  double c;
  double d;
  int i;
  int j;
  int k;
  int order = 4;
  double quad;
  double result;
  double t;
  double u;
  double v;
  double volume;
  double w;
  double weight1[4];
  double weight2[4] = {
    0.1355069134, 0.2034645680, 0.1298475476, 0.0311809709 };
  double weight3[4] = {
    0.1108884156, 0.1434587898, 0.0686338872, 0.0103522407 };
  double xtab1[4];
  double xtab2[4] = {
    0.0571041961, 0.2768430136, 0.5835904324, 0.8602401357 };
  double xtab3[4] = {
    0.0485005495, 0.2386007376, 0.5170472951, 0.7958514179 };
  double xval;
  double yval;
  double zval;
//
//  Get the Gauss-Legendre weights and abscissas for [-1,1].
//
  legendre_set ( order, xtab1, weight1 );
//
//  Adjust the rule for the interval [0,1].
//
  a = -1.0;
  b = +1.0;

  c =  0.0;
  d =  1.0;

  rule_adjust ( a, b, c, d, order, xtab1, weight1 );
//
//  Carry out the quadrature.
//
  quad = 0.0;

  for ( i = 0; i < order; i++ )
  {
    for ( j = 0; j < order; j++ )
    {
      for ( k = 0; k < order; k++ )
      {
//
//  Compute the barycentric coordinates of the point in the unit triangle.
//
        t =                                         xtab3[k];
        u =                    xtab2[j]   * ( 1.0 - xtab3[k] );
        v = xtab1[i] * ( 1.0 - xtab2[j] ) * ( 1.0 - xtab3[k] );
        w = 1.0 - t - u - v;
//
//  Compute the corresponding point in the triangle.
//
        xval = t * x[0] + u * x[1] + v * x[2] + w * x[3];
        yval = t * y[0] + u * y[1] + v * y[2] + w * y[3];
        zval = t * z[0] + u * z[1] + v * z[2] + w * z[3];

        quad = quad + 6.0 * weight1[i] * weight2[j] * weight3[k] 
          * func ( xval, yval, zval );
      }
    }
  }
  volume = tetra_volume ( x, y, z );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double tetra_sum ( double func ( double x, double y, double z ), double x[4], 
  double y[4], double z[4], int order, double xtab[], double ytab[], 
  double ztab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_SUM carries out a quadrature rule in a tetrahedron in 3D.
//
//  Integration region:
//
//    A tetrahedron whose vertices are specified.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the
//    user supplied function which is to be integrated.
//
//    Input, double X[4], Y[4], Z[4], the vertices.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], ZTAB[ORDER], the
//    abscissas.
//
//    Input, double WEIGHT[ORDER], the weights.
//
//    Output, double TETRA_SUM, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;
  double xval;
  double yval;
  double zval;

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    xval =         xtab[i]                       * x[0] 
                           + ytab[i]             * x[1] 
                                     + ztab[i]   * x[2] 
         + ( 1.0 - xtab[i] - ytab[i] - ztab[i] ) * x[3];

    yval =         xtab[i]                       * y[0] 
                           + ytab[i]             * y[1] 
                                     + ztab[i]   * y[2] 
         + ( 1.0 - xtab[i] - ytab[i] - ztab[i] ) * y[3];

    zval =         xtab[i]                       * z[0] 
                           + ytab[i]             * z[1] 
                                     + ztab[i]   * z[2] 
         + ( 1.0 - xtab[i] - ytab[i] - ztab[i] ) * z[3];

    quad = quad + weight[i] * func ( xval, yval, zval );
  }
  volume = tetra_volume ( x, y, z );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double tetra_tproduct ( double func ( double x, double y, double z ), 
  int order, double x[4], double y[4], double z[4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_TPRODUCT approximates an integral in a tetrahedron in 3D.
//
//  Discussion:
//
//    Integration is carried out over the points inside an arbitrary
//    tetrahedron whose four vertices are given.
//
//    An ORDER**3 point (2*ORDER-1)-th degree triangular product
//    Gauss-Legendre rule is used.
//
//    With ORDER = 8, this routine is equivalent to the routine TETR15
//    in the reference, page 367.
//
//    Thanks to Joerg Behrens, jbehren@gwdg.de, for numerous suggestions
//    and corrections.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int ORDER, the order of the basic quadrature rules.
//    ORDER should be between 1 and 9.
//
//    Input, double X[4], Y[4], Z[4], the vertices
//    of the tetrahedron.
//
//    Output, double TETRA_TPRODUCT, the approximate integral of the function.
//
{
  double a;
  double b;
  double c;
  double d;
  int i;
  int j;
  int k;
  double quad;
  double result;
  double volume;
  double *weight0;
  double *weight1;
  double *weight2;
  double *xtab0;
  double *xtab1;
  double *xtab2;
  double xval;
  double yval;
  double zval;

  if ( order < 1 || 9 < order )
  {
    cerr << "\n";
    cerr << "TETRA_TPRODUCT - Fatal error!\n";
    cerr << "  The quadrature rule orders must be between 1 and 9.\n";
    cerr << "  The input value was ORDER = " << order << "\n";
    exit ( 1 );
  }
//
//  Get the Gauss-Legendre ORDER point rules on [-1,1] for integrating
//    F(X),
//    X * F(X),
//    X * X * F(X).
//
  xtab0 = new double[order];
  xtab1 = new double[order];
  xtab2 = new double[order];
  weight0 = new double[order];
  weight1 = new double[order];
  weight2 = new double[order];

  legendre_set ( order, xtab0, weight0 );
  legendre_set_x1 ( order, xtab1, weight1 );
  legendre_set_x2 ( order, xtab2, weight2 );
//
//  Adjust the rules from [-1,1] to [0,1].
//
  a = -1.0;
  b = +1.0;
  c =  0.0;
  d =  1.0;

  rule_adjust ( a, b, c, d, order, xtab0, weight0 );

  rule_adjust ( a, b, c, d, order, xtab1, weight1 );

  rule_adjust ( a, b, c, d, order, xtab2, weight2 );
//
//  For rules with a weight function that is not 1, the weight vectors
//  require further adjustment.
//
  for ( i = 0; i < order; i++ )
  {
    weight1[i] = weight1[i] / 2.0;
  }
  for ( i = 0; i < order; i++ )
  {
    weight2[i] = weight2[i] / 4.0;
  }
//
//  Carry out the quadrature.
//
  quad = 0.0;

  for ( k = 0; k < order; k++ )
  {
    for ( j = 0; j < order; j++ )
    {
      for ( i = 0; i < order; i++ )
      {
        xval = x[0] + ( ( ( x[3] - x[2] )   * xtab0[i] 
                        + ( x[2] - x[1] ) ) * xtab1[j] 
                        + ( x[1] - x[0] ) ) * xtab2[k];

        yval = y[0] + ( ( ( y[3] - y[2] )   * xtab0[i] 
                        + ( y[2] - y[1] ) ) * xtab1[j] 
                        + ( y[1] - y[0] ) ) * xtab2[k];

        zval = z[0] + ( ( ( z[3] - z[2] )   * xtab0[i] 
                        + ( z[2] - z[1] ) ) * xtab1[j] 
                        + ( z[1] - z[0] ) ) * xtab2[k];

        quad = quad + 6.0 * weight0[i] * weight1[j] * weight2[k] 
          * func ( xval, yval, zval );
      }
    }
  }
//
//  Compute the volume of the tetrahedron.
//
  volume = tetra_volume ( x, y, z );
  result = quad * volume;

  delete [] xtab0;
  delete [] xtab1;
  delete [] xtab2;
  delete [] weight0;
  delete [] weight1;
  delete [] weight2;

  return result;
}
//****************************************************************************80

void tetra_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double ztab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_UNIT_SET sets quadrature weights and abscissas in the unit tetrahedron.
//
//  Integration region:
//
//      0 <= X,
//    and
//      0 <= Y,
//    and
//      0 <= Z, 
//    and
//      X + Y + Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hermann Engels,
//    Numerical Quadrature and Cubature,
//    Academic Press, 1980,
//    ISBN: 012238850X,
//    LC: QA299.3E5.
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//    ISBN: 0750663200,
//    LC: TA640.2.Z54
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
//     2, order 4, precision 1, Newton Cotes formula #1.
//     3, order 4, precision 2, Zienkiewicz #2.
//     4, order 10, precision 2, Newton Cotes formula #2
//     5, order 5, precision 3, Zienkiewicz #3.
//     6, order 8, precision 3, Newton Cotes formula #3.
//     7, order 35, precision 4, Newton Cotes formula #4.
//     8, order 11, precision 4, a Keast rule.
//
//    Input, int ORDER, the order of the rule.
//
//    Output, double XTAB[ORDER], YTAB[ORDER], ZTAB[ORDER],
//    the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  double g;
  double h;
  double z;
//
//  Newton Cotes #0.
//
  if ( rule == 1 )
  {
    xtab[0] = 0.25;
    ytab[0] = 0.25;
    ztab[0] = 0.25;
    weight[0] = 1.0;
  }
//
//  Newton Cotes #1.
//
  else if ( rule == 2 )
  {
    a = 1.0;
    b = 1.0 / 4.0;
    z = 0.0;

    xtab[0] = z;
    xtab[1] = a;
    xtab[2] = z;
    xtab[3] = z;

    ytab[0] = z;
    ytab[1] = z;
    ytab[2] = a;
    ytab[3] = z;

    ztab[0] = z;
    ztab[1] = z;
    ztab[2] = z;
    ztab[3] = a;

    weight[0] = b;
    weight[1] = b;
    weight[2] = b;
    weight[3] = b;
  }
//
//  Zienkiewicz #2.
//
  else if ( rule == 3 )
  {
    a =  0.5854101966249685;
    b =  0.1381966011250105;
    c =  0.25;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = b;
    xtab[3] = b;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = b;

    ztab[0] = b;
    ztab[1] = b;
    ztab[2] = a;
    ztab[3] = b;

    weight[0] = c;
    weight[1] = c;
    weight[2] = c;
    weight[3] = c;
  }
//
//  Newton Cotes #2.
//
  else if ( rule == 4 )
  {
    a =  1.0;
    b =  0.5;
    c = -1.0 / 20.0;
    d =  4.0 / 20.0;
    z =  0.0;

    xtab[0] = z;
    xtab[1] = a;
    xtab[2] = z;
    xtab[3] = z;
    xtab[4] = b;
    xtab[5] = z;
    xtab[6] = z;
    xtab[7] = b;
    xtab[8] = b;
    xtab[9] = z;

    ytab[0] = z;
    ytab[1] = z;
    ytab[2] = a;
    ytab[3] = z;
    ytab[4] = z;
    ytab[5] = b;
    ytab[6] = z;
    ytab[7] = b;
    ytab[8] = z;
    ytab[9] = b;

    ztab[0] = z;
    ztab[1] = z;
    ztab[2] = z;
    ztab[3] = a;
    ztab[4] = z;
    ztab[5] = z;
    ztab[6] = b;
    ztab[7] = z;
    ztab[8] = b;
    ztab[9] = b;

    weight[0] = c;
    weight[1] = c;
    weight[2] = c;
    weight[3] = c;
    weight[4] = d;
    weight[5] = d;
    weight[6] = d;
    weight[7] = d;
    weight[8] = d;
    weight[9] = d;
  }
//
//  Zienkiewicz #3.
//
  else if ( rule == 5 )
  {
    a =  1.0 / 6.0;
    b =  0.25;
    c =  0.5;
    d = -0.8;
    e =  0.45;

    xtab[0] = b;
    xtab[1] = c;
    xtab[2] = a;
    xtab[3] = a;
    xtab[4] = a;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = c;
    ytab[3] = a;
    ytab[4] = a;

    ztab[0] = b;
    ztab[1] = a;
    ztab[2] = a;
    ztab[3] = c;
    ztab[4] = a;

    weight[0] = d;
    weight[1] = e;
    weight[2] = e;
    weight[3] = e;
    weight[4] = e;
  }
//
//  Newton Cotes #3.
//  (This is actually formally a 20 point rule, but with 12 zero coefficients.)
//
  else if ( rule == 6 )
  {
    a = 1.0;
    b = 1.0 / 40.0;
    c = 1.0 /  3.0;
    d = 9.0 / 40.0;
    z = 0.0;

    xtab[0] = z;
    xtab[1] = a;
    xtab[2] = z;
    xtab[3] = z;
    xtab[4] = c;
    xtab[5] = c;
    xtab[6] = z;
    xtab[7] = c;

    ytab[0] = z;
    ytab[1] = z;
    ytab[2] = a;
    ytab[3] = z;
    ytab[4] = c;
    ytab[5] = z;
    ytab[6] = c;
    ytab[7] = c;

    ztab[0] = z;
    ztab[1] = z;
    ztab[2] = z;
    ztab[3] = a;
    ztab[4] = z;
    ztab[5] = c;
    ztab[6] = c;
    ztab[7] = c;

    weight[0] = b;
    weight[1] = b;
    weight[2] = b;
    weight[3] = b;
    weight[4] = d;
    weight[5] = d;
    weight[6] = d;
    weight[7] = d;
  }
//
//  Newton Cotes #4.
//
  else if ( rule == 7 )
  {
    a =   0.25;
    b =   0.50;
    c =   0.75;
    d =   1.00;
    e =  -5.0 / 420.0;
    f = -12.0 / 420.0;
    g =  16.0 / 420.0;
    h = 128.0 / 420.0;
    z =   0.0;

    xtab[0] = z; 
    xtab[1] = d; 
    xtab[2] = z; 
    xtab[3] = z; 
    xtab[4] = a; 
    xtab[5] = z; 
    xtab[6] = z; 
    xtab[7] = c; 
    xtab[8] = c; 
    xtab[9] = c;
    xtab[10] = z; 
    xtab[11] = a; 
    xtab[12] = z; 
    xtab[13] = z; 
    xtab[14] = a; 
    xtab[15] = z; 
    xtab[16] = b; 
    xtab[17] = z; 
    xtab[18] = z;            
    xtab[19] = b; 
    xtab[20] = b; 
    xtab[21] = z; 
    xtab[22] = a; 
    xtab[23] = b; 
    xtab[24] = a; 
    xtab[25] = a; 
    xtab[26] = b; 
    xtab[27] = z; 
    xtab[28] = b; 
    xtab[29] = z; 
    xtab[30] = a; 
    xtab[31] = a; 
    xtab[32] = z; 
    xtab[33] = a;
    xtab[34] = a;

    ytab[0] = z;
    ytab[1] = z;
    ytab[2] = d;
    ytab[3] = z;
    ytab[4] = z;
    ytab[5] = a;
    ytab[6] = z;
    ytab[7] = z;
    ytab[8] = a;
    ytab[9] = z;
    ytab[10] = c;
    ytab[11] = c;
    ytab[12] = c;
    ytab[13] = z;
    ytab[14] = z;
    ytab[15] = a;
    ytab[16] = z;
    ytab[17] = b;
    ytab[18] = z;
    ytab[19] = b;
    ytab[20] = z;
    ytab[21] = b;
    ytab[22] = a;
    ytab[23] = a;
    ytab[24] = b;
    ytab[25] = z;
    ytab[26] = z;
    ytab[27] = a;
    ytab[28] = a;
    ytab[29] = b;
    ytab[30] = b;
    ytab[31] = z;
    ytab[32] = a;
    ytab[33] = a;
    ytab[34] = a;

    ztab[0] = z;
    ztab[1] = z;
    ztab[2] = z;
    ztab[3] = d;
    ztab[4] = z;
    ztab[5] = z;
    ztab[6] = a;
    ztab[7] = z;
    ztab[8] = z;
    ztab[9] = a;
    ztab[10] = z;
    ztab[11] = z;
    ztab[12] = a;
    ztab[13] = c;
    ztab[14] = c;
    ztab[15] = c;
    ztab[16] = z;
    ztab[17] = z;
    ztab[18] = b;
    ztab[19] = z;
    ztab[20] = b;
    ztab[21] = b;
    ztab[22] = z;
    ztab[23] = z;
    ztab[24] = z;
    ztab[25] = a;
    ztab[26] = a;
    ztab[27] = a;
    ztab[28] = a;
    ztab[29] = a;
    ztab[30] = a;
    ztab[31] = b;
    ztab[32] = b;
    ztab[33] = b;
    ztab[34] = a;

    weight[0] = e;
    weight[1] = e;
    weight[2] = e;
    weight[3] = e;
    weight[4] = g;
    weight[5] = g;
    weight[6] = g;
    weight[7] = g;
    weight[8] = g;
    weight[9] = g;
    weight[10] = g;
    weight[11] = g;
    weight[12] = g;
    weight[13] = g;
    weight[14] = g;
    weight[15] = g;
    weight[16] = f;
    weight[17] = f;
    weight[18] = f;
    weight[19] = f;
    weight[20] = f;
    weight[21] = f;
    weight[22] = g;
    weight[23] = g;
    weight[24] = g;
    weight[25] = g;
    weight[26] = g;
    weight[27] = g;
    weight[28] = g;
    weight[29] = g;
    weight[30] = g;
    weight[31] = g;
    weight[32] = g;
    weight[33] = g;
    weight[34] = h;
  }
//
//  Keast Rule of order 11
//
  else if ( rule == 8 )
  {
    a =  0.25;
    b =  11.0 /    14.0;
    c =   1.0 /    14.0;
    d =  0.25 * ( 1.0 + sqrt ( 5.0 / 14.0 ) );
    e =  0.25 * ( 1.0 - sqrt ( 5.0 / 14.0 ) );
    f = -74.0 /  5625.0;
    g = 343.0 / 45000.0;
    h =  56.0 /  2250.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = c;
    xtab[5] = d;
    xtab[6] = d;
    xtab[7] = d;
    xtab[8] = e;
    xtab[9] = e;
    xtab[10] = e;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = c;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = e;
    ytab[8] = d;
    ytab[9] = d;
    ytab[10] = e;

    ztab[0] = a;
    ztab[1] = c;
    ztab[2] = c;
    ztab[3] = b;
    ztab[4] = c;
    ztab[5] = e;
    ztab[6] = d;
    ztab[7] = e;
    ztab[8] = d;
    ztab[9] = e;
    ztab[10] = d;

    weight[0] = f;
    weight[1] = g;
    weight[2] = g;
    weight[3] = g;
    weight[4] = g;
    weight[5] = h;
    weight[6] = h;
    weight[7] = h;
    weight[8] = h;
    weight[9] = h;
    weight[10] = h;
  }
  else
  {
    cerr << "\n";
    cerr << "TETRA_UNIT_SET - Fatal error!\n";
    cerr << "  Illegal value of RULE = " << rule << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int tetra_unit_size ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_UNIT_SIZE sizes quadrature weights and abscissas in the unit tetrahedron.
//
//  Integration region:
//
//      0 <= X,
//    and
//      0 <= Y,
//    and
//      0 <= Z, 
//    and
//      X + Y + Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hermann Engels,
//    Numerical Quadrature and Cubature,
//    Academic Press, 1980,
//    ISBN: 012238850X,
//    LC: QA299.3E5.
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//    ISBN: 0750663200,
//    LC: TA640.2.Z54
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
//     2, order 4, precision 1, Newton Cotes formula #1.
//     3, order 4, precision 2, Zienkiewicz #2.
//     4, order 10, precision 2, Newton Cotes formula #2
//     5, order 5, precision 3, Zienkiewicz #3.
//     6, order 8, precision 3, Newton Cotes formula #3.
//     7, order 35, precision 4, Newton Cotes formula #4.
//     8, order 11, precision 4, a Keast rule.
//
//    Output, int TETRA_UNIT_SET, the order of the rule.
//
{
  int order;
//
//  Newton Cotes #0.
//
  if ( rule == 1 )
  {
    order = 1;
  }
//
//  Newton Cotes #1.
//
  else if ( rule == 2 )
  {
    order = 4;
  }
//
//  Zienkiewicz #2.
//
  else if ( rule == 3 )
  {
    order = 4;
  }
//
//  Newton Cotes #2.
//
  else if ( rule == 4 )
  {
    order = 10;
  }
//
//  Zienkiewicz #3.
//
  else if ( rule == 5 )
  {
    order = 5;
  }
//
//  Newton Cotes #3.
//  (This is actually formally a 20 point rule, but with 12 zero coefficients//)
//
  else if ( rule == 6 )
  {
    order = 8;
  }
//
//  Newton Cotes #4.
//
  else if ( rule == 7 )
  {
    order = 35;
  }
//
//  Keast Rule of order 11
//
  else if ( rule == 8 )
  {
    order = 11;
  }
  else
  {
    cerr << "\n";
    cerr << "TETRA_UNIT_SIZE - Fatal error!\n";
    cerr << "  Illegal value of RULE = " << rule << "\n";
    exit ( 1 );
  }

  return order;
}
//****************************************************************************80

double tetra_unit_sum ( double func ( double x, double y, double z ), 
  int order, double xtab[], double ytab[], double ztab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_UNIT_SUM carries out a quadrature rule in the unit tetrahedron in 3D.
//
//  Integration region:
//
//      0 <= X,
//    and
//      0 <= Y,
//    and
//      0 <= Z, 
//    and
//      X + Y + Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], ZTAB[ORDER], the
//    abscissas.
//
//    Input, double WEIGHT[ORDER], the weights.
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;

  quad = 0.0;

  for ( i = 0; i < order; i++ )
  {
    quad = quad + weight[i] * func ( xtab[i], ytab[i], ztab[i] );
  }

  volume = tetra_unit_volume ( );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double tetra_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_UNIT_VOLUME returns the volume of the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X,
//      0 <= Y,
//      0 <= Z, 
//      X + Y + Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TETRA_UNIT_VOLUME, the volume.
//
{
  double volume;

  volume = 1.0 / 6.0;

  return volume;
}
//****************************************************************************80

double tetra_volume ( double x[4], double y[4], double z[4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRA_VOLUME computes the volume of a tetrahedron.
//
//  Integration region:
//
//    Points inside a tetrahedron whose four vertices are given.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X[4], Y[4], Z[4], the vertices.
//
//    Output, double TETRA_VOLUME, the volume of the tetrahedron.
//
{
  double volume;

  volume = parallelipiped_volume_3d ( x, y, z );

  volume = volume * tetra_unit_volume ( );

  return volume;
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
//    03 October 2003
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

double torus_1 ( double func ( double x, double y, double z ), double r1, 
  double r2, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_1 approximates an integral on the surface of a torus in 3D.
//
//  Integration region:
//
//    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2 * R2.
//
//  Discussion:
//
//    An (N+1)*(N+2) point N-th degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Input, int N, defines the degree of the formula
//    used to approximate the integral.
//
//    Output, double TORUS_1, the approximate integral of the function.
//
{
  double angle;
  double ct1;
  int i;
  int j;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double st1;
  double u;
  double volume;
  double w;
  double x;
  double y;
  double z;

  w = 1.0 / ( r1 * ( double ) ( ( n + 1 ) * ( n + 2 ) ) );
  quad = 0.0;

  for ( i = 0; i < n + 1; i++ )
  {
    angle = 2.0 * pi * ( double ) ( i + 1 ) / ( double ) ( n + 1 );
    ct1 = cos ( angle );
    st1 = sin ( angle );

    for ( j = 0; j < n + 2; j++ )
    {
      angle = 2.0 * pi * ( double ) ( j + 1 ) / ( double ) ( n + 2 );
      u = r1 + r2 * cos ( angle );
      x = u * ct1;
      y = u * st1;
      z = r2 * sin ( angle );

      quad = quad + w * u * func ( x, y, z );
    }
  }
  volume = torus_area_3d ( r1, r2 );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double torus_14s ( double func ( double x, double y, double z ), double r1, 
  double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_14S approximates an integral inside a torus in 3D.
//
//  Integration region:
//
//    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
//
//  Discussion:
//
//    A 960 point 14-th degree formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the
//    user function which is to be integrated.
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Output, double TORUS_14S, the approximate integral of the function.
//
{
  double angle;
  double ct;
  double cth;
  int i;
  int j;
  int n;
  int order = 4;
  double pi = 3.141592653589793;
  double quad;
  double r[4] = {
    0.263499230, 0.574464514, 0.818529487, 0.964659606 };
  double result;
  double st;
  double sth;
  double u;
  double volume;
  double weight[4] = {
    0.086963711, 0.163036289, 0.163036289, 0.086963711 };
  double x;
  double y;
  double z;

  quad = 0.0;

  for ( n = 1; n <= 15; n++ )
  {
    angle = 2.0 * pi * ( double ) ( n ) / 15.0;
    cth = cos ( angle );
    sth = sin ( angle );

    for ( i = 1; i <= 16; i++ )
    {
      angle = 2.0 * pi * ( double ) ( i ) / 16.0;
      ct = cos ( angle );
      st = sin ( angle );

      for ( j = 0; j < order; j++ )
      {
        u = r1 + r[j] * ct * r2;
        x = u * cth;
        y = u * sth;
        z = r[j] * st * r2;
        quad = quad + u * weight[j] * func ( x, y, z ) / ( 120.0 * r1 );
      }
    }
  }

  volume = torus_volume_3d ( r1, r2 );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double torus_5s2 ( double func ( double x, double y, double z ), double r1, 
  double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_5S2 approximates an integral inside a torus in 3D.
//
//  Integration region:
//
//    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
//
//  Discussion:
//
//    A 24 point, 5-th degree formula is used, Stroud number TOR3-S2:5-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of 
//    the user supplied function to be integrated.
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Output, double TORUS_5S2, the approximate integral of the function.
//
{
  double angle;
  double cs;
  int i;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double sn;
  double u1;
  double u2;
  double u3;
  double volume;
  double w;
  double x;
  double y;
  double z;

  w = 1.0 / 24.0;

  quad = 0.0;

  u1 = sqrt ( r1 * r1 + 0.5 * r2 * r2 );
  u2 = sqrt ( r1 * r1 + sqrt ( 2.0 ) * r1 * r2 + r2 * r2 );
  u3 = sqrt ( r1 * r1 - sqrt ( 2.0 ) * r1 * r2 + r2 * r2 );

  for ( i = 1; i <= 6; i++ )
  {
    angle = 2.0 * pi * ( double ) ( i ) / 6.0;
    cs = cos ( angle );
    sn = sin ( angle );

    x = u1 * cs;
    y = u1 * sn;
    z = r2 / sqrt ( 2.0 );
    quad = quad + w * func ( x, y, z );

    x = u1 * cs;
    y = u1 * sn;
    z = -r2 / sqrt ( 2.0 );
    quad = quad + w * func ( x, y, z );

    x = u2 * cs;
    y = u2 * sn;
    z = 0.0;
    quad = quad + w * func ( x, y, z );

    x = u3 * cs;
    y = u3 * sn;
    z = 0.0;
    quad = quad + w * func ( x, y, z );
  }

  volume = torus_volume_3d ( r1, r2 );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double torus_6s2 ( double func ( double x, double y, double z ), double r1, 
  double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_6S2 approximates an integral inside a torus in 3D.
//
//  Integration region:
//
//    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z <= R2 * R2.
//
//  Discussion:
//
//    An 84 point 6-th degree formula is used, Stroud number TOR3-S2:6-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the
//    user defined function to be integrated.
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Output, double TORUS_6S2, the approximate integral of the function.
//
{
  double cth;
  int i;
  int j;
  int k;
  int n;
  int order = 2;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double s[2] = { 0.322914992, 0.644171310 };
  double sth;
  double u;
  double v;
  double volume;
  double w;
  double weight[2] = { 0.387077796, 0.165609800 };
  double x;
  double y;
  double z;

  w = 1.0 / ( 7.0 * r1 * pi );

  quad = 0.0;

  for ( n = 1; n <= 7; n++ )
  {
    u = 0.5 * sqrt ( 3.0 ) * r2;
    cth = cos ( 2.0 * pi * ( double ) ( n ) / 7.0 );
    sth = sin ( 2.0 * pi * ( double ) ( n ) / 7.0 );

    for ( i = 1; i <= 2; i++ )
    {
      u = -u;

      x = ( r1 + u ) * cth;
      y = ( r1 + u ) * sth;
      z = 0.0;
      quad = quad + 0.232710567 * w * ( r1 + u ) * func ( x, y, z );

      x = r1 * cth;
      y = r1 * sth;
      z = u;
      quad = quad + 0.232710567 * w * r1 * func ( x, y, z );

    }
    for ( k = 0; k < order; k++ )
    {
      u = s[k] * r2;
      v = u;

      for ( i = 1; i <= 2; i++ )
      {
        u = -u;

        for ( j = 1; j <= 2; j++ )
        {
          v = -v;

          x = ( r1 + u ) * cth;
          y = ( r1 + u ) * sth;
          z = v;
          quad = quad + weight[k] * w * ( r1 + u ) * func ( x, y, z );
        }
      }
    }
  }
  volume = torus_volume_3d ( r1, r2 );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double torus_area_3d ( double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_AREA_3D returns the area of a torus in 3D.
//
//  Integration region:
//
//    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z*Z = R2*R2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Output, double TORUS_AREA_3D, the area of the torus.
//
{
  double area;
  double pi = 3.141592653589793;

  area = 4.0 * pi * pi * r1 * r2;

  return area;
}
//****************************************************************************80

double torus_square_14c ( double func ( double x, double y, double z ), 
  double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_SQUARE_14C approximates an integral in a "square" torus in 3D.
//
//  Discussion:
//
//    A 14-th degree 960 point formula is used.
//
//  Integration region:
//
//      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
//    and
//      -R2 <= Z <= R2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double R1, R2, the radii that define the torus.
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  double angle;
  double cth;
  int i;
  int j;
  int n;
  int order = 8;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double *rtab;
  double sth;
  double u;
  double volume;
  double w;
  double *weight;
  double x;
  double y;
  double z;

  rtab = new double[order];
  weight = new double[order];

  legendre_set ( order, rtab, weight );

  w = 1.0 / ( 60.0 * r1 );
  quad = 0.0;

  for ( n = 1; n <= 15; n++ )
  {
    angle = 2.0 * pi * ( double ) ( n ) / 15.0;
    cth = cos ( angle );
    sth = sin ( angle );

    for ( i = 0; i < order; i++ )
    {
      u = r1 + rtab[i] * r2;
      x = u * cth;
      y = u * sth;

      for ( j = 0; j < order; j++ )
      {
        z = rtab[j] * r2;
        quad = quad + u * w * weight[i] * weight[j] * func ( x, y, z );
      }
    }
  }

  volume = torus_square_volume_3d ( r1, r2 );
  result = quad * volume;

  delete [] rtab;
  delete [] weight;

  return result;
}
//****************************************************************************80

double torus_square_5c2 ( double func ( double x, double y, double z ), 
  double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_SQUARE_5C2 approximates an integral in a "square" torus in 3D.
//
//  Integration region:
//
//      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
//    and
//      -R2 <= Z <= R2.
//
//  Discussion:
//
//    A 24 point 5-th degree formula is used, Stroud number TOR3-C2:5-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double R1, the primary radius of the torus.
//
//    Input, double R2, one-half the length of a side of the
//    square cross-section.
//
//    Output, double TORUS_SQUARE_5C2, the approximate integral of the function.
//
{
  double b1 = 5.0 / 108.0;
  double b2 = 4.0 / 108.0;
  double cs;
  int i;
  double pi = 3.141592653589793;
  double quad;
  double result;
  double sn;
  double u1;
  double u2;
  double u3;
  double v;
  double volume;
  double x;
  double y;
  double z;

  quad = 0.0;

  u1 = sqrt ( r1 * r1 + r2 * r2 );

  v = r2 * sqrt ( 0.6 );

  u2 = sqrt ( r1 * r1 - sqrt ( 3.0 ) * r1 * r2 + r2 * r2 );

  u3 = sqrt ( r1 * r1 + sqrt ( 3.0 ) * r1 * r2 + r2 * r2 );

  for ( i = 1; i <= 6; i++ )
  {
    cs = cos ( ( double ) ( i ) * pi / 3.0 );
    sn = sin ( ( double ) ( i ) * pi / 3.0 );

    x = u1 * cs;
    y = u1 * sn;
    z = v;
    quad = quad + b1 * func ( x, y, z );

    z = -v;
    quad = quad + b1 * func ( x, y, z );

    x = u2 * cs;
    y = u2 * sn;
    z = 0.0;
    quad = quad + b2 * func ( x, y, z );

    x = u3 * cs;
    y = u3 * sn;
    z = 0.0;
    quad = quad + b2 * func ( x, y, z );
  }

  volume = torus_square_volume_3d ( r1, r2 );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double torus_square_area_3d ( double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_SQUARE_AREA_3D returns the area of a square torus in 3D.
//
//  Integration region:
//
//      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
//    and
//      -R2 <= Z <= R2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Output, double TORUS_SQUARE_AREA_3D, the area of the torus.
//
{
  double area;
  double pi = 3.141592653589793;

  area = 16.0 * pi * r1 * r2;

  return area;
}
//****************************************************************************80

double torus_square_volume_3d ( double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_SQUARE_VOLUME_3D returns the volume of a square torus in 3D.
//
//  Integration region:
//
//      R1 - R2 <= SQRT ( X*X + Y*Y ) <= R1 + R2,
//    and
//      -R2 <= Z <= R2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Output, double TORUS_SQUARE_VOLUME_3D, the volume of the torus.
//
{
  double pi = 3.141592653589793;
  double volume;

  volume = 8.0 * pi * r1 * r2 * r2;

  return volume;
}
//****************************************************************************80

double torus_volume_3d ( double r1, double r2 )

//****************************************************************************80
//
//  Purpose:
//
//    TORUS_VOLUME_3D returns the volume of a torus in 3D.
//
//  Integration region:
//
//    ( SQRT ( X*X + Y*Y ) - R1 )^2 + Z * Z = R2 * R2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R1, R2, the two radii that define the torus.
//
//    Output, double TORUS_VOLUME_3D, the volume of the torus.
//
{
  double pi = 3.141592653589793;
  double volume;

  volume = 2.0 * pi * pi * r1 * r2 * r2;

  return volume;
}
//****************************************************************************80

void triangle_rule_adjust ( double xval[3], double yval[3], int order, 
  double xtab[], double ytab[], double weight[], double xtab2[], double ytab2[], 
  double weight2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_RULE_ADJUST adjusts a unit quadrature rule to an arbitrary triangle.
//
//  Integration region:
//
//      (X,Y) = ALPHA * (X1,Y1) + BETA * (X2,Y2) + ( 1 - ALPHA - BETA ) * (X3,Y3)
//    and
//      0 <= ALPHA <= 1 - BETA
//    and
//      0 <= BETA <= 1 - ALPHA
//
//  Discussion:
//
//    This routine accepts as input abscissas and weights appropriate for
//    quadrature in the unit triangle, and returns abscissas and weights
//    appropriate for quadrature in a given triangle.
//
//    Once this routine has been called, an integral over the given triangle
//    can be approximated as:
//
//      QUAD = sum ( 1 <= I <= ORDER ) WTAB2(I) * FUNC ( XTAB2(I), YTAB2(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XVAL[3], YVAL[3], the coordinates of the nodes.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas for
//    the unit triangle.
//
//    Input, double WEIGHT[ORDER], the weights for the unit triangle.
//
//    Output, double XTAB2[ORDER], YTAB2[ORDER], the adjusted
//    abscissas.
//
//    Output, double WEIGHT2[ORDER], the adjusted weights.
//
{
  int i;
  double volume;

  volume = triangle_volume ( xval, yval );

  for ( i = 0; i < order; i++ )
  {
    xtab2[i] =         xtab[i]             * xval[0] 
             +                   ytab[i]   * xval[1] 
             + ( 1.0 - xtab[i] - ytab[i] ) * xval[2];

    ytab2[i] =         xtab[i]             * yval[0] 
             +                   ytab[i]   * yval[1] 
             + ( 1.0 - xtab[i] - ytab[i] ) * yval[2];

    weight2[i] = weight[i] * 2.0 * volume;
  }

  return;
}
//****************************************************************************80

double triangle_sub ( double func ( double x, double y ), double xval[3], 
  double yval[], int nsub, int order, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SUB carries out quadrature over subdivisions of a triangular region.
//
//  Integration region:
//
//      (X,Y) =       ALPHA          * ( XVAL[0], YVAL[0] )
//            +               BETA   * ( XVAL[1], YVAL[1] )
//            + ( 1 - ALPHA - BETA ) * ( XVAL[2], YVAL[2] )
//    and
//      0 <= ALPHA <= 1 - BETA
//    and
//      0 <= BETA <= 1 - ALPHA
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the user supplied
//    function to be integrated.
//
//    Input, double XVAL[3], YVAL[3], the coordinates of the triangle vertices.
//
//    Input, int NSUB, the number of subdivisions of each side of the
//    input triangle to be made.  NSUB = 1 means no subdivisions are made.
//    NSUB = 3 means that each side of the triangle is subdivided into
//    three portions, and that the original triangle is subdivided into
//    NSUB * NSUB triangles.  NSUB must be at least 1.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double TRIANGLE_SUB, the approximate integral of the function.
//
{
  int i;
  int j;
  int k;
  double quad;
  double result;
  double temp1;
  double temp2;
  double volume;
  double x;
  double x1;
  double x2;
  double x3;
  double y;
  double y1;
  double y2;
  double y3;
//
//  Initialize RESULT, the approximate integral.
//
  result = 0.0;
//
//  NSUB must be positive.
//
  if ( nsub <= 0 )
  {
    return result;
  }
//
//  Initialize QUAD, the quadrature sum.
//
  quad = 0.0;
//
//  The sub-triangles can be grouped into NSUB strips.
//
  for ( i = 1; i <= nsub; i++ )
  {
    temp1 = 0.0;
    temp2 = ( double ) ( i ) / ( double ) ( nsub );

    x2 = xval[1] + temp1 * ( xval[2] - xval[1] ) 
                 + temp2 * ( xval[0] - xval[1] );

    y2 = yval[1] + temp1 * ( yval[2] - yval[1] ) 
                 + temp2 * ( yval[0] - yval[1] );

    temp1 = 0.0;
    temp2 = ( double ) ( i - 1 ) / ( double ) ( nsub );

    x3 = xval[1] + temp1 * ( xval[2] - xval[1] ) 
                 + temp2 * ( xval[0] - xval[1] );

    y3 = yval[1] + temp1 * ( yval[2] - yval[1] ) 
                 + temp2 * ( yval[0] - yval[1] );
//
//  There are 2*I-1 triangles in strip number I.
//  The next triangle in the strip shares two nodes with the previous one.
//  Compute its corners, (X1,Y1), (X2,Y2), (X3,Y3).
//
    for ( j = 1; j <= 2 * i - 1; j++ )
    {
      x1 = x2;
      y1 = y2;
      x2 = x3;
      y2 = y3;
      temp1 = ( double ) ( ( ( j + 1 ) / 2 ) ) / ( double ) ( nsub );
      temp2 = ( double ) ( ( i - 1 - ( j / 2 ) ) ) / ( double ) ( nsub );

      x3 = xval[1] + temp1 * ( xval[2] - xval[1] ) 
                   + temp2 * ( xval[0] - xval[1] );

      y3 = yval[1] + temp1 * ( yval[2] - yval[1] )
                   + temp2 * ( yval[0] - yval[1] );
//
//  Now integrate over the triangle, mapping the points ( XTAB(K), YTAB(K) )
//  into the triangle.
//
      for ( k = 0; k < order; k++ )
      {
        x = x2 + xtab[k] * ( x3 - x2 ) + ytab[k] * ( x1 - x2 );
        y = y2 + xtab[k] * ( y3 - y2 ) + ytab[k] * ( y1 - y2 );
        quad = quad + weight[k] * func ( x, y );
       }
    }
  }

  volume = triangle_volume ( xval, yval ) / ( double ) ( nsub * nsub );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double triangle_sum ( double func ( double x, double y ), 
  double xval[3], double yval[3], int order, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SUM carries out a unit quadrature rule in an arbitrary triangle.
//
//  Integration region:
//
//      (X,Y) =       ALPHA          * (X1,Y1) 
//            +               BETA   * (X2,Y2) 
//            + ( 1 - ALPHA - BETA ) * (X3,Y3)
//    and
//      0 <= ALPHA <= 1 - BETA
//    and
//      0 <= BETA <= 1 - ALPHA
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y, double z ), the name of the 
//    user supplied function to be integrated.
//
//    Input, double XVAL[3], YVAL[3], the coordinates of the nodes.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;
  double x;
  double y;

  quad = 0.0;

  for ( i = 0; i < order; i++ )
  {
    x =         xtab[i]             * xval[0] 
      +                   ytab[i]   * xval[1]
      + ( 1.0 - xtab[i] - ytab[i] ) * xval[2];

    y =         xtab[i]             * yval[0] 
      +                   ytab[i]   * yval[1]
      + ( 1.0 - xtab[i] - ytab[i] ) * yval[2];

    quad = quad + weight[i] * func ( x, y );
  }

  volume = triangle_volume ( xval, yval );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double triangle_sum_adjusted ( double func ( double x, double y ), 
  int order, double xtab[], double ytab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SUM_ADJUSTED carries out an adjusted quadrature rule in a triangle.
//
//  Integration region:
//
//      (X,Y) =       ALPHA          * (X1,Y1) 
//                          + BETA   * (X2,Y2) 
//            + ( 1 - ALPHA - BETA ) * (X3,Y3)
//    and
//      0 <= ALPHA <= 1 - BETA
//    and
//      0 <= BETA <= 1 - ALPHA
//
//  Discussion:
//
//    It is assumed that a quadrature rule approprate for the unit triangle
//    was generated, and then adjusted to a particular triangle by calling
//    TRIANGLE_RULE_ADJUST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the 
//    user supplied function to be integrated.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double TRIANGLE_SUM_ADJUSTED, the approximate integral 
//    of the function.
//
{
  int i;
  double result;

  result = 0.0;

  for ( i = 0; i < order; i++ )
  {
    result = result + weight[i] * func ( xtab[i], ytab[i] );
  }

  return result;
}
//****************************************************************************80

void triangle_unit_product_set ( int rule, int order, double xtab[], 
  double ytab[], double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_PRODUCT_SET sets a product rule on the unit triangle.
//
//  Discussion:
//
//    For a given order of accuracy, a product rule on a triangle usually
//    uses more points than necessary.  That is, there is usually a rule
//    of the same order that uses fewer points.
//
//    However, one advantage of product rules is that a rule of any
//    desired order can be generated automatically.
//   
//    The integration region is:
//
//      0 <= X,
//    and
//      0 <= Y, 
//    and
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the order of the 1D rule.
//
//    Input, int ORDER, the order of the rule.
//
//    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights of the rule.
//
{
  double a;
  double b;
  double c;
  double d;
  int i;
  int j;
  int k;
  int order0;
  int order1;
  double *weight0;
  double *weight1;
  double *xtab0;
  double *xtab1;

  weight0 = new double[rule];
  weight1 = new double[rule];
  xtab0 = new double[rule];
  xtab1 = new double[rule];

  a = -1.0;
  b = +1.0;
  c =  0.0;
  d = +1.0;

  order0 = rule;
  legendre_set ( order0, xtab0, weight0 );
  rule_adjust ( a, b, c, d, order0, xtab0, weight0 );

  order1 = rule;
  legendre_set_x1 ( order1, xtab1, weight1 );
  rule_adjust ( a, b, c, d, order1, xtab1, weight1 );

  k = 0;
  for ( j = 0; j < order1; j++ )
  {
    for ( i = 0; i < order0; i++ )
    {
      xtab[k] = 1.0 - xtab1[j];
      ytab[k] = xtab0[i] * xtab1[j];
      weight[k] = weight0[i] * weight1[j];
      k = k + 1;
    }
  }

  delete [] weight0;
  delete [] weight1;
  delete [] xtab0;
  delete [] xtab1;

  return;
}
//****************************************************************************80

int triangle_unit_product_size ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_PRODUCT_SIZE sizes a product rule on the unit triangle.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X,
//    and
//      0 <= Y, 
//    and
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, int RULE, the order of the 1D rule.
//
//    Input, int TRIANGLE_UNIT_PRODUCT_SIZE, the order of the rule.
//
{
  int order;

  order = rule * rule;

  return order;
}
//****************************************************************************80

void triangle_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_SET sets a quadrature rule in the unit triangle.
//
//  Discussion:
//
//    The user is responsible for determining the value of ORDER,
//    and appropriately dimensioning the arrays XTAB, YTAB and
//    WEIGHT so that they can accommodate the data.
//
//    The value of ORDER for each rule can be found by invoking
//    the function TRIANGLE_RULE_SIZE.
//
//  Integration region:
//
//      0 <= X,
//    and
//      0 <= Y, 
//    and
//      X + Y <= 1.
//
//  Graph:
//
//      ^
//    1 | *
//      | |\
//    Y | | \
//      | |  \
//    0 | *---*
//      +------->
//        0 X 1
//
//   The rules are accessed by an index number, RULE.  The indices,
//   and the descriptions of the corresponding rules, are:
//
//     1, ORDER =  1, precision 1, Zienkiewicz #1.
//     2, ORDER =  2, precision 1, (the "vertex rule").
//     3, ORDER =  3, precision 2, Strang and Fix formula #1.
//     4, ORDER =  3, precision 2, Strang and Fix formula #2,
//                                 Zienkiewicz #2.
//     5, ORDER =  4, precision 3, Strang and Fix formula #3,
//                                 Zienkiewicz #3.
//     6, ORDER =  6, precision 3, Strang and Fix formula #4.
//     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
//     8, ORDER =  6, precision 4, Strang and Fix formula #5.
//     9, ORDER =  7, precision 4, Strang and Fix formula #6.
//    10, ORDER =  7, precision 5, Strang and Fix formula #7,
//                                 Stroud formula T2:5-1, 
//                                 Zienkiewicz #4, 
//                                 Schwarz Table 2.2.
//    11, ORDER =  9, precision 6, Strang and Fix formula #8.
//    12, ORDER = 12, precision 6, Strang and Fix formula #9.
//    13, ORDER = 13, precision 7, Strang and Fix formula #10.
//        Note that there is a typographical error in Strang and Fix
//        which lists the value of the XSI(3) component of the
//        last generator point as 0.4869... when it should be 0.04869...
//    14, ORDER =  7, precision 3.
//    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
//    16, ORDER = 64, precision 15, triangular product Gauss rule.
//    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
//    18, ORDER = 19, precision 9, from TRIEX, ACM TOMS #612.
//    19, ORDER = 28, precision 11, from TRIEX, ACM TOMS #612.
//    20, ORDER = 37, precision 13, from ACM TOMS #706.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jarle Berntsen, Terje Espelid,
//    Algorithm 706,
//    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
//    ACM Transactions on Mathematical Software,
//    Volume 18, Number 3, September 1992, pages 329-342.
//
//    Elise deDoncker, Ian Robinson,
//    Algorithm 612:
//    Integration over a Triangle Using Nonlinear Extrapolation,
//    ACM Transactions on Mathematical Software,
//    Volume 10, Number 1, March 1984, pages 17-22.
//
//    Dirk Laurie,
//    Algorithm 584,
//    CUBTRI, Automatic Cubature Over a Triangle,
//    ACM Transactions on Mathematical Software,
//    Volume 8, Number 2, 1982, pages 210-218.
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//    Hans Rudolf Schwarz,
//    Finite Element Methods,
//    Academic Press, 1988,
//    ISBN: 0126330107,
//    LC: TA347.F5.S3313.
//
//    Gilbert Strang, George Fix,
//    An Analysis of the Finite Element Method,
//    Cambridge, 1973,
//    ISBN: 096140888X,
//    LC: TA335.S77.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//    ISBN: 0750663200,
//    LC: TA640.2.Z54
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int ORDER, the order of the rule.
//
//    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights of the rule.
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  double g;
  double h;
  int i;
  int j;
  int k;
  int order2;
  double p;
  double q;
  double r;
  double s;
  double t;
  double u;
  double v;
  double w;
  double w1;
  double w2;
  double w3;
  double w4;
  double w5;
  double w6;
  double w7;
  double w8;
  double w9;
  double weight1[8];
  double weight2[8];
  double wx;
  double x;
  double xtab1[8];
  double xtab2[8];
  double y;
  double z;

//
//  1 point, precision 1.
//
  if ( rule == 1 )
  {
    xtab[0]   = 0.33333333333333333333;

    ytab[0]   = 0.33333333333333333333;

    weight[0] = 1.00000000000000000000;
  }
//
//  3 points, precision 1, the "vertex rule".
//
  else if ( rule == 2 )
  {
    xtab[0] =   1.00000000000000000000;
    xtab[1] =   0.00000000000000000000;
    xtab[2] =   0.00000000000000000000;

    ytab[0] =   0.00000000000000000000;
    ytab[1] =   1.00000000000000000000;
    ytab[2] =   0.00000000000000000000;

    weight[0] = 0.33333333333333333333;
    weight[1] = 0.33333333333333333333;
    weight[2] = 0.33333333333333333333;
  }
//
//  3 points, precision 2, Strang and Fix formula #1.
//
  else if ( rule == 3 )
  {
    xtab[0]   = 0.66666666666666666667;
    xtab[1]   = 0.16666666666666666667;
    xtab[2]   = 0.16666666666666666667;

    ytab[0]   = 0.16666666666666666667;
    ytab[1]   = 0.66666666666666666667;
    ytab[2]   = 0.16666666666666666667;

    weight[0] = 0.33333333333333333333;
    weight[1] = 0.33333333333333333333;
    weight[2] = 0.33333333333333333333;
  }
//
//  3 points, precision 2, Strang and Fix formula #2.
//
  else if ( rule == 4 )
  {
    xtab[0]   = 0.50000000000000000000;
    xtab[1]   = 0.50000000000000000000;
    xtab[2]   = 0.00000000000000000000;

    ytab[0]   = 0.00000000000000000000;
    ytab[1]   = 0.50000000000000000000;
    ytab[2]   = 0.50000000000000000000;

    weight[0] = 0.33333333333333333333;
    weight[1] = 0.33333333333333333333;
    weight[2] = 0.33333333333333333333;
  }
//
//  4 points, precision 3, Strang and Fix formula #3.
//
  else if ( rule == 5 )
  {
    a =   6.0 / 30.0;
    b =  10.0 / 30.0;
    c =  18.0 / 30.0;

    d =  25.0 / 48.0;
    e = -27.0 / 48.0;

    xtab[0] = b;
    xtab[1] = c;
    xtab[2] = a;
    xtab[3] = a;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = c;
    ytab[3] = a;

    weight[0] = e;
    weight[1] = d;
    weight[2] = d;
    weight[3] = d;
  }
//
//  6 points, precision 3, Strang and Fix formula #4.
//
  else if ( rule == 6 )
  {
    a = 0.659027622374092;
    b = 0.231933368553031;
    c = 0.109039009072877;

    xtab[0] = a;
    xtab[1] = a;
    xtab[2] = b;
    xtab[3] = b;
    xtab[4] = c;
    xtab[5] = c;

    ytab[0] = b;
    ytab[1] = c;
    ytab[2] = a;
    ytab[3] = c;
    ytab[4] = a;
    ytab[5] = b;

    weight[0] = 0.16666666666666666667;
    weight[1] = 0.16666666666666666667;
    weight[2] = 0.16666666666666666667;
    weight[3] = 0.16666666666666666667;
    weight[4] = 0.16666666666666666667;
    weight[5] = 0.16666666666666666667;
  }
//
//  6 points, precision 3, Stroud T2:3-1.
//
  else if ( rule == 7 )
  {
    a = 0.0;
    b = 0.5;
    c = 2.0 /  3.0;
    d = 1.0 /  6.0;
    v = 1.0 / 30.0;
    w = 3.0 / 10.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = d;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = c;
    ytab[5] = d;

    weight[0] = v;
    weight[1] = v;
    weight[2] = v;
    weight[3] = w;
    weight[4] = w;
    weight[5] = w;
  }
//
//  6 points, precision 4, Strang and Fix, formula #5.
//
  else if ( rule == 8 )
  {
    a = 0.816847572980459;
    b = 0.091576213509771;
    c = 0.108103018168070;
    d = 0.445948490915965;
    v = 0.109951743655322;
    w = 0.223381589678011;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = b;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = d;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = c;
    ytab[5] = d;

    weight[0] = v;
    weight[1] = v;
    weight[2] = v;
    weight[3] = w;
    weight[4] = w;
    weight[5] = w;
  }
//
//  7 points, precision 4, Strang and Fix formula #6.
//
  else if ( rule == 9 )
  {
    a = 1.0 / 3.0;
    c = 0.736712498968435;
    d = 0.237932366472434;
    e = 0.025355134551932;
    v = 0.375000000000000;
    w = 0.104166666666667;

    xtab[0] = a;
    xtab[1] = c;
    xtab[2] = c;
    xtab[3] = d;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;

    ytab[0] = a;
    ytab[1] = d;
    ytab[2] = e;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = c;
    ytab[6] = d;

    weight[0] = v;
    weight[1] = w;
    weight[2] = w;
    weight[3] = w;
    weight[4] = w;
    weight[5] = w;
    weight[6] = w;
  }
//
//  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
//
  else if ( rule == 10 )
  {
    a = 1.0 / 3.0;
    b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
    d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
    u = 0.225;
    v = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
    w = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;

    weight[0] = u;
    weight[1] = v;
    weight[2] = v;
    weight[3] = v;
    weight[4] = w;
    weight[5] = w;
    weight[6] = w;
  }
//
//  9 points, precision 6, Strang and Fix formula #8.
//
  else if ( rule == 11 )
  {
    a = 0.124949503233232;
    b = 0.437525248383384;
    c = 0.797112651860071;
    d = 0.165409927389841;
    e = 0.037477420750088;

    u = 0.205950504760887;
    v = 0.063691414286223;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = b;
    xtab[3] = c;
    xtab[4] = c;
    xtab[5] = d;
    xtab[6] = d;
    xtab[7] = e;
    xtab[8] = e;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = e;
    ytab[5] = c;
    ytab[6] = e;
    ytab[7] = c;
    ytab[8] = d;

    weight[0] = u;
    weight[1] = u;
    weight[2] = u;
    weight[3] = v;
    weight[4] = v;
    weight[5] = v;
    weight[6] = v;
    weight[7] = v;
    weight[8] = v;
  }
//
//  12 points, precision 6, Strang and Fix, formula #9.
//
  else if ( rule == 12 )
  {
    a = 0.873821971016996;
    b = 0.063089014491502;
    c = 0.501426509658179;
    d = 0.249286745170910;
    e = 0.636502499121399;
    f = 0.310352451033785;
    g = 0.053145049844816;

    u = 0.050844906370207;
    v = 0.116786275726379;
    w = 0.082851075618374;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = b;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = d;
    xtab[6] = e;
    xtab[7] = e;
    xtab[8] = f;
    xtab[9] = f;
    xtab[10] = g;
    xtab[11] = g;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = c;
    ytab[5] = d;
    ytab[6] = f;
    ytab[7] = g;
    ytab[8] = e;
    ytab[9] = g;
    ytab[10] = e;
    ytab[11] = f;

    weight[0] = u;
    weight[1] = u;
    weight[2] = u;
    weight[3] = v;
    weight[4] = v;
    weight[5] = v;
    weight[6] = w;
    weight[7] = w;
    weight[8] = w;
    weight[9] = w;
    weight[10] = w;
    weight[11] = w;
  }
//
//  13 points, precision 7, Strang and Fix, formula #10.
//
//  Note that there is a typographical error in Strang and Fix
//  which lists the value of the XSI[2] component of the
//  last generator point as 0.4869... when it should be 0.04869...
//
  else if ( rule == 13 )
  {
    h = 1.0 / 3.0;
    a = 0.479308067841923;
    b = 0.260345966079038;
    c = 0.869739794195568;
    d = 0.065130102902216;
    e = 0.638444188569809;
    f = 0.312865496004875;
    g = 0.048690315425316;

    w = -0.149570044467670;
    t =  0.175615257433204;
    u =  0.053347235608839;
    v =  0.077113760890257;

    xtab[0] = h;
    xtab[1] = a;
    xtab[2] = b;
    xtab[3] = b;
    xtab[4] = c;
    xtab[5] = d;
    xtab[6] = d;
    xtab[7] = e;
    xtab[8] = e;
    xtab[9] = f;
    xtab[10] = f;
    xtab[11] = g;
    xtab[12] = g;

    ytab[0] = h;
    ytab[1] = b;
    ytab[2] = a;
    ytab[3] = b;
    ytab[4] = d;
    ytab[5] = c;
    ytab[6] = d;
    ytab[7] = f;
    ytab[8] = g;
    ytab[9] = e;
    ytab[10] = g;
    ytab[11] = e;
    ytab[12] = f;

    weight[0] = w;
    weight[1] = t;
    weight[2] = t;
    weight[3] = t;
    weight[4] = u;
    weight[5] = u;
    weight[6] = u;
    weight[7] = v;
    weight[8] = v;
    weight[9] = v;
    weight[10] = v;
    weight[11] = v;
    weight[12] = v;
  }
//
//  7 points, precision 3.
//
  else if ( rule == 14 )
  {
    a = 1.0 / 3.0;
    b = 1.0;
    c = 0.5;
    z = 0.0;

    u = 27.0 / 60.0;
    v =  3.0 / 60.0;
    w =  8.0 / 60.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = z;
    xtab[3] = z;
    xtab[4] = z;
    xtab[5] = c;
    xtab[6] = c;

    ytab[0] = a;
    ytab[1] = z;
    ytab[2] = b;
    ytab[3] = z;
    ytab[4] = c;
    ytab[5] = z;
    ytab[6] = c;

    weight[0] = u;
    weight[1] = v;
    weight[2] = v;
    weight[3] = v;
    weight[4] = w;
    weight[5] = w;
    weight[6] = w;
  }
//
//  16 points, precision 5, Stroud T2:7-1.
//
  else if ( rule == 15 )
  {
//
//  Legendre rule of order 4.
//
    order2 = 4;

    xtab[0] = - 0.861136311594052575223946488893;
    xtab[1] = - 0.339981043584856264802665759103;
    xtab[2] =   0.339981043584856264802665759103;
    xtab[3] =   0.861136311594052575223946488893;

    weight1[0] = 0.347854845137453857373063949222;
    weight1[1] = 0.652145154862546142626936050778;
    weight1[2] = 0.652145154862546142626936050778;
    weight1[3] = 0.347854845137453857373063949222;

    for ( i = 0; i < order2; i++ )
    {
      xtab1[i] = 0.5 * ( xtab1[i] + 1.0 );
    }

    weight2[0] = 0.1355069134;
    weight2[1] = 0.2034645680;
    weight2[2] = 0.1298475476;
    weight2[3] = 0.0311809709;

    xtab2[0] = 0.0571041961;
    xtab2[1] = 0.2768430136;
    xtab2[2] = 0.5835904324;
    xtab2[3] = 0.8602401357;

    k = 0;
    for ( i = 0; i < order2; i++ )
    {
      for ( j = 0; j < order2; j++ )
      {
        xtab[k] = xtab2[j];
        ytab[k] = xtab1[i] * ( 1.0 - xtab2[j] );
        weight[k] = weight1[i] * weight2[j];
        k = k + 1;
      }
    }
  }
//
//  64 points, precision 15.
//
  else if ( rule == 16 )
  {
//
//  Legendre rule of order 8.
//
    order2 = 8;

    xtab1[0] = -0.960289856497536231683560868569;
    xtab1[1] = -0.796666477413626739591553936476;
    xtab1[2] = -0.525532409916328985817739049189;
    xtab1[3] = -0.183434642495649804939476142360;
    xtab1[4] =  0.183434642495649804939476142360;
    xtab1[5] =  0.525532409916328985817739049189;
    xtab1[6] =  0.796666477413626739591553936476;
    xtab1[7] =  0.960289856497536231683560868569;

    weight1[0] = 0.101228536290376259152531354310;
    weight1[1] = 0.222381034453374470544355994426;
    weight1[2] = 0.313706645877887287337962201987;
    weight1[3] = 0.362683783378361982965150449277;
    weight1[4] = 0.362683783378361982965150449277;
    weight1[5] = 0.313706645877887287337962201987;
    weight1[6] = 0.222381034453374470544355994426;
    weight1[7] = 0.101228536290376259152531354310;

    weight2[0] = 0.00329519144;
    weight2[1] = 0.01784290266;
    weight2[2] = 0.04543931950;
    weight2[3] = 0.07919959949;
    weight2[4] = 0.10604735944;
    weight2[5] = 0.11250579947;
    weight2[6] = 0.09111902364;
    weight2[7] = 0.04455080436;

    xtab2[0] = 0.04463395529;
    xtab2[1] = 0.14436625704;
    xtab2[2] = 0.28682475714;
    xtab2[3] = 0.45481331520;
    xtab2[4] = 0.62806783542;
    xtab2[5] = 0.78569152060;
    xtab2[6] = 0.90867639210;
    xtab2[7] = 0.98222008485;

    k = 0;
    for ( i = 0; i < order2; i++ )
    {
      for ( j = 0; j < order2; j++ )
      {
        xtab[k] = 1.0 - xtab2[j];
        ytab[k] = 0.5 * ( 1.0 + xtab1[i] ) * xtab2[j];
        weight[k] = weight1[i] * weight2[j];
        k = k + 1;
      }
    }
  }
//
//  19 points, precision 8, from CUBTRI.
//
  else if ( rule == 17 )
  {
    a = 1.0 / 3.0;
    b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
    d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
    f = ( 40.0 - 10.0 * sqrt ( 15.0 ) 
      + 10.0 * sqrt ( 7.0 ) + 2.0 * sqrt ( 105.0 ) ) / 90.0;
    g = ( 25.0 +  5.0 * sqrt ( 15.0 ) 
      -  5.0 * sqrt ( 7.0 ) - sqrt ( 105.0 ) ) / 90.0;
    p = ( 40.0 + 10.0 * sqrt ( 15.0 ) 
      + 10.0 * sqrt ( 7.0 ) - 2.0 * sqrt ( 105.0 ) ) / 90.0;
    q = ( 25.0 -  5.0 * sqrt ( 15.0 ) 
      -  5.0 * sqrt ( 7.0 ) + sqrt ( 105.0 ) ) / 90.0;
    r = ( 40.0 + 10.0 * sqrt ( 7.0 ) ) / 90.0;
    s = ( 25.0 +  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) 
      - sqrt ( 105.0 ) ) / 90.0;
    t = ( 25.0 -  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) 
      + sqrt ( 105.0 ) ) / 90.0;

    w1 = ( 7137.0 - 1800.0 * sqrt ( 7.0 ) ) / 62720.0;
    w2 = -9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 ) 
      / 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0 
      + 198763.0 * sqrt ( 105.0 ) / 939008.0;
    w2 = w2 / 3.0;
    w3 = -9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 ) 
      / 23475200.0 
      + 764885.0 * sqrt ( 7.0 ) / 939008.0 
      - 198763.0 * sqrt ( 105.0 ) / 939008.0;
    w3 = w3 / 3.0;
    w4 = ( 102791225.0 - 23876225.0 * sqrt ( 15.0 ) 
      - 34500875.0 * sqrt ( 7.0 ) 
      + 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0;
    w4 = w4 / 3.0;
    w5 = ( 102791225.0 + 23876225.0 * sqrt ( 15.0 ) 
      - 34500875.0 * sqrt ( 7.0 ) 
      - 9914825 * sqrt ( 105.0 ) ) / 59157504.0;
    w5 = w5 / 3.0;
    w6 = ( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0;
    w6 = w6 / 6.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = p;
    xtab[11] = q;
    xtab[12] = q;
    xtab[13] = r;
    xtab[14] = r;
    xtab[15] = s;
    xtab[16] = s;
    xtab[17] = t;
    xtab[18] = t;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;
    ytab[10] = q;
    ytab[11] = p;
    ytab[12] = q;
    ytab[13] = s;
    ytab[14] = t;
    ytab[15] = r;
    ytab[16] = t;
    ytab[17] = r;
    ytab[18] = s;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w6;
    weight[17] = w6;
    weight[18] = w6;
  }
//
//  19 points, precision 9.
//  Lyness and Jesperson.
//
  else if ( rule == 18 )
  {
    a = 1.0 / 3.0;
    b =  0.02063496160252593;
    c =  0.4896825191987370;
    d =  0.1258208170141290;
    e =  0.4370895914929355;
    f =  0.6235929287619356;
    g =  0.1882035356190322;
    r =  0.9105409732110941;
    s =  0.04472951339445297;
    t =  0.7411985987844980;
    u =  0.03683841205473626;
    v =  0.22196298916076574;

    w1 = 0.09713579628279610;
    w2 = 0.03133470022713983;
    w3 = 0.07782754100477543;
    w4 = 0.07964773892720910;
    w5 = 0.02557767565869810;
    w6 = 0.04328353937728940;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = r;
    xtab[11] = s;
    xtab[12] = s;
    xtab[13] = t;
    xtab[14] = t;
    xtab[15] = u;
    xtab[16] = u;
    xtab[17] = v;
    xtab[18] = v;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = s;
    ytab[10] = r;
    ytab[11] = s;
    ytab[12] = u;
    ytab[13] = v;
    ytab[14] = t;
    ytab[15] = t;
    ytab[16] = v;
    ytab[17] = t;
    ytab[18] = u;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w6;
    weight[17] = w6;
    weight[18] = w6;
  }
//
//  28 points, precision 11.
//  Lyness and Jesperson.
//
  else if ( rule == 19 )
  {
    a = 1.0 / 3.0;
    b = 0.9480217181434233;
    c = 0.02598914092828833;
    d = 0.8114249947041546;
    e = 0.09428750264792270;
    f = 0.01072644996557060;
    g = 0.4946367750172147;
    p = 0.5853132347709715;
    q = 0.2073433826145142;
    r = 0.1221843885990187;
    s = 0.4389078057004907;
    t = 0.6779376548825902;
    u = 0.04484167758913055;
    v = 0.27722066752827925;
    w = 0.8588702812826364;
    x = 0.0;
    y = 0.1411297187173636;

    w1 = 0.08797730116222190;
    w2 = 0.008744311553736190;
    w3 = 0.03808157199393533;
    w4 = 0.01885544805613125;
    w5 = 0.07215969754474100;
    w6 = 0.06932913870553720;
    w7 = 0.04105631542928860;
    w8 = 0.007362383783300573;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = p;
    xtab[11] = q;
    xtab[12] = q;
    xtab[13] = r;
    xtab[14] = s;
    xtab[15] = s;
    xtab[16] = t;
    xtab[17] = t;
    xtab[18] = u;
    xtab[19] = u;
    xtab[20] = v;
    xtab[21] = v;
    xtab[22] = w;
    xtab[23] = w;
    xtab[24] = x;
    xtab[25] = x;
    xtab[26] = y;
    xtab[27] = y;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;
    ytab[10] = q;
    ytab[11] = p;
    ytab[12] = q;
    ytab[13] = s;
    ytab[14] = r;
    ytab[15] = s;
    ytab[16] = u;
    ytab[17] = v;
    ytab[18] = t;
    ytab[19] = v;
    ytab[20] = t;
    ytab[21] = u;
    ytab[22] = x;
    ytab[23] = y;
    ytab[24] = w;
    ytab[25] = y;
    ytab[26] = w;
    ytab[27] = x;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w7;
    weight[17] = w7;
    weight[18] = w7;
    weight[19] = w7;
    weight[20] = w7;
    weight[21] = w7;
    weight[22] = w8;
    weight[23] = w8;
    weight[24] = w8;
    weight[25] = w8;
    weight[26] = w8;
    weight[27] = w8;
  }
//
//  37 points, precision 13.
//
  else if ( rule == 20 )
  {
    a = 1.0 / 3.0;
    b = 0.950275662924105565450352089520;
    c = 0.024862168537947217274823955239;
    d = 0.171614914923835347556304795551;
    e = 0.414192542538082326221847602214;
    f = 0.539412243677190440263092985511;
    g = 0.230293878161404779868453507244;

    w1 = 0.051739766065744133555179145422;
    w2 = 0.008007799555564801597804123460;
    w3 = 0.046868898981821644823226732071;
    w4 = 0.046590940183976487960361770070;
    w5 = 0.031016943313796381407646220131;
    w6 = 0.010791612736631273623178240136;
    w7 = 0.032195534242431618819414482205;
    w8 = 0.015445834210701583817692900053;
    w9 = 0.017822989923178661888748319485;
    wx = 0.037038683681384627918546472190;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w7;
    weight[17] = w7;
    weight[18] = w7;
    weight[19] = w8;
    weight[20] = w8;
    weight[21] = w8;
    weight[22] = w8;
    weight[23] = w8;
    weight[24] = w8;
    weight[25] = w9;
    weight[26] = w9;
    weight[27] = w9;
    weight[28] = w9;
    weight[29] = w9;
    weight[30] = w9;
    weight[31] = wx;
    weight[32] = wx;
    weight[33] = wx;
    weight[34] = wx;
    weight[35] = wx;
    weight[36] = wx;

    a = 0.772160036676532561750285570113;
    b = 0.113919981661733719124857214943;

    xtab[10] = a;
    ytab[10] = b;

    xtab[11] = b;
    ytab[11] = a;

    xtab[12] = b;
    ytab[12] = b;

    a = 0.009085399949835353883572964740;
    b = 0.495457300025082323058213517632;

    xtab[13] = a;
    ytab[13] = b;

    xtab[14] = b;
    ytab[14] = a;

    xtab[15] = b;
    ytab[15] = b;

    a = 0.062277290305886993497083640527;
    b = 0.468861354847056503251458179727;

    xtab[16] = a;
    ytab[16] = b;

    xtab[17] = b;
    ytab[17] = a;

    xtab[18] = b;
    ytab[18] = b;

    a = 0.022076289653624405142446876931;
    b = 0.851306504174348550389457672223;
    c = 0.126617206172027096933163647918263;

    xtab[19] = a;
    ytab[19] = b;

    xtab[20] = a;
    ytab[20] = c;

    xtab[21] = b;
    ytab[21] = a;

    xtab[22] = b;
    ytab[22] = c;

    xtab[23] = c;
    ytab[23] = a;

    xtab[24] = c;
    ytab[24] = b;

    a = 0.018620522802520968955913511549;
    b = 0.689441970728591295496647976487;
    c = 0.291937506468887771754472382212953;

    xtab[25] = a;
    ytab[25] = b;

    xtab[26] = a;
    ytab[26] = c;

    xtab[27] = b;
    ytab[27] = a;

    xtab[28] = b;
    ytab[28] = c;

    xtab[29] = c;
    ytab[29] = a;

    xtab[30] = c;
    ytab[30] = b;

    a = 0.096506481292159228736516560903;
    b = 0.635867859433872768286976979827;
    c = 0.267625659273967961282458816185681;

    xtab[31] = a;
    ytab[31] = b;

    xtab[32] = a;
    ytab[32] = c;

    xtab[33] = b;
    ytab[33] = a;

    xtab[34] = b;
    ytab[34] = c;

    xtab[35] = c;
    ytab[35] = a;

    xtab[36] = c;
    ytab[36] = b;
  }
  else
  {
    cerr << "\n";
    cerr << "TRIANGLE_UNIT_SET - Fatal error!\n";
    cerr << "  Illegal value of RULE = " << rule << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

int triangle_unit_size ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jarle Berntsen, Terje Espelid,
//    Algorithm 706,
//    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
//    ACM Transactions on Mathematical Software,
//    Volume 18, Number 3, September 1992, pages 329-342.
//
//    Elise deDoncker, Ian Robinson,
//    Algorithm 612:
//    Integration over a Triangle Using Nonlinear Extrapolation,
//    ACM Transactions on Mathematical Software,
//    Volume 10, Number 1, March 1984, pages 17-22.
//
//    DP Laurie,
//    Algorithm 584,
//    CUBTRI, Automatic Cubature Over a Triangle,
//    ACM Transactions on Mathematical Software,
//    Volume 8, Number 2, 1982, pages 210-218.
//
//    James Lyness, Dennis Jespersen,
//    Moderate Degree Symmetric Quadrature Rules for the Triangle,
//    Journal of the Institute of Mathematics and its Applications,
//    Volume 15, Number 1, February 1975, pages 19-32.
//
//    Hans Rudolf Schwarz,
//    Methode der Finiten Elemente,
//    Teubner Studienbuecher, 1980,
//    ISBN: 3-519-02349-0.
//
//    Gilbert Strang, George Fix,
//    An Analysis of the Finite Element Method,
//    Prentice Hall, 1973, page 184,
//    ISBN: 096140888X,
//    LC: TA335.S77.
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//    ISBN: 0750663200,
//    TA640.2.Z54
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//     1, ORDER =  1, precision 1, Zienkiewicz #1.
//     2, ORDER =  2, precision 1, (the "vertex rule").
//     3, ORDER =  3, precision 2, Strang and Fix formula #1.
//     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
//     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
//     6, ORDER =  6, precision 3, Strang and Fix formula #4.
//     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
//     8, ORDER =  6, precision 4, Strang and Fix formula #5.
//     9, ORDER =  7, precision 4, Strang and Fix formula #6.
//    10, ORDER =  7, precision 5, Strang and Fix formula #7,
//        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
//    11, ORDER =  9, precision 6, Strang and Fix formula #8.
//    12, ORDER = 12, precision 6, Strang and Fix formula #9.
//    13, ORDER = 13, precision 7, Strang and Fix formula #10.
//    14, ORDER =  7, precision ?.
//    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
//    16, ORDER = 64, precision 15, triangular product Gauss rule.
//    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
//    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
//    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
//    20, ORDER = 37, precision 13, from ACM TOMS #706.
//
//    Output, int TRIANGLE_UNIT_SIZE, the order of the rule.
//
{
  int size;

  if ( rule == 1 )
  {
    size = 1;
  }
  else if ( rule == 2 )
  {
    size = 3;
  }
  else if ( rule == 3 )
  {
    size = 3;
  }
  else if ( rule == 4 )
  {
    size = 3;
  }
  else if ( rule == 5 )
  {
    size = 4;
  }
  else if ( rule == 6 )
  {
    size = 6;
  }
  else if ( rule == 7 )
  {
    size = 6;
  }
  else if ( rule == 8 )
  {
    size = 6;
  }
  else if ( rule == 9 )
  {
    size = 7;
  }
  else if ( rule == 10 )
  {
    size = 7;
  }
  else if ( rule == 11 )
  {
    size = 9;
  }
  else if ( rule == 12 )
  {
    size = 12;
  }
  else if ( rule == 13 )
  {
    size = 13;
  }
  else if ( rule == 14 )
  {
    size = 7;
  }
  else if ( rule == 15 )
  {
    size = 16;
  }
  else if ( rule == 16 )
  {
    size = 64;
  }
  else if ( rule == 17 )
  {
    size = 19;
  }
  else if ( rule == 18 )
  {
    size = 19;
  }
  else if ( rule == 19 )
  {
    size = 28;
  }
  else if ( rule == 20 )
  {
    size = 37;
  }
  else
  {
    size = -1;
  }

  return size;
}
//****************************************************************************80

double triangle_unit_sum ( double func ( double x, double y ), int order, 
  double xtab[], double ytab[], double weight[] )

//****************************************************************************80
//
//// TRIANGLE_UNIT_SUM carries out a quadrature rule in the unit triangle.
//
//  Integration region:
//
//      0 <= X,
//    and
//      0 <= Y, 
//    and
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FUNC ( double x, double y ), the name of the
//    user supplied function to be integrated.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
//
//    Input, double WEIGHT[ORDER], the weights of the rule.
//
//    Output, double RESULT, the approximate integral of the function.
//
{
  int i;
  double quad;
  double result;
  double volume;

  quad = 0.0;
  for ( i = 0; i < order; i++ )
  {
    quad = quad + weight[i] * func ( xtab[i], ytab[i] );
  }

  volume = triangle_unit_volume ( );
  result = quad * volume;

  return result;
}
//****************************************************************************80

double triangle_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
//
//  Integration region:
//
//      0 <= X,
//    and
//      0 <= Y, 
//    and
//      X + Y <= 1.
//
//  Discussion:
//
//    The "volume" of a triangle is usually called its area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TRIANGLE_UNIT_VOLUME, the volume of the unit
//    triangle.
//
{
  double volume;

  volume = 1.0 / 2.0;

  return volume;
}
//****************************************************************************80

double triangle_volume ( double x[3], double y[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_VOLUME returns the "volume" of a triangle in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X[3], Y[3], the vertices of the triangle.
//
//    Output, double TRIANGLE_VOLUME, the volume of the triangle.
//
{
  double value;

  value = 0.5 * r8_abs ( 
    x[0] * ( y[1] - y[2] ) + 
    x[1] * ( y[2] - y[0] ) + 
    x[2] * ( y[0] - y[1] ) );

  return value;
}
//****************************************************************************80

double *tvec_even ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN computes an evenly spaced set of angles between 0 and 2*PI.
//
//  Discussion:
//
//    The computation realizes that 0 = 2 * PI.
//
//  Example:
//
//    NT = 4
//
//    T = ( 0, PI/2, PI, 3*PI/2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Output, double TVEC[NT], the evenly spaced angles, in radians.
//
{
  int i;
  double pi = 3.141592653589793;
  double *t;

  if ( nt < 1 )
  {
    return NULL;
  }

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( double ) ( 2 * ( i - 1 ) ) * pi / ( double ) ( nt );
  }

  return t;
}
//****************************************************************************80

double *tvec_even2 ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN2 computes evenly spaced angles between 0 and 2*PI.
//
//  Discussion:
//
//    The computation realizes that 0 = 2 * PI.  The values are equally
//    spaced in the circle, do not include 0, and are symmetric about 0.
//
//  Example:
//
//    NT = 4
//
//    T = ( PI/4, 3*PI/4, 5*PI/4, 7*PI/4 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Output, double TVEC[NT], the evenly spaced angles, in radians.
//
{
  int i;
  double pi = 3.141592653589793;
  double *t;

  if ( nt < 1 )
  {
    return NULL;
  }

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( double ) ( 2 * i - 1 ) * pi / ( double ) ( nt );
  }

  return t;
}
//****************************************************************************80

double *tvec_even3 ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN3 computes an evenly spaced set of angles between 0 and 2*PI.
//
//  Discussion:
//
//    The angles begin with 0 and end with 2*PI.
//
//  Example:
//
//    NT = 4
//
//    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Output, double TVEC[NT], the evenly spaced angles, in radians.
//
{
  int i;
  double pi = 3.141592653589793;
  double *t;

  if ( nt < 1 )
  {
    return NULL;
  }

  t = new double[nt];

  if ( nt == 1 )
  {
    t[0] = pi;
  }
  else
  {
    for ( i = 1; i <= nt; i++ )
    {
      t[i-1] = ( double ) ( 2 * ( i - 1 ) ) * pi / ( double ) ( nt - 1 );
    }
  }

  return t;
}
//****************************************************************************80

double *tvec_even_bracket ( int nt, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
//
//  Example:
//
//    NT = 4
//    THETA1 = 30
//    THETA2 = 90
//
//    T = ( 30, 50, 70, 90 )
//
//  Discussion:
//
//    The interval between THETA1 and THETA2 is divided into NT-1 subintervals.
//
//    The angles returned are the breakpoints of these subintervals,
//    including THETA1 and THETA2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Input, double THETA1, THETA2, the limiting angles.
//
//    Output, double TVEC_EVEN_BRACKET[NT], the evenly spaced angles.
//
{
  int i;
  double *t;

  if ( nt < 1 ) 
  {
    return NULL;
  }

  t = new double[nt];

  if ( nt == 1 )
  {
    t[0] = ( theta1 + theta2 ) / 2.0;
  }
  else
  {
    for ( i = 1; i <= nt; i++ )
    {
      t[i-1] = ( ( double ) ( nt - i     ) * theta1   
               + ( double ) (      i - 1 ) * theta2 ) 
               / ( double ) ( nt     - 1 );
    }
  }

  return t;
}
//****************************************************************************80

double *tvec_even_bracket2 ( int nt, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
//
//  Discussion:
//
//    The interval between THETA1 and THETA2 is divided into NT+1 subintervals.
//
//    The angles returned are the internal NT breakpoints of the subintervals.
//
//  Example:
//
//    NT = 5
//    THETA1 = 30
//    THETA2 = 90
//
//    T = ( 40, 50, 60, 70, 80 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Input, double THETA1, THETA2, the limiting angles.
//
//    Output, double TVEC_EVEN_BRACKET2[NT], the evenly spaced angles.
//
{
  int i;
  double *t;

  if ( nt < 1 ) 
  {
    return NULL;
  }

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( ( double ) ( nt + 1 - i ) * theta1   
             + ( double ) (          i ) * theta2 ) 
             / ( double ) ( nt + 1     );
  }

  return t;
}
//****************************************************************************80

double *tvec_even_bracket3 ( int nt, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN_BRACKET3 computes evenly spaced angles from THETA1 to THETA2.
//
//  Discussion:
//
//    The interval between THETA1 and THETA2 is divided into NT subintervals.
//
//    The angles returned are the midpoints of each subinterval.
//
//  Example:
//
//    NT = 3
//    THETA1 = 30
//    THETA2 = 90
//
//    T = ( 40, 60, 80 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Input, double THETA1, THETA2, the limiting angles.
//
//    Output, double TVEC_EVEN_BRACKET3[NT], the evenly spaced angles.
//
{
  int i;
  double *t;

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( ( double ) ( 2 * nt - 2 * i + 1 ) * theta1   
             + ( double ) (          2 * i - 1 ) * theta2 ) 
             / ( double ) ( 2 * nt             );
  }

  return t;
}
//****************************************************************************80

void vec_lex_next ( int dim_num, int base, int a[], bool *more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_LEX_NEXT generates vectors in lex order.
//
//  Discussion:
//
//    The vectors are produced in lexical order, starting with
//    (0,0,...,0),
//    (0,0,...,1), 
//    ...
//    (BASE-1,BASE-1,...,BASE-1).
//
//  Examples:
//
//    DIM_NUM = 2,
//    BASE = 3
//
//    0   0
//    0   1
//    0   2
//    1   0
//    1   1
//    1   2
//    2   0
//    2   1
//    2   2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the size of the vectors to be used.
//
//    Input, int BASE, the base to be used.  BASE = 2 will
//    give vectors of 0's and 1's, for instance.
//
//    Output, int A[DIM_NUM], the next vector.  
//
//    Input/output, bool *MORE.  Set this variable false before
//    the first call.  On return, MORE is TRUE if another vector has
//    been computed.  If MORE is returned FALSE, ignore the output 
//    vector and stop calling the routine.
//
{
  int i;

  if ( !( *more ) )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 0;
    }
    *more = true;
  }
  else
  {
    for ( i = dim_num - 1; 0 <= i; i-- )
    {
      a[i] = a[i] + 1;

      if ( a[i] < base )
      {
        return;
      }
      a[i] = 0;
    }
    *more = false;
  }

  return;
}
