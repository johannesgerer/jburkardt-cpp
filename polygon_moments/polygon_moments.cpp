# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "polygon_moments.hpp"

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

double moment ( int n, double x[], double y[], int p, int q )

//****************************************************************************80
//
//  Purpose:
//
//    MOMENT computes an unnormalized moment of a polygon.
//
//  Discussion:
//
//    Nu(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carsten Steger,
//    On the calculation of arbitrary moments of polygons,
//    Technical Report FGBV-96-05,
//    Forschungsgruppe Bildverstehen, Informatik IX,
//    Technische Universitaet Muenchen, October 1996.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//
//    Input, double X[N], Y[N], the vertex coordinates.
//
//    Input, int P, Q, the indices of the moment.
//
//    Output, double MOMENT, the unnormalized moment Nu(P,Q).
//
{
  int i;
  int k;
  int l;
  double nu_pq;
  double s_pq;
  double xi;
  double xj;
  double yi;
  double yj;

  nu_pq = 0.0;

  xj = x[n-1];
  yj = y[n-1];

  for ( i = 0; i < n; i++ )
  {
    xi = x[i];
    yi = y[i];

    s_pq = 0.0;
    for ( k = 0; k <= p; k++ )
    {
      for ( l = 0; l <= q; l++ )
      {
        s_pq = s_pq 
          + r8_choose ( k + l, l ) * r8_choose ( p + q - k - l, q - l ) 
          * pow ( xi, k ) * pow ( xj, p - k ) 
          * pow ( yi, l ) * pow ( yj, q - l );
      }
    }

    nu_pq = nu_pq + ( xj * yi - xi * yj ) * s_pq;

    xj = xi;
    yj = yi;
  }

  nu_pq = nu_pq / ( double ) ( p + q + 2 ) / ( double ) ( p + q + 1 ) 
    / r8_choose ( p + q, p );

  return nu_pq;
}
//****************************************************************************80

double moment_central ( int n, double x[], double y[], int p, int q )

//****************************************************************************80
//
//  Purpose:
//
//    MOMENT_CENTRAL computes central moments of a polygon.
//
//  Discussion:
//
//    The central moment Mu(P,Q) is defined by
//
//      Mu(P,Q) = Integral ( polygon ) (x-Alpha(1,0))^p (y-Alpha(0,1))^q dx dy
//              / Area ( polygon )
//
//    where 
//
//      Alpha(1,0) = Integral ( polygon ) x dx dy / Area ( polygon )
//      Alpha(0,1) = Integral ( polygon ) y dx dy / Area ( polygon )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carsten Steger,
//    On the calculation of arbitrary moments of polygons,
//    Technical Report FGBV-96-05,
//    Forschungsgruppe Bildverstehen, Informatik IX,
//    Technische Universitaet Muenchen, October 1996.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//
//    Input, double X[N], Y[N], the vertex coordinates.
//
//    Input, int P, Q, the indices of the moment.
//
//    Output, double MOMENT_CENTRAL, the unnormalized moment Mu(P,Q).
//
{
  double alpha_01;
  double alpha_10;
  double alpha_ij;
  int i;
  int j;
  double mu_pq;

  alpha_10 = moment_normalized ( n, x, y, 1, 0 );
  alpha_01 = moment_normalized ( n, x, y, 0, 1 );

  mu_pq = 0.0;

  for ( i = 0; i <= p; i++ )
  {
    for ( j = 0; j <= q; j++ )
    {
      alpha_ij = moment_normalized ( n, x, y, i, j );

      mu_pq = mu_pq + r8_mop ( p + q - i - j ) 
        * r8_choose ( p, i ) * r8_choose ( q, j ) 
        * pow ( alpha_10, p - i ) * pow ( alpha_01, q - j ) * alpha_ij;
    }
  }

  return mu_pq;
}
//****************************************************************************80

double moment_normalized ( int n, double x[], double y[], int p, int q )

//****************************************************************************80
//
//  Purpose:
//
//    MOMENT_NORMALIZED computes a normalized moment of a polygon.
//
//  Discussion:
//
//    Alpha(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy / Area ( polygon )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carsten Steger,
//    On the calculation of arbitrary moments of polygons,
//    Technical Report FGBV-96-05,
//    Forschungsgruppe Bildverstehen, Informatik IX,
//    Technische Universitaet Muenchen, October 1996.
//
//  Parameters:
//
//    Input, int N, the number of vertices of the polygon.
//
//    Input, double X[N], Y[N], the vertex coordinates.
//
//    Input, int P, Q, the indices of the moment.
//
//    Output, double MOMENT_NORMALIZED, the normalized moment Alpha(P,Q).
//
{
  double alpha_pq;
  double nu_00;
  double nu_pq;

  nu_pq = moment ( n, x, y, p, q );
  nu_00 = moment ( n, x, y, 0, 0 );

  alpha_pq = nu_pq / nu_00;

  return alpha_pq;
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
//    29 July 2011
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
  double value;

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
