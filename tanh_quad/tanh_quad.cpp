# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "tanh_quad.hpp"

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
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
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
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
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

void midpoint_rule ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    MIDPOINT_RULE computes a midpoint quadrature rule.
//
//  Discussion:
//
//    2*N+1 equally spaced points in [-1,1] are defined, but the endpoints
//    are not used in the rule.
//
//    This rule is useful in cases where there is a singularity at one
//    or both endpoints of the interval.
//
//    This rule will nest, with N = 0, 1, 3, 7, 15, 31, ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the quadrature order.
//
//    Output, double X[-N:N], the abscissas.
//
//    Output, double W[-N:N], the weights.
//
{
  int i;
  int offset;

  offset = n;

  for ( i = -n; i <= n; i++ )
  {
    x[i+offset] = ( double ) ( i ) / ( double ) ( n + 1 );
  }

  for ( i = -n; i <= n; i++ )
  {
    w[i+offset] = 2.0 / ( double ) ( 2 * ( n + 1 ) );
  }
//
//  If I correct W(-N) and W(N) in this slightly awkward way,
//  I correctly include the special case where N = 0.
//
  w[-n+offset] = w[-n+offset] + 1.0 / ( double ) ( 2 * ( n + 1 ) );
  w[ n+offset] = w[ n+offset] + 1.0 / ( double ) ( 2 * ( n + 1 ) );

  return;
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
    value = - x;
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

double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
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
//    Output, double R8VEC_DOT, the dot product of the vectors.
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
//    If the weight function W(X) is not 1, then the weight vector W will
//    require further adjustment by the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2009
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
           / ( b        - a );
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = ( ( d - c ) / ( b - a ) ) * w[i];
  }
  return;
}
//****************************************************************************80

int tanh_h_to_n ( double h, double tol )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_H_TO_N computes N as a function of H and TOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double H, the spacing.
//
//    Input, double TOL, the tolerance.
//
//    Output, int N, the corresponding quadrature order.
//
{
  double ct;
  int n;
  double pi = 3.141592653589793;
  double t;
  double w;

  n = 0;

  for ( ; ; )
  {
    t = ( double ) ( n ) * h / 2.0;

    ct = cosh ( t );

    w = 0.5 * h / ct / ct;

    if ( w <= tol )
    {
      break;
    }
    n = n + 1;
  }
  return n;
}
//****************************************************************************80

double tanh_m_to_h ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_M_TO_H computes H as a function of M.
//
//  Discussion:
//
//    H = 2^(-M).
//
//    This is simply an orderly way to index a family of decreasing values of H.
//
//  Example:
//
//     M      H
//    --  -----
//    -2      4
//    -1      2
//     0      1
//     1     1/2
//     2     1/4
//     3     1/8
//     4     1/16
//     5     1/32
//   ...     ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the level.
//
//    Output, double H, the spacing.
//
{
  double h;
  int i;

  h = 1.0;

  for ( i = -1; m <= i; i-- )
  {
    h = h * 2.0;
  }

  for ( i = 1; i <= m; i++ )
  {
    h = h / 2.0;
  }

  return h;
}
//****************************************************************************80

double tanh_n_to_h ( int n  )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_N_TO_H computes N as a function of H.
//
//  Discussion:
//
//    This formula for N(H) is suggested in Kahaner, Moler and Nash.
//    Note, however, that using this formula means that it is not possible
//    to make families of nested rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int N, the quadrature order.
//
//    Output, double H, the spacing.
//
{
  double h;
  double pi = 3.141592653589793;

  h = pi * sqrt ( 2.0 / ( double ) ( n ) ) - 1.0 / ( double ) ( n );

  return h;
}
//****************************************************************************80

void tanh_rule ( int n, double h, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_RULE computes a tanh-sinh quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the quadrature order.
//
//    Input, double H, the spacing.
//
//    Output, double X[-N:N], the abscissas.
//
//    Output, double W[-N:N], the weights.
//
{
  double ct;
  int i;
  int offset;
  double pi = 3.141592653589793;
  double t;

  offset = n;

  for ( i = -n; i <= n; i++ )
  {
    t = ( double ) ( i ) * h / 2.0;

    ct = cosh ( t );

    x[i+offset] = tanh ( t );
    w[i+offset] = 0.5 * h / ct / ct;
  }

  return;
}
//****************************************************************************80

int tanh_sinh_h_to_n ( double h, double tol )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_SINH_H_TO_N computes N as a function of H and TOL for the tanh-sinh rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double H, the spacing.
//
//    Input, double TOL, the tolerance.
//
//    Output, int TANH_SINH_H_TO_N, the corresponding quadrature order.
//
{
  double ct;
  double ct2;
  int n;
  double pi = 3.141592653589793;
  double st;
  double t;
  double w;

  n = 0;

  for ( ; ; )
  {
    t = ( double ) ( n ) * h;

    ct = cosh ( t );
    st = sinh ( t );
    ct2 = cosh ( 0.5 * pi * st );

    w = 0.5 * pi * h * ct / ct2 / ct2;

    if ( w <= tol )
    {
      break;
    }

    n = n + 1;
  }

  return n;
}
//****************************************************************************80

void tanh_sinh_rule ( int n, double h, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_SINH_RULE computes a tanh-sinh quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the quadrature order.
//
//    Input, double H, the spacing.
//
//    Output, double X[-N:N], the abscissas.
//
//    Output, double W[-N:N], the weights.
//
{
  double ct;
  double ct2;
  int i;
  int offset;
  double pi = 3.141592653589793;
  double st;
  double t;

  offset = n;

  for ( i = -n; i <= n; i++ )
  {
    t = ( double ) ( i ) * h;

    ct = cosh ( t );
    st = sinh ( t );
    ct2 = cosh ( 0.5 * pi * st );

    x[i+offset] = tanh ( 0.5 * pi * st );

    w[i+offset] = 0.5 * pi * h * ct / ct2 / ct2;
  }
  return;
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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
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

void trap_rule ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRAP_RULE computes a trapezoid quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the quadrature order.
//
//    Output, double X[-N:N], the abscissas.
//
//    Output, double W[-N:N], the weights.
//
{
  int i;
  int offset;

  offset = n;

  for ( i = -n; i <= n; i++ )
  {
    x[i+offset] = ( double ) ( i ) / ( double ) ( n );
  }

  w[-n+offset] = 1.0;
  for ( i = -n+1; i <= n-1; i++ )
  {
    w[i+offset] = 2.0;
  }
  w[+n+offset] = 1.0;

  for ( i = -n; i <= n; i++ )
  {
    w[i+offset] = w[i+offset] / ( double ) ( 2 * n );
  }
  return;
}
