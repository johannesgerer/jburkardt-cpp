# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "piecewise_linear_product_integral.hpp"

//****************************************************************************80

double piecewise_linear_product_integral ( double a, double b, int f_num, 
  double f_x[], double f_v[], int g_num, double g_x[], double g_v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PIECEWISE_LINEAR_PRODUCT_INTEGRAL: piecewise linear product integral.
//
//  Discussion:
//
//    We are given two piecewise linear functions F(X) and G(X) and we wish
//    to compute the exact value of the integral
//
//      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
//
//    The functions F(X) and G(X) are defined as tables of coordinates X and
//    values V.  A piecewise linear function is evaluated at a point X by 
//    evaluating the interpolant to the data at the endpoints of the interval 
//    containing X.  
//
//    It must be the case that A <= B.
//
//    It must be the case that the node coordinates F_X(*) and G_X(*) are
//    given in ascending order.
//
//    It must be the case that:
//
//      F_X(1) <= A and B <= F_X(F_NUM)
//      G_X(1) <= A and B <= G_X(G_NUM)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of integration.
//
//    Input, int F_NUM, the number of nodes for F.
//
//    Input, double F_X[F_NUM], the node coordinates for F.
//
//    Input, double F_V[F_NUM], the nodal values for F.
//
//    Input, int G_NUM, the number of nodes for G.
//
//    Input, double G_X[G_NUM], the node coordinates for G.
//
//    Input, double G_V[G_NUM], the nodal values for G.
//
//    Output, double INTEGRAL, the integral of F(X) * G(X)
//    from A to B.
//
{
  double bit;
  int f_left;
  double f0;
  double f1;
  double fl;
  double fr;
  int g_left;
  double g0;
  double g1;
  double gl;
  double gr;
  double h0;
  double h1;
  double h2;
  int i;
  double integral;
  double xl;
  double xr;
  double xr_max;

  integral = 0.0;

  if ( f_x[f_num-1] <= a || g_x[g_num-1] <= a )
  {
    return integral;
  }

  if ( f_num < 2 || g_num < 2 )
  {
    return integral;
  }

  xr = a;

  f_left = 0;
  r8vec_bracket3 ( f_num, f_x, xr, &f_left );
  fr = f_v[f_left] + ( xr - f_x[f_left] ) * ( f_v[f_left+1] - f_v[f_left] ) 
    / ( f_x[f_left+1] - f_x[f_left] );

  g_left = 0;
  r8vec_bracket3 ( g_num, g_x, xr, &g_left );
  gr = g_v[g_left] + ( xr - g_x[g_left] ) * ( g_v[g_left+1] - g_v[g_left] ) 
    / ( g_x[g_left+1] - g_x[g_left] );

  xr_max = b;
  xr_max = r8_min ( xr_max, f_x[f_num-1] );
  xr_max = r8_min ( xr_max, g_x[g_num-1] );

  while ( xr < xr_max )
  {
//
//  Shift right values to left.
//
    xl = xr;
    fl = fr;
    gl = gr;
//
//  Determine the new right values.
//  The hard part is figuring out how to advance XR some, but not too much.
//
    xr = xr_max;

    for ( i = 1; i <= 2; i++ )
    {
      if ( f_left + i <= f_num - 1 )
      {
        if ( xl < f_x[f_left+i] && f_x[f_left+i] < xr )
        {
          xr = f_x[f_left+i];
          break;
        }
      }
    }

    for ( i = 1; i <= 2; i++ )
    {
      if ( g_left + i <= g_num - 1 )
      {
        if ( xl < g_x[g_left+i] && g_x[g_left+i] < xr )
        {
          xr = g_x[g_left+i];
          break;
        }
      }
    }

    r8vec_bracket3 ( f_num, f_x, xr, &f_left );
    fr = f_v[f_left] + ( xr - f_x[f_left] ) * ( f_v[f_left+1] - f_v[f_left] ) 
      / ( f_x[f_left+1] - f_x[f_left] );

    r8vec_bracket3 ( g_num, g_x, xr, &g_left );
    gr = g_v[g_left] + ( xr - g_x[g_left] ) * ( g_v[g_left+1] - g_v[g_left] ) 
      / ( g_x[g_left+1] - g_x[g_left] );
//
//  Form the linear polynomials for F(X) and G(X) over [XL,XR],
//  then the product H(X), integrate H(X) and add to the running total.
//
    if ( r8_epsilon ( ) <= r8_abs ( xr - xl ) )
    {
      f1 = fl - fr;
      f0 = fr * xl - fl * xr;

      g1 = gl - gr;
      g0 = gr * xl - gl * xr;

      h2 = f1 * g1;
      h1 = f1 * g0 + f0 * g1;
      h0 = f0 * g0;

      h2 = h2 / 3.0;
      h1 = h1 / 2.0;

      bit = ( ( h2 * xr + h1 ) * xr + h0 ) * xr 
          - ( ( h2 * xl + h1 ) * xl + h0 ) * xl;

      integral = integral + bit / ( xr - xl ) / ( xr - xl );
    }
  }

  return integral;
}
//****************************************************************************80

double piecewise_linear_product_quad ( double a, double b, int f_num, 
  double f_x[], double f_v[], int g_num, double g_x[], double g_v[], 
  int quad_num )

//****************************************************************************80
//
//  Purpose:
//
//    PIECEWISE_LINEAR_PRODUCT_QUAD: estimate piecewise linear product integral.
//
//  Discussion:
//
//    We are given two piecewise linear functions F(X) and G(X) and we wish
//    to estimate the value of the integral
//
//      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
//
//    The functions F(X) and G(X) are defined as tables of coordinates X and
//    values V.  A piecewise linear function is evaluated at a point X by 
//    evaluating the interpolant to the data at the endpoints of the interval 
//    containing X.  
//
//    It must be the case that A <= B.
//
//    It must be the case that the node coordinates F_X(*) and G_X(*) are
//    given in ascending order.
//
//    It must be the case that:
//
//      F_X(1) <= A and B <= F_X(F_NUM)
//      G_X(1) <= A and B <= G_X(G_NUM)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of integration.
//
//    Input, int F_NUM, the number of nodes for F.
//
//    Input, double F_X[F_NUM], the node coordinates for F.
//
//    Input, double F_V[F_NUM], the nodal values for F.
//
//    Input, int G_NUM, the number of nodes for G.
//
//    Input, double G_X[G_NUM], the node coordinates for G.
//
//    Input, double G_V[G_NUM], the nodal values for G.
//
//    Input, int QUAD_NUM, the number of quadrature points.
//
//    Output, double PIECEWISE_LINEAR_PRODUCT_QUAD, an estimate for the integral 
//    of F(X) * G(X) from A to B.
//
{
  double a2;
  double b2;
  int f_left;
  double fq;
  int g_left;
  double gq;
  int i;
  double quad;
  double xq;

  quad = 0.0;

  f_left = 0;
  g_left = 0;

  a2 = a;
  a2 = r8_max ( a2, f_x[0] );
  a2 = r8_max ( a2, g_x[0] );

  b2 = b;
  b2 = r8_min ( b2, f_x[f_num-1] );
  b2 = r8_min ( b2, g_x[g_num-1] );

  for ( i = 1; i <= quad_num; i++ )
  {
    xq =  ( ( double ) (                2 * i - 1 ) * b2 
          + ( double ) ( 2 * quad_num - 2 * i + 1 ) * a2 )  
          / ( double ) ( 2 * quad_num             );

    r8vec_bracket3 ( f_num, f_x, xq, &f_left );

    fq = f_v[f_left] + ( xq - f_x[f_left] ) * ( f_v[f_left+1] - f_v[f_left] )
      / ( f_x[f_left+1] - f_x[f_left] );

    r8vec_bracket3 ( g_num, g_x, xq, &g_left );

    gq = g_v[g_left] + ( xq - g_x[g_left] ) * ( g_v[g_left+1] - g_v[g_left] ) 
      / ( g_x[g_left+1] - g_x[g_left] );

    quad = quad + fq * gq;
  }

  quad = quad * ( b - a ) / ( double ) ( quad_num );

  return quad;
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

void r8vec_bracket3 ( int n, double t[], double tval, int *left )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET3 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[N-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    This version of the function has been revised so that the value of
//    LEFT that is returned uses the 0-based indexing natural to C++.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, double T[N], an array that has been sorted into ascending order.
//
//    Input, double TVAL, a value to be bracketed by entries of T.
//
//    Input/output, int *LEFT.
//    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
//    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
//    is searched first, followed by the appropriate interval to the left
//    or right.  After that, a binary search is used.
//    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
//    is the closest to TVAL; it either contains TVAL, or else TVAL
//    lies outside the interval [ T[0], T[N-1] ].
//
{
  int high;
  int low;
  int mid;
//  
//  Check the input data.
//
  if ( n < 2 ) 
  {
    cout << "\n";
    cout << "R8VEC_BRACKET3 - Fatal error//\n";
    cout << "  N must be at least 2.\n";
    exit ( 1 );
  }
//
//  If *LEFT is not between 0 and N-2, set it to the middle value.
//
  if ( *left < 0 || n - 2 < *left ) 
  {
    *left = ( n - 1 ) / 2;
  }
//
//  CASE 1: TVAL < T[*LEFT]:
//  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
//
  if ( tval < t[*left] ) 
  {
    if ( *left == 0 ) 
    {
      return;
    }
    else if ( *left == 1 ) 
    {
      *left = 0;
      return;
    }
    else if ( t[*left-1] <= tval )
    {
      *left = *left - 1;
      return;
    }
    else if ( tval <= t[1] ) 
    {
      *left = 0;
      return;
    }
// 
//  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
//
    low = 1;
    high = *left - 2;

    for ( ; ; )
    {
      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval ) 
      {
        low = mid;
      }
      else 
      {
        high = mid - 1;
      }
    }
  }
// 
//  CASE 2: T[*LEFT+1] < TVAL:
//  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
//
  else if ( t[*left+1] < tval ) 
  {
    if ( *left == n - 2 ) 
    {
      return;
    }
    else if ( *left == n - 3 ) 
    {
      *left = *left + 1;
      return;
    }
    else if ( tval <= t[*left+2] )
    {
      *left = *left + 1;
      return;
    }
    else if ( t[n-2] <= tval ) 
    {
      *left = n - 2;
      return;
    }
// 
//  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
//
    low = *left + 2;
    high = n - 3;

    for ( ; ; ) 
    {

      if ( low == high ) 
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval ) 
      {
        low = mid;
      }
      else 
      {
        high = mid - 1;
      }
    }
  }
//
//  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
//  T is just where the user said it might be.
//
  else 
  {
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
