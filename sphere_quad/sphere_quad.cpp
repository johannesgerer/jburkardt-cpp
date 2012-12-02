# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sphere_quad.hpp"

//****************************************************************************80

double arc_cosine ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_COSINE computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double ARC_COSINE, an angle whose cosine is C.
//
{
  double angle;
  double pi = 3.141592653589793;

  if ( c <= -1.0 )
  {
    angle = pi;
  } 
  else if ( 1.0 <= c )
  {
    angle = 0.0;
  }
  else
  {
    angle = acos ( c );
  }
  return angle;
}
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
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the argument, the sine of an angle.
//
//    Output, double ARC_SINE, an angle whose sine is S.
//
{
  double angle;
  double pi = 3.141592653589793;

  if ( s <= -1.0 )
  {
    angle = - pi / 2.0;
  } 
  else if ( 1.0 <= s )
  {
    angle = pi / 2.0;
  }
  else
  {
    angle = asin ( s );
  }
  return angle;
}
//****************************************************************************80

double atan4 ( double y, double x )

//****************************************************************************80
//
//  Purpose:
//
//    ATAN4 computes the inverse tangent of the ratio Y / X.
//
//  Discussion:
//
//    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
//    the built in functions ATAN and ATAN2 already do.
//
//    However:
//
//    * ATAN4 always returns a positive angle, between 0 and 2 PI,
//      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
//      and [-PI,+PI] respectively;
//
//    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
//     function by contrast always returns an angle in the first or fourth
//     quadrants.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double Y, X, two quantities which represent the tangent of
//    an angle.  If Y is not zero, then the tangent is (Y/X).
//
//    Output, double ATAN4, an angle between 0 and 2 * PI, whose tangent is
//    (Y/X), and which lies in the appropriate quadrant so that the signs
//    of its cosine and sine match those of X and Y.
//
{
  double pi = 3.141592653589793;
//
//  Special cases:
//
  if ( x == 0.0 )
  {
    if ( 0.0 < y )
    {
      return ( pi / 2.0 );
    } 
    else if ( y < 0.0 )
    {
      return ( 3.0 * pi / 2.0 );
    } 
    else if ( y == 0.0 )
    {
      return ( 0.0 );
    }
  } 
  else if ( y == 0.0 )
  {
    if ( 0.0 < x )
    {
      return 0.0;
    } 
    else if ( x < 0.0 )
    {
      return pi;
    }
  }
//
//  We assume that ATAN2 is reliable when both arguments are positive.
//
  if        ( 0.0 < x && 0.0 < y )
  {
    return                  atan2 (  y,  x );
  }
  else if ( x < 0.0 && 0.0 < y )
  {
    return (           pi - atan2 (  y, -x ) );
  }
  else if ( x < 0.0 && y < 0.0 )
  {
    return (           pi + atan2 ( -y, -x ) );
  }
  else if ( 0.0 < x && y < 0.0 )
  {
    return ( 2.0 * pi - atan2 ( -y,  x ) );
  }

  return 0.0;
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
//    I4_MIN returns the smaller of two I4's.
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

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Input, int A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

void icos_shape ( int point_num, int edge_num, int face_num, 
  int face_order_max, double point_coord[], int edge_point[], int face_order[],
  int face_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    ICOS_SHAPE describes a icosahedron.
//
//  Discussion:
//
//    The input data required for this routine can be retrieved from
//    ICOS_SIZE.
//
//    The vertices lie on the unit sphere.
//
//    The dual of an icosahedron is the dodecahedron.
//
//    The data has been rearranged from a previous assignment.  
//    The STRIPACK program refuses to triangulate data if the first
//    three nodes are "collinear" on the sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points (12).
//
//    Input, int EDGE_NUM, the number of edges (30).
//
//    Input, int FACE_NUM, the number of faces (20).
//
//    Input, int FACE_ORDER_MAX, the maximum number of vertices 
//    per face (3).
//
//    Output, double POINT_COORD[3*POINT_NUM], the point coordinates.
//
//    Output, int EDGE_POINT[2*EDGE_NUM], the points that make up each 
//    edge, listed in ascending order of their indexes.
//
//    Output, int FACE_ORDER[FACE_NUM], the number of vertices per face.
//
//    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
//    contains the index of the I-th point in the J-th face.  The
//    points are listed in the counter clockwise direction defined
//    by the outward normal at the face.  The nodes of each face are 
//    ordered so that the lowest index occurs first.  The faces are 
//    then sorted by nodes.
//
{
# define DIM_NUM 3
# define EDGE_NUM 30
# define EDGE_ORDER 2
# define FACE_NUM 20
# define POINT_NUM 12

  double phi = 0.5 * ( sqrt ( 5.0 ) + 1.0 );

  double a = phi / sqrt ( 1.0 + phi * phi );
  double b = 1.0 / sqrt ( 1.0 + phi * phi );
  double z = 0.0;

  static int edge_point_save[EDGE_ORDER*EDGE_NUM] = {
     1,  2,
     1,  3,
     1,  4,
     1,  5,
     1,  6,
     2,  3,
     2,  4,
     2,  7,
     2,  8,
     3,  5,
     3,  7,
     3,  9,
     4,  6,
     4,  8,
     4, 10,
     5,  6,
     5,  9,
     5, 11,
     6, 10,
     6, 11,
     7,  8,
     7,  9,
     7, 12,
     8, 10,
     8, 12,
     9, 11,
     9, 12,
    10, 11,
    10, 12,
    11, 12 };
  static int face_order_save[FACE_NUM] = {
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };
  static int face_point_save[3*FACE_NUM] = {
     1,  2,  4,
     1,  3,  2,
     1,  4,  6,
     1,  5,  3,
     1,  6,  5,
     2,  3,  7,
     2,  7,  8,
     2,  8,  4,
     3,  5,  9,
     3,  9,  7,
     4,  8, 10,
     4, 10,  6,
     5,  6, 11,
     5, 11,  9,
     6, 10, 11,
     7,  9, 12,
     7, 12,  8,
     8, 12, 10,
     9, 11, 12,
    10, 12, 11 };
  int i;
  int j;
  static double point_coord_save[DIM_NUM*POINT_NUM] = {
      a,  b,  z,
      a, -b,  z,
      b,  z,  a,
      b,  z, -a,
      z,  a,  b,
      z,  a, -b,
      z, -a,  b,
      z, -a, -b,
     -b,  z,  a,
     -b,  z, -a,
     -a,  b,  z,
     -a, -b,  z };

  r8vec_copy ( DIM_NUM * point_num,       point_coord_save, point_coord );
  i4vec_copy ( EDGE_ORDER * edge_num,     edge_point_save,  edge_point );
  i4vec_copy ( face_num,                  face_order_save,  face_order );
  i4vec_copy ( face_order_max * face_num, face_point_save,  face_point );
//
//  Rebase at 0.
//
  for ( j = 0; j < edge_num; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      edge_point[i+j*2] = edge_point[i+j*2] - 1;
    }
  }

  for ( j = 0; j < face_num; j++ )
  {
    for ( i = 0; i < face_order_max; i++ )
    {
      face_point[i+j*face_order_max] = face_point[i+j*face_order_max] - 1;
    }
  }

  return;
# undef DIM_NUM
# undef EDGE_NUM
# undef EDGE_ORDER
# undef FACE_NUM
# undef POINT_NUM
}
//****************************************************************************80

void icos_size ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max )

//****************************************************************************80
//
//  Purpose:
//
//    ICOS_SIZE gives "sizes" for an icosahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *POINT_NUM, the number of points.
//
//    Output, int *EDGE_NUM, the number of edges.
//
//    Output, int *FACE_NUM, the number of faces.
//
//    Output, int *FACE_ORDER_MAX, the maximum order of any face.
//
{
  *point_num = 12;
  *edge_num = 30;
  *face_num = 20;
  *face_order_max = 3;

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
//    02 April 2005
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

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    The C MATH library includes a function GAMMA ( X ) which should be
//    invoked instead of this function.
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

double r8_uniform_01 ( int *seed )

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
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    Input, double A1[N], the vector to be copied.
//
//    Input, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
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

double r8vec_norm ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM returns the L2 norm of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double R8VEC_NORM, the L2 norm of A.
//
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
//****************************************************************************80

void r8vec_polarize ( int n, double a[], double p[], double a_normal[], 
  double a_parallel[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The (nonzero) vector P defines a direction.
//
//    The vector A can be written as the sum
//
//      A = A_normal + A_parallel
//
//    where A_parallel is a linear multiple of P, and A_normal
//    is perpendicular to P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], the vector to be polarized.
//
//    Input, double P[N], the polarizing direction.
//
//    Output, double A_NORMAL[N], A_PARALLEL[N], the normal
//    and parallel components of A.
//
{
  double a_dot_p;
  int i;
  double p_norm;

  p_norm = r8vec_norm ( n, p );

  if ( p_norm == 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      a_normal[i] = a[i];
    }
    for ( i = 0; i < n; i++ )
    {
      a_parallel[i] = 0.0;
    }
    return;
  }
  a_dot_p = r8vec_dot_product ( n, a, p ) / p_norm;

  for ( i = 0; i < n; i++ )
  {
    a_parallel[i] = a_dot_p * p[i] / p_norm;
  }

  for ( i = 0; i < n; i++ )
  {
    a_normal[i] = a[i] - a_parallel[i];
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

double sphere01_distance_xyz ( double xyz1[3], double xyz2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_DISTANCE_XYZ computes great circle distances on a unit sphere.
//
//  Discussion:
//
//    XYZ coordinates are used.
//
//    We assume the points XYZ1 and XYZ2 lie on the unit sphere.
//
//    This computation is a special form of the Vincenty formula.
//    It should be less sensitive to errors associated with very small 
//    or very large angular separations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    "Great-circle distance",
//    Wikipedia.
//
//  Parameters:
//
//    Input, double XYZ1[3], the coordinates of the first point.
//
//    Input, double XYZ2[3], the coordinates of the second point.
//
//    Output, double DIST, the great circle distance between the points.
//
{
  double bot;
  double dist;
  double lat1;
  double lat2;
  double lon1;
  double lon2;
  double r;
  double top;

  lat1 = arc_sine ( xyz1[2] );
  lon1 = atan4 ( xyz1[1], xyz1[0] );

  lat2 = arc_sine ( xyz2[2] );
  lon2 = atan4 ( xyz2[1], xyz2[0] );

  top = pow ( cos ( lat2 ) * sin ( lon1 - lon2 ), 2 )
      + pow ( cos ( lat1 ) * sin ( lat2 ) 
            - sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ), 2 );

  top = sqrt ( top );

  bot = sin ( lat1 ) * sin ( lat2 ) 
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 );

  dist = atan2 ( top, bot );

  return dist;
}
//****************************************************************************80

double sphere01_monomial_integral ( int e[3] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit sphere.
//
//  Discussion:
//
//    The integration region is 
//
//      X^2 + Y^2 + Z^2 = 1.
//
//    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
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
//    Academic Press, 1984, page 263.
//
//  Parameters:
//
//    Input, int E[3], the exponents of X, Y and Z in the 
//    monomial.  Each exponent must be nonnegative.
//
//    Output, double SPHERE01_MONOMIAL_INTEGRAL, the integral.
//
{
  int i;
  double integral;
  double pi = 3.141592653589793;

  if ( e[0] < 0 || e[1] < 0 || e[2] < 0 )
  {
    cout << "\n";
    cout << "SPHERE01_MONOMIAL_INTEGRAL - Fatal error!\n";
    cout << "  All exponents must be nonnegative.\n";
    cout << "  E[0] = " << e[0] << "\n";
    cout << "  E[1] = " << e[1] << "\n";
    cout << "  E[2] = " << e[2] << "\n";
    exit ( 1 );
  }

  if ( e[0] == 0 && e[1] == 0 && e[2] == 0 )
  {
    integral = 2.0 * sqrt ( pi * pi * pi ) / r8_gamma ( 1.5 );
  }
  else if ( ( e[0] % 2 ) == 1 ||
            ( e[1] % 2 ) == 1 ||
            ( e[2] % 2 ) == 1 )
  {
    integral = 0.0;
  }
  else
  {
    integral = 2.0;

    for ( i = 0; i < 3; i++ )
    {
      integral = integral * r8_gamma ( 0.5 * ( double ) ( e[i] + 1 ) );
    }

    integral = integral 
      / r8_gamma ( 0.5 * ( double ) ( e[0] + e[1] + e[2] + 3 ) );

  }
  return integral;
}
//****************************************************************************80

double sphere01_quad_icos1c ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_ICOS1C: centroid rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over the surface of the unit sphere.
//
//    This function sets up an icosahedral grid, and subdivides each
//    edge of the icosahedron into FACTOR subedges.  These edges define a grid
//    within each triangular icosahedral face.  The centroids of these
//    triangles can be determined.  All of these calculations are done,
//    essentially, on the FLAT faces of the icosahedron.  Only then are
//    the triangle vertices and centroids projected to the sphere.  
//
//    The resulting grid of spherical triangles and projected centroids
//    is used to apply a centroid quadrature rule over the surface of
//    the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_QUAD_ICOS1C, the estimated integral.
//
{
  int a;
  double a_xyz[3];
  double *a2_xyz;
  double area;
  double area_total;
  int b;
  double b_xyz[3];
  double *b2_xyz;
  int c;
  double c_xyz[3];
  double *c2_xyz;
  int edge_num;
  int *edge_point;
  int f1;
  int f2;
  int f3;
  int face;
  int face_num;
  int *face_order;
  int *face_point;
  int face_order_max;
  int i;
  double *node_xyz;
  double pi = 3.141592653589793;
  double *point_coord;
  int point_num;
  double result;
  double v[1];
//
//  Size the icosahedron.
//
  icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
//
//  Set the icosahedron.
//
  point_coord = new double[3*point_num];
  edge_point = new int[2*edge_num];
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];

  icos_shape ( point_num, edge_num, face_num, face_order_max, 
    point_coord, edge_point, face_order, face_point );
//
//  Initialize the integral data.
//
  result = 0.0;
  area_total = 0.0;
  *node_num = 0;
//
//  Pick a face of the icosahedron, and identify its vertices as A, B, C.
//
  for ( face = 0; face < face_num; face++ )
  {
    a = face_point[0+face*3];
    b = face_point[1+face*3];
    c = face_point[2+face*3];

    for ( i = 0; i < 3; i++ )
    {
      a_xyz[i] = point_coord[i+a*3];
      b_xyz[i] = point_coord[i+b*3];
      c_xyz[i] = point_coord[i+c*3];
    }
//
//  Some subtriangles will have the same direction as the face.
//  Generate each in turn, by determining the barycentric coordinates
//  of the centroid (F1,F2,F3), from which we can also work out the barycentric
//  coordinates of the vertices of the subtriangle.
//
    for ( f3 = 1; f3 <= 3 * factor - 2; f3 = f3 + 3 )
    {
      for ( f2 = 1; f2 <= 3 * factor - f3 - 1; f2 = f2 + 3 )
      {
        f1 = 3 * factor - f3 - f2;

        node_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3 );

        a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1 );
        b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1 );
        c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        fun ( 1, node_xyz, v );

        *node_num = *node_num + 1;
        result = result + area * v[0];
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] b2_xyz;
        delete [] c2_xyz;
        delete [] node_xyz;
      }
    }
//
//  The other subtriangles have the opposite direction from the face.
//  Generate each in turn, by determining the barycentric coordinates
//  of the centroid (F1,F2,F3), from which we can also work out the barycentric
//  coordinates of the vertices of the subtriangle.
//
    for ( f3 = 2; f3 <= 3 * factor - 4; f3 = f3 + 3 )
    {
      for ( f2 = 2; f2 <= 3 * factor - f3 - 2; f2 = f2 + 3 )
      {
        f1 = 3 * factor - f3 - f2;

        node_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3 );

        a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1 );
        b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1 );
        c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        fun ( 1, node_xyz, v );

        *node_num = *node_num + 1;
        result = result + area * v[0];
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] b2_xyz;
        delete [] c2_xyz;
        delete [] node_xyz;
      }
    }
  }
//
//  Discard allocated memory.
//
  delete [] edge_point;
  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return result;
}
//****************************************************************************80

double sphere01_quad_icos1m ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_ICOS1M: midside rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over the surface of the unit sphere.
//
//    This function sets up an icosahedral grid, and subdivides each
//    edge of the icosahedron into FACTOR subedges.  These edges define a grid
//    within each triangular icosahedral face.  The midsides of these
//    triangles can be determined.  All of these calculations are done,
//    essentially, on the FLAT faces of the icosahedron.  Only then are
//    the triangle vertices and midsides projected to the sphere.  
//
//    The resulting grid of spherical triangles and projected midsides
//    is used to apply a midside quadrature rule over the surface of
//    the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_QUAD_ICOS1M, the estimated integral.
//
{
  int a;
  double a_xyz[3];
  double *a2_xyz;
  double *a3_xyz;
  double area;
  double area_total;
  int b;
  double b_xyz[3];
  double *b2_xyz;
  double *b3_xyz;
  int c;
  double c_xyz[3];
  double *c2_xyz;
  double *c3_xyz;
  int edge;
  int edge_num;
  int *edge_point;
  int f;
  int f1;
  int f2;
  int f3;
  int face;
  int face_num;
  int *face_order;
  int *face_point;
  int face_order_max;
  int i;
  int j;
  int node;
  double node_norm;
  double *node_xyz;
  double pi = 3.141592653589793;
  double *point_coord;
  int point_num;
  double result;
  double va[1];
  double vb[1];
  double vc[1];
//
//  Size the icosahedron.
//
  icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
//
//  Set the icosahedron.
//
  point_coord = new double[3*point_num];
  edge_point = new int[2*edge_num];
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];

  icos_shape ( point_num, edge_num, face_num, face_order_max, 
    point_coord, edge_point, face_order, face_point );
//
//  Initialize the integral data.
//
  result = 0.0;
  *node_num = 0;
  area_total = 0.0;
//
//  Pick a face of the icosahedron, and identify its vertices as A, B, C.
//
  for ( face = 0; face < face_num; face++ )
  {
    a = face_point[0+face*3];
    b = face_point[1+face*3];
    c = face_point[2+face*3];

    for ( i = 0; i < 3; i++ )
    {
      a_xyz[i] = point_coord[i+a*3];
      b_xyz[i] = point_coord[i+b*3];
      c_xyz[i] = point_coord[i+c*3];
    }
//
//  Deal with subtriangles that have same orientation as face.
//
    for ( f1 = 0; f1 <= factor - 1; f1++ )
    {
      for ( f2 = 0; f2 <= factor - f1 - 1; f2++ )
      {
        f3 = factor - f1 - f2;

        a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1 );
        b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1 );
        c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        a3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2 + 1, 2 * f3 - 2 );
        b3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 + 1, 2 * f3 - 1 );
        c3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2,     2 * f3 - 1 );

        *node_num = *node_num + 3;
        fun ( 1, a3_xyz, va );
        fun ( 1, b3_xyz, vb );  
        fun ( 1, c3_xyz, vc ); 
        result = result + area * ( va[0] + vb[0] + vc[0] ) / 3.0;
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] a3_xyz;
        delete [] b2_xyz;
        delete [] b3_xyz;
        delete [] c2_xyz;
        delete [] c3_xyz;
      }
    }
//
//  Deal with subtriangles that have opposite orientation as face.
//
    for ( f3 = 0; f3 <= factor - 2; f3++ )
    {
      for ( f2 = 1; f2 <= factor - f3 - 1; f2++ )
      {
        f1 = factor - f2 - f3;

        a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1 );
        b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1 );
        c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        a3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2 - 1, 2 * f3 + 2 );
        b3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 - 1, 2 * f3 + 1 );
        c3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2,     2 * f3 + 1 );

        *node_num = *node_num + 3;
        fun ( 1, a3_xyz, va );
        fun ( 1, b3_xyz, vb );  
        fun ( 1, c3_xyz, vc );
        result = result + area * ( va[0] + vb[0] + vc[0] ) / 3.0;
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] a3_xyz;
        delete [] b2_xyz;
        delete [] b3_xyz;
        delete [] c2_xyz;
        delete [] c3_xyz;
      }
    }
  }
//
//  Discard allocated memory.
//
  delete [] edge_point;
  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return result;
}
//****************************************************************************80

double sphere01_quad_icos1v ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_ICOS1V: vertex rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over the surface of the unit sphere.
//
//    This function sets up an icosahedral grid, and subdivides each
//    edge of the icosahedron into FACTOR subedges.  These edges define a grid
//    within each triangular icosahedral face.  The vertices of these
//    triangles can be determined.  All of these calculations are done,
//    essentially, on the FLAT faces of the icosahedron.  Only then are
//    the triangle vertices projected to the sphere.  
//
//    The resulting grid of spherical triangles is used to apply a vertex
//    quadrature rule over the surface of the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_QUAD_ICOS2M, the estimated integral.
//
{
  int a;
  double a_xyz[3];
  double *a2_xyz;
  double area;
  double area_total;
  int b;
  double b_xyz[3];
  double *b2_xyz;
  int c;
  double c_xyz[3];
  double *c2_xyz;
  int edge;
  int edge_num;
  int *edge_point;
  int f;
  int f1;
  int f2;
  int f3;
  int face;
  int face_num;
  int *face_order;
  int *face_point;
  int face_order_max;
  int i;
  int j;
  int node;
  double pi = 3.141592653589793;
  double *point_coord;
  int point_num;
  double result;
  double va[1];
  double vb[1];
  double vc[1];
//
//  Size the icosahedron.
//
  icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
//
//  Set the icosahedron.
//
  point_coord = new double[3*point_num];
  edge_point = new int[2*edge_num];
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];

  icos_shape ( point_num, edge_num, face_num, face_order_max, 
    point_coord, edge_point, face_order, face_point );
//
//  Initialize the integral data.
//
  result = 0.0;
  *node_num = 0;
  area_total = 0.0;
//
//  Pick a face of the icosahedron, and identify its vertices as A, B, C.
//
  for ( face = 0; face < face_num; face++ )
  {
    a = face_point[0+face*3];
    b = face_point[1+face*3];
    c = face_point[2+face*3];

    for ( i = 0; i < 3; i++ )
    {
      a_xyz[i] = point_coord[i+a*3];
      b_xyz[i] = point_coord[i+b*3];
      c_xyz[i] = point_coord[i+c*3];
    }
//
//  Deal with subtriangles that have same orientation as face.
//
    for ( f1 = 0; f1 <= factor - 1; f1++ )
    {
      for ( f2 = 0; f2 <= factor - f1 - 1; f2++ )
      {
        f3 = factor - f1 - f2;

        a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1 );
        b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1 );
        c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        *node_num = *node_num + 3;
        fun ( 1, a2_xyz, va );
        fun ( 1, b2_xyz, vb );  
        fun ( 1, c2_xyz, vc ); 
        result = result + area * ( va[0] + vb[0] + vc[0] ) / 3.0;
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] b2_xyz;
        delete [] c2_xyz;
      }
    }
//
//  Deal with subtriangles that have opposite orientation as face.
//
    for ( f3 = 0; f3 <= factor - 2; f3++ )
    {
      for ( f2 = 1; f2 <= factor - f3 - 1; f2++ )
      {
        f1 = factor - f2 - f3;

        a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1 );
        b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1 );
        c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        *node_num = *node_num + 3;
        fun ( 1, a2_xyz, va );
        fun ( 1, b2_xyz, vb );  
        fun ( 1, c2_xyz, vc ); 
        result = result + area * ( va[0] + vb[0] + vc[0] ) / 3.0;
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] b2_xyz;
        delete [] c2_xyz;
      }
    }
  }
//
//  Discard allocated memory.
//
  delete [] edge_point;
  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return result;
}
//****************************************************************************80

double sphere01_quad_icos2v ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_ICOS2V: vertex rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over the surface of the unit sphere.
//
//    This function sets up an icosahedral grid, and subdivides each
//    edge of the icosahedron into FACTOR subedges.  These edges define a grid
//    within each triangular icosahedral face.  The vertices of these
//    triangles can be determined.  All of these calculations are done,
//    essentially, on the FLAT faces of the icosahedron.  Only then are
//    the triangle vertices projected to the sphere.  
//
//    The resulting grid of spherical triangles is used to apply a vertex
//    quadrature rule over the surface of the unit sphere.
//
//    This is a revision of SPHERE01_QUAD_ICOS2V that attempted to use a more
//    sophisticated scheme to map points from the planar triangle to the surface
//    of the unit sphere.  Very little improvement to the estimated integral
//    was observed, so development of this scheme has been set aside for now.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, void FUN ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_QUAD_ICOS2V, the estimated integral.
//
{
  int a;
  double a_xyz[3];
  double *a2_xyz;
  double area;
  double area_total;
  int b;
  double b_xyz[3];
  double *b2_xyz;
  int c;
  double c_xyz[3];
  double *c2_xyz;
  int edge;
  int edge_num;
  int *edge_point;
  int f;
  int f1;
  int f2;
  int f3;
  int face;
  int face_num;
  int *face_order;
  int *face_point;
  int face_order_max;
  int i;
  int j;
  int node;
  double pi = 3.141592653589793;
  double *point_coord;
  int point_num;
  double result;
  double va[1];
  double vb[1];
  double vc[1];
//
//  Size the icosahedron.
//
  icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
//
//  Set the icosahedron.
//
  point_coord = new double[3*point_num];
  edge_point = new int[2*edge_num];
  face_order = new int[face_num];
  face_point = new int[face_order_max*face_num];

  icos_shape ( point_num, edge_num, face_num, face_order_max, 
    point_coord, edge_point, face_order, face_point );
//
//  Initialize the integral data.
//
  result = 0.0;
  *node_num = 0;
  area_total = 0.0;
//
//  Pick a face of the icosahedron, and identify its vertices as A, B, C.
//
  for ( face = 0; face < face_num; face++ )
  {
    a = face_point[0+face*3];
    b = face_point[1+face*3];
    c = face_point[2+face*3];

    for ( i = 0; i < 3; i++ )
    {
      a_xyz[i] = point_coord[i+a*3];
      b_xyz[i] = point_coord[i+b*3];
      c_xyz[i] = point_coord[i+c*3];
    }
//
//  Deal with subtriangles that have same orientation as face.
//
    for ( f1 = 0; f1 <= factor - 1; f1++ )
    {
      for ( f2 = 0; f2 <= factor - f1 - 1; f2++ )
      {
        f3 = factor - f1 - f2;

        a2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1 );
        b2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1 );
        c2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        *node_num = *node_num + 3;
        fun ( 1, a2_xyz, va );
        fun ( 1, b2_xyz, vb );  
        fun ( 1, c2_xyz, vc ); 
        result = result + area * ( va[0] + vb[0] + vc[0] ) / 3.0;
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] b2_xyz;
        delete [] c2_xyz;
      }
    }
//
//  Deal with subtriangles that have opposite orientation as face.
//
    for ( f3 = 0; f3 <= factor - 2; f3++ )
    {
      for ( f2 = 1; f2 <= factor - f3 - 1; f2++ )
      {
        f1 = factor - f2 - f3;

        a2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1 );
        b2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1 );
        c2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

        area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

        *node_num = *node_num + 3;
        fun ( 1, a2_xyz, va );
        fun ( 1, b2_xyz, vb );  
        fun ( 1, c2_xyz, vc ); 
        result = result + area * ( va[0] + vb[0] + vc[0] ) / 3.0;
        area_total = area_total + area;

        delete [] a2_xyz;
        delete [] b2_xyz;
        delete [] c2_xyz;
      }
    }
  }
//
//  Discard allocated memory.
//
  delete [] edge_point;
  delete [] face_order;
  delete [] face_point;
  delete [] point_coord;

  return result;
}
//****************************************************************************80

double sphere01_quad_llc ( void f ( int n, double x[], double v[] ), double h, 
  int *n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_LLC: Longitude/Latitude grid with centroid rule.
//
//  Discussion:
//
//    The sphere is broken up into spherical triangles, whose sides
//    do not exceed the length H.  Then a centroid rule is used on
//    each spherical triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, void F ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Input, double H, the maximum length of a side of the spherical
//    quadrilaterals.
//
//    Output, int *N, the number of points used.
//
//    Output, double SPHERE01_QUAD_LLC, the approximate integral.
//
{
  double area;
  int i;
  int j;
  double phi;
  int phi_num;
  double phi1;
  double phi2;
  double pi = 3.141592653589793;
  double result;
  double sector_area;
  double sphere_area;
  double theta;
  int theta_num;
  double theta1;
  double theta2;
  double v[1];
  double *x;
  double *x1;
  double *x11;
  double *x12;
  double *x2;
  double *x21;
  double *x22;
//
//  Choose PHI and THETA counts that make short sides.
//
  phi_num = ( int ) ( pi / h );

  if ( h * ( double ) ( phi_num ) < pi )
  {
    phi_num = phi_num + 1;
  }

  theta_num = ( int ) ( 2.0 * pi / h );

  if ( h * ( double ) ( theta_num ) < pi )
  {
    theta_num = theta_num + 1;
  }

  *n = 0;
  result = 0.0;
//
//  Only one THETA (and hence, only one PHI.)
//
  if ( theta_num == 1 )
  {
    sphere_area = 4.0 * pi;

    theta = 0.0;
    phi = pi / 2.0;
    x = tp_to_xyz ( theta, phi );

    f ( 1, x, v );
    *n = *n + 1;
    result = sphere_area * v[0];
    delete [] x;
  }
//
//  Several THETA''s, one PHI.
//
  else if ( phi_num == 1 )
  {
    sphere_area = 4.0 * pi;
    sector_area = sphere_area / ( double ) ( theta_num );

    result = 0.0;

    for ( j = 1; j <= theta_num; j++ )
    {
      theta = ( double ) ( ( j - 1 ) * 2 ) * pi / ( double ) ( theta_num );
      phi = pi / 2.0;
      x = tp_to_xyz ( theta, phi );
      f ( 1, x, v );
      *n = *n + 1;
      result = result + sector_area * v[0];
      delete [] x;
    }
  }
//
//  At least two PHI''s.
//
  else
  {
    result = 0.0;
//
//  Picture in top row, with V1 = north pole:
//
//        V1
//       /  \
//      /    \
//    V12----V22
//
    phi1 = 0.0;
    phi2 = pi / ( double ) ( phi_num );

    for ( j = 1; j <= theta_num; j++ )
    {
      theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
      theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

      x1 = tp_to_xyz ( theta1, phi1 );
      x12 = tp_to_xyz ( theta1, phi2 );
      x22 = tp_to_xyz ( theta2, phi2 );

      area = sphere01_triangle_vertices_to_area ( x1, x12, x22 );
      x = sphere01_triangle_vertices_to_centroid ( x1, x12, x22 );
      f ( 1, x, v );
      *n = *n + 1;
      result = result + area * v[0];
      delete [] x;
      delete [] x1;
      delete [] x12;
      delete [] x22;
    }
//
//  Picture in all intermediate rows:
//
//    V11--V21
//     | \  |
//     |  \ |
//    V12--V22
//
    for ( i = 2; i <= phi_num - 1; i++ )
    {
      phi1 = ( double ) ( i - 1 ) * pi / ( double ) ( phi_num );
      phi2 = ( double ) ( i ) * pi / ( double ) ( phi_num );

      for ( j = 1; j <= theta_num; j++ )
      {
        theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
        theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

        x11 = tp_to_xyz ( theta1, phi1 );
        x21 = tp_to_xyz ( theta2, phi1 );
        x12 = tp_to_xyz ( theta1, phi2 );
        x22 = tp_to_xyz ( theta2, phi2 );

        area = sphere01_triangle_vertices_to_area ( x11, x12, x22 );
        x = sphere01_triangle_vertices_to_centroid ( x11, x12, x22 );
        f ( 1, x, v );
        *n = *n + 1;
        result = result + area * v[0];
        delete [] x;

        area = sphere01_triangle_vertices_to_area ( x22, x21, x11 );
        x = sphere01_triangle_vertices_to_centroid ( x22, x21, x11 );
        f ( 1, x, v );
        *n = *n + 1;
        result = result + area * v[0];
        delete [] x;
        delete [] x11;
        delete [] x12;
        delete [] x21;
        delete [] x22;
      }
    }
//
//  Picture in last row, with V2 = south pole:
//
//    V11----V21
//      \    /
//       \  /
//        V2
//
    phi1 = ( double ) ( phi_num - 1 ) * pi / ( double ) ( phi_num );
    phi2 = pi;

    for ( j = 1; j <= theta_num; j++ )
    {
      theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
      theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

      x11 = tp_to_xyz ( theta1, phi1 );
      x21 = tp_to_xyz ( theta2, phi1 );
      x2 = tp_to_xyz ( theta2, phi2 );

      area = sphere01_triangle_vertices_to_area ( x11, x2, x21 );
      x = sphere01_triangle_vertices_to_centroid ( x11, x2, x21 );
      f ( 1, x, v );
      *n = *n + 1;
      result = result + area * v[0];
      delete [] x;
      delete [] x11;
      delete [] x2;
      delete [] x21;
    }
  }
  return result;
}
//****************************************************************************80

double sphere01_quad_llm ( void f ( int n, double x[], double v[] ), double h, 
  int *n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_LLM: longitude/latitude grid plus midside rule.
//
//  Discussion:
//
//    The sphere is broken up into spherical triangles, whose sides
//    do not exceed the length H.  Then the function is evaluated
//    at the midsides, and the average is multiplied by the area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, void F ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Input, double H, the maximum length of a side of the spherical
//    quadrilaterals.
//
//    Output, int *N, the number of points used.
//
//    Output, double SPHERE01_QUAD_LLM, the approximate integral.
//
{
  double area;
  int i;
  int j;
  double m1[3];
  double m2[3];
  double m3[3];
  double phi;
  int phi_num;
  double phi1;
  double phi2;
  double pi = 3.141592653589793;
  double result;
  double sector_area;
  double sphere_area;
  double theta;
  int theta_num;
  double theta1;
  double theta2;
  double v[1];
  double *x;
  double *x1;
  double *x11;
  double *x12;
  double *x2;
  double *x21;
  double *x22;
//
//  Choose PHI and THETA counts that make short sides.
//
  phi_num = ( int ) ( pi / h );

  if ( h * ( double ) ( phi_num ) < pi )
  {
    phi_num = phi_num + 1;
  }

  theta_num = ( int ) ( 2.0 * pi / h );

  if ( h * ( double ) ( theta_num ) < pi )
  {
    theta_num = theta_num + 1;
  }

  *n = 0;
  result = 0.0;
//
//  Only one THETA (and hence, only one PHI.)
//
  if ( theta_num == 1 )
  {
    sphere_area = 4.0 * pi;

    theta = 0.0;
    phi = pi / 2.0;
    x = tp_to_xyz ( theta, phi );
    f ( 1, x, v );
    *n = *n + 1;
    result = sphere_area * v[0];
    delete [] x;
  }
//
//  Several THETA''s, one PHI.
//
  else if ( phi_num == 1 )
  {
    sphere_area = 4.0 * pi;
    sector_area = sphere_area / ( double ) ( theta_num );

    result = 0.0;

    for ( j = 1; j <= theta_num; j++ )
    {
      theta = ( double ) ( ( j - 1 ) * 2 ) * pi / ( double ) ( theta_num );
      phi = pi / 2.0;
      x = tp_to_xyz ( theta, phi );
      f ( 1, x, v );
      *n = *n + 1;
      result = result + sector_area * v[0];
      delete [] x;
    }
  }
//
//  At least two PHI''s.
//
  else
  {
    result = 0.0;
//
//  Picture:
//
//        V1
//       /  \
//      /    \
//    V12----V22
//
    phi1 = 0.0;
    phi2 = pi / ( double ) ( phi_num );

    for ( j = 1; j <= theta_num; j++ )
    {
      theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
      theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

      x1 = tp_to_xyz ( theta1, phi1 );
      x12 = tp_to_xyz ( theta1, phi2 );
      x22 = tp_to_xyz ( theta2, phi2 );

      area = sphere01_triangle_vertices_to_area ( x1, x12, x22 );

      sphere01_triangle_vertices_to_midpoints ( x1, x12, x22, m1, m2, m3 );

      f ( 1, m1, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, m2, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, m3, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      delete [] x1;
      delete [] x12;
      delete [] x22;
    }
//
//  Picture:
//
//    V11--V21
//     | \  |
//     |  \ |
//    V12--V22
//
    for ( i = 2; i <= phi_num - 1; i++ )
    {
      phi1 = ( double ) ( i - 1 ) * pi / ( double ) ( phi_num );
      phi2 = ( double ) ( i ) * pi / ( double ) ( phi_num );

      for ( j = 1; j <= theta_num; j++ )
      {
        theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
        theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

        x11 = tp_to_xyz ( theta1, phi1 );
        x21 = tp_to_xyz ( theta2, phi1 );
        x12 = tp_to_xyz ( theta1, phi2 );
        x22 = tp_to_xyz ( theta2, phi2 );

        area = sphere01_triangle_vertices_to_area ( x11, x12, x22 );

        sphere01_triangle_vertices_to_midpoints ( x11, x12, x22, m1, m2, m3 );

        f ( 1, m1, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, m2, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, m3, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;

        area = sphere01_triangle_vertices_to_area ( x22, x21, x11 );

        sphere01_triangle_vertices_to_midpoints ( x22, x21, x11, m1, m2, m3 );

        f ( 1, m1, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, m2, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, m3, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        delete [] x11;
        delete [] x12;
        delete [] x21;
        delete [] x22;
      }
    }
//
//  Picture:
//
//    V11----V21
//      \    /
//       \  /
//        V2
//
    phi1 = ( double ) ( phi_num - 1 ) * pi / ( double ) ( phi_num );
    phi2 = pi;

    for ( j = 1; j <= theta_num; j++ )
    {
      theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
      theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

      x11 = tp_to_xyz ( theta1, phi1 );
      x21 = tp_to_xyz ( theta2, phi1 );
      x2 = tp_to_xyz ( theta2, phi2 );

      area = sphere01_triangle_vertices_to_area ( x11, x2, x21 );

      sphere01_triangle_vertices_to_midpoints ( x11, x2, x21, m1, m2, m3 );

      f ( 1, m1, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, m2, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, m3, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      delete [] x11;
      delete [] x2;
      delete [] x21;
    }
  }

  return result;
}
//****************************************************************************80

double sphere01_quad_llv ( void f ( int n, double x[], double v[] ), double h, 
  int *n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_LLV: longitude/latitude grid with vertex rule.
//
//  Discussion:
//
//    The sphere is broken up into spherical triangles, whose sides
//    do not exceed the length H.  Then the function is evaluated
//    at the vertices, and the average is multiplied by the area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, void F ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Input, double H, the maximum length of a side of the spherical
//    quadrilaterals.
//
//    Output, int *N, the number of points used.
//
//    Output, double SPHERE01_QUAD_LLV, the approximate integral.
//
{
  double area;
  int i;
  int j;
  double phi;
  int phi_num;
  double phi1;
  double phi2;
  double pi = 3.141592653589793;
  double result;
  double sector_area;
  double sphere_area;
  double theta;
  int theta_num;
  double theta1;
  double theta2;
  double v[1];
  double *x;
  double *x1;
  double *x11;
  double *x12;
  double *x2;
  double *x21;
  double *x22;
//
//  Choose PHI and THETA counts that make short sides.
//
  phi_num = ( int ) ( pi / h );

  if ( h * ( double ) ( phi_num ) < pi )
  {
    phi_num = phi_num + 1;
  }

  theta_num = ( int ) ( 2.0 * pi / h );

  if ( h * ( double ) ( theta_num ) < pi )
  {
    theta_num = theta_num + 1;
  }

  *n = 0;
  result = 0.0;
//
//  Only one THETA (and hence, only one PHI.)
//
  if ( theta_num == 1 )
  {
    sphere_area = 4.0 * pi;

    theta = 0.0;
    phi = pi / 2.0;
    x = tp_to_xyz ( theta, phi );
    f ( 1, x, v );
    result = sphere_area * v[0];
    delete [] x;
  }
//
//  Several THETA''s, one PHI.
//
  else if ( phi_num == 1 )
  {
    sphere_area = 4.0 * pi;
    sector_area = sphere_area / ( double ) ( theta_num );

    result = 0.0;

    for ( j = 1; j <= theta_num; j++ )
    {
      theta = ( double ) ( ( j - 1 ) * 2 ) * pi / ( double ) ( theta_num );
      phi = pi / 2.0;
      x = tp_to_xyz ( theta, phi );
      f ( 1, x, v );
      *n = *n + 1;
      result = result + sector_area * v[0];
      delete [] x;
    }
  }
//
//  At least two PHI''s.
//
  else
  {
    result = 0.0;
//
//  Picture:
//
//        V1
//       /  \
//      /    \
//    V12----V22
//
    phi1 = 0.0;
    phi2 = pi / ( double ) ( phi_num );

    for ( j = 1; j <= theta_num; j++ )
    {
      theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
      theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

      x1 = tp_to_xyz ( theta1, phi1 );
      x12 = tp_to_xyz ( theta1, phi2 );
      x22 = tp_to_xyz ( theta2, phi2 );

      area = sphere01_triangle_vertices_to_area ( x1, x12, x22 );

      f ( 1, x1, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, x12, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, x22, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      delete [] x1;
      delete [] x12;
      delete [] x22;
    }
//
//  Picture:
//
//    V11--V21
//     | \  |
//     |  \ |
//    V12--V22
//
    for ( i = 2; i <= phi_num - 1; i++ )
    {
      phi1 = ( double ) ( i - 1 ) * pi / ( double ) ( phi_num );
      phi2 = ( double ) ( i ) * pi / ( double ) ( phi_num );

      for ( j = 1; j <= theta_num; j++ )
      {
        theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
        theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

        x11 = tp_to_xyz ( theta1, phi1 );
        x21 = tp_to_xyz ( theta2, phi1 );
        x12 = tp_to_xyz ( theta1, phi2 );
        x22 = tp_to_xyz ( theta2, phi2 );

        area = sphere01_triangle_vertices_to_area ( x11, x12, x22 );

        f ( 1, x11, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, x12, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, x22, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;

        area = sphere01_triangle_vertices_to_area ( x22, x21, x11 );

        f ( 1, x22, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, x21, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;
        f ( 1, x11, v );
        *n = *n + 1;
        result = result + area * v[0] / 3.0;

        delete [] x11;
        delete [] x12;
        delete [] x21;
        delete [] x22;
      }
    }
//
//  Picture:
//
//    V11----V21
//      \    /
//       \  /
//        V2
//
    phi1 = ( double ) ( phi_num - 1 ) * pi / ( double ) ( phi_num );
    phi2 = pi;

    for ( j = 1; j <= theta_num; j++ )
    {
      theta1 = ( double ) ( j - 1 ) * 2.0 * pi / ( double ) ( theta_num );
      theta2 = ( double ) ( j ) * 2.0 * pi / ( double ) ( theta_num );

      x11 = tp_to_xyz ( theta1, phi1 );
      x21 = tp_to_xyz ( theta2, phi1 );
      x2 = tp_to_xyz ( theta2, phi2 );

      area = sphere01_triangle_vertices_to_area ( x11, x2, x21 );

      f ( 1, x11, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, x2, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;
      f ( 1, x21, v );
      *n = *n + 1;
      result = result + area * v[0] / 3.0;

      delete [] x11;
      delete [] x2;
      delete [] x21;
    }
  }
  return result;
}
//****************************************************************************80

double sphere01_quad_mc ( void f ( int n, double x[], double v[] ), double h, 
  int *seed, int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_MC uses the Monte Carlo rule for sphere quadrature.
//
//  Discussion:
//
//    A number of points N are chosen at random on the sphere, with N
//    being determined so that, if the points were laid out on a regular
//    grid, the average spacing would be no more than H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, void F ( int n, double x[], double v[] ), evaluates the 
//    integrand.
//
//    Input, double H, the maximum length of a side of the spherical
//    quadrilaterals.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Input, int N, the number of points used.
//
//    Output, double SPHERE01_QUAD_MC, the approximate integral.
//
{
  double pi = 3.141592653589793;
  double result;
  double sphere_area;
  double *v;
  double *x;

  sphere_area = 4.0 * pi;

  x = sphere01_sample ( n, seed );

  v = new double[n];

  f ( n, x, v );

  result = sphere_area * r8vec_sum ( n, v ) / ( double ) ( n );

  delete [] v;
  delete [] x;

  return result;
}
//****************************************************************************80

int sphere01_quad_mc_size ( double h )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_QUAD_MC_SIZE sizes a Monte Carlo rule for sphere quadrature.
//
//  Discussion:
//
//    A number of points N are chosen at random on the sphere, with N
//    being determined so that, if the points were laid out on a regular
//    grid, the average spacing would be no more than H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double H, the maximum length of a side of the spherical
//    quadrilaterals.
//
//    Output, int SPHERE01_QUAD_MC_SIZE, the number of points to use.
//
{
  int n;
  double pi = 3.141592653589793;
  double sphere_area;
//
//  The sphere's area is 4 * PI.
//  Choose N so that we divide this area into N subareas of PI * H * H.
//
  sphere_area = 4.0 * pi;

  n = ( int ) ( sphere_area / h / h );
  n = i4_max ( n, 1 );

  return n;
}
//****************************************************************************80

double *sphere01_sample ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_SAMPLE picks random points on a unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of samples.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double SPHERE01_SAMPLE[3*N], the sample points.
//
{
  int j;
  double phi;
  double pi = 3.141592653589793;
  double theta;
  double vdot;
  double *x;

  x = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Pick a uniformly random VDOT, which must be between -1 and 1.
//  This represents the dot product of the random vector with the Z unit vector.
//
//  Note: this works because the surface area of the sphere between
//  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
//  a patch of area uniformly.
//
    vdot = r8_uniform_01 ( seed );
    vdot = 2.0 * vdot - 1.0;

    phi = acos ( vdot );
//
//  Pick a uniformly random rotation between 0 and 2 Pi around the
//  axis of the Z vector.
//
    theta = r8_uniform_01 ( seed );
    theta = 2.0 * pi * theta;

    x[0+j*3] = cos ( theta ) * sin ( phi );
    x[1+j*3] = sin ( theta ) * sin ( phi );
    x[2+j*3] =                 cos ( phi );
  }

  return x;
}
//****************************************************************************80

double sphere01_triangle_angles_to_area ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_ANGLES_TO_AREA: area of a spherical triangle on the unit sphere.
//
//  Discussion:
//
//    A unit sphere centered at 0 in 3D satisfies the equation:
//
//      X^2 + Y^2 + Z^2 = 1
//
//    A spherical triangle is specified by three points on the surface
//    of the sphere.
//
//    The area formula is known as Girard's formula.
//
//    The area of a spherical triangle on the unit sphere is:
//
//      AREA = ( A + B + C - PI )
//
//    where A, B and C are the (surface) angles of the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the angles of the triangle.
//
//    Output, double STRI_ANGLES_TO_AREA, the area of the spherical triangle.
//
{
  double area;
  double pi = 3.141592653589793;

  area = a + b + c - pi;

  return area;
}
//****************************************************************************80

double *sphere01_triangle_project ( double a_xyz[3], double b_xyz[3], double c_xyz[3], 
  int f1, int f2, int f3 )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_PROJECT projects from a plane triangle to a spherical triangle.
//
//  Discussion:
//
//    We assume that points A, B and C lie on the unit sphere, and they
//    thus define a spherical triangle.
//
//    They also, of course, define a planar triangle.
//
//    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
//    planar triangle.
//
//    This function determines the coordinates of the point in the planar
//    triangle identified by the barycentric coordinates, and returns the
//    coordinates of the projection of that point onto the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the coordinates
//    of the points A, B, and C.
//
//    Input, int F1, F2, F3, the barycentric coordinates
//    of a point in the triangle ABC.  Normally, these coordinates would
//    be real numbers, and would sum to 1.  For convenience, we allow these
//    to be integers which must be divided by F1+F2+F3.
//
//    Output, double NODE_XYZ[3], the coordinates of the 
//    point on the unit sphere which is the projection of the point on the plane
//    whose barycentric coordinates with respect to A, B, and C is
//    (F1,F2,F3)/(F1+F2+F3).
//
{
  int i;
  double *node_xyz;
  double norm;

  node_xyz = new double[3];

  for ( i = 0; i < 3; i++ )
  {
    node_xyz[i] = 
      ( ( double ) ( f1           ) * a_xyz[i]   
      + ( double ) (      f2      ) * b_xyz[i]   
      + ( double ) (           f3 ) * c_xyz[i] ) 
      / ( double ) ( f1 + f2 + f3 );
  }
  norm = r8vec_norm ( 3, node_xyz );

  for ( i = 0; i < 3; i++ )
  {
    node_xyz[i] = node_xyz[i] / norm;
  }
  return node_xyz;
}
//****************************************************************************80

double *sphere01_triangle_project2 ( double a_xyz[3], double b_xyz[3], double c_xyz[3], 
  int f1, int f2, int f3 )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_PROJECT2 projects from a plane triangle to a spherical triangle.
//
//  Discussion:
//
//    We assume that points A, B and C lie on the unit sphere, and they
//    thus define a spherical triangle.
//
//    They also, of course, define a planar triangle.
//
//    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
//    planar triangle.
//
//    This function determines the coordinates of the point in the planar
//    triangle identified by the barycentric coordinates, and returns the
//    coordinates of the projection of that point onto the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
//    of the points A, B, and C.
//
//    Input, int F1, F2, F3, the barycentric coordinates
//    of a point in the triangle ABC.  Normally, these coordinates would
//    be real numbers, and would sum to 1.  For convenience, we allow these
//    to be integers which must be divided by F1+F2+F3.
//
//    Output, double SPHERE01_TRIANGLE_PROJECT2[3], the coordinates of the 
//    point on the unit sphere which is the projection of the point on the 
//    plane whose barycentric coordinates with respect to A, B, and C is
//    (F1,F2,F3)/(F1+F2+F3).
//
{
  double ab[3];
  double ac[3];
  double acn[3];
  double acp[3];
  double angle;
  double bn[3];
  double bp[3];
  double cn[3];
  double cp[3];
  int i;
  double *node_xyz;
  double norm;
  double theta_ab;
  double theta_ac;
  double theta_bc;

  node_xyz = new double[3];
//
//  This check avoids 0/0 calculations later.
//
  if ( f2 == 0 && f3 == 0 )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_xyz[i] = a_xyz[i];
    }
    return node_xyz;
  }
  else if ( f1 == 0 && f3 == 0 )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_xyz[i] = b_xyz[i];
    }
    return node_xyz;
  }
  else if ( f1 == 0 && f2 == 0 )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_xyz[i] = c_xyz[i];
    }
    return node_xyz;
  }
//
//  Determine the angular distances (A,B) and (A,C).
//
  theta_ab = sphere01_distance_xyz ( a_xyz, b_xyz );

  theta_ac = sphere01_distance_xyz ( a_xyz, c_xyz );
//
//  Polarize B = BP + BN
//  Normalize BN, 
//  Same for C.
//
  r8vec_polarize ( 3, b_xyz, a_xyz, bn, bp );
  norm = r8vec_norm ( 3, bn );
  for ( i = 0; i < 3; i++ )
  {
    bn[i] = bn[i] / norm;
  }
  r8vec_polarize ( 3, c_xyz, a_xyz, cn, cp );
  norm = r8vec_norm ( 3, cn );
  for ( i = 0; i < 3; i++ )
  {
    cn[i] = cn[i] / norm;
  }
//
//  Determine AB and AC that use cos ( ( F2 + F3 ) / ( F1 + F2 + F3 ) ) of A
//  and cos ( F1 / ( F1 + F2 + F3 ) ) of B or C.
//
  angle = ( ( double ) ( f2 + f3 ) * theta_ab ) / ( double ) ( f1 + f2 + f3 );
  for ( i = 0; i < 3; i++ )
  {
    ab[i] = cos ( angle ) * a_xyz[i] + sin ( angle ) * bn[i];
  }
  angle = ( ( double ) ( f2 + f3 ) * theta_ac ) / ( double ) ( f1 + f2 + f3 );
  for ( i = 0; i < 3; i++ )
  {
    ac[i] = cos ( angle ) * a_xyz[i] + sin ( angle ) * cn[i];
  }
//
//  Determine the angular distance between AB and AC.
//
  theta_bc = sphere01_distance_xyz ( ab, ac );
//
//  Polarize AC = ACP + ACN, normalize ACN.
//
  r8vec_polarize ( 3, ac, ab, acn, acp );
  norm = r8vec_norm ( 3, acn );
  for ( i = 0; i < 3; i++ )
  {
    acn[i] = acn[i] / norm;
  }
//
//  The interval between AB and AC is marked by F2+F3+1 vertices 0 through F2+F3.
//
  angle = ( ( double ) ( f3 ) * theta_bc ) / ( double ) ( f2 + f3 );

  for ( i = 0; i < 3; i++ )
  {
    node_xyz[i] = cos ( angle ) * ab[i] + sin ( angle ) * acn[i];
  }
  return node_xyz;
}
//****************************************************************************80

double *sphere01_triangle_sample ( int n, double v1[3], double v2[3], double v3[3], 
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_SAMPLE: sample points from triangle on unit sphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    A spherical triangle on the surface of the unit sphere contains those 
//    points with radius R = 1, bounded by the vertices V1, V2, V3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Arvo,
//    Stratified sampling of spherical triangles,
//    Computer Graphics Proceedings, Annual Conference Series, 
//    ACM SIGGRAPH '95, pages 437-438, 1995.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double V1[3], V2[3], V3[3], the XYZ coordinates of
//    the vertices of the spherical triangle.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double SPHERE01_TRIANGLE_SAMPLE[3*N], the XYZ coordinates of the 
//    sample points.
//
{
  double a;
  double alpha;
  double area;
  double area_hat;
  double b;
  double beta;
  double c;
  double gamma;
  int i;
  int j;
  double norm;
  double q;
  double s;
  double t;
  double u;
  double v;
  double v31[3];
  double v4[3];
  double v42[3];
  double w;
  double *x;
  double xsi1;
  double xsi2;
  double z;

  sphere01_triangle_vertices_to_sides ( v1, v2, v3, &a, &b, &c );

  sphere01_triangle_sides_to_angles ( a, b, c, &alpha, &beta, &gamma );

  area = sphere01_triangle_angles_to_area ( alpha, beta, gamma );

  x = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Select the new area.
//
    xsi1 = r8_uniform_01 ( seed );
    area_hat = xsi1 * area;
//
//  Compute the sine and cosine of the angle phi.
//
    s = sin ( area_hat - alpha );
    t = cos ( area_hat - alpha );
//
//  Compute the pair that determines beta_hat.
//
    u = t - cos ( alpha );
    v = s + sin ( alpha ) * cos ( c );
//
//  Q is the cosine of the new edge length b_hat.
//
    q = ( ( v * t - u * s ) * cos ( alpha ) - v ) 
      / ( ( v * s + u * t ) * sin ( alpha ) );
//
//  We very occasionally get a Q value out of bounds.
//
    q = r8_max ( q, - 1.0 );
    q = r8_min ( q, + 1.0 );
//
//  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
//
    w = r8vec_dot_product ( 3, v3, v1 );
    for ( i = 0; i < 3; i++ )
    {
      v31[i] = v3[i] - w * v1[i];
    }
    norm = r8vec_norm ( 3, v31 );
    for ( i = 0; i < 3; i++ )
    {
      v31[i] = v31[i] / norm;
    }
//
//  V4 is the third vertex of the subtriangle V1, V2, V4.
//
    for ( i = 0; i < 3; i++ )
    {
      v4[i] = q * v1[i] + sqrt ( 1.0 - q * q ) * v31[i];
    }
//
//  Select cos theta, which will sample along the edge from V2 to V4.
//
    xsi2 = r8_uniform_01 ( seed );
    z = 1.0 - xsi2 * ( 1.0 - r8vec_dot_product ( 3, v4, v2 ) );
//
//  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
//
    w = r8vec_dot_product ( 3, v4, v2 );
    for ( i = 0; i < 3; i++ )
    {
      v42[i] = v4[i] - w * v2[i];
    }
    norm = r8vec_norm ( 3, v42 );
    for ( i = 0; i < 3; i++ )
    {
      v42[i] = v42[i] / norm;
    }
//
//  Construct the point.
//
    for ( i = 0; i < 3; i++ )
    {
      x[i+j*3] = z * v2[i] + sqrt ( 1.0 - z * z ) * v42[i];
    }
  }
  return x;
}
//****************************************************************************80

void sphere01_triangle_sides_to_angles ( double as, double bs, double cs, 
  double *a, double *b, double *c )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_SIDES_TO_ANGLES: angles of spherical triangle on unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double AS, BS, CS, the (geodesic) length of the sides of the
//    triangle.
//
//    Output, double *A, *B, *C, the spherical angles of the triangle.
//    Angle A is opposite the side of length AS, and so on.
//
{
  double asu;
  double bsu;
  double csu;
  double ssu;
  double tan_a2;
  double tan_b2;
  double tan_c2;

  asu = as;
  bsu = bs;
  csu = cs;
  ssu = ( asu + bsu + csu ) / 2.0;

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / 
                  ( sin ( ssu ) * sin ( ssu - asu )     ) );

  *a = 2.0 * atan ( tan_a2 );

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / 
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) );

  *b = 2.0 * atan ( tan_b2 );

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / 
                  ( sin ( ssu ) * sin ( ssu - csu )     ) );

  *c = 2.0 * atan ( tan_c2 );

  return;
}
//****************************************************************************80

void sphere01_triangle_vertices_to_angles ( double v1[3], double v2[3], 
  double v3[3], double *a, double *b, double *c )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_VERTICES_TO_ANGLES: angles of spherical triangle on unit sphere.
//
//  Discussion:
//
//    A unit sphere centered at 0 in 3D satisfies the equation:
//
//      X*X + Y*Y + Z*Z = 1
//
//    A spherical triangle is specified by three points on the surface
//    of the sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
//
//    Output, double *A, *B, *C, the angles of the spherical triangle.
{
  double as;
  double bs;
  double cs;
//
//  Compute the lengths of the sides of the spherical triangle.
//
  sphere01_triangle_vertices_to_sides ( v1, v2, v3, &as, &bs, &cs );
//
//  Get the spherical angles.
//
  sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c );

  return;
}
//****************************************************************************80

double sphere01_triangle_vertices_to_area ( double v1[3], double v2[3], double v3[3] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_VERTICES_TO_AREA: area of a spherical triangle on unit sphere.
//
//  Discussion:
//
//    A unit sphere centered at 0 in 3D satisfies the equation:
//
//      X*X + Y*Y + Z*Z = 1
//
//    A spherical triangle is specified by three points on the surface
//    of the sphere.
//
//    The area formula is known as Girard's formula.
//
//    The area of a spherical triangle on the unit sphere is:
//
//      AREA = A + B + C - PI
//
//    where A, B and C are the (surface) angles of the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
//
//    Output, double STRI_VERTICES_TO_AREA, the area of the 
//    spherical triangle.
{
  double area;
  double a;
  double as;
  double b;
  double bs;
  double c;
  double cs;
//
//  Compute the lengths of the sides of the spherical triangle.
//
  sphere01_triangle_vertices_to_sides ( v1, v2, v3, &as, &bs, &cs );
//
//  Get the spherical angles.
//
  sphere01_triangle_sides_to_angles ( as, bs, cs, &a, &b, &c );
//
//  Get the area
//
  area = sphere01_triangle_angles_to_area ( a, b, c );

  return area;
}
//****************************************************************************80

double *sphere01_triangle_vertices_to_centroid ( double v1[3], double v2[3], double v3[3] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_VERTICES_TO_CENTROID: centroid of spherical triangle on unit sphere.
//
//  Discussion:
//
//    A sphere centered at 0 in 3D satisfies the equation:
//
//      X*X + Y*Y + Z*Z = 1
//
//    A spherical triangle is specified by three points on the sphere.
//
//    The (true) centroid of a spherical triangle is the point
//
//      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
//
//    Note that the true centroid does NOT, in general, lie on the sphere.  
//
//    The "flat" centroid VF is the centroid of the planar triangle defined by
//    the vertices of the spherical triangle.
//
//    The "spherical" centroid VS of a spherical triangle is computed by
//    the intersection of the geodesic bisectors of the triangle angles.
//    The spherical centroid lies on the sphere.
//
//    VF, VT and VS lie on a line through the center of the sphere.  We can
//    easily calculate VF by averaging the vertices, and from this determine
//    VS by normalizing.
//
//    (Of course, we still will not have actually computed VT, which lies
//    somewhere between VF and VS!)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
//
//    Output, double SPHERE01_TRIANGLE_VERTICES_TO_CENTROID[3], the coordinates of the 
//    "spherical centroid" of the spherical triangle.
//
{
# define DIM_NUM 3

  int i;
  double norm;
  double *vs;

  vs = new double[3];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    vs[i] = ( v1[i] + v2[i] + v3[i] ) / 3.0;
  }

  norm = r8vec_norm ( DIM_NUM, vs );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    vs[i] = vs[i] / norm;
  }

  return vs;
# undef DIM_NUM
}
//****************************************************************************80

void sphere01_triangle_vertices_to_midpoints ( double v1[3], double v2[3], double v3[3], 
  double m1[3], double m2[3], double m3[3] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS gets the midsides of a spherical triangle.
//
//  Discussion:
//
//    The points are assumed to lie on the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
//
//    Output, double M1[3], M2[3], M3[3], the coordinates of 
//    the midpoints of the sides of the spherical triangle.
//
{
  int i;
  double norm;

  for ( i = 0; i < 3; i++ )
  {
    m1[i] = ( v1[i] + v2[i] ) / 2.0;
  }
  norm = r8vec_norm ( 3, m1 );
  for ( i = 0; i < 3; i++ )
  {
    m1[i] = m1[i] / norm;
  }

  for ( i = 0; i < 3; i++ )
  {
    m2[i] = ( v2[i] + v3[i] ) / 2.0;
  }
  norm = r8vec_norm ( 3, m2 );
  for ( i = 0; i < 3; i++ )
  {
    m2[i] = m2[i] / norm;
  }

  for ( i = 0; i < 3; i++ )
  {
    m3[i] = ( v3[i] + v1[i] ) / 2.0;
  }
  norm = r8vec_norm ( 3, m3 );
  for ( i = 0; i < 3; i++ )
  {
    m3[i] = m3[i] / norm;
  }

  return;
}
//****************************************************************************80

void sphere01_triangle_vertices_to_sides ( double v1[3], double v2[3], 
  double v3[3], double *as, double *bs, double *cs )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_VERTICES_TO_SIDES_3D: sides of spherical triangle on unit sphere.
//
//  Discussion:
//
//    We can use the ACOS system call here, but the ARC_COSINE routine
//    will automatically take care of cases where the input argument is
//    (usually slightly) out of bounds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], V3[3], the vertices of the spherical
//    triangle.
//
//    Output, double *AS, *BS, *CS, the (geodesic) length of the sides of the
//    triangle.
//
{
  *as = arc_cosine ( r8vec_dot_product ( 3, v2, v3 ) );
  *bs = arc_cosine ( r8vec_dot_product ( 3, v3, v1 ) );
  *cs = arc_cosine ( r8vec_dot_product ( 3, v1, v2 ) );

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

double *tp_to_xyz ( double theta, double phi )

//****************************************************************************80
//
//  Purpose:
//
//    TP_TO_XYZ converts unit spherical TP coordinates to XYZ coordinates.
//
//  Discussion:
//
//    The point is assume to lie on the unit sphere centered at the origin.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double THETA, PHI, the angular coordinates of a point
//    on the unit sphere.
//
//    Output, double TP_TO_XYZ[3], the XYZ coordinates.
//
{
  double *v;

  v = new double[3];

  v[0] = cos ( theta ) * sin ( phi );
  v[1] = sin ( theta ) * sin ( phi );
  v[2] =                 cos ( phi );

  return v;
}
