# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sphere_triangle_quad.hpp"

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
    value = + x;
  } 
  else
  {
    value = - x;
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
//    SPHERE01_TRIANGLE_ANGLES_TO_AREA: area of spherical triangle on unit sphere.
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

double *sphere01_triangle_project ( double a_xyz[3], double b_xyz[3], 
  double c_xyz[3], int f1, int f2, int f3 )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_PROJECT projects from plane to spherical triangle.
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

double *sphere01_triangle_project2 ( double a_xyz[3], double b_xyz[3], 
  double c_xyz[3], int f1, int f2, int f3 )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_PROJECT2 projects from plane to spherical triangle.
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

double sphere01_triangle_quad_00 ( int n, double v1[3], double v2[3], 
  double v3[3], double f ( double x[] ), int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_00: quadrature over a triangle on the unit sphere.
//
//  Discussion:
//
//    This is a Monte Carlo approach.
//
//    The integral is approximated by averaging the values at N random points,
//    multiplied by the area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of sample points.
//
//    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
//    the vertices of the triangle.
//
//    Input, double F ( double x[] ), evaluates the integrand.
//
//    Input/output, integer *SEED, a seed for the random
//    number generator.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_00, the approximate integral.
//
{
  double area;
  int j;
  double quad;
  double result;
  double *vc;

  area = sphere01_triangle_vertices_to_area ( v1, v2, v3 );

  vc = sphere01_triangle_sample ( n, v1, v2, v3, seed );

  quad = 0.0;
  for ( j = 0; j < n; j++ )
  {
    quad = quad + f ( vc+3*j );
  }

  result = quad * area / ( double ) ( n );

  delete [] vc;

  return result;
}
//****************************************************************************80

double sphere01_triangle_quad_01 ( double v1[3], double v2[3], double v3[3], 
  double f ( double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_01: quadrature over a triangle on the unit sphere.
//
//  Discussion:
//
//    The integral is approximated by the value at the centroid,
//    multiplied by the area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
//    the vertices of the triangle.
//
//    Input, double F ( double x[] ), evaluates the integrand.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_01, the approximate integral.
//
{
  double area;
  double quad;
  double result;
  double *vc;

  area = sphere01_triangle_vertices_to_area ( v1, v2, v3 );

  vc = sphere01_triangle_vertices_to_centroid ( v1, v2, v3 );

  quad = f ( vc );
  result = quad * area;

  delete [] vc;

  return result;
}
//****************************************************************************80

double sphere01_triangle_quad_02 ( double v1[3], double v2[3], double v3[3], 
  double f ( double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_02: quadrature over a triangle on the unit sphere.
//
//  Discussion:
//
//    The integral is approximated by the average of the vertex values,
//    multiplied by the area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
//    the vertices of the triangle.
//
//    Input, double F ( double x[] ), evaluates the integrand.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_02, the approximate integral.
//
{
  double area;
  double quad;
  double result;

  area = sphere01_triangle_vertices_to_area ( v1, v2, v3 );

  quad = ( f ( v1 ) + f ( v2 ) + f ( v3 ) ) / 3.0;

  result = quad * area;

  return result;
}
//****************************************************************************80

double sphere01_triangle_quad_03 ( double v1[3], double v2[3], double v3[3], 
  double f ( double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_03: quadrature over a triangle on the unit sphere.
//
//  Discussion:
//
//    The integral is approximated by the average of the midside values,
//    multiplied by the area.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real V1[3], V2[3], V3[3], the XYZ coordinates of
//    the vertices of the triangle.
//
//    Input, double F ( double x[] ), evaluates the integrand.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_03, the approximate integral.
//
{
  double area;
  double quad;
  double result;
  double v4[3];
  double v5[3];
  double v6[3];

  area = sphere01_triangle_vertices_to_area ( v1, v2, v3 );

  sphere01_triangle_vertices_to_midpoints ( v1, v2, v3, v4, v5, v6 );

  quad = ( f ( v4 ) + f ( v5 ) + f ( v6 ) ) / 3.0;

  result = quad * area;

  return result;
}
//****************************************************************************80

double sphere01_triangle_quad_icos1c ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_ICOS1C: centroid rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over a spherical triangle on the
//    unit sphere.
//
//    This function subdivides each edge of the triangle into FACTOR subedges.  
//    These edges define a grid within the triangle.  The centroids of these
//    triangles can be determined.  All of these calculations are done,
//    essentially, on the FLAT faces of the planar triangle.  Only then are
//    the triangle vertices and centroids projected to the sphere.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
//    of the spherical triangle.
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, double FUN ( double x[] ), evaluates the integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_ICOS1C, the estimated integral.
//
{
  int a;
  double *a2_xyz;
  double area;
  double area_total;
  int b;
  double *b2_xyz;
  int c;
  double *c2_xyz;
  int f1;
  int f2;
  int f3;
  int i;
  double *node_xyz;
  double pi = 3.141592653589793;
  double result;
  double v;
//
//  Initialize the integral data.
//
  result = 0.0;
  area_total = 0.0;
  *node_num = 0;
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

      a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 2, 
        f2 - 1, f3 - 1 );
      b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, 
        f2 + 2, f3 - 1 );
      c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, 
        f2 - 1, f3 + 2 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      v = fun ( node_xyz );

      *node_num = *node_num + 1;
      result = result + area * v;
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

      a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 2, 
        f2 + 1, f3 + 1 );
      b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, 
        f2 - 2, f3 + 1 );
      c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, 
        f2 + 1, f3 - 2 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      v = fun ( node_xyz );

      *node_num = *node_num + 1;
      result = result + area * v;
      area_total = area_total + area;

      delete [] a2_xyz;
      delete [] b2_xyz;
      delete [] c2_xyz;
      delete [] node_xyz;
    }
  }

  return result;
}
//****************************************************************************80

double sphere01_triangle_quad_icos1m ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_ICOS1M: midside rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over a spherical triangle on the
//    unit sphere.
//
//    This function subdivides each edge of the triangle into FACTOR subedges.  
//    These edges define a grid within the triangle.  The midsides of the
//    edges of these triangles can be determined.  All of these calculations 
//    are done, essentially, on the FLAT faces of the planar triangle.  Only
//    then are the triangle vertices and midsides projected to the sphere.  
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
//  Parameters:
//
//    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
//    of the spherical triangle.
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, double FUN ( double x[] ), evaluates the integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_ICOS1M, the estimated integral.
//
{
  int a;
  double *a2_xyz;
  double *a3_xyz;
  double area;
  double area_total;
  int b;
  double *b2_xyz;
  double *b3_xyz;
  int c;
  double *c2_xyz;
  double *c3_xyz;
  int f1;
  int f2;
  int f3;
  int i;
  double pi = 3.141592653589793;
  double result;
  double va;
  double vb;
  double vc;
//
//  Initialize the integral data.
//
  result = 0.0;
  area_total = 0.0;
  *node_num = 0;
//
//  Some subtriangles will have the same direction as the face.
//
  for ( f1 = 0; f1 <= factor - 1; f1++ )
  {
    for ( f2 = 0; f2 <= factor - f1 - 1; f2++ )
    {
      f3 = factor - f1 - f2;

      a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, 
        f2,     f3 - 1 );
      b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,
        f2 + 1, f3 - 1 );
      c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, 
        f2,     f3 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      a3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 
        2 * f2 + 1, 2 * f3 - 2 );
      b3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1,
        2 * f2 + 1, 2 * f3 - 1 );
      c3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 
        2 * f2,     2 * f3 - 1 );

      *node_num = *node_num + 3;
      va = fun ( a3_xyz );
      vb = fun ( b3_xyz );
      vc = fun ( c3_xyz );
      result = result + area * ( va + vb + vc ) / 3.0;
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
//  The other subtriangles have the opposite direction from the face.
//
  for ( f3 = 0; f3 <= factor - 2; f3++ )
  {
    for ( f2 = 1; f2 <= factor - f3 - 1; f2++ )
    {
      f1 = factor - f2 - f3;

      a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, 
        f2,     f3 + 1 );
      b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,
        f2 - 1, f3 + 1 );
      c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, 
        f2,     f3 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      a3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 
        2 * f2 - 1, 2 * f3 + 2 );
      b3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1,
        2 * f2 - 1, 2 * f3 + 1 );
      c3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 
        2 * f2,     2 * f3 + 1 );

      *node_num = *node_num + 3;
      va = fun ( a3_xyz );
      vb = fun ( b3_xyz );
      vc = fun ( c3_xyz );
      result = result + area * ( va + vb + vc ) / 3.0;
      area_total = area_total + area;

      delete [] a2_xyz;
      delete [] a3_xyz;
      delete [] b2_xyz;
      delete [] b3_xyz;
      delete [] c2_xyz;
      delete [] c3_xyz;
    }
  }

  return result;
}
//****************************************************************************80

double sphere01_triangle_quad_icos1v ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_ICOS1V: vertex rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over a spherical triangle on the
//    unit sphere.
//
//    This function subdivides each edge of the triangle into FACTOR subedges.  
//    These edges define a grid within the triangle.  The vertices of these
//    triangles can be determined.  All of these calculations 
//    are done, essentially, on the FLAT faces of the planar triangle.  Only
//    then are the triangle vertices projected to the sphere.  
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
//  Parameters:
//
//    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
//    of the spherical triangle.
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, double FUN ( double x[] ), evaluates the integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_ICOS1V, the estimated integral.
//
{
  int a;
  double *a2_xyz;
  double area;
  double area_total;
  int b;
  double *b2_xyz;
  int c;
  double *c2_xyz;
  int f1;
  int f2;
  int f3;
  int i;
  double pi = 3.141592653589793;
  double result;
  double va;
  double vb;
  double vc;
//
//  Initialize the integral data.
//
  result = 0.0;
  area_total = 0.0;
  *node_num = 0;
//
//  Some subtriangles will have the same direction as the face.
//
  for ( f1 = 0; f1 <= factor - 1; f1++ )
  {
    for ( f2 = 0; f2 <= factor - f1 - 1; f2++ )
    {
      f3 = factor - f1 - f2;

      a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, 
        f2,     f3 - 1 );
      b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,
        f2 + 1, f3 - 1 );
      c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, 
        f2,     f3 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      *node_num = *node_num + 1;
      va = fun ( a2_xyz );
      vb = fun ( b2_xyz );
      vc = fun ( c2_xyz );
      result = result + area * ( va + vb + vc ) / 3.0;
      area_total = area_total + area;

      delete [] a2_xyz;
      delete [] b2_xyz;
      delete [] c2_xyz;
    }
  }
//
//  The other subtriangles have the opposite direction from the face.
//
  for ( f3 = 0; f3 <= factor - 2; f3++ )
  {
    for ( f2 = 1; f2 <= factor - f3 - 1; f2++ )
    {
      f1 = factor - f2 - f3;

      a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, 
        f2,     f3 + 1 );
      b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,
        f2 - 1, f3 + 1 );
      c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, 
        f2,     f3 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      *node_num = *node_num + 1;
      va = fun ( a2_xyz );
      vb = fun ( b2_xyz );
      vc = fun ( c2_xyz );
      result = result + area * ( va + vb + vc ) / 3.0;
      area_total = area_total + area;

      delete [] a2_xyz;
      delete [] b2_xyz;
      delete [] c2_xyz;
    }
  }

  return result;
}
//****************************************************************************80

double sphere01_triangle_quad_icos2v ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_QUAD_ICOS2V: vertex rule, subdivide then project.
//
//  Discussion:
//
//    This function estimates an integral over a spherical triangle on the
//    unit sphere.
//
//    This function subdivides each edge of the triangle into FACTOR subedges.  
//    These edges define a grid within the triangle.  The vertices of these
//    triangles can be determined.  All of these calculations 
//    are done, essentially, on the FLAT faces of the planar triangle.  Only
//    then are the triangle vertices projected to the sphere.  
//
//    This function uses a more sophisticated projection scheme than that
//    used by SPHERE01_TRIANGLE_QUAD_ICOS1V.
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
//  Parameters:
//
//    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the vertices
//    of the spherical triangle.
//
//    Input, int FACTOR, the subdivision factor, which must
//    be at least 1.
//
//    Input, double FUN ( double x[] ), evaluates the integrand.
//
//    Output, int *NODE_NUM, the number of evaluation points.
//
//    Output, double SPHERE01_TRIANGLE_QUAD_ICOS2V, the estimated integral.
//
{
  int a;
  double *a2_xyz;
  double area;
  double area_total;
  int b;
  double *b2_xyz;
  int c;
  double *c2_xyz;
  int f1;
  int f2;
  int f3;
  int i;
  double pi = 3.141592653589793;
  double result;
  double va;
  double vb;
  double vc;
//
//  Initialize the integral data.
//
  result = 0.0;
  area_total = 0.0;
  *node_num = 0;
//
//  Some subtriangles will have the same direction as the face.
//
  for ( f1 = 0; f1 <= factor - 1; f1++ )
  {
    for ( f2 = 0; f2 <= factor - f1 - 1; f2++ )
    {
      f3 = factor - f1 - f2;

      a2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1 + 1, 
        f2,     f3 - 1 );
      b2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,
        f2 + 1, f3 - 1 );
      c2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1, 
        f2,     f3 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      *node_num = *node_num + 1;
      va = fun ( a2_xyz );
      vb = fun ( b2_xyz );
      vc = fun ( c2_xyz );
      result = result + area * ( va + vb + vc ) / 3.0;
      area_total = area_total + area;

      delete [] a2_xyz;
      delete [] b2_xyz;
      delete [] c2_xyz;
    }
  }
//
//  The other subtriangles have the opposite direction from the face.
//
  for ( f3 = 0; f3 <= factor - 2; f3++ )
  {
    for ( f2 = 1; f2 <= factor - f3 - 1; f2++ )
    {
      f1 = factor - f3 - f2;

      a2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1 - 1, 
        f2,     f3 + 1 );
      b2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,
        f2 - 1, f3 + 1 );
      c2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1, 
        f2,     f3 );

      area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

      *node_num = *node_num + 1;
      va = fun ( a2_xyz );
      vb = fun ( b2_xyz );
      vc = fun ( c2_xyz );
      result = result + area * ( va + vb + vc ) / 3.0;
      area_total = area_total + area;

      delete [] a2_xyz;
      delete [] b2_xyz;
      delete [] c2_xyz;
    }
  }

  return result;
}
//****************************************************************************80

double *sphere01_triangle_sample ( int n, double v1[3], double v2[3], 
  double v3[3], int *seed )

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

double sphere01_triangle_vertices_to_area ( double v1[3], double v2[3], 
  double v3[3] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_VERTICES_TO_AREA: area of spherical triangle on unit sphere.
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

double *sphere01_triangle_vertices_to_centroid ( double v1[3], double v2[3], 
  double v3[3] )

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
//    Output, double SPHERE01_TRIANGLE_VERTICES_TO_CENTROID[3], the coordinates 
//    of the "spherical centroid" of the spherical triangle.
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

void sphere01_triangle_vertices_to_midpoints ( double v1[3], double v2[3], 
  double v3[3], double m1[3], double m2[3], double m3[3] )

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
