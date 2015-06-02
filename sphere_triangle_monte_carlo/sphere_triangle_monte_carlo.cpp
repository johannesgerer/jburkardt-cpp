# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sphere_triangle_monte_carlo.hpp"

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

double *monomial_value ( int m, int n, int e[], double x[] )

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
//      product ( 1 <= i <= m ) x(i)^e(i)
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
//    Input, int M, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, int E[M], the exponents.
//
//    Input, double X[M*N], the point coordinates.
//
//    Output, double MONOMIAL_VALUE[N], the value of the monomial.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    if ( 0 != e[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], e[i] );
      }
    }
  }

  return v;
}
//****************************************************************************80

double r8_acos ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOS computes the arc cosine function, with argument truncation.
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
//    Output, double R8_ACOS, an angle whose cosine is C.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  if ( c <= -1.0 )
  {
    value = r8_pi;
  }
  else if ( 1.0 <= c )
  {
    value = 0.0;
  }
  else
  {
    value = acos ( c );
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
  const int i4_huge = 2147483647;
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
//      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
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

void r8vec_normalize ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMALIZE normalizes an R8VEC.
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
//    11 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, double A[N], the vector to be normalized.
//    On output, A should have unit Euclidean norm.
//
{
  int i;
  double norm;

  norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    norm = norm + a[i] * a[i];
  }
  norm = sqrt ( norm );

  if ( norm == 0.0 )
  {
    cerr << "\n";
    cerr << "R8VEC_NORMALIZE - Fatal error!\n";
    cerr << "  The vector norm is 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / norm;
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

void r8vec_transpose_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Example:
//
//    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
//    TITLE = 'My vector:  '
//
//    My vector:
//        1.0    2.1    3.2    4.3    5.4
//        6.5    7.6    8.7    9.8   10.9
//       11.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
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
  int ihi;
  int ilo;

  cout << "\n";
  cout << title << "\n";

  if ( n <= 0 )
  {
    cout << "  (Empty)\n";
    return;
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      cout << "  " << setw(12) << a[i];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

double *sphere01_sample ( int n, int &seed )

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
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double SPHERE01_SAMPLE[3*N], the sample points.
//
{
  int j;
  double phi;
  const double r8_pi = 3.141592653589793;
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

    phi = r8_acos ( vdot );
//
//  Pick a uniformly random rotation between 0 and 2 Pi around the
//  axis of the Z vector.
//
    theta = r8_uniform_01 ( seed );
    theta = 2.0 * r8_pi * theta;

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
//    Output, double SPHERE_TRIANGLE_ANGLES_TO_AREA, the area of the
//    spherical triangle.
//
{
  double area;
  const double r8_pi = 3.141592653589793;

  area = a + b + c - r8_pi;

  return area;
}
//****************************************************************************80

void sphere01_triangle_sides_to_angles ( double as, double bs, double cs,
  double &a, double &b, double &c )

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
//    Output, double &A, &B, &C, the spherical angles of the triangle.
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

  a = 2.0 * atan ( tan_a2 );

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) /
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) );

  b = 2.0 * atan ( tan_b2 );

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) /
                  ( sin ( ssu ) * sin ( ssu - csu )     ) );

  c = 2.0 * atan ( tan_c2 );

  return;
}
//****************************************************************************80

double *sphere01_triangle_sample ( int n, double v1[3], double v2[3], 
  double v3[3], int &seed )

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
//    Input/output, int &SEED, a seed for the random number generator.
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

  sphere01_triangle_vertices_to_sides ( v1, v2, v3, a, b, c );

  sphere01_triangle_sides_to_angles ( a, b, c, alpha, beta, gamma );

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
  sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs );
//
//  Get the spherical angles.
//
  sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c );
//
//  Get the area
//
  area = sphere01_triangle_angles_to_area ( a, b, c );

  return area;
}
//****************************************************************************80

void sphere01_triangle_vertices_to_sides ( double v1[3], double v2[3], 
  double v3[3], double &as, double &bs, double &cs )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE01_TRIANGLE_VERTICES_TO_SIDES_3D: sides of spherical triangle on unit sphere.
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
  as = r8_acos ( r8vec_dot_product ( 3, v2, v3 ) );
  bs = r8_acos ( r8vec_dot_product ( 3, v3, v1 ) );
  cs = r8_acos ( r8vec_dot_product ( 3, v1, v2 ) );

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
