# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sphere_grid.hpp"

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

string i4_to_string ( int i4, string format )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  char i4_char[80];
  string i4_string;

  sprintf ( i4_char, format.c_str ( ), i4 );

  i4_string = string ( i4_char );

  return i4_string;
}
//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 January 2005
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
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i - 1 << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j - 1 << ":";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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
//    The input data required for this routine can be retrieved from ICOS_NUM.
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
//    22 July 2007
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
     0,  1,
     0,  2,
     0,  3,
     0,  4,
     0,  5,
     1,  2,
     1,  3,
     1,  6,
     1,  7,
     2,  4,
     2,  6,
     2,  8,
     3,  5,
     3,  7,
     3,  9,
     4,  5,
     4,  8,
     4, 10,
     5,  9,
     5, 10,
     6,  7,
     6,  8,
     6, 11,
     7,  9,
     7, 11,
     8, 10,
     8, 11,
     9, 10,
     9, 11,
    10, 11 };
  static int face_order_save[FACE_NUM] = {
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };
  static int face_point_save[3*FACE_NUM] = {
     0,  1,  3,
     0,  2,  1,
     0,  3,  5,
     0,  4,  2,
     0,  5,  4,
     1,  2,  6,
     1,  6,  7,
     1,  7,  3,
     2,  4,  8,
     2,  8,  6,
     3,  7,  9,
     3,  9,  5,
     4,  5, 10,
     4, 10,  8,
     5,  9, 10,
     6,  8, 11,
     6, 11,  7,
     7, 11,  9,
     8, 10, 11,
     9, 11, 10 };
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

  return;
# undef DIM_NUM
# undef EDGE_NUM
# undef EDGE_ORDER
# undef FACE_NUM
# undef POINT_NUM
}
//****************************************************************************80

void icos_num ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max )

//****************************************************************************80
//
//  Purpose:
//
//    ICOS_NUM gives "sizes" for an icosahedron.
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

double r8_modp ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MODP returns the nonnegative remainder of R8 division.
//
//  Discussion:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360.0) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
//
//    If
//      REM = R8_MODP ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM is always nonnegative.
//
//  Example:
//
//        I         J     MOD  R8_MODP   R8_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number to be divided.
//
//    Input, double Y, the number that divides X.
//
//    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
//
{
  double value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_MODP - Fatal error!\n";
    cerr << "  R8_MODP ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( value < 0.0 )
  {
    value = value + r8_abs ( y );
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

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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

double r8vec_diff_norm ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
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
//    Input, double A[N], B[N], the vectors.
//
//    Output, double R8VEC_DIFF_NORM, the L2 norm of A - B.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] + b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
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

bool r8vec_eq ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EQ is true if every pair of entries in two R8VEC's is equal.
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
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], two vectors to compare.
//
//    Output, bool R8VEC_EQ, is TRUE if every pair of elements A1(I) 
//    and A2(I) are equal, and FALSE otherwise.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      return false;
    }
  }
  return true;
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

double *sphere_cubed_ijk_to_xyz ( int n, int i, int j, int k )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CUBED_IJK_TO_XYZ: cubed sphere IJK to XYZ coordinates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of sections into which each 
//    face of the cube is to be divided.
//
//    Input, int I, J, K, indices between 0 and N.  Normally,
//    at least one of the indices should have the value 0 or N.
//
//    Output, double SPHERE_CUBED_IJK_TO_XYZ[3], coordinates of the point.
//
{
  double pi = 3.141592653589793;
  double xc;
  double *xyz;
  double xyzn;
  double yc;
  double zc;

  xyz = new double[3];

  if ( i == 0 )
  {
    xc = -1.0;
  }
  else if ( i == n )
  {
    xc = +1.0;
  }
  else
  {
    xc = tan ( ( double ) ( 2 * i - n ) * 0.25 * pi / ( double ) ( n ) );
  }

  if ( j == 0 )
  {
    yc = -1.0;
  }
  else if ( j == n )
  {
    yc = +1.0;
  }
  else
  {
    yc = tan ( ( double ) ( 2 * j - n ) * 0.25 * pi / ( double ) ( n ) );
  }

  if ( k == 0 )
  {
    zc = -1.0;
  }
  else if ( k == n )
  {
    zc = +1.0;
  }
  else
  {
    zc = tan ( ( double ) ( 2 * k - n ) * 0.25 * pi / ( double ) ( n ) );
  }

  xyzn = sqrt ( xc * xc + yc * yc + zc * zc );

  xyz[0] = xc / xyzn;
  xyz[1] = yc / xyzn;
  xyz[2] = zc / xyzn;

  return xyz;
}
//****************************************************************************80

int sphere_cubed_line_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CUBED_LINE_NUM counts lines on a cubed sphere grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of sections into which each 
//    face of the cube is to be divided.
//
//    Output, int LINE_NUM, the number of lines.
//
{
  int line_num;

  line_num = 0;
//
//  If N = 1, the corners form 12 lines.
//
  if ( n == 1 )
  {
    line_num = 12;
    return line_num;
  }
//
//  If 1 < N, each of 8 corners connects to three neighboring edges.
//
  else
  {
    line_num = line_num + 8 * 3;
  }
//
//  If 2 < N, then each of the 12 edges includes lines.
//
  if ( 2 < n )
  {
    line_num = line_num + 12 * ( n - 2 );
  }
//
//  Lines that belong to one of the six faces.
//
  if ( 1 < n )
  {
    line_num = line_num + 6 * 2 * n * ( n - 1 );
  }

  return line_num;
}
//****************************************************************************80

double *sphere_cubed_lines ( int n, int line_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CUBED_LINES computes the lines on a cubed sphere grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of sections into which 
//    each face of the cube is to be divided.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Output, double SPHERE_CUBED_LINES[3*2*LINE_NUM], distinct points 
//    on the unit sphere generated by a cubed sphere grid.
//
{
  int i;
  int j;
  int l;
  double *line_data;

  line_data = new double[3*2*line_num];
  l = 0;
//
//  If N = 1, the corners form 12 lines.
//
  if ( n == 1 )
  {
    sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, n, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, n, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, n, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+1*3+l*6 );
    l = l + 1;
    return line_data;
  }
//
//  If 1 < N, each of 8 corners connects to three neighboring edges.
//
  else
  {
    sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 1, 0, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 1, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 0, 1, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n,   0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n-1, 0, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 1, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, 0, 1, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n,   n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n-1, n, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, n,   0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, n-1, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, n, 1, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 1, n, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n,   0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, n-1, 0, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, n, 1, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 1, 0, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 1, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, 0, n,   line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, 0, n-1, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n,   0, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n-1, 0, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, 1, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, 0, n,   line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, 0, n-1, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n,   n, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n-1, n, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, n,   n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, n-1, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, n, n, n,   line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, n, n, n-1, line_data+0+1*3+l*6 );

    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 1, n, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n,   n, line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, n-1, n, line_data+0+1*3+l*6 );
    l = l + 1;
    sphere_cubed_ijk_to_xyz ( n, 0, n, n,   line_data+0+0*3+l*6 );
    sphere_cubed_ijk_to_xyz ( n, 0, n, n-1, line_data+0+1*3+l*6 );
    l = l + 1;
  }
//
//  If 2 < N, then each of the 12 edges includes lines.
//
  if ( 2 < n )
  {
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, i,   0, 0, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, i+1, 0, 0, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, n,   i, 0, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, n, i+1, 0, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, n-i,   n, 0, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, n-i-1, n, 0, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, 0, n-i,   0, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, 0, n-i-1, 0, line_data+0+1*3+l*6 );
      l = l + 1;
    }

    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, i,   0, n, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, i+1, 0, n, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, n,   i, n, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, n, i+1, n, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, n-i,   n, n, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, n-i-1, n, n, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, 0, n-i,   n, line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, 0, n-i-1, n, line_data+0+1*3+l*6 );
      l = l + 1;
    }

    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, 0, 0, i,   line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, 0, 0, i+1, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, n, 0, i,   line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, n, 0, i+1, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, n, n, i,   line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, n, n, i+1, line_data+0+1*3+l*6 );
      l = l + 1;
    }
    for ( i = 1; i <= n - 2; i++ )
    {
      sphere_cubed_ijk_to_xyz ( n, 0, n, i,   line_data+0+0*3+l*6 );
      sphere_cubed_ijk_to_xyz ( n, 0, n, i+1, line_data+0+1*3+l*6 );
      l = l + 1;
    }
  }
//
//  Lines that belong to one of the six faces.
//
  if ( 1 < n )
  {
//
//  000, nn0
//
    for ( i = 1; i <= n - 1; i++ )
    {
      for ( j = 0; j <= n - 1; j++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i, j,   0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i, j+1, 0, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( i = 0; i <= n - 1; i++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i,   j, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i+1, j, 0, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
//
//  00n, nnn
//
    for ( i = 1; i <= n - 1; i++ )
    {
      for ( j = 0; j <= n - 1; j++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i, j,   n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i, j+1, n, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( i = 0; i <= n - 1; i++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i,   j, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i+1, j, n, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
//
//  000:n0n
//
    for ( i = 1; i <= n - 1; i++ )
    {
      for ( j = 0; j <= n - 1; j++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i, 0, j,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i, 0, j+1, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( i = 0; i <= n - 1; i++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i,   0, j, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i+1, 0, j, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
//
//  0n0:nnn
//
    for ( i = 1; i <= n - 1; i++ )
    {
      for ( j = 0; j <= n - 1; j++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i, n, j,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i, n, j+1, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( i = 0; i <= n - 1; i++ )
      {
        sphere_cubed_ijk_to_xyz ( n, i,   n, j, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, i+1, n, j, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
//
//  000:0nn
//
    for ( i = 1; i <= n - 1; i++ )
    {
      for ( j = 0; j <= n - 1; j++ )
      {
        sphere_cubed_ijk_to_xyz ( n, 0, i, j,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, i, j+1, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( i = 0; i <= n - 1; i++ )
      {
        sphere_cubed_ijk_to_xyz ( n, 0, i,   j, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, i+1, j, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
//
//  n00:nnn
//
    for ( i = 1; i <= n - 1; i++ )
    {
      for ( j = 0; j <= n - 1; j++ )
      {
        sphere_cubed_ijk_to_xyz ( n, n, i, j,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, i, j+1, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( i = 0; i <= n - 1; i++ )
      {
        sphere_cubed_ijk_to_xyz ( n, n, i,   j, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, i+1, j, line_data+0+1*3+l*6 );
        l = l + 1;
      }
    }

  }

  if ( l != line_num )
  {
    cerr << "\n";
    cerr << "SPHERE_CUBED_LINES - Fatal error!\n";
    cerr << "  LINE_NUM = " << line_num << "\n";
    cerr << "  L = " << l << "n";
    exit ( 1 );
  }

  return line_data;
}
//****************************************************************************80

double *sphere_cubed_points ( int n, int ns )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CUBED_POINTS computes the points on a cubed sphere grid.
//
//  Discussion:
//
//    For a value of N = 3, for instance, each of the 6 cube faces will
//    be divided into 3 sections, so that a single cube face will have
//    (3+1)x(3+1) points:
//
//      X---X---X---X
//      | 1 | 4 | 7 |
//      X---X---X---X
//      | 2 | 5 | 8 |
//      X---X---X---X
//      | 3 | 6 | 9 |
//      X---X---X---X
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of sections into which each 
//    face of the cube is to be divided.
//
//    Input, int NS, the number of points.
//
//    Output, double SPHERE_CUBED_POINTS[3*NS], distinct points on the 
//    unit sphere generated by a cubed sphere grid.
//
{
  int ns2;
  double *xyz;

  xyz = new double[3*ns];

  ns2 = 0;
//
//  Bottom full.
//
  sphere_cubed_points_face ( n, 0, 0, 0, n, n, 0, ns2, xyz );
//
//  To avoid repetition, draw the middles as grids of n-2 x n-1 points.
//
  sphere_cubed_points_face ( n, 0, 0, 1, 0,   n-1, n-1, ns2, xyz );
  sphere_cubed_points_face ( n, 0, n, 1, n-1, n,   n-1, ns2, xyz );
  sphere_cubed_points_face ( n, n, 1, 1, n,   n,   n-1, ns2, xyz );
  sphere_cubed_points_face ( n, 1, 0, 1, n,   0,   n-1, ns2, xyz );
//
//  Top full.
//
  sphere_cubed_points_face ( n, 0, 0, n, n, n, n, ns2, xyz );
  
  if ( ns2 != ns )
  {
    cerr << "\n";
    cerr << "SPHERE_CUBED_POINTS - Fatal error\n";
    cerr << "  Expected to generated NS = " << ns << " points.\n";
    cerr << "  Generated " << ns2 << " points.\n";
    exit ( 1 );
  }
  return xyz;
}
//****************************************************************************80

void sphere_cubed_points_face ( int n, int i1, int j1, int k1, int i2, int j2, 
  int k2, int &ns, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CUBED_POINTS_FACE: points on one face of a cubed sphere grid.
//
//  Discussion:
//
//    This routine starts with NS = 0, and is called repeatedly to
//    add points for another face.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of sections into which each face 
//    of the cube is to be divided.
//
//    Input, int I1, J1, K1, I2, J2, K2, the logical indices, 
//    between 0 and N, of two corners of the face grid.  It is guaranteed that 
//    I1 <= I2, J1 <= J2, and K1 <= K2.  
//
//    Input/output, int &NS, the number of points.
//
//    Input/output, double XYZ[3*NS], distinct points on the unit sphere
//    generated by a cubed sphere grid.
//
{
  int i;
  int j;
  int k;
  double pi = 3.141592653589793;
  double xyzn;
  double xc;
  double yc;
  double zc;

  for ( i = i1; i <= i2; i++ )
  {
    if ( i1 < i2 )
    {
      xc = tan ( ( double ) ( 2 * i - n ) * 0.25 * pi / ( double ) ( n ) );
    }
    else if ( i1 == 0 )
    {
      xc = -1.0;
    }
    else if ( i1 == n )
    {
      xc = +1.0;
    }
    else
    {
      xc = 0.0;
    }

    for ( j = j1; j <= j2; j++ )
    {
      if ( j1 < j2 )
      {
        yc = tan ( ( double ) ( 2 * j - n ) * 0.25 * pi / ( double ) ( n ) );
      }
      else if ( j1 == 0 )
      {
        yc = -1.0;
      }
      else if ( j1 == n )
      {
        yc = +1.0;
      }
      else
      {
        yc = 0.0;
      }

      for ( k = k1; k <= k2; k++ )
      {
        if ( k1 < k2 )
        {
          zc = tan ( ( double ) ( 2 * k - n ) * 0.25 * pi / ( double ) ( n ) );
        }
        else if ( k1 == 0 )
        {
          zc = -1.0;
        }
        else if ( k1 == n )
        {
          zc = +1.0;
        }
        else
        {
          zc = 0.0;
        }

        xyzn = sqrt ( xc * xc + yc * yc + zc * zc );

        xyz[0+ns*3] = xc / xyzn;
        xyz[1+ns*3] = yc / xyzn;
        xyz[2+ns*3] = zc / xyzn;
        ns = ns + 1;
      }
    }
  }

  return;
}
//****************************************************************************80

int sphere_cubed_point_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_CUBED_POINT_NUM counts the points on a cubed sphere grid.
//
//  Discussion:
//
//    For a value of N = 3, for instance, each of the 6 cube faces will
//    be divided into 3 sections, so that a single cube face will have
//    (3+1)x(3+1) points:
//
//      X---X---X---X
//      | 1 | 4 | 7 |
//      X---X---X---X
//      | 2 | 5 | 8 |
//      X---X---X---X
//      | 3 | 6 | 9 |
//      X---X---X---X
//
//    The number of points is simply (N+1)^3 - (N-1)^3.
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
//  Parameters:
//
//    Input, int N, the number of sections into which 
//    each face of the cube is to be divided.
//
//    Output, int SPHERE_CUBED_POINT_NUM, the number of points.
//
{
  int ns;

  ns = i4_power ( n + 1, 3 ) - i4_power ( n - 1, 3 );

  return ns;
}
//****************************************************************************80

double sphere_distance_xyz ( double xyz1[3], double xyz2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_DISTANCE_XYZ computes great circle distances on a sphere.
//
//  Discussion:
//
//    XYZ coordinates are used.
//
//    We assume the points XYZ1 and XYZ2 lie on the same sphere.
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
//    29 August 2010
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

  r = r8vec_norm ( 3, xyz1 );

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

  dist = r * atan2 ( top, bot );

  return dist;
}
//****************************************************************************80

int *sphere_grid_q4 ( int lat_num, int long_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q4: rectangular grid on a sphere.
//
//  Discussion:
//
//    The point numbering system is the same used in SPHERE_GRIDPOINTS,
//    and that routine may be used to compute the coordinates of the points.
//
//    A sphere in 3D satisfies the equation:
//
//      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT_NUM, the number of "rows" of rectangles to
//    be created.  LAT_NUM must be at least 2. 
//
//    Input, int LONG_NUM, the number of "columns" of 
//    rectangles to be created.
//
//    Output, int SPHERE_GRID_Q4[4*LAT_NUM*LONG_NUM], 
//    the indices of the nodes that make up each rectangle.
//
{
  int i;
  int j;
  int n;
  int n_max;
  int n_min;
  int ne;
  int nw;
  int s;
  int s_max;
  int s_min;
  int se;
  int sw;
  int *rectangle_node;
  int rectangle_num;

  rectangle_node = new int[4*(lat_num*long_num)];
  rectangle_num = 0;
//
//  The first row.
//
  n = 0;

  sw = 1;
  se = sw + 1;

  s_min = 1;
  s_max = long_num;

  for ( j = 1; j <= long_num; j++ )
  {
    rectangle_node[0+rectangle_num*4] = sw;
    rectangle_node[1+rectangle_num*4] = se;
    rectangle_node[2+rectangle_num*4] = n;
    rectangle_node[3+rectangle_num*4] = n;
    rectangle_num = rectangle_num + 1;

    sw = se;

    if ( se == s_max )
    {
      se = s_min;
    }
    else
    {
      se = se + 1;
    }
  }
//
//  The intermediate rows.
//
  for ( i = 2; i < lat_num; i++ )
  {
    n_max = s_max;
    n_min = s_min;

    s_max = s_max + long_num;
    s_min = s_min + long_num;

    nw = n_min;
    ne = nw + 1;
    sw = s_min;
    se = sw + 1;

    for ( j = 1; j <= long_num; j++ )
    {

      rectangle_node[0+rectangle_num*4] = sw;
      rectangle_node[1+rectangle_num*4] = se;
      rectangle_node[2+rectangle_num*4] = ne;
      rectangle_node[3+rectangle_num*4] = nw;
      rectangle_num = rectangle_num + 1;

      sw = se;
      nw = ne;
      if ( se == s_max )
      {
        se = s_min;
      }
      else
      {
        se = se + 1;
      }
      if ( ne == n_max )
      {
        ne = n_min;
      }
      else
      {
        ne = ne + 1;
      }
    }
  }
//
//  The last row.
//
  n_max = s_max;
  n_min = s_min;

  s = n_max + 1;

  nw = n_min;
  ne = nw + 1;

  for ( j = 1; j <= long_num; j++ )
  {
    rectangle_node[0+rectangle_num*4] = ne;
    rectangle_node[1+rectangle_num*4] = nw;
    rectangle_node[2+rectangle_num*4] = s;
    rectangle_node[3+rectangle_num*4] = s;
    rectangle_num = rectangle_num + 1;

    nw = ne;
    if ( ne == n_max )
    {
      ne = n_min;
    }
    else
    {
      ne = ne + 1;
    }
  }
  return rectangle_node;
}
//****************************************************************************80

int *sphere_grid_t3 ( int lat_num, int long_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T3 produces a triangle grid on a sphere.
//
//  Discussion:
//
//    I want to convert this to 0-based indexing.  However, my first attempt
//    to modify the internal coding to compute the 0-based indexes directly
//    was a failure.  So I reverted to computing 1-based node indices,
//    and then subtracting 1.  You would think this would be easy to fix.
//
//    The point numbering system is the same used in SPHERE_GRIDPOINTS,
//    and that routine may be used to compute the coordinates of the points.
//
//    A sphere in 3D satisfies the equation:
//
//      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT_NUM, LONG_NUM, the number of latitude 
//    and longitude lines to draw.  The latitudes do not include the North 
//    and South poles, which will be included automatically, so LAT_NUM = 5, 
//    for instance, will result in points along 7 lines of latitude.
//
//    Output, int SPHERE_GRID_T3[3*2*(LAT_NUM+1)*LONG_NUM], the
//    triangle vertices.
//
{
  int i;
  int j;
  int n;
  int n_max;
  int n_min;
  int ne;
  int nw;
  int s;
  int s_max;
  int s_min;
  int se;
  int sw;
  int *triangle_node;
  int triangle_num;

  triangle_node = new int[3*2*(lat_num+1)*long_num];
  triangle_num = 0;
//
//  The first row.
//
//n = 0;
//sw = 1;
//se = sw + 1;
//s_min = 1;
//s_max = long_num;

  n = 1;
  sw = 2;
  se = sw + 1;
  s_min = 2;
  s_max = long_num + 1;

  for ( j = 0; j < long_num; j++ )
  {
    triangle_node[0+triangle_num*3] = sw - 1;
    triangle_node[1+triangle_num*3] = se - 1;
    triangle_node[2+triangle_num*3] = n - 1;
    triangle_num = triangle_num + 1;

    sw = se;
    if ( se == s_max )
    {
      se = s_min;
    }
    else
    {
      se = se + 1;
    }
  }
//
//  The intermediate rows.
//
  for ( i = 1; i <= lat_num; i++ )
  {
    n_max = s_max;
    n_min = s_min;

    s_max = s_max + long_num;
    s_min = s_min + long_num;

    nw = n_min;
    ne = nw + 1;
    sw = s_min;
    se = sw + 1;

    for ( j = 0; j < long_num; j++ )
    {
      triangle_node[0+triangle_num*3] = sw - 1;
      triangle_node[1+triangle_num*3] = se - 1;
      triangle_node[2+triangle_num*3] = nw - 1;
      triangle_num = triangle_num + 1;

      triangle_node[0+triangle_num*3] = ne - 1;
      triangle_node[1+triangle_num*3] = nw - 1;
      triangle_node[2+triangle_num*3] = se - 1;
      triangle_num = triangle_num + 1;

      sw = se;
      nw = ne;
      if ( se == s_max )
      {
        se = s_min;
      }
      else
      {
        se = se + 1;
      }

      if ( ne == n_max )
      {
        ne = n_min;
      }
      else
      {
        ne = ne + 1;
      }
    }
  }
//
//  The last row.
//
  n_max = s_max;
  n_min = s_min;

  s = n_max + 1;

  nw = n_min;
  ne = nw + 1;

  for ( j = 0; j < long_num; j++ )
  {
    triangle_node[0+triangle_num*3] = ne - 1;
    triangle_node[1+triangle_num*3] = nw - 1;
    triangle_node[2+triangle_num*3] = s - 1;
    triangle_num = triangle_num + 1;

    nw = ne;

    if ( ne == n_max )
    {
      ne = n_min;
    }
    else
    {
      ne = ne + 1;
    }
  }
  return triangle_node;
}
//****************************************************************************80

int sphere_icos_edge_num ( int factor )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_ICOS_EDGE_NUM sizes an icosahedral grid on a sphere.
//
//  Discussion:
//
//    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
//
//    With FACTOR = 2, each triangle of the icosahedron is subdivided into
//    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
//    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
//
//    With FACTOR = 3, each triangle of the icosahedron is subdivided into
//    3x3 subtriangles, resulting in 180 faces, 270 edges and 
//    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
//
//    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
//    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
//      12 
//    + 20 * 3          * (FACTOR-1) / 2 
//    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
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
//    Output, int SPHERE_ICOS_EDGE_NUM, the number of edges.
//
{
  int edge_num;

  edge_num = 30 * factor * factor;

  return edge_num;
}
//****************************************************************************80

int sphere_icos_face_num ( int factor )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_ICOS_FACE_NUM sizes an icosahedral grid on a sphere.
//
//  Discussion:
//
//    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
//
//    With FACTOR = 2, each triangle of the icosahedron is subdivided into
//    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
//    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
//
//    With FACTOR = 3, each triangle of the icosahedron is subdivided into
//    3x3 subtriangles, resulting in 180 faces, 270 edges and 
//    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
//
//    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
//    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
//      12 
//    + 20 * 3          * (FACTOR-1) / 2 
//    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
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
//    Output, int SPHERE_ICOS_FACE_NUM, the number of faces.
//
{
  int face_num;

  face_num = 20 * factor * factor;

  return face_num;
}
//****************************************************************************80

int sphere_icos_point_num ( int factor )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_ICOS_POINT_NUM sizes an icosahedral grid on a sphere.
//
//  Discussion:
//
//    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
//
//    With FACTOR = 2, each triangle of the icosahedron is subdivided into
//    2x2 subtriangles, resulting in 80 faces, 120 edges, and 
//    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
//
//    With FACTOR = 3, each triangle of the icosahedron is subdivided into
//    3x3 subtriangles, resulting in 180 faces, 270 edges and 
//    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
//
//    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
//    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
//      12 
//    + 20 * 3          * (FACTOR-1) / 2 
//    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
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
//    Output, int SPHERE_ICOS_POINT_NUM, the number of points.
//
{
  int point_num;

  point_num = 12                                   
            + 10 * 3              * ( factor - 1 ) 
            + 10 * ( factor - 2 ) * ( factor - 1 );

  return point_num;
}
//****************************************************************************80

double *sphere_icos1_points ( int factor, int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_ICOS1_POINTS returns icosahedral grid points on a sphere.
//
//  Discussion:
//
//    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
//
//    With FACTOR = 2, each triangle of the icosahedron is subdivided into
//    2x2 subtriangles, resulting in 80 faces and 
//    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
//
//    With FACTOR = 3, each triangle of the icosahedron is subdivided into
//    3x3 subtriangles, resulting in 180 faces and 
//    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
//
//    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
//    resulting in 20 * FACTOR * FACTOR faces and
//      12 
//    + 20 * 3          * (FACTOR-1) / 2 
//    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2007
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
//    Input, int NODE_NUM, the number of nodes, as reported
//    by SPHERE_GRID_ICOS_NUM.
//
//    Output, double SPHERE_ICOS1_POINTS[3*NODE_NUM], the node coordinates.
//
//  Local Parameters:
//
//    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
//    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
//    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
//    We need to refer to this data to generate the grid.
//
//    NODE counts the number of nodes we have generated so far.  At the
//    end of the routine, it should be equal to NODE_NUM.
//
{
  int a;
  int b;
  int c;
  int dim;
  int edge;
  int edge_num;
  int *edge_point;
  int f;
  int f1;
  int f2;
  int face;
  int face_num;
  int *face_order;
  int *face_point;
  int face_order_max;
  int node;
  double node_norm;
  double *node_xyz;
  int point;
  double *point_coord;
  int point_num;

  node_xyz = new double[3*node_num];
//
//  Size the icosahedron.
//
  icos_num ( &point_num, &edge_num, &face_num, &face_order_max );
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
//  Generate the point coordinates.
//
//  A.  Points that are the icosahedral vertices.
//
  node = 0;
  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < 3; dim++ )
    {
      node_xyz[dim+node*3] = point_coord[dim+point*3];
    }
    node = node + 1;
  }
//
//  B. Points in the icosahedral edges, at 
//  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
//
  for ( edge = 0; edge < edge_num; edge++ )
  {
    a = edge_point[0+edge*2];
    b = edge_point[1+edge*2];

    for ( f = 1; f < factor; f++ )
    {
      for ( dim = 0; dim < 3; dim++ )
      {
        node_xyz[dim+node*3] = 
          ( ( double ) ( factor - f ) * point_coord[dim+a*3]
          + ( double ) (          f ) * point_coord[dim+b*3] ) 
          / ( double ) ( factor     );
      }

      node_norm = r8vec_norm ( 3, node_xyz+node*3 );

      for ( dim = 0; dim < 3; dim++ )
      {
        node_xyz[dim+node*3] = node_xyz[dim+node*3] / node_norm;
      }
      node = node + 1;
    }
  }
//
//  C.  Points in the icosahedral faces.
//
  for ( face = 0; face < face_num; face++ )
  {
    a = face_point[0+face*3];
    b = face_point[1+face*3];
    c = face_point[2+face*3];

    for ( f1 = 1; f1 < factor; f1++ )
    {
      for ( f2 = 1; f2 < factor - f1; f2++ )
      {
        for ( dim = 0; dim < 3; dim++ )
        {
          node_xyz[dim+node*3] = 
            ( ( double ) ( factor - f1 - f2 ) * point_coord[dim+a*3]  
            + ( double ) (          f1      ) * point_coord[dim+b*3] 
            + ( double ) (               f2 ) * point_coord[dim+c*3] ) 
            / ( double ) ( factor           );
        }
        node_norm = r8vec_norm ( 3, node_xyz+node*3 );

        for ( dim = 0; dim < 3; dim++ )
        {
          node_xyz[dim+node*3] = node_xyz[dim+node*3] / node_norm;
        }
        node = node + 1;
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

  return node_xyz;
}
//****************************************************************************80

double *sphere_icos2_points ( int factor, int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_ICOS2_POINTS returns icosahedral grid points on a sphere.
//
//  Discussion:
//
//    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
//
//    With FACTOR = 2, each triangle of the icosahedron is subdivided into
//    2x2 subtriangles, resulting in 80 faces and 
//    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
//
//    With FACTOR = 3, each triangle of the icosahedron is subdivided into
//    3x3 subtriangles, resulting in 180 faces and 
//    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
//
//    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
//    resulting in 20 * FACTOR * FACTOR faces and
//      12 
//    + 20 * 3          * (FACTOR-1) / 2 
//    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
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
//    Input, int NODE_NUM, the number of nodes, as reported
//    by SPHERE_IMP_GRID_ICOS_NUM.
//
//    Output, double SPHERE_ICOS2_POINTS[3*NODE_NUM], the node coordinates.
//
//  Local Parameters:
//
//    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
//    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
//    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
//    We need to refer to this data to generate the grid.
//
//    NODE counts the number of nodes we have generated so far.  At the
//    end of the routine, it should be equal to NODE_NUM.
//
{
  int a;
  double angle;
  double ab[3];
  double ac[3];
  double acn[3];
  double acn_norm;
  double acp[3];
  int b;
  double bn[3];
  double bn_norm;
  double bp[3];
  int c;
  double cn[3];
  double cn_norm;
  double cp[3];
  int edge;
  int edge_num;
  int *edge_point;
  int f;
  int fa;
  int fbc;
  int face;
  int face_num;
  int *face_order;
  int *face_point;
  int face_order_max;
  int i;
  int j;
  int node;
  double *node_xyz;
  double p;
  int point;
  double *point_coord;
  int point_num;
  double theta;
  double theta_ab;
  double theta_ac;
  double theta_bc;

  node_xyz = new double[3*node_num];
//
//  Size the icosahedron.
//
  icos_num ( &point_num, &edge_num, &face_num, &face_order_max );
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
//  Generate the point coordinates.
//
//  A.  Points that are the icosahedral vertices.
//
  node = 0;
  for ( j = 0; j < point_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_xyz[i+j*3] = point_coord[i+j*3];
    }
    node = node + 1;
  }
//
//  B. Points in the icosahedral edges, at 
//  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
//
  for ( edge = 0; edge < edge_num; edge++ )
  {
    a = edge_point[0+edge*2];

    b = edge_point[1+edge*2];
//
//  Determine the "distance" = angle between points A and B.
//
    theta = sphere_distance_xyz ( point_coord+a*3, point_coord+b*3 );
//
//  Polarize B into BP + BN and normalize BN.
//
    r8vec_polarize ( 3, point_coord+b*3, point_coord+a*3, bn, bp );
    bn_norm = r8vec_norm ( 3, bn );
    for ( i = 0; i < 3; i++ )
    {
      bn[i] = bn[i] / bn_norm;
    }
//
//  March from A to B, by taking equally spaced angles from 0 to THETA.
//  F = 0      => ANGLE = 0     => A
//  F = FACTOR => ANGLE = THETA => B
//
    for ( f = 1; f < factor; f++ )
    {
      angle = ( ( double ) ( f ) * theta ) / ( double ) ( factor );
     
      for ( i = 0; i < 3; i++ )
      {
        node_xyz[i+node*3] = cos ( angle ) * point_coord[i+a*3]
                           + sin ( angle ) * bn[i];
      }
      node = node + 1;
    }
  }
//
//  C.  Points in the icosahedral faces.
//
  for ( face = 0; face < face_num; face++ )
  {
    a = face_point[0+face*3];
    b = face_point[1+face*3];
    c = face_point[2+face*3];
//
//  Determine the "distance" = angle between points A and B, A and C.
//
    theta_ab = sphere_distance_xyz ( point_coord+a*3, point_coord+b*3 );
    theta_ac = sphere_distance_xyz ( point_coord+a*3, point_coord+c*3 );
//
//  Polarize B = BP + BN and normalize BN, C = CP + CN, and normalize CN.
//
    r8vec_polarize ( 3, point_coord+b*3, point_coord+a*3, bn, bp );
    bn_norm = r8vec_norm ( 3, bn );
    for ( i = 0; i < 3; i++ )
    {
      bn[i] = bn[i] / bn_norm;
    }

    r8vec_polarize ( 3, point_coord+c*3, point_coord+a*3, cn, cp );
    cn_norm = r8vec_norm ( 3, cn );
    for ( i = 0; i < 3; i++ )
    {
      cn[i] = cn[i] / cn_norm;
    }
//
//  March AB from A to B:
//    FA = 0      => ANGLE = 0        => AB = A
//    FA = FACTOR => ANGLE = THETA_AB => AB = B
//
//  March AC from A to C:
//    FA = 0      => ANGLE = 0        => AC = A
//    FA = FACTOR => ANGLE = THETA_AC => AC = C
//
    for ( fa = 2; fa < factor; fa++ )
    {
//
//  Determine points AB and AC that use cos ( FA / FACTOR ) of A 
//  and cos ( ( FACTOR - FA ) / FACTOR ) of B or C.
//
      angle = ( (double ) ( fa ) * theta_ab ) / ( double ) ( factor );
      for ( i = 0; i < 3; i++ )
      {
        ab[i] = cos ( angle ) * point_coord[i+a*3] + sin ( angle ) * bn[i];
      }

      angle = ( ( double ) ( fa ) * theta_ac ) / ( double ) ( factor );
      for ( i = 0; i < 3; i++ )
      {
        ac[i] = cos ( angle ) * point_coord[i+a*3] + sin ( angle ) * cn[i];
      }
//
//  Determine the "distance" = angle between points AB and AC.
//
      theta_bc = sphere_distance_xyz ( ab, ac );
//
//  Polarize AC into ACP + ACN and normalize ACN.
//
      r8vec_polarize ( 3, ac, ab, acn, acp );
      acn_norm = r8vec_norm ( 3, acn );
      for ( i = 0; i < 3; i++ )
      {
        acn[i] = acn[i] / acn_norm;
      }
//
//  The interval between AB and AC is broken into FA intervals.
//  Go from 1 to FA - 1.
//
      for ( fbc = 1; fbc < fa; fbc++ )
      {
        angle = ( ( double ) ( fbc ) * theta_bc ) / ( double ) ( fa );
        for ( i = 0; i < 3; i++ )
        {
          node_xyz[i+node*3] = cos ( angle ) * ab[i] + sin ( angle ) * acn[i];
        }
        node = node + 1;
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

  return node_xyz;
}
//****************************************************************************80

int sphere_line_project ( double r, double pc[3], int n, double p[], 
  int maxpnt2, double pp[], double thetamin, double thetamax )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LINE_PROJECT projects a line onto an implicit sphere.
//
//  Discussion:
//
//    The implicit form of a sphere in 3D is:
//
//        pow ( P[0] - PC[0], 2 ) 
//      + pow ( P[1] - PC[1], 2 ) 
//      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
//
//    The line to be projected is specified as a sequence of points.
//    If two successive points subtend a small angle, then the second
//    point is essentially dropped.  If two successive points subtend
//    a large angle, then intermediate points are inserted, so that
//    the projected line stays closer to the sphere.
//
//    Note that if any P coincides with the center of the sphere, then
//    its projection is mathematically undefined.  P will
//    be returned as PC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.  If R is
//    zero, PP will be returned as PC, and if R is
//    negative, points will end up diametrically opposite from where
//    you would expect them for a positive R.
//
//    Input, double PC[3], the coordinates of the center of the sphere.
//
//    Input, int N, the number of points on the line that is
//    to be projected.
//
//    Input, double P[3*N], the coordinates of the points
//    on the line that is to be projected.
//
//    Input, int MAXPNT2, the maximum number of points on the projected
//    line.  Even if the routine thinks that more points are needed,
//    no more than MAXPNT2 will be generated.
//
//    Output, double PP[3*N2], the coordinates of the
//    points representing the projected line.  The value N2 is returned
//    as the function value of this routine.  These points lie on the
//    sphere.  Successive points are separated by at least THETAMIN
//    radians, and by no more than THETAMAX radians.
//
//    Input, double THETAMIN, THETAMAX, the minimum and maximum angular
//    projections allowed between successive projected points.
//    If two successive points on the original line have projections
//    separated by more than THETAMAX radians, then intermediate points
//    will be inserted, in an attempt to keep the line closer to the
//    sphere.  If two successive points are separated by less than
//    THETAMIN radians, then the second point is dropped, and the
//    line from the first point to the next point is considered.
//
//    Output, int SPHERE_LINE_PROJECT, the number of points on 
//    the projected line.  This value can be zero, if the line has an 
//    angular projection of less than THETAMIN radians.
//
{
# define DIM_NUM 3

  double alpha;
  double ang3d;
  double dot;
  int i;
  int j;
  int nfill;
  int n2;
  double tnorm;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double pi[DIM_NUM];
//
//  Check the input.
//
  if ( r == 0.0 )
  {
    n2 = 0;
    return n2;
  }
  r8vec_copy ( DIM_NUM, pc, p1 );
  r8vec_copy ( DIM_NUM, pc, p2 );

  n2 = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( r8vec_eq ( DIM_NUM, p, pc ) )
    {
    }
    else
    {
      r8vec_copy ( DIM_NUM, p2, p1 );

      alpha = sqrt ( pow ( p[0+i*3] - pc[0], 2 )
                   + pow ( p[1+i*3] - pc[1], 2 )
                   + pow ( p[2+i*3] - pc[2], 2 ) );

      p2[0] = pc[0] + r * ( p[0+i*3] - pc[0] ) / alpha;
      p2[1] = pc[1] + r * ( p[1+i*3] - pc[1] ) / alpha;
      p2[2] = pc[2] + r * ( p[2+i*3] - pc[2] ) / alpha;
//
//  If we haven't gotten any points yet, take this point as our start.
//
      if ( n2 == 0 )
      {
        pp[0+n2*3] = p2[0];
        pp[1+n2*3] = p2[1];
        pp[2+n2*3] = p2[2];
        n2 = n2 + 1;
      }
//
//  Compute the angular projection of P1 to P2.
//
      else if ( 1 <= n2 )
      {
        dot = ( p1[0] - pc[0] ) * ( p2[0] - pc[0] ) 
            + ( p1[1] - pc[1] ) * ( p2[1] - pc[1] ) 
            + ( p1[2] - pc[2] ) * ( p2[2] - pc[2] );
        ang3d = arc_cosine (  dot / ( r * r ) );
//
//  If the angle is at least THETAMIN, (or it's the last point),
//  then we will draw a line segment.
//
        if ( thetamin < r8_abs ( ang3d ) || i == n )
        {
//      
//  Now we check to see if the line segment is too long.
//
          if ( thetamax < r8_abs ( ang3d ) )
          {
            nfill = ( int ) ( r8_abs ( ang3d ) / thetamax );

            for ( j = 1; j < nfill; j++ )
            {
              pi[0] = ( ( double ) ( nfill - j ) * ( p1[0] - pc[0] ) 
                      + ( double ) (         j ) * ( p2[0] - pc[0] ) );
              pi[1] = ( ( double ) ( nfill - j ) * ( p1[1] - pc[1] ) 
                      + ( double ) (         j ) * ( p2[1] - pc[1] ) );
              pi[2] = ( ( double ) ( nfill - j ) * ( p1[2] - pc[2] ) 
                      + ( double ) (         j ) * ( p2[2] - pc[2] ) );

              tnorm = r8vec_norm ( DIM_NUM, pi );

              if ( tnorm != 0.0 )
              {
                pi[0] = pc[0] + r * pi[0] / tnorm;
                pi[1] = pc[1] + r * pi[1] / tnorm;
                pi[2] = pc[2] + r * pi[2] / tnorm;
                pp[0+n2*3] = pi[0];
                pp[1+n2*3] = pi[1];
                pp[2+n2*3] = pi[2];
                n2 = n2 + 1;
              }
            }
          }
//
//  Now tack on the projection of point 2.
//
          pp[0+n2*3] = p2[0];
          pp[1+n2*3] = p2[1];
          pp[2+n2*3] = p2[2];
          n2 = n2 + 1;
        }
      }
    }
  }
  return n2;
# undef DIM_NUM
}
//****************************************************************************80

int *sphere_ll_lines ( int nlat, int nlong, int line_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LL_LINES produces lines from a latitude/longitude grid.
//
//  Discussion:
//
//    The point numbering system is the same used in SPHERE_LL_POINTS,
//    and that routine may be used to compute the coordinates of the points.
//
//    The implicit form of a sphere in 3D is:
//
//        pow ( P[0] - PC[0], 2 ) 
//      + pow ( P[1] - PC[1], 2 ) 
//      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LINE_MAX, the maximum number of gridlines.
//
//    Input, int NLAT, NLONG, the number of latitude and longitude
//    lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so NLAT = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Input, int LINE_NUM, the number of grid lines.
//
//    Output, int SPHERE_LL_LINE[2*LINE_NUM], contains pairs of point indices for
//    line segments that make up the grid.
//
{
  int i;
  int j;
  int l;
  int *line;
  int next;
  int newcol;
  int old;

  line = ( int * ) malloc ( 2 * line_num * sizeof ( int ) );
  l = 0;
//
//  "Vertical" lines.
//
  for ( j = 0; j <= nlong - 1; j++ )
  {
    old = 1;
    next = j + 2;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;

    for ( i = 1; i <= nlat-1; i++ )
    {
      old = next;
      next = old + nlong;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }
    old = next;
    line[0+l*2] = old;
    line[1+l*2] = 1 + nlat * nlong + 1;
    l = l + 1;
  }
//
//  "Horizontal" lines.
//
  for ( i = 1; i <= nlat; i++ )
  {
    next = 1 + ( i - 1 ) * nlong + 1;

    for ( j = 0; j <= nlong-2; j++ )
    {
      old = next;
      next = old + 1;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }

    old = next;
    next = 1 + ( i - 1 ) * nlong + 1;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;
  }
//
//  "Diagonal" lines.
//
  for ( j = 0; j <= nlong-1; j++ )
  {
    old = 1;
    next = j + 2;
    newcol = j;

    for ( i = 1; i <= nlat - 1; i++ )
    {
      old = next;
      next = old + nlong + 1;

      newcol = newcol + 1;
      if ( nlong - 1 < newcol )
      {
        newcol = 0;
        next = next - nlong;
      }
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }
  }
  return line;
}
//****************************************************************************80

int sphere_ll_line_num ( int lat_num, int long_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LL_LINE_NUM counts lines for a latitude/longitude grid.
//
//  Discussion:
//
//    The number returned is the number of pairs of points to be connected.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT_NUM, LONG_NUM, the number of latitude and
//    longitude lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so LAT_NUM = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Output, int SPHERE_LL_LINE_NUM, the number of grid lines.
//
{
  int line_num;

  line_num = long_num * ( lat_num + 1 ) 
           + lat_num * long_num 
           + long_num * ( lat_num - 1 );

  return line_num;
}
//****************************************************************************80

double *sphere_ll_points ( double r, double pc[3], int lat_num, int lon_num,
  int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LL_POINTS produces points on a latitude/longitude grid.
//
//  Discussion:
//
//    The implicit form of a sphere in 3D is:
//
//        pow ( P[0] - PC[0], 2 ) 
//      + pow ( P[1] - PC[1], 2 ) 
//      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Input, double PC[3], the coordinates of the center of the sphere.
//
//    Input, int LAT_NUM, LON_NUM, the number of latitude and longitude
//    lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so LAT_NUM = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Input, int POINT_NUM, the number of points.
//
//    Output, double SPHERE_LL_POINTS[3*POINT_NUM], the coordinates 
//    of the grid points.
//
{
  int lat;
  int lon;
  int n;
  double *p;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  p = new double[3*point_num];
  n = 0;
//
//  The north pole.
//
  theta = 0.0;
  phi = 0.0;

  p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
  p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
  p[2+n*3] = pc[2] + r * cos ( phi );
  n = n + 1;
//
//  Do each intermediate ring of latitude.
//
  for ( lat = 1; lat <= lat_num; lat++ )
  {
    phi = ( double ) ( lat ) * pi / ( double ) ( lat_num + 1 );
//
//  Along that ring of latitude, compute points at various longitudes.
//
    for ( lon = 0; lon < lon_num; lon++ )
    {
      theta = ( double ) ( lon ) * 2.0 * pi / ( double ) ( lon_num );

      p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
      p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
      p[2+n*3] = pc[2] + r * cos ( phi );
      n = n + 1;
    }
  }
//
//  The south pole.
//
  theta = 0.0;
  phi = pi;
  p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
  p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
  p[2+n*3] = pc[2] + r * cos ( phi );
  n = n + 1;

  return p;
}
//****************************************************************************80

int sphere_ll_point_num ( int lat_num, int long_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LL_POINT_NUM counts points for a latitude/longitude grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT_NUM, LONG_NUM, the number of latitude 
//    and longitude lines to draw.  The latitudes do not include the North and 
//    South poles, which will be included automatically, so LAT_NUM = 5, for 
//    instance, will result in points along 7 lines of latitude.
//
//    Output, int SPHERE_LL_POINT_NUM, the number of grid points.
//
{
  int point_num;

  point_num = 2 + lat_num * long_num;

  return point_num;
}
//****************************************************************************80

int *sphere_llq_lines ( int nlat, int nlong, int line_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLQ_LINES: latitude/longitude quadrilateral grid lines.
//
//  Discussion:
//
//    The point numbering system is the same used in SPHERE_LL_POINTS,
//    and that routine may be used to compute the coordinates of the points.
//
//    The implicit form of a sphere in 3D is:
//
//        pow ( P[0] - PC[0], 2 ) 
//      + pow ( P[1] - PC[1], 2 ) 
//      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LINE_MAX, the maximum number of gridlines.
//
//    Input, int NLAT, NLONG, the number of latitude and longitude
//    lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so NLAT = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Input, int LINE_NUM, the number of grid lines.
//
//    Output, int SPHERE_LLQ_LINE[2*LINE_NUM], contains pairs of point indices for
//    line segments that make up the grid.
//
{
  int i;
  int j;
  int l;
  int *line;
  int next;
  int newcol;
  int old;

  line = ( int * ) malloc ( 2 * line_num * sizeof ( int ) );
  l = 0;
//
//  "Vertical" lines.
//
  for ( j = 0; j <= nlong - 1; j++ )
  {
    old = 1;
    next = j + 2;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;

    for ( i = 1; i <= nlat-1; i++ )
    {
      old = next;
      next = old + nlong;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }
    old = next;
    line[0+l*2] = old;
    line[1+l*2] = 1 + nlat * nlong + 1;
    l = l + 1;
  }
//
//  "Horizontal" lines.
//
  for ( i = 1; i <= nlat; i++ )
  {
    next = 1 + ( i - 1 ) * nlong + 1;

    for ( j = 0; j <= nlong-2; j++ )
    {
      old = next;
      next = old + 1;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }

    old = next;
    next = 1 + ( i - 1 ) * nlong + 1;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;
  }
  return line;
}
//****************************************************************************80

int sphere_llq_line_num ( int lat_num, int long_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLQ_LINE_NUM counts lines for a latitude/longitude quadrilateral grid.
//
//  Discussion:
//
//    The number returned is the number of pairs of points to be connected.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT_NUM, LONG_NUM, the number of latitude and
//    longitude lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so LAT_NUM = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Output, int SPHERE_LLQ_LINE_NUM, the number of grid lines.
//
{
  int line_num;

  line_num = long_num * ( lat_num + 1 ) 
           + lat_num * long_num;

  return line_num;
}
//****************************************************************************80

double *sphere_spiralpoints ( double r, double pc[3], int n )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_SPIRALPOINTS produces spiral points on an implicit sphere.
//
//  Discussion:
//
//    The points should be arranged on the sphere in a pleasing design.
//
//    The implicit form of a sphere in 3D is:
//
//        pow ( P[0] - PC[0], 2 ) 
//      + pow ( P[1] - PC[1], 2 ) 
//      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
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
//  Reference:
//
//    Edward Saff, Arno Kuijlaars,
//    Distributing Many Points on a Sphere,
//    The Mathematical Intelligencer,
//    Volume 19, Number 1, 1997, pages 5-11.
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Input, double PC[3], the coordinates of the center of the sphere.
//
//    Input, int N, the number of points to create.
//
//    Output, double SPHERE_SPIRALPOINTS[3*N], the coordinates of the grid points.
//
{
  double cosphi = 0.0;
  int i;
  double *p;
  double pi = 3.141592653589793;
  double sinphi = 0.0;
  double theta = 0.0;

  p = new double[3*n];

  for ( i = 0; i < n; i++ )
  {
    cosphi = ( ( double ) ( n - i - 1 ) * ( -1.0 ) 
             + ( double ) (     i     ) * (  1.0 ) ) 
             / ( double ) ( n     - 1 );

    sinphi = sqrt ( 1.0 - cosphi * cosphi );

    if ( i == 0 || i == n - 1 )
    {
      theta = 0.0;
    }
    else
    {
      theta = theta + 3.6 / ( sinphi * sqrt ( ( double ) n ) );
      theta = r8_modp ( theta, 2.0 * pi );
    }
    p[0+i*3] = pc[0] + r * sinphi * cos ( theta );
    p[1+i*3] = pc[1] + r * sinphi * sin ( theta );
    p[2+i*3] = pc[2] + r * cosphi;
  }

  return p;
}
//****************************************************************************80

double *sphere_unit_sample ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_SAMPLE picks a random point on the unit sphere.
//
//  Discussion:
//
//    The unit sphere in 3D satisfies the equation:
//
//      X * X + Y * Y + Z * Z = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points to generate.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double SPHERE_UNIT_SAMPLE[3*N], the sample point.
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
//   this works because the surface area of the sphere between
//  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
//  a patch of area uniformly.
//
    vdot = 2.0 * r8_uniform_01 ( seed ) - 1.0;

    phi = arc_cosine ( vdot );
//
//  Pick a uniformly random rotation between 0 and 2 Pi around the
//  axis of the Z vector.
//
    theta = 2.0 * pi * r8_uniform_01 ( seed );

    x[0+j*3] = cos ( theta ) * sin ( phi );
    x[1+j*3] = sin ( theta ) * sin ( phi );
    x[2+j*3] = cos ( phi );
  }

  return x;
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
# define TIME_NUM 40

  static char time_buffer[TIME_NUM];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_NUM, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_NUM
}
