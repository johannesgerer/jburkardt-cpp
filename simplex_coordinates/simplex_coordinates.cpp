# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "simplex_coordinates.hpp"

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

double r8mat_det ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET computes the determinant of an R8MAT.
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
//    08 October 2005
//
//  Author:
//
//    Original FORTRAN77 version by Helmut Spaeth
//    C++ version by John Burkardt
//
//  Reference:
//
//    Helmut Spaeth,
//    Cluster Analysis Algorithms
//    for Data Reduction and Classification of Objects,
//    Ellis Horwood, 1980, page 125-127.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET, the determinant of the matrix.
//
{
  double *b;
  double det;
  int i;
  int j;
  int k;
  int kk;
  int m;
  double temp;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  det = 1.0;

  for ( k = 1; k <= n; k++ )
  {
    m = k;
    for ( kk = k+1; kk <= n; kk++ )
    {
      if ( r8_abs ( b[m-1+(k-1)*n] ) < r8_abs ( b[kk-1+(k-1)*n] ) )
      {
        m = kk;
      }
    }

    if ( m != k )
    {
      det = -det;

      temp = b[m-1+(k-1)*n];
      b[m-1+(k-1)*n] = b[k-1+(k-1)*n];
      b[k-1+(k-1)*n] = temp;
    }

    det = det * b[k-1+(k-1)*n];

    if ( b[k-1+(k-1)*n] != 0.0 )
    {
      for ( i = k+1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / b[k-1+(k-1)*n];
      }

      for ( j = k+1; j <= n; j++ )
      {
        if ( m != k )
        {
          temp = b[m-1+(j-1)*n];
          b[m-1+(j-1)*n] = b[k-1+(j-1)*n];
          b[k-1+(j-1)*n] = temp;
        }
        for ( i = k+1; i <= n; i++ )
        {
          b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * b[k-1+(j-1)*n];
        }
      }
    }
  }

  delete [] b;

  return det;
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

double *r8mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ZERO_NEW returns a new zeroed R8MAT.
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
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
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

double *simplex_coordinates1 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_COORDINATES1 computes the Cartesian coordinates of simplex vertices.
//
//  Discussion:
//
//    The simplex will have its centroid at 0;
//
//    The sum of the vertices will be zero.
//
//    The distance of each vertex from the origin will be 1.
//
//    The length of each edge will be constant.
//
//    The dot product of the vectors defining any two vertices will be - 1 / N.
//    This also means the angle subtended by the vectors from the origin
//    to any two distinct vertices will be arccos ( - 1 / N ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, double SIMPLEX_COORDINATES1[N*(N+1)], the coordinates of the vertices
//    of a simplex in N dimensions.  
//
{
  int i;
  int ii;
  int j;
  double s;
  double *x;

  x = r8mat_zero_new ( n, n + 1 );

  for ( i = 0; i < n; i++ )
  {
//
//  Set X(I,I) so that sum ( X(1:I,I)**2 ) = 1.
//
    s = 0.0;
    for ( ii = 0; ii < i; ii++ )
    {
      s = s + x[ii+i*n] * x[ii+i*n];
    }
    x[i+i*n] = sqrt ( 1.0 - s );
//
//  Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N 
//
    for ( j = i + 1; j < n + 1; j++ )
    {
      s = 0.0;
      for ( ii = 0; ii < i; ii++ )
      {
        s = s + x[ii+i*n] * x[ii+j*n];
      }
      x[i+j*n] = ( - 1.0 / ( double ) ( n ) - s ) / x[i+i*n];
    }
  }

  return x;
}
//****************************************************************************80

double *simplex_coordinates2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_COORDINATES2 computes the Cartesian coordinates of simplex vertices.
//
//  Discussion:
//
//    This routine uses a simple approach to determining the coordinates of
//    the vertices of a regular simplex in n dimensions.
//
//    We want the vertices of the simplex to satisfy the following conditions:
//
//    1) The centroid, or average of the vertices, is 0.
//    2) The distance of each vertex from the centroid is 1.
//       By 1), this is equivalent to requiring that the sum of the squares
//       of the coordinates of any vertex be 1.
//    3) The distance between any pair of vertices is equal (and is not zero.)
//    4) The dot product of any two coordinate vectors for distinct vertices
//       is -1/N; equivalently, the angle subtended by two distinct vertices
//       from the centroid is arccos ( -1/N).
//
//    Note that if we choose the first N vertices to be the columns of the
//    NxN identity matrix, we are almost there.  By symmetry, the last column
//    must have all entries equal to some value A.  Because the square of the
//    distance between the last column and any other column must be 2 (because
//    that's the distance between any pair of columns), we deduce that
//    (A-1)^2 + (N-1)*A^2 = 2, hence A = (1-sqrt(1+N))/N.  Now compute the 
//    centroid C of the vertices, and subtract that, to center the simplex 
//    around the origin.  Finally, compute the norm of one column, and rescale 
//    the matrix of coordinates so each vertex has unit distance from the origin.
//
//    This approach devised by John Burkardt, 19 September 2010.  What,
//    I'm not the first?
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Output, double SIMPLEX_COORDINATES2[N*(N+1)], the coordinates of the vertices
//    of a simplex in N dimensions.  
//
{
  double a;
  double c;
  int i;
  int j;
  double s;
  double *x;

  x = r8mat_zero_new ( n, n + 1 );

  for ( i = 0; i < n; i++ )
  {
    x[i+i*n] = 1.0;
  }

  a = ( 1.0 - sqrt ( 1.0 + ( double ) ( n ) ) ) / ( double ) ( n );

  for ( i = 0; i < n; i++ )
  {
    x[i+n*n] = a;
  }
//
//  Now adjust coordinates so the centroid is at zero.
//
  for ( i = 0; i < n; i++ )
  {
    c = 0.0;
    for ( j = 0; j < n + 1; j++ )
    {
      c = c + x[i+j*n];
    }
    c = c / ( double ) ( n + 1 );
    for ( j = 0; j < n + 1; j++ )
    {
      x[i+j*n] = x[i+j*n] - c;
    }
  }
//
//  Now scale so each column has norm 1.
//
  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + x[i+0*n] * x[i+0*n];
  }
  s = sqrt ( s );

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = x[i+j*n] / s;
    }
  }
  return x;
}
//****************************************************************************80

double simplex_volume ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_VOLUME computes the volume of a simplex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double X[N*(N+1)], the coordinates of the vertices
//    of a simplex in N dimensions.  
//
//    Output, double SIMPLEX_VOLUME, the volume of the simplex.
//
{
  double *a;
  double det;
  int i;
  int j;
  double volume;

  a = new double[n*n];
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = x[i+j*n];
    }
  }
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = a[i+j*n] - x[i+n*n];
    }
  }

  det = r8mat_det ( n, a );

  volume = r8_abs ( det );
  for ( i = 1; i <= n; i++ )
  {
    volume = volume / ( double ) ( i );
  }

  delete [] a;

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
