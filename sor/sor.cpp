# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "sor.hpp"

//****************************************************************************80

double *dif2 ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    DIF2 returns the DIF2 matrix.
//
//  Example:
//
//    N = 5
//
//    2 -1  .  .  .
//   -1  2 -1  .  .
//    . -1  2 -1  .
//    .  . -1  2 -1
//    .  .  . -1  2
//
//  Properties:
//
//    A is banded, with bandwidth 3.
//
//    A is tridiagonal.
//
//    Because A is tridiagonal, it has property A (bipartite).
//
//    A is a special case of the TRIS or tridiagonal scalar matrix.
//
//    A is integral, therefore det ( A ) is integral, and 
//    det ( A ) * inverse ( A ) is integral.
//
//    A is Toeplitz: constant along diagonals.
//
//    A is symmetric: A' = A.
//
//    Because A is symmetric, it is normal.
//
//    Because A is normal, it is diagonalizable.
//
//    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
//
//    A is positive definite.
//
//    A is an M matrix.
//
//    A is weakly diagonally dominant, but not strictly diagonally dominant.
//
//    A has an LU factorization A = L * U, without pivoting.
//
//      The matrix L is lower bidiagonal with subdiagonal elements:
//
//        L(I+1,I) = -I/(I+1)
//
//      The matrix U is upper bidiagonal, with diagonal elements
//
//        U(I,I) = (I+1)/I
//
//      and superdiagonal elements which are all -1.
//
//    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
//
//      L(I,I) =    sqrt ( (I+1) / I )
//      L(I,I-1) = -sqrt ( (I-1) / I )
//
//    The eigenvalues are
//
//      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
//                = 4 SIN^2(I*PI/(2*N+2))
//
//    The corresponding eigenvector X(I) has entries
//
//       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
//
//    Simple linear systems:
//
//      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
//
//      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
//
//    det ( A ) = N + 1.
//
//    The value of the determinant can be seen by induction,
//    and expanding the determinant across the first row:
//
//      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
//                = 2 * N - (N-1)
//                = N + 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Gregory, David Karney,
//    A Collection of Matrices for Testing Computational Algorithms,
//    Wiley, 1969,
//    ISBN: 0882756494,
//    LC: QA263.68
//
//    Morris Newman, John Todd,
//    Example A8,
//    The evaluation of matrix inversion programs,
//    Journal of the Society for Industrial and Applied Mathematics,
//    Volume 6, Number 4, pages 466-476, 1958.
//
//    John Todd,
//    Basic Numerical Mathematics,
//    Volume 2: Numerical Algebra,
//    Birkhauser, 1980,
//    ISBN: 0817608117,
//    LC: QA297.T58.
//
//    Joan Westlake,
//    A Handbook of Numerical Matrix Inversion and Solution of 
//    Linear Equations,
//    John Wiley, 1968,
//    ISBN13: 978-0471936756,
//    LC: QA263.W47.
//
//  Parameters:
//
//    Input, int M, N, the order of the matrix.
//
//    Output, double DIF2[M*N], the matrix.
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
      if ( j == i - 1 )
      {
        a[i+j*m] = -1.0;
      }
      else if ( j == i )
      {
        a[i+j*m] = 2.0;
      }
      else if ( j == i + 1 )
      {
        a[i+j*m] = -1.0;
      }
      else
      {
        a[i+j*m] = 0.0;
      }
    }
  }

  return a;
}
//****************************************************************************80

double *sor1 ( int n, double a[], double b[], double x[], double w )

//****************************************************************************80
{
  int i;
  int j;
  double *x_new;

  x_new = new double[n];
//
//  Do the Gauss-Seidel computation.
//
  for ( i = 0; i < n; i++ )
  {
    x_new[i] = b[i];
    for ( j = 0; j < i; j++ )
    {
      x_new[i] = x_new[i] - a[i+j*n] * x_new[j];
    }
    for ( j = i + 1; j < n; j++ )
    {
      x_new[i] = x_new[i] - a[i+j*n] * x[j];
    }
    x_new[i] = x_new[i] / a[i+i*n];
  }
//
//  Use W to blend the Gauss-Seidel update with the old solution.
//
  for ( i = 0; i < n; i++ )
  {
    x_new[i] = ( 1.0 - w ) * x[i] + w * x_new[i];
  }

  return x_new;
}
//****************************************************************************80

double *r8mat_mv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV multiplies a matrix times a vector.
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
//    Output, double R8MAT_MV[M], the product A*X.
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

double r8mat_residual_norm ( int m, int n, double a[], double x[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_RESIDUAL_NORM returns the norm of A*x-b.
//
//  Discussion:
//
//    A is an MxN R8MAT, a matrix of R8's.
//
//    X is an N R8VEC, and B is an M R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2011
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
//    Input, double B[M], the right hand side vector.
//
//    Output, double R8MAT_RESIDUAL_NORM, the norm of A*x-b.
//
{
  int i;
  int j;
  double *r;
  double r_norm;

  r = new double[m];

  for ( i = 0; i < m; i++ )
  {
    r[i] = - b[i];
    for ( j = 0; j < n; j++ )
    {
      r[i] = r[i] + a[i+j*m] * x[j];
    }
  }

  r_norm = 0.0;
  for ( i = 0; i < m; i++ )
  {
    r_norm = r_norm + r[i] * r[i];
  }
  r_norm = sqrt ( r_norm );

  delete [] r;

  return r_norm;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
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

double r8vec_diff_norm_squared ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The square of the L2 norm of the difference of A and B is:
//
//      R8VEC_DIFF_NORM_SQUARED = sum ( 1 <= I <= N ) ( A[I] - B[I] )^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2011
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
//    Output, double R8VEC_DIFF_NORM_SQUARED, the square of the L2 norm of A - B.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
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
