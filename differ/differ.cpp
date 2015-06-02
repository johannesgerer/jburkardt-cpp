# include <cstdlib>
# include <iostream>
# include <sstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "differ.hpp"

//****************************************************************************80

void differ_backward ( double h, int o, int p, double c[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFER_BACKWARD computes backward difference coefficients.
//
//  Discussion:
//
//    We determine coefficients C to approximate the derivative at X0
//    of order O and precision P, using equally spaced backward
//    differences, so that 
//
//      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x-ih) + O(h^(p))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double H, the spacing.  0 < H.
//
//    Input, int O, the order of the derivative to be 
//    approximated.  1 <= O.
//
//    Input, int P, the order of the error, as a power of H.
//
//    Output, double C[O+P], the coefficients.
//
//    Output, double X[O+P], the evaluation points.
//
{
  double *b;
  int i;
  int info;
  int job;
  int n;
  double t;

  n = o + p;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 - n ) * h;
  }

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }
  b[o] = 1.0;

  job = 0;
  r8vm_sl ( n, x, b, job, c, info );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "DIFFER_BACKWARD - Fatal error!\n";
    cerr << "  Vandermonde linear system is singular.\n";
    exit ( 1 );
  }

  t = r8_factorial ( o );
  for ( i = 0; i < n; i++ )
  {
    c[i] = c[i] * t;
  }

  delete [] b;

  return;
}
//****************************************************************************80

void differ_central ( double h, int  o, int p, double c[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFER_CENTRAL computes central difference coefficients.
//
//  Discussion:
//
//    We determine coefficients C to approximate the derivative at X0
//    of order O and precision P, using equally spaced central
//    differences, so that 
//
//      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+(2*i-o-p+1)*h/2) 
//        + O(h^(p))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double H, the spacing.  0 < H.
//
//    Input, int O, the order of the derivative to 
//    be approximated.  1 <= O.
//
//    Input, int P, the order of the error, as a power of H.
//
//    Output, double C[O+P], the coefficients.
//
//    Output, double X[O+P], the evaluation points.
//
{
  double *b;
  int i;
  int info;
  int job;
  int n;
  double t;

  n = o + p;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) * h / 2.0;
  }

  b = new double[n];
  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }
  b[o] = 1.0;

  job = 0;
  r8vm_sl ( n, x, b, job, c, info );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "DIFFER_CENTRAL - Fatal error!\n";
    cerr << "  Vandermonde linear system is singular.\n";
    exit ( 1 );
  }

  t = r8_factorial ( o );
  for ( i = 0; i < n; i++ )
  {
    c[i] = c[i] * t;
  }

  delete [] b;

  return;
}
//****************************************************************************80

void differ_forward ( double h, int o, int p, double c[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFER_FORWARD computes forward difference coefficients.
//
//  Discussion:
//
//    We determine coefficients C to approximate the derivative at X0
//    of order O and precision P, using equally spaced forward
//    differences, so that 
//
//      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+ih) + O(h^(p))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real H, the spacing.  0 < H.
//
//    Input, integer O, the order of the derivative to be approximated.
//    1 <= O.
//
//    Input, integer P, the order of the error, as a power of H.
//
//    Output, real C[O+P], the coefficients.
//
//    Output, real X[O+P], the evaluation points.
//
{
  double *b;
  int i;
  int info;
  int job;
  int n;
  double t;

  n = o + p;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i ) * h;
  }

  b = new double[n];
  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }
  b[o] = 1.0;

  job = 0;
  r8vm_sl ( n, x, b, job, c, info );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "DIFFER_FORWARD - Fatal error!\n";
    cerr << "  Vandermonde linear system is singular.\n";
    exit ( 1 );
  }

  t = r8_factorial ( o );
  for ( i = 0; i < n; i++ )
  {
    c[i] = c[i] * t;
  }

  delete [] b;

  return;
}
//****************************************************************************80

double *differ_inverse ( int n, double stencil[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFER_INVERSE returns the inverse of the DIFFER matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double STENCIL[N], the values that define A.
//
//    Output, double DIFFER_INVERSE[N*N], the matrix.
//
{
  double *a;
  int i;
  int indx;
  int j;
  int k;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( j == 0 )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  for ( i = 0; i < n; i++ )
  {
    indx = 0;

    for ( k = 0; k < n; k++ )
    {
      if ( k != i )
      {
        for ( j = indx + 1; 0 <= j; j-- )
        {
          a[i+j*n] = - stencil[k] * a[i+j*n] / ( stencil[i] - stencil[k] );

          if ( 0 < j )
          {
            a[i+j*n] = a[i+j*n] + a[i+(j-1)*n] / ( stencil[i] - stencil[k] );
          }
        }
        indx = indx + 1;
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = a[i+j*n] / stencil[i];
    }
  }

  return a;
}
//****************************************************************************80

double *differ_matrix ( int n, double stencil[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFER_MATRIX computes the stencil matrix from the stencil vector.
//
//  Discussion:
//
//    If N = 4, and STENCIL = ( -3, -2, -1, 1 ), then A will be
//
//    -3  -2  -1  1
//     9   4   1  1
//   -27  -8  -1  1
//    81  16   1  1
//
//    This matrix is a generalized form of a Vandermonde matrix A2:
//
//     1   1   1  1
//    -3  -2  -1  1
//     9   4   1  1
//   -27  -8  -1  1    
//
//    and if A * x = b, the A2 * x2 = b, where x2(i) = x(i) * stencil(i)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of stencil points.
//
//    Input, double STENCIL[N], the stencil vector.
//    The entries in this vector must be distinct.
//    No entry of STENCIL may be 0.
//
//    Output, double DIFFER_MATRIX[N*N], the stencil matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    a[0+j*n] = stencil[j];
    for ( i = 1; i < n; i++ )
    {
      a[i+j*n] = a[i-1+j*n] * stencil[j];
    }
  }

  return a;
}
//****************************************************************************80

double *differ_solve ( int n, double stencil[], int order )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFER_SOLVE solves for finite difference coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of stencil points.
//
//    Input, double STENCIL[N], the stencil vector.
//    The entries in this vector must be distinct.
//    No entry of STENCIL may be 0.
//
//    Input, int ORDER, the order of the derivative to
//    be approximated.  1 <= ORDER <= N.
//
//    Output, double DIFFER_SOLVE[N], the coefficients to be used
//    to multiply U(STENCIL(I))-U(0), so that the sum forms an
//    approximation to the derivative of order ORDER, with error 
//    of order H^(N+1-ORDER).
//
{ 
  double *a;
  double *b;
  double *c;
  int i;

  a = differ_matrix ( n, stencil );

  b = new double[n];
  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }
  b[order-1] = 1.0;
//
//  Solve A * C = B.
//
  c = r8mat_fs_new ( n, a, b );

  delete [] a;
  delete [] b;

  return c;
}
//****************************************************************************80

void differ_stencil ( double x0, int o, int p, double x[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFER_STENCIL computes finite difference coefficients.
//
//  Discussion:
//
//    We determine coefficients C to approximate the derivative at X0
//    of order O and precision P, using finite differences, so that 
//
//      d^o f(x)/dx^o (x0) = sum ( 0 <= i <= o+p-1 ) c(i) f(x(i)) 
//        + O(h^(p))
//
//    where H is the maximum spacing between X0 and any X(I).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X0, the point where the derivative is to 
//    be approximated.
//
//    Input, int O, the order of the derivative to be 
//    approximated.  1 <= O.
//
//    Input, int P, the order of the error, as a power of H.
//
//    Input, double X[O+P], the evaluation points.
//
//    Output, double C[O+P], the coefficients.
//
{
  double *b;
  double *dx;
  int i;
  int info;
  int job;
  int n;
  double t;

  n = o + p;

  dx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    dx[i] = x[i] - x0;
  }

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }
  b[o] = 1.0;

  job = 0;
  r8vm_sl ( n, dx, b, job, c, info );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "DIFFER_STENCIL - Fatal error!\n";
    cerr << "  Vandermonde linear system is singular.\n";
    exit ( 1 );
  }

  t = r8_factorial ( o );
  for ( i = 0; i < n; i++ )
  {
    c[i] = c[i] * t;
  }

  delete [] b;
  delete [] dx;

  return;
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

string i4_to_string ( int i4 )

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
//    16 January 2013
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
  ostringstream fred;
  string value;

  fred << i4;

  value = fred.str ( );

  return value;
}
//****************************************************************************80

double inverse_error ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_ERROR determines the error in an inverse matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Input, double B[N*N], the inverse.
//
//    Output, double ERROR_FROBENIUS, the Frobenius norm
//    of (A*B-I) + (B*A-I).
//
{
  double *c;
  int j;
  double value;

  c = r8mat_mm_new ( n, n, n, a, b );

  for ( j = 0; j < n; j++ )
  {
    c[j+j*n] = c[j+j*n] - 1.0;
  }

  value = r8mat_norm_fro ( n, n, c );

  delete [] c;

  c = r8mat_mm_new ( n, n, n, b, a );

  for ( j = 0; j < n; j++ )
  {
    c[j+j*n] = c[j+j*n] - 1.0;
  }

  value = value + r8mat_norm_fro ( n, n, c );

  delete [] c;

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

double *r8mat_fs_new ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_FS_NEW factors and solves a system with one right hand side.
//
//  Discussion:
//
//    This routine differs from R8MAT_FSS_NEW in two ways:
//    * only one right hand side is allowed;
//    * the input matrix A is not modified.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double X[N], the solution of the linear system.
//
{
  double *a2;
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  a2 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a2[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a2[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a2[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cout << "\n";
      cout << "R8MAT_FS_NEW - Fatal error!\n";
      cout << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                  = a2[jcol-1+(j-1)*n];
        a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
        a2[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
//
//  Scale the pivot row.
//
    t = a2[jcol-1+(jcol-1)*n];
    a2[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a2[jcol-1+(j-1)*n] = a2[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a2[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a2[i-1+(jcol-1)*n];
        a2[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a2[i-1+(j-1)*n] = a2[i-1+(j-1)*n] + t * a2[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a2[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  delete [] a2;

  return x;
}
//****************************************************************************80

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
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
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM_NEW[N1*N3], the product matrix C = A * B.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
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
//    Output, double R8MAT_MV_NEW[M], the product A*X.
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

double r8mat_norm_fro ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
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
//    Input, double A[M*N], the matrix whose Frobenius
//    norm is desired.
//
//    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
//
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
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
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
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
//    26 June 2013
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
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }
    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8mat_sub_new ( int m, int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SUB_NEW computes C = A - B.
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
//    30 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the order of the matrices.
//
//    Input, double A[M*N], double B[M*N], the matrices.
//
//    Output, double R8MAT_SUB_NEW[M*N], the value of A-B.
//
{
  double *c;
  int i;
  int j;

  c = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*m] = a[i+j*m] - b[i+j*m];
    }
  }

  return c;
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

double *r8vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
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
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void r8vec2_print ( int n, double a1[], double a2[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_PRINT prints an R8VEC2.
//
//  Discussion:
//
//    An R8VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A1[N], double A2[N], the vectors to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(6)  << i
         << ": " << setw(14) << a1[i]
         << "  " << setw(14) << a2[i] << "\n";
  }

  return;
}
//****************************************************************************80

void r8vm_sl ( int n, double a[], double b[], int job, double x[], int &info )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_SL solves a R8VM system.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//    Vandermonde systems are very close to singularity.  The singularity
//    gets worse as N increases, and as any pair of values defining
//    the matrix get close.  Even a system as small as N = 10 will
//    involve the 9th power of the defining values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Golub, VanLoan.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Gene Golub, Charles Van Loan,
//    Matrix Computations,
//    Third Edition,
//    Johns Hopkins, 1996.
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, int &INFO.
//    0, no error.
//    nonzero, at least two of the values in A are equal.
//
//    Output, double X[N], the solution of the linear system.
//
{
  int i;
  int j;
//
//  Check for explicit singularity.
//
  info = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j+1; i < n; i++ )
    {
      if ( a[i] == a[j] )
      {
        info = 1;
        return;
      }
    }
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = x[i-1] - a[j-1] * x[i-2];
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] / ( a[i-1] - a[i-j-1] );
      }

      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i];
      }
    }
  }
  else
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = ( x[i-1] - x[i-2] ) / ( a[i-1] - a[i-j-1] );
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i] * a[j-1];
      }
    }

  }

  return;
}
//****************************************************************************80

double *r8vm_sl_new ( int n, double a[], double b[], int job, int &info )

//****************************************************************************80
//
//  Purpose:
//
//    R8VM_SL_NEW solves a R8VM system.
//
//  Discussion:
//
//    The R8VM storage format is used for an M by N Vandermonde matrix.
//    An M by N Vandermonde matrix is defined by the values in its second
//    row, which will be written here as X(1:N).  The matrix has a first 
//    row of 1's, a second row equal to X(1:N), a third row whose entries
//    are the squares of the X values, up to the M-th row whose entries
//    are the (M-1)th powers of the X values.  The matrix can be stored
//    compactly by listing just the values X(1:N).
//
//    Vandermonde systems are very close to singularity.  The singularity
//    gets worse as N increases, and as any pair of values defining
//    the matrix get close.  Even a system as small as N = 10 will
//    involve the 9th power of the defining values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Golub, VanLoan.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Gene Golub, Charles Van Loan,
//    Matrix Computations,
//    Third Edition,
//    Johns Hopkins, 1996.
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N], the R8VM matrix.
//
//    Input, double B[N], the right hand side.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, int &INFO.
//    0, no error.
//    nonzero, at least two of the values in A are equal.
//
//    Output, double R8VM_SL_NEW[N], the solution of the linear system.
//
{
  int i;
  int j;
  double *x;
//
//  Check for explicit singularity.
//
  info = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j+1; i < n; i++ )
    {
      if ( a[i] == a[j] )
      {
        info = 1;
        return NULL;
      }
    }
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = x[i-1] - a[j-1] * x[i-2];
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] / ( a[i-1] - a[i-j-1] );
      }

      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i];
      }
    }
  }
  else
  {
    for ( j = 1; j <= n-1; j++ )
    {
      for ( i = n; j+1 <= i; i-- )
      {
        x[i-1] = ( x[i-1] - x[i-2] ) / ( a[i-1] - a[i-j-1] );
      }
    }

    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j; i <= n-1; i++ )
      {
        x[i-1] = x[i-1] - x[i] * a[j-1];
      }
    }

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
