# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>
# include <string>

using namespace std;

# include "cyclic_reduction.hpp"

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

double *r83_cr_fa ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_CR_FA decomposes a real tridiagonal matrix using cyclic reduction.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used to solve
//    linear systems A * x = b.
//
//    R83_CR_FA does not employ pivoting.  Hence, the results can be more
//    sensitive to ill-conditioning than standard Gauss elimination.  In
//    particular, R83_CR_FA will fail if any diagonal element of the matrix
//    is zero.  Other matrices may also cause R83_CR_FA to fail.
//
//    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
//    diagonally dominant, that is, if the absolute value of the diagonal
//    element is strictly greater than the sum of the absolute values of
//    the offdiagonal elements, for each equation.
//
//    The algorithm may be illustrated by the following figures:
//
//    The initial matrix is given by:
//
//          D1 U1
//          L1 D2 U2
//             L2 D3 U3
//                L3 D4 U4
//                   L4 D U5
//                      L5 D6
//
//    Rows and columns are permuted in an odd/even way to yield:
//
//          D1       U1
//             D3    L2 U3
//                D5    L4 U5
//          L1 U2    D2
//             L3 U4    D4
//                L5       D6
//
//    A block LU decomposition is performed to yield:
//
//          D1      |U1
//             D3   |L2 U3
//                D5|   L4 U5
//          --------+--------
//                  |D2'F3
//                  |F1 D4'F4
//                  |   F2 D6'
//
//    For large systems, this reduction is repeated on the lower right hand
//    tridiagonal subsystem until a completely upper triangular system
//    is obtained.  The system has now been factored into the product of a
//    lower triangular system and an upper triangular one, and the information
//    defining this factorization may be used by R83_CR_SL to solve linear
//    systems.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Roger Hockney,
//    A fast direct solution of Poisson's equation using Fourier Analysis,
//    Journal of the ACM,
//    Volume 12, Number 1, pages 95-113, January 1965.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Output, double R83_CR_FA[3*(2*N+1)], factorization information 
//    needed by R83_CR_SL.
//
{
  double *a_cr;
  int iful;
  int ifulp;
  int ihaf;
  int il;
  int ilp;
  int inc;
  int incr;
  int ipnt;
  int ipntp;
  int j;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "R83_CR_FA - Fatal error!\n";
    cout << "  Nonpositive N = " << n << "\n";
    return NULL;
  }

  a_cr = new double[3*(2*n+1)];

  if ( n == 1 )
  {
    a_cr[0+0*3] = 0.0;
    a_cr[0+1*3] = 0.0;
    a_cr[0+2*3] = 0.0;
    a_cr[1+0*3] = 0.0;
    a_cr[1+1*3] = 1.0 / a[1+0*3];
    a_cr[1+2*3] = 0.0;
    a_cr[2+0*3] = 0.0;
    a_cr[2+1*3] = 0.0;
    a_cr[2+2*3] = 0.0;

    return a_cr;
  }
//
//  Zero out the workspace entries.
//
  a_cr[0+0*3] = 0.0;
  for ( j = 1; j <= n-1; j++ )
  {
    a_cr[0+j*3] = a[0+j*3];
  }
  for ( j = n; j <= 2*n; j++ )
  {
    a_cr[0+j*3] = 0.0;
  }

  a_cr[1+0*3] = 0.0;
  for ( j = 1; j <= n; j++ )
  {
    a_cr[1+j*3] = a[1+(j-1)*3];
  }
  for ( j = n+1; j <= 2*n; j++ )
  {
    a_cr[1+j*3] = 0.0;
  }
  a_cr[2+0*3] = 0.0;
  for ( j = 1; j <= n-1; j++ )
  {
    a_cr[2+j*3] = a[2+(j-1)*3];
  }
  for ( j = n; j <= 2*n; j++ )
  {
    a_cr[2+j*3] = 0.0;
  }

  il = n;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    if ( ( il % 2 ) == 1 )
    {
      inc = il + 1;
    }
    else
    {
      inc = il;
    }

    incr = inc / 2;
    il = il / 2;
    ihaf = ipntp + incr + 1;
    ifulp = ipnt + inc + 2;

    for ( ilp = incr; 1 <= ilp; ilp-- )
    {
      ifulp = ifulp - 2;
      iful = ifulp - 1;
      ihaf = ihaf - 1;

      a_cr[1+iful*3] = 1.0 / a_cr[1+iful*3];
      a_cr[2+iful*3]  = a_cr[2+iful*3]  * a_cr[1+iful*3];
      a_cr[0+ifulp*3] = a_cr[0+ifulp*3] * a_cr[1+(ifulp+1)*3];
      a_cr[1+ihaf*3]  = a_cr[1+ifulp*3] 
        - a_cr[0+iful*3]  * a_cr[2+iful*3]
        - a_cr[0+ifulp*3] * a_cr[2+ifulp*3];
      a_cr[2+ihaf*3] = -a_cr[2+ifulp*3] * a_cr[2+(ifulp+1)*3];
      a_cr[0+ihaf*3] = -a_cr[0+ifulp*3] * a_cr[0+(ifulp+1)*3];
    }
  }

  a_cr[1+(ipntp+1)*3] = 1.0 / a_cr[1+(ipntp+1)*3];

  return a_cr;
}
//****************************************************************************80

double *r83_cr_sl ( int n, double a_cr[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_CR_SL solves a real linear system factored by R83_CR_FA.
//
//  Discussion:
//
//    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
//    LU factors of A.  It does so using a form of cyclic reduction.  If
//    the factors computed by R83_CR_FA are passed to R83_CR_SL, then one or 
//    many linear systems involving the matrix A may be solved.
//
//    Note that R83_CR_FA does not perform pivoting, and so the solution 
//    produced by R83_CR_SL may be less accurate than a solution produced 
//    by a standard Gauss algorithm.  However, such problems can be 
//    guaranteed not to occur if the matrix A is strictly diagonally 
//    dominant, that is, if the absolute value of the diagonal coefficient 
//    is greater than the sum of the absolute values of the two off diagonal 
//    coefficients, for each row of the matrix.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 January 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Roger Hockney,
//    A fast direct solution of Poisson's equation using Fourier Analysis,
//    Journal of the ACM,
//    Volume 12, Number 1, pages 95-113, January 1965.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R83_CR_SL[N], the solution.
//
{
  int i;
  int iful;
  int ifulm;
  int ihaf;
  int il;
  int ipnt;
  int ipntp;
  int ndiv;
  double *rhs;
  double *x;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "R83_CR_SL - Fatal error!\n";
    cout << "  Nonpositive N = " << n << "\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x = new double[1];
    x[0] = a_cr[1+1*3] * b[0];
    return x;
  }
//
//  Set up RHS.
//
  rhs = new double[2*n+1];

  rhs[0] = 0.0;
  for ( i = 1; i <= n; i++ )
  {
    rhs[i] = b[i-1];
  }
  for ( i = n+1; i <= 2*n; i++ )
  {
    rhs[i] = 0.0;
  }

  il = n;
  ndiv = 1;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    il = il / 2;
    ndiv = ndiv * 2;
    ihaf = ipntp;

    for ( iful = ipnt + 2; iful <= ipntp; iful = iful + 2 )
    {
      ihaf = ihaf + 1;
      rhs[ihaf] = rhs[iful] 
        - a_cr[2+(iful-1)*3] * rhs[iful-1]
        - a_cr[0+iful*3]     * rhs[iful+1];
    }
  }

  rhs[ihaf] = rhs[ihaf] * a_cr[1+ihaf*3];
  ipnt = ipntp;

  while ( 0 < ipnt )
  {
    ipntp = ipnt;
    ndiv = ndiv / 2;
    il = n / ndiv;
    ipnt = ipnt - il;
    ihaf = ipntp;

    for ( ifulm = ipnt + 1; ifulm <= ipntp; ifulm = ifulm + 2 )
    {
      iful = ifulm + 1;
      ihaf = ihaf + 1;
      rhs[iful] = rhs[ihaf];
      rhs[ifulm] = a_cr[1+ifulm*3] * ( 
                                rhs[ifulm] 
        - a_cr[2+(ifulm-1)*3] * rhs[ifulm-1] 
        - a_cr[0+ifulm*3]     * rhs[iful] );
    }
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = rhs[i+1];
  }

  delete [] rhs;

  return x;
}
//****************************************************************************80

double *r83_cr_sls ( int n, double a_cr[], int nb, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_CR_SLS solves several real linear systems factored by R83_CR_FA.
//
//  Discussion:
//
//    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
//    LU factors of A.  It does so using a form of cyclic reduction.  If
//    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or 
//    many linear systems involving the matrix A may be solved.
//
//    Note that R83_CR_FA does not perform pivoting, and so the solutions
//    produced by R83_CR_SLS may be less accurate than a solution produced 
//    by a standard Gauss algorithm.  However, such problems can be 
//    guaranteed not to occur if the matrix A is strictly diagonally 
//    dominant, that is, if the absolute value of the diagonal coefficient 
//    is greater than the sum of the absolute values of the two off diagonal 
//    coefficients, for each row of the matrix.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Roger Hockney,
//    A fast direct solution of Poisson's equation using Fourier Analysis,
//    Journal of the ACM,
//    Volume 12, Number 1, pages 95-113, January 1965.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.
//
//    Input, int NB, the number of systems.
//
//    Input, double B[N*NB], the right hand sides.
//
//    Output, double R83_CR_SL[N*NB], the solutions.
//
{
  int i;
  int iful;
  int ifulm;
  int ihaf;
  int il;
  int ipnt;
  int ipntp;
  int j;
  int ndiv;
  double *rhs;
  double *x;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "R83_CR_SLS - Fatal error!\n";
    cout << "  Nonpositive N = " << n << "\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x = new double[1*nb];
    for ( j = 0; j < nb; j++ )
    {
      x[0+j*n] = a_cr[1+1*3] * b[0+j*n];
    }
    return x;
  }
//
//  Set up RHS.
//
  rhs = new double[(2*n+1)*nb];

  for ( j = 0; j < nb; j++ )
  {
    rhs[0+j*(2*n+1)] = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      rhs[i+j*(2*n+1)] = b[i-1+j*n];
    }
    for ( i = n + 1; i <= 2 * n; i++ )
    {
      rhs[i+j*(2*n+1)] = 0.0;
    }
  }

  il = n;
  ndiv = 1;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    il = il / 2;
    ndiv = ndiv * 2;

    for ( j = 0; j < nb; j++ )
    {
      ihaf = ipntp;
      for ( iful = ipnt + 2; iful <= ipntp; iful = iful + 2 )
      {
        ihaf = ihaf + 1;
        rhs[ihaf+j*(2*n+1)] = rhs[iful+j*(2*n+1)] 
          - a_cr[2+(iful-1)*3] * rhs[iful-1+j*(2*n+1)]
          - a_cr[0+iful*3]     * rhs[iful+1+j*(2*n+1)];
      }
    }
  }

  for ( j = 0; j < nb; j++ )
  {
    rhs[ihaf+j*(2*n+1)] = rhs[ihaf+j*(2*n+1)] * a_cr[1+ihaf*3];
  }

  ipnt = ipntp;

  while ( 0 < ipnt )
  {
    ipntp = ipnt;
    ndiv = ndiv / 2;
    il = n / ndiv;
    ipnt = ipnt - il;

    for ( j = 0; j < nb; j++ )
    {
      ihaf = ipntp;
      for ( ifulm = ipnt + 1; ifulm <= ipntp; ifulm = ifulm + 2 )
      {
        iful = ifulm + 1;
        ihaf = ihaf + 1;
        rhs[iful+j*(2*n+1)] = rhs[ihaf+j*(2*n+1)];
        rhs[ifulm+j*(2*n+1)] = a_cr[1+ifulm*3] * ( 
                                  rhs[ifulm+j*(2*n+1)] 
          - a_cr[2+(ifulm-1)*3] * rhs[ifulm-1+j*(2*n+1)] 
          - a_cr[0+ifulm*3]     * rhs[iful+j*(2*n+1)] );
      }
    }
  }

  x = new double[n*nb];

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = rhs[i+1+j*(2*n+1)];
    }
  }

  delete [] rhs;

  return x;
}
//****************************************************************************80

void r83_gs_sl ( int n, double a[], double b[], double x[], int it_max, 
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    R83_GS_SL solves a R83 system using Gauss-Seidel iteration.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    This routine simply applies a given number of steps of the
//    iteration to an input approximate solution.  On first call, you can
//    simply pass in the zero vector as an approximate solution.  If
//    the returned value is not acceptable, you may call again, using
//    it as the starting point for additional iterations.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be at least 2.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input/output, double X[N], an approximate solution to the system.
//
//    Input, int IT_MAX, the maximum number of iterations to take.
//
//    Input, int JOB, specifies the system to solve.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
{
  int i;
  int it_num;
//
//  No diagonal matrix entry can be zero.
//
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      cout << "\n";
      cout << "R83_GS_SL - Fatal error!\n";
      cout << "  Zero diagonal entry, index = " << i << "\n";
      exit ( 1 );
    }
  }

  if ( job == 0 )
  {
    for ( it_num = 1; it_num <= it_max; it_num++ )
    {
      x[0] =   ( b[0]                   - a[2+0*3] * x[1]     ) / a[1+0*3];
      for ( i = 1; i < n-1; i++ )
      {
        x[i] = ( b[i] - a[0+i*3] * x[i-1] - a[2+i*3] * x[i+1] ) / a[1+i*3];
      }
      x[n-1] =   ( b[n-1] - a[0+(n-1)*3] * x[n-2]             ) / a[1+(n-1)*3];
    }
  }
  else
  {
    for ( it_num = 1; it_num <= it_max; it_num++ )
    {
      x[0] =   ( b[0]                     - a[0+1*3] * x[1]     ) 
           / a[1+0*3];
      for ( i = 1; i < n-1; i++ )
      {
        x[i] = ( b[i] - a[2+(i-1)*3] * x[i-1] - a[0+(i+1)*3] * x[i+1] ) 
             / a[1+i*3];
      }
      x[n-1] =   ( b[n-1] - a[2+(n-2)*3] * x[n-2]                     ) 
             / a[1+(n-1)*3];
    }
  }
  return;
}
//****************************************************************************80

double *r83_mxv_new ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_MXV_NEW multiplies a R83 matrix times a vector.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[3*N], the R83 matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R83_MXV_NEW[N], the product A * x.
//
{
  double *b;
  int i;

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] =        a[1+i*3] * x[i];
  }
  for ( i = 0; i < n-1; i++ )
  {
    b[i] = b[i] + a[0+(i+1)*3] * x[i+1];
  }
  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + a[2+(i-1)*3] * x[i-1];
  }

  return b;
}
//****************************************************************************80

void r83_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83_PRINT prints a R83 matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[3*N], the R83 matrix.
//
//    Input, string TITLE, a title.
//
{
  r83_print_some ( n, a, 1, 1, n, n, title );

  return;
}
//****************************************************************************80

void r83_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83_PRINT_SOME prints some of a R83 matrix.
//
//  Discussion:
//
//    The R83 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a R83 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2006
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
//    Input, double A[3*N], the R83 matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column, to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      cout << setw(7) << j << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - 1 );

    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + 1 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(6) << i << ": ";

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( 1 < i - j || 1 < j - i )
        {
          cout << "              ";
        }
        else if ( j == i + 1 )
        {
          cout << setw(12) << a[0+(j-1)*3] << "  ";
        }
        else if ( j == i )
        {
          cout << setw(12) << a[1+(j-1)*3] << "  ";
        }
        else if ( j == i - 1 )
        {
          cout << setw(12) << a[2+(j-1)*3] << "  ";
        }

      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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
//    10 September 2009
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
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << ": ";
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

double *r8vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDICATOR_NEW sets an R8VEC to the indicator vector {1,2,3...}.
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
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, double R8VEC_INDICATOR_NEW[N], the indicator array.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i <= n-1; i++ ) 
  {
    a[i] = ( double ) ( i + 1 );
  }

  return a;
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

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
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
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, integer I_LO, I_HI, the first and last indices to print.
//    The routine expects 1 <= I_LO <= I_HI <= N.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i       
         << ": " << setw(14) << a[i-1]  << "\n";
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
