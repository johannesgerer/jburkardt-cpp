# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "vandermonde.hpp"

//****************************************************************************80

double *bivand1 ( int n, double alpha[], double beta[] )

//****************************************************************************80
//
//  Purpose:
//
//    BIVAND1 returns a bidimensional Vandermonde1 matrix.
//
//  Discussion:
//
//    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
//
//    (x,y)   | (1,10)  (2,10)  (3,10)  (1,20)  (2,20)  (1,30)
//    --------+-----------------------------------------------
//    1       |     1       1       1       1       1       1  
//    x       |     1       2       3       1       2       1
//       y    |    10      10      10      20      20      30
//    x^2     |     1       4       9       1       4       1
//    x  y    |    10      20      30      20      40      30
//    x^2y^2  |   100     100     100     400     400     900
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the data vectors.
//
//    Input, double ALPHA[N], BETA[N], the values that define A.
//
//    Output, double BIVAND1[((N+1)*N)/2*((N+1)*N)/2], the matrix.
//
{
  double *a;
  int e;
  int e1;
  int e2;
  int i1;
  int i2;
  int ii;
  int j1;
  int j2;
  int jj;
  int n2;

  n2 = ( n * ( n + 1 ) ) / 2;
  a = new double[n2*n2];

  e1 = 0;
  e2 = 0;
  e = 0;

  for ( ii = 0; ii < n2; ii++ )
  {
    j1 = 0;
    j2 = 0;
    for ( jj = 0; jj < n2; jj++ )
    {
      if ( ii == 0 )
      {
        a[ii+jj*n2] = 1.0;
      }
      else
      {
        a[ii+jj*n2] = pow ( alpha[j1], e1 ) * pow ( beta[j2], e2 );
      }

      if ( j1 + j2 < n - 1 )
      {
        j1 = j1 + 1;
      }
      else
      {
        j1 = 0;
        j2 = j2 + 1;
      }
    }

    if ( e2 < e )
    {
      e1 = e1 - 1;
      e2 = e2 + 1;
    }
    else
    {
      e = e + 1;
      e1 = e;
      e2 = 0;
    }
  }

  return a;
}
//****************************************************************************80

double *bivand2 ( int n, double alpha[], double beta[] )

//****************************************************************************80
//
//  Purpose:
//
//    BIVAND2 returns a bidimensional Vandermonde1 matrix.
//
//  Discussion:
//
//    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
//
//    (x,y)   | (1,10) (2,10) (3,10) (1,20) (2,20) (3,20) (1,30) (2,30) (3,30)
//    --------+---------------------------------------------------------------
//    1       |     1      1      1      1      1      1      1      1      1  
//    x       |     1      2      3      1      2      1      1      2      3
//    x^2     |     1      4      9      1      4      1      1      4      9
//       y    |    10     10     10     20     20     20     30     30     30
//    x  y    |    10     20     30     20     40     60     30     60     90
//    x^2y    |    10     40     90     20     80    180     30    120    270
//       y^2  |   100    100    100    400    400    400    900    900    900
//    x  y^2  |   100    200    300    400    800   1200    900   1800   2700
//    x^2y^2  |   100    400    900    400   1600   3600    900   3600   8100
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the data vectors.
//
//    Input, double ALPHA[N], BETA[N], the values that define A.
//
//    Output, double BIVAND1[(N*N)*(N*N)], the matrix.
//
{
  double *a;
  int i;
  int ix;
  int iy;
  int j;
  int jx;
  int jy;

  a = new double[n*n*n*n];

  i = 0;
  for ( iy = 0; iy < n; iy++ )
  {
    for ( ix = 0; ix < n; ix++ )
    {
      j = 0;
      for ( jy = 0; jy < n; jy++ )
      {
        for ( jx = 0; jx < n; jx++ )
        {
          a[i+j*n*n] = pow ( alpha[jx], ix ) * pow ( beta[jy], iy );
          j = j + 1;
        }
      }
      i = i + 1;
    }
  }

  return a;
}
//****************************************************************************80

double *dvand ( int n, double alpha[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DVAND solves a Vandermonde system A' * x = b.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ake Bjorck, Victor Pereyra,
//    Solution of Vandermonde Systems of Equations,
//    Mathematics of Computation,
//    Volume 24, Number 112, October 1970, pages 893-903.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double ALPHA[N], the parameters that define the matrix.
//    The values should be distinct.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double DVAND[N], the solution of the linear system.
//
{
  int j;
  int k;
  double *x;

  x = r8vec_copy_new ( n, b );

  for ( k = 0; k < n - 1; k++ )
  {
    for ( j = n - 1; k < j; j-- )
    {
      x[j] = ( x[j] - x[j-1] ) / ( alpha[j] - alpha[j-k-1] );
    }
  }

  for ( k = n - 2; 0 <= k; k-- )
  {
    for ( j = k; j < n - 1; j++ )
    {
      x[j] = x[j] - alpha[k] * x[j+1];
    }
  }

  return x;
}
//****************************************************************************80

void dvandprg ( int n, double alpha[], double b[], double x[], double c[], 
  double m[] )

//****************************************************************************80
//
//  Purpose:
//
//    DVANDPRG solves a Vandermonde system A' * x = f progressively.
//
//  Discussion:
//
//    This function receives the solution to the system of equations A' * x = f
//    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
//    and new values alpha(n) and f(n).  It updates the solution.
//
//    To solve a system of Nbig equations, this function may be called 
//    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
//    current subsystem is returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ake Bjorck, Victor Pereyra,
//    Solution of Vandermonde Systems of Equations,
//    Mathematics of Computation,
//    Volume 24, Number 112, October 1970, pages 893-903.
//
//  Parameters:
//
//    Input, int N, the new order of the matrix, which is 1 
//    larger than on the previous call.  For the first call, N must be 1.
//
//    Input, double ALPHA[N], the parameters that define the matrix.
//    The values should be distinct.  The value ALPHA(N) has just been
//    added to the system.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input/output, double X[N].  On input, the first N-1 entries 
//    contain the solution of the N-1xN-1 linear system.  On output, the 
//    solution to the NxN linear system.
//
//    Input/output, double C[N], M[N].  On input, the first N-1 
//    entries contain factorization data for the N-1xN-1 linear system.  On 
//    output, factorization data for the NxN linear system.
//
{
  double cn;
  int j;
 
  c[n-1] = b[n-1];
  for ( j = n - 1; 1 <= j; j-- )
  {
    c[j-1] = ( c[j] - c[j-1] ) / ( alpha[n-1] - alpha[j-1] );
  }

  if ( n == 1 )
  {
    m[n-1] = 1.0;
  }
  else
  {
    m[n-1] = 0.0;
  }

  cn = c[0];
  x[n-1] = c[0];

  for ( j = n - 1; 1 <= j; j-- )
  {
    m[j] = m[j] - alpha[n-2] * m[j-1];
    x[n-j-1] = x[n-j-1] + m[j] * cn;
  }

  return;
}
//****************************************************************************80

double *pvand ( int n, double alpha[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    PVAND solves a Vandermonde system A * x = b.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ake Bjorck, Victor Pereyra,
//    Solution of Vandermonde Systems of Equations,
//    Mathematics of Computation,
//    Volume 24, Number 112, October 1970, pages 893-903.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double ALPHA[N], the parameters that define the matrix.
//    The values should be distinct.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double PVAND[N], the solution of the linear system.
//
{
  int j;
  int k;
  double *x;

  x = r8vec_copy_new ( n, b );

  for ( k = 0; k < n - 1; k++ )
  {
    for ( j = n - 1; k < j; j-- )
    {
      x[j] = x[j] - alpha[k] * x[j-1];
    }
  }

  for ( k = n - 2; 0 <= k; k-- )
  {
    for ( j = k + 1; j < n; j++ )
    {
      x[j] = x[j] / ( alpha[j] - alpha[j-k-1] );
    }
    for ( j = k; j < n - 1; j++ )
    {
      x[j] = x[j] - x[j+1];
    }
  }

  return x;
}
//****************************************************************************80

void pvandprg ( int n, double alpha[], double b[], double x[], double d[], 
  double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    PVANDPRG solves a Vandermonde system A * x = f progressively.
//
//  Discussion:
//
//    This function receives the solution to the system of equations A * x = f
//    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
//    and new values alpha(n) and f(n).  It updates the solution.
//
//    To solve a system of Nbig equations, this function may be called 
//    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
//    current subsystem is returned.
//
//    Note that the reference, which lists an Algol version of this algorithm, 
//    omits a minus sign, writing
//      u[j] := u[j] x delta;
//    where
//      u[j] := - u[j] x delta;
//    is actually necessary.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ake Bjorck, Victor Pereyra,
//    Solution of Vandermonde Systems of Equations,
//    Mathematics of Computation,
//    Volume 24, Number 112, October 1970, pages 893-903.
//
//  Parameters:
//
//    Input, int N, the new order of the matrix, which is 1 
//    larger than on the previous call.  For the first call, N must be 1.
//
//    Input, double ALPHA[N], the parameters that define the matrix.
//    The values should be distinct.  The value ALPHA(N) has just been
//    added to the system.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Input/output, double X[N]; on input, the solution of the 
//    N-1xN-1 linear system.  On output, the solution of the NxN linear system.
//
//    Input/output, double D[N], U[N]; on input, factorization data 
//    for the N-1xN-1 linear system.  On output, factorization data for the
//    NxN linear system.
//
{
  double delta;
  double dn;
  int j;

  d[n-1] = b[n-1];
  for ( j = n - 1; 1 <= j; j-- )
  {
    d[j-1] = d[j] - alpha[n-j-1] * d[j-1];
  }

  dn = d[0];
  u[n-1] = 1.0;

  for ( j = 1; j <= n - 1; j++ )
  {
    delta = alpha[n-1] - alpha[j-1];
    u[j-1] = - u[j-1] * delta;
    u[n-1] = u[n-1] * delta;
    x[j-1] = x[j-1] + dn / u[j-1];
  }

  x[n-1] = dn / u[n-1];

  return;
}
//****************************************************************************80

double *r8mat_mtv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MTV_NEW multiplies a transposed matrix times a vector.
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
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8MAT_MTV_NEW[N], the product A'*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
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

double *r8vec_copy_new ( int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY_NEW copies an R8VEC to a new R8VEC.
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
//    03 July 2008
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
//    Output, double R8VEC_COPY_NEW[N], the copy of A1.
//
{
  double *a2;
  int i;

  a2 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
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

double *vand1 ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    VAND1 returns the Vandermonde1 matrix A with 1's on the first row.
//
//  Formula:
//
//    A(I,J) = X(J)^(I-1)
//
//  Example:
//
//    N = 5, X = ( 2, 3, 4, 5, 6 )
//
//    1  1   1   1   1
//    2  3   4   5   6
//    4  9  16  25  36
//    8 27  64 125  216
//   16 81 256 625 1296
//
//  Properties:
//
//    A is generally not symmetric: A' /= A.
//
//    A is nonsingular if, and only if, the X values are distinct.
//
//    det ( A ) = product ( 1 <= I <= N ) ( 1 <= J .lt. I ) ( X(I) - X(J) ).
//             = product ( 1 <= J <= N ) X(J)
//             * product ( 1 <= I .lt. J ) ( X(J) - X(I) ).
//
//    A is generally ill-conditioned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Gregory, David Karney,
//    A Collection of Matrices for Testing Computational Algorithms,
//    Wiley, 1969, page 27,
//    LC: QA263.G68.
//
//    Nicholas Higham,
//    Stability analysis of algorithms for solving confluent
//    Vandermonde-like systems,
//    SIAM Journal on Matrix Analysis and Applications,
//    Volume 11, 1990, pages 23-41.
//
//  Parameters:
//
//    Input, int N, the order of the matrix desired.
//
//    Input, double X[N], the values that define A.
//
//    Output, double VAND1[N*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == 0 && x[j] == 0.0 )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = pow ( x[j], i );
      }
    }
  }

  return a;
}
