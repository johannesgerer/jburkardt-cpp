# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "ellipsoid_monte_carlo.hpp"

//****************************************************************************80

double *ellipsoid_sample ( int m, int n, double a[], double v[], double r, 
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSOID_SAMPLE samples points uniformly from an ellipsoid.
//
//  Discussion:
//
//    The points X in the ellipsoid are described by a M by M
//    positive definite symmetric matrix A, a "center" V, and 
//    a "radius" R, such that
//
//      (X-V)' * A * (X-V) <= R * R
//
//    The algorithm computes the Cholesky factorization of A:
//
//      A = U' * U.
//
//    A set of uniformly random points Y is generated, satisfying:
//
//      Y' * Y <= R * R.
//
//    The appropriate points in the ellipsoid are found by solving
//
//      U * X = Y
//      X = X + V
//
//    Thanks to Dr Karl-Heinz Keil for pointing out that the original
//    coding was actually correct only if A was replaced by its inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double A[M*M], the matrix that describes
//    the ellipsoid.
//
//    Input, double V[M], the "center" of the ellipsoid.
//
//    Input, double R, the "radius" of the ellipsoid.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double ELLIPSE_SAMPLE[M*N], the points.
//
{
  int i;
  int j;
  double *t;
  double *u;
  double *x;
//
//  Get the Cholesky factor U.
//
  u = r8po_fa ( m, a );

  if ( !u )
  {
    cerr << "\n";
    cerr << "ELLIPSOID_SAMPLE - Fatal error!\n";
    cerr << "  R8PO_FA reports that the matrix A\n";
    cerr << "  is not positive definite symmetric.\n";
    exit ( 1 );
  }
//
//  Get the points Y that satisfy Y' * Y <= 1.
//
  x = uniform_in_sphere01_map ( m, n, seed );
//
//  Get the points Y that satisfy Y' * Y <= R * R.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = r * x[i+j*m];
    }
  }
//
//  Solve U * X = Y.
//
  for ( j = 0; j < n; j++ )
  {    
    t = r8po_sl ( m, u, x + j * m );
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = t[i];
    }
    delete [] t;
  }
//
//  X = X + V.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+j*m] + v[i];
    }
  }

  delete [] u;

  return x;
}
//****************************************************************************80

double ellipsoid_volume ( int m, double a[], double v[], double r )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSOID_VOLUME returns the volume of an ellipsoid.
//
//  Discussion:
//
//    The points X in the ellipsoid are described by an M by M
//    positive definite symmetric matrix A, an M-dimensional point V,
//    and a "radius" R, such that
//      (X-V)' * A * (X-V) <= R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, double A[M*M], the matrix that describes
//    the ellipsoid.  A must be symmetric and positive definite.
//
//    Input, double V[M], the "center" of the ellipse.
//    The value of V is not actually needed by this function.
//
//    Input, double R, the "radius" of the ellipse.
//
//    Output, double ELLIPSOID_VOLUME, the volume of the ellipsoid.
//
{
  int i;
  int info;
  double sqrt_det;
  double *u;
  double volume;

  u = r8po_fa ( m, a );
  
  sqrt_det = 1.0;
  for ( i = 0; i < m; i++ )
  {
    sqrt_det = sqrt_det * u[i+i*m];
  }

  volume = pow ( r, m ) * hypersphere_unit_volume ( m ) / sqrt_det;

  delete [] u;

  return volume;
}
//****************************************************************************80

double hypersphere_unit_volume ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_UNIT_VOLUME: volume of a unit hypersphere in M dimensions.
//
//  Discussion:
//
//    The unit sphere in M dimensions satisfies the equation:
//
//      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
//
//     M  Volume
//
//     1    2
//     2    1        * PI
//     3  ( 4 /   3) * PI
//     4  ( 1 /   2) * PI^2
//     5  ( 8 /  15) * PI^2
//     6  ( 1 /   6) * PI^3
//     7  (16 / 105) * PI^3
//     8  ( 1 /  24) * PI^4
//     9  (32 / 945) * PI^4
//    10  ( 1 / 120) * PI^5
//
//    For the unit sphere, Volume(M) = 2 * PI * Volume(M-2)/ M
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Output, double HYPERSPHERE_UNIT_VOLUME, the volume of the sphere.
//
{
  int i;
  int m2;
  const double r8_pi = 3.141592653589793;
  double volume;

  if ( m % 2== 0 )
  {
    m2 = m / 2;
    volume = 1.0;
    for ( i = 1; i <= m2; i++ )
    {
      volume = volume * r8_pi / ( ( double ) i );
    }
  }
  else
  {
    m2 = ( m - 1 ) / 2;
    volume = pow ( r8_pi, m2 ) * pow ( 2.0, m );
    for ( i = m2 + 1; i <= 2 * m2 + 1; i++ )
    {
      volume = volume / ( ( double ) i );
    }
  }

  return volume;
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
//    08 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points at which the
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

void r8_print ( double r, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PRINT prints an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the value to print.
//
//    Input, string TITLE, a title.
//
{
  cout << title << "  "
       << r << "\n";

  return;
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
//    07 April 2014
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
  int i2lo_hi;
  int i2lo_lo;
  int inc;
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

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

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

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }

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

double *r8po_fa ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_FA factors a R8PO matrix.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//    The positive definite symmetric matrix A has a Cholesky factorization
//    of the form:
//
//      A = R' * R
//
//    where R is an upper triangular matrix with positive elements on
//    its diagonal.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix in R8PO storage.
//
//    Output, double R8PO_FA[N*N], the Cholesky factor in SGE
//    storage, or NULL if there was an error.
//
{
  double *b;
  int i;
  int j;
  int k;
  double s;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = 0; k <= j-1; k++ )
    {
      for ( i = 0; i <= k-1; i++ )
      {
        b[k+j*n] = b[k+j*n] - b[i+k*n] * b[i+j*n];
      }
      b[k+j*n] = b[k+j*n] / b[k+k*n];
    }

    s = b[j+j*n];
    for ( i = 0; i <= j-1; i++ )
    {
      s = s - b[i+j*n] * b[i+j*n];
    }

    if ( s <= 0.0 )
    {
      delete [] b;
      return NULL;
    }

    b[j+j*n] = sqrt ( s );
  }
//
//  Since the Cholesky factor is in R8GE format, zero out the lower triangle.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      b[i+j*n] = 0.0;
    }
  }

  return b;
}
//****************************************************************************80

double *r8po_sl ( int n, double a_lu[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_SL solves a linear system that has been factored by R8PO_FA.
//
//  Discussion:
//
//    The R8PO storage format is appropriate for a symmetric positive definite 
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A_LU[N*N], the Cholesky factor from R8PO_FA.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R8PO_SL[N], the solution vector.
//
{
  int i;
  int k;
  double *x;

  x = new double[n];

  for ( k = 0; k < n; k++ )
  {
    x[k] = b[k];
  }
//
//  Solve R' * y = b.
//
  for ( k = 0; k < n; k++ )
  {
    for ( i = 0; i < k; i++ )
    {
      x[k] = x[k] - x[i] * a_lu[i+k*n];
    }
    x[k] = x[k] / a_lu[k+k*n];
  }
//
//  Solve R * x = y.
//
  for ( k = n-1; 0 <= k; k-- )
  {
    x[k] = x[k] / a_lu[k+k*n];
    for ( i = 0; i < k; i++ )
    {
      x[i] = x[i] - a_lu[i+k*n] * x[k];
    }
  }

  return x;
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

double *r8vec_normal_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.
//
{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  double *x;
  int x_hi;
  int x_lo;

  x = new double[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  If we need just one new value, do that here to avoid null arrays.
//
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    delete [] r;
  }

  return x;
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
  const int i4_huge = 2147483647;
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

double *uniform_in_sphere01_map ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_SPHERE01_MAP maps uniform points into the unit sphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    We first generate a point ON the sphere, and then distribute it
//    IN the sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russell Cheng,
//    Random Variate Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998, pages 168.
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double X[M*N], the points.
//
{
  double exponent;
  int i;
  int j;
  double norm;
  double r;
  double *v;
  double *x;

  exponent = 1.0 / ( double ) ( m );

  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    v = r8vec_normal_01_new ( m, seed );
//
//  Compute the length of the vector.
//
    norm = r8vec_norm ( m, v );
//
//  Normalize the vector.
//
    for ( i = 0; i < m; i++ )
    {
      v[i] = v[i] / norm;
    }
//
//  Now compute a value to map the point ON the sphere INTO the sphere.
//
    r = r8_uniform_01 ( seed );
    r = pow ( r, exponent );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = r * v[i];
    }

    delete [] v;
  }

  return x;
}
