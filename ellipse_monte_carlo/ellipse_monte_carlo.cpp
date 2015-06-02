# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "ellipse_monte_carlo.hpp"

//****************************************************************************80

void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Jack Dongarra.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Lawson, Hanson, Kincaid, Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( da == 0.0 )
  {
    return;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dy[iy] + da * dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      dy[i] = dy[i] + da * dx[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      dy[i  ] = dy[i  ] + da * dx[i  ];
      dy[i+1] = dy[i+1] + da * dx[i+1];
      dy[i+2] = dy[i+2] + da * dx[i+2];
      dy[i+3] = dy[i+3] + da * dx[i+3];
    }

  }

  return;
}
//****************************************************************************80

double ddot ( int n, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Jack Dongarra.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Lawson, Hanson, Kincaid, Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  double dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
  {
    return dtemp;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dtemp = dtemp + dx[ix] * dy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      dtemp = dtemp + dx[i] * dy[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      dtemp = dtemp + dx[i  ] * dy[i  ]
                    + dx[i+1] * dy[i+1]
                    + dx[i+2] * dy[i+2]
                    + dx[i+3] * dy[i+3]
                    + dx[i+4] * dy[i+4];
    }

  }

  return dtemp;
}
//****************************************************************************80

double ellipse_area1 ( double a[], double r )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_AREA1 returns the area of an ellipse defined by a matrix.
//
//  Discussion:
//
//    The points X in the ellipse are described by a 2 by 2
//    positive definite symmetric matrix A, and a "radius" R, such that
//      X' * A * X <= R * R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[2*2], the matrix that describes
//    the ellipse.  A must be symmetric and positive definite.
//
//    Input, double R, the "radius" of the ellipse.
//
//    Output, double ELLIPSE_AREA1, the area of the ellipse.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  value = r * r * r8_pi / sqrt ( a[0+0*2] * a[1+1*2] - a[1+0*2] * a[0+1*2] );

  return value;
}
//****************************************************************************80

double ellipse_area2 ( double a, double b, double c, double d )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_AREA2 returns the area of an ellipse defined by an equation.
//
//  Discussion:
//
//    The ellipse is described by the formula
//      a x^2 + b xy + c y^2 = d
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, coefficients on the left hand side.
//
//    Input, double D, the right hand side.
//
//    Output, double ELLIPSE_AREA2, the area of the ellipse.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  value = d * r8_pi / ( 4.0 * a * c - b * b );

  return value;
}
//****************************************************************************80

double *ellipse_sample ( int n, double a[], double r, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_SAMPLE samples points in an ellipse.
//
//  Discussion:
//
//    The points X in the ellipsoid are described by a 2 by 2 positive
//    definite symmetric matrix A, and a "radius" R, such that
//      X' * A * X <= R * R
//    The algorithm computes the Cholesky factorization of A:
//      A = U' * U.
//    A set of uniformly random points Y is generated, satisfying:
//      Y' * Y <= R * R.
//    The appropriate points in the ellipsoid are found by solving
//      U * X = Y
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
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
//    Input, int N, the number of points.
//
//    Input, double A[2*2], the matrix that describes the ellipse.
//
//    Input, double R, the right hand side of the ellipse equation.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double ELLIPSE_SAMPLE[2*N], the points.
//
{
  int i;
  int info;
  int j;
  int k;
  static int m = 2;
  double *u;
  double *x;
//
//  Get the upper triangular Cholesky factor U of A.
//
  u = new double[m*m];

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u[i+j*m] = a[i+j*m];
    }
  }

  info = r8po_fa ( u, m, m );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "ELLIPSE_SAMPLE - Fatal error!\n";
    cerr << "  R8PO_FA reports that the matrix A\n";
    cerr << "  is not positive definite symmetric.\n";
    exit ( 1 );
  }
//
//  Get the points Y that satisfy Y' * Y <= R * R.
//
  x = uniform_in_sphere01_map ( m, n, seed );

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
    r8po_sl ( u, m, m, x+j*m );
  }

  delete u;

  return x;
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

int r8po_fa ( double a[], int lda, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_FA factors a real symmetric positive definite matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch and Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the symmetric matrix
//    to be  factored.  Only the diagonal and upper triangle are used.
//    On output, an upper triangular matrix R so that A = R'*R
//    where R' is the transpose.  The strict lower triangle is unaltered.
//    If INFO /= 0, the factorization is not complete.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int R8PO_FA, error flag.
//    0, for normal return.
//    K, signals an error condition.  The leading minor of order K is not
//    positive definite.
//
{
  int info;
  int j;
  int k;
  double s;
  double t;

  for ( j = 1; j <= n; j++ )
  {
    s = 0.0;

    for ( k = 1; k <= j-1; k++ )
    {
      t = a[k-1+(j-1)*lda] - ddot ( k-1, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
      t = t / a[k-1+(k-1)*lda];
      a[k-1+(j-1)*lda] = t;
      s = s + t * t;
    }

    s = a[j-1+(j-1)*lda] - s;

    if ( s <= 0.0 )
    {
      info = j;
      return info;
    }

    a[j-1+(j-1)*lda] = sqrt ( s );
  }

  info = 0;

  return info;
}
//****************************************************************************80

void r8po_sl ( double a[], int lda, int n, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_SL solves a linear system factored by DPOCO or R8PO_FA.
//
//  Discussion:
//
//    To compute inverse(A) * C where C is a matrix with P columns:
//
//      call dpoco ( a, lda, n, rcond, z, info )
//
//      if ( rcond is not too small .and. info == 0 ) then
//        do j = 1, p
//          call r8po_sl ( a, lda, n, c(1,j) )
//        end do
//      end if
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch and Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double A[LDA*N], the output from DPOCO or R8PO_FA.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  int k;
  double t;
//
//  Solve R' * Y = B.
//
  for ( k = 1; k <= n; k++ )
  {
    t = ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
    b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
  }
//
//  Solve R * X = Y.
//
  for ( k = n; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
    t = -b[k-1];
    daxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
  }

  return;
}
//****************************************************************************80

void r8vec_normal_01 ( int n, int &seed, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
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
//    Output, double X[N], a sample of the standard normal PDF.
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
  const double pi = 3.141592653589793;
  double *r;
  int x_hi;
  int x_lo;
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

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );

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
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
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
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );

    delete [] r;
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

double *uniform_in_sphere01_map ( int dim_num, int n, int &seed )

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
//    17 August 2004
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
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double X[DIM_NUM*N], the points.
//
{
  double exponent;
  int i;
  int j;
  double norm;
  double r;
  double *v;
  double *x;
//
  exponent = 1.0 / ( double ) ( dim_num );

  v = new double[dim_num];
  x = new double[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    r8vec_normal_01 ( dim_num, seed, v );
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      norm = norm + pow ( v[i], 2 );
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < dim_num; i++ )
    {
      v[i] = v[i] / norm;
    }
//
//  Now compute a value to map the point ON the sphere INTO the sphere.
//
    r = r8_uniform_01 ( seed );
    r = pow ( r, exponent );

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = r * v[i];
    }
  }

  delete [] v;

  return x;
}
