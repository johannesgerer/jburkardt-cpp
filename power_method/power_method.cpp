# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "power_method.hpp"

//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

double *fibonacci2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI2 returns the FIBONACCI2 matrix.
//
//  Example:
//
//    N = 5
//
//    0 1 0 0 0
//    1 1 0 0 0
//    0 1 1 0 0
//    0 0 1 1 0
//    0 0 0 1 1
//
//  Properties:
//
//    A is generally not symmetric: A' /= A.
//
//    A is tridiagonal.
//
//    Because A is tridiagonal, it has property A (bipartite).
//
//    A is banded, with bandwidth 3.
//
//    A is integral, therefore det ( A ) is integral, and 
//    det ( A ) * inverse ( A ) is integral.
//
//    A is a zero/one matrix.
//
//    If N = 1 then
//      det ( A ) = 0
//    else
//      det ( A ) = -1
//
//    If 1 < N, then A is unimodular.
//
//    When applied to a Fibonacci1 matrix B, the Fibonacci2 matrix
//    A produces the "next" Fibonacci1 matrix C = A*B.
//
//    Let PHI be the golden ratio (1+sqrt(5))/2.
//
//    For 2 <= N, the eigenvalues and eigenvectors are:
//
//    LAMBDA(1)     = PHI,     vector = (1,PHI,PHI^2,...PHI^(N-1));
//    LAMBDA(2:N-1) = 1        vector = (0,0,0,...,0,1);
//    LAMBDA(N)     = 1 - PHI. vector = ((-PHI)^(N-1),(-PHI)^(N-2),...,1)
//
//    Note that there is only one eigenvector corresponding to 1.
//    Hence, for 3 < N, the matrix is defective.  This fact means, 
//    for instance, that the convergence of the eigenvector in the power 
//    method will be very slow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double FIBONACCI2[N*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i == 1 )
      {
        if ( j == 2 )
        {
          a[i-1+(j-1)*n] = 1.0;
        }
        else
        {
          a[i-1+(j-1)*n] = 0.0;
        }
      }
      else
      {
        if ( j == i - 1 || j == i )
        {
          a[i-1+(j-1)*n] = 1.0;
        }
        else
        {
          a[i-1+(j-1)*n] = 0.0;
        }
      }
    }
  }
  return a;
}
//****************************************************************************80

void power_method ( int n, double a[], double y[], int it_max, double tol,
  double *lambda, int *it_num )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_METHOD applies several steps of the power method.
//
//  Discussion:
//
//    For a given NxN matrix A and an N vector Y, the power method produces
//    a series of estimates for LAMBDA, the largest eigenvalue, and Y,
//    the eigenvector corresponding to LAMBDA.
//
//    The iteration repeats the following steps
//
//      AY     = A * Y
//      LAMBDA = || AY ||
//      Y      = AY / LAMBDA
//
//    If the matrix A has a single real eigenvalue of maximum modulus,
//    then this iteration will generally produce a good estimate for that
//    eigenvalue and its corresponding eigenvector.
//
//    If there are multiple distinct eigenvalues of the same modulus,
//    perhaps two values of opposite sign, or complex eigenvalues, then
//    the situation is more complicated.
//
//    Separate issues:
//
//    * when estimating the value of LAMBDA, we use the Rayleigh quotient,
//    LAMBDA = ( y' * A * y ) / ( y' * y ).  Since we normalize Y, the
//    bottom of the fraction is 1.  Using this estimate allows us to
//    easily capture the sign of LAMDBA.  Using the eucldean norm
//    instead, for instance, would always give a positive value.
//
//    * If the dominant eigenvalue is negative, then the iteration 
//    as given will produce eigenvector iterates that alternate in sign.  
//   
//    * It is worth knowing whether the successive eigenvector estimates
//    are tending to some value.  Since an eigenvector is really a direction,
//    we need to normalize the vectors, and we need to somehow treat both
//    a vector and its negative as holding the same information.  This
//    means that the proper way of measuring the difference between two
//    eigenvector estimates is to normalize them both, and then compute
//    the cosine between them as y1'y2, followed by the sine, which is
//    sqrt ( 1 - ( y1'y2)^2 ).  If this sine is small, the vectors y1 and y2
//    are "close" in the sense of direction.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 July 2008
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
//    Input/output, double Y[N], the estimate for the eigenvector.
//
//    Input, int IT_MAX, the maximum number of iterations to take.
//    1 <= IT_MAX.
//
//    Input, double TOL, an error tolerance.
//
//    Output, double *LAMBDA, the estimate for the eigenvalue.
//
//    Output, int *IT_NUM, the number of iterations taken.
//
{
  double *ay;
  double cos_y1y2;
  static bool debug = false;
  int i;
  int it;
  double lambda_old;
  double norm;
  double sin_y1y2;
  double val_dif;
  double *y_old;

  ay = new double[n];
  y_old = new double[n];

  if ( debug )
  {
    cout << "\n";
    cout << "     IT      Lambda          Delta-Lambda    Delta-Y\n";
    cout << "\n";
  }
//
//  Force Y to be a vector of unit norm.
//
  norm = r8vec_norm_l2 ( n, y );

  for ( i = 0; i < n; i++ )
  {
    y[i] = y[i] / norm;
  }
  it = 0;

  for ( i = 0; i < n; i++ )
  {
    y_old[i] = y[i];
  }

  r8mat_mv ( n, n, a, y, ay );
  *lambda = r8vec_dot ( n, y, ay );
  norm = r8vec_norm_l2 ( n, ay );
  for ( i = 0; i < n; i++ )
  {
    y[i] = ay[i] / norm;
  }
  if ( *lambda < 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      y[i] = - y[i];
    }
  }

  val_dif = 0.0;
  cos_y1y2 = r8vec_dot ( n, y, y_old );
  sin_y1y2 = sqrt ( ( 1.0 - cos_y1y2 ) * ( 1.0 + cos_y1y2 ) );

  if ( debug )
  {
    cout << "  " << setw(5) << it
         << "  " << setw(14) << *lambda
         << "  " << setw(14) << val_dif
         << "  " << setw(14) << sin_y1y2 << "\n";
  }

  for ( it = 1; it <= it_max; it++ )
  {
    lambda_old = *lambda;
    for ( i = 0; i < n; i++ )
    {
      y_old[i] = y[i];
    }

    r8mat_mv ( n, n, a, y, ay );
    *lambda = r8vec_dot ( n, y, ay );
    norm = r8vec_norm_l2 ( n, ay );
    for ( i = 0; i < n; i++ )
    {
      y[i] = ay[i] / norm;
    }
    if ( *lambda < 0.0 )
    {
      for ( i = 0; i < n; i++ )
      {
        y[i] = - y[i];
      }
    }

    val_dif = r8_abs ( *lambda - lambda_old );
    cos_y1y2 = r8vec_dot ( n, y, y_old );
    sin_y1y2 = sqrt ( ( 1.0 - cos_y1y2 ) * ( 1.0 + cos_y1y2 ) );

    if ( debug )
    {
      cout << "  " << setw(5) << it
           << "  " << setw(14) << *lambda
           << "  " << setw(14) << val_dif
           << "  " << setw(14) << sin_y1y2 << "\n";
    }

    if ( val_dif <= tol )
    {
      break;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = ay[i] / *lambda;
  }

  *it_num = it;

  delete [] ay;
  delete [] y_old;

  return;
}
//****************************************************************************80

void power_method2 ( int n, double a[], double x_init[], int it_max, 
  double tol, complex <double> *lambda, complex <double> v[], int *it_num )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_METHOD2 applies the power method for possibly complex eigenvalues.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric VanDeVelde,
//    Concurrent Scientific Programming,
//    Springer, 1994,
//    ISBN: 0-387-94195-9,
//    LC: QA76.58.V35.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Input, double X_INIT[N], the initial estimate for the eigenvector.
//
//    Input, int IT_MAX, the maximum number of iterations to take.
//    1 <= IT_MAX.
//
//    Input, double TOL, an error tolerance.
//
//    Output, complex <double> *LAMBDA, the estimate for the eigenvalue.
//
//    Output, complex <double> V[N], the estimate for the eigenvector.
//
//    Output, int *IT_NUM, the number of iterations taken.
//
{
  double alpha;
  double beta;
  double gamma;
  int i;
  int it;
  double lambda_imag;
  double lambda_real;
  double pi_xx;
  double pi_xy;
  double pi_xz;
  double pi_yy;
  double pi_yz;
  double pi_zz;
  double *x;
  double *y;
  double *z;

  x = new double[n];
  y = new double[n];
  z = new double[n];

  *it_num = 0;
//
//  Compute data necessary to start the iteration.
//
  r8vec_copy ( n, x_init, x );

  pi_xx = r8vec_dot ( n, x, x );
  
  r8vec_divide ( n, x, pi_xx );

  r8mat_mv ( n, n, a, x, y );

  pi_xy = r8vec_dot ( n, x, y );
  pi_yy = r8vec_dot ( n, y, y );

  for ( it = 1; it <= it_max; it++ )
  {
    if ( pi_yy - pi_xy * pi_xy < tol * tol * pi_yy )
    {
      *lambda = complex <double> ( pi_xy, 0.0 );
      for ( i = 0; i < n; i++ )
      {
        v[i] = complex <double> ( y[i], 0.0 ) / pi_yy;
      }

      delete [] x;
      delete [] y;
      delete [] z;

      return;
    }
    r8mat_mv ( n, n, a, y, z );

    pi_xz = r8vec_dot ( n, x, z );
    pi_yz = r8vec_dot ( n, y, z );
    pi_zz = r8vec_dot ( n, z, z );

    alpha = - ( pi_yz - pi_xy * pi_xz ) / ( pi_yy - pi_xy * pi_xy );

    beta = ( pi_xy * pi_yz - pi_yy * pi_xz ) / ( pi_yy - pi_xy * pi_xy );

    gamma = pi_zz + alpha * alpha * pi_yy + beta * beta 
      + 2.0 * ( alpha * pi_yz + beta * pi_xz + alpha * beta * pi_xy );

    if ( gamma < tol * tol * pi_zz && alpha * alpha < 4.0 * beta )
    {
      lambda_real = - alpha / 2.0;
      lambda_imag = sqrt ( 4.0 * beta - alpha * alpha ) / 2.0;
      *lambda = complex <double> ( lambda_real, lambda_imag );

      for ( i = 0; i < n; i++ )
      {
        v[i] = ( *lambda * y[i] - z[i] ) 
         / sqrt ( beta * pi_yy + alpha * pi_yz + pi_zz );
      }

      delete [] x;
      delete [] y;
      delete [] z;

      return;
    }
    r8vec_copy ( n, y, x );
    r8vec_divide ( n, x, sqrt ( pi_yy ) );

    r8vec_copy ( n, z, y );
    r8vec_divide ( n, y, sqrt ( pi_yy ) );

    pi_xy = pi_yz / pi_yy;
    pi_yy = pi_zz / pi_yy;

    *it_num = it;
  }

  cout << "\n";
  cout << "POWER_METHOD2 - Fatal error!\n";
  cout << "  Convergence was not reached.\n";

  exit ( 1 );
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
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

void r8mat_mv ( int m, int n, double a[], double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV multiplies a matrix times a vector.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
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
//    Output, double Y[M], the product A*X.
//
{
  int i;
  int j;

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return;
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

void r8vec_divide ( int n, double a[], double s )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIVIDE divides an R8VEC by a nonzero scalar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, double A[N].  On input, the vector to be scaled.
//    On output, each entry has been divided by S.
//
//    Input, double S, the divisor.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / s;
  }
  return;
}
//****************************************************************************80

double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
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
//    Output, double R8VEC_DOT, the dot product of the vectors.
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

double r8vec_norm_l2 ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
//
//  Definition:
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
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double R8VEC_NORM_L2, the L2 norm of A.
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

double *r8vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

double *tris ( int m, int n, double x, double y, double z )

//****************************************************************************80
//
//  Purpose:
//
//    TRIS returns the TRIS matrix.
//
//  Discussion:
//
//    The matrix is a tridiagonal matrix defined by three scalars.
//
//    See page 155 of the Todd reference.
//
//  Formula:
//
//    if ( J = I-1 )
//      A(I,J) = X
//    else if ( J = I )
//      A(I,J) = Y
//    else if ( J = I + 1 )
//      A(I,J) = Z
//    else
//      A(I,J) = 0
//
//  Example:
//
//    M = 5, N = 5, X = 1, Y = 2, Z = 3
//
//    2 3 0 0 0
//    1 2 3 0 0
//    0 1 2 3 0
//    0 0 1 2 3
//    0 0 0 1 2
//
//  Properties:
//
//    A is generally not symmetric: A' /= A.
//
//    A is tridiagonal.
//
//    Because A is tridiagonal, it has property A (bipartite).
//
//    A is banded, with bandwidth 3.
//
//    A is Toeplitz: constant along diagonals.
//
//    If Y is not zero, then for A to be singular, it must be the case that
//
//      0.5 * Y / sqrt ( X * Z ) < 1
//
//    and
//
//      cos (K*PI/(N+1)) = - 0.5 * Y / sqrt ( X * Z ) for some 1 <= K <= N.
//
//    If Y is zero, then A is singular when N is odd, or if X or Z is zero.
//
//    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
//
//    A has eigenvalues
//
//      LAMBDA(I) = Y + 2 * sqrt(X*Z) * COS(I*PI/(N+1))
//
//    The eigenvalues will be complex if X * Z < 0.
//
//    If X = Z, the matrix is symmetric.
//
//    As long as X and Z are nonzero, the matrix is irreducible.
//
//    If X = Z = -1, and Y = 2, the matrix is a symmetric, positive
//    definite M matrix, the negative of the second difference matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Todd,
//    Basic Numerical Mathematics,
//    Volume 2: Numerical Algebra,
//    Birkhauser, 1980,
//    ISBN: 0817608117,
//    LC: QA297.T58.
//
//  Parameters:
//
//    Input, int M, N, the order of the matrix.
//
//    Input, double X, Y, Z, the scalars that define A.
//
//    Output, double TRIS[M*N], the matrix.
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
        a[i+j*m] = x;
      }
      else if ( j == i )
      {
        a[i+j*m] = y;
      }
      else if ( j == i + 1 )
      {
        a[i+j*m] = z;
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

complex <double> *tris_eigenvalues ( int n, double x, double y, double z )

//****************************************************************************80
//
//  Purpose:
//
//    TRIS_EIGENVALUES returns the eigenvalues of the TRIS matrix.
//
//  Discussion:
//
//    The eigenvalues will be complex if X * Z < 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double X, Y, Z, the scalars that define A.
//
//    Output, complex <double> TRIS_EIGENVALUES[N], the eigenvalues.
//
{
  double angle;
  complex <double> arg;
  int i;
  complex <double> *lambda;
  double pi = 3.141592653589793;

  lambda = new complex <double>[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( i + 1 ) * pi / ( double ) ( n + 1 );
    arg = complex <double> ( x * z, 0.0 );
    lambda[i] = y + 2.0 * sqrt ( arg ) * cos ( angle );
  }
  return lambda;
}
