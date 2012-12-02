# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "condition.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *combin ( double alpha, double beta, int n )

//****************************************************************************80
//
//  Purpose:
//
//    COMBIN returns the COMBIN matrix.
//
//  Discussion:
//
//    This matrix is known as the combinatorial matrix.
//
//  Formula:
//
//    If ( I = J ) then
//      A(I,J) = ALPHA + BETA
//    else
//      A(I,J) = BETA
//
//  Example:
//
//    N = 5, ALPHA = 2, BETA = 3
//
//    5 3 3 3 3
//    3 5 3 3 3
//    3 3 5 3 3
//    3 3 3 5 3
//    3 3 3 3 5
//
//  Properties:
//
//    A is symmetric: A' = A.
//
//    Because A is symmetric, it is normal.
//
//    Because A is normal, it is diagonalizable.
//
//    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
//
//    A is a circulant matrix: each row is shifted once to get the next row.
//
//    det ( A ) = ALPHA^(N-1) * ( ALPHA + N * BETA ).
//
//    A has constant row sums.
//
//    Because A has constant row sums,
//    it has an eigenvalue with this value,
//    and a (right) eigenvector of ( 1, 1, 1, ..., 1 ).
//
//    A has constant column sums.
//
//    Because A has constant column sums,
//    it has an eigenvalue with this value,
//    and a (left) eigenvector of ( 1, 1, 1, ..., 1 ).
//
//    LAMBDA(1:N-1) = ALPHA,
//    LAMBDA(N) = ALPHA + N * BETA.
//
//    The eigenvector associated with LAMBDA(N) is (1,1,1,...,1)/sqrt(N).
//
//    The other N-1 eigenvectors are simply any (orthonormal) basis
//    for the space perpendicular to (1,1,1,...,1).
//
//    A is nonsingular if ALPHA /= 0 and ALPHA + N * BETA /= 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2008
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
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 1, Fundamental Algorithms, Second Edition,
//    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
//
//  Parameters:
//
//    Input, double ALPHA, BETA, scalars that define A.
//
//    Input, int N, the order of the matrix.
//
//    Output, double COMBIN[N*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = alpha + beta;
      }
      else
      {
        a[i+j*n] = beta;
      }
    }
  }
  return a;
}
//****************************************************************************80

double *combin_inverse ( double alpha, double beta, int n )

//****************************************************************************80
//
//  Purpose:
//
//    COMBIN_INVERSE returns the inverse of the COMBIN matrix.
//
//  Formula:
//
//    if ( I = J )
//      A(I,J) = (ALPHA+(N-1)*BETA) / (ALPHA*(ALPHA+N*BETA))
//    else
//      A(I,J) =             - BETA / (ALPHA*(ALPHA+N*BETA))
//
//  Example:
//
//    N = 5, ALPHA = 2, BETA = 3
//
//           14 -3 -3 -3 -3
//           -3 14 -3 -3 -3
//   1/34 *  -3 -3 14 -3 -3
//           -3 -3 -3 14 -3
//           -3 -3 -3 -3 14
//
//  Properties:
//
//    A is symmetric: A' = A.
//
//    Because A is symmetric, it is normal.
//
//    Because A is normal, it is diagonalizable.
//
//    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
//
//    A is a circulant matrix: each row is shifted once to get the next row.
//
//    A is Toeplitz: constant along diagonals.
//
//    det ( A ) = 1 / (ALPHA^(N-1) * (ALPHA+N*BETA)).
//
//    A is well defined if ALPHA /= 0 and ALPHA+N*BETA /= 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 1, Fundamental Algorithms, Second Edition,
//    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
//
//  Parameters:
//
//    Input, double ALPHA, BETA, scalars that define the matrix.
//
//    Input, int N, the order of the matrix.
//
//    Output, double COMBIN_INVERSE[N*N], the matrix.
//
{
  double *a;
  double bot;
  int i;
  int j;

  if ( alpha == 0.0 )
  {
    cerr << "\n";
    cerr << "COMBIN_INVERSE - Fatal error!\n";
    cerr << "  The entries of the matrix are undefined\n";
    cerr << "  because ALPHA = 0.\n";
    exit ( 1 );
  }
  else if ( alpha + n * beta == 0.0 )
  {
    cerr << "\n";
    cerr << "COMBIN_INVERSE - Fatal error!\n";
    cerr << "  The entries of the matrix are undefined\n";
    cerr << "  because ALPHA+N*BETA is zero.\n";
    exit ( 1 );
  }

  a = new double[n*n];

  bot = alpha * ( alpha + ( double ) ( n ) * beta );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = ( alpha + ( double ) ( n - 1 ) * beta ) / bot;
      }
      else
      {
        a[i+j*n] = - beta / bot;
      }
    }
  }
  return a;
}
//****************************************************************************80

double condition_hager ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    CONDITION_HAGER estimates the L1 condition number of a matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Condition Estimates,
//    SIAM Journal on Scientific and Statistical Computing,
//    Volume 5, Number 2, June 1984, pages 311-316.
//
//  Parameters:
//
//    Input, int N, the dimension of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Output, double CONDITION_HAGER, the estimated L1 condition.
//
{
  double *a_lu;
  double *b;
  double c1;
  double c2;
  double cond;
  int i;
  int i1;
  int i2;
  int info;
  int job;
  int *pivot;

  i1 = -1;
  c1 = 0.0;
//
//  Factor the matrix.
//
  a_lu = r8mat_copy_new ( n, n, a );

  pivot = new int[n];

  info = r8ge_fa ( n, a_lu, pivot );

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 1.0 / ( double ) ( n );
  }

  while ( 1 )
  {
    job = 0;
    r8ge_sl ( n, a_lu, pivot, b, job );

    c2 = r8vec_norm_l1 ( n, b );

    for ( i = 0; i < n; i++ )
    {
      b[i] = r8_sign ( b[i] );
    }

    job = 1;
    r8ge_sl ( n, a_lu, pivot, b, job );

    i2 = r8vec_max_abs_index ( n, b );

    if ( 0 <= i1 )
    {
      if ( i1 == i2 || c2 <= c1 )
      {
        break;
      }
    }
    i1 = i2;
    c1 = c2;

    for ( i = 0; i < n; i++ )
    {
      b[i] = 0.0;
    }
    b[i1] = 1.0;
  }
  cond = c2 * r8mat_norm_l1 ( n, n, a );

  delete [] a_lu;
  delete [] b;
  delete [] pivot;

  return cond;
}
//****************************************************************************80

double condition_linpack ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    CONDITION_LINPACK estimates the L1 condition number.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    For the system A * X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    EPSILON * COND.
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
//    Input, int N, the order of the matrix A.
//
//    Input/output, double A[N*N].  On input, a matrix to be factored.
//    On output, the LU factorization of the matrix.
//
//    Output, double CONDITION_LINPACK, the estimated L1 condition.
//
{
  double anorm;
  double cond;
  double ek;
  int i;
  int info;
  int j;
  int k;
  int l;
  int *pivot;
  double s;
  double sm;
  double t;
  double wk;
  double wkm;
  double ynorm;
  double *z;
//
//  Compute the L1 norm of A.
//
  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    s = 0.0;
    for ( i = 0; i < n; i++ )
    {
      s = s + r8_abs ( a[i+j*n] );
    }
    anorm = r8_max ( anorm, s );
  }
//
//  Compute the LU factorization.
//
  pivot = new int[n];

  info = r8ge_fa ( n, a, pivot );

  if ( info != 0 ) 
  {
    delete [] pivot;
    cond = 0.0;
    return cond;
  }
//
//  COND = norm(A) * (estimate of norm(inverse(A)))
//
//  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
//
//  where
//    A * Z = Y
//  and
//    A' * Y = E
//
//  The components of E are chosen to cause maximum local growth in the
//  elements of W, where U'*W = E.  The vectors are frequently rescaled
//  to avoid overflow.
//
//  Solve U' * W = E.
//
  ek = 1.0;
  z = new double[n];
  for ( i = 0; i < n; i++ )
  {
    z[i] = 0.0;
  }

  for ( k = 0; k < n; k++ )
  {
    if ( z[k] != 0.0 ) 
    {
      ek = - r8_sign2 ( ek, z[k] );
    }

    if ( r8_abs ( a[k+k*n] ) < r8_abs ( ek - z[k] ) )
    {
      s = r8_abs ( a[k+k*n] ) / r8_abs ( ek - z[k] );
      for ( i = 0; i < n; i++ )
      {
        z[i] = s * z[i];
      }
      ek = s * ek;
    }

    wk = ek - z[k];
    wkm = -ek - z[k];
    s = r8_abs ( wk );
    sm = r8_abs ( wkm );

    if ( a[k+k*n] != 0.0 )
    {
      wk = wk / a[k+k*n];
      wkm = wkm / a[k+k*n];
    }
    else
    {
      wk = 1.0;
      wkm = 1.0;
    }

    if ( k + 2 <= n )
    {
      for ( j = k+1; j < n; j++ )
      {
        sm = sm + r8_abs ( z[j] + wkm * a[k+j*n] );
        z[j] = z[j] + wk * a[k+j*n];
        s = s + r8_abs ( z[j] );
      }

      if ( s < sm )
      {
        t = wkm - wk;
        wk = wkm;
        for ( j = k+1; j < n; j++ )
        {
          z[j] = z[j] + t * a[k+j*n];
        }
      }
    }
    z[k] = wk;
  }

  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / s;
  }
//
//  Solve L' * Y = W
//
  for ( k = n-1; 0 <= k; k-- )
  {
    for ( i = k+1; i < n; i++ )
    {
      z[k] = z[k] + z[i] * a[i+k*n];
    }
    t = r8_abs ( z[k] );
    if ( 1.0 < t )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i] = z[i] / t;
      }
    }

    l = pivot[k] - 1;

    t    = z[l];
    z[l] = z[k];
    z[k] = t;
  }

  t = 0.0;
  for ( i = 0; i < n; i++ )
  {
    t = t + r8_abs ( z[i] );
  }
  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / t;
  }

  ynorm = 1.0;
//
//  Solve L * V = Y.
//
  for ( k = 0; k < n; k++ )
  {
    l = pivot[k] - 1;

    t    = z[l];
    z[l] = z[k];
    z[k] = t;

    for ( i = k+1; i < n; i++ )
    {
      z[i] = z[i] + t * a[i+k*n];
    }

    t = r8_abs ( z[k] );

    if ( 1.0 < t )
    {
      ynorm = ynorm / t;
      for ( i = 0; i < n; i++ )
      {
        z[i] = z[i] / t;
      }
    }
  }
  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / s;
  }
  ynorm = ynorm / s;
//
//  Solve U * Z = V.
//
  for ( k = n-1; 0 <= k; k-- )
  {
    if ( r8_abs ( a[k+k*n] ) < r8_abs ( z[k] ) )
    {
      s = r8_abs ( a[k+k*n] ) / r8_abs ( z[k] );
      for ( i = 0; i < n; i++ )
      {
        z[i] = s * z[i];
      }
      ynorm = s * ynorm;
    }

    if ( a[k+k*n] != 0.0 )
    {
      z[k] = z[k] / a[k+k*n];
    }
    else
    {
      z[k] = 1.0;
    }

    for ( i = 0; i < k; i++ )
    {
      z[i] = z[i] - a[i+k*n] * z[k];
    }
  }
//
//  Normalize Z in the L1 norm.
//
  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }
  s = 1.0 / s;

  for ( i = 0; i < n; i++ )
  {
    z[i] = s * z[i];
  }
  ynorm = s * ynorm;

  cond = anorm / ynorm;

  delete [] pivot;
  delete [] z;

  return cond;
}
//****************************************************************************80

double condition_sample1 ( int n, double a[], int m )

//****************************************************************************80
//
//  Purpose:
//
//    CONDITION_SAMPLE1 estimates the L1 condition number of a matrix.
//
//  Discussion:
//
//    A naive sampling method is used.
//
//    Only "forward" sampling is used, that is, we only look at results
//    of the form y=A*x.
//
//    Presumably, solving systems A*y=x would give us a better idea of 
//    the inverse matrix.
//
//    Moreover, a power sequence y1 = A*x, y2 = A*y1, ... and the same for
//    the inverse might work better too.
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
//  Parameters:
//
//    Input, int N, the dimension of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Input, int M, the number of samples to use.
//
//    Output, double CONDITION_SAMPLE1, the estimated L1 condition.
//
{
  double a_norm;
  double ainv_norm;
  double *ax;
  double ax_norm;
  double cond;
  int i;
  int seed;
  double *x;
  double x_norm;

  a_norm = 0.0;
  ainv_norm = 0.0;
  seed = 123456789;

  for ( i = 1; i <= m; i++ )
  {
    x = r8vec_uniform_unit_new ( n, seed );
    x_norm = r8vec_norm_l1 ( n, x );
    ax = r8mat_mv_new ( n, n, a, x );
    ax_norm = r8vec_norm_l1 ( n, ax );

    if ( ax_norm == 0.0 )
    {
      cond = 0.0;
      return cond;
    }

    a_norm    = r8_max ( a_norm,    ax_norm / x_norm  );
    ainv_norm = r8_max ( ainv_norm, x_norm  / ax_norm );

    delete [] ax;
    delete [] x;
  }

  cond = a_norm * ainv_norm;

  return cond;
}
//****************************************************************************80

double *conex1 ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX1 returns the CONEX1 matrix.
//
//  Discussion:
//
//    The CONEX1 matrix is a counterexample to the LINPACK condition
//    number estimator RCOND available in the LINPACK routine DGECO.
//
//  Formula:
//
//    1  -1 -2*ALPHA   0
//    0   1    ALPHA    -ALPHA
//    0   1  1+ALPHA  -1-ALPHA
//    0   0  0           ALPHA
//
//  Example:
//
//    ALPHA = 100
//
//    1  -1  -200     0
//    0   1   100  -100
//    0   1   101  -101
//    0   0     0   100
//
//  Properties:
//
//    A is generally not symmetric: A' /= A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Cline, RK Rew,
//    A set of counterexamples to three condition number estimators,
//    SIAM Journal on Scientific and Statistical Computing,
//    Volume 4, 1983, pages 602-611.
//
//  Parameters:
//
//    Input, double ALPHA, the scalar defining A.  
//    A common value is 100.0.
//
//    Output, double CONEX1[4*4], the matrix.
//
{
  double *a;
  int n = 4;

  a = new double[n*n];

  a[0+0*n] = 1.0;
  a[1+0*n] = 0.0;
  a[2+0*n] = 0.0;
  a[3+0*n] = 0.0;

  a[0+1*n] = -1.0;
  a[1+1*n] = 1.0;
  a[2+1*n] = 1.0;
  a[3+1*n] = 0.0;

  a[0+2*n] = -2.0 * alpha;
  a[1+2*n] = alpha;
  a[2+2*n] = 1.0 + alpha;
  a[3+2*n] = 0.0;

  a[0+3*n] = 0.0;
  a[1+3*n] = -alpha;
  a[2+3*n] = -1.0 - alpha;
  a[3+3*n] = alpha;

  return a;
}
//****************************************************************************80

double *conex1_inverse ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX1_INVERSE returns the inverse of the CONEX1 matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the scalar defining A.  
//
//    Output, double CONEX1_INVERSE[4*4], the matrix.
//
{
  double *a;
  int n = 4;

  a = new double[n*n];

  a[0+0*n] =  1.0;
  a[1+0*n] =  0.0;
  a[2+0*n] =  0.0;
  a[3+0*n] =  0.0;

  a[0+1*n] =  1.0 - alpha;
  a[1+1*n] =  1.0 + alpha;
  a[2+1*n] = -1.0;
  a[3+1*n] =  0.0;

  a[0+2*n] =        alpha;
  a[1+2*n] =      - alpha;
  a[2+2*n] =  1.0;
  a[3+2*n] =  0.0;

  a[0+3*n] =  2.0;
  a[1+3*n] =  0.0;
  a[2+3*n] =  1.0 / alpha;
  a[3+3*n] =  1.0 / alpha;

  return a;
}
//****************************************************************************80

double *conex2 ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX2 returns the CONEX2 matrix.
//
//  Formula:
//
//    1   1-1/ALPHA^2  -2
//    0   1/ALPHA      -1/ALPHA
//    0   0             1
//
//  Example:
//
//    ALPHA = 100
//
//    1  0.9999  -2
//    0  0.01    -0.01
//    0  0        1
//
//  Properties:
//
//    A is generally not symmetric: A' /= A.
//
//    A is upper triangular.
//
//    det ( A ) = 1 / ALPHA.
//
//    LAMBDA = ( 1, 1/ALPHA, 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Cline, RK Rew,
//    A set of counterexamples to three condition number estimators,
//    SIAM Journal on Scientific and Statistical Computing,
//    Volume 4, 1983, pages 602-611.
//
//  Parameters:
//
//    Input, double ALPHA, the scalar defining A.  
//    A common value is 100.0.  ALPHA must not be zero.
//
//    Output, double CONEX2[3*3], the matrix.
//
{
  double *a;
  int n = 3;

  a = new double[n*n];

  if ( alpha == 0.0 )
  {
    cerr << "\n";
    cerr << "CONEX2 - Fatal error!\n";
    cerr << "  The input value of ALPHA was zero.\n";
    exit ( 1 );
  }

  a[0+0*n] =  1.0;
  a[1+0*n] =  0.0;
  a[2+0*n] =  0.0;

  a[0+1*n] = ( alpha - 1.0 ) * ( alpha + 1.0 ) / alpha / alpha;
  a[1+1*n] =  1.0 / alpha;
  a[2+1*n] =  0.0;

  a[0+2*n] = -2.0;
  a[1+2*n] = -1.0 / alpha;
  a[2+2*n] =  1.0;

  return a;
}
//****************************************************************************80

double *conex2_inverse ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX2_INVERSE returns the inverse of the CONEX2 matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the scalar defining A.  
//    A common value is 100.0.  ALPHA must not be zero.
//
//    Output, double CONEX2_INVERSE[3*3], the matrix.
//
{
  double *a;
  int n = 3;

  a = new double[n*n];

  if ( alpha == 0.0 )
  {
    cerr << "\n";
    cerr << "CONEX2_INVERSE - Fatal error!\n";
    cerr << "  The input value of ALPHA was zero.\n";
    exit ( 1 );
  }

  a[0+0*n] = 1.0;
  a[1+0*n] = 0.0;
  a[2+0*n] = 0.0;

  a[0+1*n] = ( 1.0 - alpha ) * ( 1.0 + alpha ) / alpha;
  a[1+1*n] = alpha;
  a[2+1*n] = 0.0;

  a[0+2*n] = ( 1.0 + alpha * alpha ) / alpha / alpha;
  a[1+2*n] = 1.0;
  a[2+2*n] = 1.0;

  return a;
}
//****************************************************************************80

double *conex3 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX3 returns the CONEX3 matrix.
//
//  Formula:
//
//    if ( I = J and I < N )
//      A(I,J) =  1.0 for 1<=I<N
//    else if ( I = J = N )
//      A(I,J) = -1.0
//    else if ( J < I )
//      A(I,J) = -1.0
//    else
//      A(I,J) =  0.0
//
//  Example:
//
//    N = 5
//
//     1  0  0  0  0
//    -1  1  0  0  0
//    -1 -1  1  0  0
//    -1 -1 -1  1  0
//    -1 -1 -1 -1 -1
//
//  Properties:
//
//    A is generally not symmetric: A' /= A.
//
//    A is integral, therefore det ( A ) is integral, and 
//    det ( A ) * inverse ( A ) is integral.
//
//    A is lower triangular.
//
//    det ( A ) = -1.
//
//    A is unimodular.
//
//    LAMBDA = ( 1, 1, 1, 1, -1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Cline, RK Rew,
//    A set of counterexamples to three condition number estimators,
//    SIAM Journal on Scientific and Statistical Computing,
//    Volume 4, 1983, pages 602-611.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double CONEX3[N*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( j < i )
      {
        a[i+j*n] = - 1.0;
      }
      else if ( j == i && i != n - 1 )
      {
        a[i+j*n] = 1.0;
      }
      else if ( j == i && i == n - 1 )
      {
        a[i+j*n] = - 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }
  return a;
}
//****************************************************************************80

double *conex3_inverse ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX3_INVERSE returns the inverse of the CONEX3 matrix.
//
//  Example:
//
//    N = 5
//
//     1  0  0  0  0
//     1  1  0  0  0
//     2  1  1  0  0
//     4  2  1  1  0
//    -8 -4 -2 -1 -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Cline, RK Rew,
//    A set of counterexamples to three condition number estimators,
//    SIAM Journal on Scientific and Statistical Computing,
//    Volume 4, 1983, pages 602-611.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double CONEX3_INVERSE[N*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i < n - 1 )
      {
        if ( j < i )
        {
          a[i+j*n] = pow ( 2.0, i - j - 1 );
        }
        else if ( i == j )
        {
          a[i+j*n] = 1.0;
        }
        else
        {
          a[i+j*n] = 0.0;
        }
	  }
      else if ( i == n - 1 )
      {
        if ( j < i )
        {
          a[i+j*n] = - pow ( 2.0, i - j - 1 );
        }
        else
        {
          a[i+j*n] = - 1.0;
        }
      }
    }
  }
  return a;
}
//****************************************************************************80

double *conex4 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX4 returns the CONEX4 matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CONEX4[4*4], the matrix.
//
{
  double *a;
  static double a_save[4*4] = {
     7.0,  6.0,  5.0,  5.0,
    10.0,  8.0,  7.0,  7.0,
     8.0, 10.0,  9.0,  6.0,
     7.0,  9.0, 10.0,  5.0 };

  a = r8mat_copy_new ( 4, 4, a_save );

  return a;
}
//****************************************************************************80

double *conex4_inverse ( )

//****************************************************************************80
//
//  Purpose:
//
//    CONEX4_INVERSE returns the inverse of the CONEX4 matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CONEX4_INVERSE[4*4], the matrix.
//
{
  double *a;
  static double a_save[4*4] = {
   -41.0,  25.0,  10.0, -6.0,
   -17.0,  10.0,   5.0, -3.0,
    10.0,  -6.0,  -3.0,  2.0,
    68.0, -41.0, -17.0, 10.0 };

  a = r8mat_copy_new ( 4, 4, a_save );

  return a;
}
//****************************************************************************80

double *kahan ( double alpha, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    KAHAN returns the KAHAN matrix.
//
//  Formula:
//
//    if ( I = J )
//      A(I,I) =  sin(ALPHA)^(I)
//    else if ( I < J )
//      A(I,J) = - sin(ALPHA)^(I) * cos(ALPHA)
//    else
//      A(I,J) = 0
//
//  Example:
//
//    ALPHA = 0.25, N = 4
//
//    S  -C*S    -C*S      -C*S
//    0     S^2  -C*S^2    -C*S^2
//    0     0       S^3    -C*S^3
//    0     0       0         S^4
//
//    where
//
//      S = sin(ALPHA), C=COS(ALPHA)
//
//  Properties:
//
//    A is upper triangular.
//
//    A = B * C, where B is a diagonal matrix and C is unit upper triangular.
//    For instance, for the case M = 3, N = 4:
//
//    A = | S 0    0    |  * | 1 -C -C  -C |
//        | 0 S^2  0    |    | 0  1 -C  -C |
//        | 0 0    S^3  |    | 0  0  1  -C |
//
//    A is generally not symmetric: A' /= A.
//
//    A has some interesting properties regarding estimation of
//    condition and rank.
//
//    det ( A ) = sin(ALPHA^(N*(N+1)/2).
//
//    LAMBDA(I) = sin ( ALPHA )^I
//
//    A is nonsingular if and only if sin ( ALPHA ) =/= 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nicholas Higham,
//    A survey of condition number estimation for triangular matrices,
//    SIAM Review,
//    Volume 9, 1987, pages 575-596.
//
//    W Kahan,
//    Numerical Linear Algebra,
//    Canadian Mathematical Bulletin,
//    Volume 9, 1966, pages 757-801.
//
//  Parameters:
//
//    Input, double ALPHA, the scalar that defines A.  A typical
//    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
//
//    Input, int M, N, the order of the matrix.
//
//    Output, double KAHAN[M*N], the matrix.
//
{
  double *a;
  double csi;
  int i;
  int j;
  double si;

  a = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    si = pow ( sin ( alpha ), i + 1 );
    csi = - cos ( alpha ) * si;
    for ( j = 0; j < n; j++ )
    {
      if ( j < i )
      {
        a[i+j*m] = 0.0;
      }
      else if ( j == i )
      {
        a[i+j*m] = si;
      }
      else
      {
        a[i+j*m] = csi;
      }
    }
  }
  return a;
}
//****************************************************************************80

double *kahan_inverse ( double alpha, int n )

//****************************************************************************80
//
//  Purpose:
//
//    KAHAN_INVERSE returns the inverse of the KAHAN matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the scalar that defines A.  A typical 
//    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
//
//    Input, int N, the order of the matrix.
//
//    Output, double KAHAN_INVERSE[N*N], the matrix.
//
{
  double *a;
  double ci;
  int i;
  int j;
  double si;

  a = new double[n*n];

  ci = cos ( alpha );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 1.0;
      }
      else if ( i == j - 1 )
      {
        a[i+j*n] = ci;
      }
      else if ( i < j )
      {
        a[i+j*n] = ci * pow ( 1.0 + ci, j - i - 1 );
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }
//
//  Scale the columns.
//
  for ( j = 0; j < n; j++ )
  {
    si = pow ( sin ( alpha ), j + 1 );
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = a[i+j*n] / si;
    }
  }
  return a;
}
//****************************************************************************80

int r8ge_fa ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FA performs a LINPACK-style PLU factorization of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
//
//    The two dimensional array is stored by columns in a one dimensional
//    array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
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
//    N must be positive.
//
//    Input/output, double A[N*N], the matrix to be factored.
//    On output, A contains an upper triangular matrix and the multipliers
//    which were used to obtain it.  The factorization can be written
//    A = L * U, where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Output, int PIVOT[N], a vector of pivot indices.
//
//    Output, int R8GE_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int i;
  int j;
  int k;
  int l;
  double t;
//
  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L, the index of the pivot row.
//
    l = k;

    for ( i = k+1; i <= n; i++ )
    {
      if ( r8_abs ( a[l-1+(k-1)*n] ) < r8_abs ( a[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
//
//  If the pivot index is zero, the algorithm has failed.
//
    if ( a[l-1+(k-1)*n] == 0.0 )
    {
      cout << "\n";
      cout << "R8GE_FA - Fatal error!\n";
      cout << "  Zero pivot on step " << k << "\n";
      exit ( 1 );
    }
//
//  Interchange rows L and K if necessary.
//
    if ( l != k )
    {
      t              = a[l-1+(k-1)*n];
      a[l-1+(k-1)*n] = a[k-1+(k-1)*n];
      a[k-1+(k-1)*n] = t;
    }
//
//  Normalize the values that lie below the pivot entry A(K,K).
//
    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        t              = a[l-1+(j-1)*n];
        a[l-1+(j-1)*n] = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = t;
      }

      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }

    }

  }

  pivot[n-1] = n;

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    cout << "\n";
    cout << "R8GE_FA - Fatal error!\n";
    cout << "  Zero pivot on step " << n << "\n";
    exit ( 1 );
  }

  return 0;
}
//****************************************************************************80

double *r8ge_inverse ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_INVERSE computes the inverse of a R8GE matrix factored by R8GE_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
//    SGEDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, double A[N*N], the factor information computed by R8GE_FA.
//
//    Input, int PIVOT(N), the pivot vector from R8GE_FA.
//
//    Output, double R8GE_INVERSE[N*N], the inverse matrix.
//
{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = new double[n*n];
//
//  Compute Inverse(U).
//
  for ( k = 1; k <= n; k++ )
  {
    for ( i = 1; i <= k-1; i++ )
    {
      b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
    b[k-1+(k-1)*n] = 1.0 / a[k-1+(k-1)*n];

    for ( j = k+1; j <= n; j++ )
    {
      b[k-1+(j-1)*n] = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }
    }
  }
//
//  Multiply Inverse(U) by Inverse(L).
//
  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      b[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = b[i-1+(k-1)*n] + b[i-1+(j-1)*n] * a[j-1+(k-1)*n];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp = b[i-1+(k-1)*n];
        b[i-1+(k-1)*n] = b[i-1+(pivot[k-1]-1)*n];
        b[i-1+(pivot[k-1]-1)*n] = temp;
      }

    }

  }

  return b;
}
//****************************************************************************80

void r8ge_sl ( int n, double a_lu[], int pivot[], double x[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_SL solves a R8GE system factored by R8GE_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    R8GE_SL is a simplified version of the LINPACK routine SGESL.
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
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A_LU[N*N], the LU factors from R8GE_FA.
//
//    Input, int PIVOT[N], the pivot vector from R8GE_FA.
//
//    Input/output, double X[N], on input, the right hand side vector.
//    On output, the solution vector.
//
//    Input, int JOB, specifies the operation.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
{
  int i;
  int k;
  int l;
  double t;
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve PL * Y = B.
//
    for ( k = 1; k <= n-1; k++ )
    {
      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }
      for ( i = k+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] + a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[k-1+(k-1)*n];
      for ( i = 1; i <= k-1; i++ )
      {
        x[i-1] = x[i-1] - a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      t = 0.0;
      for ( i = 1; i <= k-1; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = ( x[k-1] - t ) / a_lu[k-1+(k-1)*n];
    }
//
//  Solve ( PL )' * X = Y.
//
    for ( k = n-1; 1 <= k; k-- )
    {
      t = 0.0;
      for ( i = k+1; i <= n; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = x[k-1] + t;

      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }
    }
  }
  return;
}

