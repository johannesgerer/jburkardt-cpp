# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "linpack_d.hpp"
# include "blas1_d.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test29 ( );
void test30 ( );
void test31 ( );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    LINPACK_D_PRB tests the double precision real LINPACK routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "LINPACK_D_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LINPACK_D library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LINPACK_D_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests DCHDC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define LDA N

  double a[LDA*N];
  double b[LDA*N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  double work[N];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For a general matrix,\n";
  cout << "  DCHDC computes the Cholesky decomposition.\n";
  cout << "\n";
  cout << "  The number of equations is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= N; i++ )
    {
      if ( i == j-1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else if ( i == j )
      {
        a[i-1+(j-1)*LDA] = 2.0;
      }
      else if ( i == j+1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*LDA] = 0.0;
      }
    }
  }

  cout << "\n";
  cout << "  The matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }
//
//  Decompose the matrix.
//
  cout << "\n";
  cout << "  Decompose the matrix.\n";

  job = 0;

  for ( i = 0; i < N; i++ )
  {
    ipvt[i] = 0;
  }

  info = dchdc ( a, LDA, N, work, ipvt, job );

  if ( info != N )
  {
    cout << "\n";
    cout << "  DCHDC returned INFO = " << info << "\n";
    cout << "  This means the matrix is not positive definite.\n";
    return;
  }
//
//  Zero out the lower diagonal.
//
  for ( i = 2; i <= N; i++ )
  {
    for ( j = 1; j <= i-1; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }
//
//  Print the factorization.
//
  cout << "\n";
  cout << "  The Cholesky factor U:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }
//
//  Compute the Cholesky product.
//
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      b[i-1+(j-1)*LDA] = 0.0;
      for ( k = 1; k <= N; k++ )
      {
        b[i-1+(j-1)*LDA] = b[i-1+(j-1)*LDA]
          + a[k-1+(i-1)*LDA] * a[k-1+(j-1)*LDA];
      }
    }
  }
  cout << "\n";
  cout << "  The product U' * U:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << b[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests DCHEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA N
# define NZ 1

  double a[LDA*N];
  double b[LDA*N];
  double c[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int l;
  double s[N];
  int seed;
  double work[N];
  double z[N];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a general matrix,\n";
  cout << "  DCHEX can shift columns in a Cholesky factorization.\n";
  cout << "\n";
  cout << "  The number of equations is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= N; i++ )
    {
      if ( i == j-1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else if ( i == j )
      {
        a[i-1+(j-1)*LDA] = 2.0;
      }
      else if ( i == j+1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*LDA] = 0.0;
      }
    }
  }
  for ( i = 1; i <= N; i++ )
  {
    z[i-1] = ( double ) i;
  }

  cout << "\n";
  cout << "  The matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The vector Z:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout << "  " << setw(12) << z[i-1];
  }
//
//  Decompose the matrix.
//
  cout << "\n";
  cout << "  Decompose the matrix.\n";

  job = 0;

  for ( i = 0; i < N; i++ )
  {
    ipvt[i] = 0;
  }

  info = dchdc ( a, LDA, N, work, ipvt, job );

  if ( info != N )
  {
    cout << "\n";
    cout << "  DCHDC returned INFO = " << info << "\n";
    cout << "  This means the matrix is not positive definite.\n";
    return;
  }
//
//  Zero out the lower diagonal.
//
  for ( i = 2; i <= N; i++ )
  {
    for ( j = 1; j <= i-1; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }
//
//  Print the factorization.
//
  cout << "\n";
  cout << "  The Cholesky factor U:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }
//
//  Right circular shift columns L through K.
//
  k = 1;
  l = 3;

  cout << "\n";
  cout << "  Right circular shift columns K  = " << k <<
    " through L = " << l << "\n";

  job = 1;
  dchex ( a, LDA, N, k, l, z, N, NZ, c, s, job );
//
//  Left circular shift columns K+1 through L.
//
  cout << "\n";
  cout << "  Left circular shift columns K+1 = " << k+1 <<
    " through L = " << l << "\n";

  job = 2;
  dchex ( a, LDA, N, k+1, l, z, N, NZ, c, s, job );
//
//  Print the factorization.
//
  cout << "\n";
  cout << "  The shifted Cholesky factor U:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The shifted vector Z:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout << "  " << setw(12) << z[i-1];
  }
//
// Compute the Cholesky product.
//
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      b[i-1+(j-1)*LDA] = 0.0;
      for ( k = 1; k <= N; k++ )
      {
        b[i-1+(j-1)*LDA] = b[i-1+(j-1)*LDA]
          + a[k-1+(i-1)*LDA] * a[k-1+(j-1)*LDA];
      }
    }
  }

  cout << "\n";
  cout << "  The shifted product U' * U:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << b[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  return;
# undef LDA
# undef N
# undef NZ
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests DCHUD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define P 20
# define LDR P
# define NZ 1

  double b[P];
  double beta[P];
  double c[P];
  int i;
  int info;
  int j;
  int job;
  double r[LDR*P];
  double rho[NZ];
  double row[P];
  double s[P];
  int seed;
  double x[P];
  double y[NZ];
  double z[P*NZ];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a general matrix,\n";
  cout << "  DCHUD updates a Cholesky decomposition.\n";
  cout << "\n";
  cout << "  In this example, we use DCHUD to solve a\n";
  cout << "  least squares problem R * b = z.\n";
  cout << "\n";
  cout << "  The number of equations is P = " << P << "\n";
//
//  Initialize.
//
  for ( j = 1; j <= P; j++ )
  {
    for ( i = 1; i <= P; i++ )
    {
      r[i-1+(j-1)*LDR] = 0.0;
    }
  }
  for ( j = 1; j <= NZ; j++ )
  {
    for ( i = 1; i <= P; i++ )
    {
      z[i-1+(j-1)*P] = 0.0;
    }
  }

  for ( i = 1; i <= P; i++ )
  {
    x[i-1] = ( double ) i;
  }
//
//  Use DCHUD to form R, Z and RHO by adding X and Y a row at a time.
//  X is a row of the least squares matrix and Y the right hand side.
//
  seed = 123456789;

  for ( i = 1; i <= P; i++ )
  {
    r8mat_uniform_01 ( 1, P, &seed, row );
    y[0] = 0.0;
    for ( j = 1; j <= P; j++ )
    {
      y[0] = y[0] + row[j-1] * x[j-1];
    }
    rho[0] = 0.0;
    dchud ( r, LDR, P, row, z, P, NZ, y, rho, c, s );
  }
//
//  Generate the least squares solution, b = inverse ( R ) * Z.
//
  for ( j = 1; j <= NZ; j++ )
  {
    for ( i = 1; i <= P; i++ )
    {
      b[i-1] = z[i-1+(j-1)*P];
    }
    job = 01;

    info = dtrsl ( r, LDR, P, b, job );

    cout << "\n";
    cout << "  Solution vector # " << j << "\n";
    cout << "  (Should be (1,2,3...,n))\n";
    cout << "\n";

    for ( i = 1; i <= P; i++ )
    {
      if ( i <= 5 || P-5 < i )
      {
        cout << "  " << setw(6) << i
             << "  " << setw(14) << b[i-1] << "\n";
      }
      if ( i == 5 )
      {
        cout << "  ......  ..............\n";
      }
    }
  }

  return;
# undef LDR
# undef NZ
# undef P
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests DGBCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define ML 1
# define MU 1
# define LDA ( 2 * ML + MU + 1 )

  double a[LDA*N];
  int i;
  int ipivot[N];
  int j;
  int m;
  double rcond;
  double z[N];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For a general banded matrix,\n";
  cout << "  DGBCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  m = ML + MU + 1;
  cout << "  The bandwidth of the matrix is " << m << "\n";

  for ( j = 1; j <= N; j++ )
  {
    a[m-2+(j-1)*LDA] = -1.0;
    a[m-1+(j-1)*LDA] =  2.0;
    a[m  +(j-1)*LDA] = -1.0;
  }
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = dgbco ( a, LDA, N, ML, MU, ipivot, z );

  cout << "\n";
  cout << "  Estimated reciprocal condition = " << rcond << "\n";

  return;
# undef LDA
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests DGBFA and DGBSL.
//
//  Discussion:
//
//    The problem solved here is a larger version of this one:
//
//    Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
//                (-1  2 -1  0  0)                          (0)
//                ( 0 -1  2 -1  0)                          (0)
//                ( 0  0 -1  2 -1)                          (0)
//                ( 0  0  0 -1  2)                          (1)
//
//
//    Solution is   (1)
//                  (1)
//                  (1)
//                  (1)
//                  (1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Local Parameters:
//
//    N is the number of equations.
//
//    ML is the number of subdiagonals,
//    MU the number of superdiagonals.
//
//    LDA is the leading dimension of the array used to store the
//    matrix, which must be at least 2*ML+MU+1.
//
{
# define N 10
# define ML 1
# define MU 1
# define LDA ( 2 * ML + MU + 1 )

  double a[LDA*N];
  double b[N];
  int i;
  int info;
  int ipivot[N];
  int j;
  int job;
  int m;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For a general banded matrix,\n";
  cout << "  DGBFA computes the LU factors,\n";
  cout << "  DGBSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the right hand side B.
//
  b[0] = 1.0;
  for ( i = 2; i <= N-1; i++ )
  {
    b[i-1] = 0.0;
  }
  b[N-1] = 1.0;
//
//  Set the matrix A.
//
  m = ML + MU + 1;
  cout << "  The bandwidth of the matrix is " << m << "\n";

  for ( j = 1; j <= N; j++ )
  {
    a[m-2+(j-1)*LDA] = -1.0;
    a[m-1+(j-1)*LDA] =  2.0;
    a[m  +(j-1)*LDA] = -1.0;
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dgbfa ( a, LDA, N, ML, MU, ipivot );

  if ( info != 0 )
  {
    cout << "  Error!  DGBFA returns INFO = " << info << "\n";
    return;
  }
//
//  Call DGBSL to solve the linear system.  The solution
//  is returned in B, that is, it overwrites the right hand side.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  job = 0;
  dgbsl ( a, LDA, N, ML, MU, ipivot, b, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,1,1,1,1,...,1,1))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }

  return;
# undef LDA
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests DGBFA and DGBDI.
//
//  Discussion:
//
//    Matrix A is ( 2 -1  0  0  0)
//                (-1  2 -1  0  0)
//                ( 0 -1  2 -1  0)
//                ( 0  0 -1  2 -1)
//                ( 0  0  0 -1  2)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Local Parameters:
//
//    N is the number of equations.
//
//    ML is the number of subdiagonals,
//    MU the number of superdiagonals.
//
//    LDA is the leading dimension of the array used to store the
//    matrix, which must be at least 2*ML+MU+1.
//
{
# define N_MAX 128
# define ML 1
# define MU 1
# define LDA ( 2 * ML + MU + 1 )

  double a[LDA*N_MAX];
  double det[2];
  int i;
  int ihi;
  int ilo;
  int info;
  int ipivot[N_MAX];
  int j;
  int m;
  int n;
  int n_log;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a general banded matrix,\n";
  cout << "  DGBFA factors the matrix,\n";
  cout << "  DGBDI computes the determinant as\n";
  cout << "    det = MANTISSA * 10^EXPONENT\n";
  cout << "\n";
  cout << "  Find the determinant of the -1,2,-1 matrix\n";
  cout << "  for N = 2, 4, 8, 16, 32, 64, 128.\n";
  cout << "\n";
  cout << "  (For this matrix, det ( A ) = N + 1.)\n";
//
//  Set the matrix A.
//
  m = ML + MU + 1;
  cout << "  The bandwidth of the matrix is " << m << "\n";
  cout << "\n";
  cout << "       N    Mantissa       Exponent\n";
  cout << "\n";

  n = 1;

  for ( n_log = 1; n_log <= 7; n_log++ )
  {
    n = 2 * n;

    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= LDA; i++ )
      {
        a[i-1+(j-1)*LDA] = 0.0;
      }
    }

    for ( j = 1; j <= n; j++ )
    {
      i = j;
      a[i-j+ML+MU+(j-1)*LDA] = 2.0;
    }

    for ( j = 2; j <= n; j++ )
    {
      i = j - 1;
      a[i-j+ML+MU+(j-1)*LDA] = -1.0;
    }

    for ( j = 1; j <= n-1; j++ )
    {
      i = j + 1;
      a[i-j+ML+MU+(j-1)*LDA] = -1.0;
    }

    info = dgbfa ( a, LDA, n, ML, MU, ipivot );

    if ( info != 0 )
    {
      cout << "  Error!  DGBFA returns INFO = " << info << "\n";
      return;
    }

    dgbdi ( a, LDA, n, ML, MU, ipivot, det );

    cout << "  " << setw(6)  << n
         << "  " << setw(14) << det[0]
         << "  " << setw(14) << det[1] << "\n";
  }

  return;
# undef LDA
# undef ML
# undef MU
# undef N_MAX
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests DGBFA and DGBSL.
//
//  Discussion:
//
//    DGBFA and DGBSL are for general banded matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 100
# define ML 25
# define MU 25
# define LDA ( 2 * ML + MU + 1 )

  double a[LDA*N];
  double b[N];
  int i;
  int ihi;
  int ilo;
  int info;
  int ipivot[N];
  int j;
  int job;
  int m;
  double temp;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  For a general banded matrix,\n";
  cout << "  DGBFA computes the LU factors,\n";
  cout << "  DGBSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Assign values to matrix A and right hand side B.
//
//  We want to try a problem with a significant bandwidth.
//
  m = ML + MU + 1;
  cout << "  The bandwidth of the matrix is " << m << "\n";

  for ( j = 1; j <= N; j++ )
  {
    ilo = i4_max ( 1, j - MU );
    ihi = i4_min ( N, j + ML );

    temp = 0.0;
    for ( i = ilo; i <= ihi; i++ )
    {
      a[i-j+m-1+(j-1)*LDA] = -1.0;
      temp = temp - 1.0;
    }
    temp = temp + 1.0;
    a[m-1+(j-1)*LDA] = 4.0 - temp;
    b[j-1] = 4.0;
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dgbfa ( a, LDA, N, ML, MU, ipivot );

  if ( info != 0 )
  {
    cout << "  Error!  DGBFA returns INFO = " << info << "\n";
    return;
  }
//
//  Call DGBSL to solve the linear system.  The solution
//  is returned in B, that is, it overwrites the right hand side.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  job = 0;
  dgbsl ( a, LDA, N, ML, MU, ipivot, b, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,1,1,1,1,...,1,1))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }

  return;
# undef LDA
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 calls DGECO and DGESL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Local Parameters:
//
//    LDA defines the maximum matrix size we will use.
//
{
# define LDA 10

  double a[LDA*LDA];
  double b[LDA];
  int i;
  int ipvt[LDA];
  int job;
  int n;
  double rcond;
  double z[LDA];

  n = 3;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For a general matrix,\n";
  cout << "  DGECO computes the LU factors and computes\n";
  cout << "  its reciprocal condition number;\n";
  cout << "  DGESL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << n << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*LDA] = 1.0;
  a[0+1*LDA] = 2.0;
  a[0+2*LDA] = 3.0;

  a[1+0*LDA] = 4.0;
  a[1+1*LDA] = 5.0;
  a[1+2*LDA] = 6.0;

  a[2+0*LDA] = 7.0;
  a[2+1*LDA] = 8.0;
  a[2+2*LDA] = 0.0;
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  rcond = dgeco ( a, LDA, n, ipvt, z );

  cout << "  The reciprocal matrix condition number = " << rcond << "\n";

  if ( rcond + 1.0 == 1.0 )
  {
    cout << "  Error!  The matrix is nearly singular!\n";
    return;
  }
//
//  Set a right hand side.
//
  b[0] = 6.0;
  b[1] = 15.0;
  b[2] = 15.0;
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  job = 0;
  dgesl ( a, LDA, n, ipvt, b, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Solution returned by DGESL\n";
  cout << "  (Should be (1,1,1))\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(14) << b[i-1] << "\n";
  }
//
//  A second right hand side can be solved without refactoring a.
//
  cout << "\n";
  cout << "  Call DGESL for a new right hand\n";
  cout << "  side for the same, factored matrix.\n";
//
//  Set the right hand side.
//
  b[0] = 1.0;
  b[1] = 4.0;
  b[2] = 7.0;
//
//  Solve the system.
//
  cout << "\n";
  cout << "  Solve a linear system.\n";

  job = 0;
  dgesl ( a, LDA, n, ipvt, b, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Solution returned by DGESL\n";
  cout << "  (should be (1,0,0))\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(14) << b[i-1] << "\n";
  }
//
//  The transposed problem A'*X = B can be solved by DGESL
//  as well, without any refactoring.
//
  cout << "\n";
  cout << "  Call DGESL for transposed problem.\n";
//
//  Set the right hand side.
//
  b[0] = 6.0;
  b[1] = 6.0;
  b[2] = -3.0;
//
//  Solve the transposed system.
//
  cout << "\n";
  cout << "  Call DGESL to solve a transposed linear system.\n";

  job = 1;
  dgesl ( a, LDA, n, ipvt, b, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Solution returned by DGESL\n";
  cout << "  (should be (-1,0,1))\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(14) << b[i-1] << "\n";
  }

  return;
# undef LDA
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests DGEFA and DGEDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define LDA 3
//
//  Matrix listed by columns.
//
  double a[LDA*N] = {
    1.0, 4.0, 7.0,
    2.0, 5.0, 8.0,
    3.0, 6.0, 0.0 };
  double det[2];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  double work[N];

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For a general banded matrix,\n";
  cout << "  DGEFA computes the LU factors;\n";
  cout << "  DGEDI computes the inverse and determinant.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dgefa ( a, LDA, N, ipvt );

  if ( info != 0 )
  {
    cout << "  Error!  The matrix is nearly singular!\n";
    return;
  }
//
//  Get the inverse and determinant.
//
  cout << "\n";
  cout << "  Get the inverse and determinant.\n";

  job = 11;
  dgedi ( a, LDA, N, ipvt, det, work, job );

  cout << "\n";
  cout << "  The determinant = " << det[0] << " * 10^"<< det[1] << "\n";

  cout << "\n";
  cout << "  The inverse matrix:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests DGEFA and DGESL.
//
//  Discussion:
//
//    Solve A*x = b where A is a given matrix, and B a right hand side.
//
//    We will also assume that A is stored in the simplest
//    possible way.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define LDA N
//
//  The entries of the matrix A are listed by columns.
//
  double a[LDA*N] = {
    1.0, 4.0, 7.0,
    2.0, 5.0, 8.0,
    3.0, 6.0, 0.0 };
  double b[N] = { 6.0, 15.0, 15.0 };
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For a general banded matrix,\n";
  cout << "  DGEFA computes the LU factors;\n";
  cout << "  DGESL solves a factored linear system;\n";
  cout << "\n";
  cout << "  The number of equations is N = " << N << "\n";

  cout << "\n";
  cout << "  The matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(14) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The right hand side B:\n";
  cout << "\n";
  for ( i = 1; i <= N; i++ )
  {
    cout << "  " << setw(14) << b[i-1];
  }
  cout << "\n";
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dgefa ( a, LDA, N, ipvt );

  if ( info != 0 )
  {
    cout << "  DGEFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  job = 0;
  dgesl ( a, LDA, N, ipvt, b, job );

  cout << "\n";
  cout << "  DGESL returns the solution:\n";
  cout << "  (Should be (1,1,1))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout << "  " << setw(14) << b[i-1];
  }
  cout << "\n";

  return;
# undef N
# undef LDA
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests DGEFA and DGESL.
//
//  Discussion:
//
//    In this example, we solve a relatively large linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 100
# define LDA N

  double a[LDA*N];
  double b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  For a general banded matrix,\n";
  cout << "  DGEFA computes the LU factors;\n";
  cout << "  DGESL solves a factored linear system;\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Assign values to the matrix A and the right hand side B.
//
//  The problem is just an enlarged version of the
//  problem for N = 5, which is:
//
//  Matrix A is ( n -1 -1 -1 -1)    Right hand side B is  (1)
//              (-1  n -1 -1 -1)                          (1)
//              (-1 -1  n -1 -1)                          (1)
//              (-1 -1 -1  n -1)                          (1)
//              (-1 -1 -1 -1  n)                          (1)
//
//  Solution is   (1)
//                (1)
//                (1)
//                (1)
//                (1)
//
  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 1.0;
  }

  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= N; i++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*LDA] = ( double ) N;
      }
      else
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
    }
   }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dgefa ( a, LDA, N, ipvt );

  if ( info != 0 )
  {
    cout << "  DGEFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  cout << "\n";
  cout << "  Solve the factored system.\n";

  job = 0;
  dgesl ( a, LDA, N, ipvt, b, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,1,1,1,1,...,1,1))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests DGTSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 100

  double b[N];
  double c[N];
  double d[N];
  double e[N];
  int i;
  int info;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  For a general tridiagonal matrix,\n";
  cout << "  DGTSL factors and solves a linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
  cout << "\n";
//
//  Set up the linear system, by storing the values of the
//  subdiagonal, diagonal, and superdiagonal in C, D, and E,
//  and the right hand side in B.
//
  c[0] = 0.0;
  for ( i = 2; i <= N; i++ )
  {
    c[i-1] = -1.0;
  }

  for ( i = 1; i <= N; i++ )
  {
    d[i-1] = 2.0;
  }

  for ( i = 1; i <= N-1; i++ )
  {
    e[i-1] = -1.0;
  }
  e[N-1] = 0.0;

  for ( i = 1; i <= N-1; i++ )
  {
    b[i-1] = 0.0;
  }
  b[N-1] = ( double ) ( N + 1 );
//
// Factor the matrix and solve the system.
//
  cout << "\n";
  cout << "  Factor the matrix and solve the system.\n";

  info = dgtsl ( N, c, d, e, b );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  DGTSL returns nonzero INFO = " << info << "\n";
    return;
  }
//
//  Print the results.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,2,3,4,5,...,n-1,n))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }

  return;
# undef N
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests DPBCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define LDA 2

  double a[LDA*N];
  int i;
  int info;
  int j;
  int m;
  double rcond;
  double z[N];

  cout << "\n";
  cout << "TEST13\n";
  cout << "  For a positive definite symmetric banded matrix,\n";
  cout << "  DPBCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the number of nonzero diagonals.
//
  m = 1;
//
//  Set the value of the subdiagonal and diagonal.
//
  for ( j = 1; j <= N; j++ )
  {
    a[0+(j-1)*LDA] = -1.0;
    a[1+(j-1)*LDA] = 2.0;
  }

  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = dpbco ( a, LDA, N, m, z );

  cout << "\n";
  cout << "  Reciprocal condition  = " << rcond << "\n";

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests DPBDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 128
# define LDA 2

  double a[LDA*N_MAX];
  double det[2];
  int i;
  int info;
  int j;
  int m;
  int n;
  int n_log;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  For a positive definite symmetric banded matrix,\n";
  cout << "  DPBDI computes the determinant as\n";
  cout << "    det = MANTISSA * 10^EXPONENT\n";
  cout << "\n";
  cout << "  Find the determinant of the -1,2,-1 matrix\n";
  cout << "  for N = 2, 4, 8, 16, 32, 64, 128.\n";
  cout << "\n";
  cout << "  (For this matrix, det ( A ) = N + 1.)\n";
//
//  Set the number of  nonzero diagonals.
//
  m = 1;

  cout << "\n";
  cout << "       N    Mantissa       Exponent\n";
  cout << "\n";

  n = 1;

  for ( n_log = 1; n_log <= 7; n_log++ )
  {
    n = 2 * n;

    a[0+0*LDA] =  0.0;
    for ( j = 2; j <= n; j++ )
    {
      a[0+(j-1)*LDA] = -1.0;
    }
    for ( j = 1; j <= n; j++ )
    {
      a[1+(j-1)*LDA] = 2.0;
    }

    info = dpbfa ( a, LDA, n, m );

    if ( info != 0 )
    {
      cout << "  Error!  DPBFA returns INFO = " << info << "\n";
      return;
    }

    dpbdi ( a, LDA, n, m, det );

    cout << "  " << setw(6)  << n
         << "  " << setw(14) << det[0]
         << "  " << setw(14) << det[1] << "\n";
  }

  return;
# undef LDA
# undef N_MAX
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests DPBFA and DPBSL.
//
//  Discussion:
//
//    DPBFA and DPBSL are for a positive definite symmetric band matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define LDA 2

  double a[LDA*N];
  double b[N];
  int i;
  int info;
  int j;
  int m;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  For a positive definite symmetric banded matrix,\n";
  cout << "  DPBFA computes the LU factors.\n";
  cout << "  DPBSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Assign values to matrix A and right hand side B.
//
//  The problem is just an enlarged version of the
//  problem for N = 5, which is:
//
//  Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
//              (-1  2 -1  0  0)                          (0)
//              ( 0 -1  2 -1  0)                          (0)
//              ( 0  0 -1  2 -1)                          (0)
//              ( 0  0  0 -1  2)                          (1)
//
//
//  solution is   (1)
//                (1)
//                (1)
//                (1)
//                (1)
//
//  Set the right hand side.
//
  b[0] = 1.0;
  for ( i = 2; i <= N-1; i++ )
  {
    b[i-1] = 0.0;
  }
  b[N-1] = 1.0;
//
//  Set the number of nonzero diagonals.
//
  m = 1;
//
//  Set the value of the subdiagonal and diagonal.
//
  for ( j = 1; j <= N; j++ )
  {
    a[0+(j-1)*LDA] = -1.0;
    a[1+(j-1)*LDA] = 2.0;
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dpbfa ( a, LDA, N, m );

  if ( info != 0 )
  {
    cout << "  Error!  DPBFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  dpbsl ( a, LDA, N, m, b );
//
//  Print the results.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,1,1,1,1,...,1,1))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }
  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests DPOCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA N

  double a[LDA*N];
  int i;
  int info;
  int j;
  int job;
  double rcond;
  double z[N];

  cout << "\n";
  cout << "TEST16\n";
  cout << "  For a positive definite symmetric banded matrix,\n";
  cout << "  DPOCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a[i+j*LDA] = 0.0;
    }
  }

  for ( i = 1; i <= N; i++ )
  {
    a[i-1+(i-1)*LDA] = 2.0;
    if ( 1 < i )
    {
      a[i-1+(i-2)*LDA] = -1.0;
    }
    if ( i < N )
    {
      a[i-1+(i)*LDA] = -1.0;
    }
  }

  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = dpoco ( a, LDA, N, z );

  cout << "\n";
  cout << "  Reciprocal condition  = " << rcond << "\n";

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests DPOFA and DPODI.
//
//  Discussion:
//
//    DPOFA factors a positive definite symmetric matrix,
//    and DPODI can compute the determinant or the inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA N

  double a[LDA*N];
  double det[2];
  int i;
  int info;
  int j;
  int job;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  For a positive definite symmetric matrix,\n";
  cout << "  DPOFA computes the LU factors.\n";
  cout << "  DPODI computes the inverse or determinant.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a[i+j*LDA] = 0.0;
    }
  }

  for ( i = 1; i <= N; i++ )
  {
    a[i-1+(i-1)*LDA] = 2.0;
    if ( 1 < i )
    {
      a[i-1+(i-2)*LDA] = -1.0;
    }
    if ( i < N )
    {
      a[i-1+(i)*LDA] = -1.0;
    }
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dpofa ( a, LDA, N );

  if ( info != 0 )
  {
    cout << "  Error, DPOFA returns INFO = " << info << "\n";
    return;
  }
//
//  Invert the matrix.
//
  cout << "\n";
  cout << "  Get the determinant and inverse.\n";

  job = 11;
  dpodi ( a, LDA, N, det, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Determinant = " << det[0] << " * 10^" << det[1] << "\n";
//
//  DPODI produces only the 'upper half triangle' of the inverse,
//  which is actually symmetric.  Thus, the lower half could be
//  produced by copying from the upper half.  However, the first row
//  of A, as returned, is exactly the first row of the inverse.
//
  cout << "\n";
  cout << "  First row of inverse:\n";
  cout << "\n";
  for ( j = 1; j <= N; j++ )
  {
    cout << "  " << setw(12) << a[0+(j-1)*LDA];
  }
  cout << "\n";

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests DPOFA and DPOSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 20
# define LDA N

  double a[LDA*N];
  double b[N];
  int i;
  int info;
  int j;
  int job;
  double x[N];

  cout << "\n";
  cout << "TEST18\n";
  cout << "  For a positive definite symmetric matrix,\n";
  cout << "  DPOFA computes the LU factors.\n";
  cout << "  DPOSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a[i+j*LDA] = 0.0;
    }
  }

  for ( i = 1; i <= N; i++ )
  {
    a[i-1+(i-1)*LDA] = 2.0;
    if ( 1 < i )
    {
      a[i-1+(i-2)*LDA] = -1.0;
    }
    if ( i < N )
    {
      a[i-1+(i)*LDA] = -1.0;
    }
  }
//
//  Set the right hand side.
//
  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( double ) i;
  }
  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 0.0;
    for ( j = 1; j <= N; j++ )
    {
      b[i-1] = b[i-1] + a[i-1+(j-1)*LDA] * x[j-1];
    }
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dpofa ( a, LDA, N );

  if ( info != 0 )
  {
    cout << "  Error, DPOFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  dposl ( a, LDA, N, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,2,3,4,5,...,n-1,n))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests DPPCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[(N*(N+1))/2];
  int i;
  int j;
  int k;
  double rcond;
  double z[N];

  cout << "\n";
  cout << "TEST19\n";
  cout << "  For a positive definite symmetric packed matrix,\n";
  cout << "  DPPCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  k = 0;
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      k = k + 1;
      if ( i == j - 1 )
      {
        a[k-1] = -1.0;
      }
      else if ( i == j )
      {
        a[k-1] = 2.0;
      }
      else
      {
        a[k-1] = 0.0;
      }
    }
  }
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition number.\n";

  rcond = dppco ( a, N, z );

  cout << "\n";
  cout << "  Reciprocal condition number = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests DPPFA and DPPDI.
//
//  Discussion:
//
//    DPPFA factors a packed positive definite symmetric matrix,
//    and DPPDI can compute the determinant or the inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[(N*(N+1))/2];
  double b[N*N];
  double det[2];
  int i;
  int info;
  int j;
  int job;
  int k;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  For a positive definite symmetric packed matrix,\n";
  cout << "  DPPFA computes the LU factors.\n";
  cout << "  DPPDI computes the inverse or determinant.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  k = 0;
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      k = k + 1;
      if ( i == j - 1 )
      {
        a[k-1] = -1.0;
      }
      else if ( i == j )
      {
        a[k-1] = 2.0;
      }
      else
      {
        a[k-1] = 0.0;
      }
    }
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dppfa ( a, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, DPPFA returns INFO = " << info << "\n";
    return;
  }
//
//  Invert the matrix.
//
  cout << "\n";
  cout << "  Get the determinant and inverse.\n";

  job = 11;
  dppdi ( a, N, det, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Determinant = " << det[0] << " * 10^" << det[1] << "\n";
//
//  DPPDI produces only the 'upper half triangle' of the inverse,
//  which is actually symmetric.  Thus, the lower half could be
//  produced by copying from the upper half.  However, the first row
//  of A, as returned, is exactly the first row of the inverse.
//
  k = 0;
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      k = k + 1;
      b[i-1+(j-1)*N] = a[k-1];
      b[j-1+(i-1)*N] = a[k-1];
    }
  }

  cout << "\n";
  cout << "  The inverse matrix:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << b[i-1+(j-1)*N];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests DPPFA and DPPSL.
//
//  Discussion:
//
//    DPOFA factors a positive definite symmetric matrix,
//    and DPOSL can solve a factored linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double a[(N*(N+1))/2];
  double b[N];
  int i;
  int info;
  int j;
  int job;
  int k;
  double x[N];

  cout << "\n";
  cout << "TEST21\n";
  cout << "  For a positive definite symmetric packed matrix,\n";
  cout << "  DPPFA computes the LU factors.\n";
  cout << "  DPPSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( double ) i;
  }

  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 0.0;
  }
//
//  Set the matrix A.
//
  k = 0;
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      k = k + 1;
      if ( i == j - 1 )
      {
        a[k-1] = -1.0;
        b[i-1] = b[i-1] + a[k-1] * x[j-1];
        b[j-1] = b[j-1] + a[k-1] * x[i-1];
      }
      else if ( i == j )
      {
        a[k-1] = 2.0;
        b[i-1] = b[i-1] + a[k-1] * x[i-1];
      }
      else
      {
        a[k-1] = 0.0;
      }
    }
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dppfa ( a, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, DPPFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  dppsl ( a, N, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,2,3,4,5,...,n-1,n))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }
  return;
# undef N
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests DPTSL.
//
//  Discussion:
//
//    DPTSL factors and solves a positive definite symmetric tridiagonal system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double b[N];
  double d[N];
  double e[N];
  int i;
  double x[N];

  cout << "\n";
  cout << "TEST22\n";
  cout << "  For a positive definite symmetric tridiagonal matrix,\n";
  cout << "  DPTSL factors and solves a linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the matrix A.
//
  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( double ) i;
  }

  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 0.0;
  }
  for ( i = 1; i <= N; i++ )
  {
    d[i-1] = 2.0;
  }
  for ( i = 1; i <= N-1; i++ )
  {
    e[i-1] = -1.0;
  }
  e[N-1] = 0.0;

  for ( i = 1; i <= N; i++ )
  {
    if ( 1 < i )
    {
      b[i-1] = b[i-1] + e[i-2] * x[i-2];
    }
    b[i-1] = b[i-1] + d[i-1] * x[i-1];
    if ( i < N )
    {
      b[i-1] = b[i-1] + e[i-1] * x[i];
    }
  }
//
//  Factor and solve the system.
//
  cout << "\n";
  cout << "  Factor the matrix and solve the system.\n";

  dptsl ( N, d, e, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,2,3,4,5,...,n-1,n))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }

  return;
# undef N
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests DQRDC and DQRSL.
//
//  Discussion:
//
//    DQRDC and DQRSL compute the QR factorization, and use it
//    to solve linear systems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define P 3
# define LDA N

  double a[LDA*P] = {
    1.0, 1.0, 0.0,
    1.0, 0.0, 1.0,
    0.0, 1.0, 1.0 };
  double b[LDA*P];
  int i;
  int info;
  int ipvt[P];
  int j;
  int job;
  int k;
  double q[N*N];
  double qraux[P];
  double qty[N];
  double qy[N];
  double r[N*P];
  double rsd[N];
  double work[P];
  double xb[N];
  double y[N];

  cout << "\n";
  cout << "TEST23\n";
  cout << "  For a general rectangular matrix,\n";
  cout << "  DQRDC computes the QR decomposition of a\n";
  cout << "  matrix, but does not return Q and R explicitly.\n";
  cout << "\n";
  cout << "  Show how Q and R can be recovered using DQRSL.\n";
//
//  Print the matrix A.
//
  cout << "\n";
  cout << "  The matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= P; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }
//
//  Decompose the matrix.
//
  cout << "\n";
  cout << "  Decompose the matrix.\n";

  job = 0;
  for ( j = 1; j <= P; j++ )
  {
    ipvt[j] = 0;
  }

  dqrdc ( a, LDA, N, P, qraux, ipvt, work, job );
//
//  Print out what DQRDC has stored in A...
//
  cout << "\n";
  cout << "  The packed matrix A which describes Q and R:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= P; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }
//
//  ...and in QRAUX.
//
  cout << "\n";
  cout << "  The QRAUX vector, containing some additional\n";
  cout << "  information defining Q:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout << "  " << setw(12) << qraux[i-1];
  }
  cout << "\n";
//
//  Print out the resulting R factor.
//
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= P; j++ )
    {
      if ( j < i )
      {
        r[i-1+(j-1)*N] = 0.0;
      }
      else
      {
        r[i-1+(j-1)*N] = a[i-1+(j-1)*LDA];
      }
    }
  }

  cout << "\n";
  cout << "  The R factor:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= P; j++ )
    {
      cout << "  " << setw(12) << r[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }
//
//  Call DQRSL to extract the information about the Q matrix.
//  We do this, essentially, by asking DQRSL to tell us the
//  value of Q*Y, where Y is a column of the identity matrix.
//
  job = 10000;

  for ( i = 1; i <= N; i++ )
  {
//
//  Set the vector Y.
//
    for ( j = 1; j <= N; j++ )
    {
      y[j-1] = 0.0;
    }
    y[i-1] = 1.0;
//
//  Ask DQRSL to tell us what Q*Y is.
//
    info = dqrsl ( a, LDA, N, P, qraux, y, qy, qty, b, rsd, xb, job );

    if ( info != 0 )
    {
      cout << "  Error!  DQRSL returns INFO = " << info << "\n";
      return;
    }
//
//  Copy QY into the appropriate column of Q.
//
    for ( j = 1; j <= N; j++ )
    {
      q[j-1+(i-1)*N] = qy[j-1];
    }
  }
//
//  Now print out the Q matrix we have extracted.
//
  cout << "\n";
  cout << "  The Q factor:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << q[i-1+(j-1)*N];
    }
    cout << "\n";
  }
//
//  Compute Q*R to verify that it equals A.
//
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= P; j++ )
    {
      b[i-1+(j-1)*LDA] = 0.0;
      for ( k = 1; k <= N; k++ )
      {
        b[i-1+(j-1)*LDA] = b[i-1+(j-1)*LDA]
          + q[i-1+(k-1)*N] * r[k-1+(j-1)*N];
      }
    }
  }
//
//  Print the result.
//
  cout << "\n";
  cout << "  The product Q * R:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= P; j++ )
    {
      cout << "  " << setw(12) << b[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  return;
# undef LDA
# undef N
# undef P
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests DSICO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 100
# define LDA N

  double a[LDA*N];
  int i;
  int ipvt[N];
  int j;
  double rcond;
  double z[N];

  cout << "\n";
  cout << "TEST24\n";
  cout << "  For a symmetric indefinite matrix,\n";
  cout << "  DSICO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Assign values to the matrix A and the right hand side B.
//
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*LDA] = 2.0;
      }
      else if ( j == i+1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*LDA] = 0.0;
      }
    }
  }
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = dsico ( a, LDA, N, ipvt, z );

  cout << "\n";
  cout << "  Estimated reciprocal condition = " << rcond << "\n";

  return;
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests DSIFA and DSISL.
//
//  Discussion:
//
//    DSIFA and DSISL are for symmetric indefinite matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 100
# define LDA N

  double a[LDA*N];
  double b[N];
  int i;
  int info;
  int ipvt[N];
  int j;

  cout << "\n";
  cout << "TEST25\n";
  cout << "  For a symmetric indefinite matrix,\n";
  cout << "  DSIFA factor a symmetric indefinite matrix;\n";
  cout << "  DSISL solves a factored linear system,\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Assign values to the matrix A and the right hand side B.
//
  for ( i = 1; i < N; i++ )
  {
    b[i-1] = 0.0;
  }
  b[N-1] = ( double ) ( N + 1 );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*LDA] = 2.0;
      }
      else if ( j == i+1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*LDA] = 0.0;
      }
    }
  }
//
//  Factor the matrix A.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = dsifa ( a, LDA, N, ipvt );

  if ( info != 0 )
  {
    cout << "  Error!  DSIFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  dsisl ( a, LDA, N, ipvt, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,2,3,4,5,...,n-1,n))\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests DSPCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 100

  double a[(N*(N+1))/2];
  int i;
  int ipvt[N];
  int j;
  int k;
  double rcond;
  double z[N];

  cout << "\n";
  cout << "TEST26\n";
  cout << "  For a symmetric indefinite packed matrix,\n";
  cout << "  DSPCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Assign values to the matrix A.
//
  k = 0;
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      k = k + 1;
      if ( i == j )
      {
        a[k-1] = 2.0;
      }
      else if ( j == i+1 )
      {
        a[k-1] = -1.0;
      }
      else
      {
        a[k-1] = 0.0;
      }
    }
  }
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = dspco ( a, N, ipvt, z );

  cout << "\n";
  cout << "  Estimated reciprocal condition = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests DSPFA and DSPSL.
//
//  Discussion:
//
//    DSPFA and DSPSL are for packed symmetric indefinite matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 100

  double a[(N*(N+1))/2];
  double b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  For a symmetric indefinite packed matrix,\n";
  cout << "  DSPFA computes the LU factors,\n";
  cout << "  DSPSL solves a factored linear system,\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Assign values to the matrix A and the right hand side B.
//
  for ( i = 1; i <= N-1; i++ )
  {
    b[i-1] = 0.0;
  }
  b[N-1] = ( double ) ( N + 1 );

  k = 0;
  for ( j = 1; j <= N; j++ )
  {
    for ( i = 1; i <= j; i++ )
    {
      k = k + 1;
      if ( i == j )
      {
        a[k-1] = 2.0;
      }
      else if ( j == i+1 )
      {
        a[k-1] = -1.0;
      }
      else
      {
        a[k-1] = 0.0;
      }
    }
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n" << flush;

  info = dspfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "  Error!  DSPFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n" << flush;

  dspsl ( a, N, ipvt, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The first and last 5 entries of solution:\n";
  cout << "  (Should be (1,2,3,4,5,...,n-1,n))\n" << flush;
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    if ( i <= 5 || N-5 < i )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(14) << b[i-1] << "\n";
    }
    if ( i == 5 )
    {
      cout << "  ......  ..............\n";
    }
  }


  return;
# undef N
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests DSVDC.
//
//  Discussion:
//
//    DSVDC computes the singular value decomposition:
//
//      A = U * S * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 4

  double a[M*N];
  double b[M*N];
//
//  E must be dimensioned at least maximum(M+1,N).
//
  double e[M+N];
  int i;
  int info;
  int j;
  int job;
  int k;
  int lda;
  int ldu;
  int ldv;
//
//  S must be dimensioned at least maximum(M+1,N).
//
  double s[M+N];
  int seed;
  double sigma[M*N];
  double u[M*M];
  double v[N*N];
  double work[M];

  cout << "\n";
  cout << "TEST28\n";
  cout << "  For an MxN matrix A in general storage,\n";
  cout << "  DSVDC computes the singular value decomposition:\n";
  cout << "    A = U * S * V'\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
//
//  Set A.
//
  seed = 123456789;

  r8mat_uniform_01 ( M, N, &seed, a );

  cout << "\n";
  cout << "  The matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(10) << a[(i-1)+(j-1)*M];
    }
    cout << "\n";
  }
//
//  Decompose the matrix.
//
  cout << "\n";
  cout << "  Decompose the matrix.\n";

  job = 11;
  lda = M;
  ldu = M;
  ldv = N;

  info = dsvdc ( a, lda, M, N, s, e, u, ldu, v, ldv, work, job );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning:\n";
    cout << "  DSVDC returned nonzero INFO = " << info << "\n";
    return;
  }

  cout << "\n";
  cout << "  Singular values:\n";
  cout << "\n";

  for ( i = 1; i <= i4_min ( M, N ); i++ )
  {
    cout << "  "
         << setw(4)  << i+1  << "  "
         << setw(14) << s[i-1] << "\n";
  }

  cout << "\n";
  cout << "  Left Singular Vector Matrix U:\n";
  cout << "\n";

  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= M; j++ )
    {
      cout << "  " << setw(10) << u[(i-1)+(j-1)*M];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Right Singular Vector Matrix V:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(10) << v[(i-1)+(j-1)*N];
    }
    cout << "\n";
  }

  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      if ( i == j )
      {
        sigma[(i-1)+(j-1)*M] = s[i-1];
      }
      else
      {
        sigma[(i-1)+(j-1)*M] = 0.0;
      }
    }
  }
//
//  Verify that A = U * S * V'.
//
  for ( i = 1; i <= M; i++ )
  {
    for ( k = 1; k <= N; k++ )
    {
      b[(i-1)+(k-1)*M] = 0.0;
      for ( j = 1; j <= N; j++ )
      {
      b[(i-1)+(k-1)*M] = b[(i-1)+(k-1)*M] + sigma[i-1+(j-1)*M] * v[k-1+(j-1)*N];
      }
    }
  }

  for ( i = 1; i <= M; i++ )
  {
    for ( k = 1; k <= N; k++ )
    {
      a[(i-1)+(k-1)*M] = 0.0;
      for ( j = 1; j <= M; j++ )
      {
      a[(i-1)+(k-1)*M] = a[(i-1)+(k-1)*M] + u[i-1+(j-1)*M] * b[j-1+(k-1)*M];
      }
    }
  }

  cout << "\n";
  cout << "  The product U * S * V' (should equal A):\n";
  cout << "\n";

  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(10) << a[(i-1)+(j-1)*M];
    }
    cout << "\n";
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests DTRCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA N

  double a[LDA*N];
  int i;
  int j;
  int job;
  double rcond;
  int  seed = 123456789;
  double z[N];

  cout << "\n";
  cout << "TEST29\n";
  cout << "  For a triangular matrix,\n";
  cout << "  DTRCO computes the LU factors and\n";
  cout << "  computes its reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Lower triangular matrix A.
//
  r8mat_uniform_01 ( LDA, N, &seed, a );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = i+1; j <= N; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }

  cout << "\n";
  cout << "  Lower triangular matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  job = 0;
  rcond = dtrco ( a, LDA, N, z, job );

  cout << "\n";
  cout << "  The reciprocal condition number = " << rcond << "\n";
//
//  Upper triangular matrix A.
//
  r8mat_uniform_01 ( LDA, N, &seed, a );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= i-1; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }

  cout << "\n";
  cout << "  Upper triangular matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  job = 1;

  rcond = dtrco ( a, LDA, N, z, job );

  cout << "\n";
  cout << "  The reciprocal condition number = " << rcond << "\n";

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests DTRDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA N

  double a[LDA*N];
  double det[2];
  int i;
  int info;
  int j;
  int job;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST30\n";
  cout << "  For a triangular matrix,\n";
  cout << "  DTRDI computes the determinant or inverse.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Lower triangular matrix A.
//
  r8mat_uniform_01 ( N, N, &seed, a );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = i+1; j <= N; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }

  cout << "\n";
  cout << "  Lower triangular matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  job = 110;

  info = dtrdi ( a, LDA, N, det, job );

  cout << "\n";
  cout << "  The determinant = " << det[0] << " * 10^(" << det[1] << ").\n";

  cout << "\n";
  cout << "  The inverse matrix:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }
//
//  Upper triangular matrix A.
//
  r8mat_uniform_01 ( N, N, &seed, a );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= i-1; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }

  cout << "\n";
  cout << "  Upper triangular matrix A:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  job = 111;

  info = dtrdi ( a, LDA, N, det, job );

  cout << "\n";
  cout << "  The determinant = " << det[0] << " * 10^(" << det[1] << ").\n";

  cout << "\n";
  cout << "  The inverse matrix:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
    }
    cout << "\n";
  }

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests DTRSL.
//
//  Discussion:
//
//    DTRSL solves triangular linear systems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA 5

  double a[LDA*N];
  double b[N];
  int i;
  int info;
  int j;
  int job;
  int seed = 123456789;
  double x[N];

  cout << "\n";
  cout << "TEST31\n";
  cout << "  For a triangular matrix,\n";
  cout << "  DTRSL solves a linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Lower triangular matrix A.
//
  r8mat_uniform_01 ( N, N, &seed, a );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = i+1; j <= N; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }

  for ( i = 1; i <= N; i++)
  {
    x[i-1] = ( double ) ( i );
  }

  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 0.0;
    for ( j = 1; j <= N; j++ )
    {
      b[i-1] = b[i-1] + a[i-1+(j-1)*LDA] * x[j-1];
    }
  }

  cout << "\n";
  cout << "  For a lower triangular matrix A,\n";
  cout << "  solve A * x = b\n";

  job = 00;

  info = dtrsl ( a, LDA, N, b, job );

  cout << "\n";
  cout << "  The solution (should be 1,2,3,4,5):\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout                       << "  "
         << setw(6)  << i      << "  "
         << setw(14) << b[i-1] << "\n";
  }

  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 0.0;
    for ( j = 1; j <= N; j++ )
    {
      b[i-1] = b[i-1] + a[j-1+(i-1)*LDA] * x[j-1];
    }
  }

  cout << "\n";
  cout << "  For a lower triangular matrix A,\n";
  cout << "  solve A' * x = b\n";

  job = 10;

  info = dtrsl ( a, LDA, N, b, job );

  cout << "\n";
  cout << "  The solution (should be 1,2,3,4,5):\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout                       << "  "
         << setw(6)  << i      << "  "
         << setw(14) << b[i-1] << "\n";
  }
//
//  Upper triangular matrix A.
//
  r8mat_uniform_01 ( N, N, &seed, a );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= i-1; j++ )
    {
      a[i-1+(j-1)*LDA] = 0.0;
    }
  }

  for ( i = 1; i <= N; i++)
  {
    x[i-1] = ( double ) ( i );
  }

  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 0.0;
    for ( j = 1; j <= N; j++ )
    {
      b[i-1] = b[i-1] + a[i-1+(j-1)*LDA] * x[j-1];
    }
  }

  cout << "\n";
  cout << "  For an upper triangular matrix A,\n";
  cout << "  solve A * x = b\n";

  job = 01;

  info = dtrsl ( a, LDA, N, b, job );

  cout << "\n";
  cout << "  The solution (should be 1,2,3,4,5):\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout                       << "  "
         << setw(6)  << i      << "  "
         << setw(14) << b[i-1] << "\n";
  }

  for ( i = 1; i <= N; i++ )
  {
    b[i-1] = 0.0;
    for ( j = 1; j <= N; j++ )
    {
      b[i-1] = b[i-1] + a[j-1+(i-1)*LDA] * x[j-1];
    }
  }
  cout << "\n";
  cout << "  For an upper triangular matrix A,\n";
  cout << "  solve A' * x = b\n";

  job = 11;

  info = dtrsl ( a, LDA, N, b, job );

  cout << "\n";
  cout << "  The solution (should be 1,2,3,4,5):\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout                       << "  "
         << setw(6)  << i      << "  "
         << setw(14) << b[i-1] << "\n";
  }

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 fills a double precision array with unit pseudorandom values.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
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
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    P A Lewis, A S Goodman, J M Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and D_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2009
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
