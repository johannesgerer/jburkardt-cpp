# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "linpack_c.hpp"
# include "blas1_c.hpp"

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
void test32 ( );
void test33 ( );
void test34 ( );
void test345 ( );
void test35 ( );
void test36 ( );
void test37 ( );
complex <float> c4_uniform_01 ( int *seed );
complex <float> *c4mat_uniform_01 ( int m, int n, int *seed );
complex <float> *c4vec_uniform_01 ( int n, int *seed );
float r4_uniform_01 ( int *seed );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LINPACK_C_PRB.
//
//  Discussion:
//
//    LINPACK_C_PRB tests the LINPACK_C library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LINPACK_C_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LINPACK_C library.\n";

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
  test32 ( );
  test33 ( );
  test34 ( );
  test345 ( );
  test35 ( );
  test36 ( );
  test37 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LINPACK_C_PRB\n";
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
//    TEST01 tests CCHDC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  complex <float> c[N*N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int lda;
  lda = N;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  CCHDC computes the Cholesky decomposition.\n";
  cout << "\n";
  cout << "  The number of equations is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*lda] = complex <float> ( 2.5281,  0.0000 );
  a[1+0*lda] = complex <float> ( 2.1341,  0.2147 );
  a[2+0*lda] = complex <float> ( 2.4187, -0.2932 );

  a[0+1*lda] = complex <float> ( 2.1341, -0.2147 );
  a[1+1*lda] = complex <float> ( 3.0371,  0.0000 );
  a[2+1*lda] = complex <float> ( 2.0905, -1.1505 );

  a[0+2*lda] = complex <float> ( 2.4187,  0.2932 );
  a[1+2*lda] = complex <float> ( 2.0905,  1.1505 );
  a[2+2*lda] = complex <float> ( 2.7638,  0.0000 );

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
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

  info = cchdc ( a, lda, N, ipvt, job );

  if ( info != N )
  {
    cout << "\n";
    cout << "  CCHDC returned INFO = " << info << "\n";
    cout << "  The matrix is not Hermitian positive definite.\n";
    return;
  }
//
//  Zero out the lower diagonal.
//
  for ( i = 1; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      a[i+j*lda] = complex <float> ( 0.0, 0.0 );
    }
  }
//
//  Print the factorization.
//
  cout << "\n";
  cout << "  The Cholesky factor U:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Compute the Cholesky product.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = complex <float> ( 0.0, 0.0 );
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + conj ( a[k+i*lda] ) * a[k+j*lda];
      }
    }
  }

  cout << "\n";
  cout << "  The product U^H * U:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << c[i+j*N];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CCHEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define NZ 1

  complex <float> a[N*N];
  complex <float> b[N*N];
  float c[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int l;
  int lda;
  int ldz;
  complex <float> s[N];
  int seed;
  complex <float> z[N*NZ];

  lda = N;
  ldz = N;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  CCHEX can shift columns in a Cholesky factorization.\n";
  cout << "\n";
  cout << "  The number of equations is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*lda] = complex <float> ( 2.5281,  0.0000 );
  a[1+0*lda] = complex <float> ( 2.1341,  0.2147 );
  a[2+0*lda] = complex <float> ( 2.4187, -0.2932 );

  a[0+1*lda] = complex <float> ( 2.1341, -0.2147 );
  a[1+1*lda] = complex <float> ( 3.0371,  0.0000 );
  a[2+1*lda] = complex <float> ( 2.0905, -1.1505 );

  a[0+2*lda] = complex <float> ( 2.4187,  0.2932 );
  a[1+2*lda] = complex <float> ( 2.0905,  1.1505 );
  a[2+2*lda] = complex <float> ( 2.7638,  0.0000 );

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    z[i] = complex <float> ( i+1, 0.0 );
  }

  cout << "\n";
  cout << "  The vector Z:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << z[i] << "\n";
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

  info = cchdc ( a, lda, N, ipvt, job );

  if ( info != N )
  {
    cout << "\n";
    cout << "  CCHDC returned INFO = " << info << "\n";
    cout << "  This means the matrix is not positive definite.\n";
    return;
  }
//
//  Zero out the lower diagonal.
//
  for ( i = 1; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      a[i+j*lda] = complex <float> ( 0.0, 0.0 );
    }
  }
//
//  Print the factorization.
//
  cout << "\n";
  cout << "  The Cholesky factor U:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Right circular shift columns L through K.
//
  k = 1;
  l = 3;

  cout << "\n";
  cout << "  Right circular shift columns K  = " << k
       << " through L = " << l << "\n";

  job = 1;
  cchex ( a, lda, N, k, l, z, ldz, NZ, c, s, job );
//
//  Left circular shift columns K+1 through L.
//
  k = 2;
  l = 3;
  cout << "\n";
  cout << "  Left circular shift columns K = " << k
       << " through L = " << l << "\n";

  job = 2;
  cchex ( a, lda, N, k, l, z, ldz, NZ, c, s, job );
//
//  Print the factorization.
//
  cout << "\n";
  cout << "  The shifted Cholesky factor U:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The shifted vector Z:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << z[i] << "\n";
  }
//
//  Compute the Cholesky product.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      b[i+j*N] = complex <float> ( 0.0, 0.0 );
      for ( k = 0; k < N; k++ )
      {
        b[i+j*N] = b[i+j*N] + conj ( a[k+i*lda] ) * a[k+j*lda];
      }
    }
  }
  cout << "\n";
  cout << "  The shifted product U' * U:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << b[i+j*lda];
    }
    cout << "\n";
  }
  return;
# undef N
# undef NZ
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests CCHUD and CTRSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define P 20
# define NZ 1

  complex <float> b[P];
  complex <float> beta[P];
  float c[P];
  int i;
  int info;
  int j;
  int job;
  int ldr;
  int ldz;
  complex <float> r[P*P];
  float rho[NZ];
  complex <float> *row;
  complex <float> s[P];
  int seed;
  complex <float> x[P];
  complex <float> y[NZ];
  complex <float> z[P*NZ];

  ldr = P;
  ldz = P;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a complex Hermitian matrix\n";
  cout << "  CCHUD updates a Cholesky decomposition.\n";
  cout << "  CTRSL solves a triangular linear system.\n";
  cout << "\n";
  cout << "  In this example, we use CCHUD to solve a\n";
  cout << "  least squares problem R * b = z.\n";
  cout << "\n";
  cout << "  The number of equations is P = " << P << "\n";
//
//  Initialize.
//
  for ( j = 0; j < P; j++ )
  {
    for ( i = 0; i < P; i++ )
    {
      r[i+j*ldr] = complex <float> ( 0.0, 0.0 );
    }
  }
  for ( j = 0; j < NZ; j++ )
  {
    for ( i = 0; i < P; i++ )
    {
      z[i+j*ldz] = complex <float> ( 0.0, 0.0 );
    }
  }

  for ( i = 1; i <= P; i++ )
  {
    x[i-1] = complex <float> ( i, ( i % 2 ) );
  }
//
//  Use CCHUD to form R, Z and RHO by adding X and Y a row at a time.
//  X is a row of the least squares matrix and Y the right hand side.
//
  seed = 123456789;

  for ( i = 1; i <= P; i++ )
  {
    row = c4vec_uniform_01 ( P, &seed );
    y[0] = cdotu ( P, row, 1, x, 1 );
    rho[0] = 0.0;
    cchud ( r, ldr, P, row, z, ldz, NZ, y, rho, c, s );
    delete [] row;
  }
//
//  Generate the least squares solution, b = inverse ( R ) * Z.
//
  for ( j = 1; j <= NZ; j++ )
  {
    for ( i = 1; i <= P; i++ )
    {
      b[i-1] = z[i-1+(j-1)*ldz];
    }
    job = 1;

    info = ctrsl ( r, ldr, P, b, job );

    cout << "\n";
    cout << "  Solution vector # " << j << "\n";
    cout << "  (Should be (1,1) (2,0), (3,1) (4,0) ...)\n";
    cout << "\n";

    for ( i = 1; i <= P; i++ )
    {
      if ( i <= 5 || P-5 < i )
      {
        cout << "  " << setw(8) << i
             << "  " << setw(20) << b[i-1] << "\n";
      }
      if ( i == 5 )
      {
        cout << "  ......  ..............\n";
      }
    }
  }

  return;
# undef NZ
# undef P
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests CGBCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define ML 1
# define MU 1

  complex <float> a[(2*ML+MU+1)*N];
  complex <float> a_save[N*N];
  int i;
  int i1;
  int i2;
  int ipvt[N];
  int j;
  int k;
  int lda;
  int m;
  float rcond;
  int seed;

  lda = 2 * ML + MU + 1;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For a complex general band storage matrix:\n";
  cout << "  CGBCO factors the matrix and estimates the\n";
  cout << "  reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
  cout << "  The lower band is ML =  " << ML << "\n";
  cout << "  The upper band is MU =  " << MU << "\n";
//
//  Set the values of the matrix A.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a_save[i+j*N] = complex <float> ( 0.0, 0.0 );
    }
  }

  m = ML + MU + 1;

  seed = 123456789;

  for ( j = 1; j <= N; j++ )
  {
    i1 = i4_max ( 1, j - MU );
    i2 = i4_min ( N, j + ML );
    for ( i = i1; i <= i2; i++ )
    {
      k = i - j + m;
      a[k-1+(j-1)*lda] = c4_uniform_01 ( &seed );
      a_save[i-1+(j-1)*N] = a[k-1+(j-1)*lda];
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  rcond = cgbco ( a, lda, N, ML, MU, ipvt );

  cout << "\n";
  cout << "  Estimated reciprocal condition RCOND = " << rcond << "\n";

  return;
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
//    TEST05 tests CGBFA and CGBSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define ML 1
# define MU 1

  complex <float> a[(2*ML+MU+1)*N];
  complex <float> a_save[N*N];
  complex <float> b[N];
  int i;
  int i1;
  int i2;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int lda;
  int m;
  int seed;
  complex <float> *x;

  lda = 2 * ML + MU + 1;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For a complex general band storage matrix:\n";
  cout << "  CGBFA factors the matrix;\n";
  cout << "  CGBSL solves a factored linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N  << "\n";
  cout << "  The lower band is ML =  " << ML << "\n";
  cout << "  The upper band is MU =  " << MU << "\n";
//
//  Set the values of the matrix A.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a_save[i+j*N] = complex <float> ( 0.0, 0.0 );
    }
  }

  m = ML + MU + 1;

  seed = 123456789;

  for ( j = 1; j <= N; j++ )
  {
    i1 = i4_max ( 1, j - MU );
    i2 = i4_min ( N, j + ML );
    for ( i = i1; i <= i2; i++ )
    {
      k = i - j + m;
      a[k-1+(j-1)*lda] = c4_uniform_01 ( &seed );
      a_save[i-1+(j-1)*N] = a[k-1+(j-1)*lda];
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Set the values of the right hand side vector B.
//
  x = c4vec_uniform_01 ( N, &seed );

  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a_save[i+j*N] * x[j];
    }
  }
  cout << "\n";
  cout << "  The right hand side:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }
//
//  Factor the matrix A.
//
  info = cgbfa ( a, lda, N, ML, MU, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CGBFA returned INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  job = 0;
  cgbsl ( a, lda, N, ML, MU, ipvt, b, job );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i]
         << "  " << setw(20) << x[i] << "\n";;
  }

  delete [] x;

  return;
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
//    TEST06 tests CGBFA and CGBDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define ML 1
# define MU 1

  complex <float> a[(2*ML+MU+1)*N];
  complex <float> a_save[N*N];
  complex <float> det[2];
  int i;
  int i1;
  int i2;
  int info;
  int ipvt[N];
  int j;
  int k;
  int lda;
  int m;
  int seed;

  lda = 2 * ML + MU + 1;
  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a complex general band storage matrix:\n";
  cout << "  CGBFA factors the matrix.\n";
  cout << "  CGBDI computes the determinant.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N  << "\n";
  cout << "  The lower band is ML =  " << ML << "\n";
  cout << "  The upper band is MU =  " << MU << "\n";
//
//  Set the values of the matrix A.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a_save[i+j*N] = complex <float> ( 0.0, 0.0 );
    }
  }

  m = ML + MU + 1;

  seed = 123456789;

  for ( j = 1; j <= N; j++ )
  {
    i1 = i4_max ( 1, j - MU );
    i2 = i4_min ( N, j + ML );
    for ( i = i1; i <= i2; i++ )
    {
      k = i - j + m;
      a[k-1+(j-1)*lda] = c4_uniform_01 ( &seed );
      a_save[i-1+(j-1)*N] = a[k-1+(j-1)*lda];
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  info = cgbfa ( a, lda, N, ML, MU, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CGBFA returned INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  cgbdi ( a, lda, N, ML, MU, ipvt, det );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << "  * 10^ " << det[1] << "\n";

  return;
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests CGECO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> *a;
  int i;
  int ipvt[N];
  int j;
  int lda;
  float rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  For a complex general storage matrix:\n";
  cout << "  CGECO factors the matrix and estimates the\n";
  cout << "  reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c4mat_uniform_01 ( N, N, &seed );

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  rcond = cgeco ( a, lda, N, ipvt );

  cout << "\n";
  cout << "  Estimated reciprocal condition RCOND = " << rcond << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests CGEFA and CGESL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> *a;
  complex <float> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int lda;
  int seed;
  complex <float> *x;

  lda = N;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For a complex general storage matrix:\n";
  cout << "  CGEFA factors the matrix.\n";
  cout << "  CGESL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c4mat_uniform_01 ( N, N, &seed );

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*N];
    }
    cout << "\n";
  }
//
//  Set the values of the right hand side vector B.
//
  x = c4vec_uniform_01 ( N, &seed );

  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[i+j*N] * x[j];
    }
  }
  cout << "\n";
  cout << "  The right hand side:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }
//
//  Factor the matrix A.
//
  info = cgefa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CGEFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  job = 0;
  cgesl ( a, lda, N, ipvt, b, job );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i]
         << "  " << setw(20) << x[i] << "\n";;
  }

  delete [] a;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests CGEFA and CGEDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> *a;
  complex <float> a_save[N*N];
  complex <float> c[N*N];
  complex <float> det[2];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int lda;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For a complex general storage matrix:\n";
  cout << "  CGEFA factors the matrix.\n";
  cout << "  CGEDI computes the determinant or inverse.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c4mat_uniform_01 ( N, N, &seed );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_save[i+j*N] = a[i+j*N];
    }
  }
  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  info = cgefa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CGEFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  cgedi ( a, lda, N, ipvt, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inverse.
//
  job = 01;
  cgedi ( a, lda, N, ipvt, det, job );

  for ( i = 0; i < N; i++ )
  {
    for ( k = 0; k < N; k++ )
    {
      c[i+k*N] = 0.0;
      for ( j = 0; j < N; j++ )
      {
        c[i+k*N] = c[i+k*N] + a[i+j*N] * a_save[j+k*N];
      }
    }
  }
  cout << "\n";
  cout << "  The product inv(A) * A is\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << c[i+j*N];
    }
    cout << "\n";
  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests CGTSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  complex <float> b[N];
  complex <float> *c;
  complex <float> d[N];
  complex <float> *e;
  int i;
  int info;
  int seed;
  complex <float> x[N];

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For a complex tridiagonal matrix:\n";
  cout << "  CGTSL solves a linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  seed = 123456789;

  c = c4vec_uniform_01 ( N, &seed );
  c[0] = complex <float> ( 0.0, 0.0 );

  e = c4vec_uniform_01 ( N, &seed );
  e[N-1] = complex <float> ( 0.0, 0.0 );

  for ( i = 1; i <= N; i++ )
  {
    d[i-1] = 0.0;
    if ( i < N )
    {
      d[i-1] = d[i-1] - complex <float> ( 2.0 ) * e[i-1];
    }
    if ( 1 < i )
    {
      d[i-1] = d[i-1] - complex <float> ( 2.0 ) * c[i-1];
    }
  }
//
//  Set the desired solution
//
  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = complex <float> ( i, 10 * i );
  }
//
//  Compute the corresponding right hand side.
//
  b[0] = d[0] * x[0] + e[0] * x[1];
  for ( i = 2; i <= N-1; i++ )
  {
    b[i-1] = c[i-1] * x[i-2] + d[i-1] * x[i-1] + e[i-1] * x[i];
  }
  b[N-1] = c[N-1] * x[N-2] + d[N-1] * x[N-1];
//
//  Solve the tridiagonal system.
//
  info = cgtsl ( N, c, d, e, b );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(26) << b[i]
         << "  " << setw(26) << x[i] << "\n";
  }

  delete [] c;
  delete [] e;

  return;
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests CHICO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  int i;
  int ipvt[N];
  int j;
  int lda;
  float rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  For a complex Hermitian matrix:\n";
  cout << "  CHICO factors the matrix and estimates\n";
  cout << "  the reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = complex <float> ( r4_uniform_01 ( &seed ), 0.0 );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
      a[j+i*lda] = conj ( a[i+j*lda] );
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  rcond = chico ( a, lda, N, ipvt );

  cout << "\n";
  cout << "  Estimated reciprocal condition RCOND = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests CHIFA and CHISL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  complex <float> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int lda;
  int seed;
  complex <float> *x;

  lda = N;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  For a complex Hermitian matrix:\n";
  cout << "  CHIFA factors the matrix.\n";
  cout << "  CHISL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = complex <float> ( r4_uniform_01 ( &seed ), 0.0 );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
      a[j+i*lda] = conj ( a[i+j*lda] );
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*N];
    }
    cout << "\n";
  }
//
//  Set the values of the right hand side vector B.
//
  x = c4vec_uniform_01 ( N, &seed );

  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[i+j*N] * x[j];
    }
  }
  cout << "\n";
  cout << "  The right hand side:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }
//
//  Factor the matrix A.
//
  info = chifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CHIFA returned an error flag INFO = " << info << "\n";
    delete [] x;
    return;
  }
//
//  Solve the system.
//
  chisl ( a, lda, N, ipvt, b );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(26) << b[i]
         << "  " << setw(26) << x[i] << "\n";
  }

  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests CHIFA and CHIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  complex <float> a_save[N*N];
  complex <float> c[N*N];
  float det[2];
  int i;
  int inert[3];
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int lda;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  For a complex hermitian matrix:\n";
  cout << "  CHIFA factors the matrix.\n";
  cout << "  CHIDI computes the determinant, inverse,\n";
  cout << "  or inertia.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = complex <float> ( r4_uniform_01 ( &seed ), 0.0 );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
      a[j+i*lda] = conj ( a[i+j*lda] );
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_save[i+j*lda] = a[i+j*lda];
    }
  }
  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  info = chifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CHIFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  chidi ( a, lda, N, ipvt, det, inert, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inertia.
//
  job = 100;
  chidi ( a, lda, N, ipvt, det, inert, job );

  cout << "\n";
  cout << "  The inertia:\n";
  cout << "\n";

  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << inert[i] << "\n";
  }
//
//  Get the inverse.
//
  job = 1;
  chidi ( a, lda, N, ipvt, det, inert, job );
//
//  Only the upper triangle is set, so the user must set up the
//  lower triangle:
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      a[i+j*lda] = conj ( a[j+i*lda] );
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( k = 0; k < N; k++ )
    {
      c[i+k*lda] = 0.0;
      for ( j = 0; j < N; j++ )
      {
        c[i+k*lda] = c[i+k*lda] + a[i+j*lda] * a_save[j+k*lda];
      }
    }
  }
  cout << "\n";
  cout << "  The product inv(A) * A is\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << c[i+j*lda];
    }
    cout << "\n";
  }
  return;
# undef N
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests CHPCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> a_save[N*N];
  int i;
  int ipvt[N];
  int j;
  int k;
  float rcond;
  int seed;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  For a complex Hermitian matrix\n";
  cout << "  using packed storage,\n";
  cout << "  CHPCO factors the matrix and estimates\n";
  cout << "  the reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  k = 0;
  seed = 123456789;

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[k] = c4_uniform_01 ( &seed );
      a_save[i+j*N] = a[k];
      a_save[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
    a[k] = complex <float> ( r4_uniform_01 ( &seed ), 0.0 );
    a_save[j+j*N] = a[k];
    k = k + 1;
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  rcond = chpco ( a, N, ipvt );

  cout << "\n";
  cout << "  Estimated reciprocal condition RCOND = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests CHPFA and CHPSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> a_save[N*N];
  complex <float> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;
  int seed;
  complex <float> *x;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  For a complex Hermitian matrix,\n";
  cout << "  using packed storage,\n";
  cout << "  CHPFA factors the matrix.\n";
  cout << "  CHPSL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  k = 0;
  seed = 123456789;

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[k] = c4_uniform_01 ( &seed );
      a_save[i+j*N] = a[k];
      a_save[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
    a[k] = complex <float> ( r4_uniform_01 ( &seed ), 0.0 );
    a_save[j+j*N] = a[k];
    k = k + 1;
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Set the values of the right hand side vector B.
//
  x = c4vec_uniform_01 ( N, &seed );

  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a_save[i+j*N] * x[j];
    }
  }
  cout << "\n";
  cout << "  The right hand side:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }
//
//  Factor the matrix A.
//
  info = chpfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CHPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  chpsl ( a, N, ipvt, b );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(26) << b[i]
         << "  " << setw(26) << x[i] << "\n";
  }

  delete [] x;

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests CHPFA and CHPDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> a_save[N*N];
  complex <float> b[N*N];
  complex <float> c[N*N];
  float det[2];
  int i;
  int inert[3];
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int seed;

  cout << "\n";
  cout << "TEST16\n";
  cout << "  For a complex hermitian matrix,\n";
  cout << "  using packed storage,\n";
  cout << "  CHPFA factors the matrix.\n";
  cout << "  CHPDI computes the determinant, inverse,\n";
  cout << "  or inertia.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  k = 0;
  seed = 123456789;

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[k] = c4_uniform_01 ( &seed );
      a_save[i+j*N] = a[k];
      a_save[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
    a[k] = complex <float> ( r4_uniform_01 ( &seed ), 0.0 );
    a_save[j+j*N] = a[k];
    k = k + 1;
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  info = chpfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CHPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  chpdi ( a, N, ipvt, det, inert, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inertia.
//
  job = 100;
  chpdi ( a, N, ipvt, det, inert, job );

  cout << "\n";
  cout << "  The inertia:\n";
  cout << "\n";

  for ( i = 0; i < 3; i++ )
  {
    cout << "  " << setw(8) << inert[i] << "\n";
  }
//
//  Get the inverse.
//
  job = 1;
  chpdi ( a, N, ipvt, det, inert, job );
//
//  Only the upper triangle is set, so the user must set up the
//  lower triangle:
//
  k = 0;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      b[i+j*N] = a[k];
      b[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
    b[j+j*N] = a[k];
    k = k + 1;
  }

  for ( i = 0; i < N; i++ )
  {
    for ( k = 0; k < N; k++ )
    {
      c[i+k*N] = 0.0;
      for ( j = 0; j < N; j++ )
      {
        c[i+k*N] = c[i+k*N] + b[i+j*N] * a_save[j+k*N];
      }
    }
  }
  cout << "\n";
  cout << "  The product inv(A) * A is\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << c[i+j*N];
    }
    cout << "\n";
  }
  return;
# undef N
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests CPBCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define M 1

  complex <float> a[(M+1)*N];
  int i;
  int info;
  int lda;
  float rcond;

  lda = M + 1;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  For a complex positive definite hermitian band matrix,\n";
  cout << "  CPBCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  a[0+0*lda] = complex <float> ( 0.0000,  0.0000 );
  a[0+1*lda] = complex <float> ( 2.1341, -0.2147 );
  a[0+2*lda] = complex <float> ( 2.0905,  1.1505 );

  a[1+0*lda] = complex <float> ( 4.5281,  0.0000 );
  a[1+1*lda] = complex <float> ( 5.0371,  0.0000 );
  a[1+2*lda] = complex <float> ( 4.7638,  0.0000 );
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = cpbco ( a, lda, N, M, &info );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CPBCO returned INFO = " << info << "\n";
    cout << "  The factorization was not completed.\n";
    return;
  }

  cout << "\n";
  cout << "  Reciprocal condition  = " << rcond << "\n";

  return;
# undef N
# undef M
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests CPBDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define M 1

  complex <float> a[(M+1)*N];
  float det[2];
  int i;
  int info;
  int lda;

  lda = M + 1;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  For a complex positive definite hermitian band matrix,\n";
  cout << "  CPBDI computes the determinant as\n";
  cout << "    det = MANTISSA * 10^EXPONENT\n";
  cout << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  a[0+0*lda] = complex <float> ( 0.0000,  0.0000 );
  a[0+1*lda] = complex <float> ( 2.1341, -0.2147 );
  a[0+2*lda] = complex <float> ( 2.0905,  1.1505 );

  a[1+0*lda] = complex <float> ( 4.5281,  0.0000 );
  a[1+1*lda] = complex <float> ( 5.0371,  0.0000 );
  a[1+2*lda] = complex <float> ( 4.7638,  0.0000 );

  info = cpbfa ( a, lda, N, M );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error!  CPBFA returns INFO = " << info << "\n";
    return;
  }

  cpbdi ( a, lda, N, M, det );

  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests CPBFA and CPBSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define M 1

  complex <float> a[(M+1)*N];
  complex <float> b[N];
  int i;
  int info;
  int lda;

  lda = M + 1;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  For a complex positive definite hermitian band matrix,\n";
  cout << "  CPBFA computes the LU factors.\n";
  cout << "  CPBSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  a[0+0*lda] = complex <float> ( 0.0000,  0.0000 );
  a[0+1*lda] = complex <float> ( 2.1341, -0.2147 );
  a[0+2*lda] = complex <float> ( 2.0905,  1.1505 );

  a[1+0*lda] = complex <float> ( 4.5281,  0.0000 );
  a[1+1*lda] = complex <float> ( 5.0371,  0.0000 );
  a[1+2*lda] = complex <float> ( 4.7638,  0.0000 );
//
//  Set the right hand side.
//
  b[0] = complex <float> (  8.7963, -0.4294 );
  b[1] = complex <float> ( 18.4798,  3.6662 );
  b[2] = complex <float> ( 18.4724, -2.3010 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = cpbfa ( a, lda, N, M );

  if ( info != 0 )
  {
    cout << "  Error!  CPBFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  cpbsl ( a, lda, N, M, b );
//
//  Print the results.
//
  cout << "\n";
  cout << "  The solution:\n";
  cout << "  (Should be roughly (1,2,3)):\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }

  return;
# undef N
# undef M
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests CPOCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  int i;
  int info;
  int lda;
  float rcond;

  lda = N;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  CPOCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*3] = complex <float> ( 2.5281,  0.0000 );
  a[1+0*3] = complex <float> ( 2.1341,  0.2147 );
  a[2+0*3] = complex <float> ( 2.4187, -0.2932 );

  a[0+1*3] = complex <float> ( 2.1341, -0.2147 );
  a[1+1*3] = complex <float> ( 3.0371,  0.0000 );
  a[2+1*3] = complex <float> ( 2.0905, -1.1505 );

  a[0+2*3] = complex <float> ( 2.4187,  0.2932 );
  a[1+2*3] = complex <float> ( 2.0905,  1.1505 );
  a[2+2*3] = complex <float> ( 2.7638,  0.0000 );
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = cpoco ( a, lda, N, &info );

  cout << "\n";
  cout << "  Reciprocal condition  = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests CPOFA and CPODI.
//
//  Discussion:
//
//    CPOFA factors a positive definite symmetric matrix,
//    and CPODI can compute the determinant or the inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  float det[2];
  int i;
  int info;
  int j;
  int job;
  int lda;

  lda = N;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  CPOFA computes the LU factors,\n";
  cout << "  CPODI computes the inverse or determinant.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*3] = complex <float> ( 2.5281,  0.0000 );
  a[1+0*3] = complex <float> ( 2.1341,  0.2147 );
  a[2+0*3] = complex <float> ( 2.4187, -0.2932 );

  a[0+1*3] = complex <float> ( 2.1341, -0.2147 );
  a[1+1*3] = complex <float> ( 3.0371,  0.0000 );
  a[2+1*3] = complex <float> ( 2.0905, -1.1505 );

  a[0+2*3] = complex <float> ( 2.4187,  0.2932 );
  a[1+2*3] = complex <float> ( 2.0905,  1.1505 );
  a[2+2*3] = complex <float> ( 2.7638,  0.0000 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = cpofa ( a, lda, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, CPOFA returns INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant and inverse.
//
  cout << "\n";
  cout << "  Get the determinant and inverse.\n";

  job = 11;
  cpodi ( a, lda, N, det, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Determinant  = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  CPODI produces only the 'upper half triangle' of the inverse,
//  which is actually symmetric.  Thus, the lower half could be
//  produced by copying from the upper half.  However, the first row
//  of A, as returned, is exactly the first row of the inverse.
//
  cout << "\n";
  cout << "  First row of inverse:\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(20) << a[0+j*N];
  }
  cout << "\n";

  return;
# undef N
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests CPOFA and CPOSL.
//
//  Discussion:
//
//    CPOFA factors a Hermitian positive definite matrix,
//    and CPOSL can solve a factored linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  complex <float> b[N];
  int i;
  int info;
  int j;
  int lda;
  complex <float> x[N];

  lda = N;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  CPOFA computes the LU factors.\n";
  cout << "  CPOSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*3] = complex <float> ( 2.5281,  0.0000 );
  a[1+0*3] = complex <float> ( 2.1341,  0.2147 );
  a[2+0*3] = complex <float> ( 2.4187, -0.2932 );

  a[0+1*3] = complex <float> ( 2.1341, -0.2147 );
  a[1+1*3] = complex <float> ( 3.0371,  0.0000 );
  a[2+1*3] = complex <float> ( 2.0905, -1.1505 );

  a[0+2*3] = complex <float> ( 2.4187,  0.2932 );
  a[1+2*3] = complex <float> ( 2.0905,  1.1505 );
  a[2+2*3] = complex <float> ( 2.7638,  0.0000 );
//
//  Set the right hand side.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( 2 * i + 1, 2 * i + 2 );
  }
  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[i+j*N] * x[j];
    }
  }
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = cpofa ( a, lda, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, CPOFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  cposl ( a, lda, N, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The solution:\n";
  cout << "  (Should be (1+2i),(3+4i),(5+6i):\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
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
//    TEST23 tests CPPCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  int info;
  float rcond;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite packed matrix,\n";
  cout << "  CPPCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0] = complex <float> ( 2.5281,  0.0000 );

  a[1] = complex <float> ( 2.1341, -0.2147 );
  a[2] = complex <float> ( 3.0371,  0.0000 );

  a[3] = complex <float> ( 2.4187,  0.2932 );
  a[4] = complex <float> ( 2.0905,  1.1505 );
  a[5] = complex <float> ( 2.7638,  0.0000 );
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition number.\n";

  rcond = cppco ( a, N, &info );

  cout << "\n";
  cout << "  Reciprocal condition number = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests CPPFA and CPPDI.
//
//  Discussion:
//
//    CPPFA factors a Hermitian positive definite packed matrix,
//    and CPPDI can compute the determinant or the inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> b[N*N];
  float det[2];
  int i;
  int info;
  int j;
  int job;
  int k;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite packed matrix,\n";
  cout << "  CPPFA factors the matrix.\n";
  cout << "  CPPDI computes the inverse or determinant.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0] = complex <float> ( 2.5281,  0.0000 );

  a[1] = complex <float> ( 2.1341, -0.2147 );
  a[2] = complex <float> ( 3.0371,  0.0000 );

  a[3] = complex <float> ( 2.4187,  0.2932 );
  a[4] = complex <float> ( 2.0905,  1.1505 );
  a[5] = complex <float> ( 2.7638,  0.0000 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = cppfa ( a, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, CPPFA returns INFO = " << info << "\n";
    return;
  }
//
//  Invert the matrix.
//
  cout << "\n";
  cout << "  Get the determinant and inverse.\n";

  job = 11;
  cppdi ( a, N, det, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Determinant  = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  CPPDI produces only the 'upper half triangle' of the inverse,
//  which is actually symmetric.  Thus, the lower half could be
//  produced by copying from the upper half.
//
  k = 0;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      b[i+j*N] = a[k];
      b[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
  }

  cout << "\n";
  cout << "  Inverse:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << b[i+j*N];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests CPPFA and CPPSL.
//
//  Discussion:
//
//    CPOFA factors a Hermitian positive definite packed matrix,
//    and CPOSL can solve a factored linear system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> b[N];
  int i;
  int info;
  int j;
  int k;
  complex <float> x[N];

  cout << "\n";
  cout << "TEST25\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite packed matrix,\n";
  cout << "  CPPFA factors the matrix.\n";
  cout << "  CPPSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0] = complex <float> ( 2.5281,  0.0000 );

  a[1] = complex <float> ( 2.1341, -0.2147 );
  a[2] = complex <float> ( 3.0371,  0.0000 );

  a[3] = complex <float> ( 2.4187,  0.2932 );
  a[4] = complex <float> ( 2.0905,  1.1505 );
  a[5] = complex <float> ( 2.7638,  0.0000 );

  b[0] = complex <float> ( 20.12350, 28.92670 );
  b[1] = complex <float> ( 14.36550, 34.92680 );
  b[2] = complex <float> ( 27.69760, 26.03750 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = cppfa ( a, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, CPPFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  cppsl ( a, N, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The solution:\n";
  cout << "  (Should be (1+2i),(3+4i),(5+6i):\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests CPTSL.
//
//  Discussion:
//
//    CPTSL factors and solves a Hermitian positive definite
//    tridiagonal system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> b[N];
  complex <float> d[N];
  complex <float> e[N];
  int i;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite tridiagonal matrix,\n";
  cout << "  CPTSL factors and solves a linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  e[0] = complex <float> ( 2.1341, -0.2147 );
  e[1] = complex <float> ( 2.0905,  1.1505 );
  e[2] = complex <float> ( 0.0000,  0.0000 );

  d[0] = complex <float> ( 4.5281,  0.0000 );
  d[1] = complex <float> ( 5.0371,  0.0000 );
  d[2] = complex <float> ( 4.7638,  0.0000 );
//
//  Set the right hand side.
//
  b[0] = complex <float> (  8.7963, -0.4294 );
  b[1] = complex <float> ( 18.4798,  3.6662 );
  b[2] = complex <float> ( 18.4724, -2.3010 );
//
//  Factor and solve the system.
//
  cout << "\n";
  cout << "  Factor the matrix and solve the system.\n";

  cptsl ( N, d, e, b );
//
//  Print the result.
//
  cout << "\n";
  cout << "  The solution:\n";
  cout << "  (Should be roughly (1,2,3)):\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(26) << b[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests CQRDC and CQRSL.
//
//  Discussion:
//
//    CQRDC and CQRSL compute the QR factorization, and use it
//    to solve linear systems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define P 3

  complex <float> *a;
  complex <float> b[N*P];
  int i;
  int info;
  int ipvt[P];
  int j;
  int job;
  int k;
  int lda;
  complex <float> q[N*N];
  complex <float> qraux[P];
  complex <float> qty[N];
  complex <float> qy[N];
  complex <float> r[N*P];
  complex <float> rsd[N];
  int seed;
  complex <float> xb[N];
  complex <float> y[N];

  lda = N;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  For a complex general matrix,\n";
  cout << "  CQRDC computes the QR decomposition of a\n";
  cout << "  matrix, but does not return Q and R explicitly.\n";
  cout << "\n";
  cout << "  Show how Q and R can be recovered using CQRSL.\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c4mat_uniform_01 ( N, P, &seed );

  cout << "\n";
  cout << "  The matrix A is\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Decompose the matrix.
//
  cout << "\n";
  cout << "  Decompose the matrix.\n";

  job = 0;
  for ( i = 0; i < P; i++ )
  {
    ipvt[i] = 0;
  }

  cqrdc ( a, lda, N, P, qraux, ipvt, job );
//
//  Print out what CQRDC has stored in A...
//
  cout << "\n";
  cout << "  The packed matrix A which describes Q and R:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
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

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << qraux[i] << "\n";
  }
//
//  Print out the resulting R factor.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      if ( j < i )
      {
        r[i+j*lda] = complex <float> ( 0.0, 0.0 );
      }
      else
      {
        r[i+j*lda] = a[i+j*lda];
      }
    }
  }
  cout << "\n";
  cout << "  The R factor:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      cout << "  " << setw(20) << r[i+j*N];
    }
    cout << "\n";
  }
//
//  Call CQRSL to extract the information about the Q matrix.
//  We do this, essentially, by asking CQRSL to tell us the
//  value of Q*Y, where Y is a column of the identity matrix.
//
  job = 10000;

  for ( j = 0; j < N; j++ )
  {
//
//  Set the vector Y.
//
    for ( i = 0; i < N; i++ )
    {
      y[i] = complex <float> ( 0.0, 0.0 );
    }
    y[j] = complex <float> ( 1.0, 0.0 );
//
//  Ask CQRSL to tell us what Q*Y is.
//
    info = cqrsl ( a, lda, N, P, qraux, y, qy, qty, b, rsd, xb, job );

    if ( info != 0 )
    {
      cout << "  Error!  CQRSL returns INFO = " << info << "\n";
      return;
    }
//
//  Copy QY into the appropriate column of Q.
//
    for ( i = 0; i < N; i++ )
    {
      q[i+j*N] = qy[i];
    }
  }
//
//  Now print out the Q matrix we have extracted.
//
  cout << "\n";
  cout << "  The Q factor:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << q[i+j*N];
    }
    cout << "\n";
  }
//
//  Compute Q*R to verify that it equals A.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      b[i+j*N] = complex <float> ( 0.0, 0.0 );
      for ( k = 0; k < N; k++ )
      {
        b[i+j*N] = b[i+j*N] + q[i+k*N] * r[k+j*N];
      }
    }
    cout << "\n";
  }
//
//  Print the result.
//
  cout << "\n";
  cout << "  The product Q * R:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < P; j++ )
    {
      cout << "  " << setw(20) << b[i+j*N];
    }
    cout << "\n";
  }

  delete [] a;

  return;
# undef N
# undef P
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests CSICO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  int i;
  int ipvt[N];
  int j;
  int lda;
  float rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST28\n";
  cout << "  For a complex symmetric matrix:\n";
  cout << "  CSICO factors the matrix and estimates\n";
  cout << "  the reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = c4_uniform_01 ( &seed );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
      a[j+i*lda] = a[i+j*lda];
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  rcond = csico ( a, lda, N, ipvt );

  cout << "\n";
  cout << "  Estimated reciprocal condition RCOND = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests CSIFA and CSISL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  complex <float> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int lda;
  int seed;
  complex <float> x[N];

  lda = N;

  cout << "\n";
  cout << "TEST29\n";
  cout << "  For a complex symmetric matrix:\n";
  cout << "  CSIFA factors the matrix.\n";
  cout << "  CSISL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = c4_uniform_01 ( &seed );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
      a[j+i*lda] = a[i+j*lda];
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Set the values of the right hand side vector B.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = c4_uniform_01 ( &seed );
  }

  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[i+j*lda] * x[j];
    }
  }
  cout << "\n";
  cout << "  The right hand side:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }
//
//  Factor the matrix A.
//
  info = csifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CSIFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  csisl ( a, lda, N, ipvt, b );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i]
         << "  " << setw(20) << x[i] << "\n";;
  }

  return;
# undef N
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests CSIFA and CSIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  complex <float> a_save[N*N];
  complex <float> c[N*N];
  complex <float> det[2];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int lda;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST30\n";
  cout << "  For a complex symmetric matrix:\n";
  cout << "  CSIFA factors the matrix.\n";
  cout << "  CSIDI computes the determinant or inverse.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = c4_uniform_01 ( &seed );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
      a[j+i*lda] = a[i+j*lda];
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_save[i+j*lda] = a[i+j*lda];
    }
  }
  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  info = csifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CSIFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  csidi ( a, lda, N, ipvt, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inverse.
//
  job = 1;
  csidi ( a, lda, N, ipvt, det, job );
//
//  Only the upper triangle is set, so the user must set up the
//  lower triangle:
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      a[i+j*lda] = a[j+i*lda];
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( k = 0; k < N; k++ )
    {
      c[i+k*N] = 0.0;
      for ( j = 0; j < N; j++ )
      {
        c[i+k*N] = c[i+k*N] + a[i+j*N] * a_save[j+k*N];
      }
    }
  }
  cout << "\n";
  cout << "  The product inv(A) * A is\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << c[i+j*N];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests CSPCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> a_save[N*N];
  int i;
  int ipvt[N];
  int j;
  int k;
  float rcond;
  int seed;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  For a complex symmetric matrix\n";
  cout << "  in packed storage,\n";
  cout << "  CSPCO factors the matrix and estimates\n";
  cout << "  the reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the packed matrix A.
//
  k = 0;
  seed = 123456789;

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[k] = c4_uniform_01 ( &seed );
      k = k + 1;
    }
    a[k] = c4_uniform_01 ( &seed );
    k = k + 1;
  }
//
//  Copy the packed matrix into a "normal" matrix.
//
  k = 0;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a_save[i+j*N] = a[k];
      k = k + 1;
    }
  }

  for ( j = 0; j < N; j++ )
  {
    for ( i = j+1; i < N; i++ )
    {
      a_save[i+j*N] = a_save[j+i*N];
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  rcond = cspco ( a, N, ipvt );

  cout << "\n";
  cout << "  Estimated reciprocal condition RCOND = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests CSPFA and CSPSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> a_save[N*N];
  complex <float> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;
  int seed;
  complex <float> *x;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  For a complex symmetric matrix\n";
  cout << "  in packed storage,\n";
  cout << "  CSPFA factors the matrix.\n";
  cout << "  CSPSL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the packed matrix A.
//
  k = 0;
  seed = 123456789;

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[k] = c4_uniform_01 ( &seed );
      k = k + 1;
    }
    a[k] = c4_uniform_01 ( &seed );
    k = k + 1;
  }
//
//  Copy the packed matrix into a "normal" matrix.
//
  k = 0;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a_save[i+j*N] = a[k];
      k = k + 1;
    }
  }

  for ( j = 0; j < N; j++ )
  {
    for ( i = j+1; i < N; i++ )
    {
      a_save[i+j*N] = a_save[j+i*N];
    }
  }

  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Set the values of the right hand side vector B.
//
  x = c4vec_uniform_01 ( N, &seed );

  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a_save[i+j*N] * x[j];
    }
  }
  cout << "\n";
  cout << "  The right hand side:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i] << "\n";
  }
//
//  Factor the matrix A.
//
  info = cspfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CSPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  cspsl ( a, N, ipvt, b );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(20) << b[i]
         << "  " << setw(20) << x[i] << "\n";;
  }

  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests CSPFA and CSPDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[(N*(N+1))/2];
  complex <float> a_save[N*N];
  complex <float> b_save[N*N];
  complex <float> c[N*N];
  complex <float> det[2];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int seed;

  cout << "\n";
  cout << "TEST33\n";
  cout << "  For a complex symmetric matrix\n";
  cout << "  in packed storage,\n";
  cout << "  CSPFA factors the matrix.\n";
  cout << "  CSPDI computes the determinant or inverse.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the packed matrix A.
//
  k = 0;
  seed = 123456789;

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[k] = c4_uniform_01 ( &seed );
      k = k + 1;
    }
    a[k] = c4_uniform_01 ( &seed );
    k = k + 1;
  }
//
//  Copy the packed matrix into a "normal" matrix.
//
  k = 0;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a_save[i+j*N] = a[k];
      k = k + 1;
    }
  }

  for ( j = 0; j < N; j++ )
  {
    for ( i = j+1; i < N; i++ )
    {
      a_save[i+j*N] = a_save[j+i*N];
    }
  }
  cout << "\n";
  cout << "  The matrix:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a_save[i+j*N];
    }
    cout << "\n";
  }
//
//  Factor the matrix A.
//
  info = cspfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  CSPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  cspdi ( a, N, ipvt, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inverse.
//
  job = 1;
  cspdi ( a, N, ipvt, det, job );
//
//  Copy the packed matrix into a "normal" matrix.
//
  k = 0;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      b_save[i+j*N] = a[k];
      k = k + 1;
    }
  }
  for ( j = 0; j < N; j++ )
  {
    for ( i = j+1; i < N; i++ )
    {
      b_save[i+j*N] = b_save[j+i*N];
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( k = 0; k < N; k++ )
    {
      c[i+k*N] = 0.0;
      for ( j = 0; j < N; j++ )
      {
        c[i+k*N] = c[i+k*N] + a_save[i+j*N] * b_save[j+k*N];
      }
    }
  }
  cout << "\n";
  cout << "  The product inv(A) * A is\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << c[i+j*N];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests CSVDC.
//
//  Discussion:
//
//    CSVDC computes the singular value decomposition:
//
//      A = U * S * conjg-transpose ( V )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 3

  complex <float> *a;
//
//  E must be dimensioned at least the maximum of M+1 and N.
//
  complex <float> e[M+N];
  int i;
  int info;
  int j;
  int lda;
  int ldu;
  int ldv;
  int job;
  int k;
//
//  S must be dimensioned at least the maximum of M+1 and N.
//
  complex <float> s[M+N];
  int seed;
  complex <float> sigma[M*N];
  complex <float> u[M*M];
  complex <float> v[N*N];

  cout << "\n";
  cout << "TEST34\n";
  cout << "  For an MxN matrix A in complex general storage,\n";
  cout << "  CSVDC computes the singular value decomposition:\n";
  cout << "    A = U * S * V^H\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
//
//  Set A.
//
  seed = 123456789;
  lda = M;

  a = c4mat_uniform_01 ( M, N, &seed );

  cout << "\n";
  cout << "  The matrix A:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Decompose the matrix.
//
  cout << "\n";
  cout << "  Decompose the matrix.\n";

  job = 11;
  ldu = M;
  ldv = N;

  info = csvdc ( a, lda, M, N, s, e, u, ldu, v, ldv, job );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning:\n";
    cout << "  CSVDC returned nonzero INFO = " << info << "\n";
    return;
  }

  cout << "\n";
  cout << "  Singular values:\n";
  cout << "\n";

  for ( i = 0; i < i4_min ( M, N ); i++ )
  {
    cout << "  " << setw(8) << i + 1
         << setw(20) << s[i] << "\n";
  }

  cout << "\n";
  cout << "  Left Singular Vector Matrix U:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < M; j++ )
    {
      cout << "  " << setw(20) << u[i+j*ldu];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Right Singular Vector Matrix V:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << v[i+j*ldv];
    }
    cout << "\n";
  }
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*lda] = complex <float> ( 0.0, 0.0 );
      for ( k = 0; k < i4_min ( M, N ); k++ )
      {
        a[i+j*lda] =  a[i+j*lda] + u[i+k*ldu] * s[k] * conj ( v[j+k*ldv] );
      }
    }
  }

  cout << "\n";
  cout << "  The product U * S * V^H (should equal A):\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test345 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST345 tests CSVDC.
//
//  Discussion:
//
//    CSVDC computes the singular value decomposition:
//
//      A = U * S * conjg-transpose ( V )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2011
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 4

  complex <float> *a;
//
//  E must be dimensioned at least the maximum of M+1 and N.
//
  complex <float> e[M+N];
  complex <float> I;
  int i;
  int info;
  int j;
  int lda;
  int ldu;
  int ldv;
  int job;
  int k;
//
//  S must be dimensioned at least the maximum of M+1 and N.
//
  complex <float> s[M+N];
  int seed;
  complex <float> sigma[M*N];
  complex <float> u[M*M];
  complex <float> v[N*N];

  cout << "\n";
  cout << "TEST345\n";
  cout << "  For an MxN matrix A in complex general storage,\n";
  cout << "  CSVDC computes the singular value decomposition:\n";
  cout << "    A = U * S * V^H\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
//
//  Set A.
//
  I = complex <float> ( 0.0, 1.0 );

  lda = M;
  a = new complex <float>[M*N];

  a[0+0*M] =   1.0;
  a[1+0*M] = - I;
  a[2+0*M] = - 1.0;
  a[3+0*M] =   I;

  a[0+1*M] =   1.0;
  a[1+1*M] = - 1.0;
  a[2+1*M] = - 1.0;
  a[3+1*M] =   1.0;

  a[0+2*M] =   1.0;
  a[1+2*M] =   1.0;
  a[2+2*M] =   1.0;
  a[3+2*M] =   1.0;

  a[0+3*M] =   1.0;
  a[1+3*M] =   I;
  a[2+3*M] = - 1.0;
  a[3+3*M] = - I;

  cout << "\n";
  cout << "  The matrix A:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }
//
//  Decompose the matrix.
//
  cout << "\n";
  cout << "  Decompose the matrix.\n";

  job = 11;
  ldu = M;
  ldv = N;

  info = csvdc ( a, lda, M, N, s, e, u, ldu, v, ldv, job );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning:\n";
    cout << "  CSVDC returned nonzero INFO = " << info << "\n";
    return;
  }

  cout << "\n";
  cout << "  Singular values:\n";
  cout << "\n";

  for ( i = 0; i < i4_min ( M, N ); i++ )
  {
    cout << "  " << setw(8) << i + 1
         << setw(20) << s[i] << "\n";
  }

  cout << "\n";
  cout << "  Left Singular Vector Matrix U:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < M; j++ )
    {
      cout << "  " << setw(20) << u[i+j*ldu];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Right Singular Vector Matrix V:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << v[i+j*ldv];
    }
    cout << "\n";
  }
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*lda] = complex <float> ( 0.0, 0.0 );
      for ( k = 0; k < i4_min ( M, N ); k++ )
      {
        a[i+j*lda] =  a[i+j*lda] + u[i+k*ldu] * s[k] * conj ( v[j+k*ldv] );
      }
    }
  }

  cout << "\n";
  cout << "  The product U * S * V^H (should equal A):\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << fixed << setw(20) << a[i+j*lda];
    }
    cout << "\n";
  }

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 tests CTRCO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  int i;
  int j;
  int job;
  int lda;
  float rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST35\n";
  cout << "  For a complex triangular matrix,\n";
  cout << "  CTRCO estimates the condition.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
    }
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = complex <float> ( 0.0, 0.0 );
    }
  }
//
//  Get the condition of the lower triangular matrix.
//
  job = 0;
  rcond = ctrco ( a, lda, N, job );

  cout << "\n";
  cout << "  Estimated reciprocal condition RCOND = " << rcond << "\n";

  return;
# undef N
}
//****************************************************************************80

void test36 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST36 tests CTRDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  complex <float> a[N*N];
  complex <float> a_save[N*N];
  complex <float> c[N*N];
  complex <float> det[2];
  int i;
  int info;
  int j;
  int job;
  int k;
  int lda;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST36\n";
  cout << "  For a complex triangular matrix,\n";
  cout << "  CTRDI computes the determinant or inverse.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      a[i+j*lda] = c4_uniform_01 ( &seed );
    }
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = complex <float> ( 0.0, 0.0 );
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_save[i+j*lda] = a[i+j*lda];
    }
  }
//
//  Get the determinant of the lower triangular matrix.
//
  job = 100;
  info = ctrdi ( a, lda, N, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << real ( det[1] ) << "\n";
//
//  Get the inverse of the lower triangular matrix.
//
  job = 10;
  info = ctrdi ( a, lda, N, det, job );

  for ( i = 0; i < N; i++ )
  {
    for ( k = 0; k < N; k++ )
    {
      c[i+k*lda] = 0.0;
      for ( j = 0; j < N; j++ )
      {
        c[i+k*lda] = c[i+k*lda] + a[i+j*lda] * a_save[j+k*lda];
      }
    }
  }
  cout << "\n";
  cout << "  The product inv(A) * A is\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(20) << c[i+j*lda];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test37 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST37 tests CTRSL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  complex <float> a[N*N];
  complex <float> b[N];
  int i;
  int info;
  int j;
  int job;
  int k;
  int lda;
  int seed;
  complex <float> x[N];

  lda = N;

  cout << "\n";
  cout << "TEST37\n";
  cout << "  For a complex triangular matrix,\n";
  cout << "  CTRSL solves a linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  seed = 123456789;

  k = 0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      k = k + 1;
      a[i+j*lda] = c4_uniform_01 ( &seed );
    }
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = complex <float> ( 0.0, 0.0 );
    }
  }
//
//  Set the desired solution
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <float> ( i + 1, 10 * ( i + 1 ) );
  }
//
//  Compute the corresponding right hand side.
//
  for ( i = 0; i < N; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < N; j++ )
    {
      b[i] = b[i] + a[i+j*N] * x[j];
    }
  }
//
//  Solve the lower triangular system.
//
  job = 0;
  info = ctrsl ( a, lda, N, b, job );

  cout << "\n";
  cout << "  Computed                     Exact\n";
  cout << "  Solution                     Solution\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(26) << b[i]
         << "  " << setw(26) << x[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

complex <float> c4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4_UNIFORM_01 returns a unit complex pseudorandom number.
//
//  Discussion:
//
//    The angle should be uniformly distributed between 0 and 2 * PI,
//    the square root of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4_UNIFORM_01, a pseudorandom complex value.
//
{
  float r;
  int k;
  float pi = 3.1415926E+00;
  float theta;
  complex <float> value;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = sqrt ( ( float ) ( ( double ) ( *seed ) * 4.656612875E-10 ) );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  theta = 2.0 * pi * ( float )
    ( ( double ) ( *seed ) * 4.656612875E-10 );

  value = complex <float> ( r * cos ( theta ), r * sin ( theta ) );

  return value;
}
//****************************************************************************80

complex <float> *c4mat_uniform_01 ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01 returns a unit complex pseudorandom matrix.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex C4MAT_UNIFORM_01[M*N], the pseudorandom complex matrix.
//
{
  complex <float> *c;
  int i;
  int j;
  float r;
  int k;
  float pi = 3.1415926;
  float theta;

  c = new complex <float>[m*n];

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

      r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * complex <float> ( cos ( theta ), sin ( theta ) );
    }
  }

  return c;
}
//****************************************************************************80

complex <float> *c4vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01 returns a unit complex pseudorandom vector.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values to compute.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4VEC_UNIFORM_01[N], the pseudorandom
//    complex vector.
//
{
  complex <float> *c;
  int i;
  float r;
  int k;
  float pi = 3.1415926;
  float theta;

  c = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return c;
}
//****************************************************************************80

float r4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a real pseudorandom number.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r4_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R4_UNIFORM_01
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
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
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
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float r;

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
  r = ( float ) ( *seed ) * 4.656612875E-10;

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
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
