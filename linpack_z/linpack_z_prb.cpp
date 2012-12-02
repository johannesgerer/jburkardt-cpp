# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "linpack_z.hpp"
# include "blas1_z.hpp"

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
complex <double> c8_uniform_01 ( int *seed );
complex <double> *c8mat_uniform_01 ( int m, int n, int *seed );
complex <double> *c8vec_uniform_01 ( int n, int *seed );
double r8_uniform_01 ( int *seed );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    LINPACK_Z_PRB tests the double precision complex LINPACK routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "LINPACK_Z_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LINPACK_Z library.\n";

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
  cout << "LINPACK_Z_PRB\n";
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
//    TEST01 tests ZCHDC.
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

  complex <double> a[N*N];
  complex <double> c[N*N];
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
  cout << "  ZCHDC computes the Cholesky decomposition.\n";
  cout << "\n";
  cout << "  The number of equations is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*lda] = complex <double> ( 2.5281,  0.0000 );
  a[1+0*lda] = complex <double> ( 2.1341,  0.2147 );
  a[2+0*lda] = complex <double> ( 2.4187, -0.2932 );

  a[0+1*lda] = complex <double> ( 2.1341, -0.2147 );
  a[1+1*lda] = complex <double> ( 3.0371,  0.0000 );
  a[2+1*lda] = complex <double> ( 2.0905, -1.1505 );

  a[0+2*lda] = complex <double> ( 2.4187,  0.2932 );
  a[1+2*lda] = complex <double> ( 2.0905,  1.1505 );
  a[2+2*lda] = complex <double> ( 2.7638,  0.0000 );

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

  info = zchdc ( a, lda, N, ipvt, job );

  if ( info != N )
  {
    cout << "\n";
    cout << "  ZCHDC returned INFO = " << info << "\n";
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
      a[i+j*lda] = complex <double> ( 0.0, 0.0 );
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
      c[i+j*N] = complex <double> ( 0.0, 0.0 );
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
//    TEST02 tests ZCHEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define NZ 1

  complex <double> a[N*N];
  complex <double> b[N*N];
  double c[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int k;
  int l;
  int lda;
  int ldz;
  complex <double> s[N];
  int seed;
  complex <double> z[N*NZ];

  lda = N;
  ldz = N;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  ZCHEX can shift rows and columns in a Cholesky factorization.\n";
  cout << "\n";
  cout << "  The number of equations is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*lda] = complex <double> ( 2.5281,  0.0000 );
  a[1+0*lda] = complex <double> ( 2.1341,  0.2147 );
  a[2+0*lda] = complex <double> ( 2.4187, -0.2932 );

  a[0+1*lda] = complex <double> ( 2.1341, -0.2147 );
  a[1+1*lda] = complex <double> ( 3.0371,  0.0000 );
  a[2+1*lda] = complex <double> ( 2.0905, -1.1505 );

  a[0+2*lda] = complex <double> ( 2.4187,  0.2932 );
  a[1+2*lda] = complex <double> ( 2.0905,  1.1505 );
  a[2+2*lda] = complex <double> ( 2.7638,  0.0000 );

  cout << "\n";
  cout << "  The matrix A:\n";
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
    z[i] = complex <double> ( i+1, 0.0 );
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

  info = zchdc ( a, lda, N, ipvt, job );

  if ( info != N )
  {
    cout << "\n";
    cout << "  ZCHDC returned INFO = " << info << "\n";
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
      a[i+j*lda] = complex <double> ( 0.0, 0.0 );
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
  cout << "  Right circular shift rows and columns K  = " << k
       << " through L = " << l << "\n";
  cout << "\n";
  cout << "  Logical matrix is now:\n";
  cout << "\n";
  cout << "  33 31 32\n";
  cout << "  13 11 12\n";
  cout << "  23 21 22\n";

  job = 1;
  zchex ( a, lda, N, k, l, z, ldz, NZ, c, s, job );
//
//  Left circular shift columns K+1 through L.
//
  cout << "\n";
  cout << "  Left circular shift rows and columns K+1 = " << k + 1
       << " through L = " << l << "\n";
  cout << "\n";
  cout << "  Logical matrix is now:\n";
  cout << "\n";
  cout << "  33 32 31\n";
  cout << "  23 22 21\n";
  cout << "  13 12 11\n";

  job = 2;
  zchex ( a, lda, N, k+1, l, z, ldz, NZ, c, s, job );
//
//  Print the factorization.
//
  cout << "\n";
  cout << "  The shifted Cholesky factor UU:\n";
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
  cout << "  The shifted vector ZZ:\n";
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
      b[i+j*N] = complex <double> ( 0.0, 0.0 );
      for ( k = 0; k < N; k++ )
      {
        b[i+j*N] = b[i+j*N] + conj ( a[k+i*lda] ) * a[k+j*lda];
      }
    }
  }
  cout << "\n";
  cout << "  The shifted product AA = UU' * UU:\n";
  cout << "  The rows and columns of the original matrix A reappear,\n";
  cout << "  but in reverse order.\n";
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
//    TEST03 tests ZCHUD and ZTRSL.
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
# define NZ 1

  complex <double> b[P];
  complex <double> beta[P];
  double c[P];
  int i;
  int info;
  int j;
  int job;
  int ldr;
  int ldz;
  complex <double> r[P*P];
  double rho[NZ];
  complex <double> *row;
  complex <double> s[P];
  int seed;
  complex <double> x[P];
  complex <double> y[NZ];
  complex <double> z[P*NZ];

  ldr = P;
  ldz = P;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a complex Hermitian matrix\n";
  cout << "  ZCHUD updates a Cholesky decomposition.\n";
  cout << "  ZTRSL solves a triangular linear system.\n";
  cout << "\n";
  cout << "  In this example, we use ZCHUD to solve a\n";
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
      r[i+j*ldr] = complex <double> ( 0.0, 0.0 );
    }
  }
  for ( j = 0; j < NZ; j++ )
  {
    for ( i = 0; i < P; i++ )
    {
      z[i+j*ldz] = complex <double> ( 0.0, 0.0 );
    }
  }

  for ( i = 1; i <= P; i++ )
  {
    x[i-1] = complex <double> ( i, ( i % 2 ) );
  }
//
//  Use ZCHUD to form R, Z and RHO by adding X and Y a row at a time.
//  X is a row of the least squares matrix and Y the right hand side.
//
  seed = 123456789;

  for ( i = 1; i <= P; i++ )
  {
    row = c8vec_uniform_01 ( P, &seed );
    y[0] = zdotu ( P, row, 1, x, 1 );
    rho[0] = 0.0;
    zchud ( r, ldr, P, row, z, ldz, NZ, y, rho, c, s );
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

    info = ztrsl ( r, ldr, P, b, job );

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
//    TEST04 tests ZGBCO.
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
# define ML 1
# define MU 1

  complex <double> a[(2*ML+MU+1)*N];
  complex <double> a_save[N*N];
  int i;
  int i1;
  int i2;
  int ipvt[N];
  int j;
  int k;
  int lda;
  int m;
  double rcond;
  int seed;

  lda = 2 * ML + MU + 1;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For a complex general band storage matrix:\n";
  cout << "  ZGBCO factors the matrix and estimates the\n";
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
      a_save[i+j*N] = complex <double> ( 0.0, 0.0 );
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
      a[k-1+(j-1)*lda] = c8_uniform_01 ( &seed );
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
  rcond = zgbco ( a, lda, N, ML, MU, ipvt );

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
//    TEST05 tests ZGBFA and ZGBSL.
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
# define ML 1
# define MU 1

  complex <double> a[(2*ML+MU+1)*N];
  complex <double> a_save[N*N];
  complex <double> b[N];
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
  complex <double> *x;

  lda = 2 * ML + MU + 1;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For a complex general band storage matrix:\n";
  cout << "  ZGBFA factors the matrix;\n";
  cout << "  ZGBSL solves a factored linear system.\n";
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
      a_save[i+j*N] = complex <double> ( 0.0, 0.0 );
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
      a[k-1+(j-1)*lda] = c8_uniform_01 ( &seed );
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
  x = c8vec_uniform_01 ( N, &seed );

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
  info = zgbfa ( a, lda, N, ML, MU, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZGBFA returned INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  job = 0;
  zgbsl ( a, lda, N, ML, MU, ipvt, b, job );

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
//    TEST06 tests ZGBFA and ZGBDI.
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
# define ML 1
# define MU 1

  complex <double> a[(2*ML+MU+1)*N];
  complex <double> a_save[N*N];
  complex <double> det[2];
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
  cout << "  ZGBFA factors the matrix.\n";
  cout << "  ZGBDI computes the determinant.\n";
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
      a_save[i+j*N] = complex <double> ( 0.0, 0.0 );
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
      a[k-1+(j-1)*lda] = c8_uniform_01 ( &seed );
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
  info = zgbfa ( a, lda, N, ML, MU, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZGBFA returned INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  zgbdi ( a, lda, N, ML, MU, ipvt, det );

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
//    TEST07 tests ZGECO.
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

  complex <double> *a;
  int i;
  int ipvt[N];
  int j;
  int lda;
  double rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  For a complex general storage matrix:\n";
  cout << "  ZGECO factors the matrix and estimates the\n";
  cout << "  reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c8mat_uniform_01 ( N, N, &seed );

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
  rcond = zgeco ( a, lda, N, ipvt );

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
//    TEST08 tests ZGEFA and ZGESL.
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

  complex <double> *a;
  complex <double> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int job;
  int lda;
  int seed;
  complex <double> *x;

  lda = N;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For a complex general storage matrix:\n";
  cout << "  ZGEFA factors the matrix.\n";
  cout << "  ZGESL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c8mat_uniform_01 ( N, N, &seed );

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
  x = c8vec_uniform_01 ( N, &seed );

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
  info = zgefa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZGEFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  job = 0;
  zgesl ( a, lda, N, ipvt, b, job );

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
//    TEST09 tests ZGEFA and ZGEDI.
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

  complex <double> *a;
  complex <double> a_save[N*N];
  complex <double> c[N*N];
  complex <double> det[2];
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
  cout << "  ZGEFA factors the matrix.\n";
  cout << "  ZGEDI computes the determinant or inverse.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c8mat_uniform_01 ( N, N, &seed );

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
  info = zgefa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZGEFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  zgedi ( a, lda, N, ipvt, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inverse.
//
  job = 01;
  zgedi ( a, lda, N, ipvt, det, job );

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
//    TEST10 tests ZGTSL.
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

  complex <double> b[N];
  complex <double> *c;
  complex <double> d[N];
  complex <double> *e;
  int i;
  int info;
  int seed;
  complex <double> x[N];

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For a complex tridiagonal matrix:\n";
  cout << "  ZGTSL solves a linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  seed = 123456789;

  c = c8vec_uniform_01 ( N, &seed );
  c[0] = complex <double> ( 0.0, 0.0 );

  e = c8vec_uniform_01 ( N, &seed );
  e[N-1] = complex <double> ( 0.0, 0.0 );

  for ( i = 1; i <= N; i++ )
  {
    d[i-1] = 0.0;
    if ( i < N )
    {
      d[i-1] = d[i-1] - complex <double> ( 2.0 ) * e[i-1];
    }
    if ( 1 < i )
    {
      d[i-1] = d[i-1] - complex <double> ( 2.0 ) * c[i-1];
    }
  }
//
//  Set the desired solution
//
  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = complex <double> ( i, 10 * i );
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
  info = zgtsl ( N, c, d, e, b );

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
//    TEST11 tests ZHICO.
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

  complex <double> a[N*N];
  int i;
  int ipvt[N];
  int j;
  int lda;
  double rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  For a complex Hermitian matrix:\n";
  cout << "  ZHICO factors the matrix and estimates\n";
  cout << "  the reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = complex <double> ( r8_uniform_01 ( &seed ), 0.0 );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c8_uniform_01 ( &seed );
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
  rcond = zhico ( a, lda, N, ipvt );

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
//    TEST12 tests ZHIFA and ZHISL.
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

  complex <double> a[N*N];
  complex <double> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int lda;
  int seed;
  complex <double> *x;

  lda = N;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  For a complex Hermitian matrix:\n";
  cout << "  ZHIFA factors the matrix.\n";
  cout << "  ZHISL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = complex <double> ( r8_uniform_01 ( &seed ), 0.0 );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c8_uniform_01 ( &seed );
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
  x = c8vec_uniform_01 ( N, &seed );

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
  info = zhifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZHIFA returned an error flag INFO = " << info << "\n";
    delete [] x;
    return;
  }
//
//  Solve the system.
//
  zhisl ( a, lda, N, ipvt, b );

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
//    TEST13 tests ZHIFA and ZHIDI.
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

  complex <double> a[N*N];
  complex <double> a_save[N*N];
  complex <double> c[N*N];
  double det[2];
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
  cout << "  ZHIFA factors the matrix.\n";
  cout << "  ZHIDI computes the determinant, inverse,\n";
  cout << "  or inertia.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = complex <double> ( r8_uniform_01 ( &seed ), 0.0 );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c8_uniform_01 ( &seed );
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
  info = zhifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZHIFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  zhidi ( a, lda, N, ipvt, det, inert, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inertia.
//
  job = 100;
  zhidi ( a, lda, N, ipvt, det, inert, job );

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
  zhidi ( a, lda, N, ipvt, det, inert, job );
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
//    TEST14 tests ZHPCO.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> a_save[N*N];
  int i;
  int ipvt[N];
  int j;
  int k;
  double rcond;
  int seed;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  For a complex Hermitian matrix\n";
  cout << "  using packed storage,\n";
  cout << "  ZHPCO factors the matrix and estimates\n";
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
      a[k] = c8_uniform_01 ( &seed );
      a_save[i+j*N] = a[k];
      a_save[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
    a[k] = complex <double> ( r8_uniform_01 ( &seed ), 0.0 );
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
  rcond = zhpco ( a, N, ipvt );

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
//    TEST15 tests ZHPFA and ZHPSL.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> a_save[N*N];
  complex <double> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;
  int seed;
  complex <double> *x;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  For a complex Hermitian matrix,\n";
  cout << "  using packed storage,\n";
  cout << "  ZHPFA factors the matrix.\n";
  cout << "  ZHPSL solves a linear system.\n";
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
      a[k] = c8_uniform_01 ( &seed );
      a_save[i+j*N] = a[k];
      a_save[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
    a[k] = complex <double> ( r8_uniform_01 ( &seed ), 0.0 );
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
  x = c8vec_uniform_01 ( N, &seed );

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
  info = zhpfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZHPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  zhpsl ( a, N, ipvt, b );

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
//    TEST16 tests ZHPFA and ZHPDI.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> a_save[N*N];
  complex <double> b[N*N];
  complex <double> c[N*N];
  double det[2];
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
  cout << "  ZHPFA factors the matrix.\n";
  cout << "  ZHPDI computes the determinant, inverse,\n";
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
      a[k] = c8_uniform_01 ( &seed );
      a_save[i+j*N] = a[k];
      a_save[j+i*N] = conj ( a[k] );
      k = k + 1;
    }
    a[k] = complex <double> ( r8_uniform_01 ( &seed ), 0.0 );
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
  info = zhpfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZHPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  zhpdi ( a, N, ipvt, det, inert, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inertia.
//
  job = 100;
  zhpdi ( a, N, ipvt, det, inert, job );

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
  zhpdi ( a, N, ipvt, det, inert, job );
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
//    TEST17 tests ZPBCO.
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
# define M 1

  complex <double> a[(M+1)*N];
  int i;
  int info;
  int lda;
  double rcond;

  lda = M + 1;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  For a complex positive definite hermitian band matrix,\n";
  cout << "  ZPBCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  a[0+0*lda] = complex <double> ( 0.0000,  0.0000 );
  a[0+1*lda] = complex <double> ( 2.1341, -0.2147 );
  a[0+2*lda] = complex <double> ( 2.0905,  1.1505 );

  a[1+0*lda] = complex <double> ( 4.5281,  0.0000 );
  a[1+1*lda] = complex <double> ( 5.0371,  0.0000 );
  a[1+2*lda] = complex <double> ( 4.7638,  0.0000 );
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = zpbco ( a, lda, N, M, &info );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZPBCO returned INFO = " << info << "\n";
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
//    TEST18 tests ZPBDI.
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
# define M 1

  complex <double> a[(M+1)*N];
  double det[2];
  int i;
  int info;
  int lda;

  lda = M + 1;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  For a complex positive definite hermitian band matrix,\n";
  cout << "  ZPBDI computes the determinant as\n";
  cout << "    det = MANTISSA * 10^EXPONENT\n";
  cout << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  a[0+0*lda] = complex <double> ( 0.0000,  0.0000 );
  a[0+1*lda] = complex <double> ( 2.1341, -0.2147 );
  a[0+2*lda] = complex <double> ( 2.0905,  1.1505 );

  a[1+0*lda] = complex <double> ( 4.5281,  0.0000 );
  a[1+1*lda] = complex <double> ( 5.0371,  0.0000 );
  a[1+2*lda] = complex <double> ( 4.7638,  0.0000 );

  info = zpbfa ( a, lda, N, M );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error!  ZPBFA returns INFO = " << info << "\n";
    return;
  }

  zpbdi ( a, lda, N, M, det );

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
//    TEST19 tests ZPBFA and ZPBSL.
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
# define M 1

  complex <double> a[(M+1)*N];
  complex <double> b[N];
  int i;
  int info;
  int lda;

  lda = M + 1;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  For a complex positive definite hermitian band matrix,\n";
  cout << "  ZPBFA computes the LU factors.\n";
  cout << "  ZPBSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  a[0+0*lda] = complex <double> ( 0.0000,  0.0000 );
  a[0+1*lda] = complex <double> ( 2.1341, -0.2147 );
  a[0+2*lda] = complex <double> ( 2.0905,  1.1505 );

  a[1+0*lda] = complex <double> ( 4.5281,  0.0000 );
  a[1+1*lda] = complex <double> ( 5.0371,  0.0000 );
  a[1+2*lda] = complex <double> ( 4.7638,  0.0000 );
//
//  Set the right hand side.
//
  b[0] = complex <double> (  8.7963, -0.4294 );
  b[1] = complex <double> ( 18.4798,  3.6662 );
  b[2] = complex <double> ( 18.4724, -2.3010 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = zpbfa ( a, lda, N, M );

  if ( info != 0 )
  {
    cout << "  Error!  ZPBFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  zpbsl ( a, lda, N, M, b );
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
//    TEST20 tests ZPOCO.
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

  complex <double> a[N*N];
  int i;
  int info;
  int lda;
  double rcond;

  lda = N;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  ZPOCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*3] = complex <double> ( 2.5281,  0.0000 );
  a[1+0*3] = complex <double> ( 2.1341,  0.2147 );
  a[2+0*3] = complex <double> ( 2.4187, -0.2932 );

  a[0+1*3] = complex <double> ( 2.1341, -0.2147 );
  a[1+1*3] = complex <double> ( 3.0371,  0.0000 );
  a[2+1*3] = complex <double> ( 2.0905, -1.1505 );

  a[0+2*3] = complex <double> ( 2.4187,  0.2932 );
  a[1+2*3] = complex <double> ( 2.0905,  1.1505 );
  a[2+2*3] = complex <double> ( 2.7638,  0.0000 );
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition.\n";

  rcond = zpoco ( a, lda, N, &info );

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
//    TEST21 tests ZPOFA and ZPODI.
//
//  Discussion:
//
//    ZPOFA factors a positive definite symmetric matrix,
//    and ZPODI can compute the determinant or the inverse.
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

  complex <double> a[N*N];
  double det[2];
  int i;
  int info;
  int j;
  int job;
  int lda;

  lda = N;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  ZPOFA computes the LU factors,\n";
  cout << "  ZPODI computes the inverse or determinant.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*3] = complex <double> ( 2.5281,  0.0000 );
  a[1+0*3] = complex <double> ( 2.1341,  0.2147 );
  a[2+0*3] = complex <double> ( 2.4187, -0.2932 );

  a[0+1*3] = complex <double> ( 2.1341, -0.2147 );
  a[1+1*3] = complex <double> ( 3.0371,  0.0000 );
  a[2+1*3] = complex <double> ( 2.0905, -1.1505 );

  a[0+2*3] = complex <double> ( 2.4187,  0.2932 );
  a[1+2*3] = complex <double> ( 2.0905,  1.1505 );
  a[2+2*3] = complex <double> ( 2.7638,  0.0000 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = zpofa ( a, lda, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, ZPOFA returns INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant and inverse.
//
  cout << "\n";
  cout << "  Get the determinant and inverse.\n";

  job = 11;
  zpodi ( a, lda, N, det, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Determinant  = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  ZPODI produces only the 'upper half triangle' of the inverse,
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
//    TEST22 tests ZPOFA and ZPOSL.
//
//  Discussion:
//
//    ZPOFA factors a Hermitian positive definite matrix,
//    and ZPOSL can solve a factored linear system.
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

  complex <double> a[N*N];
  complex <double> b[N];
  int i;
  int info;
  int j;
  int lda;
  complex <double> x[N];

  lda = N;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  For a complex Hermitian positive definite matrix,\n";
  cout << "  ZPOFA computes the LU factors.\n";
  cout << "  ZPOSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0+0*3] = complex <double> ( 2.5281,  0.0000 );
  a[1+0*3] = complex <double> ( 2.1341,  0.2147 );
  a[2+0*3] = complex <double> ( 2.4187, -0.2932 );

  a[0+1*3] = complex <double> ( 2.1341, -0.2147 );
  a[1+1*3] = complex <double> ( 3.0371,  0.0000 );
  a[2+1*3] = complex <double> ( 2.0905, -1.1505 );

  a[0+2*3] = complex <double> ( 2.4187,  0.2932 );
  a[1+2*3] = complex <double> ( 2.0905,  1.1505 );
  a[2+2*3] = complex <double> ( 2.7638,  0.0000 );
//
//  Set the right hand side.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( 2 * i + 1, 2 * i + 2 );
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

  info = zpofa ( a, lda, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, ZPOFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  zposl ( a, lda, N, b );
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
//    TEST23 tests ZPPCO.
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

  complex <double> a[(N*(N+1))/2];
  int info;
  double rcond;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite packed matrix,\n";
  cout << "  ZPPCO estimates the reciprocal condition number.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0] = complex <double> ( 2.5281,  0.0000 );

  a[1] = complex <double> ( 2.1341, -0.2147 );
  a[2] = complex <double> ( 3.0371,  0.0000 );

  a[3] = complex <double> ( 2.4187,  0.2932 );
  a[4] = complex <double> ( 2.0905,  1.1505 );
  a[5] = complex <double> ( 2.7638,  0.0000 );
//
//  Estimate the condition.
//
  cout << "\n";
  cout << "  Estimate the condition number.\n";

  rcond = zppco ( a, N, &info );

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
//    TEST24 tests ZPPFA and ZPPDI.
//
//  Discussion:
//
//    ZPPFA factors a Hermitian positive definite packed matrix,
//    and ZPPDI can compute the determinant or the inverse.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> b[N*N];
  double det[2];
  int i;
  int info;
  int j;
  int job;
  int k;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite packed matrix,\n";
  cout << "  ZPPFA factors the matrix.\n";
  cout << "  ZPPDI computes the inverse or determinant.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0] = complex <double> ( 2.5281,  0.0000 );

  a[1] = complex <double> ( 2.1341, -0.2147 );
  a[2] = complex <double> ( 3.0371,  0.0000 );

  a[3] = complex <double> ( 2.4187,  0.2932 );
  a[4] = complex <double> ( 2.0905,  1.1505 );
  a[5] = complex <double> ( 2.7638,  0.0000 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = zppfa ( a, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, ZPPFA returns INFO = " << info << "\n";
    return;
  }
//
//  Invert the matrix.
//
  cout << "\n";
  cout << "  Get the determinant and inverse.\n";

  job = 11;
  zppdi ( a, N, det, job );
//
//  Print the results.
//
  cout << "\n";
  cout << "  Determinant  = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  ZPPDI produces only the 'upper half triangle' of the inverse,
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
//    TEST25 tests ZPPFA and ZPPSL.
//
//  Discussion:
//
//    ZPOFA factors a Hermitian positive definite packed matrix,
//    and ZPOSL can solve a factored linear system.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> b[N];
  int i;
  int info;
  int j;
  int k;
  complex <double> x[N];

  cout << "\n";
  cout << "TEST25\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite packed matrix,\n";
  cout << "  ZPPFA factors the matrix.\n";
  cout << "  ZPPSL solves a factored linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  a[0] = complex <double> ( 2.5281,  0.0000 );

  a[1] = complex <double> ( 2.1341, -0.2147 );
  a[2] = complex <double> ( 3.0371,  0.0000 );

  a[3] = complex <double> ( 2.4187,  0.2932 );
  a[4] = complex <double> ( 2.0905,  1.1505 );
  a[5] = complex <double> ( 2.7638,  0.0000 );

  b[0] = complex <double> ( 20.12350, 28.92670 );
  b[1] = complex <double> ( 14.36550, 34.92680 );
  b[2] = complex <double> ( 27.69760, 26.03750 );
//
//  Factor the matrix.
//
  cout << "\n";
  cout << "  Factor the matrix.\n";

  info = zppfa ( a, N );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Error, ZPPFA returns INFO = " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  cout << "\n";
  cout << "  Solve the linear system.\n";

  zppsl ( a, N, b );
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
//    TEST26 tests ZPTSL.
//
//  Discussion:
//
//    ZPTSL factors and solves a Hermitian positive definite
//    tridiagonal system.
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

  complex <double> b[N];
  complex <double> d[N];
  complex <double> e[N];
  int i;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  For a complex Hermitian\n";
  cout << "  positive definite tridiagonal matrix,\n";
  cout << "  ZPTSL factors and solves a linear system.\n";
  cout << "  The matrix size is N = " << N << "\n";
//
//  Set the value of the superdiagonal and diagonal.
//
  e[0] = complex <double> ( 2.1341, -0.2147 );
  e[1] = complex <double> ( 2.0905,  1.1505 );
  e[2] = complex <double> ( 0.0000,  0.0000 );

  d[0] = complex <double> ( 4.5281,  0.0000 );
  d[1] = complex <double> ( 5.0371,  0.0000 );
  d[2] = complex <double> ( 4.7638,  0.0000 );
//
//  Set the right hand side.
//
  b[0] = complex <double> (  8.7963, -0.4294 );
  b[1] = complex <double> ( 18.4798,  3.6662 );
  b[2] = complex <double> ( 18.4724, -2.3010 );
//
//  Factor and solve the system.
//
  cout << "\n";
  cout << "  Factor the matrix and solve the system.\n";

  zptsl ( N, d, e, b );
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
//    TEST27 tests ZQRDC and ZQRSL.
//
//  Discussion:
//
//    ZQRDC and ZQRSL compute the QR factorization, and use it
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

  complex <double> *a;
  complex <double> b[N*P];
  int i;
  int info;
  int ipvt[P];
  int j;
  int job;
  int k;
  int lda;
  complex <double> q[N*N];
  complex <double> qraux[P];
  complex <double> qty[N];
  complex <double> qy[N];
  complex <double> r[N*P];
  complex <double> rsd[N];
  int seed;
  complex <double> xb[N];
  complex <double> y[N];

  lda = N;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  For a complex general matrix,\n";
  cout << "  ZQRDC computes the QR decomposition of a\n";
  cout << "  matrix, but does not return Q and R explicitly.\n";
  cout << "\n";
  cout << "  Show how Q and R can be recovered using ZQRSL.\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  a = c8mat_uniform_01 ( N, P, &seed );

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

  zqrdc ( a, lda, N, P, qraux, ipvt, job );
//
//  Print out what ZQRDC has stored in A...
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
        r[i+j*lda] = complex <double> ( 0.0, 0.0 );
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
//  Call ZQRSL to extract the information about the Q matrix.
//  We do this, essentially, by asking ZQRSL to tell us the
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
      y[i] = complex <double> ( 0.0, 0.0 );
    }
    y[j] = complex <double> ( 1.0, 0.0 );
//
//  Ask ZQRSL to tell us what Q*Y is.
//
    info = zqrsl ( a, lda, N, P, qraux, y, qy, qty, b, rsd, xb, job );

    if ( info != 0 )
    {
      cout << "  Error!  ZQRSL returns INFO = " << info << "\n";
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
      b[i+j*N] = complex <double> ( 0.0, 0.0 );
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
//    TEST28 tests ZSICO.
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

  complex <double> a[N*N];
  int i;
  int ipvt[N];
  int j;
  int lda;
  double rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST28\n";
  cout << "  For a complex symmetric matrix:\n";
  cout << "  ZSICO factors the matrix and estimates\n";
  cout << "  the reciprocal condition number.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = c8_uniform_01 ( &seed );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c8_uniform_01 ( &seed );
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
  rcond = zsico ( a, lda, N, ipvt );

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
//    TEST29 tests ZSIFA and ZSISL.
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

  complex <double> a[N*N];
  complex <double> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int lda;
  int seed;
  complex <double> x[N];

  lda = N;

  cout << "\n";
  cout << "TEST29\n";
  cout << "  For a complex symmetric matrix:\n";
  cout << "  ZSIFA factors the matrix.\n";
  cout << "  ZSISL solves a linear system.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = c8_uniform_01 ( &seed );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c8_uniform_01 ( &seed );
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
    x[i] = c8_uniform_01 ( &seed );
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
  info = zsifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZSIFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  zsisl ( a, lda, N, ipvt, b );

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
//    TEST30 tests ZSIFA and ZSIDI.
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

  complex <double> a[N*N];
  complex <double> a_save[N*N];
  complex <double> c[N*N];
  complex <double> det[2];
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
  cout << "  ZSIFA factors the matrix.\n";
  cout << "  ZSIDI computes the determinant or inverse.\n";
  cout << "\n";
  cout << "  The matrix order is N = " << N << "\n";
//
//  Set the values of the matrix A.
//
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i+i*lda] = c8_uniform_01 ( &seed );
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = c8_uniform_01 ( &seed );
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
  info = zsifa ( a, lda, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZSIFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  zsidi ( a, lda, N, ipvt, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inverse.
//
  job = 1;
  zsidi ( a, lda, N, ipvt, det, job );
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
//    TEST31 tests ZSPCO.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> a_save[N*N];
  int i;
  int ipvt[N];
  int j;
  int k;
  double rcond;
  int seed;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  For a complex symmetric matrix\n";
  cout << "  in packed storage,\n";
  cout << "  ZSPCO factors the matrix and estimates\n";
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
      a[k] = c8_uniform_01 ( &seed );
      k = k + 1;
    }
    a[k] = c8_uniform_01 ( &seed );
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
  rcond = zspco ( a, N, ipvt );

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
//    TEST32 tests ZSPFA and ZSPSL.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> a_save[N*N];
  complex <double> b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;
  int seed;
  complex <double> *x;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  For a complex symmetric matrix\n";
  cout << "  in packed storage,\n";
  cout << "  ZSPFA factors the matrix.\n";
  cout << "  ZSPSL solves a linear system.\n";
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
      a[k] = c8_uniform_01 ( &seed );
      k = k + 1;
    }
    a[k] = c8_uniform_01 ( &seed );
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
  x = c8vec_uniform_01 ( N, &seed );

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
  info = zspfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZSPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  zspsl ( a, N, ipvt, b );

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
//    TEST33 tests ZSPFA and ZSPDI.
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

  complex <double> a[(N*(N+1))/2];
  complex <double> a_save[N*N];
  complex <double> b_save[N*N];
  complex <double> c[N*N];
  complex <double> det[2];
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
  cout << "  ZSPFA factors the matrix.\n";
  cout << "  ZSPDI computes the determinant or inverse.\n";
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
      a[k] = c8_uniform_01 ( &seed );
      k = k + 1;
    }
    a[k] = c8_uniform_01 ( &seed );
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
  info = zspfa ( a, N, ipvt );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  ZSPFA returned an error flag INFO = " << info << "\n";
    return;
  }
//
//  Get the determinant.
//
  job = 10;
  zspdi ( a, N, ipvt, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << det[1] << "\n";
//
//  Get the inverse.
//
  job = 1;
  zspdi ( a, N, ipvt, det, job );
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
//**********************************************************************

void test34 ( )

//**********************************************************************
//
//  Purpose:
//
//    TEST34 tests ZSVDC.
//
//  Discussion:
//
//    ZSVDC computes the singular value decomposition:
//
//      A = U * S * conjg-transpose ( V )
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
# define M 4
# define N 3

  complex <double> *a;
//
//  E should be dimensioned at least the maximum of M+1 and N.
//
  complex <double> e[M+N];
  int i;
  int info;
  int j;
  int lda;
  int ldu;
  int ldv;
  int job;
  int k;
//
//  S should be dimensioned at least the maximum of M+1 and N.
//
  complex <double> s[M+N];
  int seed;
  complex <double> sigma[M*N];
  complex <double> u[M*M];
  complex <double> v[N*N];

  cout << "\n";
  cout << "TEST34\n";
  cout << "  For an MxN matrix A in complex general storage,\n";
  cout << "  ZSVDC computes the singular value decomposition:\n";
  cout << "    A = U * S * V^H\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
//
//  Set A.
//
  seed = 123456789;
  lda = M;

  a = c8mat_uniform_01 ( M, N, &seed );

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

  info = zsvdc ( a, lda, M, N, s, e, u, ldu, v, ldv, job );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning:\n";
    cout << "  ZSVDC returned nonzero INFO = " << info << "\n";
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
      a[i+j*lda] = complex <double> ( 0.0, 0.0 );
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
//    TEST345 tests ZSVDC.
//
//  Discussion:
//
//    ZSVDC computes the singular value decomposition:
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

  complex <double> *a;
//
//  E must be dimensioned at least the maximum of M+1 and N.
//
  complex <double> e[M+N];
  complex <double> I;
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
  complex <double> s[M+N];
  int seed;
  complex <double> sigma[M*N];
  complex <double> u[M*M];
  complex <double> v[N*N];

  cout << "\n";
  cout << "TEST345\n";
  cout << "  For an MxN matrix A in double complex general storage,\n";
  cout << "  ZSVDC computes the singular value decomposition:\n";
  cout << "    A = U * S * V^H\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
//
//  Set A.
//
  I = complex <double> ( 0.0, 1.0 );

  lda = M;
  a = new complex <double>[M*N];

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

  info = zsvdc ( a, lda, M, N, s, e, u, ldu, v, ldv, job );

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
      a[i+j*lda] = complex <double> ( 0.0, 0.0 );
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
//    TEST35 tests ZTRCO.
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

  complex <double> a[N*N];
  int i;
  int j;
  int job;
  int lda;
  double rcond;
  int seed;

  lda = N;

  cout << "\n";
  cout << "TEST35\n";
  cout << "  For a complex triangular matrix,\n";
  cout << "  ZTRCO estimates the condition.\n";
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
      a[i+j*lda] = c8_uniform_01 ( &seed );
    }
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = complex <double> ( 0.0, 0.0 );
    }
  }
//
//  Get the condition of the lower triangular matrix.
//
  job = 0;
  rcond = ztrco ( a, lda, N, job );

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
//    TEST36 tests ZTRDI.
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

  complex <double> a[N*N];
  complex <double> a_save[N*N];
  complex <double> c[N*N];
  complex <double> det[2];
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
  cout << "  ZTRDI computes the determinant or inverse.\n";
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
      a[i+j*lda] = c8_uniform_01 ( &seed );
    }
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = complex <double> ( 0.0, 0.0 );
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
  info = ztrdi ( a, lda, N, det, job );

  cout << "\n";
  cout << "  Determinant = " << det[0]
       << " * 10^ " << real ( det[1] ) << "\n";
//
//  Get the inverse of the lower triangular matrix.
//
  job = 10;
  info = ztrdi ( a, lda, N, det, job );

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
//    TEST37 tests ZTRSL.
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

  complex <double> a[N*N];
  complex <double> b[N];
  int i;
  int info;
  int j;
  int job;
  int k;
  int lda;
  int seed;
  complex <double> x[N];

  lda = N;

  cout << "\n";
  cout << "TEST37\n";
  cout << "  For a complex triangular matrix,\n";
  cout << "  ZTRSL solves a linear system.\n";
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
      a[i+j*lda] = c8_uniform_01 ( &seed );
    }
    for ( j = i+1; j < N; j++ )
    {
      a[i+j*lda] = complex <double> ( 0.0, 0.0 );
    }
  }
//
//  Set the desired solution
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = complex <double> ( i + 1, 10 * ( i + 1 ) );
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
  info = ztrsl ( a, lda, N, b, job );

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

complex <double> c8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8_UNIFORM_01 returns a unit double complex pseudorandom number.
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
//    23 June 2009
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
//    Output, complex <double> C8_UNIFORM_01, a pseudorandom complex value.
//
{
  double r;
  int k;
  double pi = 3.14159265358979300;
  double theta;
  complex <double> value;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = sqrt ( ( ( double ) ( *seed ) * 4.656612875E-10 ) );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

  value = complex <double> ( r * cos ( theta ), r * sin ( theta ) );

  return value;
}
//****************************************************************************80

complex <double> *c8mat_uniform_01 ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_UNIFORM_01 returns a unit complex pseudorandom matrix.
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
//    23 June 2009
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
//    Output, complex C8MAT_UNIFORM_01[M*N], the pseudorandom complex matrix.
//
{
  complex <double> *c;
  int i;
  int j;
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;

  c = new complex <double>[m*n];

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

      c[i+j*m] = r * complex <double> ( cos ( theta ), sin ( theta ) );
    }
  }

  return c;
}
//****************************************************************************80

complex <double> *c8vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNIFORM_01 returns a unit complex pseudorandom vector.
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
//    23 June 2009
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
//    Output, complex <double> C8VEC_UNIFORM_01[N], the pseudorandom
//    complex vector.
//
{
  complex <double> *c;
  int i;
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;

  c = new complex <double> [n];

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

    c[i] = r * complex <double> ( cos ( theta ), sin ( theta ) );
  }

  return c;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a real pseudorandom number.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
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
//    23 June 2009
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
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

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
  r = ( double ) ( *seed ) * 4.656612875E-10;

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
