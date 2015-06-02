# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "blas0.hpp"
# include "blas3.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS3_PRB.
//
//  Discussion:
//
//    BLAS3_PRB tests the BLAS3 library.
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
{
  timestamp ( );
  cout << "\n";
  cout << "BLAS3_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BLAS3 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BLAS3_PRB\n";
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
//    TEST01 tests DGEMM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double *b;
  double beta;
  double *c;
  int k;
  int lda;
  int ldb;
  int ldc;
  int m;
  int n;
  char transa;
  char transb;
  char transc;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  DGEMM carries out matrix multiplications\n";
  cout << "  for double precision real matrices.\n";
  cout << "\n";
  cout << "  1: C = alpha * A  * B  + beta * C;\n";
  cout << "  2: C = alpha * A' * B  + beta * C;\n";
  cout << "  3: C = alpha * A  * B' + beta * C;\n";
  cout << "  4: C = alpha * A' * B' + beta * C;\n";
  cout << "\n";
  cout << "  We carry out all four calculations, but in each case,\n";
  cout << "  we choose our input matrices so that we get the same result.\n";
//
//  C = alpha * A * B + beta * C.
//
  transa = 'N';
  transb = 'N';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( transa, lda, m, k );
  ldb = k;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A * B + beta * C:" );

  delete [] a;
  delete [] b;
  delete [] c;
//
//  C = alpha * A' * B + beta * C.
//
  transa = 'T';
  transb = 'N';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = k;
  a = r8mat_test ( transa, lda, m, k );
  ldb = k;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb,  m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A' * B + beta * C:" );

  delete [] a;
  delete [] b;
  delete [] c;
//
//  C = alpha * A * B' + beta * C.
//
  transa = 'N';
  transb = 'T';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( transa, lda, m, k );
  ldb = n;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A * B' + beta * C:" );

  delete [] a;
  delete [] b;
  delete [] c;
//
//  C = alpha * A' * B' + beta * C.
//
  transa = 'T';
  transb = 'T';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = k;
  a = r8mat_test ( transa, lda, m, k );
  ldb = n;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb,  m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A' * B' + beta * C:" );

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests DTRMM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double *b;
  char diag;
  int i;
  int j;
  int lda;
  int ldb;
  int m;
  int n;
  char side;
  char transa;
  char transb;
  char uplo;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  DTRMM multiplies a triangular matrix A and a\n";
  cout << "  rectangular matrix B\n";
  cout << "\n";
  cout << "  1: B = alpha * A  * B;\n";
  cout << "  2: B = alpha * A' * B;\n";
//
//  B = alpha * A * B.
//
  side = 'L';
  uplo = 'U';
  transa = 'N';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = new double[lda*m];
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrmm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  B = alpha * A * B:" );

  delete [] a;
  delete [] b;
//
//  B = alpha * A' * B.
//
  side = 'L';
  uplo = 'U';
  transa = 'T';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = new double[lda*m];
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrmm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  B = alpha * A' * B:" );

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests DTRSM.
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
{
  double *a;
  double alpha;
  double *b;
  char diag;
  int i;
  int j;
  int lda;
  int ldb;
  int m;
  int n;
  char side;
  char transa;
  char transb;
  char uplo;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  DTRSM solves a linear system involving a triangular\n";
  cout << "  matrix A and a rectangular matrix B.\n";
  cout << "\n";
  cout << "  1: Solve A  * X  = alpha * B;\n";
  cout << "  2: Solve A' * X  = alpha * B;\n";
  cout << "  3: Solve X  * A  = alpha * B;\n";
  cout << "  4: Solve X  * A' = alpha * B;\n";
//
//  Solve A * X = alpha * B.
//
  side = 'L';
  uplo = 'U';
  transa = 'N';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = new double[lda*m];

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = inv ( A ) * alpha * B:" );

  delete [] a;
  delete [] b;
//
//  Solve A' * X = alpha * B.
//
  side = 'L';
  uplo = 'U';
  transa = 'T';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = new double[lda*m];
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = inv ( A' ) * alpha * B:" );

  delete [] a;
  delete [] b;
//
//  Solve X * A = alpha * B.
//
  side = 'R';
  uplo = 'U';
  transa = 'N';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = n;
  ldb = m;

  a = new double[lda*n];
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = alpha * B * inv ( A ):" );

  delete [] a;
  delete [] b;
//
//  Solve X * A'' = alpha * B.
//
  side = 'R';
  uplo = 'U';
  transa = 'T';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = n;
  ldb = m;

  a = new double[lda*n];
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = alpha * B * inv ( A' ):" );

  delete [] a;
  delete [] b;

  return;
}

