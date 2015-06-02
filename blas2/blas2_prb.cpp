# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "blas0.hpp"
# include "blas2.hpp"

int main ( );
void test01 ( );
void test05 ( );
void test06 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS2_PRB.
//
//  Discussion:
//
//    BLAS2_PRB tests the BLAS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BLAS2_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BLAS library.\n";

  test01 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BLAS2_PRB\n";
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
//    TEST01 tests DGEMV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double beta;
  int i;
  int incx;
  int incy;
  int j;
  int lda;
  int m;
  int n;
  char trans;
  double *x;
  double *y;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For a general matrix A,\n";
  cout << "  DGEMV computes y := alpha * A * x + beta * y\n";
  cout << "  or             y := alpha * A'' * x + beta * y.\n";
//
//  y = alpha * A * x + beta * y
//
  trans = 'N';
  m = 5;
  n = 4;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( trans, lda, m, n );
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  incx = 1;
  beta = 3.0;
  y = new double[m];
  for ( i = 0; i < m; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }
  incy = 1;

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( n, x, "  Vector X:" );
  r8vec_print ( m, y, "  Vector Y:" );

  dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy );

  r8vec_print ( m, y, "  Result Y = alpha * A  * x + beta * y" );

  delete [] a;
  delete [] x;
  delete [] y;
//
//  y = alpha * A' * x + beta * y
//
  trans = 'T';
  m = 5;
  n = 4;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( trans, lda, n, m );
  x = new double[m];
  for ( i = 0; i < m; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  incx = 1;
  beta = 3.0;
  y = new double[n];
  for ( i = 0; i < n; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }
  incy = 1;

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( m, x, "  Vector X:" );
  r8vec_print ( n, y, "  Vector Y:" );

  dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy );

  r8vec_print ( n, y, "  Result Y = alpha * A  * x + beta * y" );

  delete [] a;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests DGER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  int i;
  int incx;
  int incy;
  int lda;
  int m;
  int n;
  char trans;
  double *x;
  double *y;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For a general matrix A,\n";
  cout << "  DGER computes A := A + alpha * x * y'\n";

  m = 5;
  n = 4;
  alpha = 2.0;
  trans = 'N';
  lda = m;
  a = r8mat_test ( trans, lda, m, n );

  x = new double[m];
  for ( i = 0; i < m; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  incx = 1;

  y = new double[n];
  for ( i = 0; i < n; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }
  incy = 1;

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( m, x, "  Vector X:" );
  r8vec_print ( n, y, "  Vector Y:" );

  dger ( m, n, alpha, x, incx, y, incy, a, lda );

  r8mat_print ( m, n, a, "  Result A = A + alpha * x * y" );

  delete [] a;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests DTRMV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  char diag;
  int i;
  int incx;
  int j;
  int lda = 5;
  int m = 5;
  int n = 5;
  int test;
  char trans;
  char uplo;
  double *x;

  a = new double[lda*n];
  x = new double[n];

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a triangular matrix A,\n";
  cout << "  DTRMV computes y := A * x or y := A' * x\n";

  for ( test = 1; test <= 2; test++ )
  {
    uplo = 'U';

    if ( test == 1 )
    {
      trans = 'N';
    }
    else
    {
      trans = 'T';
    }

    diag = 'N';

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i <= j; i++ )
      {
        a[i+j*m] = ( double ) ( i + j + 2 );
      }
      for ( i = j + 1; i < m; i++ )
      {
        a[i+j*m] = 0.0;
      }
    }

    incx = 1;
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( double ) ( i + 1 );
    }

    dtrmv ( uplo, trans, diag, n, a, lda, x, incx );

    if ( trans == 'N' )
    {
      r8vec_print ( n, x, "  Result y = A * x" );
    }
    else
    {
      r8vec_print ( n, x, "  Result y = A' * x" );
    }
  }

  delete [] a;
  delete [] x;

  return;
}