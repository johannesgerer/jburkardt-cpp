# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "eispack.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test065 ( );
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

//****************************************************************************80

int main ()

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for EISPACK_PRB1.
//
//  Discussion:
//
//    EISPACK_PRB1 calls the EISPACK sample programs.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "EISPACK_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the EISPACK library.\n";
//
//  test01 ( );
//  test02 ( );
//  test03 ( );
//  test04 ( );
//  test05 ( );
//
  test06 ( );
  test065 ( );
  test07 ( );
//  test08 ( );
//  test09 ( );
//  test10 ( );

//  test11 ( );
//  test12 ( );
//  test13 ( );
//  test14 ( );
//  test15 ( );
//  test16 ( );

//
//  Terminate.
//
  cout << "\n";
  cout << "EISPACK_PRB1\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests RS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a[4*4] = {
    5.0, 4.0, 1.0, 1.0,
    4.0, 5.0, 1.0, 1.0,
    1.0, 1.0, 4.0, 2.0,
    1.0, 1.0, 2.0, 4.0 };
  double a2[4*4];
  int i;
  int ierr;
  int j;
  int k;
  int matz;
  int n = 4;
  double *r;
  double w[4];
  double x[4*4];
//
//  Save a copy of the matrix.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }
  cout << "\n";
  cout << "TEST06\n";
  cout << "  RS computes the eigenvalues and eigenvectors\n";
  cout << "  of a real symmetric matrix.\n";
  cout << "\n";
  cout << "  Matrix order = " << n << "\n";

  r8mat_print ( n, n, a, "  The matrix A:" );

  matz = 1;

  ierr = rs ( n, a, w, matz, x );

  if ( ierr != 0 )
  {
    cout << "\n";
    cout << "TEST06 - Warning!\n";
    cout << "  The error return flag IERR = " << ierr << "\n";
    return;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  if ( matz != 0 )
  {
    r8mat_print ( n, n, x, "  The eigenvector matrix:" );

    r = r8mat_mm_new ( n, n, n, a2, x );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i+j*n] = r[i+j*n] - w[j] * x[i+j*n];
      }
    }

    r8mat_print ( n, n, r, "  The residual (A-Lambda*I)*X:" );
  }

  delete [] r;

  return;
}
//****************************************************************************80

void test065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 tests RS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double a2[3*3];
  int i;
  int ierr;
  int j;
  int k;
  int matz;
  int n = 3;
  double *r;
  int seed;
  double t;
  double w[3];
  double x[3*3];

  cout << "\n";
  cout << "TEST065\n";
  cout << "  RS computes the eigenvalues and eigenvectors\n";
  cout << "  of a real symmetric matrix.\n";
  cout << "\n";
  cout << "  Matrix order = " << n << "\n";

  seed = 123456789;
  a = r8mat_uniform_01_new ( n, n, seed );

  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      t = 0.5 * ( a[i+j*n] + a[j+i*n] );
      a[i+j*n] = t;
      a[j+i*n] = t;
    }
  }
//
//  Save a copy of the matrix.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }

  r8mat_print ( n, n, a, "  The matrix A:" );

  matz = 1;

  ierr = rs ( n, a, w, matz, x );

  if ( ierr != 0 )
  {
    cout << "\n";
    cout << "TEST065 - Warning!\n";
    cout << "  The error return flag IERR = " << ierr << "\n";
    return;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  if ( matz != 0 )
  {
    r8mat_print ( n, n, x, "  The eigenvector matrix:" );

    r = r8mat_mm_new ( n, n, n, a2, x );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i+j*n] = r[i+j*n] - w[j] * x[i+j*n];
      }
    }

    r8mat_print ( n, n, r, "  The residual (A-Lambda*I)*X:" );

    delete [] r;
  }

  delete [] a;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests RSB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a2;
  int i;
  int ierr;
  int j;
  int matz;
  int mb = 2;
  int n = 5;
  double *r;
  double *w;
  double *x;

  a = new double[n*mb];

  for ( j = 0; j < mb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }
  j = mb - 1;
  for ( i = 0; i < n; i++ )
  {
    a[i+j*n] = 2.0;
  }
  j = 0;
  for ( i = 1; i < n; i++ )
  {
    a[i+j*n] = -1.0;
  }
  a2 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a2[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a2[i+j*n] = - 1.0;
      }
      else
      {
        a2[i+j*n] = 0.0;
      }
    }
  }

  cout << "\n";
  cout << "TEST07\n";
  cout << "  RSB computes the eigenvalues and eigenvectors\n";
  cout << "  of a real symmetric band matrix.\n";
  cout << "\n";
  cout << "  Matrix order = " << n << "\n";

  r8mat_print ( n, n, a2, "  The matrix A:" );

  w = new double[n];
  x = new double[n*n];
  matz = 1;

  ierr = rsb ( n, mb, a, w, matz, x );

  if ( ierr != 0 )
  {
    cout << "\n";
    cout << "TEST07 - Warning!\n";
    cout << "  The error return flag IERR = " << ierr << "\n";
    return;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  if ( matz != 0 )
  {
    r8mat_print ( n, n, x, "  The eigenvector matrix X:" );

    r = r8mat_mm_new ( n, n, n, a2, x );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i+j*n] = r[i+j*n] - x[i+j*n] * w[j];
      }
    }
    r8mat_print ( n, n, r, "  The residual (A-Lambda*I)*X:" );

    delete [] r;
  }

  delete [] a;
  delete [] a2;
  delete [] w;
  delete [] x;

  return;
}

