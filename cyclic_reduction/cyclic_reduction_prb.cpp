# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "cyclic_reduction.hpp"

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
//    CYCLIC_REDUCTION_PRB calls the CYCLIC_REDUCTION tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CYCLIC_REDUCTION_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CYCLIC_REDUCTION library.\n";

  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CYCLIC_REDUCTION_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R83_CR_FA, R83_CR_SLS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_cr;
  double *b;
  bool debug = true;
  int i;
  int j;
  int n = 5;
  int nb = 2;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R83_CR_FA factors a real tridiagonal matrix;\n";
  cout << "  R83_CR_SLS solves 1 or more systems.\n";
  cout << "\n";
  cout << "  Matrix order N = " << n << "\n";
  cout << "  Demonstrate multiple system solution method.\n";
//
//  Set the matrix.
//
  a = new double[3*n];

  a[0+0*3] = 0.0;
  for ( j = 1; j < n; j++ )
  {
    a[0+j*3] = - 1.0;
  }
  for ( j = 0; j < n; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < n - 1; j++ )
  {
    a[2+j*3] = - 1.0;
  }
  a[2+(n-1)*3] = 0.0;

  if ( debug )
  {
    r83_print ( n, a, "  Input matrix:" );
  }
//
//  Factor the matrix once.
//
  a_cr = r83_cr_fa ( n, a );

  if ( debug )
  {
    r83_print ( 2 * n + 1, a_cr, "  Cyclic reduction factor information:" );
  }
//
//  Solve 2 systems simultaneously.
//
  b = new double[n*nb];

  for ( i = 0; i < n - 1; i++ )
  {
    b[i+0*n] = 0.0;
  }
  b[n-1+0*n] = ( double ) ( n + 1 );

  b[0+1*n] = 1.0;
  for ( i = 1; i < n - 1; i++ )
  {
    b[i+1*n] = 0.0;
  }
  b[n-1+1*n] = 1.0;
//
//  Solve the linear systems.
//
  x = r83_cr_sls ( n, a_cr, nb, b );

  r8mat_print_some ( n, nb, x, 1, 1, 10, nb, "  Solutions:" );

  delete [] a;
  delete [] a_cr;
  delete [] b;
  delete [] x;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R83_CR_FA, R83_CR_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_cr;
  double *b;
  bool debug = false;
  int i;
  int j;
  int n = 10;
  double *x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a real tridiagonal matrix,\n";
  cout << "  using CYCLIC REDUCTION,\n";
  cout << "  R83_CR_FA factors;\n";
  cout << "  R83_CR_SL solves a system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << n << "\n";
  cout << "  The matrix is NOT symmetric.\n";
//
//  Set the matrix values.
//
  a = new double[3*n];

  a[0+0*3] = 0.0;
  for ( j = 2; j <= n; j++ )
  {
    a[0+(j-1)*3] = ( double ) ( j );
  }
  for ( j = 1; j <= n; j++ )
  {
    a[1+(j-1)*3] = 4.0 * ( double ) ( j );
  }
  for ( j = 1; j <= n - 1; j++ )
  {
    a[2+(j-1)*3] = ( double ) ( j );
  }
  a[2+(n-1)*3] = 0.0;

  if ( debug )
  {
    r83_print ( n, a, "  The matrix:" );
  }
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( n );
//
//  Compute the corresponding right hand side.
//
  b = r83_mxv_new ( n, a, x );

  if ( debug )
  {
    r8vec_print  ( n, b, "  The right hand side:" );
  }
  delete [] x;
//
//  Factor the matrix.
//
  a_cr = r83_cr_fa ( n, a );

  if ( debug )
  {
    r83_print ( 2 * n + 1, a_cr, "  The factor information:" );
  }
//
//  Solve the linear system.
//
  x = r83_cr_sl ( n, a_cr, b );

  r8vec_print  ( n, x, "  The solution:" );

  delete [] a;
  delete [] a_cr;
  delete [] b;
  delete [] x;

  return;
}
