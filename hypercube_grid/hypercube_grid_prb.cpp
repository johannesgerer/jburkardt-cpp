# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hypercube_grid.hpp"

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
//    MAIN is the main program for HYPERCUBE_GRID_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HYPERCUBE_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the HYPERCUBE_GRID library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HYPERCUBE_GRID_PRB:\n";
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
//    TEST01 tests HYPERCUBE_GRID on a two dimensional example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 2

  double a[M] = { 0.0, 0.0 };
  double b[M] = { 1.0, 10.0 };
  int c[M] = { 2, 4 };
  int i;
  int m = M;
  int n;
  int ns[M] = { 4, 5 };
  double *x;

  n = i4vec_product ( m, ns );

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Create a grid using HYPERCUBE_GRID.\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Number of grid points N = " << n << "\n";
  cout << "\n";
  cout << "     I    NS     C      A         B\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout 
      << setw(6) << i << "  "
      << setw(6) << ns[i] << "  "
      << setw(6) << c[i] << "  "
      << setw(8) << a[i] << "  "
      << setw(8) << b[i] << "\n";
  }

  x = hypercube_grid ( m, n, ns, a, b, c );
  r8mat_transpose_print ( m, n, x, "  Grid points:" );
  delete [] x;

  return;
# undef M
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HYPERCUBE_GRID on a five dimensional example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 5

  double a[M] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  double b[M] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  int c[M] = { 1, 2, 3, 4, 5 };
  int i;
  int m = M;
  int n;
  int ns[M] = { 2, 2, 2, 2, 2 };
  double *x;

  n = i4vec_product ( m, ns );

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Create a grid using HYPERCUBE_GRID.\n";
  cout << "  Use a two point grid in each dimension.\n";
  cout << "  Use a different centering option in each dimension.\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Number of grid points N = " << n << "\n";
  cout << "\n";
  cout << "     I    NS     C      A         B\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout 
      << setw(6) << i << "  "
      << setw(6) << ns[i] << "  "
      << setw(6) << c[i] << "  "
      << setw(8) << a[i] << "  "
      << setw(8) << b[i] << "\n";
  }

  x = hypercube_grid ( m, n, ns, a, b, c );
  r8mat_transpose_print ( m, n, x, "  Grid points:" );
  delete [] x;

  return;
# undef M
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests HYPERCUBE_GRID on a three dimensional example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 3

  double a[M] = { -1.0, -1.0, -1.0 };
  double b[M] = { +1.0, +1.0, +1.0 };
  int c[M] = { 1, 1, 1 };
  int i;
  int m = M;
  int n;
  int ns[M] = { 3, 3, 3 };
  double *x;

  n = i4vec_product ( m, ns );

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Create a grid using HYPERCUBE_GRID.\n";
  cout << "  Use the same parameters in every dimension.\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Number of grid points N = " << n << "\n";
  cout << "\n";
  cout << "     I    NS     C      A         B\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout 
      << setw(6) << i << "  "
      << setw(6) << ns[i] << "  "
      << setw(6) << c[i] << "  "
      << setw(8) << a[i] << "  "
      << setw(8) << b[i] << "\n";
  }

  x = hypercube_grid ( m, n, ns, a, b, c );
  r8mat_transpose_print ( m, n, x, "  Grid points:" );
  delete [] x;

  return;
# undef M
}
