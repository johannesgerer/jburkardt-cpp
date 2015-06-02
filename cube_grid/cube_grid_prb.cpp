# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "cube_grid.hpp"

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
//    MAIN is the main program for CUBE_GRID_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CUBE_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the CUBE_GRID library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CUBE_GRID_PRB:\n";
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
//    TEST01 tests CUBE_GRID using the same parameters for all dimensions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[3] = { -1.0, -1.0, -1.0 };
  double b[3] = { +1.0, +1.0, +1.0 };
  int c[3] = { 1, 1, 1 };
  int i;
  int n;
  int ns[3] = { 3, 3, 3 };
  double *x;

  n = ns[0] * ns[1] * ns[2];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Create a grid using CUBE_GRID.\n";
  cout << "  Use the same parameters in every dimension.\n";
  cout << "  Number of grid points N = " << n << "\n";
  cout << "\n";
  cout << "     I    NS     C      A         B\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout 
      << setw(6) << i << "  "
      << setw(4) << ns[i] << "  "
      << setw(4) << c[i] << "  "
      << setw(8) << a[i] << "  "
      << setw(8) << b[i] << "\n";
  }

  x = cube_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 3, n, x, "  Grid points:" );
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses a different number of points in each coordinate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[3] = { 0.0, 0.0, 0.0 };
  double b[3] = { 1.0, 1.0, 1.0  };
  int c[3] = { 2, 2, 2 };
  int i;
  int n;
  int ns[3] = { 4, 2, 3 };
  double *x;

  n = ns[0] * ns[1] * ns[2];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Create a grid using CUBE_GRID.\n";
  cout << "  se a different number of points in each dimension..\n";
  cout << "  Number of grid points N = " << n << "\n";
  cout << "\n";
  cout << "     I    NS     C      A         B\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout 
      << setw(6) << i << "  "
      << setw(4) << ns[i] << "  "
      << setw(4) << c[i] << "  "
      << setw(8) << a[i] << "  "
      << setw(8) << b[i] << "\n";
  }

  x = cube_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 3, n, x, "  Grid points:" );
  delete [] x;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses a cube with different sizes in each dimension.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[3] = {   0.0, -2.0, 50.0 };
  double b[3] = { +10.0, +2.0, 51.0 };
  int c[3] = { 3, 4, 5 };
  int i;
  int n;
  int ns[3] = { 3, 3, 3 };
  double *x;

  n = ns[0] * ns[1] * ns[2];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Create a grid using CUBE_GRID.\n";
  cout << "  Use a different physical size in each dimension.\n";
  cout << "  Number of grid points N = " << n << "\n";
  cout << "\n";
  cout << "     I    NS     C      A         B\n";
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout 
      << setw(6) << i << "  "
      << setw(4) << ns[i] << "  "
      << setw(4) << c[i] << "  "
      << setw(8) << a[i] << "  "
      << setw(8) << b[i] << "\n";
  }

  x = cube_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 3, n, x, "  Grid points:" );
  delete [] x;

  return;
}
