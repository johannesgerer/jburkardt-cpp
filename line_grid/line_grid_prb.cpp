# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "line_grid.hpp"

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
//    MAIN is the main program for LINE_GRID_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LINE_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LINE_GRID library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LINE_GRID_PRB:\n";
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
//    TEST01 tests LINE_GRID using simple parameters.
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
  double a = -1.0;
  double b = +1.0;
  int c = 1;
  int n = 11;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Create a grid using LINE_GRID.\n";
  cout << "  Use simple parameters.\n";
  cout << "\n";
  cout << "     N     C      A         B\n";
  cout << "\n";
  cout 
    << setw(4) << n << "  "
    << setw(4) << c << "  "
    << setw(8) << a << "  "
    << setw(8) << b << "\n";

  x = line_grid ( n, a, b, c );
  r8vec_print ( n, x, "  Grid points:" );
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 changes the number of points.
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
  double a = 0.0;
  double b = 1.0;
  int c = 2;
  int n;
  int test;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Create a grid using LINE_GRID.\n";
  cout << "  Try an increasing number of points.\n";

  n = 4;

  for ( test = 1; test <= 3; test++ )
  {
    n = 2 * n + 1;

    cout << "\n";
    cout << "     N     C      A         B\n";
    cout << "\n";
    cout 
      << setw(4) << n << "  "
      << setw(4) << c << "  "
      << setw(8) << a << "  "
      << setw(8) << b << "\n";

    x = line_grid ( n, a, b, c );
    r8vec_print ( n, x, "  Grid points:" );
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tries all the centering options.
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
  double a = 0.0;
  double b = 100.0;
  int c;
  int n;
  double *x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Create a grid using LINE_GRID.\n";
  cout << "  Try the different centering options.\n";

  for ( c = 1; c <= 5; c++ )
  {
    cout << "\n";
    cout << "     N     C      A         B\n";
    cout << "\n";
    cout 
      << setw(4) << n << "  "
      << setw(4) << c << "  "
      << setw(8) << a << "  "
      << setw(8) << b << "\n";

    x = line_grid ( n, a, b, c );
    r8vec_print ( n, x, "  Grid points:" );
    delete [] x;
  }

  return;
}
