# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "circle_arc_grid.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CIRCLE_ARC_GRID_PRB.
//
//  Discussion:
//
//    CIRCLE_ARC_GRID_PRB tests CIRCLE_ARC_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CIRCLE_ARC_GRID_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test CIRCLE_ARC_GRID.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CIRCLE_ARC_GRID_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 demonstrates the use of CIRCLE_ARC_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a[2];
  double c[2];
  string filename = "arc.txt";
  int n;
  double r;
  double *xy;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compute points along a 90 degree arc\n";

  r = 2.0;
  c[0] = 5.0;
  c[1] = 5.0;
  a[0] = 0.0;
  a[1] = 90.0;
  n = 10;
//
//  Echo the input.
//
  cout << "\n";
  cout << "  Radius =           " << r << "\n";
  cout << "  Center =           " << c[0] << "  " << c[1] << "\n";
  cout << "  Angle 1 =          " << a[0]  << "\n";
  cout << "  Angle 2 =          " << a[1]  << "\n";;
  cout << "  Number of points = " << n  << "\n";
//
//  Compute the data.
//
  xy = circle_arc_grid ( r, c, a, n );
//
//  Print a little of the data.
//
  r82vec_print_part ( n, xy, 5, "  A few of the points:" );
//
//  Write the data.
//
  r8mat_write ( filename, 2, n, xy );
  cout << "\n";
  cout << "  Data written to \"" << filename << "\".\n";
//
//  Free memory.
//
  delete [] xy;

  return;
}
