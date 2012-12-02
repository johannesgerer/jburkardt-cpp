# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "ball_grid.hpp"

int main ( );
void ball_grid_test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_GRID_TEST tests BALL_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BALL_GRID_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the BALL_GRID library.\n";

  ball_grid_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BALL_GRID_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void ball_grid_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_GRID_TEST01 tests BALL_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *bg;
  double c[3];
  string filename = "ball_grid_test01.xyz";
  int n;
  int ng;
  double r;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  BALL_GRID can define a grid of points\n";
  cout << "  with N+1 points on a horizontal or vertical radius,\n";
  cout << "  based on any ball.\n";

  n = 10;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;
  c[2] = 2.0;

  cout << "\n";
  cout << "  We use N = " << n << "\n";
  cout << "  Radius R = " << r << "\n";
  cout << "  Center C = (" << c[0] << "," << c[1] << "," << c[2] << ")\n";

  ng = ball_grid_count ( n, r, c );

  cout << "\n";
  cout << "  Number of grid points will be " << ng << "\n";

  bg = ball_grid ( n, r, c, ng );

  r83vec_print_part ( ng, bg, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 3, ng, bg );

  cout << "\n";
  cout << "  Data written to the file \"" << filename << "\".\n";

  return;
}
