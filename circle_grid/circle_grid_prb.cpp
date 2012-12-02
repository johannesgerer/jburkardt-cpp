# include <cstdlib>
# include <iostream>
# include <cmath>

using namespace std;

# include "circle_grid.hpp"

int main ( );
void circle_grid_test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_GRID_TEST tests CIRCLE_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CIRCLE_GRID_TEST:\n";
  cout << "  C version\n";
  cout << "  Test the CIRCLE_GRID library.\n";

  circle_grid_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CIRCLE_GRID_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void circle_grid_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_GRID_TEST01 tests CIRCLE_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  double c[2];
  double *cg;
  string filename = "circle_grid_test01.xy";
  int n;
  int ng;
  double r;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  CIRCLE_GRID can define a grid of points\n";
  cout << "  with N+1 points on a horizontal or vertical radius,\n";
  cout << "  based on any circle.\n";

  n = 20;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;

  cout << "\n";
  cout << "  We use N = " << n << "\n";
  cout << "  Radius R = " << r << "\n";;
  cout << "  Center C = (" << c[0] << "," << c[1] << ")\n";

  ng = circle_grid_count ( n, r, c );

  cout << "\n";
  cout << "  Number of grid points will be " << ng << "\n";

  cg = circle_grid ( n, r, c, ng );

  r82vec_print_part ( ng, cg, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 2, ng, cg );

  cout << "\n";
  cout << "  Data written to the file \"" << filename << "\"\n";

  delete [] cg;

  return;
}
