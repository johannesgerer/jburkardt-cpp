# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "ellipse_grid.hpp"

int main ( );
void ellipse_grid_test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ELLIPSE_GRID_PRB.
//
//  Discussion:
//
//    ELLIPSE_GRID_PRB tests the ELLIPSE_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ELLIPSE_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ELLIPSE_GRID library.\n";

  ellipse_grid_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ELLIPSE_GRID_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void ellipse_grid_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_GRID_TEST01 tests ELLIPSE_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  double c[2];
  string filename = "ellipse_grid_test01.xy";
  int n;
  int ng;
  double r[2];
  double *xy;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  ELLIPSE_GRID can define a grid of points\n";
  cout << "  with N+1 points on the minor half axis,\n";
  cout << "  based on any ellipse.\n";

  n = 8;
  r[0] = 2.0;
  r[1] = 1.0;
  c[0] = 1.0;
  c[1] = 2.0;

  cout << "\n";
  cout << "  We use N = " << n << "\n";
  cout << "  Radius R = (" << r[0] << "," << r[1] << ")\n";
  cout << "  Center C = (" << c[0] << "," << c[1] << ")\n";

  ng = ellipse_grid_count ( n, r, c );

  cout << "\n";
  cout << "  Number of grid points will be " << ng << "\n";

  xy = ellipse_grid ( n, r, c, ng );

  r82vec_print_part ( ng, xy, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 2, ng, xy );

  cout << "\n";
  cout << "  Data written to the file \"" << filename << "\".\n";

  delete [] xy;

  return;
}
