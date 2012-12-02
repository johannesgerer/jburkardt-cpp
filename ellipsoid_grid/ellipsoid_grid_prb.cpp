# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "ellipsoid_grid.hpp"

int main ( );
void ellipsoid_grid_test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSOID_GRID_TEST tests ELLIPSOID_GRID.
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
  cout << "ELLIPSOID_GRID_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the ELLIPSOID_GRID library.\n";

  ellipsoid_grid_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ELLIPSOID_GRID_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void ellipsoid_grid_test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSOID_GRID_TEST01 tests ELLIPSOID_GRID.
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
  double c[3];
  string filename = "ellipsoid_grid_test01.xyz";
  int n;
  int ng;
  double r[3];
  double *xyz;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  ELLIPSOID_GRID can define a grid of points\n";
  cout << "  with N+1 points on the minor half axis,\n";
  cout << "  based on any ellipsoid.\n";

  n = 4;
  r[0] = 2.0;
  r[1] = 1.0;
  r[2] = 1.5;
  c[0] = 1.0;
  c[1] = 2.0;
  c[2] = 1.5;

  cout << "\n";
  cout << "  We use N = " << n << "\n";
  cout << "  Radius R = (" << r[0] << "," << r[1] << "," << r[2] << ")\n";
  cout << "  Center C = (" << c[0] << "," << c[1] << "," << c[2] << ")\n";

  ng = ellipsoid_grid_count ( n, r, c );

  cout << "\n";
  cout << "  Number of grid points will be " << ng << "\n";

  xyz = ellipsoid_grid ( n, r, c, ng );

  r83vec_print_part ( ng, xyz, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 3, ng, xyz );

  cout << "\n";
  cout << "  Data written to the file \"" << filename << "\".\n";

  delete [] xyz;

  return;
}
