# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "sphere_fibonacci_grid.hpp"

int main ( );
void sphere_fibonacci_grid_points_test ( );
void sphere_fibonacci_grid_display_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_FIBONACCI_GRID_TEST tests the SPHERE_FIBONACCI_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPHERE_FIBONACCI_GRID_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_FIBONACCI_GRID library.\n";

  sphere_fibonacci_grid_points_test ( );
  sphere_fibonacci_grid_display_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_FIBONACCI_GRID_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void sphere_fibonacci_grid_points_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_FIBONACCI_GRID_POINTS_TEST tests SPHERE_FIBONACCI_GRID_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int ng;
  double *xg;

  cout << "\n";
  cout << "SPHERE_FIBONACCI_GRID_POINTS_TEST\n";
  cout << "  SPHERE_FIBONACCI_GRID_POINTS computes points on a sphere\n";
  cout << "  that lie on a Fibonacci spiral.\n";

  ng = 1000;
  cout << "\n";
  cout << "  Number of points NG = " << ng << "\n";

  xg = sphere_fibonacci_grid_points ( ng );

  r8mat_transpose_print_some ( 3, ng, xg, 1, 1, 3, 10, 
    "  Part of the grid array:" );
//
//  Write the nodes to a file.
//
  filename = "sphere_fibonacci_grid_n1000.xyz";

  r8mat_write ( filename, 3, ng, xg );

  delete [] xg;

  return;
}
//****************************************************************************80

void sphere_fibonacci_grid_display_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_FIBONACCI_GRID_DISPLAY_TEST tests SPHERE_FIBONACCI_GRID_DISPLAY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int ng;
  string prefix;
  double *xg;

  cout << "\n";
  cout << "SPHERE_FIBONACCI_GRID_DISPLAY_TEST\n";
  cout << "  SPHERE_FIBONACCI_GRID_DISPLAY displays points\n";
  cout << "  on a sphere that lie on a Fibonacci spiral.\n";

  ng = 1000;
  cout << "\n";
  cout << "  Number of points NG = " << ng << "\n";

  xg = sphere_fibonacci_grid_points ( ng );
//
//  Display the nodes on a sphere.
//
  prefix = "sphere_fibonacci_grid_n1000";

  sphere_fibonacci_grid_display ( ng, xg, prefix );

  delete [] xg;

  return;
}
