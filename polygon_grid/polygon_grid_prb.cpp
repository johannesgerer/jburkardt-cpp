# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "polygon_grid.hpp"

int main ( );
void polygon_grid_count_test ( );
void polygon_grid_display_test ( );
void polygon_grid_points_test01 ( );
void polygon_grid_points_test02 ( );
void polygon_grid_points_test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_TEST tests the POLYGON_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POLYGON_GRID_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the POLYGON_GRID library.\n";

  polygon_grid_count_test ( );

  polygon_grid_points_test01 ( );
  polygon_grid_points_test02 ( );
  polygon_grid_points_test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLYGON_GRID_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void polygon_grid_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_COUNT_TEST tests POLYGON_GRID_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int ng;
  int nv;

  cout << "\n";
  cout << "POLYGON_GRID_COUNT_TEST:\n";
  cout << "  POLYGON_GRID_COUNT counts NG, the number of points in\n";
  cout << "  a grid defined with N+1 points on each side of a\n";
  cout << "  polygon of NV vertices.\n";

  for ( nv = 3; nv <= 5; nv++ )
  {
    cout << "\n";
    cout << "  Polygonal vertex count NV = " << nv << "\n";
    cout << "\n";
    cout << "   N     NG\n";
    cout << "\n";
    for ( n = 0; n <= 5; n++ )
    {
      ng = polygon_grid_count ( n, nv );
      cout << "  " << setw(2) << n
           << "  " << setw(5) << ng << "\n";
    }
  }

  return;
}
//****************************************************************************80

void polygon_grid_points_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_POINTS_TEST01 tests POLYGON_GRID_POINTS
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int n;
  int ng;
  int nv = 3;
  string prefix;
  double v[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.5, 0.86602540378443860 };
  double *xg;

  cout << "\n";
  cout << "POLYGON_GRID_POINTS_TEST01:\n";
  cout << "  POLYGON_GRID_POINTS returns grid points for a polygon\n";
  cout << "  of NV vertices, with N+1 points on a side\n";
  cout << "\n";
  cout << "  For this test, the polygon is a triangle.\n";

  r8mat_transpose_print ( 2, nv, v, "  Polygon vertices:" );
//
//  Count the grid points.
//
  n = 5;
  ng = polygon_grid_count ( n, nv );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  Number of grid points will be NG = " << ng << "\n";
//
//  Compute the grid points.
//
  xg = polygon_grid_points ( n, nv, v, ng );

  r8mat_transpose_print ( 2, ng, xg, "  The grid point array:" );
//
//  Display the points.
//
  prefix = "triangle";

  polygon_grid_display ( n, nv, v, ng, xg, prefix );
//
//  Write the points to a file.
//
  filename = "triangle.xy";

  r8mat_write ( filename, 2, ng, xg );

  cout << "\n";
  cout << "  Data written to the file '" << filename << "'\n";

  delete [] xg;

  return;
}
//****************************************************************************80

void polygon_grid_points_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_POINTS_TEST02 tests POLYGON_GRID_POINTS
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int n;
  int ng;
  int nv = 4;
  string prefix;
  double v[2*4] = {
    1.0, 1.0, 
    2.0, 0.0, 
    4.0, 3.0, 
    0.0, 5.0 };
  double *xg;

  cout << "\n";
  cout << "POLYGON_GRID_POINTS_TEST02:\n";
  cout << "  POLYGON_GRID_POINTS returns grid points for a polygon\n";
  cout << "  of NV vertices, with N+1 points on a side\n";
  cout << "\n";
  cout << "  For this test, the polygon is a convex quadrilateral\n";
  cout << "  with sides of varying length.\n";
//
//  Define the polygon.
//
  r8mat_transpose_print ( 2, nv, v, "  Polygon vertices:" );
//
//  Count the grid points.
//
  n = 7;
  ng = polygon_grid_count ( n, nv );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  Number of grid points will be NG = " << ng << "\n";
//
//  Compute the grid points.
//
  xg = polygon_grid_points ( n, nv, v, ng );

  r8mat_transpose_print ( 2, ng, xg, "  The grid point array:" );
//
//  Display the points.
//
  prefix = "quad";

  polygon_grid_display ( n, nv, v, ng, xg, prefix );
//
//  Write the points to a file.
//
  filename = "quad.xy";

  r8mat_write ( filename, 2, ng, xg );

  cout << "\n";
  cout << "  Data written to the file '" << filename << "'\n";

  delete [] xg;

  return;
}
//****************************************************************************80

void polygon_grid_points_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_POINTS_TEST03 tests POLYGON_GRID_POINTS
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int n;
  int ng;
  int nv = 6;
  string prefix;
  double v[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    2.0, 1.0, 
    1.0, 1.0, 
    1.0, 2.0, 
    0.0, 2.0 };
  double *xg;

  cout << "\n";
  cout << "POLYGON_GRID_POINTS_TEST03:\n";
  cout << "  POLYGON_GRID_POINTS returns grid points for a polygon\n";
  cout << "  of NV vertices, with N+1 points on a side\n";
  cout << "\n";
  cout << "  For this test, the polygon is nonconvex and six sided.\n";
  cout << "  Two degenerate triangles are created, and some grid points\n";
  cout << "  are generated several times.\n";
//
//  Define the polygon.
//
  r8mat_transpose_print ( 2, nv, v, "  Polygon vertices:" );
//
//  Count the grid points.
//
  n = 5;
  ng = polygon_grid_count ( n, nv );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  Number of grid points will be NG = " << ng << "\n";
//
//  Compute the grid points.
//
  xg = polygon_grid_points ( n, nv, v, ng );

  r8mat_transpose_print ( 2, ng, xg, "  The grid point array:" );
//
//  Display the points.
//
  prefix = "ell";

  polygon_grid_display ( n, nv, v, ng, xg, prefix );
//
//  Write the points to a file.
//
  filename = "ell.xy";

  r8mat_write ( filename, 2, ng, xg );

  cout << "\n";
  cout << "  Data written to the file '" << filename << "'\n";

  delete [] xg;

  return;
}

