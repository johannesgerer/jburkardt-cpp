# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "sphere_llt_grid.hpp"

int main ( );
void sphere_llt_grid_point_count_test ( );
void sphere_llt_grid_points_test ( );
void sphere_llt_grid_line_count_test ( );
void sphere_llt_grid_lines_test ( );
void sphere_llt_grid_display_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_LLT_GRID_PRB.
//
//  Discussion:
//
//    SPHERE_LLT_GRID_PRB tests the SPHERE_LLT_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPHERE_LLT_GRID_TEST\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_LLT_GRID library.\n";

  sphere_llt_grid_point_count_test ( );
  sphere_llt_grid_points_test ( );
  sphere_llt_grid_line_count_test ( );
  sphere_llt_grid_lines_test ( );
  sphere_llt_grid_display_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_LLT_GRID_TEST\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void sphere_llt_grid_point_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_POINT_COUNT_TEST tests SPHERE_LLT_GRID_POINT_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num;
  int long_log;
  int long_num;
  int point_num;

  cout << "\n";
  cout << "SPHERE_LLT_GRID_POINT_COUNT_TEST\n";
  cout << "  SPHERE_LLT_GRID_POINT_COUNT counts the points used for a\n";
  cout << "  grid based on triangles defined by latitude and longitude\n";
  cout << "  lines on a sphere in 3D.\n";
  cout << "\n";
  cout << "     LAT_NUM    LONG_NUM   POINT_NUM\n";

  for ( lat_num = 1; lat_num <= 17; lat_num = lat_num + 2 );
  {
    cout << "\n";
    long_num = 1;
    for ( long_log = 1; long_log <= 4; long_log++ )
    {
      long_num = long_num * 2;
      point_num = sphere_llt_grid_point_count ( lat_num, long_num );
      cout << "  " << setw(8) << lat_num
           << "  " << setw(8) << long_num
           << "  " << setw(8) << point_num << "\n";
    }
  }

  return;
}
//****************************************************************************80

void sphere_llt_grid_points_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_POINTS_TEST tests SPHERE_LLT_GRID_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int lat_num;
  int long_num;
  int node_num;
  double *node_xyz;
  double pc[3] = { 0.0, 0.0, 0.0 };
  double r;

  lat_num = 3;
  long_num = 4;

  r = 10.0;

  cout << "\n";
  cout << "SPHERE_LLT_GRID_POINTS_TEST\n";
  cout << "  SPHERE_LLT_POINTS produces latitude/longitude\n";
  cout << "  points on a sphere in 3D.\n";

  cout << "\n";
  cout << "  Radius = " << r << "\n";

  r8vec_print ( 3, pc, "  Center:" );

  cout << "\n";
  cout << "  Number of latitudes is  " << lat_num << "\n";
  cout << "  Number of longitudes is " << long_num << "\n";

  node_num = sphere_llt_grid_point_count ( lat_num, long_num );

  cout << "\n";
  cout << "  The number of grid points is " << node_num << "\n";

  node_xyz = sphere_llt_grid_points ( r, pc, lat_num, long_num, node_num );

  cout << "\n";

  k = 0;
  cout << "  " << setw(8) << k
       << "  " << setw(14) << node_xyz[0+k*3]
       << "  " << setw(14) << node_xyz[1+k*3]
       << "  " << setw(14) << node_xyz[2+k*3] << "\n";
  k = k + 1;

  cout << "\n";

  for ( i = 1; i <= lat_num; i++ )
  {
    cout << "\n";
    for ( j = 0; j < long_num; j++ )
    {
      cout << "  " << setw(8) << k
           << "  " << setw(14) << node_xyz[0+k*3]
           << "  " << setw(14) << node_xyz[1+k*3]
           << "  " << setw(14) << node_xyz[2+k*3] << "\n";
      k = k + 1;
      cout << "\n";
    }
  }

  cout << "\n";

  cout << "  " << setw(8) << k
       << "  " << setw(14) << node_xyz[0+k*3]
       << "  " << setw(14) << node_xyz[1+k*3]
       << "  " << setw(14) << node_xyz[2+k*3] << "\n";
  k = k + 1;
  cout << "\n";

  delete [] node_xyz;

  return;
}
//****************************************************************************80

void sphere_llt_grid_line_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_LINE_COUNT_TEST tests SPHERE_LLT_GRID_LINE_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num;
  int line_num;
  int long_log;
  int long_num;

  lat_num = 3;
  long_num = 4;

  cout << "\n";
  cout << "SPHERE_LLT_GRID_LINE_COUNT_TEST\n";
  cout << "  SPHERE_LLT_GRID_LINE_COUNT counts the lines used for a\n";
  cout << "  grid based on triangles defined by latitude and longitude\n";
  cout << "  lines on a sphere in 3D.\n";
  cout << "\n";
  cout << "     LAT_NUM    LONG_NUM   LINE_NUM\n";

  for ( lat_num = 1; lat_num <= 17; lat_num = lat_num + 2 )
  {
    cout << "\n";
    long_num = 1;
    for ( long_log = 1; long_log <= 4; long_log++ )
    {
      long_num = long_num * 2;
      line_num = sphere_llt_grid_line_count ( lat_num, long_num );
      cout << "  " << setw(8) << lat_num
           << "  " << setw(8) << long_num
           << "  " << setw(8) << line_num << "\n";
    }
  }

  return;
}
//****************************************************************************80

void sphere_llt_grid_lines_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_LINES_TEST tests SPHERE_LLT_GRID_LINES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num;
  int *line_data;
  int line_num;
  int long_num;

  lat_num = 3;
  long_num = 4;

  cout << "\n";
  cout << "SPHERE_LLT_GRID_LINES_TEST\n";
  cout << "  SPHERE_LLT_GRID_LINES computes grid lines\n";
  cout << "  on a sphere in 3D.\n";
  cout << "\n";
  cout << "  Number of latitudes is  " << lat_num << "\n";
  cout << "  Number of longitudes is " << long_num << "\n";

  line_num = sphere_llt_grid_line_count ( lat_num, long_num );

  cout << "\n";
  cout << "  Number of line segments is " << line_num << "\n";

  line_data = sphere_llt_grid_lines ( lat_num, long_num, line_num );

  i4mat_transpose_print ( 2, line_num, line_data, 
    "  Grid line vertices:" );

  delete [] line_data;

  return;
}
//****************************************************************************80

void sphere_llt_grid_display_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_DISPLAY_TEST tests SPHERE_LLT_GRID_DISPLAY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num;
  int *line_data;
  int line_num;
  int long_num;
  int node_num;
  double *node_xyz;
  double pc[3] = { 0.0, 0.0, 0.0 };
  string prefix;
  double r;

  lat_num = 10;
  long_num = 12;

  r = 10.0;

  cout << "\n";
  cout << "SPHERE_LLT_GRID_DISPLAY_TEST\n";
  cout << "  SPHERE_LLT_GRID_DISPLAY displays an LLT grid on a sphere.\n";
  cout << "\n";
  cout << "  Number of latitudes is  " << lat_num << "\n";
  cout << "  Number of longitudes is " << long_num << "\n";
//
//  Get points.
//
  node_num = sphere_llt_grid_point_count ( lat_num, long_num );

  cout << "\n";
  cout << "  The number of grid points is " << node_num << "\n";

  node_xyz = sphere_llt_grid_points ( r, pc, lat_num, long_num, node_num );
//
//  Get lines.
//
  line_num = sphere_llt_grid_line_count ( lat_num, long_num );

  cout << "\n";
  cout << "  Number of line segments is " << line_num << "\n";

  line_data = sphere_llt_grid_lines ( lat_num, long_num, line_num );

  prefix = "sphere_llt_grid";

  sphere_llt_grid_display ( node_num, node_xyz, line_num, line_data, prefix );

  delete [] line_data;
  delete [] node_xyz;

  return;
}
