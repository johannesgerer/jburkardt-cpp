# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>

using namespace std;

# include "sphere_grid.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_GRID_PRB.
//
//  Discussion:
//
//    SPHERE_GRID_PRB tests routines from the SPHERE_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "SPHERE_GRID_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_GRID library.\n";
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_GRID_PRB\n";
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
//    TEST01 tests SPHERE_ICOS_POINT_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int edge_num;
  int factor;
  int factor_log;
  int node_num;
  int triangle_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SPHERE_ICOS_POINT_NUM determines the size\n";
  cout << "  (number of vertices, edges and faces) in a grid\n";
  cout << "  on a sphere, made by subdividing an initial\n";
  cout << "  projected icosahedron.\n";
  cout << "\n";
  cout << "  N determines the number of subdivisions of each\n";
  cout << "  edge of the icosahedral faces.\n";
  cout << "\n";
  cout << "         N         V         E         F\n";
  cout << "  --------  --------  --------  --------\n";
  cout << "\n";

  for ( factor = 1; factor <= 20; factor++ )
  {
    node_num = sphere_icos_point_num ( factor );
    edge_num = sphere_icos_edge_num ( factor );
    triangle_num = sphere_icos_face_num ( factor );
    cout << "  " << setw(8) << factor
         << "  " << setw(8) << node_num
         << "  " << setw(8) << edge_num
         << "  " << setw(8) << triangle_num << "\n";
  }

  cout << "\n";
  cout << "  Repeat, but using N constrained by doubling:\n";
  cout << "\n";
  cout << "         N         V         E         F\n";
  cout << "  --------  --------  --------  --------\n";
  cout << "\n";

  factor = 1;
  for ( factor_log = 0; factor_log <= 10; factor_log++ )
  {
    node_num = sphere_icos_point_num ( factor );
    edge_num = sphere_icos_edge_num ( factor );
    triangle_num = sphere_icos_face_num ( factor );
    cout << "  " << setw(8) << factor
         << "  " << setw(8) << node_num
         << "  " << setw(8) << edge_num
         << "  " << setw(8) << triangle_num << "\n";
    factor = factor * 2;
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests SPHERE_ICOS1_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int factor;
  string filename;
  int node;
  int node_num;
  double *node_xyz;
  ofstream output;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  SPHERE_ICOS_POINT_NUM \"sizes\" a grid generated\n";
  cout << "  on an icosahedron and projected to a sphere.\n";
  cout << "  SPHERE_ICOS1_POINTS creates the grid points.\n";

  factor = 3;

  cout << "\n";
  cout << "  Sizing factor =       " << factor << "\n";

  node_num = sphere_icos_point_num ( factor );

  cout << "\n";
  cout << "  Number of vertices =  " << node_num << "\n";

  node_xyz = sphere_icos1_points ( factor, node_num );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 20, 
    "  Initial part of NODE_XYZ array:" );
//
//  Write the nodes to a file.
//
  if ( true )
  {
    filename = "sphere_icos1_points_f" + i4_to_string ( factor, "%d" ) + ".xyz";

    output.open ( filename.c_str ( ) );
    for ( node = 0; node < node_num; node++ )
    {
      output << "  " << node_xyz[0+node*3]
             << "  " << node_xyz[1+node*3]
             << "  " << node_xyz[2+node*3] << "\n";
    }
    output.close ( );

    cout << "\n";
    cout << "  Wrote data to \"" << filename << "\"\n";
  }

  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests SPHERE_ICOS2_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int factor;
  string filename;
  int node;
  int node_num;
  double *node_xyz;
  ofstream output;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  SPHERE_ICOS_POINT_NUM \"sizes\" a grid generated\n";
  cout << "  on an icosahedron and projected to a sphere.\n";
  cout << "  SPHERE_ICOS2_POINTS creates the grid.\n";

  factor = 3;

  cout << "\n";
  cout << "  Sizing factor FACTOR = " << factor << "\n";

  node_num = sphere_icos_point_num ( factor );

  cout << "\n";
  cout << "  Number of nodes =     " << node_num << "\n";

  node_xyz = sphere_icos2_points ( factor, node_num );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 20, 
    "  Initial part of NODE_XYZ array:" );
//
//  Write the nodes to a file.
//
  if ( true )
  {
    filename = "sphere_icos2_points_f" + i4_to_string ( factor, "%d" ) + ".xyz";

    output.open ( filename.c_str ( ) );
    for ( node = 0; node < node_num; node++ )
    {
      output << "  " << node_xyz[0+node*3]
             << "  " << node_xyz[1+node*3]
             << "  " << node_xyz[2+node*3] << "\n";
    }
    output.close ( );

    cout << "\n";
    cout << "  Wrote data to \"" << filename << "\"\n";
  }

  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests SPHERE_LL_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num = 3;
  int lon_num = 4;

  double pc[3] = { 0.0, 0.0, 0.0 };
  int i;
  int j;
  int k;
  int node_num;
  double *node_xyz;
  double r = 10.0;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  SPHERE_LL_POINTS produces latitude/longitude\n";
  cout << "  points on a sphere in 3D.\n";

  cout << "\n";
  cout << "  Radius = " << r << "\n";

  r8vec_print ( 3, pc, "  Center:" );

  cout << "\n";
  cout << "  The number of latitudes =  " << lat_num << "\n";
  cout << "  The number of longitudes = " << lon_num << "\n";

  node_num = sphere_ll_point_num ( lat_num, lon_num );
  cout << "\n";
  cout << "  The number of grid points is " << node_num << "\n";

  node_xyz = sphere_ll_points ( r, pc, lat_num, lon_num, node_num );

  cout << "\n";

  k = 0;
  cout << "  " << setw(8) << k
       << "  " << setw(12) << node_xyz[0+k*3]
       << "  " << setw(12) << node_xyz[1+k*3]
       << "  " << setw(12) << node_xyz[2+k*3] << "\n";

  for ( i = 1; i <= lat_num; i++ )
  {
    cout << "\n";
    for ( j = 0; j < lon_num; j++ )
    {
      k = k + 1;
      cout << "  " << setw(8) << k
           << "  " << setw(12) << node_xyz[0+k*3]
           << "  " << setw(12) << node_xyz[1+k*3]
           << "  " << setw(12) << node_xyz[2+k*3] << "\n";
    }
  }

  cout << "\n";
  k = k + 1;
  cout << "  " << setw(8) << k
       << "  " << setw(12) << node_xyz[0+k*3]
       << "  " << setw(12) << node_xyz[1+k*3]
       << "  " << setw(12) << node_xyz[2+k*3] << "\n";

  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests SPHERE_SPIRALPOINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  double center_xyz[3] = { 0.0, 0.0, 0.0 };
  int node;
  int node_num = 500;
  double *node_xyz;
  ofstream output;
  double r = 1.0;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  SPHERE_SPIRALPOINTS produces a spiral of\n";
  cout << "  points on an implicit sphere in 3D.\n";
  cout << "\n";
  cout << "  Radius = " << r << "\n";

  r8vec_print ( 3, center_xyz, "  Center:" );

  cout << "\n";
  cout << "  The number of spiral points is " << node_num << "\n";

  node_xyz = sphere_spiralpoints ( r, center_xyz, node_num );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 10, 
    "  The spiral points:" );
//
//  Write the nodes to a file.
//
  if ( true )
  {
    filename = "sphere_grid_spiral_n" + i4_to_string ( node_num, "%d" ) + ".xyz";

    output.open ( filename.c_str ( ) );
    for ( node = 0; node < node_num; node++ )
    {
      output << "  " << node_xyz[0+node*3]
             << "  " << node_xyz[1+node*3]
             << "  " << node_xyz[2+node*3] << "\n";
    }
    output.close ( );

    cout << "\n";
    cout << "  Wrote data to \"" << filename << "\"\n";
  }

  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests SPHERE_LL_LINES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num = 3;
  int *line;
  int line_num;
  int long_num = 4;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  SPHERE_LL_LINES computes latitude/longitude\n";
  cout << "  lines on a sphere in 3D.\n";
  cout << "\n";
  cout << "  Number of latitudes is  " << lat_num << "\n";
  cout << "  Number of longitudes is " << long_num << "\n";

  line_num = sphere_ll_line_num ( lat_num, long_num );

  cout << "\n";
  cout << "  Number of line segments is " << line_num << "\n";

  line = sphere_ll_lines ( lat_num, long_num, line_num );

  i4mat_transpose_print ( 2, line_num, line, "  Grid line vertices:" );

  free ( line );

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests SPHERE_GRID_Q4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num = 3;
  int long_num = 4;
  int rectangle_num = lat_num * long_num;
  int *rectangle_node;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  SPHERE_GRID_Q4 computes a grid\n";
  cout << "  of Q4 rectangular elements on a sphere in 3D.\n";
  cout << "\n";
  cout << "  Number of latitudes is      " << lat_num << "\n";
  cout << "  Number of longitudes is     " << long_num << "\n";
  cout << "  The number of rectangles is " << rectangle_num << "\n";

  rectangle_node = sphere_grid_q4 ( lat_num, long_num );

  i4mat_transpose_print ( 4, rectangle_num, rectangle_node, 
    "  Rectangle vertices:" );

  delete [] rectangle_node;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests SPHERE_GRID_T3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int lat_num = 3;
  int lon_num = 4;

  int triangle_num;
  int *triangle_node;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  SPHERE_GRID_T3 computes a grid\n";
  cout << "  of T3 triangular elements on a sphere in 3D.\n";
  cout << "\n";
  cout << "  Number of latitudes is  " << lat_num << "\n";
  cout << "  Number of longitudes is " << lon_num << "\n";

  triangle_node = sphere_grid_t3 ( lat_num, lon_num );

  triangle_num = 2 * ( lat_num + 1 ) * lon_num;

  cout << "\n";
  cout << "  The number of triangles is " << triangle_num << "\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_node, 
    "  Triangle vertices:" );

  delete [] triangle_node;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests SPHERE_UNIT_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int node;
  int node_num;
  double *node_xyz;
  ofstream output;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For the unit sphere in 3 dimensions:\n";
  cout << "  SPHERE_UNIT_SAMPLE does a random sampling.\n";

  node_num = 1000;

  node_xyz = sphere_unit_sample ( node_num, &seed );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 10, 
    "  First 10 values:" );
//
//  Write the nodes to a file.
//
  if ( true )
  {
    filename = "sphere_sample_n" + i4_to_string ( node_num, "%d" ) + ".xyz";

    output.open ( filename.c_str ( ) );
    for ( node = 0; node < node_num; node++ )
    {
      output << "  " << node_xyz[0+node*3]
             << "  " << node_xyz[1+node*3]
             << "  " << node_xyz[2+node*3] << "\n";
    }
    output.close ( );

    cout << "\n";
    cout << "  Wrote data to \"" << filename << "\"\n";
  }

  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests SPHERE_CUBED_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int j;
  int n;
  int ns;
  ofstream output;
  double *xyz;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  SPHERE_CUBED_POINTS computes points on a cubed sphere grid.\n";

  n = 10;
  cout << "\n";
  cout << "  Number of divisions on each face = " << n << "\n";

  ns = sphere_cubed_point_num ( n );
  cout << "  Total number of points = " << ns << "\n";

  xyz = sphere_cubed_points ( n, ns );

  r8mat_transpose_print_some ( 3, ns, xyz, 1, 1, 3, 20, "  Initial part of XYZ array:" );
/*
  Write the nodes to a file.
*/
  if ( 1 )
  {
    filename = "sphere_cubed_f" + i4_to_string ( n, "%d" ) + ".xyz";

    output.open ( filename.c_str ( ) );
    for ( j = 0; j < n; j++ )
    {
      output << "  " << xyz[0+j*3]
             << "  " << xyz[1+j*3]
             << "  " << xyz[2+j*3] << "\n";
    }
    output.close ( );

    cout << "\n";
    cout << "  Wrote data to \"" << filename << "\"\n";
  }

  delete [] xyz;

  return;
}
