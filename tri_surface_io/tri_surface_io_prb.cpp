# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "tri_surface_io.hpp"

int main ( );
void test01 ( string node_file_name, string triangle_file_name );
void test02 ( string node_file_name, string triangle_file_name );
void test03 ( string node_file_name, string triangle_file_name );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRI_SURFACE_IO_PRB runs the TRI_SURFACE_IO tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TRI_SURFACE_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the TRI_SURFACE_IO library.\n";

  test01 ( "sphere_nodes.txt", "sphere_triangles.txt" );
  test02 ( "sphere_nodes.txt", "sphere_triangles.txt" );
  test03 ( "cube_nodes.txt", "cube_triangles.txt" );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRI_SURFACE_IO_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( string node_file_name, string triangle_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests TRI_SURFACE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2008
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  int node_num;
  int order_num;
  int triangle_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  TRI_SURFACE_SIZE determines the size of various objects\n";
  cout << "  in a TRI_SURFACE file.\n";

  tri_surface_size ( node_file_name, triangle_file_name, &dim_num, &node_num, 
    &order_num, &triangle_num );

  tri_surface_size_print ( node_file_name, triangle_file_name, dim_num, 
    node_num, order_num, triangle_num );

  return;
}
//****************************************************************************80

void test02 ( string node_file_name, string triangle_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests TRI_SURFACE_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2008
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  int node_num;
  double *node_xyz;
  int order_num;
  int triangle_num;
  int *triangle_node;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  TRI_SURFACE_READ reads data from a TRI_SURFACE file.\n";

  tri_surface_size ( node_file_name, triangle_file_name, &dim_num, &node_num, 
    &order_num, &triangle_num );

  tri_surface_size_print ( node_file_name, triangle_file_name, dim_num, 
    node_num, order_num, triangle_num );

  tri_surface_read ( node_file_name, triangle_file_name, dim_num, node_num, 
    order_num, triangle_num, &node_xyz, &triangle_node );

  tri_surface_print ( node_file_name, triangle_file_name, dim_num, node_num, 
    order_num, triangle_num, node_xyz, triangle_node );

  delete [] node_xyz;
  delete [] triangle_node;

  return;
}
//****************************************************************************80

void test03 ( string node_file_name, string triangle_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests TRI_SURFACE_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define NODE_NUM 8
# define ORDER_NUM 3
# define TRIANGLE_NUM 12

  int dim_num = DIM_NUM;
  int node_num = NODE_NUM;
  int order_num = 3;
  int triangle_num = 12;

  double node_xyz[DIM_NUM*NODE_NUM] = {
    0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    1.0, 1.0, 0.0, 
    0.0, 0.0, 1.0, 
    1.0, 0.0, 1.0, 
    0.0, 1.0, 1.0, 
    1.0, 1.0, 1.0 };
  int triangle_node[ORDER_NUM*TRIANGLE_NUM] = {
    1, 3, 2, 
    2, 3, 4, 
    1, 6, 5, 
    1, 2, 6, 
    3, 7, 4, 
    4, 7, 8, 
    5, 6, 8, 
    5, 8, 7, 
    1, 5, 7, 
    1, 7, 3, 
    2, 4, 6, 
    6, 4, 8 };

  cout << "\n";
  cout << "TEST03\n";
  cout << "  TRI_SURFACE_WRITE writes TRI_SURFACE data to two files.\n";

  tri_surface_write ( node_file_name, triangle_file_name, dim_num, node_num, 
    order_num, triangle_num, node_xyz, triangle_node );

  cout << "\n";
  cout << "  Graphics data was written to:\n";
  cout << "    Node file:     \"" << node_file_name     << "\".\n";
  cout << "    Triangle file: \"" << triangle_file_name << "\".\n";

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef ORDER_NUM
# undef TRIANGLE_NUM
}

