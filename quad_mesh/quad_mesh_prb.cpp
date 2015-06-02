# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "quad_mesh.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test062 ( );
void test07 ( );
void test072 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test105 ( );
void test11 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUAD_MESH_PRB.
//
//  Discussion:
//
//    QUAD_MESH_PRB tests the QUAD_MESH library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "QUAD_MESH_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the QUAD_MESH library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test062 ( );
  test07 ( );
  test072 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  test105 ( );
  test11 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QUAD_MESH_PRB\n";
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
//    TEST01 tests ADJ_SIZE_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int adj_num;
  int *adj_col;
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  int k;
  int node;
  int node_num;
  double *node_xy;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  ADJ_SIZE counts the node adjacencies.\n";
//
//  Get the sizes of the example.
//
  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );
//
//  Get the count of the node adjacencies.
//
  adj_col = new int[node_num+1];

  adj_num = adj_size_q4_mesh ( node_num, element_num, element_node,
    element_neighbor, adj_col );

  cout << "\n";
  cout << "  Number of adjacency entries is " << adj_num << "\n";

  cout << "\n";
  cout << "  Adjacency pointers:\n";
  cout << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(8) << node
         << "  " << setw(8) << adj_col[node]
         << "  " << setw(8) << adj_col[node+1] - 1 << "\n";
  }

  delete [] adj_col;
  delete [] node_xy;
  delete [] element_neighbor;
  delete [] element_node;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests AREA_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double *element_area;
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  double mesh_area;
  int node_num;
  double *node_xy;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  AREA_Q4_MESH computes the area of each element\n";
  cout << "  in a Q4 mesh.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  element_area = new double[element_num];

  area_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_area, &mesh_area );

  r8vec_print ( element_num, element_area, "  Element areas:" );

  cout << "\n";
  cout << "   Mesh =   " << mesh_area << "\n";

  delete [] element_area;
  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 demonstrates AREA_QUAD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double quad_xy[2*4] = {
    1.0, 2.0,
    5.0, 2.0,
    5.0, 3.0,
    4.0, 4.0 };

  cout << "\n";
  cout << "TEST03\n";
  cout << "  AREA_QUAD computes the area of a quadrilateral.\n";

  area = area_quad ( quad_xy );

  cout << "\n";
  cout << "  Area = " << area << "\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests BOUNDARY_EDGE_COUNT_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int boundary_edge_num;
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  int node_num;
  double *node_xy;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  BOUNDARY_EDGE_COUNT_Q4_MESH counts the\n";
  cout << "    boundary edges.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  boundary_edge_num = boundary_edge_count_q4_mesh ( element_num, element_node );

  cout << "\n";
  cout << "  Number of boundary edges = " << boundary_edge_num << "\n";
  cout << "  Correct number =           " << 22 << "\n";

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests BOUNDARY_EDGE_COUNT_EULER_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int boundary_edge_num;
  int element_num;
  int hole_num;
  int node_num;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  BOUNDARY_EDGE_COUNT_EULER_Q4_MESH counts the\n";
  cout << "    boundary edges using Euler's formula.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  boundary_edge_num = boundary_edge_count_euler_q4_mesh ( node_num,
    element_num, hole_num );

  cout << "\n";
  cout << "  Number of boundary edges = " << boundary_edge_num << "\n";
  cout << "  Correct number =           " << 22 << "\n";

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests EXAMPLE1_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  int node_num;
  double *node_xy;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  cout << "\n";
  cout << "  Number of nodes =    " << node_num << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Number of holes =    " << hole_num << "\n";

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  r8mat_transpose_print ( 2, node_num, node_xy, "  Node coordinates:" );
  i4mat_transpose_print ( 4, element_num, element_node, "  Elements:" );
  i4mat_transpose_print ( 4, element_num, element_neighbor,
    "  Element neighbors" );

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST062 tests EXAMPLE1_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_node;
  int element_num;
  int element_show;
  int hole_num;
  int node_num;
  int node_show;
  double *node_xy;
  string output_filename;

  cout << "\n";
  cout << "TEST062\n";
  cout << "  EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.\n";
  cout << "  PLOT_Q4_MESH plots it.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  cout << "\n";
  cout << "  Number of nodes =    " << node_num << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Number of holes =    " << hole_num << "\n";

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  node_show = 2;
  element_show = 2;
  output_filename = "q4_mesh_ex1.eps";

  plot_q4_mesh ( node_num, element_num, node_xy, element_node,
    node_show, element_show, output_filename );

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests EXAMPLE2_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  int node_num;
  double *node_xy;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.\n";

  example2_q4_mesh_size ( &node_num, &element_num, &hole_num );

  cout << "\n";
  cout << "  Number of nodes =    " << node_num << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Number of holes =    " << hole_num << "\n";

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example2_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  r8mat_transpose_print ( 2, node_num, node_xy, "  Node coordinates:" );
  i4mat_transpose_print ( 4, element_num, element_node, "  Elements:" );
  i4mat_transpose_print ( 4, element_num, element_neighbor,
    "  Element neighbors" );

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test072 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST072 tests EXAMPLE2_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_node;
  int element_num;
  int element_show;
  int hole_num;
  int node_num;
  int node_show;
  double *node_xy;
  string output_filename;

  cout << "\n";
  cout << "TEST072\n";
  cout << "  EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.\n";
  cout << "  PLOT_Q4_MESH plots it.\n";

  example2_q4_mesh_size ( &node_num, &element_num, &hole_num );

  cout << "\n";
  cout << "  Number of nodes =    " << node_num << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Number of holes =    " << hole_num << "\n";

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example2_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  node_show = 2;
  element_show = 2;
  output_filename = "q4_mesh_ex2.eps";

  plot_q4_mesh ( node_num, element_num, node_xy, element_node,
    node_show, element_show, output_filename );

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests NEIGHBOR_ELEMENTS_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_neighbor2;
  int *element_node;
  int element_num;
  int hole_num;
  int node_num;
  double *node_xy;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  NEIGHBOR_ELEMENTS_Q4_MESH determines the\n";
  cout << "  adjacency relationships between elements.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  i4mat_transpose_print ( 4, element_num, element_neighbor,
    "  Element neighbors as reported by EXAMPLE1_Q4_MESH:" );

  element_neighbor2 = neighbor_elements_q4_mesh ( element_num, element_node );

  i4mat_transpose_print ( 4, element_num, element_neighbor2,
    "  Element neighbors computed by NEIGHBOR_ELEMENTS_Q4_MESH:" );

  delete [] element_neighbor;
  delete [] element_neighbor2;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 writes data to files.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  int node_num;
  double *node_xy;
  string output_filename;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  Write Q4 Mesh Example #1 to files.\n";

  example2_q4_mesh_size ( &node_num, &element_num, &hole_num );

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example2_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  output_filename = "q4_mesh_ex2_element_neighbors.txt";
  i4mat_write ( output_filename, 4, element_num, element_neighbor );
  cout << "\n";
  cout << "  Element neighbors written to \"" << output_filename << "\".\n";

  output_filename = "q4_mesh_ex2_elements.txt";
  i4mat_write ( output_filename, 4, element_num, element_node );
  cout << "  Elements written to \"" << output_filename << "\".\n";

  output_filename = "q4_mesh_ex2_xy.txt";
  r8mat_write ( output_filename, 2, node_num, node_xy );
  cout << "  Node coordinates written to \"" << output_filename << "\".\n";

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests NODE_ORDER_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  int node_num;
  int *node_order;
  double *node_xy;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  NODE_ORDER_4_MESH computes the order\n";
  cout << "  of the nodes in a Q4 mesh.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  node_order = node_order_q4_mesh ( element_num, element_node, node_num );

  i4vec_print ( node_num, node_order, "      NODE         ORDER" );

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_order;
  delete [] node_xy;

  return;
}
//****************************************************************************80

void test105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST105 tests SAMPLE_Q4_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int *element_node;
  int element_num;
  int hole_num;
  int node_num;
  double *node_xy;
  int sample;
  int sample_num;
  int *sample_element;
  double *sample_xy;
  int seed;

  cout << "\n";
  cout << "TEST105\n";
  cout << "  SAMPLE_Q4_MESH returns uniform sample points from\n";
  cout << "  a Q4 mesh.\n";

  example1_q4_mesh_size ( &node_num, &element_num, &hole_num );

  element_neighbor = new int[4*element_num];
  element_node = new int[4*element_num];
  node_xy = new double[2*node_num];

  example1_q4_mesh ( node_num, element_num, node_xy, element_node,
    element_neighbor );

  sample_num = 20;

  sample_xy = new double[2*sample_num];
  sample_element = new int[sample_num];

  seed = 123456789;

  sample_q4_mesh ( node_num, node_xy, element_num, element_node,
    sample_num, &seed, sample_xy, sample_element );

  cout << "\n";
  cout << "             X        Y     Element\n";
  cout << "\n";

  for ( sample = 0; sample < sample_num; sample++ )
  {
    cout << "  " << setw(8) << sample
         << "  " << setw(8) << sample_xy[0+sample*2]
         << "  " << setw(8) << sample_xy[1+sample*2]
         << "  " << setw(8) << sample_element[sample] << "\n";
  }

  delete [] element_neighbor;
  delete [] element_node;
  delete [] node_xy;
  delete [] sample_element;
  delete [] sample_xy;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 demonstrates SAMPLE_QUAD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  string output_filename = "sample_quad.txt";
  double quad_xy[2*4] = {
    1.0, 2.0,
    5.0, 2.0,
    5.0, 3.0,
    4.0, 4.0 };
  int seed;
  double *xy;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  SAMPLE_QUAD computes N random points in a quadrilateral.\n";
  cout << "  Write them to a file.\n";

  n = 5000;

  seed = 123456789;

  xy = sample_quad_new ( quad_xy, n, &seed );

  r8mat_write ( output_filename, 2, n, xy );

  cout << "\n";
  cout << "  Point coordinates written to \"" << output_filename << "\".\n";

  delete [] xy;

  return;
}
