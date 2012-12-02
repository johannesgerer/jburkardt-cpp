# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "triangle_io.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_IO_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE_IO library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_IO_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 gets the example node data and writes it to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *node_att;
  int node_att_num;
  double *node_coord;
  int node_dim;
  string node_file = "example.node";
  int *node_marker;
  int node_marker_num;
  int node_num;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Get example node data, write to a node file.\n";
//
//  Get node example size.
//
  node_size_example ( &node_num, &node_dim, &node_att_num, &node_marker_num );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "  Spatial dimension =" << node_dim << "\n";
  cout << "  Number of node attributes = " << node_att_num << "\n";
  cout << "  Number of node markers = " << node_marker_num << "\n";
//
//  Allocate memory for node data.
//
  node_coord = new double[ node_dim * node_num ];
  node_att = new double[ node_att_num * node_num ];
  node_marker = new int[ node_marker_num * node_num ];
//
//  Get the node data.
//
  node_data_example ( node_num, node_dim, node_att_num, node_marker_num,
    node_coord, node_att, node_marker );
//
//  Print some of the data.
//
  r8mat_transpose_print_some ( node_dim, node_num, node_coord, 
    1, 1, node_dim, 10, "  Coordinates for first 10 nodes:" );

  r8mat_transpose_print_some ( node_att_num, node_num, node_att,
    1, 1, node_att_num, 10, "  Attributes for first 10 nodes:" );

  i4mat_transpose_print_some ( node_marker_num, node_num, node_marker,
    1, 1, node_marker_num, 10, "  Markers for first 10 nodes:" ); 
//
//  Write the node information to node file.
//
  node_write ( node_file, node_num, node_dim, node_att_num, node_marker_num, 
    node_coord, node_att, node_marker );

  cout << "\n";
  cout << "  Node data written to file \"" << node_file << "\"\n";
//
//  Clean up.
//
  delete [] node_att;
  delete [] node_coord;
  delete [] node_marker;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 gets the example element data and writes it to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *element_att;
  int element_att_num;
  string element_file = "example.ele";
  int *element_node;
  int element_num;
  int element_order;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Get example element data, write to an element file.\n";
//
//  Get element example size.
//
  element_size_example ( &element_num, &element_order, &element_att_num );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Order of elements = " << element_order << "\n";
  cout << "  Number of element attributes = " << element_att_num << "\n";
//
//  Allocate memory.
//
  element_node = new int[ element_order * element_num ];
  element_att = new double[ element_att_num * element_num ];
//
//  Get the data.
//
  element_data_example ( element_num, element_order, element_att_num, 
    element_node, element_att );
//
//  Print some of the data.
//
  i4mat_transpose_print_some ( element_order, element_num, element_node,
    1, 1, element_order, 10, "  Node connectivity of first 10 elements:" );

  r8mat_transpose_print_some ( element_att_num, element_num, element_att,
    1, 1, element_att_num, 10, "  Attributes for first 10 elements:" ); 
//
//  Write the node information to node file.
//
  element_write ( element_file, element_num, element_order, element_att_num, 
    element_node, element_att );

  cout << "\n";
  cout << "  Element data written to file \"" << element_file << "\"\n";
//
//  Clean up.
//
  delete [] element_att;
  delete [] element_node;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 reads the example node data from a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *node_att;
  int node_att_num;
  double *node_coord;
  int node_dim;
  string node_file = "example.node";
  int *node_marker;
  int node_marker_num;
  int node_num;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Read node data from a node file.\n";
//
//  Get the data size.
//
  node_size_read ( node_file, &node_num, &node_dim, &node_att_num, 
    &node_marker_num );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Node data read from file \"" << node_file << "\"\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "  Spatial dimension = " << node_dim << "\n";
  cout << "  Number of node attributes = " << node_att_num << "\n";
  cout << "  Number of node markers = " << node_marker_num << "\n";
//
//  Allocate memory.
//
  node_coord = new double[ node_dim * node_num ];
  node_att = new double[ node_att_num * node_num ];
  node_marker = new int[ node_marker_num * node_num ];
//
//  Get the data.
//
  node_data_read ( node_file, node_num, node_dim, node_att_num, node_marker_num,
    node_coord, node_att, node_marker );
//
//  Print some of the data.
//
  r8mat_transpose_print_some ( node_dim, node_num, node_coord, 
    1, 1, node_dim, 10, "  Coordinates for first 10 nodes:" );

  r8mat_transpose_print_some ( node_att_num, node_num, node_att,
    1, 1, node_att_num, 10, "  Attributes for first 10 nodes:" );

  i4mat_transpose_print_some ( node_marker_num, node_num, node_marker,
    1, 1, node_marker_num, 10, "  Markers for first 10 nodes:" ); 
//
//  Clean up.
//
  delete [] node_att;
  delete [] node_coord;
  delete [] node_marker;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 reads the example element data from a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *element_att;
  int element_att_num;
  string element_file = "example.ele";
  int *element_node;
  int element_num;
  int element_order;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Read element data from an element file.\n";
//
//  Get data size.
//
  element_size_read ( element_file, &element_num, &element_order, 
    &element_att_num );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Element data read from file \"" << element_file << "\"\n";
  cout << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Element order = " << element_order << "\n";
  cout << "  Number of element attributes = " << element_att_num << "\n";
//
//  Allocate memory.
//
  element_node = new int[ element_order * element_num ];
  element_att = new double[ element_att_num * element_num ];  
//
//  Get the data.
//
  element_data_read ( element_file, element_num, element_order, element_att_num, 
    element_node, element_att );
//
//  Print some of the data.
//
  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Connectivity for first 10 elements:" );

  r8mat_transpose_print_some ( element_att_num, element_num, element_att,
    1, 1, element_att_num, 10, "  Attributes for first 10 elements:" );
//
//  Clean up.
//
  delete [] element_att;
  delete [] element_node;

  return;
}
