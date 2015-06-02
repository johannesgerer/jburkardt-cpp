# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "gmsh_io.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GMSH_IO_PRB.
//
//  Discussion:
//
//    GMSH_IO_PRB tests the GMSH_IO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "GMSH_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the GMSH_IO library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "GMSH_IO_PRB\n";
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
//    TEST01 gets the example 2D data and writes it to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *element_node;
  int element_num;
  int element_order;
  string gmsh_filename = "example_2d.msh";
  int m;
  int node_num;
  double *node_x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Get example 2D data, write to a file.\n";
//
//  Get sizes.
//
  gmsh_mesh2d_node_size_example ( node_num, m );

  gmsh_mesh2d_element_size_example ( element_num, element_order );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "  Spatial dimension = " << m << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Order of elements = " << element_order << "\n";
//
//  Get the data.
//
  node_x = gmsh_mesh2d_node_data_example ( node_num, m );

  element_node = gmsh_mesh2d_element_data_example ( element_num, element_order );
//
//  Print some of the data.
//
  r8mat_transpose_print_some ( m, node_num, node_x, 
    1, 1, m, 10, "  Coordinates for first 10 nodes:" );

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Node connectivity of first 10 elements:" );
//
//  Write the GMSH file.
//
  gmsh_mesh2d_write ( gmsh_filename, m, node_num, node_x, 
    element_order, element_num, element_node );

  cout << "\n";
  cout << "  Wrote example data to file \"" << gmsh_filename << "\"\n";
//
//  Clean up.
//
  delete [] element_node;
  delete [] node_x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 reads the example data from a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//   John Burkardt
//
{
  int *element_node;
  int element_num;
  int element_order;
  string gmsh_filename = "example_2d.msh";
  int m;
  int node_num;
  double *node_x;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Read data from a file.\n";
//
//  Get the data size.
//
  gmsh_size_read ( gmsh_filename, node_num, m, element_num, 
    element_order );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Node data read from file \"" << gmsh_filename << "\"\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "  Spatial dimension = " << m << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Element order = " << element_order << "\n";
//
//  Allocate memory.
//
  node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
  element_node = ( int * ) 
    malloc ( element_order * element_num * sizeof ( int ) );
//
//  Get the data.
//
  gmsh_data_read ( gmsh_filename, m, node_num, node_x, 
    element_order, element_num, element_node );
//
//  Print some of the data.
//
  r8mat_transpose_print_some ( m, node_num, node_x, 
    1, 1, m, 10, "  Coordinates for first 10 nodes:" );

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Connectivity for first 10 elements:" );
//
//  Clean up.
//
  delete [] element_node;
  delete [] node_x;

  return;
}
