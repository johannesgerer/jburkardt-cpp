# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "medit_io.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( string filename );
void test04 ( string filename );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MEDIT_IO_PRB.
//
//  Discussion:
//
//    MEDIT_IO_PRB tests the MEDIT_IO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  string filename;

  timestamp ( );
  cout << "\n";
  cout << "MEDIT_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the MEDIT_IO library.\n";
//
//  Create the file hexahexa_2x2x2.mesh
//
  test01 ( );
//
//  Read and print the sizes of file hexahexa_2x2x2.mesh.
//
  filename = "hexahexa_2x2x2.mesh";
  test03 ( filename );
//
//  Create the file cyl248.mesh
//
  test02 ( );
//
//  Read and print the sizes of file cyl248.mesh.
//
  filename = "cyl248.mesh";
  test03 ( filename );
//
//  Read and print the data in file cyl248.mesh.
//
  filename = "cyl248.mesh";
  test04 ( filename );
//
//  Terminate.
//
  cout << "\n";
  cout << "MEDIT_IO_PRB\n";
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
//    TEST01 creates a MESH dataset and writes it to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  string filename;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int *quadrilateral_label;
  int *quadrilateral_vertex;
  int quadrilaterals;
  int *tetrahedron_label;
  int *tetrahedron_vertex;
  int tetrahedrons;
  int *triangle_label;
  int *triangle_vertex;
  int triangles;
  double *vertex_coordinate;
  int *vertex_label;
  int vertices;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Create a hexahedral mesh and write it to a file.\n";
//
//  Get sizes.
//
  hexahexa_2x2x2_size ( &dim, &vertices, &edges, &triangles, &quadrilaterals, 
    &tetrahedrons, &hexahedrons );
//
//  Allocate memory.
//
  edge_label = new int[edges];
  edge_vertex = new int[2*edges];
  hexahedron_label = new int[hexahedrons];
  hexahedron_vertex = new int[8*hexahedrons];
  quadrilateral_label = new int[quadrilaterals];
  quadrilateral_vertex = new int[4*quadrilaterals];
  tetrahedron_label = new int[tetrahedrons];
  tetrahedron_vertex = new int[4*tetrahedrons];
  triangle_label = new int[triangles];
  triangle_vertex = new int[3*triangles];
  vertex_coordinate = new double[dim*vertices];
  vertex_label = new int[vertices];
//
//  Get the data.
//
  hexahexa_2x2x2_data ( dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, 
    edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
//
//  Write the data.
//
  filename = "hexahexa_2x2x2.mesh";

  mesh_write ( filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );

  cout << "\n";
  cout << "  Created the file \"" << filename << "\".\n";
//
//  Deallocate memory.
//
  delete [] edge_label;
  delete [] edge_vertex;
  delete [] hexahedron_label;
  delete [] hexahedron_vertex;
  delete [] quadrilateral_label;
  delete [] quadrilateral_vertex;
  delete [] tetrahedron_label;
  delete [] tetrahedron_vertex;
  delete [] triangle_label;
  delete [] triangle_vertex;
  delete [] vertex_coordinate;
  delete [] vertex_label;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 creates a MESH dataset and writes it to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  string filename;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int *quadrilateral_label;
  int *quadrilateral_vertex;
  int quadrilaterals;
  int *tetrahedron_label;
  int *tetrahedron_vertex;
  int tetrahedrons;
  int *triangle_label;
  int *triangle_vertex;
  int triangles;
  double *vertex_coordinate;
  int *vertex_label;
  int vertices;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Create a tetrahedral mesh and write it to a file.\n";
//
//  Get sizes.
//
  cyl248_size ( &dim, &vertices, &edges, &triangles, &quadrilaterals, 
    &tetrahedrons, &hexahedrons );
//
//  Allocate memory.
//
  edge_label = new int[edges];
  edge_vertex = new int[2*edges];
  hexahedron_label = new int[hexahedrons];
  hexahedron_vertex = new int[8*hexahedrons];
  quadrilateral_label = new int[quadrilaterals];
  quadrilateral_vertex = new int[4*quadrilaterals];
  tetrahedron_label = new int[tetrahedrons];
  tetrahedron_vertex = new int[4*tetrahedrons];
  triangle_label = new int[triangles];
  triangle_vertex = new int[3*triangles];
  vertex_coordinate = new double[dim*vertices];
  vertex_label = new int[vertices];
//
//  Get the data.
//
  cyl248_data ( dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, 
    edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
//
//  Write the data.
//
  filename = "cyl248.mesh";

  mesh_write ( filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );

  cout << "\n";
  cout << "  Created the file \"" << filename << "\".\n";
//
//  Deallocate memory.
//
  delete [] edge_label;
  delete [] edge_vertex;
  delete [] hexahedron_label;
  delete [] hexahedron_vertex;
  delete [] quadrilateral_label;
  delete [] quadrilateral_vertex;
  delete [] tetrahedron_label;
  delete [] tetrahedron_vertex;
  delete [] triangle_label;
  delete [] triangle_vertex;
  delete [] vertex_coordinate;
  delete [] vertex_label;

  return;
}
//****************************************************************************80

void test03 ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_IO_TEST03 reads and prints the sizes in a MESH dataset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int edges;
  int hexahedrons;
  int quadrilaterals;
  int tetrahedrons;
  int triangles;
  int vertices;

  cout << "\n";
  cout << "MESH_IO_TEST03\n";
  cout << "  Read a mesh file and print its sizes.\n";
//
//  Read sizes.
//
  mesh_size_read ( filename, &dim, &vertices, &edges, &triangles, 
    &quadrilaterals, &tetrahedrons, &hexahedrons );
//
//  Print sizes.
//
  cout << "\n";
  cout << "  Header information for \"" << filename << "\"\n";

  mesh_size_print ( dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons );

  return;
}
//****************************************************************************80

void test04 ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_IO_TEST04 reads a MESH dataset and prints its data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int *quadrilateral_label;
  int *quadrilateral_vertex;
  int quadrilaterals;
  int *tetrahedron_label;
  int *tetrahedron_vertex;
  int tetrahedrons;
  int *triangle_label;
  int *triangle_vertex;
  int triangles;
  double *vertex_coordinate;
  int *vertex_label;
  int vertices;

  cout << "\n";
  cout << "MESH_IO_TEST04\n";
  cout << "  Read a mesh file and print its data.\n";
//
//  Read sizes.
//
  mesh_size_read ( filename, &dim, &vertices, &edges, &triangles, 
    &quadrilaterals, &tetrahedrons, &hexahedrons );
//
//  Allocate memory.
//
  edge_label = new int[edges];
  edge_vertex = new int[2*edges];
  hexahedron_label = new int[hexahedrons];
  hexahedron_vertex = new int[8*hexahedrons];
  quadrilateral_label = new int[quadrilaterals];
  quadrilateral_vertex = new int[4*quadrilaterals];
  tetrahedron_label = new int[tetrahedrons];
  tetrahedron_vertex = new int[4*tetrahedrons];
  triangle_label = new int[triangles];
  triangle_vertex = new int[3*triangles];
  vertex_coordinate = new double[dim*vertices];
  vertex_label = new int[vertices];
//
//  Read the data.
//
  mesh_data_read ( filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
//
//  Print the data.
//
  cout << "\n";
  cout << "  Data for file \"" << filename << "\".\n";

  mesh_data_print ( dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,  edge_vertex, 
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, 
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, 
    hexahedron_vertex, hexahedron_label );
//
//  Deallocate memory.
//
  delete [] edge_label;
  delete [] edge_vertex;
  delete [] hexahedron_label;
  delete [] hexahedron_vertex;
  delete [] quadrilateral_label;
  delete [] quadrilateral_vertex;
  delete [] tetrahedron_label;
  delete [] tetrahedron_vertex;
  delete [] triangle_label;
  delete [] triangle_vertex;
  delete [] vertex_coordinate;
  delete [] vertex_label;

  return;
}
