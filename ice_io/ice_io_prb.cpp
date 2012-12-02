# include <cstdlib>
# include <iostream>
# include <cmath>
# include <cstring>

using namespace std;

# include "ice_io.h"
# include "netcdf.hpp"

int main ( void );
void test01 (
  void ice_size ( int *dim, int *vertices, int *edges, int *triangles,
    int *quadrilaterals, int *tetrahedrons, int *hexahedrons ),
  void ice_data ( int dim, int vertices, int edges, int triangles,
    int quadrilaterals, int tetrahedrons, int hexahedrons,
    double vertex_coordinate[], int vertex_label[], int edge_vertex[],
    int edge_label[], int triangle_vertex[], int triangle_label[],
    int quadrilateral_vertex[], int quadrilateral_label[], int tetrahedron_vertex[],
    int tetrahedron_label[], int hexahedron_vertex[], int hexahedron_label[] ),
  string filename );
void test02 ( string filename );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    ICE_IO_PRB tests the ICE_IO library.
//
//  Discussion:
//
//    We begin by creating a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
{
  string filename_hexa = "hexahexa_2x2x2.nc";
  string filename_cyl = "cyl248.nc";

  timestamp ( );
  cout << "\n";
  cout << "ICE_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the ICE_IO library.\n";
//
//  Create "hexahexa_2x2x2.nc"
//
  test01 ( hexahexa_2x2x2_size, hexahexa_2x2x2_data, filename_hexa );
//
//  Read "hexahexa_2x2x2.nc"
//
  test02 ( filename_hexa );
//
//  Create "cyl248.nc"
//
  test01 ( cyl248_size, cyl248_data, filename_cyl );
//
//  Read "cyl248.nc"
//
  test02 ( filename_cyl );
//
//  Terminate.
//
  cout << "\n";
  cout << "ICE_IO_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 (
  void ice_size ( int *dim, int *vertices, int *edges,
    int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons ),
  void ice_data ( int dim, int vertices, int edges, int triangles,
    int quadrilaterals, int tetrahedrons, int hexahedrons,
    double vertex_coordinate[], int vertex_label[], int edge_vertex[],
    int edge_label[], int triangle_vertex[], int triangle_label[],
    int quadrilateral_vertex[], int quadrilateral_label[],
    int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
    int hexahedron_label[] ),
  string filename )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 writes an ICE grid dataset to a NETCDF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, void ICE_SIZE(), a function which determines sizes.
//
//    Input, void ICE_DATA(), a function which determines data.
//
//    Input, string FILENAME, the name of the file to be created.
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
  cout << "TEST01:\n";
  cout << "  Create an ICE grid dataset, print it,\n";
  cout << "  and write it to an NETCDF file.\n";
//
//  Get sizes.
//
  ice_size ( &dim, &vertices, &edges, &triangles,
    &quadrilaterals, &tetrahedrons, &hexahedrons );
//
//  Print sizes;
//
  size_print ( dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons );
//
//  Allocate memory.
//
  vertex_coordinate = new double [ 3 * vertices ];
  vertex_label = new int [ vertices ];
  edge_vertex = new int [ 2 * edges ];
  edge_label = new int [ edges ];
  triangle_vertex = new int [ 3 * triangles ];
  triangle_label = new int [ triangles ];
  quadrilateral_vertex = new int [ 4 * quadrilaterals ];
  quadrilateral_label = new int [ quadrilaterals ];
  tetrahedron_vertex = new int [ 4 * tetrahedrons ];
  tetrahedron_label = new int [ tetrahedrons ];
  hexahedron_vertex = new int [ 8 * hexahedrons ];
  hexahedron_label = new int [ hexahedrons ];
//
//  Get data.
//
  ice_data ( dim, vertices, edges, triangles, quadrilaterals, tetrahedrons,
    hexahedrons, vertex_coordinate, vertex_label, edge_vertex, edge_label,
    triangle_vertex, triangle_label, quadrilateral_vertex, quadrilateral_label,
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, hexahedron_label );
//
//  Print the data.
//
  cout << "\n";
  cout << "  Data to be written to \"" << filename << "\":\n";

  data_print ( dim, vertices, edges, triangles, quadrilaterals, tetrahedrons,
    hexahedrons, vertex_coordinate, vertex_label, edge_vertex, edge_label,
    triangle_vertex, triangle_label, quadrilateral_vertex, quadrilateral_label,
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, hexahedron_label );
//
//  Create the file.
//
  ice_write ( filename, dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex,
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex,
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
    hexahedron_vertex, hexahedron_label );

  cout << "\n";
  cout << "  Created the file \"" << filename << "\"\n";
//
//  Free memory.
//
  delete [] vertex_coordinate;
  delete [] vertex_label;
  delete [] edge_vertex;
  delete [] edge_label;
  delete [] triangle_vertex;
  delete [] triangle_label;
  delete [] quadrilateral_vertex;
  delete [] quadrilateral_label;
  delete [] tetrahedron_vertex;
  delete [] tetrahedron_label;
  delete [] hexahedron_vertex;
  delete [] hexahedron_label;

  return;
}
//****************************************************************************80

void test02 ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 reads an ICE grid dataset from a NETCDF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be read.
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
  cout << "TEST02:\n";
  cout << "  Read an ICE grid dataset from a NETCDF file,\n";
  cout << "  and print the data.\n";
//
//  Read sizes;
//
  size_read ( filename, &dim, &vertices, &edges, &triangles, &quadrilaterals,
    &tetrahedrons, &hexahedrons );
//
//  Print sizes;
//
  size_print ( dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons );
//
//  Allocate memory.
//
  vertex_coordinate = new double [ 3 * vertices ];
  vertex_label = new int [ vertices ];
  edge_vertex = new int [ 2 * edges ];
  edge_label = new int [ edges ];
  triangle_vertex = new int [ 3 * triangles ];
  triangle_label = new int [ triangles ];
  quadrilateral_vertex = new int [ 4 * quadrilaterals ];
  quadrilateral_label = new int [ quadrilaterals ];
  tetrahedron_vertex = new int [ 4 * tetrahedrons ];
  tetrahedron_label = new int [ tetrahedrons ];
  hexahedron_vertex = new int [ 8 * hexahedrons ];
  hexahedron_label = new int [ hexahedrons ];
//
//  Read the file
//
  data_read ( filename, dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex,
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex,
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
    hexahedron_vertex, hexahedron_label );
//
//  Print the data.
//
  cout << "\n";
  cout << "  Data from file \"" << filename << "\"\n";

  data_print ( dim, vertices, edges, triangles, quadrilaterals, tetrahedrons,
    hexahedrons, vertex_coordinate, vertex_label, edge_vertex, edge_label,
    triangle_vertex, triangle_label, quadrilateral_vertex, quadrilateral_label,
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, hexahedron_label );
//
//  Free memory.
//
  delete [] vertex_coordinate;
  delete [] vertex_label;
  delete [] edge_vertex;
  delete [] edge_label;
  delete [] triangle_vertex;
  delete [] triangle_label;
  delete [] quadrilateral_vertex;
  delete [] quadrilateral_label;
  delete [] tetrahedron_vertex;
  delete [] tetrahedron_label;
  delete [] hexahedron_vertex;
  delete [] hexahedron_label;

  return;
}
