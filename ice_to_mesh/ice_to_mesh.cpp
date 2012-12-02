# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstring>
# include <fstream>

using namespace std;

# include "netcdfcpp.h"

int main ( int argc, char **argv );
void data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void data_read ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] );
void mesh_write ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons );
void size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char **argv )

//****************************************************************************80
//
//  Purpose:
//
//    ICE_TO_MESH reads ICE data from a NETCDF file and writes to a MESH file.
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
//  Parameters:
//
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  string filename_mesh;
  string filename_nc;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  string prefix;
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

  timestamp ( );
  cout << "\n";
  cout << "ICE_TO_MESH:\n";
  cout << "  C++ version\n";
  cout << "  Read ICE data from NETCDF file, write to MESH file.\n";
//
//  Check the input argument.
//
  if ( 2 <= argc )
  {
    prefix = argv[1];
  }
  else
  {
    cout << "\n";
    cout << "  Enter the filename prefix:\n";
    cin >> prefix;
  }
//
//  Create the file names.
//
  filename_nc = prefix + ".nc";
  filename_mesh = prefix + ".mesh";
//
//  Read sizes;
//
  size_read ( filename_nc, &dim, &vertices, &edges, &triangles,
    &quadrilaterals, &tetrahedrons, &hexahedrons );
//
//  Print sizes.
//
  size_print ( dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons );
//
//  Allocate memory.
//
  vertex_coordinate = new double [ dim * vertices ];
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
//  Read the data.
//
  cout << "\n";
  cout << "  Reading \"" << filename_nc << "\".\n";

  data_read ( filename_nc, dim, vertices, edges, triangles,
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
    vertex_label, edge_vertex, edge_label, triangle_vertex,
    triangle_label, quadrilateral_vertex, quadrilateral_label,
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
    hexahedron_label );
//
//  Print the data.
//
  if ( vertices < 250 )
  {
    data_print ( dim, vertices, edges, triangles, quadrilaterals,
      tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,
      edge_vertex, edge_label, triangle_vertex, triangle_label,
      quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex,
      tetrahedron_label, hexahedron_vertex, hexahedron_label );
  }
//
//  Write the data.
//
  cout << "\n";
  cout << "  Writing \"" << filename_mesh << "\".\n";

  mesh_write ( filename_mesh, dim, vertices, edges, triangles,
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
    vertex_label, edge_vertex, edge_label,  triangle_vertex, triangle_label,
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex,
    tetrahedron_label, hexahedron_vertex, hexahedron_label );

  cout << "  Conversion completed.\n";
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
//
//  Terminate.
//
  cout << "\n";
  cout << "ICE_TO_MESH:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_PRINT prints the data of a MESH dataset.
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
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
//    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Input, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
  int i;
  int j;

  cout << "\n";
  cout << "  Vertices:\n";
  cout << "\n";
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < dim; i++ )
    {
      cout << "  " << vertex_coordinate[i+j*dim];
    }
    cout << "  (" << vertex_label[j] << ")\n";
  }

  if ( 0 < edges )
  {
    cout << "\n";
    cout << "  Edges:\n";
    cout << "\n";
    for ( j = 0; j < edges; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        cout << "  " << edge_vertex[i+j*2];
    }
    cout << "  (" << edge_label[j] << ")\n";
    }
  }

  if ( 0 < triangles )
  {
    cout << "\n";
    cout << "  Triangles:\n";
    cout << "\n";
    for ( j = 0; j < triangles; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        cout << "  " << triangle_vertex[i+j*3];
      }
      cout << "  (" << triangle_label[j] << ")\n";
    }
  }

  if ( 0 < quadrilaterals )
  {
    cout << "\n";
    cout << "  Quadrilaterals:\n";
    cout << "\n";
    for ( j = 0; j < quadrilaterals; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        cout << "  " << quadrilateral_vertex[i+j*4];
      }
      cout << "  (" << quadrilateral_label[j] << ")\n";
    }
  }

  if ( 0 < tetrahedrons )
  {
    cout << "\n";
    cout << "  Tetrahedrons:\n";
    cout << "\n";
    for ( j = 0; j < tetrahedrons; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        cout << "  " << tetrahedron_vertex[i+j*4];
      }
      cout << "  (" << tetrahedron_label[j] << ")\n";
    }
  }

  if ( 0 < hexahedrons )
  {
    cout << "\n";
    cout << "  Hexahedrons:\n";
    cout << "\n";
    for ( j = 0; j < hexahedrons; j++ )
    {
      for ( i = 0; i < 8; i++ )
      {
        cout << "  " << hexahedron_vertex[i+j*8];
      }
      cout << "  (" << hexahedron_label[j] << ")\n";
    }
  }
  return;
}
//****************************************************************************80

void data_read ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_READ reads ICE data from a NETCDF file.
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
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//    Russ Rew,
//    The NetCDF C++ Interface Guide,
//    Unidata Program Center, August 2008.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be read.
//    Ordinarily, the name should include the extension ".nc".
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
//    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Output, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
//
//  Open the file in "read only" mode.
//
  NcFile dataFile ( filename.c_str ( ), NcFile::ReadOnly );

  if ( !dataFile.is_valid ( ) )
  {
    cout << "\n";
    cout << "DATA_READ: Fatal error!\n";
    cout << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Vertices.
//
  NcVar *var_vertex_coordinate = dataFile.get_var ( "Vertex_Coordinate" );
  var_vertex_coordinate->get ( &vertex_coordinate[0], dim, vertices );

  NcVar *var_vertex_label = dataFile.get_var ( "Vertex_Label" );
  var_vertex_label->get ( &vertex_label[0], vertices );
//
//  Edges.
//
  if ( 0 < edges )
  {
    NcVar *var_edge_vertex = dataFile.get_var ( "Edge_Vertex" );
    var_edge_vertex->get ( &edge_vertex[0], 2, edges );

    NcVar *var_edge_label = dataFile.get_var ( "Edge_Label" );
    var_edge_label->get ( &edge_label[0], edges );
  }
//
//  Triangles.
//
  if ( 0 < triangles )
  {
    NcVar *var_triangle_vertex = dataFile.get_var ( "Triangle_Vertex" );
    var_triangle_vertex->get ( &triangle_vertex[0], 3, triangles );

    NcVar *var_triangle_label = dataFile.get_var ( "Triangle_Label" );
    var_triangle_label->get ( &triangle_label[0], triangles );
  }
//
//  Quadrilaterals.
//
  if ( 0 < quadrilaterals )
  {
    NcVar *var_quadrilateral_vertex = dataFile.get_var ( "Quadrilateral_Vertex" );
    var_quadrilateral_vertex->get ( &quadrilateral_vertex[0], 4, quadrilaterals );

    NcVar *var_quadrilateral_label = dataFile.get_var ( "Quadrilateral_Label" );
    var_quadrilateral_label->get ( &quadrilateral_label[0], quadrilaterals );
  }
//
//  Tetrahedrons.
//
  if ( 0 < tetrahedrons )
  {
    NcVar *var_tetrahedron_vertex = dataFile.get_var ( "Tetrahedron_Vertex" );
    var_tetrahedron_vertex->get ( &tetrahedron_vertex[0], 4, tetrahedrons );

    NcVar *var_tetrahedron_label = dataFile.get_var ( "Tetrahedron_Label" );
    var_tetrahedron_label->get ( &tetrahedron_label[0], tetrahedrons );
  }
//
//  Hexahedrons.
//
  if ( 0 < hexahedrons )
  {
    NcVar *var_hexahedron_vertex = dataFile.get_var ( "Hexahedron_Vertex" );
    var_hexahedron_vertex->get ( &hexahedron_vertex[0], 8, hexahedrons );

    NcVar *var_hexahedron_label = dataFile.get_var ( "Hexahedron_Label" );
    var_hexahedron_label->get ( &hexahedron_label[0], hexahedrons );
  }
//
//  Close the file.
//
  dataFile.close ( );

  return;
}
//****************************************************************************80

void mesh_write ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_WRITE writes mesh data to a MESH file.
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
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be created.
//    Ordinarily, the name should include the extension ".mesh".
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
//    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Input, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
  int i;
  int j;
  ofstream output;

  output.open ( filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "MESH_WRITE - Fatal error!\n";
    cerr << "  Unable to open output file.\n";
    exit ( 1 );
  }

  output << "MeshVersionFormatted 1\n";
  output << "#  Created by mesh_write.C\n";
//
//  Vertices.
//
  output << "\n";
  output << "Vertices\n";
  output << vertices << "\n";
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < dim; i++ )
    {
      output << "  " << vertex_coordinate[i+j*dim];
    }
    output << "  " << vertex_label[j] << "\n";
  }
//
//  Edges.
//
  if ( 0 < edges )
  {
    output << "\n";
    output << "Edges\n";
    output << edges << "\n";
    for ( j = 0; j < edges; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        output << "  " << edge_vertex[i+j*2];
    }
    output << "  " << edge_label[j] << "\n";
    }
  }
//
//  Triangles.
//
  if ( 0 < triangles )
  {
    output << "\n";
    output << "Triangles\n";
    output << triangles << "\n";
    for ( j = 0; j < triangles; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        output << "  " << triangle_vertex[i+j*3];
      }
      output << "  " << triangle_label[j] << "\n";
    }
  }
//
//  Quadrilaterals.
//
  if ( 0 < quadrilaterals )
  {
    output << "\n";
    output << "Quadrilaterals\n";
    output << quadrilaterals << "\n";
    for ( j = 0; j < quadrilaterals; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        output << "  " << quadrilateral_vertex[i+j*4];
      }
      output << "  " << quadrilateral_label[j] << "\n";
    }
  }
//
//  Tetrahedra.
//
  if ( 0 < tetrahedrons )
  {
    output << "\n";
    output << "Tetrahedra\n";
    output << tetrahedrons << "\n";
    for ( j = 0; j < tetrahedrons; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        output << "  " << tetrahedron_vertex[i+j*4];
      }
      output << "  " << tetrahedron_label[j] << "\n";
    }
  }
//
//  Hexahedra.
//
  if ( 0 < hexahedrons )
  {
    output << "\n";
    output << "Hexahedra\n";
    output << hexahedrons << "\n";
    for ( j = 0; j < hexahedrons; j++ )
    {
      for ( i = 0; i < 8; i++ )
      {
        output << "  " << hexahedron_vertex[i+j*8];
      }
      output << "  " << hexahedron_label[j] << "\n";
    }
  }
//
//  End
//
  output << "\n";
  output << "End\n";

  output.close ( );

  return;
}
//****************************************************************************80

void size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    SIZE_PRINT prints the sizes of an ICE dataset.
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
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
{
  cout << "\n";
  cout << "  Number of dimensions = " << dim << "\n";
  cout << "  Number of vertices = " << vertices << "\n";
  cout << "  Number of edges = " << edges << "\n";
  cout << "  Number of triangles = " << triangles << "\n";
  cout << "  Number of quadrilaterals = " << quadrilaterals << "\n";
  cout << "  Number of tetrahedrons = " << tetrahedrons << "\n";
  cout << "  Number of hexahedrons = " << hexahedrons << "\n";

  return;
}
//*****************************************************************************80

void size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

//*****************************************************************************80
//
//  Purpose:
//
//    SIZE_READ reads ICE sizes from a NETCDF file.
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
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be read.
//    Ordinarily, the name should include the extension ".nc".
//
//    Output, int *DIM, the spatial dimension, which should be 2 or 3.
//
//    Output, int *VERTICES, the number of vertices.
//
//    Output, int *EDGES, the number of edges (may be 0).
//
//    Output, int *TRIANGLES, the number of triangles (may be 0).
//
//    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
//
{
  NcDim *dim_dimension;
  NcDim *dim_edges;
  NcDim *dim_eight;
  NcDim *dim_four;
  NcDim *dim_hexahedrons;
  NcToken dim_name;
  int dim_num;
  NcDim *dim_quadrilaterals;
  NcDim *dim_tetrahedrons;
  NcDim *dim_three;
  NcDim *dim_triangles;
  NcDim *dim_two;
  NcDim *dim_vertices;
  NcDim *dim_pointer;
  int i;
//
//  Initialize everything to nothing.
//
  *dim = 0;
  *vertices = 0;
  *edges = 0;
  *triangles = 0;
  *quadrilaterals = 0;
  *tetrahedrons = 0;
  *hexahedrons = 0;
//
//  Open the file in "read only" mode.
//
  NcFile dataFile ( filename.c_str ( ), NcFile::ReadOnly );

  if ( !dataFile.is_valid ( ) )
  {
    cout << "\n";
    cout << "SIZE_READ: Fatal error!\n";
    cout << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Get the dimension information.
//
//  I would much prefer to write "0" as the size of certain dimensions, but I am not
//  allowed to, so I simply omit them from the file.
//
//  Therefore, when I open the file and try to determine dimensions, some dimensions
//  are "missing", which I would have presumed I could discover painlessly by asking
//  for pointers to them, and getting NULLs back.  But that doesn't seem to work either.
//
//  So my bonehead backup is simply to read all the dimensions by index, retrieve
//  their names, and see what I find.
//
  dim_num = dataFile.num_dims ( );

  for ( i = 0; i < dim_num; i++ )
  {
    dim_pointer = dataFile.get_dim ( i );
    dim_name = dim_pointer->name ( );

    if ( !strcmp ( dim_name, "Dimension" ) )
    {
      *dim = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Vertices" ) )
    {
      *vertices = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Edges" ) )
    {
      *edges = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Triangles" ) )
    {
      *triangles = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Quadrilaterals" ) )
    {
      *quadrilaterals = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Tetrahedrons" ) )
    {
      *tetrahedrons = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Hexahedrons" ) )
    {
      *hexahedrons = dim_pointer->size ( );
    }
    else
    {
      cout << "  Ignoring information about dimension \"" << dim_name << "\".\n";
    }
  }
//
//  Close the file.
//
  dataFile.close ( );

  return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
