# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void element_data_read ( string element_file, int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] );
void element_size_read ( string element_file, int *element_num, 
  int *element_order, int *element_att_num );
void mesh_write ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void node_data_read ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] );
void node_size_read ( string node_file, int *node_num, int *node_dim, 
  int *node_att_num, int *node_marker_num );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_TO_MEDIT.
//
//  Discussion:
//
//    The TRIANGLE program creates "node" and "element" files that define
//    a triangular mesh.  A typical pair of such files might have the names
//    "suv.node" and "suv.ele".
//
//    This program reads this pair of files and creates a MEDIT mesh file, whose
//    name might be "suv.mesh".
//
//  Usage:
//
//    triangle_to_medit prefix
//
//    where 'prefix' is the common filename prefix so that:
//
//    * prefix.node contains the coordinates of the nodes;
//    * prefix.ele contains the indices of nodes forming each element.
//    * prefix.mesh will be the MESH file created by the program.
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
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  double *element_att;
  int element_att_num;
  string element_filename;
  int *element_node;
  int element_num;
  int element_order;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int i;
  int j;
  string mesh_filename;
  double *node_att;
  int node_att_num;
  string node_filename;
  int *node_marker;
  double *node_coord;
  int node_dim;
  int node_marker_num;
  int node_num;
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
  cout << "TRIANGLE_TO_MEDIT:\n";
  cout << "  C++ version\n";
  cout << "  Read a pair of NODE and ELE files created by TRIANGLE.\n";
  cout << "  Write a corresponding MEDIT mesh file.\n";
//
//  Get the filename prefix.
//
  if ( 1 <= argc )
  {
    prefix = argv[1];
  }
  else
  {
    cout << "\n";
    cout << "  Please enter the filename prefix:\n";

    cin >> prefix;
  }
//
//  Create the file names.
//
  node_filename = prefix + ".node";
  element_filename = prefix + ".ele";
  mesh_filename = prefix + ".mesh";

  cout << "\n";
  cout << "  Read Node file \"" << node_filename << "\"\n";
  cout << "    and Element file \"" << element_filename << "\".\n";
  cout << "  Create mesh file \"" << mesh_filename << "\".\n";
//
//  Read the TRIANGLE NODE data.
//
  node_size_read ( node_filename, &node_num, &node_dim, &node_att_num, 
    &node_marker_num );

  node_coord = new double[2*node_num];
  node_att = new double[node_att_num*node_num];
  node_marker = new int[node_num];

  node_data_read ( node_filename, node_num, node_dim, node_att_num, 
    node_marker_num, node_coord, node_att, node_marker );
//
//  Read the TRIANGLE ELE data.
//
  element_size_read ( element_filename, &element_num, &element_order, 
    &element_att_num );

  element_node = new int[element_order*element_num];
  element_att = new double[element_att_num*element_num];

  element_data_read ( element_filename, element_num, element_order, 
    element_att_num, element_node, element_att );
//
//  Write the MESH data.
//
  dim = 2;
  vertices = node_num;
  edges = 0;
  triangles = element_num;
  quadrilaterals = 0;
  tetrahedrons = 0;
  hexahedrons = 0;
  vertex_coordinate = new double[2*vertices];
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      vertex_coordinate[i+j*2] = node_coord[i+j*2];
    }
  }
  vertex_label = new int[vertices];
  for ( j = 0; j < vertices; j++ )
  {
    vertex_label[j] = node_marker[j];
  }
  edge_vertex = NULL;
  edge_label = NULL;
  triangle_vertex = new int[3*triangles];
  for ( j = 0; j < triangles; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_vertex[i+j*3] = element_node[i+j*3];
    }
  }
  triangle_label = new int[triangles];
  for ( j = 0; j < triangles; j++ )
  {
    triangle_label[j] = 0;
  }
  quadrilateral_vertex = NULL;
  quadrilateral_label = NULL;
  tetrahedron_vertex = NULL;
  tetrahedron_label = NULL;
  hexahedron_vertex = NULL;
  hexahedron_label = NULL;

  mesh_write ( mesh_filename, dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex,
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, 
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
    hexahedron_vertex, hexahedron_label );
//
//  Free memory.
//
  delete [] element_att;
  delete [] element_node;
  delete [] node_att;
  delete [] node_coord;
  delete [] node_marker;
  delete [] triangle_label;
  delete [] triangle_vertex;
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_TO_MEDIT:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void element_data_read ( string element_file, int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] )

//*****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_DATA_READ reads the header information from an element file.
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
//  Parameters:
//
//    Input, string ELEMENT_FILE, the name of the file to be read.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_ATT_NUM, number of element attributes listed on each 
//    node record.
//
//    Output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
//    nodes that make up each element.
//
//    Output, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
//    of each element.
//
{
  int element;
  int i;
  int i1;
  int i2;
  int i3;
  ifstream input;
  int ival;
  double value;

  element = -1;

  input.open ( element_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "ELEMENT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
/*
  Read, but ignore, dimension line.
*/
    if ( element == -1 )
    {
      input >> i1 >> i2 >> i3;
    }
    else
    {
      input >> ival;

      for ( i = 0; i < element_order; i++ )
      {
        input >> ival;
        element_node[i+element*element_order] = ival;
      }
      for ( i = 0; i < element_att_num; i++ )
      {
        input >> value;
        element_att[i+element*element_att_num] = value;
      }
    }

    element = element + 1;

    if ( element_num <= element )
    {
      break;
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void element_size_read ( string element_file, int *element_num, 
  int *element_order, int *element_att_num )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_SIZE_READ reads the header information from an element file.
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
//  Parameters:
//
//    Input, string ELEMENT_FILE, the name of the file to be read.
//
//    Output, int *ELEMENT_NUM, the number of elements.
//
//    Output, int *ELEMENT_ORDER, the order of the elements.
//
//    Output, int *ELEMENT_ATT_NUM, the number of element attributes.
//
{
  ifstream input;

  input.open ( element_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "ELEMENT_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  input >> *element_num >> *element_order >> *element_att_num;

  input.close ( );

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
  output << "#  Created by mesh_write.cpp\n";
//
//  Dimension.
//
  output << "\n";
  output << "Dimension\n";
  output << dim << "\n";
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

void node_data_read ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_HEADER_READ reads the header information from a node file.
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
//  Parameters:
//
//    Input, string NODE_FILE, the name of the node file to be read.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, int NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Input, int NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
//    Output, double NODE_COORD[NODE_DIM*NODE_NUM], the nodal coordinates.
//
//    Output, double NODE_ATT[NODE_ATT_NUM*NODE_NUM], the nodal attributes.
//
//    Output, int NODE_MARKER[NODE_MARKER_NUM*NODE_NUM], the node markers.
//
{
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  ifstream input;
  int ival;
  int node;
  double value;

  node = -1;

  input.open ( node_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "NODE_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read, but ignore, dimension line.
//
    if ( node == -1 )
    {
      input >> i1 >> i2 >> i3 >> i4;
    }
    else
    {
      input >> ival;

      for ( i = 0; i < node_dim; i++ )
      {
        input >> value;
        node_coord[i+node*node_dim] = value;
      }
      for ( i = 0; i < node_att_num; i++ )
      {
        input >> value;
        node_att[i+node*node_att_num] = value;
      }
      for ( i = 0; i < node_marker_num; i++ )
      {
        input >> ival;
        node_marker[i+node*node_marker_num] = ival;
      }
    }

    node = node + 1;

    if ( node_num <= node )
    {
      break;
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void node_size_read ( string node_file, int *node_num, int *node_dim, 
  int *node_att_num, int *node_marker_num )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_SIZE_READ reads the header information from a node file.
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
//  Parameters:
//
//    Input, string NODE_FILE, the name of the node file to be read.
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *NODE_DIM, the spatial dimension.
//
//    Output, int *NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Output, int *NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
{
  ifstream input;

  input.open ( node_file.c_str ( ) );

  input >> *node_num
        >> *node_dim
        >> *node_att_num
        >> *node_marker_num;

  input.close ( );

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
