# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <string>

using namespace std;

# include "netcdf.hpp"

int main ( int argc, char **argv );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void i4vec_copy ( int n, int a1[], int a2[] );
void ice_write ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void mesh_data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void mesh_data_read ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void mesh_size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons );
void mesh_size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void r8vec_zero ( int n, double a[] );
bool s_begin ( string s1, string s2 );
bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
string s_newline_to_null ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char **argv )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_TO_ICE reads "ICE" data from a MESH file and writes it to a NETCDF file.
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
  cout << "MESH_TO_ICE:\n";
  cout << "  C++ version\n";
  cout << "  Read ICE data from a MESH file, write to a NETCDF file.\n";
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
  mesh_size_read ( filename_mesh, &dim, &vertices, &edges, &triangles,
    &quadrilaterals, &tetrahedrons, &hexahedrons );
//
//  Print sizes.
//
  mesh_size_print ( dim, vertices, edges, triangles, quadrilaterals,
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
  cout << "  Reading \"" << filename_mesh << "\".\n";

  mesh_data_read ( filename_mesh, dim, vertices, edges, triangles,
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
    mesh_data_print ( dim, vertices, edges, triangles, quadrilaterals,
      tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,
      edge_vertex, edge_label, triangle_vertex, triangle_label,
      quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex,
      tetrahedron_label, hexahedron_vertex, hexahedron_label );
  }
//
//  Write the data.
//
  cout << "\n";
  cout << "  Writing \"" << filename_nc << "\".\n";

  ice_write ( filename_nc, dim, vertices, edges, triangles,
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
  cout << "MESH_TO_ICE:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

char ch_cap ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  bool value;

  if ( 97 <= ch1 && ch1 <= 122 )
  {
    ch1 = ch1 - 32;
  }
  if ( 97 <= ch2 && ch2 <= 122 )
  {
    ch2 = ch2 - 32;
  }

  value = ( ch1 == ch2 );

  return value;
}
//****************************************************************************80

int ch_to_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     CH  DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the
//    character was 'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

void i4vec_zero ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
//****************************************************************************80

void ice_write ( string filename, int dim, int vertices, int edges,
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
//    ICE_WRITE writes 3D ICE sizes and data to a NETCDF file.
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
//    Input, string FILENAME, the name of the file to be created.
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
  NcDim *dim_dimension;
  NcDim *dim_edges;
  NcDim *dim_eight;
  NcDim *dim_four;
  NcDim *dim_hexahedrons;
  NcDim *dim_quadrilaterals;
  NcDim *dim_tetrahedrons;
  NcDim *dim_three;
  NcDim *dim_triangles;
  NcDim *dim_two;
  NcDim *dim_vertices;
  NcVar *var_edge_vertex;
  NcVar *var_edge_label;
  NcVar *var_hexahedron_vertex;
  NcVar *var_hexahedron_label;
  NcVar *var_quadrilateral_vertex;
  NcVar *var_quadrilateral_label;
  NcVar *var_tetrahedron_vertex;
  NcVar *var_tetrahedron_label;
  NcVar *var_triangle_vertex;
  NcVar *var_triangle_label;
  NcVar *var_vertex_coordinate;
  NcVar *var_vertex_label;
//
//  Create the file.
//
  NcFile dataFile ( filename.c_str ( ), NcFile::Replace );

  if ( !dataFile.is_valid ( ) )
  {
    cout << "\n";
    cout << "ICE_WRITE - Fatal error!\n";
    cout << "  Could not open the file.\n";
    exit ( 1 );
  }
//
//  Dimension information.
//
  dim_dimension = dataFile.add_dim ( "Dimension", dim );
  dim_vertices = dataFile.add_dim ( "Vertices", vertices );

  if ( 0 < edges )
  {
    dim_edges = dataFile.add_dim ( "Edges", edges );
  }

  if ( 0 < triangles )
  {
    dim_triangles = dataFile.add_dim ( "Triangles", triangles );
  }

  if ( 0 < quadrilaterals )
  {
    dim_quadrilaterals = dataFile.add_dim ( "Quadrilaterals", quadrilaterals );
  }

  if ( 0 < tetrahedrons )
  {
    dim_tetrahedrons = dataFile.add_dim ( "Tetrahedrons", tetrahedrons );
  }

  if ( 0 < hexahedrons )
  {
    dim_hexahedrons = dataFile.add_dim ( "Hexahedrons", hexahedrons );
  }

  dim_two = dataFile.add_dim ( "Two", 2 );
  dim_three = dataFile.add_dim ( "Three", 3 );
  dim_four = dataFile.add_dim ( "Four", 4 );
  dim_eight = dataFile.add_dim ( "Eight", 8 );
//
//  Define variables.
//
  if ( dim == 2 )
  {
    var_vertex_coordinate      = dataFile.add_var ( "Vertex_Coordinate",    ncDouble, dim_two, dim_vertices );
  }
  else if ( dim == 3 )
  {
    var_vertex_coordinate      = dataFile.add_var ( "Vertex_Coordinate",    ncDouble, dim_three, dim_vertices );
  }
  var_vertex_label           = dataFile.add_var ( "Vertex_Label",         ncInt,               dim_vertices );

  if ( 0 < edges )
  {
    var_edge_vertex          = dataFile.add_var ( "Edge_Vertex",          ncInt,    dim_two,   dim_edges );
    var_edge_label           = dataFile.add_var ( "Edge_Label",           ncInt,               dim_edges );
  }

  if ( 0 < triangles )
  {
    var_triangle_vertex      = dataFile.add_var ( "Triangle_Vertex",      ncInt, dim_three, dim_triangles );
    var_triangle_label       = dataFile.add_var ( "Triangle_Label",       ncInt,            dim_triangles );
  }

  if ( 0 < quadrilaterals )
  {
    var_quadrilateral_vertex = dataFile.add_var ( "Quadrilateral_Vertex", ncInt, dim_four, dim_quadrilaterals );
    var_quadrilateral_label  = dataFile.add_var ( "Quadrilateral_Label",  ncInt,           dim_quadrilaterals );
  }

  if ( 0 < tetrahedrons )
  {
    var_tetrahedron_vertex   = dataFile.add_var ( "Tetrahedron_Vertex",   ncInt, dim_four,  dim_tetrahedrons );
    var_tetrahedron_label    = dataFile.add_var ( "Tetrahedron_Label",    ncInt,            dim_tetrahedrons );
  }

  if ( 0 < hexahedrons )
  {
    var_hexahedron_vertex    = dataFile.add_var ( "Hexahedron_Vertex",    ncInt, dim_eight, dim_hexahedrons );
    var_hexahedron_label     = dataFile.add_var ( "Hexahedron_Label",     ncInt,            dim_hexahedrons );
  }
//
//  Write the data.
//
  var_vertex_coordinate->put      ( &vertex_coordinate[0], dim, vertices );
  var_vertex_label->put           ( &vertex_label[0],            vertices );
  if ( 0 < edges )
  {
    var_edge_vertex->put          ( &edge_vertex[0],          2, edges );
    var_edge_label->put           ( &edge_label[0],              edges );
  }
  if ( 0 < triangles )
  {
    var_triangle_vertex->put      ( &triangle_vertex[0],      3, triangles );
    var_triangle_label->put       ( &triangle_label[0],          triangles );
  }
  if ( 0 < quadrilaterals )
  {
    var_quadrilateral_vertex->put ( &quadrilateral_vertex[0], 4, quadrilaterals );
    var_quadrilateral_label->put  ( &quadrilateral_label[0],     quadrilaterals );
  }
  if ( 0 < tetrahedrons )
  {
    var_tetrahedron_vertex->put   ( &tetrahedron_vertex[0],   4, tetrahedrons );
    var_tetrahedron_label->put    ( &tetrahedron_label[0],       tetrahedrons );
  }
  if ( 0 < hexahedrons )
  {
    var_hexahedron_vertex->put    ( &hexahedron_vertex[0],      8, hexahedrons );
    var_hexahedron_label->put     ( &hexahedron_label[0],          hexahedrons );
  }
//
//  Close the file.
//
  dataFile.close ( );

  return;
}
//****************************************************************************80

void mesh_data_print ( int dim, int vertices, int edges, int triangles,
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
//    MESH_DATA_PRINT prints the data of a MESH dataset.
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

void mesh_data_read ( string filename, int dim, int vertices, int edges,
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
//    MESH_DATA_READ reads data from a MESH file.
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
//    Input, string FILENAME, the name of the MESH file.
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals
//    (may be 0).
//
//    Input, int TETRAHEDRAONS, the number of tetrahedrons
//    (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons
//    (may be 0).
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
  int dim2;
  int edge;
  int edges2;
  int hexahedron;
  int hexahedrons2;
  int i;
  int i4vec[9];
  int ierror;
  ifstream input;
  string keyword;
  int length;
  int line_num;
  int quadrilateral;
  int quadrilaterals2;
  double r8vec[9];
  int tetrahedron;
  int tetrahedrons2;
  string text;
  int triangle;
  int triangles2;
  int vertex;
  int vertices2;
//
//  Initialize everything to nothing.
//
  i4vec_zero ( edges, edge_label );
  i4vec_zero ( 2 * edges, edge_vertex );
  i4vec_zero ( hexahedrons, hexahedron_label );
  i4vec_zero ( 8 * hexahedrons, hexahedron_vertex );
  i4vec_zero ( quadrilaterals, quadrilateral_label );;
  i4vec_zero ( 4 * quadrilaterals, quadrilateral_vertex );
  i4vec_zero ( tetrahedrons, tetrahedron_label );
  i4vec_zero ( 4 * tetrahedrons, tetrahedron_vertex );
  i4vec_zero ( triangles, triangle_label );
  i4vec_zero ( 3 * triangles, triangle_vertex );
  r8vec_zero ( dim * vertices, vertex_coordinate );
  i4vec_zero ( vertices, vertex_label );
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MESH_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Read lines til you get alphanumerics and determine a "mode"
//
  line_num = 0;
  keyword = "NONE";

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    line_num = line_num + 1;

    if ( s_len_trim ( text ) == 0 )
    {
      keyword = "NONE";
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
//
//  Remove initial blanks.
//

//
//  Expecting a keyword.
//
        if ( s_eqi ( text, "CORNERS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "DIMENSION" ) )
    {
      keyword = "DIMENSION";
    }
    else if ( s_eqi ( text, "EDGES" ) )
    {
      keyword = "EDGES";
    }
    else if ( s_eqi ( text, "END" ) )
    {
      cout << "\n";
      cout << "  END statement encountered.\n";
      break;
    }
    else if ( s_eqi ( text, "HEXAHEDRA" ) ||
              s_eqi ( text, "HEXAHEDRONS" ) )
    {
      keyword = "HEXAHEDRONS";
    }
    else if ( s_begin ( text, "MESHVERSIONFORMATTED" ) )
    {
    }
    else if ( s_eqi ( text, "NORMALATQUADRILATERALVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATTRIANGLEVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "QUADRILATERALS" ) )
    {
      keyword = "QUADRILATERALS";
    }
    else if ( s_eqi ( text, "REQUIREDEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "REQUIREDVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "RIDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTATEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TETRAHEDRA" ) ||
              s_eqi ( text, "TETRAHEDRONS" ) )
    {
      keyword = "TETRAHEDRONS";
    }
    else if ( s_eqi ( text, "TRIANGLES" ) )
    {
      keyword = "TRIANGLES";
    }
    else if ( s_eqi ( text, "VERTICES" ) )
    {
      keyword = "VERTICES";
    }
//
//  Presumably, numeric data to be processed by keyword.
//
    else if ( s_eqi ( keyword, "DIMENSION" ) )
    {
      dim2 = atoi ( text.c_str ( ) );
      keyword = "NONE";
    }
    else if ( s_eqi ( keyword, "EDGES" ) )
    {
      edges2 = atoi ( text.c_str ( ) );
      keyword = "EDGE_VERTEX";
      edge = 0;
    }
    else if ( s_eqi ( keyword, "EDGE_VERTEX" ) )
    {
      s_to_i4vec ( text, 3, i4vec );
      for ( i = 0; i < 2; i++ )
      {
        edge_vertex[i+edge*2] = i4vec[i];
      }
      edge_label[edge] = i4vec[2];
      edge = edge + 1;
    }
    else if ( s_eqi ( keyword, "HEXAHEDRONS" ) )
    {
      hexahedrons2 = atoi ( text.c_str ( ) );
      keyword = "HEXAHEDRON_VERTEX";
      hexahedron = 0;
    }
    else if ( s_eqi ( keyword, "HEXAHEDRON_VERTEX" ) )
    {
      s_to_i4vec ( text, 9, i4vec );
      for ( i = 0; i < 8; i++ )
      {
        hexahedron_vertex[i+hexahedron*8] = i4vec[i];
      }
      hexahedron_label[hexahedron] = i4vec[8];
      hexahedron = hexahedron + 1;
    }
    else if ( s_eqi ( keyword, "QUADRILATERALS" ) )
    {
      quadrilaterals2 = atoi ( text.c_str ( ) );
      keyword = "QUADRILATERAL_VERTEX";
      quadrilateral = 0;
    }
    else if ( s_eqi ( keyword, "QUADRILATERAL_VERTEX" ) )
    {
      s_to_i4vec ( text, 5, i4vec );
      for ( i = 0; i < 4; i++ )
      {
        quadrilateral_vertex[i+quadrilateral*4] = i4vec[i];
      }
      quadrilateral_label[quadrilateral] = i4vec[4];
      quadrilateral = quadrilateral + 1;
    }
    else if ( s_eqi ( keyword, "TETRAHEDRONS" ) )
    {
      tetrahedrons2 = atoi ( text.c_str ( ) );
      keyword = "TETRAHEDRON_VERTEX";
      tetrahedron = 0;
    }
    else if ( s_eqi ( keyword, "TETRAHEDRON_VERTEX" ) )
    {
      s_to_i4vec ( text, 5, i4vec );
      for ( i = 0; i < 4; i++ )
      {
        tetrahedron_vertex[i+tetrahedron*4] = i4vec[i];
      }
      tetrahedron_label[tetrahedron] = i4vec[4];
      tetrahedron = tetrahedron + 1;
    }
    else if ( s_eqi ( keyword, "TRIANGLES" ) )
    {
      triangles2 = atoi ( text.c_str ( ) );
      keyword = "TRIANGLE_VERTEX";
      triangle = 0;
    }
    else if ( s_eqi ( keyword, "TRIANGLE_VERTEX" ) )
    {
      s_to_i4vec ( text, 4, i4vec );
      for ( i = 0; i < 3; i++ )
      {
        triangle_vertex[i+triangle*3] = i4vec[i];
      }
      triangle_label[triangle] = i4vec[3];
      triangle = triangle + 1;
    }
    else if ( s_eqi ( keyword, "VERTICES" ) )
    {
      vertices2 = atoi ( text.c_str ( ) );
      keyword = "VERTEX_COORDINATE";
      vertex = 0;
    }
    else if ( s_eqi ( keyword, "VERTEX_COORDINATE" ) )
    {
      s_to_r8vec ( text, dim + 1, r8vec );
      for ( i = 0; i < dim; i++ )
      {
        vertex_coordinate[i+vertex*dim] = r8vec[i];
      }
      vertex_label[vertex] = ( int ) r8vec[dim];
      vertex = vertex + 1;
    }
    else if ( s_eqi ( keyword, "SKIP" ) )
    {
    }
    else
    {
      cerr << "\n";
      cerr << "MESH_DATA_READ - Fatal error!\n";
      cerr << "  Could not find keyword while reading line "
           << line_num << "\n";
      cerr << "\"" << text << "\".\n";
      exit ( 1 );
    }
  }
//
//  Close the file.
//
  input.close ( );

  cout << "\n";
  cout << "  Read " << line_num << " lines from \"" << filename << "\".\n";

  return;
}
//****************************************************************************80

void mesh_size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_SIZE_PRINT prints the sizes of an ICE dataset.
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
//****************************************************************************80

void mesh_size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_SIZE_READ reads sizes from a MESH file.
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
//    Input, string FILENAME, the name of the MESH file.
//
//    Output, int *DIM, the spatial dimension, which should be 2 or 3.
//
//    Output, int *VERTICES, the number of vertices.
//
//    Output, int *EDGES, the number of edges (may be 0).
//
//    Output, int *TRIANGLES, the number of triangles (may be 0).
//
//    Output, int *QUADRILATERALS, the number of quadrilaterals
//    (may be 0).
//
//    Output, int *TETRAHEDRAONS, the number of tetrahedrons
//    (may be 0).
//
//    Output, int *HEXAHEDRONS, the number of hexahedrons
//    (may be 0).
//
{
  int ierror;
  ifstream input;
  string keyword;
  int length;
  int line_num;
  string text;
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
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MESH_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Read lines til you get alphanumerics and determine a "mode"
//
  line_num = 0;
  keyword = "NONE";

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    line_num = line_num + 1;

    if ( s_len_trim ( text ) == 0 )
    {
      keyword = "NONE";
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
//
//  Remove initial blanks.
//

//
//  Expecting a keyword.
//
        if ( s_eqi ( text, "CORNERS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "DIMENSION" ) )
    {
      keyword = "DIMENSION";
    }
    else if ( s_eqi ( text, "EDGES" ) )
    {
      keyword = "EDGES";
    }
    else if ( s_eqi ( text, "END" ) )
    {
      cout << "\n";
      cout << "  END statement encountered.\n";
      break;
    }
    else if ( s_eqi ( text, "HEXAHEDRA" ) ||
              s_eqi ( text, "HEXAHEDRONS" ) )
    {
      keyword = "HEXAHEDRONS";
    }
    else if ( s_begin ( text, "MESHVERSIONFORMATTED" ) )
    {
    }
    else if ( s_eqi ( text, "NORMALATQUADRILATERALVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATTRIANGLEVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "QUADRILATERALS" ) )
    {
      keyword = "QUADRILATERALS";
    }
    else if ( s_eqi ( text, "REQUIREDEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "REQUIREDVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "RIDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTATEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TETRAHEDRA" ) ||
              s_eqi ( text, "TETRAHEDRONS" ) )
    {
      keyword = "TETRAHEDRONS";
    }
    else if ( s_eqi ( text, "TRIANGLES" ) )
    {
      keyword = "TRIANGLES";
    }
    else if ( s_eqi ( text, "VERTICES" ) )
    {
      keyword = "VERTICES";
    }
//
//  Presumably, numeric data to be processed by keyword.
//
    else if ( s_eqi ( keyword, "DIMENSION" ) )
    {
      *dim = atoi ( text.c_str ( ) );
      keyword = "NONE";
    }
    else if ( s_eqi ( keyword, "EDGES" ) )
    {
      *edges = atoi ( text.c_str ( ) );
      keyword = "EDGE_VERTEX";
    }
    else if ( s_eqi ( keyword, "EDGE_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "HEXAHEDRONS" ) )
    {
      *hexahedrons = atoi ( text.c_str ( ) );
      keyword = "HEXAHEDRON_VERTEX";
    }
    else if ( s_eqi ( keyword, "HEXAHEDRON_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "QUADRILATERALS" ) )
    {
      *quadrilaterals = atoi ( text.c_str ( ) );
      keyword = "QUADRILATERAL_VERTEX";
    }
    else if ( s_eqi ( keyword, "QUADRILATERAL_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "TETRAHEDRONS" ) )
    {
      *tetrahedrons = atoi ( text.c_str ( ) );
      keyword = "TETRAHEDRON_VERTEX";
    }
    else if ( s_eqi ( keyword, "TETRAHEDRON_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "TRIANGLES" ) )
    {
      *triangles = atoi ( text.c_str ( ) );
      keyword = "TRIANGLE_VERTEX";
    }
    else if ( s_eqi ( keyword, "TRIANGLE_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "VERTICES" ) )
    {
      *vertices = atoi ( text.c_str ( ) );
      keyword = "VERTEX_COORDINATE";
    }
    else if ( s_eqi ( keyword, "VERTEX_COORDINATE" ) )
    {
    }
    else if ( s_eqi ( keyword, "SKIP" ) )
    {
    }
    else
    {
      cerr << "\n";
      cerr << "MESH_SIZE_READ - Fatal error!\n";
      cerr << "  Could not find keyword while reading line "
           << line_num << "\n";
      cerr << "\"" << text << "\".\n";
      exit ( 1 );
    }
  }
//
//  Close the file.
//
  input.close ( );

  cout << "\n";
  cout << "  Read " << line_num << " lines from \"" << filename << "\".\n";

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

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO zeroes an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
//****************************************************************************80

bool s_begin ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_BEGIN reports whether string 1 begins with string 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, string S2, two strings.
//
//    Output, bool S_BEGIN, is true if S1 is the same as S2 up to
//    the end of S2, and false otherwise.
//
{
  int i;
  int n;
  int n1;
  int n2;

  n1 = s1.length ( );
  n2 = s2.length ( );

  if ( n1 < n2 )
  {
    return false;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal.
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ )
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length )
  {
    for ( i = nchar; i < s1_length; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return false;
      }
    }
  }
  else if ( nchar < s2_length )
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' )
      {
        return false;
      }
    }
  }

  return true;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

string s_newline_to_null ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_NEWLINE_TO_NULL replaces carriage returns or newlines by nulls.
//
//  Discussion:
//
//    The function FGETS will read a string containing a line of text read from
//    input.  However, the string will include the linefeed character '/n', or,
//    for a PC-formatted file, the carriage return and linefeed pair '/r' + '/n'.
//
//    It may be desirable that the string not contain these characters.  The
//    easiest way to deal with this is simply to replace the first instance of
//    '/r' or '/n' by a null character, which terminates the string.
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
//    Input, string S, the string to be modified.
//
//    Output, string S_NEWLINE_TO_NULL, the modified string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
//
//  Handle carriage return;
//
    if ( s[i] == '\r' )
    {
      s2[i] = '\0';
      break;
    }
//
//  Handle linefeed.
//
    if ( s[i] == '\n' )
    {
      s2[i] = '\0';
      break;
    }
  }
  return s2;
}
//****************************************************************************80

int s_to_i4 ( string s, int *last, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an I4 from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be examined.
//
//    Output, int *LAST, the last character of S used to make IVAL.
//
//    Output, bool *ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  for ( ; ; )
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = true;
    *last = 0;
  }

  return ival;
}
//****************************************************************************80

bool s_to_i4vec ( string s, int n, int ivec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4VEC reads an I4VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int IVEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

double s_to_r8 ( string s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( string s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
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
