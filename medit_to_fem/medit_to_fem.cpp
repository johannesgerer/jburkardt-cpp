# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void i4mat_write ( string output_filename, int m, int n, int table[] );
void i4vec_zero ( int n, int a[] );
void mesh_data_read ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void mesh_size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void r8mat_write ( string output_filename, int m, int n, double table[] );
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

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MEDIT_TO_FEM.
//
//  Discussion:
//
//    MEDIT_TO_FEM converts mesh data from MEDIT to FEM format.
//
//  Usage:
//
//    medit_to_fem prefix
//
//    where 'prefix' is the common filename prefix:
//
//    * 'prefix'.mesh is the MEDIT mesh file.
//    * 'prefix'_nodes.txt will contain the node coordinates.
//    * 'prefix'_elements.txt will contain the element node connectivity.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 November 2014
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
  int element_num;
  int element_order;
  string fem_element_filename;
  string fem_node_filename;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int m;
  string medit_filename;
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
  cout << "MEDIT_TO_FEM\n";
  cout << "  C++ version:\n";
  cout << "  Read a mesh description created by the MEDIT program:\n";
  cout << "  * 'prefix'.mesh, the MEDIT mesh file.\n";
  cout << "  Write out two simple FEM files listing nodes and elements.\n";
  cout << "  * 'prefix'_nodes.txt, node coordinates.\n";
  cout << "  * 'prefix'_elements.txt, element connectivity.\n";
//
//  Get the filename prefix.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "  Please enter the filename prefix.\n";

    cin >> prefix;
  }
  else 
  {
    prefix = argv[1];
  }
//
//  Create the filenames.
//
  medit_filename = prefix + ".mesh";
  fem_node_filename = prefix + "_nodes.txt";
  fem_element_filename = prefix + "_elements.txt";
//
//  Read MEDIT sizes.
//
  mesh_size_read ( medit_filename, &dim, &vertices, &edges, &triangles, 
    &quadrilaterals, &tetrahedrons, &hexahedrons );
//
//  Report sizes.
//
  cout << "\n";
  cout << "  Number of dimensions = " << dim << "\n";
  cout << "  Number of vertices = " << vertices << "\n";
  cout << "  Number of edges = " << edges << "\n";
  cout << "  Number of triangles = " << triangles << "\n";
  cout << "  Number of quadrilaterals = " << quadrilaterals << "\n";
  cout << "  Number of tetrahedrons = " << tetrahedrons << "\n";
  cout << "  Number of hexahedrons = " << hexahedrons << "\n";
//
//  Allocate memory.
//
  edge_label = new int[edges];
  edge_vertex = new int[2 * edges];
  hexahedron_label = new int[hexahedrons];
  hexahedron_vertex = new int[8 * hexahedrons];
  quadrilateral_label = new int[quadrilaterals];
  quadrilateral_vertex = new int[4 * quadrilaterals];
  tetrahedron_label = new int[tetrahedrons];
  tetrahedron_vertex = new int[4 * tetrahedrons];
  triangle_label = new int[triangles];
  triangle_vertex = new int[3 * triangles];
  vertex_coordinate = new double[dim * vertices];
  vertex_label = new int[vertices];
//
//  Read MEDIT data.
//
  mesh_data_read ( medit_filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
//
//  Choose the FEM data.
//
//  We need to assume that there is only one element type.
//  If there are elements of multiple dimension, take the highest.
//
  m = dim;
  node_num = vertices;
  r8mat_write ( fem_node_filename, dim, vertices, vertex_coordinate );

  cout << "\n";
  cout << "  Created node coordinate file '" << fem_node_filename << "'\n";

  if ( 0 < hexahedrons && dim == 3 )
  {
    element_order = 8;
    element_num = hexahedrons;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      hexahedron_vertex );
  }
  else if ( 0 < tetrahedrons && dim == 3 )
  {
    element_order = 4;
    element_num = tetrahedrons;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      tetrahedron_vertex );
  }
  else if ( 0 < quadrilaterals && dim == 2 )
  {
    element_order = 4;
    element_num = quadrilaterals;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      quadrilateral_vertex );
  }
  else if ( 0 < triangles && dim == 2 )
  {
    element_order = 3;
    element_num = triangles;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      triangle_vertex );
  }
  else
  {
    cerr << "\n";
    cerr << "MEDIT_TO_FEM - Fatal error!\n";
    cerr << "  Unexpected values for spatial dimension\n";
    cerr << "  and number of nonzero objects.\n";
    exit ( 1 );
  }

  cout << "  Created element connectivity file '" << fem_element_filename << "'\n";
//
//  Free memory.
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
//
//  Terminate.
//
  cout << "\n";
  cout << "MEDIT_TO_FEM:\n";
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

void i4mat_write ( string output_filename, int m, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_WRITE writes an I4MAT file with no header.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "I4MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
//    Input, int TETRAHEDRONS, the number of tetrahedrons
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
//    Output, int *TETRAHEDRONS, the number of tetrahedrons
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

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
//    10 October 2014
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
    if ( s[n-1] != ' ' && s[n-1] != '\n' )
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

