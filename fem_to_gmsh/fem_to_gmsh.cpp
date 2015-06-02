# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
int file_column_count ( string filename );
int file_row_count ( string filename );
void gmsh_mesh1d_write ( string gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh2d_write ( string gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh3d_write ( string gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
int i4mat_max ( int m, int n, int a[] );
int i4mat_min ( int m, int n, int a[] );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM_TO_GMSH.
//
//  Discussion:
//
//    FEM_TO_GMSH converts a 1D, 2D or 3D mesh from FEM to GMSH format.
//
//  Usage:
//
//    fem_to_gmsh prefix
//
//    where 'prefix' is the common filename prefix:
//
//    * 'prefix'_nodes.txt contains the node coordinates,
//    * 'prefix'_elements.txt contains the element node connectivity.
//    * 'prefix'.msh will contain the Gmsh version of the data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int element;
  string element_filename;
  int *element_node;
  int element_num;
  int element_order;
  string gmsh_filename;
  int i;
  int m;
  int node;
  string node_filename;
  int node_num;
  double *node_x;
  string prefix;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "FEM_TO_GMSH\n";
  cout << "  C++ version:\n";
  cout << "  Convert a 1D, 2D or 3D mesh from FEM to GMSH format.\n";
  cout << "\n";
  cout << "  Read \"prefix\"_nodes.txt, node coordinates.\n";
  cout << "  Read \"prefix\"_elements.txt, element node connectivity.\n";
  cout << "\n";
  cout << "  Create \"prefix\".msh, a corresponding Gmsh mesh file.\n";
//
//  Get the filename prefix.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "FEM_TO_GMSH:\n";
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
  node_filename = prefix + "_nodes.txt";
  element_filename = prefix + "_elements.txt";
  gmsh_filename = prefix + ".msh";
//
//  Read the node data.
//
  r8mat_header_read ( node_filename, &m, &node_num );

  cout << "\n";
  cout << "  Read the header of \"" << node_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension = " << m << "\n";
  cout << "  Number of nodes  = " << node_num << "\n";

  node_x = r8mat_data_read ( node_filename, m, node_num );

  cout << "\n";
  cout << "  Read the data in \"" << node_filename << "\".\n";

  r8mat_transpose_print_some ( m, node_num, node_x, 1, 1, m, 5, 
    "  Portion of node coordinate data:" );
//
//  Read the element data.
//
  i4mat_header_read ( element_filename, &element_order, &element_num );

  if ( m == 1 )
  {
    if ( element_order == 2 )
    {
    }
    else
    {
    }
  }
  else if ( m == 2 )
  {
    if ( element_order == 3 )
    {
    }
    else if ( element_order == 6 )
    {
    }
    else
    {
      cerr << "\n";
      cerr << "FEM_TO_GMSH - Fatal error!\n";
      cerr << "  2D mesh data must use 3 or 6 nodes.\n";
      exit ( 1 );
    }
  }
  else if ( m == 3 )
  {
    if ( element_order == 4 )
    {
    }
    else if ( element_order == 10 )
    {
    }
    else if ( element_order == 20 )
    {
    }
    else
    {
      cerr << "\n";
      cerr << "FEM_TO_GMSH - Fatal error!\n";
      cerr << "  3D mesh data must use 4, 10, or 20 nodes.\n";
      exit ( 1 );
    }
  }
  cout << "\n";
  cout << "  Read the header of \"" << element_filename << "\".\n";
  cout << "\n";
  cout << "  Element order = " << element_order << "\n";
  cout << "  Number of elements  = " << element_num << "\n";

  element_node = i4mat_data_read ( element_filename, element_order,
    element_num );

  cout << "\n";
  cout << "  Read the data in \"" << element_filename << "\".\n";

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Initial portion of element data:" );
//
//  Write out the Gmsh version of the data.
//
  if ( m == 1 )
  {
    gmsh_mesh1d_write ( gmsh_filename, m, node_num, node_x, element_order,
      element_num, element_node );
  }
  else if ( m == 2 )
  {
    gmsh_mesh2d_write ( gmsh_filename, m, node_num, node_x, element_order,
      element_num, element_node );
  }
  else if ( m == 3 )
  {
    gmsh_mesh3d_write ( gmsh_filename, m, node_num, node_x, element_order,
      element_num, element_node );
  }

  cout << "\n";
  cout << "  Created the GMSH file \"" << gmsh_filename << "\".\n";
//
//  Free memory.
//
  delete [] element_node;
  delete [] node_x;
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM_TO_GMSH:\n";
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
//    13 June 2003
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
  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     

  return ( ch1 == ch2 );
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

int file_column_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file are presumed to consist of COLUMN_NUM words, 
//    separated by spaces.  There may also be some blank lines, and some 
//    comment lines, which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  char text[255];
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    column_num = -1;
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    exit ( 1 );
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    input.getline ( text, sizeof ( text ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) == 0 )
    {
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
    got_one = true;
    break;
  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( filename.c_str ( ) );

    for ( ; ; )
    {
      input.getline ( text, sizeof ( text ) );

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( text ) == 0 )
      {
        continue;
      }
      got_one = true;
      break;
    }
  }

  input.close ( );

  if ( !got_one )
  {
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  The file does not seem to contain any data.\n";
    exit ( 1 );
  }

  column_num = s_word_count ( text );

  return column_num;
}
//****************************************************************************80

int file_row_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  int record_num;
  int row_num;
  char text[255];

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the file: \"" << filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input.getline ( text, sizeof ( text ) );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( text[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( text ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;

  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

void gmsh_mesh1d_write ( string gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH1D_WRITE writes 1d mesh data as a Gmsh mesh file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Christophe Geuzaine, Jean-Francois Remacle,
//    Gmsh: a three-dimensional finite element mesh generator with
//    built-in pre- and post-processing facilities,
//    International Journal for Numerical Methods in Engineering,
//    Volume 79, Number 11, pages 1309-1331, 2009.
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the name of the Gmsh file.
//
//    Input, int M, the spatial dimension.
//
//    Input, inte NODE_NUM, the number of nodes.
//
//    Input, double NODE_X[M*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the nodes that make up each element.
//
{
  int element;
  int element_type;
  ofstream gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
//
//  Detect and correct 0-based node indexing.
//
  mesh_base_one ( node_num, element_order, element_num, element_node );
//
//  Open the file.
//
  gmsh.open ( gmsh_filename.c_str ( ) );
//
//  Write the data.
//
  gmsh << "$MeshFormat\n";
  gmsh << "2.2 0 8\n";
  gmsh << "$EndMeshFormat\n";

  gmsh << "$Nodes\n";
  gmsh << node_num << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    gmsh << node + 1 
         << "  " << node_x[0+node*m]
         << "  0.0  0.0\n";
  }
  gmsh << "$EndNodes\n";

  element_type = 1;

  tag_num = 2;
  tag1 = 0;
  gmsh << "$Elements\n";
  gmsh << element_num << "\n";
  for ( element = 0; element < element_num; element++ )
  {
    gmsh << element + 1
         << "  " << element_type
         << "  " << tag_num
         << "  " << tag1
         << "  " << element + 1;
    for ( i = 0; i < element_order; i++ )
    {
      gmsh << "  " << element_node[i+element*element_order];
    }
    gmsh << "\n";
  }
  gmsh << "$EndElements\n";

  gmsh.close ( );

  return;
}
//****************************************************************************80

void gmsh_mesh2d_write ( string gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH2D_WRITE writes 2d mesh data as a Gmsh mesh file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Christophe Geuzaine, Jean-Francois Remacle,
//    Gmsh: a three-dimensional finite element mesh generator with
//    built-in pre- and post-processing facilities,
//    International Journal for Numerical Methods in Engineering,
//    Volume 79, Number 11, pages 1309-1331, 2009.
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the name of the Gmsh file.
//
//    Input, int M, the spatial dimension.
//
//    Input, inte NODE_NUM, the number of nodes.
//
//    Input, double NODE_X[M*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the nodes that make up each element.
//
{
  int element;
  int element_type;
  ofstream gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
//
//  Detect and correct 0-based node indexing.
//
  mesh_base_one ( node_num, element_order, element_num, element_node );
//
//  Open the file.
//
  gmsh.open ( gmsh_filename.c_str ( ) );
//
//  Write the data.
//
  gmsh << "$MeshFormat\n";
  gmsh << "2.2 0 8\n";
  gmsh << "$EndMeshFormat\n";

  gmsh << "$Nodes\n";
  gmsh << node_num << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    gmsh << node + 1 
         << "  " << node_x[0+node*m]
         << "  " << node_x[1+node*m]
         << "  0.0\n";
  }
  gmsh << "$EndNodes\n";

  if ( element_order == 3 )
  {
    element_type = 2;
  }
  else if ( element_order == 6 )
  {
    element_type = 9;
  }

  tag_num = 2;
  tag1 = 0;
  gmsh << "$Elements\n";
  gmsh << element_num << "\n";
  for ( element = 0; element < element_num; element++ )
  {
    gmsh << element + 1
         << "  " << element_type
         << "  " << tag_num
         << "  " << tag1
         << "  " << element + 1;
    for ( i = 0; i < element_order; i++ )
    {
      gmsh << "  " << element_node[i+element*element_order];
    }
    gmsh << "\n";
  }
  gmsh << "$EndElements\n";

  gmsh.close ( );

  return;
}
//****************************************************************************80

void gmsh_mesh3d_write ( string gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH3D_WRITE writes 3D mesh data as a Gmsh mesh file.
//
//  Discussion:
//
//    The node ordering for the 20 node element is not standard.
//
//    Assuming the vertices are A, B, C and D, Gmsh uses the following ordering:
//
//    1:    a
//    2:        b
//    3:            c
//    4:                d
//    5: (2*a  +b        )/3
//    6: (  a+2*b        )/3
//    7: (    2*b+  c    )/3
//    8: (      b+2*c    )/3
//    9: (  a    +2*c    )/3
//   10: (2*a    +  c    )/3
//   11: (2*a        +  d)/3
//   12: (  a        +2*d)/3
//   13: (     b     +2*d)/3
//   14: (   2*b     +  d)/3
//   15: (       +  c+2*d)/3
//   16: (       +2*c+  d)/3
//   17: (  a+  b+  c    )/3
//   18: (  a+  b    +  d)/3
//   19: (      b+  c+  d)/3
//   20: (  a+      c+  d)/3
//
//    Leo Rebholz used the following ordering:
//
//    1:    a
//    2:        b
//    3:            c
//    4:                d
//    5: (2*a  +b        )/3
//    6: (2*a    +  c    )/3
//    7: (  a+2*b        )/3
//    8: (  a    +2*c    )/3
//    9: (  a+  b+  c    )/3
//   10: (    2*b+  c    )/3
//   11: (      b+2*c    )/3
//   12: (2*a        +  d)/3
//   13: (   2*b     +  d)/3
//   14: (       +2*c+  d)/3
//   15: (  a+  b    +  d)/3
//   16: (      b+  c+  d)/3
//   17: (  a+      c+  d)/3
//   18: (  a        +2*d)/3
//   19: (     b     +2*d)/3
//   20: (       +  c+2*d)/3
//
//    Since the only 20 node data we have is from Leo, we will assume that
//    all 20 node input data is in Leo's format, and needs to be converted
//    to the Gmsh convention.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Christophe Geuzaine, Jean-Francois Remacle,
//    Gmsh: a three-dimensional finite element mesh generator with
//    built-in pre- and post-processing facilities,
//    International Journal for Numerical Methods in Engineering,
//    Volume 79, Number 11, pages 1309-1331, 2009.
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the name of the Gmsh file.
//
//    Input, int M, the spatial dimension.
//
//    Input, inte NODE_NUM, the number of nodes.
//
//    Input, double NODE_X[M*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the nodes that make up each element.
//
{
  int element;
  int element_type;
  ofstream gmsh;
  int i;
  int i2;
  int leo_to_gmsh[20] = {
     0,  1,  2,  3,  4, 
     6,  9, 10,  7,  5, 
    11, 17, 18, 12, 19, 
    13,  8, 14, 15, 16 };
  int node;
  int tag_num;
  int tag1;
//
//  Detect and correct 0-based node indexing.
//
  mesh_base_one ( node_num, element_order, element_num, element_node );
//
//  Open the file.
//
  gmsh.open ( gmsh_filename.c_str ( ) );
//
//  Write the data.
//
  gmsh << "$MeshFormat\n";
  gmsh << "2.2 0 8\n";
  gmsh << "$EndMeshFormat\n";

  gmsh << "$Nodes\n";
  gmsh << node_num << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    gmsh << node + 1 
         << "  " << node_x[0+node*m]
         << "  " << node_x[1+node*m]
         << "  " << node_x[2+node*m] << "\n";
  }
  gmsh << "$EndNodes\n";

  if ( element_order == 4 )
  {
    element_type = 4;
  }
  else if ( element_order == 10 )
  {
    element_type = 11;
  }
  else if ( element_order == 20 )
  {
    element_type = 29;
  }

  tag_num = 2;
  tag1 = 0;
  gmsh << "$Elements\n";
  gmsh << element_num << "\n";
  for ( element = 0; element < element_num; element++ )
  {
    gmsh << element + 1
         << "  " << element_type
         << "  " << tag_num
         << "  " << tag1
         << "  " << element + 1;
    if ( element_order == 20 )
    {
      for ( i = 0; i < element_order; i++ )
      {
        i2 = leo_to_gmsh[i];
        gmsh << "  " << element_node[i2+element*element_order];
      }
      gmsh << "\n";
    }
    else
    {
      for ( i = 0; i < element_order; i++ )
      {
        gmsh << "  " << element_node[i+element*element_order];
      }
      gmsh << "\n";
    }
  }
  gmsh << "$EndElements\n";

  gmsh.close ( );

  return;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int *i4mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_DATA_READ reads data from an I4MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
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
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int I4MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  int *table;
  int *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "I4MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new int[m*n];

  x = new int[m];

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_i4vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80
 
void i4mat_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_HEADER_READ reads the header from an I4MAT file.
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
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points
//
{
  *m = file_column_count ( input_filename );
 
  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }
 
  *n = file_row_count ( input_filename );
 
  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }
 
  return;
}
//****************************************************************************80

int i4mat_max ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_MAX returns the maximum of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int I4MAT_MAX, the maximum entry of A.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int value;

  value = - i4_huge;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

int i4mat_min ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_MIN returns the minimum of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int I4MAT_MIN, the minimum entry of A.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int value;

  value = i4_huge;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < value )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title to be printed.
//
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title for the matrix.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j << ": ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ONE ensures that the element definition is one-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
{
  int element;
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  This will be converted to 1-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] + 1;
      }
    }
  }
  else if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ONE - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

double *r8mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DATA_READ reads the data from an R8MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
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
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double R8MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  double *table;
  double *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80
 
void r8mat_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HEADER_READ reads the header from an R8MAT file.
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
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points.
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << ": ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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
    if ( s[n-1] != ' ' && s[n-1] != '\n' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
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
//    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
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

int s_word_count ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
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
//    Input, string S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int char_count;
  int i;
  int word_count;

  word_count = 0;
  blank = true;

  char_count = s.length ( );

  for ( i = 0; i < char_count; i++ )
  {
    if ( isspace ( s[i] ) )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_count = word_count + 1;
      blank = false;
    }
  }

  return word_count;
}
//**********************************************************************

void timestamp ( )

//**********************************************************************
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
//    24 September 2003
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

