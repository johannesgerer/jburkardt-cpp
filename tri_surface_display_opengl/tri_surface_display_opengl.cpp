# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstring>

//
//  This is the include statement I need on for Mac OS X.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

int main ( int argc, char *argv[] );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void display ( );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void mesh_base_zero ( int node_num, int element_order,
  int element_num, int element_node[] );
void mouse ( int btn, int state, int x, int y );
void myinit ( );
void myReshape ( int w, int h );
double *node_normal_set ( );;
double r8_max ( double x, double y );
double *r83vec_max ( int n, double a[] );
double *r83vec_min ( int n, double a[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void spinSurface ( );
void timestamp ( );
//
//  Global data.
//
  static GLint axis = 2;
  static GLfloat theta[3] = { 0.0, 0.0, 0.0 };

  int dim_num = 0;
  int *element_node = NULL;
  int element_num = 0;
  int element_order = 0;
  int node_num = 0;
  double *node_normal = NULL;
  double *node_xyz = NULL;
  double *node_xyz_max = NULL;
  double *node_xyz_min = NULL;
  double node_xyz_range[3];
  int pixel_height = 0;
  int pixel_width = 0;
  bool spinning = true;
  double theta_speed = 0.020;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRI_SURFACE_DISPLAY_OPENGL.
//
//  Discussion:
//
//    This program reads two files defining the nodes and the faces
//    of a triangulated surface in 3D.
//
//    It displays the surface using OpenGL.
//
//    Since rotation is about the point (0,0,0), I have to compute
//    the centroid and subtract it from the node coordinates.
//
//  Usage:
//
//    tri_surface_display_opengl prefix
//
//    where 'prefix' is the common prefix for the files:
//
//    * prefix_nodes.txt,     the node coordinates;
//    * prefix_elements.txt,  the nodes that make up each element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  double *centroid;
  int element;
  string element_filename;
  int i;
  int j;
  int node;
  string node_filename;
  int node_max;
  int node_min;
  int order;
  string prefix;

  cout << "\n";
  cout << "TRI_SURFACE_DISPLAY_OPENGL:\n";
  cout << "  C++ version:\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  This program reads files defining the nodes and elements\n";
  cout << "  of a triangular mesh that forms a surface in 3D. \n";
  cout << "\n";
  cout << "  An image of the surface is displayed using OpenGL.\n";
  cout << "\n";
  cout << "  The image rotates slowly around the X, Y or Z axis.\n";
  cout << "  Click the mouse to change the axis.\n";
//
//  Get the common filename prefix.
//
  if ( argc <= 1 )
  {
    cout << "\n";
    cout << "TRI_SURFACE_DISPLAY_OPENGL:\n";
    cout << "  Please enter the common filename prefix.\n";

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

//
//  Read the node data.
//
  r8mat_header_read ( node_filename, &dim_num, &node_num );

  cout << "\n";
  cout << "  The node file has been examined.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM =  " << dim_num << "\n";
  cout << "  The number of nodes NODE_NUM =   " << node_num << "\n";

  if ( dim_num != 3 )
  {
    cout << "\n";
    cout << "TRI_SURFACE_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The spatial dimension of the data should be 3.\n";
    exit ( 1 );
  }

  node_xyz = r8mat_data_read ( node_filename, dim_num, node_num );

  node_xyz_min = r83vec_min ( node_num, node_xyz );
  node_xyz_max = r83vec_max ( node_num, node_xyz );

  node_xyz_range[0] = node_xyz_max[0] - node_xyz_min[0];
  node_xyz_range[1] = node_xyz_max[1] - node_xyz_min[1];
  node_xyz_range[2] = node_xyz_max[2] - node_xyz_min[2];

  cout << "\n";
  cout << "  Minimum: " << node_xyz_min[0]
       << "  "          << node_xyz_min[1]
       << "  "          << node_xyz_min[2] << "\n";
  cout << "  Maximum: " << node_xyz_max[0]
       << "  "          << node_xyz_max[1]
       << "  "          << node_xyz_max[2] << "\n";
  cout << "  Range:   " << node_xyz_range[0]
       << "  "          << node_xyz_range[1]
       << "  "          << node_xyz_range[2] << "\n";

  if ( node_xyz_range[0] == 0.0 )
  {
    cout << "\n";
    cout << "TRI_SURFACE_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The X data range is 0.\n";
    exit ( 1 );
  }

  if ( node_xyz_range[1] == 0.0 )
  {
    cout << "\n";
    cout << "TRI_SURFACE_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The Y data range is 0.\n";
    exit ( 1 );
  }

  if ( node_xyz_range[2] == 0.0 )
  {
    cout << "\n";
    cout << "TRI_SURFACE_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The Z data range is 0.\n";
    exit ( 1 );
  }

  r8mat_transpose_print_some ( dim_num, node_num, node_xyz, 1, 1, 3, 5,
    "  First five nodes:" );
//
//  Read the element data.
//
  i4mat_header_read ( element_filename, &element_order, &element_num );

  cout << "\n";
  cout << "  The element file has been examined.\n";
  cout << "\n";
  cout << "  The element order =  " << element_order << "\n";
  cout << "  The number of elements =  " << element_num << "\n";

  if ( element_order < 3 )
  {
    cout << "\n";
    cout << "TRI_SURFACE_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The element order must be at least 3!\n";
    exit ( 1 );
  }

  element_node = i4mat_data_read ( element_filename, element_order, element_num );

  i4mat_transpose_print_some ( element_order, element_num, element_node,
    1, 1, element_order, 5, "  First five elements:" );
//
//  Detect and correct 1-based node indexing.
//
  mesh_base_zero ( node_num, element_order, element_num, element_node );
//
//  Since the spin function works around (0,0,0), we have to compute the
//  centroid and subtract it.
//
  centroid = new double[3];
  for ( i = 0; i < 3; i++ )
  {
    centroid[i] = 0.0;
  }

  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      centroid[i] = centroid[i] + node_xyz[i+j*3];
    }
  }
  for ( i = 0; i < 3; i++ )
  {
    centroid[i] = centroid[i] / ( double ) node_num;
  }

  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_xyz[i+j*3] = node_xyz[i+j*3] - centroid[i];
    }
  }
//
//  Compute the node normals.
//
  node_normal = node_normal_set ( );

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
//
//  What is appropriate here?  What are the projection axes
//  for a 3D plot?
//
  if ( node_xyz_range[1] < node_xyz_range[0] )
  {
    pixel_width = 500;
    pixel_height = ( int )
      ( ( double ) ( 500 ) * node_xyz_range[1] / node_xyz_range[0] );
  }
  else
  {
    pixel_width = ( int )
      ( ( double ) ( 500 ) * node_xyz_range[0] / node_xyz_range[1] );
    pixel_height = 500;
  }
  cout << "  Pixels:  " << pixel_width  << "  "
                        << pixel_height << "\n";

  glutInitWindowSize ( pixel_width, pixel_height );

  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Rotating Triangulated Surface" );
  glutReshapeFunc ( myReshape );
  glutDisplayFunc ( display );
  glutIdleFunc ( spinSurface );
  glutMouseFunc ( mouse );
//
//  Enable hidden surface removal.
//
  glEnable ( GL_DEPTH_TEST );
//
//  Do "my" initializations.
//
  myinit ( );

  glutMainLoop ( );
//
//  Free memory.
//
  delete [] centroid;
  delete [] node_xyz;
  delete [] element_node;
//
//  Terminate.
//
  cout << "\n";
  cout << "TRI_SURFACE_DISPLAY_OPENGL:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );
  return 0;
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
//    Output, int CH_TO_DIGIT, the corresponding value.  If the character was
//    'illegal', then DIGIT is -1.
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

void display ( )

//****************************************************************************80
//
//  Purpose:
//
//    DISPLAY generates the graphics output.
//
//  Discussion;
//
//    Display a tetrahedral mesh as a collection of lines and nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  double color[3] = { 0.0, 1.0, 0.0 };
  int element;
  int node;
  double norm[3];
  int order;
  double p[3];
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
//
//  ?
//
  glLoadIdentity ( );

  glRotatef ( theta[0], 1.0, 0.0, 0.0 );
  glRotatef ( theta[1], 0.0, 1.0, 0.0 );
  glRotatef ( theta[2], 0.0, 0.0, 1.0 );
//
//  Specify that elements facing the front are to be filled.
//
  glPolygonMode ( GL_FRONT, GL_FILL );
//
//  Draw the elements, in GREEN.
//
  glColor3f ( 0.0, 1.0, 0.0 );

  for ( element = 0; element < element_num; element++ )
  {
    glBegin ( GL_POLYGON );
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];

      glColor3dv ( color );

      norm[0] = node_normal[0+node*3];
      norm[1] = node_normal[1+node*3];
      norm[2] = node_normal[2+node*3];
      glNormal3dv ( norm );

      p[0] = node_xyz[0+node*3];
      p[1] = node_xyz[1+node*3];
      p[2] = node_xyz[2+node*3];
      glVertex3dv ( p );
    }
    glEnd ( );
  }
//
//  Draw lines delimiting the elements, in RED.
//
  glColor3f ( 1.0, 0.0, 0.0 );

  for ( element = 0; element < element_num; element++ )
  {
    glBegin ( GL_LINE_LOOP );
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      p[0] = node_xyz[0+node*3];
      p[1] = node_xyz[1+node*3];
      p[2] = node_xyz[2+node*3];
      glVertex3dv ( p );
    }
    glEnd ( );
  }
//
//  Draw the nodes in BLUE.
//  Draw them AFTER the lines, so they are more likely to show up!
//
  glColor3f ( 0.0, 0.0, 1.0 );

  for ( node = 0; node < node_num; node++ )
  {
    glBegin ( GL_POINTS );
      p[0] = node_xyz[0+node*3];
      p[1] = node_xyz[1+node*3];
      p[2] = node_xyz[2+node*3];
      glVertex3dv ( p );
    glEnd ( );
  }
//
//  Clear all the buffers.
//
  glFlush ( );
//
//  Switch between the two buffers for fast animation.
//
  glutSwapBuffers ( );

  return;
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
//    05 July 2009
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
  string text;
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
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) <= 0 )
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
      input >> text;

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
    cerr << "FILE_COLUMN_COUNT - Warning!\n";
    cerr << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( text );

  return column_num;
}
//****************************************************************************80

int file_row_count ( string input_filename )

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
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  string line;
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return (-1);
  }

  for ( ; ; )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
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
//    I4_MIN returns the smaller of two I4's.
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
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
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
      cout << setw(5) << j << "  ";
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

void mesh_base_zero ( int node_num, int element_order,
  int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ZERO ensures that the element definition is zero-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//    Since a C++ program will naturally assume a 0-based indexing, it is
//    necessary to check a given element definition and, if it is actually
//    1-based, to convert it.
//
//    This function attempts to detect 1-based node indexing and correct it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2010
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
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = node_num + 1;
  node_max = -1;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      node_min = i4_min ( node_min, node );
      node_max = i4_max ( node_max, node );
    }
  }

  if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  This will be converted to 0-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] - 1;
      }
    }
  }
  else if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else if ( ( 1 == node_min && node_max <  node_num ) ||
            ( 1 <  node_min && node_max == node_num ) )
  {
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] - 1;
      }
    }
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  Node indices will be decremented.\n";
  }
  else if ( ( 0 == node_min && node_max <  node_num - 1 ) ||
            ( 0 <  node_min && node_max == node_num - 1 ) )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else if ( node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  Node indices will be used without change.\n";
  }
  return;
}
//****************************************************************************80

void mouse ( int btn, int state, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    MOUSE determines the response to mouse input.
//
//  Discussion:
//
//    The original routine assumed the user had a three button mouse, and
//    dedicated one axis to each.  Since Apple prefers the esthetics of a
//    one button mouse, this routine simply increments the axis by 1,
//    no matter which button is pushed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2008
//
//  Author:
//
//    Edward Angel
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  if ( btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
  {
    if ( spinning )
    {
      spinning = false;
      theta_speed = 0.0;
    }
    else
    {
      spinning = true;
      axis = axis + 1;
      theta_speed = 0.020;
    }
  }
  if ( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN )
  {
    if ( spinning )
    {
      spinning = false;
      theta_speed = 0.0;
    }
    else
    {
      spinning = true;
      axis = axis + 1;
      theta_speed = 0.020;
    }
  }
  if ( btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN )
  {
    if ( spinning )
    {
      spinning = false;
      theta_speed = 0.0;
    }
    else
    {
      spinning = true;
      axis = axis + 1;
      theta_speed = 0.020;
    }
  }
  axis = axis % 3;

  return;
}
//****************************************************************************80

void myinit ( )

//****************************************************************************80
//
//  Purpose:
//
//    MYINIT initializes OpenGL state variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  double margin;
  double x_max;
  double x_min;
  double y_max;
  double y_min;
  double z_max;
  double z_min;
//
//  Set the background to WHITE.
//
  glClearColor ( 1.0, 1.0, 1.0, 1.0 );
//
//  Make vertices bigger than the default size of 1.0.
//
  glPointSize ( 5.0 );
//
//  Set up the viewing window with origin at the lower left.
//
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
//
//  Determine an amount MARGIN by which it would be appropriate to spread the
//  data range, so that all the data is comfortably inside the picture.
//
  margin =                  node_xyz_max[0] - node_xyz_min[0];
  margin = r8_max ( margin, node_xyz_max[1] - node_xyz_min[1] );
  margin = r8_max ( margin, node_xyz_max[2] - node_xyz_min[2] );

  margin = 0.025 * margin;

  x_min = ( double ) ( node_xyz_min[0] - margin );
  x_max = ( double ) ( node_xyz_max[0] + margin );
  y_min = ( double ) ( node_xyz_min[1] - margin );
  y_max = ( double ) ( node_xyz_max[1] + margin );
  z_min = ( double ) ( node_xyz_min[2] - margin );
  z_max = ( double ) ( node_xyz_max[2] + margin );
//
//  Specify the clipping volume.
//
  glOrtho ( x_min, x_max, y_min, y_max, z_min, z_max );

  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

void myReshape ( int w, int h )

//****************************************************************************80
//
//  Purpose:
//
//    MYRESHAPE determines the window mapping.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2008
//
//  Author:
//
//    Edward Angel
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  glViewport ( 0, 0, w, h );
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );

  if ( w <= h )
  {
    glOrtho (
      -2.0, 2.0,
      -2.0 * ( GLfloat ) h / ( GLfloat ) w, 2.0 * ( GLfloat ) h / ( GLfloat ) w,
      -10.0, 10.0 );
  }
  else
  {
    glOrtho (
      -2.0 * ( GLfloat ) h / ( GLfloat ) w, 2.0 * ( GLfloat ) h / ( GLfloat ) w,
      -2.0, 2.0,
      -10.0, 10.0 );
  }

  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

double *node_normal_set ( )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_NORMAL_SET computes node normal vectors.
//
//  Discussion:
//
//    We assume we are given no normal vector information to start with.
//    We then consider each face of the surface.  We take every
//    set of three consecutive vertices, and construct the normal vector
//    to the triangle that they define.  We add this normal vector to the
//    normal vector information for the node associated with the middle of
//    the three vertices.
//
//    Once we have processed all the vertices on a face, and all the faces,
//    we average the node normal information.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  int element;
  int i;
  int n0;
  int n1;
  int n2;
  int node;
  double *node_normal;
  double triangle_normal[3];
  double norm;
  double v1[3];
  double v2[3];
  int vert0;
  int vert1;
  int vert2;
//
//  Make space.
//
  node_normal = new double[3*node_num];
//
//  Zero out the space.
//
  for ( node = 0; node < node_num; node++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_normal[i+node*3] = 0.0;
    }
  }
//
//  For each element:
//    Start with the last, first and second vertex.
//
  for ( element = 0; element < element_num; element++ )
  {
    vert0 = element_order - 2;
    vert1 = element_order - 1;
    for ( vert2 = 0; vert2 < element_order; vert2++ )
    {
      n0 = element_node[vert0+element*element_order];
      n1 = element_node[vert1+element*element_order];
      n2 = element_node[vert2+element*element_order];
//
//  Determine a unit normal vector associated with the plane of the triangle.
//
      for ( i = 0; i < 3; i++ )
      {
        v1[i] = node_xyz[i+n1*3] - node_xyz[i+n0*3];
      }
      for ( i = 0; i < 3; i++ )
      {
        v2[i] = node_xyz[i+n2*3] - node_xyz[i+n0*3];
      }

      triangle_normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
      triangle_normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
      triangle_normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

      norm = 0.0;
      for ( i = 0; i < 3; i++ )
      {
        norm = norm + pow ( triangle_normal[i], 2 );
      }
      norm = sqrt ( norm );

      for ( i = 0; i < 3; i++ )
      {
        triangle_normal[i] = triangle_normal[i] / norm;
      }

      for ( i = 0; i < 3; i++ )
      {
        node_normal[i+n1*3] = node_normal[i+n1*3] + triangle_normal[i];
      }
      vert0 = vert1;
      vert1 = vert2;
    }
  }
//
//  Renormalize.
//
  for ( node = 0; node < node_num; node++ )
  {
    norm = 0.0;
    for ( i = 0; i < 3; i++ )
    {
      norm = norm + pow ( node_normal[i+node*3], 2 );
    }

    if ( norm == 0.0 )
    {
      norm = 3.0;
      for ( i = 0; i < 3; i++ )
      {
        node_normal[i+node*3] = 1.0;
      }
    }

    norm = ( double ) sqrt ( norm );

    for ( i = 0; i < 3; i++ )
    {
      node_normal[i+node*3] = node_normal[i+node*3] / norm;
    }
  }

  return node_normal;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
//****************************************************************************80

double *r83vec_max ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83VEC_MAX returns the maximum value in an R83VEC.
//
//  Discussion:
//
//    An R83VEC is an array of triples of double precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[3*N], the array.
//
//    Output, double R83VEC_MAX[3]; the largest entries in each row.
//
{
# define DIM_NUM 3

  double *amax = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = new double[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amax[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[0+j*DIM_NUM] )
      {
        amax[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amax;
# undef DIM_NUM
}
//****************************************************************************80

double *r83vec_min ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83VEC_MIN returns the minimum value in an R83VEC.
//
//  Discussion:
//
//    An R83VEC is an array of triples of double precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[3*N], the array.
//
//    Output, double R83VEC_MIN[3]; the smallest entries in each row.
//
{
# define DIM_NUM 3

  double *amin = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = new double[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amin[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( a[0+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
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
//    10 September 2009
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
      cout << setw(5) << j << " ";
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
    if ( s[n-1] != ' ' )
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
//****************************************************************************80

void spinSurface ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPINSURFACE adjusts the angle of rotation and redisplays the picture.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2008
//
//  Author:
//
//    Edward Angel
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  theta[axis] = theta[axis] + theta_speed;
  if ( 360.0 < theta[axis] )
  {
    theta[axis] = theta[axis] - 360.0;
  }
  glutPostRedisplay ( );

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
