# include <cstdlib>
# include <ctime>
# include <cmath>
# include <fstream>
# include <iostream>

//
//  This is the include statement I need for Mac OS X.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

int main ( int argc, char *argv[] );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void display ( );
int i4vec_max ( int n, int a[] );
int i4vec_min ( int n, int a[] );
void myinit ( );
double r8_max ( double x, double y );
double *r82vec_max ( int n, double a[] );
double *r82vec_min ( int n, double a[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( );
void xy_data_read ( string input_filename, int point_num, double xy[] );
void xy_header_read ( string input_filename, int *point_num );
void xyf_data_read ( string input_filename, int face_num, int face_data_num,
  int face_pointer[], int face_data[] );
void xyf_header_read ( string input_filename, int *face_num, int *face_data_num );
//
//  Global data.
//
  int *face_data = NULL;
  int face_data_num;
  int face_num;
  int *face_pointer = NULL;
  int pixel_height;
  int pixel_width;
  int point_num = 0;
  double *xy = NULL;
  double *xy_max = NULL;
  double *xy_min = NULL;
  double xy_range[2];

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for XYF_DISPLAY_OPENGL.
//
//  Discussion:
//
//    This program reads a text file of point coordinates and another
//    text file of point indices forming faces,
//    and displays an OpenGL image of the points and faces between them.
//
//    The data is stored using the XYF format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
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
  int i;
  string prefix;
  string xy_filename;
  string xyf_filename;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "XYF_DISPLAY_OPENGL:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display a set of 2D points and faces.\n";
//
//  Assume the input files are paired, as in "FILE.XY" and "FILE.XYF".
//  Expect to receive the root name only.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "XYF_DISPLAY_OPENGL:\n";
    cout << "  Please enter the common \"first\" name of the XY and XYF files.\n";

    cin >> prefix;
  }
  else 
  {
    prefix = argv[1];
  }
//
//  Deal with the XY data.
//
  xy_filename = prefix + ".xy";

  xy_header_read ( xy_filename, &point_num );

  cout << "\n";
  cout << "  The number of points POINT_NUM = " << point_num << "\n";

  xy = new double[2*point_num];

  xy_data_read ( xy_filename, point_num, xy );

  xy_min = r82vec_min ( point_num, xy );
  xy_max = r82vec_max ( point_num, xy );

  xy_range[0] = xy_max[0] - xy_min[0];
  xy_range[1] = xy_max[1] - xy_min[1];

  cout << "\n";
  cout << "  Minimum: " << xy_min[0] << "  "
                        << xy_min[1] << "\n";
  cout << "  Maximum: " << xy_max[0] << "  "
                        << xy_max[1] << "\n";
  cout << "  Range:   " << xy_range[0] << "  "
                        << xy_range[1] << "\n";

  if ( xy_range[0] == 0.0 ) 
  {
    cout << "\n";
    cout << "XYF_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The X data range is 0.\n";
    exit ( 1 );
  }

  if ( xy_range[1] == 0.0 ) 
  {
    cout << "\n";
    cout << "XYF_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The Y data range is 0.\n";
    exit ( 1 );
  }
//
//  Deal with the XYF data.
//
  xyf_filename = prefix + ".xyf";

  xyf_header_read ( xyf_filename, &face_num, &face_data_num );

  face_pointer = new int[face_num+1];
  face_data = new int[face_data_num];

  xyf_data_read ( xyf_filename, face_num, face_data_num, face_pointer, face_data );
//
//  If the XYF data was created using 1-based indexing, we need to correct that,
//  if possible.
//
//  We can be sure that we are doing the right thing if:
//  A) the index 0 does not occur in FACE_DATA and
//  B) the index POINT_NUM does occur.
//
  if ( 0 < i4vec_min ( face_data_num, face_data ) &&
       point_num == i4vec_max ( face_data_num, face_data ) )
  {
    cout << "\n";
    cout << "  Converting face data from 1-based to 0-based indexing.\n";

    for ( i = 0; i < face_data_num; i++ )
    {
      face_data[i] = face_data[i] - 1;
    }
  }

  glutInit ( &argc, argv );

  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );

  if ( xy_range[1] < xy_range[0] )
  {
    pixel_width = 500;
    pixel_height = ( int ) 
      ( ( double ) ( 500 ) * xy_range[1] / xy_range[0] );
  }
  else
  {
    pixel_width = ( int ) 
      ( ( double ) ( 500 ) * xy_range[0] / xy_range[1] );
    pixel_height = 500;
  }
  cout << "  Pixels:  " << pixel_width  << "  "
                        << pixel_height << "\n";

  glutInitWindowSize ( pixel_width, pixel_height );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( prefix.c_str() );
  glutDisplayFunc ( display );

  myinit ( );

  glutMainLoop ( );
//
//  Free memory.
//
  delete [] face_data;
  delete [] face_pointer;
  delete [] xy;
//
//  Terminate.
//
  cout << "\n";
  cout << "XYF_DISPLAY_OPENGL:\n";
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
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the character was
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

void display ( void )

//****************************************************************************80
//
//  Purpose:
//
//    DISPLAY generates the graphics output.
//
//  Discussion;
//
//    We plot the points whose coordinates are in XY.
//
//    We also dump an ASCII PPM file of the screen.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2006
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
  int i;
  int j;
  int face;
  float p[2];
  int point;
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Draw faces in light GREEN.
//
  glColor3f ( 0.70, 1.0, 0.70 );

  for ( face = 0; face < face_num; face++ )
  {
    glBegin ( GL_POLYGON );

    for ( j = face_pointer[face]; j < face_pointer[face+1]; j++ )
    {
      point = face_data[j];
      p[0] = ( float ) xy[0+point*2];
      p[1] = ( float ) xy[1+point*2];
      glVertex2fv ( p );
    }
 
    glEnd ( );
  }
//
//  Draw lines in RED.
//  To close up the boundary, you must repeat the first point as the last one.
//
  glColor3f ( 1.0, 0.0, 0.0 );

  for ( face = 0; face < face_num; face++ )
  {
    glBegin ( GL_LINE_STRIP );

    for ( j = face_pointer[face]; j < face_pointer[face+1]; j++ )
    {
      point = face_data[j];
      p[0] = ( float ) xy[0+point*2];
      p[1] = ( float ) xy[1+point*2];
      glVertex2fv ( p );
    }

    j = face_pointer[face];
    point = face_data[j];
    p[0] = ( float ) xy[0+point*2];
    p[1] = ( float ) xy[1+point*2];
    glVertex2fv ( p );
 
    glEnd ( );
  }
//
//  Draw points in BLUE.
//
  glColor3f ( 0.0, 0.0, 1.0 );

  for ( point = 0; point < point_num; point++ )
  {
    glBegin ( GL_POINTS );

      p[0] = ( float ) xy[0+point*2];
      p[1] = ( float ) xy[1+point*2];

      glVertex2fv ( p );

    glEnd ( );
  }
//
//  Clear all the buffers.
//
  glFlush ( );

  return;
}
//****************************************************************************80

int i4vec_max ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAX returns the value of the maximum element in an I4VEC.
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
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MAX, the value of the maximum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;  
}
//****************************************************************************80

int i4vec_min ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the value of the minimum element in an I4VEC.
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
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MIN, the value of the minimum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value; 
}
//****************************************************************************80

void myinit ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MYINIT initializes OpenGL state variables dealing with viewing and attributes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2006
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
{
  double margin;
  float x_max;
  float x_min;
  float y_max;
  float y_min;
//
//  Set the background to WHITE.
//
  glClearColor ( 1.0, 1.0, 1.0, 1.0 );
//
//  The default point size is 1.0.
//
  if ( point_num <= 100 )
  {
    glPointSize ( 5.0 );
  }
  else if ( point_num <= 500 )
  {
    glPointSize ( 3.0 );
  }
  else
  {
    glPointSize ( 1.0 );
  }
//
//  Set up the viewing window with origin at the lower left.
//
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
//
//  Determine an amount MARGIN by which it would be appropriate to spread the
//  data range, so that all the data is comfortably inside the picture.
//
  margin = 0.025 * r8_max ( 
    xy_max[0] - xy_min[0],
    xy_max[1] - xy_min[1] );

  x_min = ( float ) ( xy_min[0] - margin );
  x_max = ( float ) ( xy_max[0] + margin );
  y_min = ( float ) ( xy_min[1] - margin );
  y_max = ( float ) ( xy_max[1] + margin );
//
//  Specify an orthogonal view.
//
  gluOrtho2D ( x_min, x_max, y_min, y_max );

  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8.
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

double *r82vec_max ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_MAX returns the maximum value in a R82VEC.
//
//  Discussion:
//
//    An R82VEC is an array of pairs of double precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[2*N], the array.
//
//    Output, double R82VEC_MAX[2]; the largest entries in each row.
//
{
  double *amax = NULL;
  int i;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = new double[2];

  amax[0] = a[0+0*2];
  amax[1] = a[1+0*2];

  for ( i = 1; i < n; i++ )
  {
    if ( amax[0] < a[0+i*2] )
    {
      amax[0] = a[0+i*2];
    }
    if ( amax[1] < a[1+i*2] )
    {
      amax[1] = a[1+i*2];
    }
  }
  return amax;
}
//****************************************************************************80

double *r82vec_min ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_MIN returns the minimum value in a R82VEC.
//
//  Discussion:
//
//    A R82VEC is an array of pairs of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[2*N], the array.
//
//    Output, double R82VEC_MIN[2]; the smallest entries in each row.
//
{
  double *amin = NULL;
  int i;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = new double[2];

  amin[0] = a[0+0*2];
  amin[1] = a[1+0*2];

  for ( i = 1; i < n; i++ )
  {
    if ( a[0+i*2] < amin[0] )
    {
      amin[0] = a[0+i*2];
    }
    if ( a[1+i*2] < amin[1] )
    {
      amin[1] = a[1+i*2];
    }
  }
  return amin;
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

void timestamp ( void )

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
//****************************************************************************80

void xy_data_read ( string input_filename, int point_num, double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    XY_DATA_READ reads the data in an XY file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Output, double XY[2*POINT_NUM], the point coordinates.
//
{
  bool error;
  int i;
  ifstream input;
  int j;
  string text;
  double temp[2];

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XY_DATA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  j = 0;

  while ( j < point_num )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }
    
    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }
//
//  Extract two real numbers.
//
    error = s_to_r8vec ( text, 2, temp );

    if ( error )
    {
      cout << "\n";
      cout << "XY_DATA_READ - Fatal error!\n";
      cout << "  S_TO_R8VEC returned error flag.\n";
      exit ( 1 );
    }

    xy[0+j*2] = temp[0];
    xy[1+j*2] = temp[1];
    j = j + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xy_header_read ( string input_filename, int *point_num )

//****************************************************************************80
//
//  Purpose:
//
//    XY_HEADER_READ reads the header of an XY file.
//
//  Discussion:
//
//    All we do here is count the number of lines that are not comments 
//    and not blank.  Each such line is assumed to represent a single point 
//    coordinate record.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *POINT_NUM, the number of points.
//
{
  ifstream input;
  string text;

  *point_num = 0;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XY_HEADER_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  while ( 1 )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    *point_num = *point_num + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyf_data_read ( string input_filename, int face_num, int face_data_num,
  int face_pointer[], int face_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYF_DATA_READ reads the data in an XYF file.
//
//  Discussion:
//
//    This routine assumes that the file contains exactly three kinds of
//    records:
//
//    COMMENTS which begin with a '#' character in column 1;
//    BLANKS which contain nothing but 'whitespace';
//    FACE ITEMS, which are indices of points on a face.
//
//    The routine ignores comments and blank faces and returns
//    the number of face items.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
//    Output, int FACE_POINTER[FACE_NUM+1], pointers to the
//    first face item for each face.
//
//    Output, int FACE_DATA[FACE_DATA_NUM], the face items.
//
{
  int ierror;
  int ilo;
  ifstream input;
  int n;
  int face;
  string text;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYF_DATA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  face = 0;
  face_pointer[0] = 0;

  while ( face < face_num )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      cout << "\n";
      cout << "XYF_DATA_READ - Fatal error!\n";
      cout << "  Unexpected end of information.\n";
      exit ( 1 );
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    n = s_word_count ( text );
    face_pointer[face+1] = face_pointer[face] + n;
 
    ilo = face_pointer[face];

    ierror = s_to_i4vec ( text, n, face_data+ilo );

    if ( ierror != 0 )
    {
      cout << "\n";
      cout << "XYF_DATA_READ - Fatal error!\n";
      cout << "  Error from S_TO_I4VEC.\n";
      exit ( 1 );
    }
    face = face + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyf_header_read ( string input_filename, int *face_num, 
  int *face_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYF_HEADER_READ determines the number of face items in an XYF file.
//
//  Discussion:
//
//    This routine assumes that the file contains exactly three kinds of
//    records:
//
//    COMMENTS which begin with a '#' character in column 1;
//    BLANKS which contain nothing but 'whitespace';
//    FACE ITEMS, which are indices of points on a face.
//
//    The routine ignores comments and blanks and returns
//    the number of face items.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *FACE_NUM, the number of faces.
//
//    Output, int *FACE_DATA_NUM, the number of face items.
//
{
  int i;
  int i4_val;
  int ierror;
  ifstream  input;
  int length;
  int n;
  string text;

  *face_data_num = 0;
  *face_num = 0;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYF_HEADER_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    n = s_word_count ( text );

    *face_data_num = *face_data_num + n;

    *face_num = *face_num + 1;
  }

  input.close ( );

  return;
}
