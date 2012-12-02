# include <cstdlib>
# include <ctime>
# include <cmath>
# include <fstream>
# include <iostream>
# include <iomanip>

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
void mouse ( int btn, int state, int x, int y );
void myinit ( );
void myReshape ( int w, int h );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r83vec_max ( int n, double a[] );
double *r83vec_min ( int n, double a[] );
int s_len_trim ( string s );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
void spin_image ( );
void timestamp ( );
void xyz_data_print ( int point_num, double xyz[] );
void xyz_data_read ( string input_filename, int point_num, double xyz[] );
void xyz_header_print ( int point_num );
void xyz_header_read ( string input_filename, int *point_num );
//
//  Global data.
//
  static GLint axis = 2;
  int click_num = 0;
  int dim_num = 3;
  int pixel_height;
  int pixel_width;
  int point_num = 0;
  bool spinning = true;
  static GLfloat theta[3] = { 0.0, 0.0, 0.0 };
  double theta_speed = 0.020;
  double *xyz = NULL;
  double xyz_center[3];
  double *xyz_max = NULL;
  double *xyz_min = NULL;
  double xyz_range[3];
  double xyz_scale;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_XYZ_DISPLAY_OPENGL.
//
//  Discussion:
//
//    This program reads a text file of 3D point coordinates and displays an
//    OpenGL image of the points.
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
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  int dim;
  string xyz_filename;
  int point;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "SPHERE_XYZ_DISPLAY_OPENGL:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display a set of 3D points.\n";
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "SPHERE_XYZ_DISPLAY_OPENGL:\n";
    cout << "  Please enter the name of a file of 3D point coordinates.\n";

    cin >> xyz_filename;
  }
  else 
  {
    xyz_filename = argv[1];
  }

  xyz_header_read ( xyz_filename, &point_num );

  cout << "\n";
  cout << "  The number of points POINT_NUM = " << point_num << "\n";

  xyz = new double[3*point_num];

  xyz_data_read ( xyz_filename, point_num, xyz );

  if ( false )
  {
    xyz_data_print ( point_num, xyz );
  }

  xyz_min = r83vec_min ( point_num, xyz );

  xyz_max = r83vec_max ( point_num, xyz );

  xyz_range[0] = xyz_max[0] - xyz_min[0];
  xyz_range[1] = xyz_max[1] - xyz_min[1];
  xyz_range[2] = xyz_max[2] - xyz_min[2];

  cout << "\n";
  cout << "  Minimum: " << xyz_min[0]
       << "  " << xyz_min[1]
       << "  " << xyz_min[2] << "\n";
  cout << "  Maximum: " << xyz_max[0] 
       << "  " << xyz_max[1]
       << "  " << xyz_max[2] << "\n";
  cout << "  Range:   " << xyz_range[0] 
       << "  " << xyz_range[1]
       << "  " << xyz_range[2] << "\n";

  if ( xyz_range[0] == 0.0 ) 
  {
    cout << "\n";
    cout << "SPHERE_XYZL_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The X data range is 0.\n";
    exit ( 1 );
  }

  if ( xyz_range[1] == 0.0 ) 
  {
    cout << "\n";
    cout << "SPHERE_XYZL_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The Y data range is 0.\n";
    exit ( 1 );
  }
  if ( xyz_range[2] == 0.0 ) 
  {
    cout << "\n";
    cout << "SPHERE_XYZL_DISPLAY_OPENGL - Fatal error!\n";
    cout << "  The Z data range is 0.\n";
    exit ( 1 );
  }

  if ( false )
  {
  xyz_scale = 0.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    xyz_center[dim] = ( xyz_min[dim] + xyz_max[dim] ) / 2.0;
    xyz_scale = r8_max ( xyz_scale, ( xyz_max[dim] - xyz_min[dim] ) / 2.0 );
  }
  xyz_scale = sqrt ( 3.0 ) * xyz_scale;
//
//  Translate the data so it is centered.
//  Scale the data so it fits in the unit cube.
//
  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      xyz[dim+point*dim_num] = ( xyz[dim+point*dim_num] - xyz_center[dim] ) / xyz_scale;
    }
  }

  }

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( xyz_filename.c_str() );
  glutReshapeFunc ( myReshape );
  glutDisplayFunc ( display );
  glutIdleFunc ( spin_image );
  glutMouseFunc ( mouse );
//
//  Enable hidden surface removal.
//
  glEnable ( GL_DEPTH_TEST );

  myinit ( );

  glutMainLoop ( );
//
//  Free memory.
//
  delete [] xyz;
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_XYZ_DISPLAY_OPENGL:\n";
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 December 2008
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
  float p[3];
  int point;
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glLoadIdentity ( );
  
  glRotatef ( theta[0], 1.0, 0.0, 0.0 );
  glRotatef ( theta[1], 0.0, 1.0, 0.0 );
  glRotatef ( theta[2], 0.0, 0.0, 1.0 );

  glColor3f ( 0.6, 0.6, 0.8 );

  glutSolidSphere ( 0.97, 20, 16 );
//
//  Draw the points in RED.
//  Note that OpenGL is using FLOAT's for real numbers, while we prefer DOUBLE's.
//
  glColor3f ( 1.0, 0.0, 0.0 );

  for ( point = 0; point < point_num; point++ )
  {
    glBegin ( GL_POINTS );

      p[0] = ( float ) xyz[0+point*dim_num];
      p[1] = ( float ) xyz[1+point*dim_num];
      p[2] = ( float ) xyz[2+point*dim_num];

      glVertex3fv ( p );

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
//    dedicated one axis to each.  
//
//    Since Apple prefers the esthetics of a one button mouse, we're forced
//    to live with that choice.  This routine alternately pauses rotation,
//    or increments the rotation axis by 1, no matter which button is pushed.
//
//  Modified:
//
//    30 December 2008
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
  float dist;
  float dist_min;
  int gen;
  int gen_min;
  int point;

  if ( ( btn == GLUT_LEFT_BUTTON   && state == GLUT_DOWN ) ||
       ( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN ) ||
       ( btn == GLUT_RIGHT_BUTTON  && state == GLUT_DOWN ) )
  {
    if ( spinning )
    {
      spinning = false;
      theta_speed = 0.0;
    }
    else
    {
      spinning = true;
      theta_speed = 0.020;
      axis = axis + 1;
    }

    click_num = click_num + 1;
  }

  axis = axis % 3;

  return;
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
//    30 December 2008
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

//
//  Set the background to WHITE.
//
  glClearColor ( 1.0, 1.0, 1.0, 1.0 );
//
//  The default point size is 1.0.
//
  if ( point_num <= 1000 )
  {
    glPointSize ( 5.0 );
  }
  else if ( point_num <= 4000 )
  {
    glPointSize ( 3.0 );
  }
  else
  {
    glPointSize ( 1.0 );
  }

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
//  Modified:
//
//    30 December 2008
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
      -1.05, +1.05, 
      - 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w, 
      -10.0, 10.0 );
  }
  else
  {
    glOrtho ( 
      - 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w,  
      - 1.05, + 1.05, 
      -10.0, 10.0 );
  }

  glMatrixMode ( GL_MODELVIEW );

  return;
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
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
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
//    04 January 2009
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
      if ( amax[i] < a[i+j*DIM_NUM] )
      {
        amax[i] = a[i+j*DIM_NUM];
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
//    04 January 2009
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
      if ( a[i+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[i+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
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

void spin_image ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPIN_IMAGE adjusts the angle of rotation and redisplays the picture.
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

void xyz_data_print ( int point_num, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_DATA_PRINT prints the data for an XYZ file.
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
//    Input, int POINT_NUM, the number of points.
//
//    Input, double XY[3*POINT_NUM], the arrays of coordinate data.
//
{
  int j;

  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << setw(10) << xyz[0+j*3] << "  "
         << setw(10) << xyz[1+j*3] << "  " 
         << setw(10) << xyz[2+j*3] << "\n";
    
  }
  return;
}
//****************************************************************************80

void xyz_data_read ( string input_filename, int point_num, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_DATA_READ reads the data in an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Output, double XYZ[3*POINT_NUM], the point coordinates.
//
{
  bool error;
  int i;
  ifstream input;
  int j;
  string text;
  double temp[3];

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZ_DATA_READ - Fatal error!\n";
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
    error = s_to_r8vec ( text, 3, temp );

    if ( error )
    {
      cout << "\n";
      cout << "XYZ_DATA_READ - Fatal error!\n";
      cout << "  S_TO_R8VEC returned error flag.\n";
      exit ( 1 );
    }

    xyz[0+j*3] = temp[0];
    xyz[1+j*3] = temp[1];
    xyz[2+j*3] = temp[2];
    j = j + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyz_header_print ( int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_HEADER_PRINT prints the header of an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
{
  cout << "\n";
  cout << "  Number of points = " << point_num << "\n";

  return;
}
//****************************************************************************80

void xyz_header_read ( string input_filename, int *point_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_HEADER_READ reads the header of an XYZ file.
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
//    01 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Output, int *POINT_NUM, the number of points.
//
{
  ifstream input;
  string text;

  *point_num = 0;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZ_HEADER_READ - Fatal error!\n";
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
