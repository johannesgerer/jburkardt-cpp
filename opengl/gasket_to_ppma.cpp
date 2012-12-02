# include <cstdlib>
# include <iostream>
# include <fstream>

# include <GLUT/glut.h>
//# include <GL/gl.h>
//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

typedef float point2[2];

int main ( int argc, char *argv[] );
void display ( void );
void divide_triangle ( point2 a, point2 b, point2 c, int m );
void i3mat_flip_cols ( int m, int n, int a[] );
void myinit ( void );
bool ppma_write ( char *file_out_name, int xsize, int ysize, int *rgb );
bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *rgb );
bool ppma_write_header ( ofstream &file_out, char *file_out_name, int xsize, 
  int ysize, int rgb_max );
void triangle ( point2 a, point2 b, point2 c );

//
//  This data needs to be shared by several routines.
//
point2 v[] = { 
  { -2.0, -0.8 },
  {  2.0, -0.8 },
  {  0.0,  1.5 } };

int n;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for the OpenGL Gasket example.
//
//  Discussion:
//
//    This program is especially interesting because it transfers the 
//    OpenGL pixels to an ASCII PPM file.
//
//    This program draws the Sierpinski gasket by displaying filled polygons.
//
//    The main program calls GLUT functions to set up the windows,
//    name the required callbacks and callback functions, in particular
//    the DISPLAY callback.
//
//  Modified:
//
//    21 September 2003
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
  cout << "\n";
  cout << "GASKET_TO_PPMA:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display the image of a Sierpinski gasket.\n";
  cout << "\n";
  cout << "  The gasket is depicted using 'open' (unfilled) polygons.\n";
  cout << "\n";
  cout << "  The OpenGL pixels are read from the screen into an array\n";
  cout << "  which is then written to an ASCII PPM file.\n";

  if ( 2 <= argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "GASKET_TO_PPMA:\n";
    cout << "  Please enter N, the number of recursive steps.\n";
    cout << "  A reasonable value is 4 or 5.\n";
    cin >> n;
  }

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Sierpinski Gasket (Unfilled Polygons)" );
  glutDisplayFunc ( display );

  myinit ( );

  glutMainLoop ( );

  cout << "\n";
  cout << "GASKET_TO_PPMA:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void display ( void )

//****************************************************************************80
//
//  Purpose:
//
//    DISPLAY generates the graphics output.
//
//  Modified:
//
//    21 June 2006
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
# define XPIXELS 500
# define YPIXELS 500

  char *file_out_name = "gasket_ascii.ppm";
  GLint *rgb;
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Carry out N levels of subdivision on the triangle.
//
  divide_triangle ( v[0], v[1], v[2], n );
//
//  Clear all buffers.
//
  glFlush ( );
//
//  Read the current graphics pixels to an array.
//
  rgb = new GLint[3*XPIXELS*YPIXELS];

  glReadPixels ( 0, 0, XPIXELS, YPIXELS, GL_RGB, GL_INT, rgb );
//
//  OpenGL and PPM disagree about the ordering of pixels in the vertical direction,
//  so we need to flip them.
//
  i3mat_flip_cols ( XPIXELS, YPIXELS, rgb );
//
//  Write the data to an ASCII PPM file.
//
  ppma_write ( file_out_name, XPIXELS, YPIXELS, rgb );

  delete [] rgb;
//
//  Clear all buffers.
//
  glFlush ( );

  return;
# undef XPIXELS
# undef YPIXELS
}
//****************************************************************************80

void divide_triangle ( point2 a, point2 b, point2 c, int m )

//****************************************************************************80
//
//  Purpose:
//
//    DIVIDE_TRIANGLE is called recursively to subdivide a triangle M times.
//
//  Modified:
//
//    08 September 2003
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
//  Parameters:
//
//    Input, POINT2 A, B, C, the corners of the current triangle.
//
//    Input, int M, counts the number of levels of recursion that remain to do.
//    If M = 0, then the triangle is ready to draw.
//
{
  int j;
  point2 v0;
  point2 v1;
  point2 v2;

  if ( 0 < m ) 
  {
    for ( j = 0; j < 2; j++ )
    {
      v0[j] = ( a[j] + b[j] ) / 2.0;
      v1[j] = ( a[j] + c[j] ) / 2.0;
      v2[j] = ( b[j] + c[j] ) / 2.0;
    }
    divide_triangle ( a, v0, v1, m-1 );
    divide_triangle ( c, v1, v2, m-1 );
    divide_triangle ( b, v2, v0, m-1 );
  }
//
//  If M = 0, then this triangle is the result of M subdivisions, and can now be
//  displayed.
//
  else
  {
    triangle ( a, b, c );
  }
  return;
}
//****************************************************************************80

void i3mat_flip_cols ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I3MAT_FLIP_COLS swaps the columns of an I3MAT.
//
//  Discussion:
//
//    An I3MAT is a matrix, each of whose entries is an I3, a triple of integers.
//
//    An I3MAT can be stored as a 3 x M x N array, where M counts the "columns"
//    and N counts the "rows".
//
//  Modified:
//
//    22 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[3*M*N], the matrix whose columns are to be flipped.
//
{
  int b;
  int i;
  int j;
  int k;
  
  for ( k = 0; k < ( n / 2 ); k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        b                    = a[i+j*3+     k *m*3];
        a[i+j*3+     k *m*3] = a[i+j*3+(n-1-k)*m*3];
        a[i+j*3+(n-1-k)*m*3] = b;
      }
    }
  }

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
//  Modified:
//
//    08 September 2003
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
//
//  Set the background to WHITE.
//
  glClearColor ( 1.0, 1.0, 1.0, 1.0 );
//
//  Draw in BLUE.
//
  glColor3f ( 0.0, 0.0, 1.0 );
//
//  Set up 500 by 500 viewing window with origin at the lower left.
//
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
  gluOrtho2D ( -2.0, 2.0, -2.0, 2.0);
  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

bool ppma_write ( char *file_out_name, int xsize, int ysize, int *rgb )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE writes the header and data for an ASCII portable pixel map file.
// 
//  Example:
//
//    P3
//    # feep.ppm
//    4 4
//    15
//     0  0  0    0  0  0    0  0  0   15  0 15
//     0  0  0    0 15  7    0  0  0    0  0  0
//     0  0  0    0  0  0    0 15  7    0  0  0
//    15  0 15    0  0  0    0  0  0    0  0  0
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
//    Input, char *FILE_OUT_NAME, the name of the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *RGBB, the array of size 3 by XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE, is
//    true, if an error was detected, or
//    false, if the file was written.
//
{
  bool error;
  ofstream file_out;
  int i;
  int j;
  int k;
  int *rgb_index;
  int rgb_max;
//
//  Open the output file.
//
  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << file_out_name << "\".\n";
    return true;
  }
//
//  Compute the maximum.
//
  rgb_max = 0;
  rgb_index = rgb;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      for ( k = 0; k < 3; k++ )
      {
        if ( rgb_max < *rgb_index )
        {
          rgb_max = *rgb_index;
        }
        rgb_index = rgb_index + 1;
      }
    }
  }
//
//  Write the header.
//
  error = ppma_write_header ( file_out, file_out_name, xsize, ysize, rgb_max );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_HEADER failed.\n";
    return true;
  }
//
//  Write the data.
//
  error = ppma_write_data ( file_out, xsize, ysize, rgb );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  file_out.close ( );

  return false;
}
//****************************************************************************80

bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *rgb )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE_DATA writes the data for an ASCII portable pixel map file.
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
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *RGB, the array of 3 by XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE_DATA, is
//    true, if an error was detected, or
//    false, if the data was written.
//
{
  int i;
  int j;
  int k;
  int *rgb_index;
  int rgb_num;

  rgb_index = rgb;
  rgb_num = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      for ( k = 0; k < 3; k++ )
      {
        file_out << *rgb_index << " ";
        rgb_index = rgb_index + 1;
      }
      if ( *rgb_index % 12 == 0 || i == xsize - 1 )
      {
        file_out << "\n";
      }
    }
  }
  return false;
}
//****************************************************************************80

bool ppma_write_header ( ofstream &file_out, char *file_out_name, int xsize, 
  int ysize, int rgb_max )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE_HEADER writes the header of an ASCII portable pixel map file.
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, char *FILE_OUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int RGB_MAX, the maximum RGB value.
//
//    Output, bool PPMA_WRITE_HEADER, is
//    true, if an error was detected, or
//    false, if the header was written.
//
{
  file_out << "P3\n";
  file_out << "# " << file_out_name 
           << " created by GASKET_TO_PPMA::PPMA_WRITE.C.\n";
  file_out << xsize << "  " << ysize << "\n";
  file_out << rgb_max << "\n";

  return false;
}
//****************************************************************************80

void triangle ( point2 a, point2 b, point2 c )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE displays one triangle.
//
//  Modified:
//
//    08 September 2003
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
//  Parameters:
//
//    Input, point2 A, B, C, the vertices of the triangle.
//
{
  glBegin ( GL_LINE_LOOP );
    glVertex2fv ( a );
    glVertex2fv ( b );
    glVertex2fv ( c );
  glEnd ( );

  return;
}
