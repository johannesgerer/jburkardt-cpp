# include <cstdlib>
# include <iostream>

# include <GLUT/glut.h>
//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

typedef float point2[2];

int main ( int argc, char *argv[] );
void display ( void );
void divide_triangle ( point2 a, point2 b, point2 c, int m );
void myinit ( void );
void triangle ( point2 a, point2 b, point2 c );

//
//  This data needs to be shared by several routines.
//
point2 v[] = { 
  { -1.0,  0.58 },
  {  1.0, -0.58 },
  {  0.0,  1.15 } };

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
//    This program draws the Sierpinski gasket by displaying filled polygons.
//
//    The main program calls GLUT functions to set up the windows,
//    name the required callbacks and callback functions, in particular
//    the DISPLAY callback.
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
  cout << "\n";
  cout << "GASKET_POLY_FILLED:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display the image of a Sierpinski gasket.\n";
  cout << "\n";
  cout << "  The gasket is depicted using filled polygons.\n";

  if ( 2 <= argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "GASKET_POLY_FILLED:\n";
    cout << "  Please enter N, the number of recursive steps.\n";
    cout << "  A reasonable number if 4 or 5.\n";
    cin >> n;
  }

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Sierpinski Gasket (Filled Polygons)" );
  glutDisplayFunc ( display );

  myinit ( );

  glutMainLoop ( );

  cout << "\n";
  cout << "GASKET_POLY_FILLED:\n";
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
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );

  divide_triangle ( v[0], v[1], v[2], n );
//
//  Clear all buffers.
//
  glFlush ( );

  return;
}
//****************************************************************************80

void divide_triangle ( point2 a, point2 b, point2 c, int m )

//****************************************************************************80
//
//  Purpose:
//
//    DIVIDE_TRIANGLE...
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
  else
  {
    triangle ( a, b, c );
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
{
  glBegin ( GL_TRIANGLES );
    glVertex2fv ( a );
    glVertex2fv ( b );
    glVertex2fv ( c );
  glEnd ( );

  return;
}
