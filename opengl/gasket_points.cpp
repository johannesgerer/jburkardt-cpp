# include <cstdlib>
# include <iostream>

//
//  The location of the necessary include files depends on your
//  local implementation.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

int main ( int argc, char *argv[] );
void display ( );
void myinit ( );

int point_num = 50000;

typedef GLfloat point2[2];

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
//    This program draws the Sierpinski gasket by displaying points.
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
  cout << "GASKET_POINTS:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display the image of a Sierpinski gasket.\n";
  cout << "\n";
  cout << "  The gasket is depicted using " << point_num << " points.\n";

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Sierpinski Gasket (Points)" );
  glutDisplayFunc ( display );

  myinit ( );

  glutMainLoop ( );

  cout << "\n";
  cout << "GASKET_POINTS:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void display ( )

//****************************************************************************80
//
//  Purpose:
//
//    DISPLAY generates the graphics output.
//
//  Modified:
//
//    27 August 2003
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
  int k;
  point2 p = { 75.0, 50.0 };
  point2 vertices[3] = { 
    {0.0, 0.0}, { 250.0, 500.0}, { 500.0, 0.0} };
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Compute and plot the points.
//
  for ( k = 0; k < point_num; k++ )
  {
    j = rand ( ) % 3;
//
//  Compute point halfway between current point and randomly selected vertex.
//
    p[0] = 0.5 * ( p[0] + vertices[j][0] );
    p[1] = 0.5 * ( p[1] + vertices[j][1] );
//
//  Plot the new point.
//
    glBegin ( GL_POINTS );
      glVertex2fv ( p );
    glEnd ( );
  }
//
//  Clear all buffers.
//
  glFlush ( );

  return;
}
//****************************************************************************80

void myinit ( )

//****************************************************************************80
//
//  Purpose:
//
//    MYINIT initializes OpenGL state variables dealing with viewing and attributes.
//
//  Modified:
//
//    27 August 2003
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
//  Draw in RED.
//
  glColor3f ( 1.0, 0.0, 0.0 );
//
//  Set up 500 by 500 viewing window with origin at the lower left.
//
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
  gluOrtho2D ( 0.0, 500.0, 0.0, 500.0 );
  glMatrixMode ( GL_MODELVIEW );

  return;
}
