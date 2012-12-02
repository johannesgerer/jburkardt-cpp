# include <cstdlib>
# include <iostream>

# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

int main ( int argc, char *argv[] );
void display ( void );
void myinit ( void );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for the OpenGL VORONOI_DISPLAY example.
//
//  Discussion:
//
//    This program displays a Voronoi diagram.  It shows a set of points in the
//    unit square, and the lines that form the "boundaries" of the Voronoi
//    regions.  Each region is made up of the area that is closest to a
//    particular point.  All of the calculations have been done already.
//    This program simply displays the results.
//
//    The main program calls GLUT functions to set up the windows,
//    name the required callbacks and callback functions, in particular
//    the DISPLAY callback.
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
  cout << "\n";
  cout << "VORONOI_DISPLAY:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display a set of points ('generators'),\n";
  cout << "  and their Voronoi regions.\n";

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Voronoi Diagram" );
  glutDisplayFunc ( display );

  myinit ( );

  glutMainLoop ( );

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
//  Discussion;
//
//    Here we have stored the results of a Voronoi diagram calculation.
//
//    This program is simply a demonstration.  In a real usage, the generators
//    would be read in or randomly generated, and the lines defining the boundaries
//    of the Voronoi regions would be computed as part of the computation.
//
//    Notice, also, that the semi-infinite rays that form part of the diagram
//    have not been included in the data.
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
# define GENERATOR_NUM 9
# define VERTEX_NUM 12
# define LINE_NUM 16

  int i;
  int j;
  int k;
//
//  Coordinates of the generators.
//
  float generator_xy[2*GENERATOR_NUM] = { 
    0.0, 0.0,
    0.0, 1.0,
    0.2, 0.5,
    0.3, 0.6,
    0.4, 0.5,
    0.6, 0.3,
    0.6, 0.5,
    1.0, 0.0,
    1.0, 1.0 };
//
//  Pairs of vertices that form the Voronoi boundary line segments.
//
  int line_vertex[2*LINE_NUM] = { 
     0,  2,
     0,  1,
     1,  4,
     1,  7,
     2,  3,
     2,  5,
     3,  4,
     3,  6,
     4,  8,
     5, 10,
     6,  8,
     6, 10,
     7,  9,
     8,  9,
     9, 11,
    10, 11 };
//
//  Coordinates of points that are endpoints of Voronoi boundary line segments.
//
  float vertex_xy[2*VERTEX_NUM] = { 
    -0.525,  0.500,
     0.064,  0.735,
     0.287,  0.175,
     0.300,  0.200,
     0.300,  0.500,
     0.500, -0.250,
     0.500,  0.400,
     0.500,  1.062,
     0.500,  0.700,
     0.576,  0.928,
     0.987,  0.400,
     1.112,  0.500 };
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Draw the generator points in BLUE.
//
  glColor3f ( 0.0, 0.0, 1.0 );

  for ( i = 0; i < GENERATOR_NUM; i++ )
  {
    glBegin ( GL_POINTS );
      glVertex2fv ( generator_xy+i*2 );
    glEnd ( );
  }
//
//  Draw the boundary lines, in RED.
//
  glColor3f ( 1.0, 0.0, 0.0 );

  for ( i = 0; i < LINE_NUM; i++ )
  {
    glBegin ( GL_LINES );

      j = line_vertex[0+i*2];
      glVertex2fv ( vertex_xy+j*2 );

      k = line_vertex[1+i*2];
      glVertex2fv ( vertex_xy+k*2 );

    glEnd ( );
  }
//
//  Clear all the buffers.
//
  glFlush ( );

  return;
# undef GENERATOR_NUM
# undef VERTEX_NUM
# undef LINE_NUM
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
//  Our generator points have X and Y values between 0 and 1.
//  We decide to show a little wider range, because some points are
//  right on the boundary of the unit box, and the boundary lines extend beyond that.
//
  gluOrtho2D ( -0.1, 1.1, -0.1, 1.1 );

  glMatrixMode ( GL_MODELVIEW );

  return;
}
