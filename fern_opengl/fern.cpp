# include <cstdlib>
# include <iostream>

# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

int main ( int argc, char *argv[] );
void display ( );
void myinit ( );

typedef GLfloat point2[2];

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for the OpenGL Fern example.
//
//  Discussion:
//
//    This program draws the Barnsley fractal fern by plotting points.
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
//    09 May 2011
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
//    Michael Barnsley,
//    Fractals Everywhere,
//    Academic Press, 1988,
//    ISBN: 0120790696,
//    LC: QA614.86.B37.
//
//    Cleve Moler,
//    Experiments with MATLAB,
//    ebook: http://www.mathworks.com/moler/exm/index.html
//
{
  cout << "\n";
  cout << "FERN:\n";
  cout << "  C++ version\n";
  cout << "  This OpenGL program displays the Barnsley fractal fern.\n";

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( 400, 600 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Barnsley Fractal Fern" );
  glutDisplayFunc ( display );

  myinit ( );

  glutMainLoop ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FERN:\n";
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2011
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
//    Michael Barnsley,
//    Fractals Everywhere,
//    Academic Press, 1988,
//    ISBN: 0120790696,
//    LC: QA614.86.B37.
//
//    Cleve Moler,
//    Experiments with MATLAB,
//    ebook: http://www.mathworks.com/moler/exm/index.html
//
{
  int i;
  point2 p;
  int point_num = 500000;
  double prob[4] = { 0.85, 0.92, 0.99, 1.00 };
  double r;
  double x;
  double y;
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Compute and plot the points.
//
  p[0] = drand48 ( );
  p[1] = drand48 ( );

  for ( i = 0; i < point_num; i++ )
  {
    r = drand48 ( );

    if ( r < prob[0] )
    {
      x =   0.85 * p[0] + 0.04 * p[1] + 0.0;
      y = - 0.04 * p[0] + 0.85 * p[1] + 1.6;
    }
    else if ( r < prob[1] )
    {
      x =   0.20 * p[0] - 0.26 * p[1] + 0.0;
      y =   0.23 * p[0] + 0.22 * p[1] + 1.6;
    }
    else if ( r < prob[2] )
    {
      x = - 0.15 * p[0] + 0.28 * p[1] + 0.0;
      y =   0.26 * p[0] + 0.24 * p[1] + 0.44;
    }
    else
    {
      x =   0.00 * p[0] + 0.00 * p[1] + 0.0;
      y =   0.00 * p[0] + 0.16 * p[1] + 0.0;
    }

    p[0] = x;
    p[1] = y;
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2011
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
//    Michael Barnsley,
//    Fractals Everywhere,
//    Academic Press, 1988,
//    ISBN: 0120790696,
//    LC: QA614.86.B37.
//
{
//
//  Set the background to WHITE.
//
  glClearColor ( 1.0, 1.0, 1.0, 1.0 );
//
//  Draw in FOREST GREEN.
//
  glColor3f ( 0.133, 0.545, 0.133 );
//
//  Set up a viewing window with origin at the lower left.
//
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
  gluOrtho2D ( -4.0, 4.0, -1.0, 11.0 );
  glMatrixMode ( GL_MODELVIEW );

  return;
}
