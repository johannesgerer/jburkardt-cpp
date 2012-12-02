# include <cstdlib>
# include <cmath>
# include <iostream>

//
//  This is the include statement I need on for Mac OS X.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

//
//  Get P, Q, A, B from the user, store them here, make them available to DISPLAY.
//
int a;
int b;
int p;
int q;

int main ( int argc, char *argv[] );
void display ( );
void myinit ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for the CAUSTIC_OPENGL program.
//
//  Discussion:
//
//    This program displays a caustic.  It draws Q points on a curve,
//    and then connects each point I to point I+P.
//
//    The curve is traced out by the values 
//
//      X(I) = cos ( A * 2 * I * PI / Q )
//      Y(I) = sin ( B * 2 * I * PI / Q ) 
//
//    for I = 0 to Q-1.
//
//    Note that if if A = B, the curve will be a circle.  
//
//    The values of Q, P, A and B are input from the user.
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
//    30 November 2008
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
  char title[80];

  cout << "\n";
  cout << "CAUSTIC_OPENGL:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display a caustic.\n";

  if ( 2 <= argc )
  {
    q = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "CAUSTIC_OPENGL:\n";
    cout << "  Please enter Q, the number of points on the circle.\n";
    cout << "  A reasonable number might be 100.\n";
    cin >> q;
  }

  if ( 3 <= argc )
  {
    p = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "CAUSTIC_OPENGL:\n";
    cout << "  Please enter P, the point to which point 0 is connected.\n";
    cout << "  Boring values are 0 and 1.\n";
    cin >> p;
  }

  if ( 4 <= argc )
  {
    a = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "CAUSTIC_OPENGL:\n";
    cout << "  Please enter A, the scale for x(i) = cos ( A * 2 * pi * i / Q ).\n";
    cout << "  A should be a small integer.\n";
    cin >> a;
  }

  if ( 5 <= argc )
  {
    b = atoi ( argv[4] );
  }
  else
  {
    cout << "\n";
    cout << "CAUSTIC_OPENGL:\n";
    cout << "  Please enter B, the scale for y(i) = sin ( B * 2 * pi * i / Q ).\n";
    cout << "  B should be a small integer.\n";
    cin >> b;
  }

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  sprintf ( title, "Caustic  %d  %d  %d  %d", q, p, a, b );
  glutCreateWindow ( title );
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
//    The user has specified:
//
//      Q, the number of points on the curve,
//      P, the index of the point to which point 0 is connected,
//      A, which defines x(i) = cos ( a * 2 * pi * i / q );
//      B, which defines y(i) = sin ( b * 2 * pi * i / q );
//
//    To form the caustic, the program draws the points [ X(i), Y(i) ]
//    for I = 0 to Q-1, and connects each point I to point I+P, using
//    modular arithmetic where necessary.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 November 2008
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
  int k;
  float pi = 3.14159265358979;
  float r = 1.0;
  float theta;
  float *xy;
//
//  Define the points.
//
  xy = new float[2*q];

  k = 0;
  for ( i = 0; i < q; i++ )
  {
    theta = ( float ) ( i * 2 ) * pi / ( float ) q;
    xy[k]   = r * cos ( a * theta );
    xy[k+1] = r * sin ( b * theta );
    k = k + 2;
  }
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Draw the points in BLUE.
//
  glColor3f ( 0.0, 0.0, 1.0 );

  for ( i = 0; i < q; i++ )
  {
    glBegin ( GL_POINTS );
      glVertex2fv ( xy+i*2 );
    glEnd ( );
  }
//
//  Draw the boundary lines, in GREEN.
//
  glColor3f ( 0.0, 1.0, 0.0 );

  for ( i = 0; i < q; i++ )
  {
    glBegin ( GL_LINES );

      glVertex2fv ( xy+i*2 );

      j = ( i + 1 ) % q;

      glVertex2fv ( xy+j*2 );

    glEnd ( );
  }
//
//  Draw the caustic lines, in RED.
//
  glColor3f ( 1.0, 0.0, 0.0 );

  for ( i = 0; i < q; i++ )
  {
    glBegin ( GL_LINES );

      glVertex2fv ( xy+i*2 );

      j = ( i + p ) % q;

      glVertex2fv ( xy+j*2 );

    glEnd ( );
  }
//
//  Clear all the buffers.
//
  glFlush ( );

  delete [] xy;

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
//    30 November 2008
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
//  Our generator points have X and Y values between -1 and 1.
//  We show a little wider range as a grace margin.
//
  gluOrtho2D ( -1.1, 1.1, -1.1, 1.1 );

  glMatrixMode ( GL_MODELVIEW );

  return;
}
