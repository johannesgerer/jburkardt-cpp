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
//  User input values are placed here.
//
float a1;
float a2;
float b1;
float b2;
int n;
float r1 = 1;
float r2 = 1;
float t2;

int main ( int argc, char *argv[] );
void display ( );
void myinit ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for the LISSAJOUS_OPENGL program.
//
//  Discussion:
//
//    This program displays a Lissajous figure.
//
//    The curve is traced out by N values of 
//
//      x(I) = R1 * sin ( ( A1 * t(i) + B1 ) * 2 * pi )
//      y(I) = R2 * sin ( ( A2 * t(i) + B2 ) * 2 * pi ) 
//
//    for 0 = T1 <= t(i) <= T2.
//
//    For now, the values R1 and R2 are set to 1.
//
//    The user is allowed to enter values for N, T2, A1, B1, A2, B2.
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
//    01 December 2008
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
  cout << "LISSAJOUS_OPENGL:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display a Lissajous figure of the form:\n";
  cout << "\n";
  cout << "    x(i) = sin ( ( A1 * t + B1 ) * 2 * pi ).\n";
  cout << "    y(i) = sin ( ( A2 * t + B2 ) * 2 * pi ).\n";

  if ( 2 <= argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "LISSAJOUS_OPENGL:\n";
    cout << "  Please enter N, the number of points.\n";
    cout << "  A reasonable number might be 1000.\n";
    cin >> n;
  }

  if ( 3 <= argc )
  {
    t2 = atof ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "LISSAJOUS_OPENGL:\n";
    cout << "  Please enter T2, the final value of time.\n";
    cout << "  A value of 1 is  reasonable if A1 and A2 are integers.\n";
    cin >> t2;
  }

  if ( 4 <= argc )
  {
    a1 = atof ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "LISSAJOUS_OPENGL:\n";
    cout << "  Please enter A1, for x(i) = sin ( ( A1 * t + B1 ) * 2 * pi ).\n";
    cin >> a1;
  }

  if ( 5 <= argc )
  {
    b1 = atof ( argv[4] );
  }
  else
  {
    cout << "\n";
    cout << "LISSAJOUS_OPENGL:\n";
    cout << "  Please enter B1, for x(i) = sin ( ( A1 * t + B1 ) * 2 * pi ).\n";
    cin >> b1;
  }

  if ( 6 <= argc )
  {
    a2 = atof ( argv[5] );
  }
  else
  {
    cout << "\n";
    cout << "LISSAJOUS_OPENGL:\n";
    cout << "  Please enter A2, for y(i) = sin ( ( A2 * t + B2 ) * 2 * pi ).\n";
    cin >> a2;
  }

  if ( 7 <= argc )
  {
    b2 = atof ( argv[6] );
  }
  else
  {
    cout << "\n";
    cout << "LISSAJOUS_OPENGL:\n";
    cout << "  Please enter B2, for y(i) = sin ( ( A2 * t + B2 ) * 2 * pi ).\n";
    cin >> b2;
  }

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  sprintf ( title, "Lissajous  %d  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f", 
    n, t2, a1, b1, a2, b2 );
  glutCreateWindow ( title );
  glutDisplayFunc ( display );

  myinit ( );

  glutMainLoop ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LISSAJOUS_OPENGL:\n";
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
//  Discussion;
//
//    The user has specified:
//
//      N, the number of points on the curve,
//      T2, the final value of T,
//      A1 and B1, for the X definition, and
//      A2 and B2, for the Y definition, which define
//
//        t(i) = i * t2 / ( n - 1 );
//        x(i) = sin ( ( a1 * t(i) + b1 ) * 2 * pi );
//        y(i) = sin ( ( a2 * t(i) + b2 ) * 2 * pi );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2008
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
  float t;
  float *xy;
//
//  Define the points.
//
  xy = new float[2*n];

  k = 0;
  for ( i = 0; i < n; i++ )
  {
    t = ( float ) ( i ) * t2 / ( float ) ( n - 1 );
    xy[k]   = r1 * sin ( ( a1 * t + b1 ) * 2.0 * pi );
    xy[k+1] = r2 * sin ( ( a2 * t + b2 ) * 2.0 * pi );
    k = k + 2;
  }
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Connect the points, in BLUE.
//
  glColor3f ( 0.0, 0.0, 1.0 );

  for ( i = 0; i < n - 1; i++ )
  {
    glBegin ( GL_LINES );

      glVertex2fv ( xy+i*2 );
      glVertex2fv ( xy+(i+1)*2 );

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
