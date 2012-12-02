# include <cstdlib>
# include <iostream>

# include <GLUT/glut.h>
//# include <GL/gl.h>
//# include <GL/glu.h>
//# include <GL/glx.h>
//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;
//
// Define a 2D point data type.
//

typedef struct { float x, y;} point_2d;

//
// Cosine and sine of 95 degrees.
//
float c = -0.087155;
float s =  0.996195;

float r = 0.0;
float g = 0.2;
float b = 0.0;

float window_xmin = 0.0;
float window_xmax = 600.0;

float window_ymin = 0.0;
float window_ymax = 400.0;

point_2d mid = { 250.0, 250.0 };
point_2d next = { 250.0, 250.0 };
point_2d old = { 250.0, 250.0 };

point_2d direction = { 1.0, 0.0 };

int main ( int argc, char* argv[] );
void clear ( void );
void display ( void );
void mouse ( int btn, int state, int x, int y );

//****************************************************************************80

int main ( int argc, char* argv[] )

//****************************************************************************80
// 
//  Purpose:
//
//    MAIN is the main routine for the TURTLE program.
//
//  Discussion:
//
//    This program illustrates a simple use of the mouse buttons to control
//    an OpenGL graphics program.
//
//    A line of a certain color is always being drawn on the screen.
//
//    Each time the LEFT button is clicked, the line turns left 95 degrees
//    and the line color is modified. 
//
//    A RIGHT button makes a right turn of 95 degrees and a color change. 
//
//    The MIDDLE button stops the program.
//
//  Modified:
//
//    23 June 2006
//
{
  cout << "\n";
  cout << "TURTLE\n";
  cout << "  C++ version\n";

  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
  glutInitWindowSize ( window_xmax - window_xmin, window_ymax - window_ymin );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Mouse Patrol" );

  glutIdleFunc ( display );
  glutMouseFunc ( mouse );  
  glClearColor ( 0.9, 0.9, 0.9, 0.0 );
  glColor3f ( r, g, b );

  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
  gluOrtho2D ( window_xmin, window_xmax, window_ymin, window_ymax );
  glMatrixMode ( GL_MODELVIEW );

  glutDisplayFunc ( clear );

  glutMainLoop();
}
//****************************************************************************80

void clear ( void )

//****************************************************************************80
//
//  Purpose:
//
//    CLEAR ???
//
//  Modified:
//
//    23 June 2006
//
{
  glClear ( GL_COLOR_BUFFER_BIT );
}
//****************************************************************************80

void display ( void )

//****************************************************************************80
// 
//  Purpose:
//
//    DISPLAY draws little line segments in the current direction.
//
//  Modified:
//
//    23 June 2006
//
{
  float delta = 0.01;
//
//  If the point goes out of the region, make it come in the other side.
//
  next.x = old.x + delta * direction.x; 
  next.y = old.y + delta * direction.y;
//
//  Take care of simple wraparound cases.
//
  if ( next.x < window_xmin )
  {
    mid.x = window_xmin;
    mid.y = old.y + ( ( mid.x - old.x ) / direction.x ) * direction.y;
    glBegin ( GL_LINES );
      glVertex2f ( old.x, old.y );
      glVertex2f ( mid.x, mid.y );
    glEnd ( );

    mid.x = window_xmax;
    next.x = next.x + window_xmax;
    glBegin ( GL_LINES );
      glVertex2f ( mid.x, mid.y );
      glVertex2f ( next.x, next.y );
    glEnd ( );

  }
  else if ( window_xmax < next.x )
  {
    mid.x = window_xmax;
    mid.y = old.y + ( ( mid.x - old.x ) / direction.x ) * direction.y;
    glBegin ( GL_LINES );
      glVertex2f ( old.x, old.y );
      glVertex2f ( mid.x, mid.y );
    glEnd ( );

    mid.x = window_xmin;
    next.x = next.x - window_xmax;
    glBegin ( GL_LINES );
      glVertex2f ( mid.x, mid.y );
      glVertex2f ( next.x, next.y );
    glEnd ( );
  }
  else if ( next.y < window_ymin )
  {
    mid.y = window_ymin;
    mid.x = old.x + ( ( mid.y - old.y ) / direction.y ) * direction.x;
    glBegin ( GL_LINES );
      glVertex2f ( old.x, old.y );
      glVertex2f ( mid.x, mid.y );
    glEnd ( );

    mid.y = window_ymax;
    next.y = next.y + window_ymax;
    glBegin ( GL_LINES );
      glVertex2f ( mid.x, mid.y );
      glVertex2f ( next.x, next.y );
    glEnd ( );

  }
  else if ( window_ymax < next.y )
  {
    mid.y = window_ymax;
    mid.x = old.x + ( ( mid.y - old.y ) / direction.y ) * direction.x;
    glBegin ( GL_LINES );
      glVertex2f ( old.x, old.y );
      glVertex2f ( mid.x, mid.y );
    glEnd ( );

    mid.y = window_ymin;
    next.y = next.y - window_ymax;
    glBegin ( GL_LINES );
      glVertex2f ( mid.x, mid.y );
      glVertex2f ( next.x, next.y );
    glEnd ( );
  }
  else
  { 
    glBegin ( GL_LINES );
      glVertex2f ( old.x, old.y );
      glVertex2f ( next.x, next.y );
    glEnd ( );
  }

  old.x = next.x;
  old.y = next.y;

  glFlush();
}
//****************************************************************************80

void mouse ( int btn, int state, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    MOUSE responds to the user mouse button input.
//
//  Modified:
//
//    23 June 2006
//
{
  float temp;

  if ( state == GLUT_DOWN )
  {
    if ( btn == GLUT_LEFT_BUTTON )
    {
      r = r / 2.0;
      b = ( 1.0 + b ) / 2.0;
      glColor3f ( r, g, b );

      temp =          c * direction.x - s * direction.y;
      direction.y = + s * direction.x + c * direction.y;
      direction.x = temp;

      glutIdleFunc ( display );
    }
    else if ( btn == GLUT_RIGHT_BUTTON )
    {
      r = ( 1.0 + r ) / 2.0;
      b = 2.0 * b;
      glColor3f ( r, g, b );

      temp =          c * direction.x + s * direction.y;
      direction.y = - s * direction.x + c * direction.y;
      direction.x = temp;

      glutIdleFunc ( display );
    }
    else if ( btn == GLUT_MIDDLE_BUTTON )
    {
      exit ( 0 );
    }
  }
}
