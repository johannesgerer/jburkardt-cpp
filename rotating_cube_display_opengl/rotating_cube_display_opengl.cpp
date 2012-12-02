# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <fstream>

//
//  This is the include statement I need for Mac OS X.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

GLfloat vertices[][3] = {
  { -1.0, -1.0, -1.0 },
  {  1.0, -1.0, -1.0 },
  {  1.0,  1.0, -1.0 },
  { -1.0,  1.0, -1.0 },
  { -1.0, -1.0,  1.0 },
  {  1.0, -1.0,  1.0 },
  {  1.0,  1.0,  1.0 },
  { -1.0,  1.0,  1.0 } };

GLfloat normals[][3] = {
  { -1.0, -1.0, -1.0 },
  {  1.0, -1.0, -1.0 },
  {  1.0,  1.0, -1.0 },
  { -1.0,  1.0, -1.0 },
  { -1.0, -1.0,  1.0 },
  {  1.0, -1.0,  1.0 },
  {  1.0,  1.0,  1.0 },
  { -1.0,  1.0,  1.0 } };

GLfloat colors[][3] = {
  { 0.0, 0.0, 0.0 },
  { 1.0, 0.0, 0.0 },
  { 1.0, 1.0, 0.0 },
  { 0.0, 1.0, 0.0 },
  { 0.0, 0.0, 1.0 },
  { 1.0, 0.0, 1.0 },
  { 1.0, 1.0, 1.0 },
  { 0.0, 1.0, 1.0 } };

static GLint axis = 2;
static GLfloat theta[3] = { 0.0, 0.0, 0.0 };

int main ( int argc, char *argv[] );
void colorcube ( );
void display ( );
void mouse ( int btn, int state, int x, int y );
void myReshape ( int w, int h );
void polygon ( int a, int b, int c, int d );
void spinCube ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ROTATING_CUBE_DISPLAY_OPENGL.
//
//  Discussion:
//
//    This program constructs a cube, each of whose vertices is given a 
//    different color, and displays the cube.  The cube rotates slowly
//    about the X, Y or Z axis.  Each time the user clicks the mouse, the
//    "next" axis is used for rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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
  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Rotating cube" );
  glutReshapeFunc ( myReshape );
  glutDisplayFunc ( display );
  glutIdleFunc ( spinCube );
  glutMouseFunc ( mouse );
//
//  Enable hidden surface removal.
//
  glEnable ( GL_DEPTH_TEST );
  glutMainLoop ( );
//
//  Terminate.
//
  return 0;
}
//****************************************************************************80

void colorcube ( )

//****************************************************************************80
//
//  Purpose:
//
//    COLORCUBE defines the 6 faces of the color cube object.
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
  polygon ( 0, 3, 2, 1 );
  polygon ( 2, 3, 7, 6 );
  polygon ( 0, 4, 7, 3 );
  polygon ( 1, 2, 6, 5 );
  polygon ( 4, 5, 6, 7 );
  polygon ( 0, 1, 5, 4 );

  return;
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

//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glLoadIdentity ( );
  
  glRotatef ( theta[0], 1.0, 0.0, 0.0 );
  glRotatef ( theta[1], 0.0, 1.0, 0.0 );
  glRotatef ( theta[2], 0.0, 0.0, 1.0 );

  colorcube ( );
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
//    dedicated one axis to each.  Since Apple prefers the esthetics of a
//    one button mouse, this routine simply increments the axis by 1,
//    no matter which button is pushed.
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
  if ( btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
  {
    axis = axis + 1;
  }
  if ( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN )
  {
    axis = axis + 1;
  }
  if ( btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN )
  {
    axis = axis + 1;
  }
  axis = axis % 3;

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
  glViewport ( 0, 0, w, h );
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );

  if ( w <= h )
  {
    glOrtho ( 
      -2.0, 2.0, 
      -2.0 * ( GLfloat ) h / ( GLfloat ) w, 2.0 * ( GLfloat ) h / ( GLfloat ) w, 
      -10.0, 10.0 );
  }
  else
  {
    glOrtho ( 
      -2.0 * ( GLfloat ) h / ( GLfloat ) w, 2.0 * ( GLfloat ) h / ( GLfloat ) w,  
      -2.0, 2.0, 
      -10.0, 10.0 );
  }

  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

void polygon ( int a, int b, int c, int d )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON defines the colors, vertices and normals for a quadrilateral.
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
  glBegin ( GL_POLYGON );

  glColor3fv ( colors[a] );
  glNormal3fv ( normals[a] );
  glVertex3fv ( vertices[a] );

  glColor3fv ( colors[b] );
  glNormal3fv ( normals[b] );
  glVertex3fv ( vertices[b] );

  glColor3fv ( colors[c] );
  glNormal3fv ( normals[c] );
  glVertex3fv ( vertices[c] );

  glColor3fv ( colors[d] );
  glNormal3fv ( normals[d] );
  glVertex3fv ( vertices[d] );

  glEnd ( );

  return;
}
//****************************************************************************80

void spinCube ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPINCUBE adjusts the angle of rotation and redisplays the picture.
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
  theta[axis] = theta[axis] + 0.020;
  if ( 360.0 < theta[axis] ) 
  {
    theta[axis] = theta[axis] - 360.0;
  }
  glutPostRedisplay ( );

  return;
}
