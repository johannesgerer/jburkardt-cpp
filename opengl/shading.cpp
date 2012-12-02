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

# include "graphicsDefs.hpp"
# include "axis.hpp"
# include "bird.hpp"
# include "gourd.hpp"
//
//  Global variables:
//
GLfloat light_position[] = { -1.0, 0.5, 6.0, 0.0 };

GLfloat mat_ambient[] = { 0.8, 0.8, 0.8, 1.0 };
GLfloat mat_diffuse[] = { 0.9, 0.7, 0.1, 1.0 };
GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess[] = { 50.0 };

GLfloat mat_ground[] = { 0.1, 1.0, 0.2, 1.0 };

float ang = 0.0;
int coloring_mode = 0;

int main ( int argc, char *argv[] );
void display ( );
void idle ( );
void drawSolidObject ( point3 vert[], int tris[][3], int numtris, 
  point3 norms[] );
void menu ( int );
void myInit ( );
void myReshape ( int, int );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    Allow user to select different types of shading for an object.
//
//  Discussion:
//
//    A right mouse click will produce a menu for selection of shading styles.
//
//  Modified:
//
//    01 June 2007
//
//  Author:
//
//    Margaret Geroch
//
{
  cout << "\n";
  cout << "SHADING:\n";
  cout << "  C++ version\n";

  glutInit ( &argc, argv );

  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize ( 500, 500 );
  glutCreateWindow ( "Various shading styles" );
  glutReshapeFunc ( myReshape );
  glutDisplayFunc ( display );
  glutIdleFunc ( idle );

  glutCreateMenu ( menu );
  glutAddMenuEntry ( "No shading", 0 );
  glutAddMenuEntry ( "Flat shading", 1 );
  glutAddMenuEntry ( "Gouraud shading", 2 );
  glutAddMenuEntry ( "Quit", 3 );
  glutAttachMenu ( GLUT_RIGHT_BUTTON );  

  myInit ( );

  glutMainLoop ( );
//
//  Terminate.
//
  return 0;
}
//****************************************************************************80

void display ( )

//****************************************************************************80
//
//  Purpose:
//
//    DISPLAY displays whenever the window is drawn.
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    Margaret Geroch
//
{
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode ( GL_MODELVIEW );
  glLoadIdentity ( );
//
//  First (last) push everything backwards into the projection volume.
//
  glTranslatef ( 0.0, 0.0, -5.0 );

  gluLookAt ( 0.1, 0.2, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 );
  glLightfv ( GL_LIGHT0, GL_POSITION, light_position );
//
//  Save the view matrix.
//
  glPushMatrix ( );
//
//  Objects:
//
  switch ( coloring_mode )
  {
    case 0:
      glColor4fv ( mat_diffuse );
      glDisable ( GL_LIGHTING );
      glDisable ( GL_LIGHT0 );
    break;
    case 1:
      glEnable ( GL_LIGHTING );
      glEnable ( GL_LIGHT0 );
      glShadeModel ( GL_FLAT );
      glMaterialfv ( GL_FRONT, GL_AMBIENT, mat_ambient );
      glMaterialfv ( GL_FRONT, GL_DIFFUSE, mat_diffuse );
      glMaterialfv ( GL_FRONT, GL_SPECULAR, mat_specular );
      glMaterialfv ( GL_FRONT, GL_SHININESS, mat_shininess );
    break;
    case 2:
      glEnable ( GL_LIGHTING );
      glEnable ( GL_LIGHT0 );
      glShadeModel ( GL_SMOOTH );
      glMaterialfv ( GL_FRONT, GL_AMBIENT, mat_ambient );
      glMaterialfv ( GL_FRONT, GL_DIFFUSE, mat_diffuse );
      glMaterialfv ( GL_FRONT, GL_SPECULAR, mat_specular );
      glMaterialfv ( GL_FRONT, GL_SHININESS, mat_shininess );
    break;
  }

  glPushMatrix ( );
  glTranslatef ( 0.7, 0.5, -0.8 );
  glRotatef ( ang, 0.0, 1.0, 0.0 );
  glutSolidTeapot ( 1.0 );
  glPopMatrix ( );

  glPushMatrix ( );
  glTranslatef ( -1.7, -0.5, 0.8 );
  glRotatef ( ang, 0.0, 1.0, 0.0 );
  glScalef ( 2.0, 2.0, 2.0 );
  glEnable ( GL_NORMALIZE );
  drawSolidObject ( BirdVerts, BirdFaces, numBirdFaces, BirdVertNorms );
  glPopMatrix ( );

  glPushMatrix ( );
  glTranslatef ( 0.0, -0.5, 0.0 );
  glRotatef ( ang, 0.0, 1.0, 0.0 );
  drawSolidObject ( GourdVerts, GourdFaces, numGourdFaces, GourdVertNorms );
  glPopMatrix ( );

  glutSwapBuffers ( );
}
//****************************************************************************80

void drawSolidObject ( point3 verts[], int tris[][3], int numtris, 
  point3 norms[] )

//****************************************************************************80
//
//  Purpose:
//
//    DRAWSOLIDOBJECT draws a triangulated solid object.
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    Margaret Geroch
//
{
  int j;

  glBegin ( GL_TRIANGLES );

  for ( j = 0; j < numtris; ++j )
  {
    glNormal3fv ( norms[tris[j][0]] );
    glVertex3fv ( verts[tris[j][0]] );
    
    glNormal3fv ( norms[tris[j][1]] );
    glVertex3fv ( verts[tris[j][1]] );

    glNormal3fv ( norms[tris[j][2]] );
    glVertex3fv ( verts[tris[j][2]] );
  }
  glEnd( );
}
//****************************************************************************80

void idle ( void )

//****************************************************************************80
//
//  Purpose:
//
//    IDLE specifies what to do when idle.
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    Margaret Geroch
//
{
  if ( 360.0 <= ang )
  {
    ang = 0.0;
  }
  else
  {
    ang = ang + 0.5;
  }

  glutPostRedisplay ( );
}
//****************************************************************************80

void menu ( int id )

//****************************************************************************80
//
//  Purpose:
//
//    MENU responds to user choice on the menu.
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    Margaret Geroch
//
{
  if ( id == 3 )
  {
    exit ( 0 );
  }
  else
  {
    coloring_mode = id;
    glutPostRedisplay ( );
  }
}
//****************************************************************************80

void myInit ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MYINIT initializes the world and camera view.
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    Margaret Geroch
//
{
//
//  Values for a default white light and light cyan material:
//

  GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };

  ang = 0.0;
//
//  Set up the default light and material.
//

//
//  glShadeModel ( GL_SMOOTH );
//
  glEnable ( GL_DEPTH_TEST );

  glLightfv ( GL_LIGHT0, GL_AMBIENT, light_ambient );
  glLightfv ( GL_LIGHT0, GL_DIFFUSE, light_diffuse );
  glLightfv ( GL_LIGHT0, GL_SPECULAR, light_specular );

  glEnable ( GL_LIGHTING );
  glEnable ( GL_LIGHT0 );
//
//  Set up the perspective frustrum.
//
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
  gluPerspective ( 60, 1, 2, 9 );
  
  glMatrixMode ( GL_MODELVIEW );
}
//****************************************************************************80

void myReshape ( int w, int h )

//****************************************************************************80
//
//  Purpose:
//
//    RESHAPE recalculates the world view if the window is reshaped.
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    Margaret Geroch
//
{
  glViewport ( 0, 0, w, h );

  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
  gluPerspective ( 60, w/h, 2, 9 );
  glMatrixMode ( GL_MODELVIEW );
  glutPostRedisplay ( );
}
