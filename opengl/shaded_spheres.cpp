# include <cstdlib>
# include <iostream>
# include <cmath>

# include <GLUT/glut.h>
//# include <GL/glx.h>
//# include <GL/gl.h>
//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
GLfloat light_diffuse[] = { 0.8, 0.7, 0.0, 1.0 };
GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };

int num_spheres;
float phi = 0.0;
float theta = 0.0;

int main ( int argc, char* argv[] );

void display ( void );
void init ( void );
void keyboard ( unsigned char key, int x, int y );
void mouse ( int btn, int state, int x, int y );
float r4_random ( float rlo, float rhi );
void reshape ( int w, int h );

//****************************************************************************80

int main ( int argc, char* argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    SHADED_SPHERES displays a set of spheres with shading.
//
//  Discussion:
//
//    Display a number of 3D spheres of random size and position.
//    Have the scene lighted from a given direction.
//    Allow the light source to move, based on user mouse input.
//    * Left mouse button increments PHI.
//    * Right mouse button increments THETA.
//
//  Modified:
//
//    04 February 2000
//
//  Author:
//
//    John Burkardt
//
{
//
//  Retrieve the desired number of spheres from the command line.
//  If no number is specified, use a random value between 1 and 10.
//
  cout << "\n";
  cout << "SHADED_SPHERES:\n";
  cout << "  C++ version\n";
  cout << "  Display a set of shaded spheres.\n";

  if ( argc <= 1 )
  {
    num_spheres = 1 + ( rand ( ) % 10 );
  }
  else
  {
    num_spheres = atoi ( argv[1] );
  }

  cout << "\n";
  cout << "  Using " << num_spheres << " spheres.\n";
//  
//  Open window with initial window size, title bar, 
//  RGBA display mode, and handle input events.
//
  glutInit ( &argc, argv );
//
//  The option GLUT_DEPTH enables depth buffering, which means that
//  hidden surfaces are properly hidden.
//
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH );

  glutInitWindowSize ( 500, 500 );

  glutCreateWindow ( "Shaded Spheres - John Burkardt\n" );

  init ( );

  glutReshapeFunc ( reshape );

  glutDisplayFunc ( display );

  glutMouseFunc ( mouse );
  
  glutKeyboardFunc ( keyboard );

  glutMainLoop ( );

  return 0; 
}
//****************************************************************************80

void display ( void )

//****************************************************************************80
//
//  Purpose:
//
//    DISPLAY is passed to the glutDisplayFunc routine.
//
//  Modified:
//
//    04 February 2000
//
{
  float angle;
  int i;
  int j;
  float r1;
  float r2;
  float x;
  float y;
  float z;
//
//  When using depth buffering, you have to clear the depth buffer bit
//  on every call to the DISPLAY function.
//
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
//
//  I'm being lazy here.  Every time I draw the picture, I'm calling the
//  random number generator to tell me what I draw and how big.  But I
//  really only want this to be random on the first call.  When I call
//  again, I want to draw the same picture (but with different lighting).
//  To force this to happen, I make the random number generator start
//  at the same spot.
//
  srand ( 1776 );

  glPushMatrix ();

  glRotatef ( 20.0, 1.0, 0.0, 0.0 );

  for ( i = 1; i <= num_spheres; i++ )
  {
    j = rand ( ) % 3;
    glPushMatrix ();

    if ( j == 0 )
    {
      angle = ( float ) ( rand ( ) % 360 );
      r1 = r4_random ( 0.125, 0.25 );
      r2 = r4_random ( 0.30, 0.6 );
      x = r4_random ( -2.5, 2.5 );
      y = r4_random ( -2.5, 2.5 );
      z = r4_random ( -2.5, 2.5 );
      glTranslatef ( x, y, z ); 
      glRotatef ( angle, 1.0, 0.0, 0.0 );
      glutSolidTorus ( r1, r2, 15, 15 );
    }
    else if ( j == 1 )
    {
      angle = (float) ( rand ( ) % 360 );
      r1 = r4_random ( 0.25, 1.0 );
      r2 = r4_random ( 0.25, 1.0 );
      x = r4_random ( -2.5, 2.5 );
      y = r4_random ( -2.5, 2.5 );
      z = r4_random ( -2.5, 2.5 );
      glTranslatef ( x, y, z ); 
      glRotatef ( angle, 1.0, 0.0, 0.0 );
      glutSolidCone ( r1, r2, 15, 15 );
    }
    else if ( j == 2 )
    {
      r1 = r4_random ( 0.25, 1.0 );
      x = r4_random ( -2.5, 2.5 );
      y = r4_random ( -2.5, 2.5 );
      z = r4_random ( -2.5, 2.5 );
      glTranslatef ( x, y, z ); 
      glutSolidSphere ( r1, 15, 15 );
    }

    glPopMatrix ();

  }

  glPopMatrix ();

  glFlush ();
}
//****************************************************************************80

void init ( void )

//****************************************************************************80
// 
//  Purpose:
// 
//    INIT initializes material properties and the light source.
//
//  Modified:
//
//    04 February 2000
//
{
  glLightfv ( GL_LIGHT0, GL_AMBIENT, light_ambient );
  glLightfv ( GL_LIGHT0, GL_DIFFUSE, light_diffuse );
  glLightfv ( GL_LIGHT0, GL_SPECULAR, light_specular );

  light_position[0] = sin ( theta ) * sin ( phi );
  light_position[1] = cos ( theta ) * sin ( phi );
  light_position[2] = cos ( phi );
  light_position[3] = 0.0;

  glLightfv ( GL_LIGHT0, GL_POSITION, light_position );
   
  glEnable ( GL_LIGHTING );
  glEnable ( GL_LIGHT0 );
  glEnable ( GL_DEPTH_TEST );
}
//****************************************************************************80

void keyboard ( unsigned char key, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    KEYBOARD terminates execution upon seeing an ESCAPE character from the keyboard.
//
//  Modified:
//
//    04 February 2000
//
{
   switch ( key )
   {
      case 27:

        exit ( 0 );
        break;
   }
}
//****************************************************************************80

void mouse ( int btn, int state, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    MOUSE responds to user mouse button input.
//
//  Modified:
//
//    04 February 2000
//
{
  if ( state == GLUT_DOWN )
  {
    if ( btn == GLUT_LEFT_BUTTON )
    {
      phi = phi + 3.14159265 / 12.0;

      light_position[0] = sin ( theta ) * sin ( phi );
      light_position[1] = cos ( theta ) * sin ( phi );

      light_position[2] = cos ( phi );
      light_position[3] = 0.0;

      glDisable ( GL_LIGHT0 );
      glDisable ( GL_LIGHTING );

      glLightfv ( GL_LIGHT0, GL_AMBIENT, light_ambient );
      glLightfv ( GL_LIGHT0, GL_DIFFUSE, light_diffuse );
      glLightfv ( GL_LIGHT0, GL_SPECULAR, light_specular );
      glLightfv ( GL_LIGHT0, GL_POSITION, light_position );

      glEnable ( GL_LIGHT0 );
      glEnable ( GL_LIGHTING );

      cout << "New light angle PHI is " << phi << "\n";

      glutPostRedisplay ( );
    }
    else if ( btn == GLUT_RIGHT_BUTTON )
    {
      theta = theta + 3.14159265 / 12.0;

      light_position[0] = sin ( theta ) * sin ( phi );
      light_position[1] = cos ( theta ) * sin ( phi );
      light_position[2] = cos ( phi );
      light_position[3] = 0.0;

      glDisable ( GL_LIGHT0 );
      glDisable ( GL_LIGHTING );

      glLightfv ( GL_LIGHT0, GL_AMBIENT, light_ambient );
      glLightfv ( GL_LIGHT0, GL_DIFFUSE, light_diffuse );
      glLightfv ( GL_LIGHT0, GL_SPECULAR, light_specular );
      glLightfv ( GL_LIGHT0, GL_POSITION, light_position );

      glEnable ( GL_LIGHT0 );
      glEnable ( GL_LIGHTING );

      cout << "New light angle THETA is " << theta << "\n";

      glutPostRedisplay ( );

    }
    else if ( btn == GLUT_MIDDLE_BUTTON )
    {
      cout << "\n";
      cout << "Bye bye!\n";
      exit ( 0 );
    }
  }
}
//****************************************************************************80

float r4_random ( float rlo, float rhi )

//****************************************************************************80
//
//  Purpose:
//
//   R4_RANDOM returns a random real in a given range.
//
//  Modified:
//
//    04 February 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_RANDOM, the randomly chosen value.
//
//    Input, float RLO, RHI, the minimum and maximum values.
//
{
  float t;

  t = ( float ) ( rand ( ) ) / 32767.0;

  return ( 1.0 - t ) * rlo + t * rhi;
}
//****************************************************************************80

void reshape ( int w, int h )

//*****************************************************************************
//
//  Purpose:
//
//    RESHAPE is passed to glutReshapeFunc to handle window reshaping requests.
//
//  Modified:
//
//    04 February 2000
//
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();

  if ( w <= h ) 
  {
    glOrtho ( -2.5, 2.5, -2.5*(GLfloat)h/(GLfloat)w, 
      2.5*(GLfloat)h/(GLfloat)w, -10.0, 10.0 );
  }
  else 
  {
    glOrtho ( -2.5*(GLfloat)w/(GLfloat)h, 
      2.5*(GLfloat)w/(GLfloat)h, -2.5, 2.5, -10.0, 10.0 );
  }

  glMatrixMode ( GL_MODELVIEW );
  glLoadIdentity ( );
}
