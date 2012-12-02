# include <cstdlib>
# include <ctime>
# include <cmath>
# include <fstream>
# include <iostream>
# include <iomanip>
# include <cstring>

//
//  This is the include statement I need for Mac OS X.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

using namespace std;

int main ( int argc, char *argv[] );
void display ( );
void mouse ( int btn, int state, int x, int y );
void myinit ( );
void myReshape ( int w, int h );
float *r4mat_uniform_01 ( int m, int n, int *seed );
float *r4mat_zero_new ( int m, int n );
void r4vec_normal_01 ( int n, int *seed, float x[] );
float *r4vec_uniform_01 ( int n, int *seed );
void spin_image ( );
void timestamp ( );
float *uniform_on_sphere01_map ( int dim_num, int n, int *seed );
//
//  Global data.
//
  static GLint axis = 2;
  int click_num = 0;
  int dim_num = 3;
  int *face_data = NULL;
  int face_data_num;
  int face_num;
  int *face_pointer = NULL;
  static GLfloat *gen_color;
  int gen_num;
  static GLfloat *gen_vec;
  int pixel_height;
  int pixel_width;
  static GLfloat *point_color;
  static int point_num = 0;
  static GLfloat *point_vec;
  int seed = 123456789;
  bool spinning = true;
  static GLfloat theta[3] = { 0.0, 0.0, 0.0 };
  double theta_speed = 0.020;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_VORONOI_DISPLAY_OPENGL.
//
//  Discussion:
//
//    This program chooses random generators on a sphere, assigns each a
//    color, and then gradually colors in other points on the sphere with
//    the color of the nearest generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2009
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
  timestamp ( );

  cout << "\n";
  cout << "SPHERE_VORONOI_DISPLAY_OPENGL:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  This is a program which uses OpenGL\n";
  cout << "  to display an approximate Voronoi diagram on a sphere.\n";
//
//  If the number of generators was not specified, enter it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "SPHERE_VORONOI_DISPLAY_OPENGL:\n";
    cout << "  Please enter the number of generators:\n";

    cin >> gen_num;
  }
  else 
  {
   gen_num = atoi ( argv[1] );
  }
//
//  Something here about the number of grid points.
//
  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize ( 500, 500 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( "Voronoi on a Sphere" );
  glutReshapeFunc ( myReshape );
  glutDisplayFunc ( display );
  glutIdleFunc ( spin_image );
  glutMouseFunc ( mouse );
//
//  Enable hidden surface removal.
//
  glEnable ( GL_DEPTH_TEST );

  myinit ( );

  glutMainLoop ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_VORONOI_DISPLAY_OPENGL:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

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
//    31 January 2009
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
  int gen;
  int hi;
  int offset;
  int point;
  GLfloat point_size;
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glLoadIdentity ( );
  
  glRotatef ( theta[0], 1.0, 0.0, 0.0 );
  glRotatef ( theta[1], 0.0, 1.0, 0.0 );
  glRotatef ( theta[2], 0.0, 0.0, 1.0 );
//
//  Draw SOME of the points.
//
  point_size = 5.0;

  glPointSize ( point_size );

  if ( 7 <= click_num )
  {

    hi = ( point_num * ( ( click_num - 6 ) / 2 ) ) / 20;

    if ( point_num < hi ) 
    {
      hi = point_num;
    }

    for ( point = 0; point < hi; point++ )
    {
      offset = point * 3;

      glColor3fv ( point_color+offset );

      glBegin ( GL_POINTS );

        glVertex3fv ( point_vec+offset );

      glEnd ( );
    }

  }
//
//  Draw the generators.
//
  point_size = 15.0;

  glPointSize ( point_size );

  for ( gen = 0; gen < gen_num; gen++ )
  { 
    offset = gen * 3;

    glColor3fv ( gen_color+offset );

    glBegin ( GL_POINTS );

      glVertex3fv ( gen_vec+offset );

    glEnd ( );
  }
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
//    dedicated one axis to each.  
//
//    Since Apple prefers the esthetics of a one button mouse, we're forced
//    to live with that choice.  This routine alternately pauses rotation,
//    or increments the rotation axis by 1, no matter which button is pushed.
//
//  Modified:
//
//    30 December 2008
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
  float dist;
  float dist_min;
  int gen;
  int gen_min;
  int point;

  if ( ( btn == GLUT_LEFT_BUTTON   && state == GLUT_DOWN ) ||
       ( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN ) ||
       ( btn == GLUT_RIGHT_BUTTON  && state == GLUT_DOWN ) )
  {
    if ( spinning )
    {
      spinning = false;
      theta_speed = 0.0;
    }
    else
    {
      spinning = true;
      theta_speed = 0.020;
      axis = axis + 1;
    }

    click_num = click_num + 1;
  }

  axis = axis % 3;

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
//    30 December 2008
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
{
  double dist;
  double dist_min;
  int gen;
  int gen_min;
  GLfloat line_width;
  int point;
  GLfloat point_size;
  GLfloat scale;
  GLfloat t;
//
//  Set the background to WHITE.
//
  glClearColor ( 1.0, 1.0, 1.0, 1.0 );
//
//  Antialiasing.
//
  glEnable ( GL_POINT_SMOOTH );
//
//  Define the generators.
//
  gen_vec = uniform_on_sphere01_map ( 3, gen_num, &seed );
//
//  Define a color for each generator.
//
  gen_color = r4mat_uniform_01 ( 3, gen_num, &seed );
//
//  Brighten up the colors!
//
  for ( gen = 0; gen < gen_num; gen++ )
  {
    t = gen_color[0+gen*3];
    if ( t < gen_color[1+gen*3] )
    {
      t = gen_color[1+gen*3];
    }
    if ( t < gen_color[2+gen*3] )
    {
      t = gen_color[2+gen*3];
    }
    gen_color[0+gen*3] = gen_color[0+gen*3] / t;
    gen_color[1+gen*3] = gen_color[1+gen*3] / t;
    gen_color[2+gen*3] = gen_color[2+gen*3] / t;
  }
//
//  Set up the points NOW.
//
  point_num = 100000;

  point_vec = uniform_on_sphere01_map ( 3, point_num, &seed );

  point_color = new float[3*point_num];

  for ( point = 0; point < point_num; point++ )
  {
    gen_min = -1;
    dist_min = 10000.0;
    for ( gen = 0; gen < gen_num; gen++ )
    {
      dist = pow ( point_vec[0+point*3] - gen_vec[0+gen*3], 2 )
           + pow ( point_vec[1+point*3] - gen_vec[1+gen*3], 2 )
           + pow ( point_vec[2+point*3] - gen_vec[2+gen*3], 2 );
      if ( dist < dist_min )
      {
        dist_min = dist;
        gen_min = gen;
      }
    }
    point_color[0+point*3] = gen_color[0+gen_min*3];
    point_color[1+point*3] = gen_color[1+gen_min*3];
    point_color[2+point*3] = gen_color[2+gen_min*3];
  }

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
//    30 December 2008
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
      -1.05, +1.05, 
      - 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w, 
      -10.0, 10.0 );
  }
  else
  {
    glOrtho ( 
      - 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w,  
      - 1.05, + 1.05, 
      -10.0, 10.0 );
  }

  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

float *r4mat_uniform_01 ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0.  On output, SEED has been updated.
//
//    Output, float R4MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new float[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
      r[i+j*m] = ( float ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

float *r4mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_ZERO_NEW returns a new zeroed R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, float R4MAT_ZERO[M*N], the new zeroed matrix.
//
{
  float *a;
  int i;
  int j;

  a = new float[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
//*****************************************************************************

void r4vec_normal_01 ( int n, int *seed, float x[] )

//*****************************************************************************
//
//  Purpose:
//
//    R4VEC_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, float X[N], a sample of the standard normal PDF.
//
{
  int i;
  int m;
  float pi = 3.141592653589793;
  float *r;
  static int made = 0;
  static int saved = 0;
  int xhi;
  int xlo;
  static float y = 0.0;
//
//  I'd like to allow the user to reset the random number seed.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return;
  }
  else if ( n == 0 )
  {
    return;
  }
//
//  Record the range of X we need to fill in.
//
  xlo = 1;
  xhi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    xlo = 2;
  }
//
//  If we don't need any more values, return.
//
  if ( xhi - xlo + 1 == 0 )
  {
    return;
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  if ( xhi - xlo + 1 == 1 )
  {
    r = r4vec_uniform_01 ( 2, seed );

    x[xhi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =        sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( ( xhi-xlo+1) % 2 ) == 0 )
  {
    m = ( xhi-xlo+1 ) / 2;

    r = r4vec_uniform_01 ( 2*m, seed );

    for ( i = 0; i < 2*m; i = i + 2 )
    {
      x[xlo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[xlo+i]   = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    made = made + xhi - xlo + 1;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    xhi = xhi - 1;

    m = ( xhi-xlo+1 ) / 2 + 1;

    r = r4vec_uniform_01 ( 2*m, seed );

    for ( i = 0; i < 2*m-2; i = i + 2 )
    {
      x[xlo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[xlo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    x[n-1] = sqrt ( -2.0 * log ( r[2*m-2] ) ) * cos ( 2.0 * pi * r[2*m-1] );
    y =      sqrt ( -2.0 * log ( r[2*m-2] ) ) * sin ( 2.0 * pi * r[2*m-1] );

    saved = 1;

    made = made + xhi - xlo + 2;

  }

  delete [] r;

  return;
}
//****************************************************************************80

float *r4vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIFORM_01 returns a unit unit pseudorandom R4VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, float R4VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new float[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( float ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void spin_image ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPIN_IMAGE adjusts the angle of rotation and redisplays the picture.
//
//  Modified:
//
//    15 December 2008
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
  theta[axis] = theta[axis] + theta_speed;

  if ( 360.0 < theta[axis] ) 
  {
    theta[axis] = theta[axis] - 360.0;
  }
  glutPostRedisplay ( );

  return;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

float *uniform_on_sphere01_map ( int dim_num, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russell Cheng,
//    Random Variate Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998, pages 168.
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, float UNIFORM_ON_SPHERE01_MAP[DIM_NUM*N], the points.
//
{
  int i;
  int j;
  float norm;
  float *u;
  float *x;

  u = new float[dim_num];
  x = new float[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    r4vec_normal_01 ( dim_num, seed, u );
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      norm = norm + u[i] * u[i];
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = u[i] / norm;
    }

  }

  delete [] u;
  return x;
}
