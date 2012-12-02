# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
//
//  This is the include statement I need for Mac OS X.
//
# include <GLUT/glut.h>
//
//  The following include statement might be more typical:
//
//# include <GL/glut.h>

using namespace std;

int main ( int argc, char *argv[] );
void box_draw ( int i, int j, bool state_ij );
void display ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform ( int a, int b, int *seed );
void my_init ( ) ;
void my_mouse ( int btn, int mouse_state, int x, int y );
float r4_abs ( float x );
int r4_nint ( float x );
void state_randomize ( int moves, int m, int n, bool state[], int *seed );
void state_update ( int m, int n, bool state[], int i, int j );
void timestamp ( );
//
//  Global data.
//
  int box_size;
  int seed = 123456789;
  bool *state;
  int m;
  int n;
  int pixel_height;
  int pixel_width;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LIGHTS_OUT_OPENGL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int moves;

  timestamp ( );

  cout << "\n";
  cout << "LIGHTS_OUT_OPENGL\n";
  cout << "  C++ version\n";
  cout << "  This program sets up a version of the \"LIGHTS OUT\" game.\n";
  cout << "\n";
  cout << "  Your goal is to turn all the lights out.\n";
  cout << "  Clicking in any square switches it from ON to OFF or vice versa.\n";
  cout << "  But it also switches the four neighbors.\n";
  cout << "\n"; 
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  Expect width N.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "LIGHTS_OUT_OPENGL:\n";
    cout << "  Enter the WIDTH of the board, (10 is a good value).\n";

    cin >> n;
  }
  else 
  {
    n = atoi ( argv[1] );
  }
//
//  Expect height M.
//
  if ( argc <= 2 ) 
  {
    cout << "\n";
    cout << "LIGHTS_OUT_OPENGL:\n";
    cout << "  Enter the HEIGHT of the board, (10 is a good value).\n";

    cin >> m;
  }
  else 
  {
    m = atoi ( argv[2] );
  }
//
//  Expect number of randomized moves MOVES.
//
  if ( argc <= 3 ) 
  {
    cout << "\n";
    cout << "LIGHTS_OUT_OPENGL:\n";
    cout << "  Enter the number of scrambling steps.\n";
    cout << "  (Use a low number like 3 for beginners.)\n";

    cin >> moves;
  }
  else 
  {
    moves = atoi ( argv[3] );
  }
//
//  Randomize the state.
//
  state = new bool[m * n];

  state_randomize ( moves, m, n, state, &seed );

  cout << "\n";
  cout << "  The board will be " << n << " boxes wide by " << m << " boxes high.\n";

  glutInit ( &argc, argv );
//
//  Use double buffering; otherwise the screen jitters when the user
//  updates it.
//
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB );
//
//  Set up the screen, using 800 pixels in the long dimension.
//
  if ( m == n )
  {
    box_size = ( 800 / m );
    pixel_width = 800;
    pixel_height = 800;
  }
  else if ( m < n )
  {
    box_size = ( 800 / n );
    pixel_width = n * box_size;
    pixel_height = m * box_size;
  }
  else if ( n < m )
  {
    box_size = ( 800 / m );
    pixel_width = n * box_size;
    pixel_height = m * box_size;
  }
  cout << "  Box size = " << box_size << "\n";
  cout << "  Pixels(WxH):  " << pixel_width  << "  "
                             << pixel_height << "\n";
//
//  Now call OpenGL to get going.
//
  glutInitWindowSize ( pixel_width, pixel_height );
  glutInitWindowPosition ( 0, 0 );

  glutCreateWindow ( "Lights Out!" );
  glutDisplayFunc ( display );
  my_init ( );
  glutMouseFunc ( my_mouse );
  glutMainLoop ( );
//
//  Free memory.
//
  delete [] state;
//
//  Terminate.
//
  cout << "\n";
  cout << "LIGHTS_OUT_OPENGL\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void box_draw ( int i, int j, bool state_ij ) 

//****************************************************************************80
//
//  Purpose:
//
//    BOX_DRAW draws one box of the Lights Out array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the row and column of the box.
//
//    Input, bool STATE_IJ, is TRUE if the box is "ON".
//
{
  float p[2];
  GLfloat a;
  GLfloat b;
  GLfloat c;
  GLfloat gray[3] = { 0.8, 0.8, 0.8 };
  GLfloat yellow[3] = { 1.0, 1.0, 0.0 };
//
//  Set the color.
//
  if ( state_ij )
  {
    glColor3fv ( yellow );
  }
  else
  {
    glColor3fv ( gray );
  }
//
//  Locate (A,B), the lower left corner of the box.
//
//    A,B+C---A+C,B+C
//     |         |
//     |         |
//    A,B-----A+C,B
//
  c = box_size;
  a =                      j   * c;
  b = ( m - 1 - i ) * c;
//
//  Fill the box with color, but leave a margin of 3 pixels.
//
  glBegin ( GL_POLYGON );
    p[0] = a + 3;
    p[1] = b + 3;
    glVertex2fv ( p );
    p[0] = a + c - 3;
    p[1] = b + 3;
    glVertex2fv ( p );
    p[0] = a + c - 3;
    p[1] = b + c - 3;
    glVertex2fv ( p );
    p[0] = a + 3;
    p[1] = b + c - 3;
    glVertex2fv ( p );
  glEnd ( );
//
//  Draw box boundaries in BLUE.
//
  glColor3f ( 0.0, 0.0, 1.0 );

  glBegin ( GL_LINE_LOOP );
    p[0] = a;
    p[1] = b;
    glVertex2fv ( p );
    p[0] = a + c;
    p[1] = b;
    glVertex2fv ( p );
    p[0] = a + c;
    p[1] = b + c;
    glVertex2fv ( p );
    p[0] = a;
    p[1] = b + c;
    glVertex2fv ( p );
  glEnd ( );
//
//  Clear all the buffers.
//
  glFlush ( );

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
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT );
//
//  Draw each box.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      box_draw ( i, j, state[i+j*m] );
    }
  }
  glFlush ( );
//
//  Time to swap buffers.
//
  glutSwapBuffers ( );

  return;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 ) 
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

void my_init ( ) 

//****************************************************************************80
//
//  Purpose:
//
//    MY_INIT initializes OpenGL state variables dealing with viewing and attributes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  glClearColor ( 1.0, 1.0, 1.0, 0.0 );

  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );
//
//  Change this to proportions for MxN
//
  gluOrtho2D ( 0.0, ( double ) pixel_width, 0.0, ( double ) pixel_height );
  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

void my_mouse ( int btn, int mouse_state, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    MY_MOUSE reacts to mouse events.
//
//  Discussion:
//
//    Right now, the only mouse event is a click on a box inside the board,
//    which causes that box and its neighbors to change state.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  int i;
  int j;
  int k;
//
//  (I,J) are the coordinates of the point where we have clicked the mouse.
//
  i = y / box_size;
  j = x / box_size;
//
//  Update the board, so that box (I,J) and its neighbors are switched.
//
  if ( btn == GLUT_LEFT_BUTTON && mouse_state == GLUT_DOWN )
  {
    state_update ( m, n, state, i, j );
  }
  else if ( btn == GLUT_MIDDLE_BUTTON && mouse_state == GLUT_DOWN )
  {
    state_update ( m, n, state, i, j );
  }
  else if ( btn == GLUT_RIGHT_BUTTON && mouse_state == GLUT_DOWN )
  {
    state_update ( m, n, state, i, j );
  }
//
//  Redisplay the screen.
//  Since this causes a jerky screen, it would be best to double buffer!
//
  display ( );

  return;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

void state_randomize ( int moves, int m, int n, bool state[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    STATE_RANDOMIZE randomizes the state.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MOVES, the number of moves to make.
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, bool STATE[M*N], the Lights Out state.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
{
  int i;
  int j;
  int k;
//
//  Start with all boxes OFF.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      state[i+j*m] = false;
    }
  }
//
//  Choose a box at random, and update it (and its four neighbors).
//
  for ( k = 0; k < moves; k++ )
  {
    i = i4_uniform ( 0, m - 1, seed );
    j = i4_uniform ( 0, n - 1, seed );
    state_update ( m, n, state, i, j );
  }

  return;
}
//****************************************************************************80

void state_update ( int m, int n, bool state[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    STATE_UPDATE updates the state after button (I,J) has been pressed.
//
//  Discussion:
//
//    Reverse the states of 
//
//               (I-1,J)
//      (I,J-1)  (I,  J)  (I,J+1)
//               (I+1,J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, bool STATE[M*N], the Lights Out state.
//
//    Input, int I, J, the row and column that were pressed.
//
{
  int center;
  int down;
  int left;
  int right;
  int up;

  if ( i < 0 || m <= i )
  {
    cerr << "\n";
    cerr << "STATE_UPDATE - Fatal error!\n";
    cerr << "  Illegal row index I.\n";
    exit ( 1 );
  }

  if ( j < 0 || n <= j )
  {
    cerr << "\n";
    cerr << "STATE_UPDATE - Fatal error!\n";
    cerr << "  Illegal column index J.\n";
    exit ( 1 );
  }
//
//  Locate the center box and its four neighbors.
//
  up     = ( i - 1 ) +   j       * m;
  down   = ( i + 1 ) +   j       * m;
  center =   i       +   j       * m;
  left   =   i       + ( j - 1 ) * m;
  right  =   i       + ( j + 1 ) * m;
//
//  Reverse the center, and reverse each neighbor, as long as
//  we don't have to go outside the legal range.
//
  if ( 0 < i )
  {
    state[up] = !state[up];
  }

  if ( 0 < j )
  {
    state[left] = !state[left];
  }

  state[center] = !state[center];

  if ( j < n - 1 )
  {
    state[right] = !state[right];
  }

  if ( i < m - 1 )
  {
    state[down] = !state[down];
  }

  return;
}
//****************************************************************************80

void timestamp ( )

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
//    04 October 2003
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
