# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
//
//  This is the include statement I need for Mac OS X.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>

using namespace std;

int main ( int argc, char *argv[] );
void box_draw ( int i, int j, bool state_ij );
void display ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_uniform ( int a, int b, int *seed );
int i4_wrap ( int ival, int ilo, int ihi );
void my_init ( );
void my_keyboard ( unsigned char key, int x, int y );
void my_mouse ( int btn, int mouse_state, int x, int y );
float r4_abs ( float x );
int r4_nint ( float x );
void state_randomize ( int moves, int m, int n, bool state[], int *seed );
void state_reset ( int m, int n, bool state[], int i, int j );
void state_update ( int m, int n, bool state[] );
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
  bool wrap;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LIFE_OPENGL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int iwrap;
  int j;

  timestamp ( );

  cout << "\n";
  cout << "LIFE_OPENGL\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  Expect width N.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "LIFE_OPENGL:\n";
    cout << "  Please the WIDTH of the board, in boxes  (10 is a good value).\n";

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
    cout << "LIFE_OPENGL:\n";
    cout << "  Please enter the HEIGHT of the board, in boxes (10 is a good value).\n";

    cin >> m;
  }
  else 
  {
    m = atoi ( argv[2] );
  }
//
//  Wrap around?
//
  if ( argc <= 3 ) 
  {
    cout << "\n";
    cout << "LIFE_OPENGL:\n";
    cout << "  Please enter 1 if the board is wraparound, 0 otherwise\n";

    cin >> iwrap;
    wrap = bool ( iwrap );
  }
  else 
  {
    iwrap = atoi ( argv[3] );
    wrap = ( bool ) iwrap;
  }

  state = new bool[m * n];
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      state[i+j*m] = false;
    }
  }

  glutInit ( &argc, argv );
//
//  Use double buffering; otherwise the screen jitters when the user
//  updates it.
//
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB );

  if ( m == n )
  {
    box_size = ( 500 / m );
    pixel_width = 500;
    pixel_height = 500;
  }
  else if ( m < n )
  {
    box_size = ( 500 / n );
    pixel_width = n * box_size;
    pixel_height = m * box_size;
  }
  else if ( n < m )
  {
    box_size = ( 500 / m );
    pixel_width = n * box_size;
    pixel_height = m * box_size;
  }
  cout << "  Box size = " << box_size << "\n";
  cout << "  Pixels(WxH):  " << pixel_width  << "  "
                             << pixel_height << "\n";

  glutInitWindowSize ( pixel_width, pixel_height );
  glutInitWindowPosition ( 0, 0 );

  glutCreateWindow ( "Life" );
  glutDisplayFunc ( display );
  my_init ( );
  glutKeyboardFunc ( my_keyboard );
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
  cout << "LIFE\n";
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
//    BOX_DRAW draws one box of the array.
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
  float p[2];
  GLfloat a;
  GLfloat b;
  GLfloat c;

  GLfloat gray[3] = { 0.8, 0.8, 0.8 };
  GLfloat yellow[3] = { 1.0, 1.0, 0.0 };
 
  if ( state_ij )
  {
    glColor3fv ( yellow );
  }
  else
  {
    glColor3fv ( gray );
  }

  c = box_size;
  a =                      j   * c;
  b = ( m - 1 - i ) * c;
//
//  Fill boxes with gray or yellow.
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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
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

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

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
  gluOrtho2D ( 0.0, ( double ) pixel_width, 0.0, (double ) pixel_height );
  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

void my_keyboard ( unsigned char key, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    MY_KEYBOARD reacts to keyboard events.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  state_update ( m, n, state );
//
//  Redisplay the screen.
//  Since this causes a jerky screen, it would be best to double buffer!
//
  display ( );

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

  i = y / box_size;
  j = x / box_size;

  if ( btn == GLUT_LEFT_BUTTON && mouse_state == GLUT_DOWN )
  {
    state_reset ( m, n, state, i, j );
  }
  else if ( btn == GLUT_MIDDLE_BUTTON && mouse_state == GLUT_DOWN )
  {
    state_reset ( m, n, state, i, j );
  }
  else if ( btn == GLUT_RIGHT_BUTTON && mouse_state == GLUT_DOWN )
  {
    state_reset ( m, n, state, i, j );
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

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      state[i+j*m] = false;
    }
  }

  for ( k = 0; k < moves; k++ )
  {
    i = i4_uniform ( 0, m - 1, seed );
    j = i4_uniform ( 0, n - 1, seed );
    state[i+j*m] = true;
  }

  return;
}
//****************************************************************************80

void state_reset ( int m, int n, bool state[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    STATE_RESET allows the user to reset the state of a cell.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, bool STATE[M*N], the current state array.
//
//    Input, int I, J, the row and column that were pressed.
//
{
  if ( i < 0 || m <= i )
  {
    cerr << "\n";
    cerr << "STATE_RESET - Fatal error!\n";
    cerr << "  Illegal row index I.\n";
    exit ( 1 );
  }

  if ( j < 0 || n <= j )
  {
    cerr << "\n";
    cerr << "STATE_RESET - Fatal error!\n";
    cerr << "  Illegal column index J.\n";
    exit ( 1 );
  }
  state[i+j*m] = !state[i+j*m];

  return;
}
//****************************************************************************80

void state_update ( int m, int n, bool state[] )

//****************************************************************************80
//
//  Purpose:
//
//    STATE_UPDATE updates the cell states.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, bool STATE[M*N], the current state.
//
{
  bool *state_new;
  int i;
  int im;
  int ip;
  int j;
  int jm;
  int jp;
  int neighbors;

  state_new = new bool[m*n];

  if ( wrap )
  {

    for ( j = 0; j < n; j++ )
    {
      jp = i4_wrap ( j + 1, 0, n - 1 );
      jm = i4_wrap ( j - 1, 0, n - 1 );
      for ( i = 0; i < m; i++ )
      {
        ip = i4_wrap ( i + 1, 0, m - 1 );
        im = i4_wrap ( i - 1, 0, m - 1 );
        neighbors = ( int ) state[im+jm*m]
                  + ( int ) state[im+j *m]
                  + ( int ) state[im+jp*m]
                  + ( int ) state[i +jm*m]
                  + ( int ) state[i +jp*m]
                  + ( int ) state[ip+jm*m]
                  + ( int ) state[ip+j *m]
                  + ( int ) state[ip+jp*m];

        if ( neighbors == 3 )
        {
          state_new[i+j*m] = true;
        }
        else if ( neighbors == 2 )
        {
          state_new[i+j*m] = state[i+j*m];
        }
        else
        {
          state_new[i+j*m] = false;
        }
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        neighbors = 0;
        if ( 0 < i && 0 < j )
        {
          neighbors = neighbors + ( int ) state[(i-1)+(j-1)*m];
        }
        if ( 0 < i )
        {
          neighbors = neighbors + ( int ) state[(i-1)+(j)*m];
        }
        if ( 0 < i && j < n - 1 )
        {
          neighbors = neighbors + ( int ) state[(i-1)+(j+1)*m];
        }
        if ( 0 < j )
        {
          neighbors = neighbors + ( int ) state[(i)+(j-1)*m];
        }
        if ( j < n - 1 )
        {
          neighbors = neighbors + ( int ) state[(i)+(j+1)*m];
        }
        if ( i < m - 1 && 0 < j )
        {
          neighbors = neighbors + ( int ) state[(i+1)+(j-1)*m];
        }
        if ( i < m - 1 )
        {
          neighbors = neighbors + ( int ) state[(i+1)+(j)*m];
        }
        if ( i < m - 1 && j < n - 1 )
        {
          neighbors = neighbors + ( int ) state[(i+1)+(j+1)*m];
        }

        if ( neighbors == 3 )
        {
          state_new[i+j*m] = true;
        }
        else if ( neighbors == 2 )
        {
          state_new[i+j*m] = state[i+j*m];
        }
        else
        {
          state_new[i+j*m] = false;
        }
      }
    }
  }
//
//  Update.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0;  i < m; i++ )
    {
      state[i+j*m] = state_new[i+j*m];
    }
  }

  delete [] state_new;

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
