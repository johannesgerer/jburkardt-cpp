# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "wandzura.hpp"

//****************************************************************************80

void file_name_inc ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the input string contains no digits, a blank string is returned.
//
//    If a blank string is input, then an error condition results.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, character *FILE_NAME, (a pointer to) the character string
//    to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = s_len_trim ( file_name );

  if ( lens <= 0 )
  {
    cout << "\n";
    cout << "FILE_NAME_INC - Fatal error!\n";
    cout << "  Input file name is blank.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens-1; 0 <= i; i-- )
  {
    c = *(file_name+i);

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      if ( c == '9' )
      {
        c = '0';
        *(file_name+i) = c;
      }
      else
      {
        c = c + 1;
        *(file_name+i) = c;
        return;
      }
    }
  }

  if ( change == 0 )
  {
    strcpy ( file_name, " " );
  }

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
//    I4_MIN returns the smaller of two I4's.
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
//  Example:
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
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80*

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80*
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

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal double precision number,
//    and is usually defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  return HUGE_VAL;
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         Value
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
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
//****************************************************************************80

void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are literally halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0) and
//    (0,1) respectively.
//
//    Input, int N, the number of objects to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] )
                 + t[i+1*2] *       + ref[0+j*2]
                 + t[i+2*2] *                    + ref[1+j*2];
    }
  }

  return;
}
//****************************************************************************80

int s_len_trim ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
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

double triangle_area ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA computes the area of a triangle.
//
//  Discussion:
//
//    If the triangle's vertices are given in counter clockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed
//    area of a triangle, it is possible to easily compute the area
//    of a nonconvex polygon as the sum of the (possibly negative)
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA, the area of the triangle.
//
{
  double area;

  area = 0.5 * (
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) +
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) +
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );

  return area;
}
//****************************************************************************80

void triangle_points_plot ( char *file_name, double node_xy[], int node_show,
  int point_num, double point_xy[], int point_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_POINTS_PLOT plots a triangle and some points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_NAME, the name of the output file.
//
//    Input, double NODE_XY[2*3], the coordinates of the nodes
//    of the triangle.
//
//    Input, int NODE_SHOW,
//   -1, do not show the triangle, or the nodes.
//    0, show the triangle, do not show the nodes;
//    1, show the triangle and the nodes;
//    2, show the triangle, the nodes and number them.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[2*POINT_NUM], the coordinates of the
//    points.
//
//    Input, int POINT_SHOW,
//    0, do not show the points;
//    1, show the points;
//    2, show the points and number them.
//
{
  int circle_size;
  int delta;
  int e;
  ofstream file_unit;
  int i;
  int node;
  int node_num = 3;
  int point;
  char string[40];
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( x_max < node_xy[0+node*2] )
    {
      x_max = node_xy[0+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( x_max < point_xy[0+point*2] )
    {
      x_max = point_xy[0+point*2];
    }
  }

  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_xy[0+node*2] < x_min )
    {
      x_min = node_xy[0+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( point_xy[0+point*2] < x_min )
    {
      x_min = point_xy[0+point*2];
    }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( y_max < node_xy[1+node*2] )
    {
      y_max = node_xy[1+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( y_max < point_xy[1+point*2] )
    {
      y_max = point_xy[1+point*2];
    }
  }

  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_xy[1+node*2] < y_min )
    {
      y_min = node_xy[1+node*2];
    }
  }
  for ( point = 0; point < point_num; point++ )
  {
    if ( point_xy[1+point*2] < y_min )
    {
      y_min = point_xy[1+point*2];
    }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max      = y_ps_max - delta;
    y_ps_min      = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "TRIANGLE_POINTS_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%//PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: triangulation_order3_plot.C\n";
  file_unit << "%%Title: " << file_name << "\n";
  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page: 1 1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "%  Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << "0.900  0.900  0.900 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "%  Set the RGB color to black.\n";
  file_unit << "%\n";
  file_unit << "0.000  0.000  0.000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size.\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title.\n";
  file_unit << "%\n";
  file_unit << "%  210  702  moveto\n";
  file_unit << "%  (Triangulation)  show\n";
  file_unit << "%\n";
  file_unit << "%  Define a clipping polygon.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( 1 <= node_show )
  {
    circle_size = 5;

    file_unit << "%\n";
    file_unit << "%  Draw filled dots at the nodes.\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to blue.\n";
    file_unit << "%\n";
    file_unit << "0.000  0.150  0.750 setrgbcolor\n";
    file_unit << "%\n";

    for ( node = 0; node < 3; node++ )
    {
      x_ps = r8_nint (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (         node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  " << x_ps
                << "  " << y_ps
                << "  " << circle_size
                << "  0 360 arc closepath fill\n";
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to darker blue.\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850 setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = r8_nint (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit  << "  " << x_ps
                 << "  " << y_ps + 5
                 << "  moveto (" << node+1 << ") show\n";
    }
  }
//
//  Draw the points.
//
  if ( point_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( point_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( point_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( point_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= point_show )
  {
    file_unit << "%\n";
    file_unit << "%  Draw filled dots at the points.\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to green.\n";
    file_unit << "%\n";
    file_unit << "0.150  0.750  0.000 setrgbcolor\n";
    file_unit << "%\n";

    for ( point = 0; point < point_num; point++ )
    {
      x_ps = r8_nint (
        ( ( x_max - point_xy[0+point*2]         ) * ( double ) ( x_ps_min )
        + (         point_xy[0+point*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                       - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - point_xy[1+point*2]         ) * ( double ) ( y_ps_min )
        + (         point_xy[1+point*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                       - y_min ) );

      file_unit << "newpath  " << x_ps
                << "  " << y_ps
                << "  " << circle_size
                << "  0 360 arc closepath fill\n";
    }
  }
//
//  Label the points.
//
  if ( 2 <= point_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the point:\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to darker green.\n";
    file_unit << "%\n";
    file_unit << "0.250  0.850  0.000 setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";

    for ( point = 0; point < point_num; point++ )
    {
      x_ps = r8_nint (
        ( ( x_max - point_xy[0+point*2]         ) * ( double ) ( x_ps_min )
        + (       + point_xy[0+point*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                       - x_min ) );

      y_ps = r8_nint (
        ( ( y_max - point_xy[1+point*2]         ) * ( double ) ( y_ps_min )
        + (         point_xy[1+point*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                       - y_min ) );

      file_unit << "  " << x_ps
                << "  " << y_ps + 5
                << "  moveto (" << point+1 << ") show\n";
    }
  }
//
//  Draw the triangle.
//
  if ( 0 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to red.\n";
    file_unit << "%\n";
    file_unit << "0.900  0.200  0.100 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "%  Draw the triangle.\n";
    file_unit << "%\n";

    file_unit << "newpath\n";

    for ( i = 0; i <= 3; i++ )
    {
      node = i4_wrap ( i, 0, 2 );

      x_ps = ( r8_nint ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (         node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( r8_nint ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      if ( i == 0 )
      {
        file_unit << x_ps << "  " << y_ps << "  moveto\n";
      }
      else
      {
        file_unit << x_ps << "  " << y_ps << "  lineto\n";
      }
    }
    file_unit << "stroke\n";
  }

  file_unit << "%\n";
  file_unit << "restore  showpage\n";
  file_unit << "%\n";
  file_unit << "%  End of page.\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
//****************************************************************************80

int wandzura_degree ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_DEGREE returns the degree of a given Wandzura rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int WANDZURA_DEGREE, the polynomial degree of exactness of
//    the rule.
//
{
  int degree;

  if ( rule == 1 )
  {
    degree = 5;
  }
  else if ( rule == 2 )
  {
    degree = 10;
  }
  else if ( rule == 3 )
  {
    degree = 15;
  }
  else if ( rule == 4 )
  {
    degree = 20;
  }
  else if ( rule == 5 )
  {
    degree = 25;
  }
  else if ( rule == 6 )
  {
    degree = 30;
  }
  else
  {
    degree = -1;
    cout << "\n";
    cout << "WANDZURA_DEGREE - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return degree;
}
//****************************************************************************80

int wandzura_order_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_ORDER_NUM returns the order of a given Wandzura rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int WANDZURA_ORDER_NUM, the order (number of points) of the rule.
//
{
  int order;
  int order_num;
  int *suborder;
  int suborder_num;

  suborder_num = wandzura_suborder_num ( rule );

  suborder = wandzura_suborder ( rule, suborder_num );

  order_num = 0;
  for ( order = 0; order < suborder_num; order++ )
  {
    order_num = order_num + suborder[order];
  }

  delete [] suborder;

  return order_num;
}
//****************************************************************************80

void wandzura_rule ( int rule, int order_num, double xy[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_RULE returns the points and weights of a Wandzura rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int ORDER_NUM, the order (number of points) of the rule.
//
//    Output, double XY[2*ORDER_NUM], the points of the rule.
//
//    Output, double W[ORDER_NUM], the weights of the rule.
//
{
  int k;
  int o;
  int s;
  int *suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
//
//  Get the suborder information.
//
  suborder_num = wandzura_suborder_num ( rule );

  suborder_xyz = new double[3*suborder_num];
  suborder_w = new double[suborder_num];

  suborder = wandzura_suborder ( rule, suborder_num );

  wandzura_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
//
//  Expand the suborder information to a full order rule.
//
  o = 0;

  for ( s = 0; s < suborder_num; s++ )
  {
    if ( suborder[s] == 1 )
    {
      xy[0+o*2] = suborder_xyz[0+s*3];
      xy[1+o*2] = suborder_xyz[1+s*3];
      w[o] = suborder_w[s];
      o = o + 1;
    }
    else if ( suborder[s] == 3 )
    {
      for ( k = 0; k < 3; k++ )
      {
        xy[0+o*2] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
        xy[1+o*2] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
        w[o] = suborder_w[s];
        o = o + 1;
      }
    }
    else if ( suborder[s] == 6 )
    {
      for ( k = 0; k < 3; k++ )
      {
        xy[0+o*2] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
        xy[1+o*2] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
        w[o] = suborder_w[s];
        o = o + 1;
      }

      for ( k = 0; k < 3; k++ )
      {
        xy[0+o*2] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
        xy[1+o*2] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
        w[o] = suborder_w[s];
        o = o + 1;
      }
    }
    else
    {
      cout << "\n";
      cout << "WANDZURA_RULE - Fatal error!\n;";
      cout << "  Illegal SUBORDER(" << s << ") = " << suborder[s] << "\n";
      exit ( 1 );
    }
  }

  delete [] suborder;
  delete [] suborder_xyz;
  delete [] suborder_w;

  return;
}
//****************************************************************************80

int wandzura_rule_num ( void )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_RULE_NUM returns the number of Wandzura rules available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Output, int WANDZURA_RULE_NUM, the number of rules available.
//
{
  int rule_num;

  rule_num = 6;

  return rule_num;
}
//****************************************************************************80

int *wandzura_suborder ( int rule, int suborder_num )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBORDER returns the suborders for a Wandzura rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int WANDZURA_SUBORDER[SUBORDER_NUM], the suborders of the rule.
//
{
  int *suborder;

  suborder = new int[suborder_num];

  if ( rule == 1 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
  }
  else if ( rule == 2 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 6;
    suborder[6] = 6;
  }
  else if ( rule == 3 )
  {
    suborder[0] = 3;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 6;
    suborder[7] = 6;
    suborder[8] = 6;
    suborder[9] = 6;
    suborder[10] = 6;
    suborder[11] = 6;
  }
  else if ( rule == 4 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 3;
    suborder[7] = 3;
    suborder[8] = 3;
    suborder[9] = 6;
    suborder[10] = 6;
    suborder[11] = 6;
    suborder[12] = 6;
    suborder[13] = 6;
    suborder[14] = 6;
    suborder[15] = 6;
    suborder[16] = 6;
    suborder[17] = 6;
    suborder[18] = 6;
  }
  else if ( rule == 5 )
  {
    suborder[0] = 3;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 3;
    suborder[7] = 3;
    suborder[8] = 3;
    suborder[9] = 3;
    suborder[10] = 6;
    suborder[11] = 6;
    suborder[12] = 6;
    suborder[13] = 6;
    suborder[14] = 6;
    suborder[15] = 6;
    suborder[16] = 6;
    suborder[17] = 6;
    suborder[18] = 6;
    suborder[19] = 6;
    suborder[20] = 6;
    suborder[21] = 6;
    suborder[22] = 6;
    suborder[23] = 6;
    suborder[24] = 6;
    suborder[25] = 6;
  }
  else if ( rule == 6 )
  {
    suborder[0] = 1;
    suborder[1] = 3;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 3;
    suborder[5] = 3;
    suborder[6] = 3;
    suborder[7] = 3;
    suborder[8] = 3;
    suborder[9] = 3;
    suborder[10] = 3;
    suborder[11] = 3;
    suborder[12] = 3;
    suborder[13] = 6;
    suborder[14] = 6;
    suborder[15] = 6;
    suborder[16] = 6;
    suborder[17] = 6;
    suborder[18] = 6;
    suborder[19] = 6;
    suborder[20] = 6;
    suborder[21] = 6;
    suborder[22] = 6;
    suborder[23] = 6;
    suborder[24] = 6;
    suborder[25] = 6;
    suborder[26] = 6;
    suborder[27] = 6;
    suborder[28] = 6;
    suborder[29] = 6;
    suborder[30] = 6;
    suborder[31] = 6;
    suborder[32] = 6;
    suborder[33] = 6;
    suborder[34] = 6;
    suborder[35] = 6;
  }
  else
  {
    cout << "\n";
    cout << "WANDZURA_SUBORDER - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder;
}
//****************************************************************************80

int wandzura_suborder_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBORDER_NUM returns the number of suborders for a Wandzura rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int WANDZURA_SUBORDER_NUM, the number of suborders of the rule.
//
{
  int suborder_num;

  if ( rule == 1 )
  {
    suborder_num = 3;
  }
  else if ( rule == 2 )
  {
    suborder_num = 7;
  }
  else if ( rule == 3 )
  {
    suborder_num = 12;
  }
  else if ( rule == 4 )
  {
    suborder_num = 19;
  }
  else if ( rule == 5 )
  {
    suborder_num = 26;
  }
  else if ( rule == 6 )
  {
    suborder_num = 36;
  }
  else
  {
    suborder_num = -1;
    cout << "\n";
    cout << "WANDZURA_SUBORDER_NUM - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder_num;
}
//****************************************************************************80

void wandzura_subrule ( int rule, int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE returns a compressed Wandzura rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  if ( rule == 1 )
  {
    wandzura_subrule_1 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 2 )
  {
    wandzura_subrule_2 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 3 )
  {
    wandzura_subrule_3 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 4 )
  {
    wandzura_subrule_4 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 5 )
  {
    wandzura_subrule_5 ( suborder_num, suborder_xyz, suborder_w );
  }
  else if ( rule == 6 )
  {
    wandzura_subrule_6 ( suborder_num, suborder_xyz, suborder_w );
  }
  else
  {
    cout << "\n";
    cout << "WANDZURA_SUBRULE - Fatal error!\n";
    cout << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void wandzura_subrule_1 ( int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE_1 returns a compressed Wandzura rule 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_1[3*3] = {
      0.33333333333333, 0.33333333333333, 0.33333333333333,
      0.05971587178977, 0.47014206410512, 0.47014206410512,
      0.79742698535309, 0.10128650732346, 0.10128650732346 };

  double suborder_w_rule_1[3] = {
      0.2250000000000000E+00,
      0.1323941527885062E+00,
      0.1259391805448271E+00 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_1[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_1[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_1[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_1[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule_2 ( int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE_2 returns a compressed Wandzura rule 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_2[3*7] = {
      0.33333333333333, 0.33333333333333, 0.33333333333333,
      0.00426913409105, 0.49786543295447, 0.49786543295447,
      0.14397510054189, 0.42801244972906, 0.42801244972906,
      0.63048717451355, 0.18475641274322, 0.18475641274322,
      0.95903756285664, 0.02048121857168, 0.02048121857168,
      0.03500298989727, 0.13657357625603, 0.82842343384669,
      0.03754907025844, 0.33274360058864, 0.62970732915292 };

  double suborder_w_rule_2[7] = {
      0.8352339980519638E-01,
      0.7229850592056743E-02,
      0.7449217792098051E-01,
      0.7864647340310853E-01,
      0.6928323087107504E-02,
      0.2951832033477940E-01,
      0.3957936719606124E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_2[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_2[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_2[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_2[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule_3 ( int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE_3 returns a compressed Wandzura rule 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_3[3*12] = {
      0.08343840726175, 0.45828079636912, 0.45828079636913,
      0.19277907084174, 0.40361046457913, 0.40361046457913,
      0.41360566417395, 0.29319716791303, 0.29319716791303,
      0.70706442611445, 0.14646778694277, 0.14646778694277,
      0.88727426466879, 0.05636286766560, 0.05636286766560,
      0.96684974628326, 0.01657512685837, 0.01657512685837,
      0.00991220330923, 0.23953455415479, 0.75055324253598,
      0.01580377063023, 0.40487880731834, 0.57931742205143,
      0.00514360881697, 0.09500211311304, 0.89985427806998,
      0.04892232575299, 0.14975310732227, 0.80132456692474,
      0.06876874863252, 0.28691961244133, 0.64431163892615,
      0.16840441812470, 0.28183566809908, 0.54975991377622 };

  double suborder_w_rule_3[12] = {
      0.3266181884880529E-01,
      0.2741281803136436E-01,
      0.2651003659870330E-01,
      0.2921596213648611E-01,
      0.1058460806624399E-01,
      0.3614643064092035E-02,
      0.8527748101709436E-02,
      0.1391617651669193E-01,
      0.4291932940734835E-02,
      0.1623532928177489E-01,
      0.2560734092126239E-01,
      0.3308819553164567E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_3[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_3[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_3[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_3[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule_4 ( int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE_4 returns a compressed Wandzura rule 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_4[3*19] = {
      0.33333333333333, 0.33333333333333, 0.33333333333333,
      0.00150064932443, 0.49924967533779, 0.49924967533779,
      0.09413975193895, 0.45293012403052, 0.45293012403052,
      0.20447212408953, 0.39776393795524, 0.39776393795524,
      0.47099959493443, 0.26450020253279, 0.26450020253279,
      0.57796207181585, 0.21101896409208, 0.21101896409208,
      0.78452878565746, 0.10773560717127, 0.10773560717127,
      0.92186182432439, 0.03906908783780, 0.03906908783780,
      0.97765124054134, 0.01117437972933, 0.01117437972933,
      0.00534961818734, 0.06354966590835, 0.93110071590431,
      0.00795481706620, 0.15710691894071, 0.83493826399309,
      0.01042239828126, 0.39564211436437, 0.59393548735436,
      0.01096441479612, 0.27316757071291, 0.71586801449097,
      0.03856671208546, 0.10178538248502, 0.85964790542952,
      0.03558050781722, 0.44665854917641, 0.51776094300637,
      0.04967081636276, 0.19901079414950, 0.75131838948773,
      0.05851972508433, 0.32426118369228, 0.61721909122339,
      0.12149778700439, 0.20853136321013, 0.66997084978547,
      0.14071084494394, 0.32317056653626, 0.53611858851980 };

  double suborder_w_rule_4[19] = {
      0.2761042699769952E-01,
      0.1779029547326740E-02,
      0.2011239811396117E-01,
      0.2681784725933157E-01,
      0.2452313380150201E-01,
      0.1639457841069539E-01,
      0.1479590739864960E-01,
      0.4579282277704251E-02,
      0.1651826515576217E-02,
      0.2349170908575584E-02,
      0.4465925754181793E-02,
      0.6099566807907972E-02,
      0.6891081327188203E-02,
      0.7997475072478163E-02,
      0.7386134285336024E-02,
      0.1279933187864826E-01,
      0.1725807117569655E-01,
      0.1867294590293547E-01,
      0.2281822405839526E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_4[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_4[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_4[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_4[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule_5 ( int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE_5 returns a compressed Wandzura rule 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_5[3*26] = {
      0.02794648307317, 0.48602675846341, 0.48602675846341,
      0.13117860132765, 0.43441069933617, 0.43441069933617,
      0.22022172951207, 0.38988913524396, 0.38988913524396,
      0.40311353196039, 0.29844323401980, 0.29844323401980,
      0.53191165532526, 0.23404417233737, 0.23404417233737,
      0.69706333078196, 0.15146833460902, 0.15146833460902,
      0.77453221290801, 0.11273389354599, 0.11273389354599,
      0.84456861581695, 0.07771569209153, 0.07771569209153,
      0.93021381277141, 0.03489309361430, 0.03489309361430,
      0.98548363075813, 0.00725818462093, 0.00725818462093,
      0.00129235270444, 0.22721445215336, 0.77149319514219,
      0.00539970127212, 0.43501055485357, 0.55958974387431,
      0.00638400303398, 0.32030959927220, 0.67330639769382,
      0.00502821150199, 0.09175032228001, 0.90322146621800,
      0.00682675862178, 0.03801083585872, 0.95516240551949,
      0.01001619963993, 0.15742521848531, 0.83255858187476,
      0.02575781317339, 0.23988965977853, 0.73435252704808,
      0.03022789811992, 0.36194311812606, 0.60782898375402,
      0.03050499010716, 0.08355196095483, 0.88594304893801,
      0.04595654736257, 0.14844322073242, 0.80560023190501,
      0.06744280054028, 0.28373970872753, 0.64881749073219,
      0.07004509141591, 0.40689937511879, 0.52305553346530,
      0.08391152464012, 0.19411398702489, 0.72197448833499,
      0.12037553567715, 0.32413434700070, 0.55549011732214,
      0.14806689915737, 0.22927748355598, 0.62265561728665,
      0.19177186586733, 0.32561812259598, 0.48261001153669 };

  double suborder_w_rule_5[26] = {
      0.8005581880020417E-02,
      0.1594707683239050E-01,
      0.1310914123079553E-01,
      0.1958300096563562E-01,
      0.1647088544153727E-01,
      0.8547279074092100E-02,
      0.8161885857226492E-02,
      0.6121146539983779E-02,
      0.2908498264936665E-02,
      0.6922752456619963E-03,
      0.1248289199277397E-02,
      0.3404752908803022E-02,
      0.3359654326064051E-02,
      0.1716156539496754E-02,
      0.1480856316715606E-02,
      0.3511312610728685E-02,
      0.7393550149706484E-02,
      0.7983087477376558E-02,
      0.4355962613158041E-02,
      0.7365056701417832E-02,
      0.1096357284641955E-01,
      0.1174996174354112E-01,
      0.1001560071379857E-01,
      0.1330964078762868E-01,
      0.1415444650522614E-01,
      0.1488137956116801E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_5[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_5[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_5[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_5[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule_6 ( int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE_6 returns a compressed Wandzura rule 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_6[3*36] = {
      0.33333333333333, 0.33333333333333, 0.33333333333333,
      0.00733011643277, 0.49633494178362, 0.49633494178362,
      0.08299567580296, 0.45850216209852, 0.45850216209852,
      0.15098095612541, 0.42450952193729, 0.42450952193729,
      0.23590585989217, 0.38204707005392, 0.38204707005392,
      0.43802430840785, 0.28098784579608, 0.28098784579608,
      0.54530204829193, 0.22734897585403, 0.22734897585403,
      0.65088177698254, 0.17455911150873, 0.17455911150873,
      0.75348314559713, 0.12325842720144, 0.12325842720144,
      0.83983154221561, 0.08008422889220, 0.08008422889220,
      0.90445106518420, 0.04777446740790, 0.04777446740790,
      0.95655897063972, 0.02172051468014, 0.02172051468014,
      0.99047064476913, 0.00476467761544, 0.00476467761544,
      0.00092537119335, 0.41529527091331, 0.58377935789334,
      0.00138592585556, 0.06118990978535, 0.93742416435909,
      0.00368241545591, 0.16490869013691, 0.83140889440718,
      0.00390322342416, 0.02503506223200, 0.97106171434384,
      0.00323324815501, 0.30606446515110, 0.69070228669389,
      0.00646743211224, 0.10707328373022, 0.88645928415754,
      0.00324747549133, 0.22995754934558, 0.76679497516308,
      0.00867509080675, 0.33703663330578, 0.65428827588746,
      0.01559702646731, 0.05625657618206, 0.92814639735063,
      0.01797672125369, 0.40245137521240, 0.57957190353391,
      0.01712424535389, 0.24365470201083, 0.73922105263528,
      0.02288340534658, 0.16538958561453, 0.81172700903888,
      0.03273759728777, 0.09930187449585, 0.86796052821639,
      0.03382101234234, 0.30847833306905, 0.65770065458860,
      0.03554761446002, 0.46066831859211, 0.50378406694787,
      0.05053979030687, 0.21881529945393, 0.73064491023920,
      0.05701471491573, 0.37920955156027, 0.56377573352399,
      0.06415280642120, 0.14296081941819, 0.79288637416061,
      0.08050114828763, 0.28373128210592, 0.63576756960645,
      0.10436706813453, 0.19673744100444, 0.69889549086103,
      0.11384489442875, 0.35588914121166, 0.53026596435959,
      0.14536348771552, 0.25981868535191, 0.59481782693256,
      0.18994565282198, 0.32192318123130, 0.48813116594672 };

  double suborder_w_rule_6[36] = {
      0.1557996020289920E-01,
      0.3177233700534134E-02,
      0.1048342663573077E-01,
      0.1320945957774363E-01,
      0.1497500696627150E-01,
      0.1498790444338419E-01,
      0.1333886474102166E-01,
      0.1088917111390201E-01,
      0.8189440660893461E-02,
      0.5575387588607785E-02,
      0.3191216473411976E-02,
      0.1296715144327045E-02,
      0.2982628261349172E-03,
      0.9989056850788964E-03,
      0.4628508491732533E-03,
      0.1234451336382413E-02,
      0.5707198522432062E-03,
      0.1126946125877624E-02,
      0.1747866949407337E-02,
      0.1182818815031657E-02,
      0.1990839294675034E-02,
      0.1900412795035980E-02,
      0.4498365808817451E-02,
      0.3478719460274719E-02,
      0.4102399036723953E-02,
      0.4021761549744162E-02,
      0.6033164660795066E-02,
      0.3946290302129598E-02,
      0.6644044537680268E-02,
      0.8254305856078458E-02,
      0.6496056633406411E-02,
      0.9252778144146602E-02,
      0.9164920726294280E-02,
      0.1156952462809767E-01,
      0.1176111646760917E-01,
      0.1382470218216540E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xyz[0+s*3] = suborder_xy_rule_6[0+s*3];
    suborder_xyz[1+s*3] = suborder_xy_rule_6[1+s*3];
    suborder_xyz[2+s*3] = suborder_xy_rule_6[2+s*3];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_6[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule2_1 ( int suborder_num, double suborder_xy[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE2_1 returns a compressed Wandzura rule 1.
//
//  Discussion:
//
//    This version of the rules uses as reference the equilateral
//    triangle whose vertices are (-1/2,-sqrt(3)/2), (1,0), (-1/2,sqrt(3)/2).
//
//    This, in fact, is the data as printed in the reference.
//
//    Currently, we don't use this routine at all.  The values of
//    X and Y here could be converted to lie XSI and ETA in the
//    standard (0,0), (1,0), (0,1) reference triangle by
//
//      XSI = ( 2/3) * X                 + 1/3
//      ETA = (-1/3) * X + sqrt(3)/3 * Y + 1/3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XY[2*SUBORDER_NUM],
//    the coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_1[2*3] = {
       0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4104261923153453E+00,   0.0000000000000000E+00,
       0.6961404780296310E+00,   0.0000000000000000E+00  };

  double suborder_w_rule_1[3] = {
      0.2250000000000000E+00,
      0.1323941527885062E+00,
      0.1259391805448271E+00 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xy[0+s*2] = suborder_xy_rule_1[0+s*2];
    suborder_xy[1+s*2] = suborder_xy_rule_1[1+s*2];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_1[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule2_2 ( int suborder_num, double suborder_xy[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE2_2 returns a compressed Wandzura rule 2.
//
//  Discussion:
//
//    This version of the rules uses as reference the equilateral
//    triangle whose vertices are (-1/2,-sqrt(3)/2), (1,0), (-1/2,sqrt(3)/2).
//
//    This, in fact, is the data as printed in the reference.
//
//    Currently, we don't use this routine at all.  The values of
//    X and Y here could be converted to lie XSI and ETA in the
//    standard (0,0), (1,0), (0,1) reference triangle by
//
//      XSI = ( 2/3) * X                 + 1/3
//      ETA = (-1/3) * X + sqrt(3)/3 * Y + 1/3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XY[2*SUBORDER_NUM],
//    the coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_2[2*7] = {
       0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4935962988634245E+00,   0.0000000000000000E+00,
      -0.2840373491871686E+00,   0.0000000000000000E+00,
       0.4457307617703263E+00,   0.0000000000000000E+00,
       0.9385563442849673E+00,   0.0000000000000000E+00,
      -0.4474955151540920E+00,  -0.5991595522781586E+00,
      -0.4436763946123360E+00,  -0.2571781329392130E+00  };

  double suborder_w_rule_2[7] = {
      0.8352339980519638E-01,
      0.7229850592056743E-02,
      0.7449217792098051E-01,
      0.7864647340310853E-01,
      0.6928323087107504E-02,
      0.2951832033477940E-01,
      0.3957936719606124E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xy[0+s*2] = suborder_xy_rule_2[0+s*2];
    suborder_xy[1+s*2] = suborder_xy_rule_2[1+s*2];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_2[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule2_3 ( int suborder_num, double suborder_xy[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE2_3 returns a compressed Wandzura rule 3.
//
//  Discussion:
//
//    This version of the rules uses as reference the equilateral
//    triangle whose vertices are (-1/2,-sqrt(3)/2), (1,0), (-1/2,sqrt(3)/2).
//
//    This, in fact, is the data as printed in the reference.
//
//    Currently, we don't use this routine at all.  The values of
//    X and Y here could be converted to lie XSI and ETA in the
//    standard (0,0), (1,0), (0,1) reference triangle by
//
//      XSI = ( 2/3) * X                 + 1/3
//      ETA = (-1/3) * X + sqrt(3)/3 * Y + 1/3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XY[2*SUBORDER_NUM],
//    the coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_3[2*12] = {
      -0.3748423891073751E+00,   0.0000000000000000E+00,
      -0.2108313937373917E+00,   0.0000000000000000E+00,
       0.1204084962609239E+00,   0.0000000000000000E+00,
       0.5605966391716812E+00,   0.0000000000000000E+00,
       0.8309113970031897E+00,   0.0000000000000000E+00,
       0.9502746194248890E+00,   0.0000000000000000E+00,
      -0.4851316950361628E+00,  -0.4425551659467111E+00,
      -0.4762943440546580E+00,  -0.1510682717598242E+00,
      -0.4922845867745440E+00,  -0.6970224211436132E+00,
      -0.4266165113705168E+00,  -0.5642774363966393E+00,
      -0.3968468770512212E+00,  -0.3095105740458471E+00,
      -0.2473933728129512E+00,  -0.2320292030461791E+00  };

  double suborder_w_rule_3[12] = {
      0.3266181884880529E-01,
      0.2741281803136436E-01,
      0.2651003659870330E-01,
      0.2921596213648611E-01,
      0.1058460806624399E-01,
      0.3614643064092035E-02,
      0.8527748101709436E-02,
      0.1391617651669193E-01,
      0.4291932940734835E-02,
      0.1623532928177489E-01,
      0.2560734092126239E-01,
      0.3308819553164567E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xy[0+s*2] = suborder_xy_rule_3[0+s*2];
    suborder_xy[1+s*2] = suborder_xy_rule_3[1+s*2];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_3[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule2_4 ( int suborder_num, double suborder_xy[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE2_4 returns a compressed Wandzura rule 4.
//
//  Discussion:
//
//    This version of the rules uses as reference the equilateral
//    triangle whose vertices are (-1/2,-sqrt(3)/2), (1,0), (-1/2,sqrt(3)/2).
//
//    This, in fact, is the data as printed in the reference.
//
//    Currently, we don't use this routine at all.  The values of
//    X and Y here could be converted to lie XSI and ETA in the
//    standard (0,0), (1,0), (0,1) reference triangle by
//
//      XSI = ( 2/3) * X                 + 1/3
//      ETA = (-1/3) * X + sqrt(3)/3 * Y + 1/3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XY[2*SUBORDER_NUM],
//    the coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_4[2*19] = {
       0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4977490260133565E+00,   0.0000000000000000E+00,
      -0.3587903720915737E+00,   0.0000000000000000E+00,
      -0.1932918138657104E+00,   0.0000000000000000E+00,
       0.2064993924016380E+00,   0.0000000000000000E+00,
       0.3669431077237697E+00,   0.0000000000000000E+00,
       0.6767931784861860E+00,   0.0000000000000000E+00,
       0.8827927364865920E+00,   0.0000000000000000E+00,
       0.9664768608120111E+00,   0.0000000000000000E+00,
      -0.4919755727189941E+00,  -0.7513212483763635E+00,
      -0.4880677744007016E+00,  -0.5870191642967427E+00,
      -0.4843664025781043E+00,  -0.1717270984114328E+00,
      -0.4835533778058150E+00,  -0.3833898305784408E+00,
      -0.4421499318718065E+00,  -0.6563281974461070E+00,
      -0.4466292382741727E+00,  -0.6157647932662624E-01,
      -0.4254937754558538E+00,  -0.4783124082660027E+00,
      -0.4122204123735024E+00,  -0.2537089901614676E+00,
      -0.3177533194934086E+00,  -0.3996183176834929E+00,
      -0.2889337325840919E+00,  -0.1844183967233982E+00  };

  double suborder_w_rule_4[19] = {
      0.2761042699769952E-01,
      0.1779029547326740E-02,
      0.2011239811396117E-01,
      0.2681784725933157E-01,
      0.2452313380150201E-01,
      0.1639457841069539E-01,
      0.1479590739864960E-01,
      0.4579282277704251E-02,
      0.1651826515576217E-02,
      0.2349170908575584E-02,
      0.4465925754181793E-02,
      0.6099566807907972E-02,
      0.6891081327188203E-02,
      0.7997475072478163E-02,
      0.7386134285336024E-02,
      0.1279933187864826E-01,
      0.1725807117569655E-01,
      0.1867294590293547E-01,
      0.2281822405839526E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xy[0+s*2] = suborder_xy_rule_4[0+s*2];
    suborder_xy[1+s*2] = suborder_xy_rule_4[1+s*2];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_4[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule2_5 ( int suborder_num, double suborder_xy[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE2_5 returns a compressed Wandzura rule 5.
//
//  Discussion:
//
//    This version of the rules uses as reference the equilateral
//    triangle whose vertices are (-1/2,-sqrt(3)/2), (1,0), (-1/2,sqrt(3)/2).
//
//    This, in fact, is the data as printed in the reference.
//
//    Currently, we don't use this routine at all.  The values of
//    X and Y here could be converted to lie XSI and ETA in the
//    standard (0,0), (1,0), (0,1) reference triangle by
//
//      XSI = ( 2/3) * X                 + 1/3
//      ETA = (-1/3) * X + sqrt(3)/3 * Y + 1/3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XY[2*SUBORDER_NUM],
//    the coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_5[2*26] = {
      -0.4580802753902387E+00,   0.0000000000000000E+00,
      -0.3032320980085228E+00,   0.0000000000000000E+00,
      -0.1696674057318916E+00,   0.0000000000000000E+00,
       0.1046702979405866E+00,   0.0000000000000000E+00,
       0.2978674829878846E+00,   0.0000000000000000E+00,
       0.5455949961729473E+00,   0.0000000000000000E+00,
       0.6617983193620190E+00,   0.0000000000000000E+00,
       0.7668529237254211E+00,   0.0000000000000000E+00,
       0.8953207191571090E+00,   0.0000000000000000E+00,
       0.9782254461372029E+00,   0.0000000000000000E+00,
      -0.4980614709433367E+00,  -0.4713592181681879E+00,
      -0.4919004480918257E+00,  -0.1078887424748246E+00,
      -0.4904239954490375E+00,  -0.3057041948876942E+00,
      -0.4924576827470104E+00,  -0.7027546250883238E+00,
      -0.4897598620673272E+00,  -0.7942765584469995E+00,
      -0.4849757005401057E+00,  -0.5846826436376921E+00,
      -0.4613632802399150E+00,  -0.4282174042835178E+00,
      -0.4546581528201263E+00,  -0.2129434060653430E+00,
      -0.4542425148392569E+00,  -0.6948910659636692E+00,
      -0.4310651789561460E+00,  -0.5691146659505208E+00,
      -0.3988357991895837E+00,  -0.3161666335733065E+00,
      -0.3949323628761341E+00,  -0.1005941839340892E+00,
      -0.3741327130398251E+00,  -0.4571406037889341E+00,
      -0.3194366964842710E+00,  -0.2003599744104858E+00,
      -0.2778996512639500E+00,  -0.3406754571040736E+00,
      -0.2123422011990124E+00,  -0.1359589640107579E+00  };

  double suborder_w_rule_5[26] = {
      0.8005581880020417E-02,
      0.1594707683239050E-01,
      0.1310914123079553E-01,
      0.1958300096563562E-01,
      0.1647088544153727E-01,
      0.8547279074092100E-02,
      0.8161885857226492E-02,
      0.6121146539983779E-02,
      0.2908498264936665E-02,
      0.6922752456619963E-03,
      0.1248289199277397E-02,
      0.3404752908803022E-02,
      0.3359654326064051E-02,
      0.1716156539496754E-02,
      0.1480856316715606E-02,
      0.3511312610728685E-02,
      0.7393550149706484E-02,
      0.7983087477376558E-02,
      0.4355962613158041E-02,
      0.7365056701417832E-02,
      0.1096357284641955E-01,
      0.1174996174354112E-01,
      0.1001560071379857E-01,
      0.1330964078762868E-01,
      0.1415444650522614E-01,
      0.1488137956116801E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xy[0+s*2] = suborder_xy_rule_5[0+s*2];
    suborder_xy[1+s*2] = suborder_xy_rule_5[1+s*2];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_5[s];
  }

  return;
}
//****************************************************************************80

void wandzura_subrule2_6 ( int suborder_num, double suborder_xy[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    WANDZURA_SUBRULE2_6 returns a compressed Wandzura rule 6.
//
//  Discussion:
//
//    This version of the rules uses as reference the equilateral
//    triangle whose vertices are (-1/2,-sqrt(3)/2), (1,0), (-1/2,sqrt(3)/2).
//
//    This, in fact, is the data as printed in the reference.
//
//    Currently, we don't use this routine at all.  The values of
//    X and Y here could be converted to lie XSI and ETA in the
//    standard (0,0), (1,0), (0,1) reference triangle by
//
//      XSI = ( 2/3) * X                 + 1/3
//      ETA = (-1/3) * X + sqrt(3)/3 * Y + 1/3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wandzura, Hong Xiao,
//    Symmetric Quadrature Rules on a Triangle,
//    Computers and Mathematics with Applications,
//    Volume 45, Number 12, June 2003, pages 1829-1840.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XY[2*SUBORDER_NUM],
//    the coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xy_rule_6[2*36] = {
       0.0000000000000000E+00,   0.0000000000000000E+00,
      -0.4890048253508517E+00,   0.0000000000000000E+00,
      -0.3755064862955532E+00,   0.0000000000000000E+00,
      -0.2735285658118844E+00,   0.0000000000000000E+00,
      -0.1461412101617502E+00,   0.0000000000000000E+00,
       0.1570364626117722E+00,   0.0000000000000000E+00,
       0.3179530724378968E+00,   0.0000000000000000E+00,
       0.4763226654738105E+00,   0.0000000000000000E+00,
       0.6302247183956902E+00,   0.0000000000000000E+00,
       0.7597473133234094E+00,   0.0000000000000000E+00,
       0.8566765977763036E+00,   0.0000000000000000E+00,
       0.9348384559595755E+00,   0.0000000000000000E+00,
       0.9857059671536891E+00,   0.0000000000000000E+00,
      -0.4986119432099803E+00,  -0.1459114994581331E+00,
      -0.4979211112166541E+00,  -0.7588411241269780E+00,
      -0.4944763768161339E+00,  -0.5772061085255766E+00,
      -0.4941451648637610E+00,  -0.8192831133859931E+00,
      -0.4951501277674842E+00,  -0.3331061247123685E+00,
      -0.4902988518316453E+00,  -0.6749680757240147E+00,
      -0.4951287867630010E+00,  -0.4649148484601980E+00,
      -0.4869873637898693E+00,  -0.2747479818680760E+00,
      -0.4766044602990292E+00,  -0.7550787344330482E+00,
      -0.4730349181194722E+00,  -0.1533908770581512E+00,
      -0.4743136319691660E+00,  -0.4291730489015232E+00,
      -0.4656748919801272E+00,  -0.5597446281020688E+00,
      -0.4508936040683500E+00,  -0.6656779209607333E+00,
      -0.4492684814864886E+00,  -0.3024354020045064E+00,
      -0.4466785783099771E+00,  -0.3733933337926417E-01,
      -0.4241903145397002E+00,  -0.4432574453491491E+00,
      -0.4144779276264017E+00,  -0.1598390022600824E+00,
      -0.4037707903681949E+00,  -0.5628520409756346E+00,
      -0.3792482775685616E+00,  -0.3048723680294163E+00,
      -0.3434493977982042E+00,  -0.4348816278906578E+00,
      -0.3292326583568731E+00,  -0.1510147586773290E+00,
      -0.2819547684267144E+00,  -0.2901177668548256E+00,
      -0.2150815207670319E+00,  -0.1439403370753732E+00    };

  double suborder_w_rule_6[36] = {
      0.1557996020289920E-01,
      0.3177233700534134E-02,
      0.1048342663573077E-01,
      0.1320945957774363E-01,
      0.1497500696627150E-01,
      0.1498790444338419E-01,
      0.1333886474102166E-01,
      0.1088917111390201E-01,
      0.8189440660893461E-02,
      0.5575387588607785E-02,
      0.3191216473411976E-02,
      0.1296715144327045E-02,
      0.2982628261349172E-03,
      0.9989056850788964E-03,
      0.4628508491732533E-03,
      0.1234451336382413E-02,
      0.5707198522432062E-03,
      0.1126946125877624E-02,
      0.1747866949407337E-02,
      0.1182818815031657E-02,
      0.1990839294675034E-02,
      0.1900412795035980E-02,
      0.4498365808817451E-02,
      0.3478719460274719E-02,
      0.4102399036723953E-02,
      0.4021761549744162E-02,
      0.6033164660795066E-02,
      0.3946290302129598E-02,
      0.6644044537680268E-02,
      0.8254305856078458E-02,
      0.6496056633406411E-02,
      0.9252778144146602E-02,
      0.9164920726294280E-02,
      0.1156952462809767E-01,
      0.1176111646760917E-01,
      0.1382470218216540E-01 };

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_xy[0+s*2] = suborder_xy_rule_6[0+s*2];
    suborder_xy[1+s*2] = suborder_xy_rule_6[1+s*2];
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = suborder_w_rule_6[s];
  }

  return;
}
