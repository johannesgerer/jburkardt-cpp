# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "triangle_ncc_rule.hpp"

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
    cerr << "\n";
    cerr << "FILE_NAME_INC - Fatal error!\n";
    cerr << "  Input file name is blank.\n";
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

int triangle_ncc_degree ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_DEGREE returns the degree of an NCC rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int TRIANGLE_NCC_DEGREE, the polynomial degree of exactness of
//    the rule.
//
{
  int degree;

  if ( 1 <= rule && rule <= 9 )
  {
    degree = rule - 1;
  }
  else
  {
    degree = -1;
    cerr << "\n";
    cerr << "TRIANGLE_NCC_DEGREE - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return degree;
}
//****************************************************************************80

int triangle_ncc_order_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_ORDER_NUM returns the order of an NCC rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int TRIANGLE_NCC_ORDER_NUM, the order (number of points)
//    of the rule.
//
{
  int order;
  int order_num;
  int *suborder;
  int suborder_num;

  suborder_num = triangle_ncc_suborder_num ( rule );

  suborder = triangle_ncc_suborder ( rule, suborder_num );

  order_num = 0;
  for ( order = 0; order < suborder_num; order++ )
  {
    order_num = order_num + suborder[order];
  }

  delete [] suborder;

  return order_num;
}
//****************************************************************************80

void triangle_ncc_rule ( int rule, int order_num, double xy[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_RULE returns the points and weights of an NCC rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
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
  suborder_num = triangle_ncc_suborder_num ( rule );

  suborder_xyz = new double[3*suborder_num];
  suborder_w = new double[suborder_num];

  suborder = triangle_ncc_suborder ( rule, suborder_num );

  triangle_ncc_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
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
      cerr << "\n";
      cerr << "TRIANGLE_NCC_RULE - Fatal error!\n;";
      cerr << "  Illegal SUBORDER(" << s << ") = " << suborder[s] << "\n";
      exit ( 1 );
    }
  }

  delete [] suborder;
  delete [] suborder_xyz;
  delete [] suborder_w;

  return;
}
//****************************************************************************80

int triangle_ncc_rule_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_RULE_NUM returns the number of NCC rules available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Output, int TRIANGLE_NCC_RULE_NUM, the number of rules available.
//
{
  int rule_num;

  rule_num = 9;

  return rule_num;
}
//****************************************************************************80

int *triangle_ncc_suborder ( int rule, int suborder_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBORDER returns the suborders for an NCC rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int TRIANGLE_NCC_SUBORDER[SUBORDER_NUM],
//    the suborders of the rule.
//
{
  int *suborder;

  suborder = new int[suborder_num];

  if ( rule == 1 )
  {
    suborder[0] = 1;
  }
  else if ( rule == 2 )
  {
    suborder[0] = 3;
  }
  else if ( rule == 3 )
  {
    suborder[0] = 3;
  }
  else if ( rule == 4 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 1;
  }
  else if ( rule == 5 )
  {
    suborder[0] = 6;
    suborder[1] = 3;
    suborder[2] = 3;
  }
  else if ( rule == 6 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 6;
    suborder[3] = 3;
    suborder[4] = 3;
  }
  else if ( rule == 7 )
  {
    suborder[0] = 6;
    suborder[1] = 6;
    suborder[2] = 3;
    suborder[3] = 3;
    suborder[4] = 6;
    suborder[5] = 1;
  }
  else if ( rule == 8 )
  {
    suborder[0] = 3;
    suborder[1] = 6;
    suborder[2] = 6;
    suborder[3] = 3;
    suborder[4] = 6;
    suborder[5] = 6;
    suborder[6] = 3;
    suborder[7] = 3;
  }
  else if ( rule == 9 )
  {
    suborder[0] = 6;
    suborder[1] = 6;
    suborder[2] = 3;
    suborder[3] = 6;
    suborder[4] = 6;
    suborder[5] = 3;
    suborder[6] = 6;
    suborder[7] = 3;
    suborder[8] = 3;
  }
  else
  {
    cerr << "\n";
    cerr << "TRIANGLE_NCC_SUBORDER - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder;
}
//****************************************************************************80

int triangle_ncc_suborder_num ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBORDER_NUM returns the number of suborders for an NCC rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int TRIANGLE_NCC_SUBORDER_NUM, the number of suborders
//    of the rule.
//
{
  int suborder_num;

  if ( rule == 1 )
  {
    suborder_num = 1;
  }
  else if ( rule == 2 )
  {
    suborder_num = 1;
  }
  else if ( rule == 3 )
  {
    suborder_num = 1;
  }
  else if ( rule == 4 )
  {
    suborder_num = 3;
  }
  else if ( rule == 5 )
  {
    suborder_num = 3;
  }
  else if ( rule == 6 )
  {
    suborder_num = 5;
  }
  else if ( rule == 7 )
  {
    suborder_num = 6;
  }
  else if ( rule == 8 )
  {
    suborder_num = 8;
  }
  else if ( rule == 9 )
  {
    suborder_num = 9;
  }
  else
  {
    suborder_num = -1;
    cerr << "\n";
    cerr << "TRIANGLE_NCC_SUBORDER_NUM - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  return suborder_num;
}
//****************************************************************************80

void triangle_ncc_subrule ( int rule, int suborder_num, double suborder_xyz[],
  double suborder_w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE returns a compressed NCC rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
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
  int i;
  int s;
  int suborder_w_d;
  int *suborder_w_n;
  int suborder_xyz_d;
  int *suborder_xyz_n;

  suborder_xyz_n = new int[3*suborder_num];
  suborder_w_n = new int[suborder_num];

  if ( rule == 1 )
  {
    triangle_ncc_subrule_01 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 2 )
  {
    triangle_ncc_subrule_02 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 3 )
  {
    triangle_ncc_subrule_03 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 4 )
  {
    triangle_ncc_subrule_04 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 5 )
  {
    triangle_ncc_subrule_05 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 6 )
  {
    triangle_ncc_subrule_06 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 7 )
  {
    triangle_ncc_subrule_07 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 8 )
  {
    triangle_ncc_subrule_08 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 9 )
  {
    triangle_ncc_subrule_09 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else
  {
    cerr << "\n";
    cerr << "TRIANGLE_NCC_SUBRULE - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit ( 1 );
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz[i+s*3] =
          ( double ) ( suborder_xyz_n[i+s*3] )
        / ( double ) ( suborder_xyz_d );
    }
  }
  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = ( double ) suborder_w_n[s] / ( double ) suborder_w_d;
  }

  delete [] suborder_w_n;
  delete [] suborder_xyz_n;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_01 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_01 returns a compressed NCC rule 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_01[3*1] = {
      1, 1, 1
  };
  int suborder_xyz_d_01 = 3;
  int suborder_w_n_01[1] = { 1 };
  int suborder_w_d_01 = 1;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_01[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_01;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_01[s];
  }
  *suborder_w_d = suborder_w_d_01;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_02 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_02 returns a compressed NCC rule 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_02[3*1] = {
      1, 0, 0
  };
  int suborder_xyz_d_02 = 1;
  int suborder_w_n_02[1] = { 1 };
  int suborder_w_d_02 = 3;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_02[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_02;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_02[s];
  }
  *suborder_w_d = suborder_w_d_02;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_03 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_03 returns a compressed NCC rule 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_03[3*1] = {
      1, 1, 0
  };
  int suborder_xyz_d_03 = 2;
  int suborder_w_n_03[1] = { 1 };
  int suborder_w_d_03 = 3;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_03[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_03;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_03[s];
  }
  *suborder_w_d = suborder_w_d_03;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_04 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_04 returns a compressed NCC rule 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_04[3*3] = {
      3, 0, 0,
      2, 1, 0,
      1, 1, 1
  };
  int suborder_xyz_d_04 = 3;
  int suborder_w_n_04[3] = { 4, 9, 54 };
  int suborder_w_d_04 = 120;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_04[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_04;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_04[s];
  }
  *suborder_w_d = suborder_w_d_04;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_05 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_05 returns a compressed NCC rule 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_05[3*3] = {
      3, 1, 0,
      2, 2, 0,
      2, 1, 1
  };
  int suborder_xyz_d_05 = 4;
  int suborder_w_n_05[3] = { 4, -1, 8 };
  int suborder_w_d_05 = 45;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_05[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_05;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_05[s];
  }
  *suborder_w_d = suborder_w_d_05;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_06 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_06 returns a compressed NCC rule 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_06[3*5] = {
      5, 0, 0,
      4, 1, 0,
      3, 2, 0,
      3, 1, 1,
      2, 2, 1
  };
  int suborder_xyz_d_06 = 5;
  int suborder_w_n_06[5] = { 11, 25, 25, 200, 25 };
  int suborder_w_d_06 = 1008;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_06[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_06;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_06[s];
  }
  *suborder_w_d = suborder_w_d_06;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_07 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_07 returns a compressed NCC rule 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_07[3*6] = {
      5, 1, 0,
      4, 2, 0,
      4, 1, 1,
      3, 3, 0,
      3, 2, 1,
      2, 2, 2
  };
  int suborder_xyz_d_07 = 6;
  int suborder_w_n_07[6] = { 36, -27, 72, 64, 72, -54 };
  int suborder_w_d_07 = 840;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_07[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_07;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_07[s];
  }
  *suborder_w_d = suborder_w_d_07;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_08 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_08 returns a compressed NCC rule 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_08[3*8] = {
      7, 0, 0,
      6, 1, 0,
      5, 2, 0,
      5, 1, 1,
      4, 3, 0,
      4, 2, 1,
      3, 3, 1,
      3, 2, 2
  };
  int suborder_xyz_d_08 = 7;
  int suborder_w_n_08[8] = { 1336, 2989, 3577, 32242, 2695, -6860, 44590, 3430 };
  int suborder_w_d_08 = 259200;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_08[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_08;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_08[s];
  }
  *suborder_w_d = suborder_w_d_08;

  return;
}
//****************************************************************************80

void triangle_ncc_subrule_09 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NCC_SUBRULE_09 returns a compressed NCC rule 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Silvester,
//    Symmetric Quadrature Formulae for Simplexes,
//    Mathematics of Computation,
//    Volume 24, Number 109, January 1970, pages 95-100.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
//    the numerators of the barycentric coordinates of the abscissas.
//
//    Output, int *SUBORDER_XYZ_D,
//    the denominator of the barycentric coordinates of the abscissas.
//
//    Output, int SUBORDER_W_N[SUBORDER_NUM],
//    the numerator of the suborder weights.
//
//    Output, int SUBORDER_W_D,
//    the denominator of the suborder weights.
//
{
  int i;
  int s;
  int suborder_xyz_n_09[3*9] = {
      7, 1, 0,
      6, 2, 0,
      6, 1, 1,
      5, 3, 0,
      5, 2, 1,
      4, 4, 0,
      4, 3, 1,
      4, 2, 2,
      3, 3, 2
  };
  int suborder_xyz_d_09 = 8;
  int suborder_w_n_09[9] = {
    368, -468, 704, 1136, 832, -1083, 672, -1448, 1472 };
  int suborder_w_d_09 = 14175;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      suborder_xyz_n[i+s*3] = suborder_xyz_n_09[i+s*3];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_09;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_09[s];
  }
  *suborder_w_d = suborder_w_d_09;

  return;
}

