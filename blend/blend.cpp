# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "blend.hpp"

//****************************************************************************80

void blend_0d1 ( double r, double x0, double x1, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_0D1 extends scalar data at endpoints to a line.
//
//  Diagram:
//
//    0-----r-----1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the coordinate where an interpolated value is desired.
//
//    Input, double X0, X1, the data values at the ends of the line.
//
//    Output, double *X, the interpolated data value at (R).
//
{
  *x = ( 1.0 - r ) * x0 + r * x1;

  return;
}
//****************************************************************************80

void blend_1d1 ( double r, double s, double x00, double x01, double x10,
  double x11, double xr0, double xr1, double x0s, double x1s, double *x )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_1D1 extends scalar data along the boundary into a square.
//
//  Diagram:
//
//    01-----r1-----11
//     |      .      |
//     |      .      |
//    0s.....rs.....1s
//     |      .      |
//     |      .      |
//    00-----r0-----10
//
//  Formula:
//
//    Written as a polynomial in R and S, the interpolation map has the form
//
//      X(R,S) =
//           1     * ( x0s + xr0 - x00 )
//         + r     * ( x00 + x1s - x0s - x10 )
//         + s     * ( x00 + xr0 - x01 - xr1 )
//         + r * s * ( x01 + x10 - x00 - x11 )
//
//    The nonlinear term ( r * s ) has an important role:
//
//    If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in a plane,
//    and the mapping is affine.  All the interpolated data will lie
//    on the plane defined by the four corner values.  In particular,
//    on any line through the square, data values at intermediate points
//    will lie between the values at the endpoints.
//
//    If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
//    not lie in a plane, and the interpolation map is nonlinear.  On
//    any line through the square, data values at intermediate points
//    may lie above or below the data values at the endpoints.  The
//    size of the coefficient of r * s will determine how severe this
//    effect is.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, S, the coordinates where an interpolated value is desired.
//
//    Input, double X00, X01, X10, X11, the data values at the corners.
//
//    Input, double XR0, XR1, X0S, X1S, the data values at points along the edges.
//
//    Output, double *X, the interpolated data value at (R,S).
//
{
  *x =
     - ( 1.0 - r ) * ( 1.0 - s ) * x00
     + ( 1.0 - r )               * x0s
     - ( 1.0 - r ) *         s   * x01
     +               ( 1.0 - s ) * xr0
     +                       s   * xr1
     -         r   * ( 1.0 - s ) * x10
     +         r                 * x1s
     -         r   *         s   * x11;

  return;
}
//****************************************************************************80

void blend_2d1 ( double r, double s, double t,
  double x000, double x001, double x010, double x011,
  double x100, double x101, double x110, double x111,
  double xr00, double xr01, double xr10, double xr11,
  double x0s0, double x0s1, double x1s0, double x1s1,
  double x00t, double x01t, double x10t, double x11t,
  double x0st, double x1st, double xr0t, double xr1t, double xrs0, double xrs1,
  double *x )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_2D1 extends scalar data along the surface into a cube.
//
//  Diagram:
//
//    010-----r10-----110    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    0s0.....rs0.....1s0    0s1.....rs1.....1s1     S
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    001-----r01-----101     +----R
//       BOTTOM                      TOP
//
//    011-----0s1-----001    111-----1s1-----101
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    01t.....0st.....00t    11t.....1st.....10t          T
//      |       .       |      |       .       |          |
//      |       .       |      |       .       |          |
//    010-----0s0-----000    110-----1s0-----100     S----+
//       LEFT                       RIGHT
//
//    001-----r01-----101    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    00t.....r0t.....100    01t.....r1t.....11t     T
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    010-----r10-----110     +----R
//       FRONT                       BACK
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, S, T, the coordinates where an interpolated value is desired.
//
//    Input, double X000, X001, X010, X011, X100, X101, X110, X111, the data
//    values at the corners.
//
//    Input, double XR00, XR01, XR10, XR11, X0S0, X0S1, X1S0, X1S1, X00T, X01T,
//    X10T, X11T, the data values at points along the edges.
//
//    Input, double X0ST, X1ST, XR0T, XR1T, XRS0, XRS1, the data values
//    at points on the faces.
//
//    Output, double *X, the interpolated data value at (R,S,T).
//
{
//
//  Interpolate the interior point.
//

 *x =    ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - t ) * x000
       - ( 1.0 - r ) * ( 1.0 - s )               * x00t
       + ( 1.0 - r ) * ( 1.0 - s ) *         t   * x001
       - ( 1.0 - r )               * ( 1.0 - t ) * x0s0
       + ( 1.0 - r )                             * x0st
       - ( 1.0 - r )               *         t   * x0s1
       + ( 1.0 - r ) *         s   * ( 1.0 - t ) * x010
       - ( 1.0 - r ) *         s                 * x01t
       + ( 1.0 - r ) *         s   *         t   * x011
       -               ( 1.0 - s ) * ( 1.0 - t ) * xr00
       +               ( 1.0 - s )               * xr0t
       -               ( 1.0 - s ) *         t   * xr01
       +                             ( 1.0 - t ) * xrs0
       +                                     t   * xrs1
       -                       s   * ( 1.0 - t ) * xr10
       +                       s                 * xr1t
       -                       s   *         t   * xr11
       +         r   * ( 1.0 - s ) * ( 1.0 - t ) * x100
       -         r   * ( 1.0 - s )               * x10t
       +         r   * ( 1.0 - s ) *         t   * x101
       -         r                 * ( 1.0 - t ) * x1s0
       +         r                               * x1st
       -         r                 *         t   * x1s1
       +         r   *         s   * ( 1.0 - t ) * x110
       -         r   *         s                 * x11t
       +         r   *         s   *         t   * x111;

  return;
}
//****************************************************************************80

void blend_i_0d1 ( double x[], int m )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_I_0D1 extends indexed scalar data at endpoints along a line.
//
//  Diagram:
//
//    ( X0, ..., ..., ..., ..., ..., X6 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input/output, double X[M].
//
//    On input, X[0] and X[M-1] contain scalar values which are to be
//    interpolated through the entries X[1] through X[M-2].  It is assumed
//    that the dependence of the data is linear in the vector index I.
//
//    On output, X[1] through X[M-2] have been assigned interpolated values.
//
//    Input, int M, the number of entries in X.
//
{
  int i;
  double r;

  for ( i = 1; i < m - 1; i++ )
  {
    r = ( double ) i  / ( double ) ( m - 1 );

    blend_0d1 ( r, x[0], x[m-1], &x[i] );

  }

  return;
}
//****************************************************************************80

void blend_ij_0d1 ( double x[], int m1, int m2 )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_IJ_0D1 extends indexed scalar data at corners into a table.
//
//  Diagram:
//
//    ( X00,  ..., ..., ..., ..., ..., X06 )
//    ( ...,  ..., ..., ..., ..., ..., ... )
//    ( ...,  ..., ..., ..., ..., ..., ... )
//    ( ...,  ..., ..., ..., ..., ..., ... )
//    ( X40,  ..., ..., ..., ..., ..., X46 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input/output, double X[M1*M2], a singly dimensioned array that
//    is "really" doubly dimensioned.  The double dimension index [I][J]
//    corresponds to the single dimension index I * M2 + J.
//
//    On input, data values have been stored in the entries
//    [0], [M2-1], [M1*M2-M2] and [M1*M2-1], which correspond to the double
//    dimension entries [0][0], [0][M2-1], [M1-1][0] and [M1-1][M2-1].
//
//    On output, all entries in X have been assigned a value.
//
//    Input, int M1, M2, the number of rows and columns in the doubly
//    dimensioned data.
//
{
  int i;
  int j;
  double r;
  double s;
//
//  Interpolate values along the edges.
//
  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );

    blend_0d1 ( r, x[0], x[m1*m2-m2], &x[i*m2] );

    blend_0d1 ( r, x[m2-1], x[m1*m2-1], &x[i*m2+m2-1] );

  }

  for ( j = 1; j < m2 - 1; j++ )
  {
    s = ( double ) j / ( double ) ( m2 - 1 );

    blend_0d1 ( s, x[0], x[m2-1], &x[j] );

    blend_0d1 ( s, x[m1*m2-m2], x[m1*m2-1], &x[(m1-1)*m2+j] );

  }
//
//  Interpolate values in the interior.
//
  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );

    for ( j = 1; j < m2 - 1; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );

      blend_1d1 ( r, s, x[0], x[m2-1], x[m1*m2-m2], x[m1*m2-1],
        x[i*m2], x[i*m2+m2-1], x[j], x[(m1-1)*m2+j], &x[i*m2+j] );

    }
  }
  return;
}
//****************************************************************************80

void blend_ij_1d1 ( double x[], int m1, int m2 )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_IJ_1D1 extends indexed scalar data along edges into a table.
//
//  Diagram:
//
//    ( X00,  X01,  X02,  X03,  X04,  X05,  X06 )
//    ( X10,  ...,  ...,  ...,  ...,  ...,  X16 )
//    ( X20,  ...,  ...,  ...,  ...,  ...,  X26 )
//    ( X30,  ...,  ...,  ...,  ...,  ...,  X36 )
//    ( X40,  X41,  X42,  X43,  X44,  X45,  X46 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input/output, double X[M1*M2], a singly dimensioned array that
//    is "really" doubly dimensioned.  The double dimension index [I][J]
//    corresponds to the single dimension index I * M2 + J.
//
//    On input, data is contained in the "edge entries"
//    X[0][J], X[I][0], X[M1-1][J] and X[I][M2-1],
//    for I = 0 to M1-1, and J = 0 to M2-1.
//
//    On output, all entries in X have been assigned a value.
//
//    Input, int M1, M2, the number of rows and columns in X.
//
{
  int i;
  int j;
  double r;
  double s;
//
//  Interpolate values in the interior.
//

  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );

    for ( j = 1; j < m2 - 1; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );

      blend_1d1 ( r, s, x[0], x[m2-1], x[m1*m2-m2], x[m1*m2-1],
        x[i*m2], x[i*m2+m2-1], x[j], x[(m1-1)*m2+j], &x[i*m2+j] );

    }
  }
  return;
}
//****************************************************************************80

void blend_ijk_0d1 ( double x[], int m1, int m2, int m3 )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_IJK_0D1 extends indexed scalar data along corners into a cubic table.
//
//  Diagram:
//
//    ( X000,   ...,  ...,  ...,  ...,  ...,  X060 )
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   First "layer"
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//    ( X400,   ...,  ...,  ...,  ...,  ...,  X460 )
//
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   Middle "layers"
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//
//    ( X003,  ...,  ...,  ...,  ...,  ...,  X063  )
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   Last "layer"
//    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
//    ( X403,  ...,  ...,  ...,  ...,  ...,  X463  )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input/output, double X[M1*M2*M3], a singly dimensioned array that
//    is "really" triply dimensioned.  The triple dimension index
//    [I][J][K] corresponds to the single dimension index
//    I * M2*M3 + J * M2 + K
//
//    On input, there is already scalar data in the entries X[I][J][K]
//    corresponding to "cornders" of the table, that is, entries for which
//    each of the three indices I, J and K is equal to their
//    minimum or maximum possible values.
//
//    Input, int M1, M2, M3, the number of rows, columns, and layers in X.
//
{
  int i;
  int j;
  int k;
  double r;
  double s;
  double t;
//
//  Interpolate values along the "edges", that is, index triplets (i,j,k)
//  with exactly two of I, J, K an "extreme" value.
//
  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    blend_0d1 ( r, x[0],              x[(m1-1)*m3*m2],
                  &x[i*m3*m2] );

    blend_0d1 ( r, x[(m2-1)*m2],      x[((m1-1)*m3+m2-1)*m2],
                  &x[(i*m3+m2-1)*m2] );

    blend_0d1 ( r, x[m3-1],           x[(m1-1)*m3*m2+m3-1],
                  &x[i*m3*m2+m3-1] );

    blend_0d1 ( r, x[(m2-1)*m2+m3-1], x[((m1-1)*m3+m2-1)*m2+m3-1],
                  &x[(i*m3+m2-1)*m2+m3-1] );
  }

  for ( j = 1; j < m2 - 1; j++ )
  {
    s = ( double ) j / ( double ) ( m2 - 1 );
    blend_0d1 ( s, x[0], x[(m2-1)*m2],
                  &x[j*m2] );

    blend_0d1 ( s, x[(m1-1)*m3*m2], x[((m1-1)*m3+m2-1)*m2],
                  &x[((m1-1)*m3+j)*m2] );

    blend_0d1 ( s, x[m3-1], x[(m2-1)*m2+m3-1],
                  &x[j*m2+m3-1] );

    blend_0d1 ( s, x[(m1-1)*m3*m2+m3-1], x[((m1-1)*m3+m2-1)*m2+m3-1],
                  &x[((m1-1)*m3+j)*m2+m3-1] );
  }

  for ( k = 1; k < m3 - 1; k++ )
  {
    t = ( double ) k / ( double ) ( m3 - 1 );
    blend_0d1 ( t, x[0], x[m3-1],
                  &x[k] );

    blend_0d1 ( t, x[(m1-1)*m3*m2], x[(m1-1)*m3*m2+m3-1],
                  &x[(m1-1)*m3*m2+k] );

    blend_0d1 ( t, x[(m2-1)*m2], x[(m2-1)*m2+m3-1],
                  &x[(m2-1)*m2+k] );

    blend_0d1 ( t, x[((m1-1)*m3+m2-1)*m2], x[((m1-1)*m3+m2-1)*m2+m3-1],
                  &x[((m1-1)*m3+m2-1)*m2+k] );
  }
//
//  Interpolate values along the "faces", that is, index triplets (i,j,k)
//  with exactly one of I, J, K is an "extreme" value.
//
  for ( j = 1; j < m2 - 1; j++ )
  {
    s = ( double ) j / ( double ) ( m2 - 1 );
    for ( k = 1; k < m3 - 1; k++ )
    {
      t = ( double ) k / ( double ) ( m3 - 1 );

      blend_1d1 ( s, t,
        x[0],                   x[m3-1],
        x[(m2-1)*m2],           x[(m2-1)*m2+m3-1],
        x[j*m2],                x[j*m2+m3-1],
        x[k],                   x[(m2-1)*m2+k],
        &x[j*m2+k] );

      blend_1d1 ( s, t,
        x[(m1-1)*m3*m2],        x[(m1-1)*m3*m2+m3-1],
        x[((m1-1)*m3+m2-1)*m2], x[((m1-1)*m3+m2-1)*m2+m3-1],
        x[((m1-1)*m3+j)*m2],    x[((m1-1)*m3+j)*m2+m3-1],
        x[(m1-1)*m3*m2+k],      x[((m1-1)*m3+m2-1)*m2+k],
        &x[((m1-1)*m3+j)*m2+k] );

    }
  }

  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( k = 1; k < m3 - 1; k++ )
    {
      t = ( double ) k / ( double ) ( m3 - 1 );

      blend_1d1 ( r, t,
        x[0],                   x[m3-1],
        x[(m1-1)*m3*m2],        x[(m1-1)*m3*m2+m3-1],
        x[i*m3*m2],             x[i*m3*m2+m3-1],
        x[k],                   x[(m1-1)*m3*m2+k],
        &x[i*m3*m2+k] );

      blend_1d1 ( r, t,
        x[(m2-1)*m2],           x[(m2-1)*m2+m3-1],
        x[((m1-1)*m3+m2-1)*m2], x[((m1-1)*m3+m2-1)*m2+m3-1],
        x[(i*m3+m2-1)*m2],      x[(i*m3+m2-1)*m2+m3-1],
        x[(m2-1)*m2+k],         x[((m1-1)*m3+m2-1)*m2+k],
        &x[(i*m3+m2-1)*m2+k] );

    }
  }

  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 1; j < m2 - 1; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );

      blend_1d1 ( r, s,
        x[0],            x[(m2-1)*m2],
        x[(m1-1)*m3*m2], x[((m1-1)*m3+m2-1)*m2],
        x[i*m3*m2],      x[(i*m3+m2-1)*m2],
        x[j*m2],         x[((m1-1)*m3+j)*m2],
        &x[(i*m3+j)*m2] );

      blend_1d1 ( r, s,
       x[m3-1],              x[(m2-1)*m2+m3-1],
       x[(m1-1)*m3*m2+m3-1], x[((m1-1)*m3+m2-1)*m2+m3-1],
       x[i*m3*m2+m3-1],      x[(i*m3+m2-1)*m2+m3-1],
       x[j*m2+m3-1],         x[((m1-1)*m3+j)*m2+m3-1],
       &x[(i*m3+j)*m2+m3-1] );

    }
  }
//
//  Interpolate values in the interior.
//
  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );

    for ( j = 1; j < m2 - 1; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );

      for ( k = 1; k < m3 - 1; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );

        blend_2d1 ( r, s, t,
          x[0],                    x[m3-1],
          x[(m2-1)*m2],            x[(m2-1)*m2+m3-1],
          x[(m1-1)*m3*m2],         x[(m1-1)*m3*m2+m3-1],
          x[((m1-1)*m3+m2-1)*m2],  x[((m1-1)*m3+m2-1)*m2+m3-1],
          x[i*m3*m2],              x[i*m3*m2+m3-1],
          x[(i*m3+m2-1)*m2],       x[(i*m3+m2-1)*m2+m3-1],
          x[j*m2],                 x[j*m2+m3-1],
          x[((m1-1)*m3+j)*m2],     x[((m1-1)*m3+j)*m2+m3-1],
          x[k],                    x[(m2-1)*m2+k],
          x[(m1-1)*m3*m2+k],       x[((m1-1)*m3+m2-1)*m2+k],
          x[j*m3+k],               x[((m1-1)*m3+j)*m2+k],
          x[i*m3*m2+k],            x[(i*m3+m2-1)*m2+k],
          x[(i*m3+j)*m2],          x[(i*m3+j)*m2+m3-1],
          &x[(i*m3+j)*m2+k] );

      }

    }

  }
  return;
}
//****************************************************************************80

void blend_ijk_1d1 ( double x[], int m1, int m2, int m3 )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_IJK_1D1 extends indexed scalar data along "edges" into a cubic table.
//
//  Diagram:
//
//    ( X000,   X010,   X020,   X030,   X040,   X050 )
//    ( X100,   ...,    ...,    ...,    ...,    X150 )
//    ( X200,   ...,    ...,    ...,    ...,    X250 )   Layer 1
//    ( X300,   ...,    ...,    ...,    ...,    X350 )
//    ( X400,   X410,   X420,   X430,   X440,   X450 )
//
//    ( X001,   ...,    ...,    ...,    ...,    X051 )
//    ( ....,   ...,    ...,    ...,    ...,    ...  )
//    ( ....,   ...,    ...,    ...,    ...,    ...  )   Layer K
//    ( ....,   ...,    ...,    ...,    ...,    ...  )   1 < K < M3
//    ( X401,   ...,    ...,    ...,    ...,    X451 )
//
//    ( X002,   X012,   X022,   X032,   X042,   X052 )
//    ( X102,   ...,    ...,    ...,    ...,    X152 )
//    ( X202,   ...,    ...,    ...,    ...,    X252 )   Layer M3
//    ( X302    ...,    ...,    ...,    ...,    X352 )
//    ( X402,   X412,   X422,   X432,   X442,   X452 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input/output, double X[M1*M2*M3], a singly dimensioned array that
//    is "really" triply dimensioned.  The triple dimension index
//    [I][J][K] corresponds to the single dimension index
//    I * M2*M3 + J * M2 + K
//
//    On input, there is already scalar data in the entries X[I][J][K]
//    corresponding to "edges" of the table, that is, entries for which
//    at least two of the three indices I, J and K are equal to their
//    minimum or maximum possible values.
//
//    Input, int M1, M2, M3, the number of rows, columns, and layers in X.
//
{
  int i;
  int j;
  int k;
  double r;
  double s;
  double t;
//
//  Interpolate values along the "faces", that is, index triplets (i,j,k)
//  with exactly one of I, J, K is an "extreme" value.
//
  for ( j = 1; j < m2 - 1; j++ )
  {
    s = ( double ) j / ( double ) ( m2 - 1 );
    for ( k = 1; k < m3 - 1; k++ )
    {
      t = ( double ) k / ( double ) ( m3 - 1 );

      blend_1d1 ( s, t,
        x[0],                   x[m3-1],
        x[(m2-1)*m2],           x[(m2-1)*m2+m3-1],
        x[j*m2],                x[j*m2+m3-1],
        x[k],                   x[(m2-1)*m2+k],
        &x[j*m2+k] );

      blend_1d1 ( s, t,
        x[(m1-1)*m3*m2],        x[(m1-1)*m3*m2+m3-1],
        x[((m1-1)*m3+m2-1)*m2], x[((m1-1)*m3+m2-1)*m2+m3-1],
        x[((m1-1)*m3+j)*m2],    x[((m1-1)*m3+j)*m2+m3-1],
        x[(m1-1)*m3*m2+k],      x[((m1-1)*m3+m2-1)*m2+k],
        &x[((m1-1)*m3+j)*m2+k] );

    }
  }

  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( k = 1; k < m3 - 1; k++ )
    {
      t = ( double ) k / ( double ) ( m3 - 1 );

      blend_1d1 ( r, t,
        x[0],                   x[m3-1],
        x[(m1-1)*m3*m2],        x[(m1-1)*m3*m2+m3-1],
        x[i*m3*m2],             x[i*m3*m2+m3-1],
        x[k],                   x[(m1-1)*m3*m2+k],
        &x[i*m3*m2+k] );

      blend_1d1 ( r, t,
        x[(m2-1)*m2],           x[(m2-1)*m2+m3-1],
        x[((m1-1)*m3+m2-1)*m2], x[((m1-1)*m3+m2-1)*m2+m3-1],
        x[(i*m3+m2-1)*m2],      x[(i*m3+m2-1)*m2+m3-1],
        x[(m2-1)*m2+k],         x[((m1-1)*m3+m2-1)*m2+k],
        &x[(i*m3+m2-1)*m2+k] );

    }
  }

  for ( i = 1; i < m1 - 1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 1; j < m2 - 1; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );

      blend_1d1 ( r, s,
        x[0],            x[(m2-1)*m2],
        x[(m1-1)*m3*m2], x[((m1-1)*m3+m2-1)*m2],
        x[i*m3*m2],      x[(i*m3+m2-1)*m2],
        x[j*m2],         x[((m1-1)*m3+j)*m2],
        &x[(i*m3+j)*m2] );

      blend_1d1 ( r, s,
       x[m3-1],              x[(m2-1)*m2+m3-1],
       x[(m1-1)*m3*m2+m3-1], x[((m1-1)*m3+m2-1)*m2+m3-1],
       x[i*m3*m2+m3-1],      x[(i*m3+m2-1)*m2+m3-1],
       x[j*m2+m3-1],         x[((m1-1)*m3+j)*m2+m3-1],
       &x[(i*m3+j)*m2+m3-1] );

    }
  }
//
//  Interpolate values in the interior.
//
  for ( i = 1; i < m1 - 1; i++ )
  {

    r = ( double ) i / ( double ) ( m1 - 1 );

    for ( j = 1; j < m2 - 1; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );

      for ( k = 1; k < m3 - 1; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );

        blend_2d1 ( r, s, t,
          x[0],                    x[m3-1],
          x[(m2-1)*m2],            x[(m2-1)*m2+m3-1],
          x[(m1-1)*m3*m2],         x[(m1-1)*m3*m2+m3-1],
          x[((m1-1)*m3+m2-1)*m2],  x[((m1-1)*m3+m2-1)*m2+m3-1],
          x[i*m3*m2],              x[i*m3*m2+m3-1],
          x[(i*m3+m2-1)*m2],       x[(i*m3+m2-1)*m2+m3-1],
          x[j*m2],                 x[j*m2+m3-1],
          x[((m1-1)*m3+j)*m2],     x[((m1-1)*m3+j)*m2+m3-1],
          x[k],                    x[(m2-1)*m2+k],
          x[(m1-1)*m3*m2+k],       x[((m1-1)*m3+m2-1)*m2+k],
          x[j*m3+k],               x[((m1-1)*m3+j)*m2+k],
          x[i*m3*m2+k],            x[(i*m3+m2-1)*m2+k],
          x[(i*m3+j)*m2],          x[(i*m3+j)*m2+m3-1],
          &x[(i*m3+j)*m2+k] );

      }

    }

  }
  return;
}
//****************************************************************************80

void blend_ijk_2d1 ( double x[], int m1, int m2, int m3 )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_IJK_2D1 extends indexed scalar data along faces into a cubic table.
//
//  Diagram:
//
//    ( X000    X010    X020    X030    X040    X050 )
//    ( X100    X110    X120    X130    X140    X150 )
//    ( X200    X210    X220    X230    X240    X250 )   Layer 1
//    ( X300    X310    X320    X330    X340    X350 )
//    ( X400    X410    X420    X430    X440    X450 )
//
//    ( X001    X011    X021    X031    X041    X051 )
//    ( X101    ...     ....    ....    ....    X151 )
//    ( X201    ...     ....    ....    ....    X251 )   Layer K
//    ( X301    ...     ....    ....    ....    X351 )   1 < K < M3
//    ( X401    X411    X421    X431    X441    X451 )
//
//    ( X002    X012    X022    X032    X042    X052 )
//    ( X102    X112    X122    X132    X142    X152 )
//    ( X202    X212    X222    X232    X242    X252 )   Layer M3
//    ( X302    X312    X322    X332    X342    X352 )
//    ( X402    X412    X422    X432    X442    X452 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input/output, double X[M1*M2*M3], a singly dimensioned array that
//    is "really" triply dimensioned.  The triple dimension index
//    [I][J][K] corresponds to the single dimension index
//    I * M2*M3 + J * M2 + K
//
//    On input, there is already scalar data in the entries X[I][J][K]
//    corresponding to "faces" of the table, that is, entries for which
//    at least one of the three indices I, J and K is equal to their
//    minimum or maximum possible values.
//
//    On output, all entries in X have been assigned a value, using the
//    table indices as independent variables.
//
//    Input, int M1, M2, M3, the number of rows, columns, and layers in X.
//
{
  int i;
  int j;
  int k;
  double r;
  double s;
  double t;
//
//  Interpolate values in the interior.
//
  for ( i = 1; i < m1 - 1; i++ )
  {

    r = ( double ) i / ( double ) ( m1 - 1 );

    for ( j = 1; j < m2 - 1; j++ )
    {

      s = ( double ) j / ( double ) ( m2 - 1 );

      for ( k = 1; k < m3 - 1; k++ )
      {

        t = ( double ) k / ( double ) ( m3 - 1 );

        blend_2d1 ( r, s, t,
          x[0],                    x[m3-1],
          x[(m2-1)*m2],            x[(m2-1)*m2+m3-1],
          x[(m1-1)*m3*m2],         x[(m1-1)*m3*m2+m3-1],
          x[((m1-1)*m3+m2-1)*m2],  x[((m1-1)*m3+m2-1)*m2+m3-1],
          x[i*m3*m2],              x[i*m3*m2+m3-1],
          x[(i*m3+m2-1)*m2],       x[(i*m3+m2-1)*m2+m3-1],
          x[j*m2],                 x[j*m2+m3-1],
          x[((m1-1)*m3+j)*m2],     x[((m1-1)*m3+j)*m2+m3-1],
          x[k],                    x[(m2-1)*m2+k],
          x[(m1-1)*m3*m2+k],       x[((m1-1)*m3+m2-1)*m2+k],
          x[j*m3+k],               x[((m1-1)*m3+j)*m2+k],
          x[i*m3*m2+k],            x[(i*m3+m2-1)*m2+k],
          x[(i*m3+j)*m2],          x[(i*m3+j)*m2+m3-1],
          &x[(i*m3+j)*m2+k] );

      }

    }

  }

  return;
}
//****************************************************************************80

void blend_r_0dn ( double r, double x[], int n,
  void ( *bound_r ) ( double r, int i, double *xi ) )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_R_0DN extends vector data at endpoints into a line.
//
//  Diagram:
//
//    0-----r-----1
//
//  Discussion:
//
//    This is simply linear interpolation.  BLEND_R_0DN is provided
//    mainly as a "base routine" which can be compared to its generalizations,
//    such as BLEND_RS_0DN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, the (R) coordinate of the point to be evaluated.
//
//    Output, double X[N], the interpolated value at the point (R).
//
//    Input, int N, the dimension of the vector space.
//
//    (*void) BOUND_R(), is a subroutine which is given (R) coordinates
//    and an component value I, and returns XI, the value of the I-th component
//    of the N-vector at that point.  BOUND_R will only be called for
//    "corners", that is, for values (R) where R is either 0.0
//    or 1.0.  BOUND_R has the form:
//
//      void bound_r ( double r, int i, double *xi )
//
{
  int i;
  double x0;
  double x1;

  for ( i = 0; i < n; i++ )
  {
//
//  Get the I-th coordinate component at the two corners.
//

    bound_r ( 0.0, i, &x0 );
    bound_r ( 1.0, i, &x1 );
//
//  Interpolate the I-th coordinate component of the interior point.
//
    blend_0d1 ( r, x0, x1, &x[i] );

  }

  return;
}
//****************************************************************************80

void blend_rs_0dn ( double r, double s, double x[], int n,
  void ( *bound_rs ) ( double r, double s, int i, double *xi ) )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_RS_0DN extends vector data at corners into a square.
//
//  Diagram:
//
//    01-----r1-----11
//     |      .      |
//     |      .      |
//    0s.....rs.....1s
//     |      .      |
//     |      .      |
//    00-----r0-----10
//
//  Discussion:
//
//    BLEND_RS_0DN should be equivalent to the use of a bilinear finite
//    element method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, S, the (R,S) coordinates of the point to be evaluated.
//
//    Output, double X[N], the interpolated value at the point (R,S).
//
//    Input, int N, the dimension of the vector space.
//
//    External, BOUND_RS, is a subroutine which is given (R,S) coordinates
//    and an component value I, and returns XI, the value of the I-th component
//    of the N-vector at that point.  BOUND_RS will only be called for
//    "corners", that is, for values (R,S) where R and S are either 0.0
//    or 1.0.  BOUND_RS has the form:
//
//      void bound_rs ( double r, double s, int i, double *xi )
//
{
  int i;
  double x00;
  double x01;
  double x10;
  double x11;
  double xr0;
  double xr1;
  double x0s;
  double x1s;

  for ( i = 0; i < n; i++ )
  {
//
//  Get the I-th coordinate component at the four corners.
//
    bound_rs ( 0.0, 0.0, i, &x00 );
    bound_rs ( 0.0, 1.0, i, &x01 );
    bound_rs ( 1.0, 0.0, i, &x10 );
    bound_rs ( 1.0, 1.0, i, &x11 );
//
//  Interpolate the I-th coordinate component at the sides.
//
    blend_0d1 ( r, x00, x10, &xr0 );
    blend_0d1 ( r, x01, x11, &xr1 );
    blend_0d1 ( s, x00, x01, &x0s );
    blend_0d1 ( s, x10, x11, &x1s );
//
//  Interpolate the I-th coordinate component of the interior point.
//
    blend_1d1 ( r, s, x00, x01, x10, x11, xr0, xr1, x0s, x1s, &x[i] );

  }

  return;
}
//****************************************************************************80

void blend_rs_1dn ( double r, double s, double x[], int n,
  void ( *bound_rs ) ( double r, double s, int i, double *xi ) )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_RS_1DN extends vector data along sides into a square.
//
//  Diagram:
//
//    01-----r1-----11
//     |      .      |
//     |      .      |
//    0s.....rs.....1s
//     |      .      |
//     |      .      |
//    00-----r0-----10
//
//  Discussion:
//
//    BLEND_RS_1DN is NOT equivalent to a bilinear finite element method,
//    since the data is sampled everywhere along the boundary lines,
//    rather than at a finite number of nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, S, the (R,S) coordinates of the point to be evaluated.
//
//    Output, double X[N], the interpolated value at the point (R,S).
//
//    Input, int N, the dimension of the vector space.
//
//    External, BOUND_RS, is a subroutine which is given (R,S) coordinates
//    and an component value I, and returns XI, the value of the I-th component
//    of the N-vector at that point.  BOUND_RS will only be called for
//    "sides", that is, for values (R,S) where at least one of R and S is
//    either 0.0 or 1.0.  BOUND_RS has the form:
//
//      void bound_rs ( double r, double s, int i, double *xi )
//
{
  int i;
  double x00;
  double x01;
  double x10;
  double x11;
  double xr0;
  double xr1;
  double x0s;
  double x1s;

  for ( i = 0; i < n; i++ )
  {
//
//  Get the I-th coordinate component at the four corners.
//
    bound_rs ( 0.0, 0.0, i, &x00 );
    bound_rs ( 0.0, 1.0, i, &x01 );
    bound_rs ( 1.0, 0.0, i, &x10 );
    bound_rs ( 1.0, 1.0, i, &x11 );
//
//  Get the I-th coordinate component at the sides.
//
    bound_rs ( r, 0.0, i, &xr0 );
    bound_rs ( r, 1.0, i, &xr1 );
    bound_rs ( 0.0, s, i, &x0s );
    bound_rs ( 1.0, s, i, &x1s );
//
//  Interpolate the I-th coordinate component of the interior point.
//
    blend_1d1 ( r, s, x00, x01, x10, x11, xr0, xr1, x0s, x1s, &x[i] );

  }

  return;
}
//****************************************************************************80

void blend_rst_0dn ( double r, double s, double t, double x[], int n,
  void ( *bound_rst ) ( double r, double s, double t, int i, double *xi ) )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_RST_0DN extends vector data at corners into a cube.
//
//  Diagram:
//
//    010-----r10-----110    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    0s0.....rs0.....1s0    0s1.....rs1.....1s1     S
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    001-----r01-----101     +----R
//       BOTTOM                      TOP
//
//    011-----0s1-----001    111-----1s1-----101
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    01t.....0st.....00t    11t.....1st.....10t          T
//      |       .       |      |       .       |          |
//      |       .       |      |       .       |          |
//    010-----0s0-----000    110-----1s0-----100     S----+
//       LEFT                       RIGHT
//
//    001-----r01-----101    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    00t.....r0t.....100    01t.....r1t.....11t     T
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    010-----r10-----110     +----R
//       FRONT                       BACK
//
//  Discussion:
//
//    BLEND_RST_0DN is equivalent to a trilinear finite element method.
//    Data along the edges, faces, and interior of the cube is interpolated
//    from the data at the corners.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, S, T, the (R,S,T) coordinates of the point to be evaluated.
//
//    Output, double X[N], the interpolated value at the point (R,S,T).
//
//    Input, int N, the dimension of the vector space.
//
//    External, BOUND_RST, is a subroutine which is given (R,S,T) coordinates
//    and an component value I, and returns XI, the value of the I-th component
//    of the N-vector at that point.  BOUND_RST will only be called for
//    "corners", that is, for values (R,S,T) where R, S and T are either 0.0
//    or 1.0.  BOUND_RST has the form:
//
//      void bound_rst ( double r, double s, double t, int i, double *xi )
//
{
  int i;
  double x000;
  double x001;
  double x010;
  double x011;
  double x100;
  double x101;
  double x110;
  double x111;
  double xr00;
  double xr01;
  double xr10;
  double xr11;
  double x0s0;
  double x0s1;
  double x1s0;
  double x1s1;
  double x00t;
  double x01t;
  double x10t;
  double x11t;
  double x0st;
  double x1st;
  double xr0t;
  double xr1t;
  double xrs0;
  double xrs1;

  for ( i = 0; i < n; i++ )
  {
//
//  Get the I-th coordinate component at the corners.
//
    bound_rst ( 0.0, 0.0, 0.0, i, &x000 );
    bound_rst ( 0.0, 0.0, 1.0, i, &x001 );
    bound_rst ( 0.0, 1.0, 0.0, i, &x010 );
    bound_rst ( 0.0, 1.0, 1.0, i, &x011 );
    bound_rst ( 1.0, 0.0, 0.0, i, &x100 );
    bound_rst ( 1.0, 0.0, 1.0, i, &x101 );
    bound_rst ( 1.0, 1.0, 0.0, i, &x110 );
    bound_rst ( 1.0, 1.0, 1.0, i, &x111 );
//
//  Interpolate the I-th coordinate component at the edges.
//
    blend_0d1 ( r, x000, x100, &xr00 );
    blend_0d1 ( r, x001, x101, &xr01 );
    blend_0d1 ( r, x010, x110, &xr10 );
    blend_0d1 ( r, x011, x111, &xr11 );

    blend_0d1 ( s, x000, x010, &x0s0 );
    blend_0d1 ( s, x001, x011, &x0s1 );
    blend_0d1 ( s, x100, x110, &x1s0 );
    blend_0d1 ( s, x101, x111, &x1s1 );

    blend_0d1 ( t, x000, x001, &x00t );
    blend_0d1 ( t, x010, x011, &x01t );
    blend_0d1 ( t, x100, x101, &x10t );
    blend_0d1 ( t, x110, x111, &x11t );
//
//  Interpolate the I-th component on the faces.
//
    blend_1d1 ( s, t, x000, x001, x010, x011, x0s0, x0s1, x00t, x01t, &x0st );

    blend_1d1 ( s, t, x100, x101, x110, x111, x1s0, x1s1, x10t, x11t, &x1st );

    blend_1d1 ( r, t, x000, x001, x100, x101, xr00, xr01, x00t, x10t, &xr0t );

    blend_1d1 ( r, t, x010, x011, x110, x111, xr10, xr11, x01t, x11t, &xr1t );

    blend_1d1 ( r, s, x000, x010, x100, x110, xr00, xr10, x0s0, x1s0, &xrs0 );

    blend_1d1 ( r, s, x001, x011, x101, x111, xr01, xr11, x0s1, x1s1, &xrs1 );
//
//  Interpolate the I-th coordinate component of the interior point.
//
    blend_2d1 ( r, s, t, x000, x001, x010, x011, x100, x101,
      x110, x111, xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1,
      x00t, x01t, x10t, x11t, x0st, x1st, xr0t, xr1t, xrs0, xrs1,
      &x[i] );

  }

  return;
}
//****************************************************************************80

void blend_rst_1dn ( double r, double s, double t, double x[], int n,
  void ( *bound_rst ) ( double r, double s, double t, int i, double *xi ) )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_RST_1DN extends vector data on edges into a cube.
//
//  Diagram:
//
//    010-----r10-----110    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    0s0.....rs0.....1s0    0s1.....rs1.....1s1     S
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    001-----r01-----101     +----R
//       BOTTOM                      TOP
//
//    011-----0s1-----001    111-----1s1-----101
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    01t.....0st.....00t    11t.....1st.....10t          T
//      |       .       |      |       .       |          |
//      |       .       |      |       .       |          |
//    010-----0s0-----000    110-----1s0-----100     S----+
//       LEFT                       RIGHT
//
//    001-----r01-----101    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    00t.....r0t.....100    01t.....r1t.....11t     T
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    010-----r10-----110     +----R
//       FRONT                       BACK
//
//  Discussion:
//
//    BLEND_RST_1D is NOT equivalent to a trilinear finite element method,
//    since the data is sampled everywhere along the corners and edges,
//    rather than at a finite number of nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, S, T, the (R,S,T) coordinates of the point to be evaluated.
//
//    Output, double X(N), the interpolated value at the point (R,S,T).
//
//    Input, int N, the dimension of the vector space.
//
//    External, BOUND_RST, is a subroutine which is given (R,S,T) coordinates
//    and an component value I, and returns XI, the value of the I-th component
//    of the N-vector at that point.  BOUND_RST will only be called for
//    "edges", that is, for values (R,S,T) where at least two of R, S and T
//    are either 0.0 or 1.0.  BOUND_RST has the form:
//
//      void bound_rst ( r, s, t, i, xi )
//
{
  int i;
  double x000;
  double x001;
  double x010;
  double x011;
  double x100;
  double x101;
  double x110;
  double x111;
  double xr00;
  double xr01;
  double xr10;
  double xr11;
  double x0s0;
  double x0s1;
  double x1s0;
  double x1s1;
  double x00t;
  double x01t;
  double x10t;
  double x11t;
  double x0st;
  double x1st;
  double xr0t;
  double xr1t;
  double xrs0;
  double xrs1;

  for ( i = 0; i < n; i++ )
  {
//
//  Get the I-th coordinate component at the corners.
//
    bound_rst ( 0.0, 0.0, 0.0, i, &x000 );
    bound_rst ( 0.0, 0.0, 1.0, i, &x001 );
    bound_rst ( 0.0, 1.0, 0.0, i, &x010 );
    bound_rst ( 0.0, 1.0, 1.0, i, &x011 );
    bound_rst ( 1.0, 0.0, 0.0, i, &x100 );
    bound_rst ( 1.0, 0.0, 1.0, i, &x101 );
    bound_rst ( 1.0, 1.0, 0.0, i, &x110 );
    bound_rst ( 1.0, 1.0, 1.0, i, &x111 );
//
//  Get the I-th coordinate component at the edges.
//
    bound_rst ( r, 0.0, 0.0, i, &xr00 );
    bound_rst ( r, 0.0, 1.0, i, &xr01 );
    bound_rst ( r, 1.0, 0.0, i, &xr10 );
    bound_rst ( r, 1.0, 1.0, i, &xr11 );

    bound_rst ( 0.0, s, 0.0, i, &x0s0 );
    bound_rst ( 0.0, s, 1.0, i, &x0s1 );
    bound_rst ( 1.0, s, 0.0, i, &x1s0 );
    bound_rst ( 1.0, s, 1.0, i, &x1s1 );

    bound_rst ( 0.0, 0.0, t, i, &x00t );
    bound_rst ( 0.0, 1.0, t, i, &x01t );
    bound_rst ( 1.0, 0.0, t, i, &x10t );
    bound_rst ( 1.0, 1.0, t, i, &x11t );
//
//  Interpolate the I-th component on the faces.
//
    blend_1d1 ( s, t, x000, x001, x010, x011, x0s0, x0s1, x00t, x01t, &x0st );

    blend_1d1 ( s, t, x100, x101, x110, x111, x1s0, x1s1, x10t, x11t, &x1st );

    blend_1d1 ( r, t, x000, x001, x100, x101, xr00, xr01, x00t, x10t, &xr0t );

    blend_1d1 ( r, t, x010, x011, x110, x111, xr10, xr11, x01t, x11t, &xr1t );

    blend_1d1 ( r, s, x000, x010, x100, x110, xr00, xr10, x0s0, x1s0, &xrs0 );

    blend_1d1 ( r, s, x001, x011, x101, x111, xr01, xr11, x0s1, x1s1, &xrs1 );
//
//  Interpolate the I-th coordinate component of the interior point.
//
    blend_2d1 ( r, s, t, x000, x001, x010, x011, x100, x101,
      x110, x111, xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1,
      x00t, x01t, x10t, x11t, x0st, x1st, xr0t, xr1t, xrs0, xrs1,
      &x[i] );

  }

  return;
}
//****************************************************************************80

void blend_rst_2dn ( double r, double s, double t, double x[], int n,
  void ( *bound_rst ) ( double r, double s, double t, int i, double *xi ) )

//****************************************************************************80
//
//  Purpose:
//
//    BLEND_RST_2DN extends vector data on faces into a cube.
//
//  Diagram:
//
//    010-----r10-----110    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    0s0.....rs0.....1s0    0s1.....rs1.....1s1     S
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    001-----r01-----101     +----R
//       BOTTOM                      TOP
//
//    011-----0s1-----001    111-----1s1-----101
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    01t.....0st.....00t    11t.....1st.....10t          T
//      |       .       |      |       .       |          |
//      |       .       |      |       .       |          |
//    010-----0s0-----000    110-----1s0-----100     S----+
//       LEFT                       RIGHT
//
//    001-----r01-----101    011-----r11-----111
//      |       .       |      |       .       |
//      |       .       |      |       .       |
//    00t.....r0t.....100    01t.....r1t.....11t     T
//      |       .       |      |       .       |     |
//      |       .       |      |       .       |     |
//    000-----r00-----100    010-----r10-----110     +----R
//       FRONT                       BACK
//
//  Discussion:
//
//    BLEND_RST_2DN is NOT equivalent to a trilinear finite element method,
//    since the data is sampled everywhere along the corners, edges, and
//    faces, rather than at a finite number of nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input, double R, S, T, the (R,S,T) coordinates of the point to be evaluated.
//
//    Output, double X[N], the interpolated value at the point (R,S,T).
//
//    Input, int N, the dimension of the vector space.
//
//    External, BOUND_RST, is a subroutine which is given (R,S,T) coordinates
//    and an component value I, and returns XI, the value of the I-th component
//    of the N-vector at that point.  BOUND_RST will only be called for
//    "faces", that is, for values (R,S,T) where at least one of R, S and T
//    is either 0.0 or 1.0.  BOUND_RST has the form:
//
//      void bound_rst ( r, s, t, i, xi )
//
{
  int i;
  double x000;
  double x001;
  double x010;
  double x011;
  double x100;
  double x101;
  double x110;
  double x111;
  double xr00;
  double xr01;
  double xr10;
  double xr11;
  double x0s0;
  double x0s1;
  double x1s0;
  double x1s1;
  double x00t;
  double x01t;
  double x10t;
  double x11t;
  double x0st;
  double x1st;
  double xr0t;
  double xr1t;
  double xrs0;
  double xrs1;

  for ( i = 0; i < n; i++ )
  {
//
//  Get the I-th coordinate component at the corners.
//
    bound_rst ( 0.0, 0.0, 0.0, i, &x000 );
    bound_rst ( 0.0, 0.0, 1.0, i, &x001 );
    bound_rst ( 0.0, 1.0, 0.0, i, &x010 );
    bound_rst ( 0.0, 1.0, 1.0, i, &x011 );
    bound_rst ( 1.0, 0.0, 0.0, i, &x100 );
    bound_rst ( 1.0, 0.0, 1.0, i, &x101 );
    bound_rst ( 1.0, 1.0, 0.0, i, &x110 );
    bound_rst ( 1.0, 1.0, 1.0, i, &x111 );
//
//  Get the I-th coordinate component at the edges.
//
    bound_rst ( r, 0.0, 0.0, i, &xr00 );
    bound_rst ( r, 0.0, 1.0, i, &xr01 );
    bound_rst ( r, 1.0, 0.0, i, &xr10 );
    bound_rst ( r, 1.0, 1.0, i, &xr11 );

    bound_rst ( 0.0, s, 0.0, i, &x0s0 );
    bound_rst ( 0.0, s, 1.0, i, &x0s1 );
    bound_rst ( 1.0, s, 0.0, i, &x1s0 );
    bound_rst ( 1.0, s, 1.0, i, &x1s1 );

    bound_rst ( 0.0, 0.0, t, i, &x00t );
    bound_rst ( 0.0, 1.0, t, i, &x01t );
    bound_rst ( 1.0, 0.0, t, i, &x10t );
    bound_rst ( 1.0, 1.0, t, i, &x11t );
//
//  Get the I-th component on the faces.
//
    bound_rst ( 0.0, s, t, i, &x0st );
    bound_rst ( 1.0, s, t, i, &x1st );
    bound_rst ( r, 0.0, t, i, &xr0t );
    bound_rst ( r, 1.0, t, i, &xr1t );
    bound_rst ( r, s, 0.0, i, &xrs0 );
    bound_rst ( r, s, 1.0, i, &xrs1 );
//
//  Interpolate the I-th coordinate component of the interior point.
//
    blend_2d1 ( r, s, t, x000, x001, x010, x011, x100, x101,
      x110, x111, xr00, xr01, xr10, xr11, x0s0, x0s1, x1s0, x1s1,
      x00t, x01t, x10t, x11t, x0st, x1st, xr0t, xr1t, xrs0, xrs1,
      &x[i] );

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

void r8block_print ( int l, int m, int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8BLOCK_PRINT prints a double precision block (a 3D matrix).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, M, N, the dimensions of the block.
//
//    Input, double A[L*M*N], the matrix to be printed.
//
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;
  int j;
  int jhi;
  int jlo;
  int k;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( k = 1; k <= n; k++ )
  {
    cout << "\n";
    cout << "  K = " << k << "\n";
    cout << "\n";
    for ( jlo = 1; jlo <= m; jlo = jlo + 5 )
    {
      jhi = i4_min ( jlo + 4, m );
      cout << "\n";
      cout << "      ";
      for ( j = jlo; j <= jhi; j++ )
      {
        cout << setw(7) << j << "       ";
      }
      cout << "\n";
      cout << "\n";
      for ( i = 1; i <= l; i++ )
      {
        cout << "  " << setw(4) << i;
        for ( j = jlo; j <= jhi; j++ )
        {
          cout << "  " << setw(12) << a[i-1+(j-1)*l+(k-1)*l*m];
        }
        cout << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, char *TITLE, a title to be printed.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, char *TITLE, a title for the matrix.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  cout << "\n";

  return;
# undef INCX
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
