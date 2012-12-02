# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

int main ( );
void boundary_3d ( double umin, double vmin, double wmin, double umax, 
  double vmax, double wmax, double u, double v, double w, double *x, 
  double *y, double *z );
void sub_box_tiler_3d ( ofstream &file_out, double umin, double vmin, 
  double wmin, double umax, double vmax, double wmax, double u0, double v0, 
  double w0, double u1, double v1, double w1 );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    TILER_3D illustrates the use of 3D blending.
//
//  Discussion:
//
//    This main program works in a 3D rectangular space we can think
//    of as "UVW" space.  We are interested in the space of data values
//    bounded by (umin,vmin,wmin) and (umax,vmax,wmax), and we plan to
//    divide this up, using indices (I,J,K), into NI by NJ by NK sub-boxes.
//
//    The code below considers each sub-box indexed by (I,J,K) and determines
//    the values (u0,v0,w0) and (u1,v1,w1) that characterize its corners.
//    The coordinates of this sub-box and the coordinates of the big box
//    are then passed to SUB_BOX_TILER_3D.
//
//    The picture would be a LOT more interesting if the boundary were
//    a bit more wiggly, there were more sub-boxes, and the object in
//    each sub-box had more parts.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon and Charles Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press, 1999.
//
{
  char *file_out_name = "tiler_3d.tri";
  ofstream file_out;
  int i;
  int j;
  int k;
  int ni = 3;
  int nj = 3;
  int nk = 3;
  double u0;
  double u1;
  double umax = 30.0;
  double umin = 150.0;
  double v0;
  double v1;
  double vmax = 5.0;
  double vmin = 1.0;
  double w0;
  double w1;
  double wmax = 30.0;
  double wmin = -30.0;

  timestamp ( );

  cout << "\n";
  cout << "TILER_3D\n";
  cout << "  C++ version\n";
  cout << "  A simple example of transfinite interpolation in 3D.\n";
//
//  Write the first line of the output file, which is the number
//  of triangles in the output 3D shape.  (In our simple example,
//  each sub-box will contain a tetrahedron.)
//
  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "TILER_3D:\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  file_out << 4 * ni * nj * nk << "\n";
//
//  Consider items with index (I,*,*):
//
  for ( i = 1; i <= ni; i++ )
  {

    u0 = ( ( double ) ( ni + 1 - i ) * umin + ( double ) ( i - 1 ) * umax ) 
      / ( double ) ni;

    u1 = ( ( double ) ( ni - i ) * umin + ( double ) ( i ) * umax ) 
      / ( double ) ni;
//
//  Consider items with index (I,J,*):
//
    for ( j = 1; j <= nj; j++ )
    {

      v0 = ( ( double ) ( nj - j + 1 ) * vmin 
           + ( double ) (      j - 1 ) * vmax ) 
           / ( double )   nj;

      v1 = ( ( double ) ( nj - j ) * vmin 
           + ( double ) (      j ) * vmax ) 
           / ( double )   nj;
//
//  Consider items with index (I,J,K):
//
      for ( k = 1; k <= nk; k++ )
      {

        w0 = ( ( double ) ( nk - k + 1 ) * wmin 
             + ( double ) (      k - 1 ) * wmax ) 
             / ( double )   nk;
  
        w1 = ( ( double ) ( nk - k ) * wmin 
             + ( double ) (      k ) * wmax ) 
             / ( double )   nk;
//
//  Fill sub-box (I,J,K) with the 3-D "tile".
//   
        sub_box_tiler_3d ( file_out, umin, vmin, wmin, umax, vmax, wmax, 
          u0, v0, w0, u1, v1, w1 );
      }
    }
  }

  file_out.close ( );

  cout << "\n";
  cout << "TILER_3D:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void boundary_3d ( double umin, double vmin, double wmin, double umax, 
  double vmax, double wmax, double u, double v, double w, double *x, 
  double *y, double *z )

//****************************************************************************80
//
//  Purpose:
//
//    BOUNDARY_3D returns the (X,Y,Z) coordinates of a point (U,V,W).
//
//  Discussion:
//
//    In this example, a single formula describes the mapping
//
//      (U,V,W) => (X,Y,Z).
//
//    It is more common (and more interesting) for the formula
//    to depend on which face of the boundary is being considered.
//
//    This routine is only called for points (U,V,W) where one of the 
//    values U, V and W is "extreme", so there are generally six cases
//    to consider, for general boundaries.  The coding has been set
//    up with this in mind.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon and Charles Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press, 1999.
//
//  Parameters:
//
//    Input, double UMIN, VMIN, WMIN, UMAX, VMAX, WMAX, the (U,V,W) coordinates
//    of the two opposite corners of the big box.
//
//    Input, double U, V, W, the (U,V,W) coordinates of a point in the big box.
//
//    Output, double *X, *Y, *Z, the (X,Y,Z) coordinates of the point.
//
{
# define PI 3.14159265
# define DEG2RAD ( PI / 180.0 )

  double psi;
  double theta;

  theta = u * DEG2RAD;
  psi = w * DEG2RAD;

       if ( u == umin )
  {
    *x = v * cos ( theta ) * cos ( psi );
    *y = v * sin ( theta ) * cos ( psi );
    *z = v                 * sin ( psi );
  }
  else if ( u == umax )
  {
    *x = v * cos ( theta ) * cos ( psi );
    *y = v * sin ( theta ) * cos ( psi );
    *z = v                 * sin ( psi );
  }
  else if ( v == vmin )
  {
    *x = v * cos ( theta ) * cos ( psi );
    *y = v * sin ( theta ) * cos ( psi );
    *z = v                 * sin ( psi );
  }
  else if ( v == vmax )
  {
    *x = v * cos ( theta ) * cos ( psi );
    *y = v * sin ( theta ) * cos ( psi );
    *z = v                 * sin ( psi );
  }
  else if ( w == wmin )
  {
    *x = v * cos ( theta ) * cos ( psi );
    *y = v * sin ( theta ) * cos ( psi );
    *z = v                 * sin ( psi );
  }
  else if ( w == wmax )
  {
    *x = v * cos ( theta ) * cos ( psi );
    *y = v * sin ( theta ) * cos ( psi );
    *z = v                 * sin ( psi );
  }
  else
  {
    cout << "\n";
    cout << "BOUNDARY_3D - Fatal error!\n";
    cout << "  Illegal input value of (U,V,W).\n";
    exit ( 1 );
  }

  return;
# undef DEG2RAD
# undef PI
}
//****************************************************************************80

void sub_box_tiler_3d ( ofstream &file_out, double umin, 
  double vmin, double wmin, double umax, double vmax, double wmax, 
  double u0, double v0, double w0, double u1, double v1, double w1 )

//****************************************************************************80
//
//  Purpose:
//
//    SUB_BOX_TILER_3D "tiles" a 3D sub-box with a given pattern.
//
//  Discussion:
//
//    This routine knows the (U,V,W) coordinates of the big box and the
//    sub box, and knows the shape of the object to be place in the sub-box.
//    It uses transfinite interpolation to put the shape in the box.
//    This requires that, for each point of the shape to be mapped, the 
//    (X,Y,Z) coordinates be evaluated for points on the surface of
//    the big box, namely, at 8 corners, 12 edges, and 6 faces.  These
//    values are then blended to give a sensible (X,Y,Z) coordinate for
//    the point belonging to the shape.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon and Charles Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press, 1999.
//
//  Parameters:
//
//    Input, ofstream &file_out, a call by reference for the output file.
//
//    Input, double UMIN, VMIN, WMIN, UMAX, VMAX, WMAX, the (U,V,W) coordinates
//    of the two opposite corners of the big box.
//
//    Input, double U0, V0, W0, U1, V1, W1, the (U,V,W) coordinates of the
//    two oppositie corners of the sub-box.
//
{
# define NPOINT 4
# define R0 0.0
# define R1 1.0
# define S0 0.0
# define S1 1.0
# define T0 0.0
# define T1 1.0

  int i;
  double r;
  double r_tab[NPOINT];
  double s;
  double s_tab[NPOINT];
  double t;
  double t_tab[NPOINT];
  double u;
  double v;
  double w;
  double x[NPOINT];
  double x000;
  double x00t;
  double x001;
  double x010;
  double x011;
  double x01t;
  double x0s0;
  double x0s1;
  double x0st;
  double x100;
  double x101;
  double x10t;
  double x110;
  double x111;
  double x11t;
  double x1s0;
  double x1s1;
  double x1st;
  double xr00;
  double xr01;
  double xr0t;
  double xr10;
  double xr11;
  double xr1t;
  double xrs0;
  double xrs1;
  double y[NPOINT];
  double y000;
  double y00t;
  double y001;
  double y010;
  double y011;
  double y01t;
  double y0s0;
  double y0s1;
  double y0st;
  double y100;
  double y101;
  double y10t;
  double y110;
  double y111;
  double y11t;
  double y1s0;
  double y1s1;
  double y1st;
  double yr00;
  double yr01;
  double yr0t;
  double yr10;
  double yr11;
  double yr1t;
  double yrs0;
  double yrs1;
  double z[NPOINT];
  double z000;
  double z00t;
  double z001;
  double z010;
  double z011;
  double z01t;
  double z0s0;
  double z0s1;
  double z0st;
  double z100;
  double z101;
  double z10t;
  double z110;
  double z111;
  double z11t;
  double z1s0;
  double z1s1;
  double z1st;
  double zr00;
  double zr01;
  double zr0t;
  double zr10;
  double zr11;
  double zr1t;
  double zrs0;
  double zrs1;
//
//  Here are the (R,S,T) coordinates of the tetrahedron that we can think
//  of as the "tile" we need to place in each sub-box.
//
//  The (R,S,T) coordinate system is assumed to range from 0 to 1.
//
  r_tab[0] =  0.2;
  s_tab[0] =  0.2;
  t_tab[0] =  0.2;

  r_tab[1] =  0.8;
  s_tab[1] =  0.2;
  t_tab[1] =  0.2;

  r_tab[2] =  0.5;
  s_tab[2] =  0.8;
  t_tab[2] =  0.2;

  r_tab[3] =  0.5;
  s_tab[3] =  0.5;
  t_tab[3] =  0.8;
//
//  Compute the (X,Y,Z) coordinates of the corners of the (U,V,W) box.
//  These really only need to be computed once ever.
//
  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umin, vmin, wmin, &x000, &y000, &z000 );

  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umax, vmin, wmin, &x100, &y100, &z100 );

  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umin, vmax, wmin, &x010, &y010, &z010 );

  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umax, vmax, wmin, &x110, &y110, &z110 );

  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umin, vmin, wmax, &x001, &y001, &z001 );

  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umax, vmin, wmax, &x101, &y101, &z101 );

  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umin, vmax, wmax, &x011, &y011, &z011 );

  boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
    umax, vmax, wmax, &x111, &y111, &z111 );
//
//  Now figure out the (X,Y,Z) coordinates of the tile point with
//  given (R,S,T) coordinates.  This depends on the positions of 
//  all sorts of points on the corners, edges, and faces of the big box.
//
  for ( i = 0; i < NPOINT; i++ )
  {
//
//  Get the (R,S,T) coordinates of point I.
//
    r = r_tab[i];
    s = s_tab[i];
    t = t_tab[i];
//
//  Get the corresponding point (U,V,W) in the rectangular space.
//
    u = ( ( R1 - r      ) * u0 
        + (      r - R0 ) * u1 ) 
        / ( R1     - R0 );

    v = ( ( S1 - s      ) * v0 
        + (      s - S0 ) * v1 ) 
        / ( S1     - S0 );

    w = ( ( T1 - t      ) * w0 
        + (      t - T0 ) * w1 ) 
        / ( T1     - T0 );
//
//  Evaluate (X,Y,Z) on the 12 edges "near" (U,V,W).
//
    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umin, vmin, w, &x00t, &y00t, &z00t );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umin, vmax, w, &x01t, &y01t, &z01t );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umax, vmin, w, &x10t, &y10t, &z10t );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umax, vmax, w, &x11t, &y11t, &z11t );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umin, v, wmin, &x0s0, &y0s0, &z0s0 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umin, v, wmax, &x0s1, &y0s1, &z0s1 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umax, v, wmin, &x1s0, &y1s0, &z1s0 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umax, v, wmax, &x1s1, &y1s1, &z1s1 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, vmin, wmin, &xr00, &yr00, &zr00 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, vmin, wmax, &xr01, &yr01, &zr01 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, vmax, wmin, &xr10, &yr10, &zr10 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, vmax, wmax, &xr11, &yr11, &zr11 );
//
//  Evaluate (X,Y,Z) on the six faces near (U,V,W).
//
    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umin, v, w, &x0st, &y0st, &z0st );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      umax, v, w, &x1st, &y1st, &z1st );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, vmin, w, &xr0t, &yr0t, &zr0t );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, vmax, w, &xr1t, &yr1t, &zr1t );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, v, wmin, &xrs0, &yrs0, &zrs0 );

    boundary_3d ( umin, vmin, wmin, umax, vmax, wmax, 
      u, v, wmax, &xrs1, &yrs1, &zrs1 );
//
//  Now figure out the location of 
//    point I => (R,S,T) => (U,V,W) => (X(I),Y(I),Z(I)).
//
    x[i] = (
        ( umax - u    ) * ( vmax - v    ) * ( wmax - w    ) * x000
      - ( umax - u    ) * ( vmax - v    ) * ( wmax - wmin ) * x00t
      + ( umax - u    ) * ( vmax - v    ) * ( w    - wmin ) * x001
      - ( umax - u    ) * ( vmax - vmin ) * ( wmax - w    ) * x0s0
      + ( umax - u    ) * ( vmax - vmin ) * ( wmax - wmin ) * x0st
      - ( umax - u    ) * ( vmax - vmin ) * ( w    - wmin ) * x0s1
      + ( umax - u    ) * ( v    - vmin ) * ( wmax - w    ) * x010
      - ( umax - u    ) * ( v    - vmin ) * ( wmax - wmin ) * x01t
      + ( umax - u    ) * ( v    - vmin ) * ( w    - wmin ) * x011
      - ( umax - umin ) * ( vmax - v    ) * ( wmax - w    ) * xr00
      + ( umax - umin ) * ( vmax - v    ) * ( wmax - wmin ) * xr0t
      - ( umax - umin ) * ( vmax - v    ) * ( w    - wmin ) * xr01
      + ( umax - umin ) * ( vmax - vmin ) * ( wmax - w    ) * xrs0
      + ( umax - umin ) * ( vmax - vmin ) * ( w    - wmin ) * xrs1
      - ( umax - umin ) * ( v    - vmin ) * ( wmax - w    ) * xr10
      + ( umax - umin ) * ( v    - vmin ) * ( wmax - wmin ) * xr1t
      - ( umax - umin ) * ( v    - vmin ) * ( w    - wmin ) * xr11
      + ( u    - umin ) * ( vmax - v    ) * ( wmax - w    ) * x100
      - ( u    - umin ) * ( vmax - v    ) * ( wmax - wmin ) * x10t
      + ( u    - umin ) * ( vmax - v    ) * ( w    - wmin ) * x101
      - ( u    - umin ) * ( vmax - vmin ) * ( wmax - w    ) * x1s0
      + ( u    - umin ) * ( vmax - vmin ) * ( wmax - wmin ) * x1st
      - ( u    - umin ) * ( vmax - vmin ) * ( w    - wmin ) * x1s1
      + ( u    - umin ) * ( v    - vmin ) * ( wmax - w    ) * x110
      - ( u    - umin ) * ( v    - vmin ) * ( wmax - wmin ) * x11t
      + ( u    - umin ) * ( v    - vmin ) * ( w    - wmin ) * x111 
  ) / ( ( umax - umin ) * ( vmax - vmin ) * ( wmax - wmin ) );

    y[i] = (
        ( umax - u    ) * ( vmax - v    ) * ( wmax - w    ) * y000
      - ( umax - u    ) * ( vmax - v    ) * ( wmax - wmin ) * y00t
      + ( umax - u    ) * ( vmax - v    ) * ( w    - wmin ) * y001
      - ( umax - u    ) * ( vmax - vmin ) * ( wmax - w    ) * y0s0
      + ( umax - u    ) * ( vmax - vmin ) * ( wmax - wmin ) * y0st
      - ( umax - u    ) * ( vmax - vmin ) * ( w    - wmin ) * y0s1
      + ( umax - u    ) * ( v    - vmin ) * ( wmax - w    ) * y010
      - ( umax - u    ) * ( v    - vmin ) * ( wmax - wmin ) * y01t
      + ( umax - u    ) * ( v    - vmin ) * ( w    - wmin ) * y011
      - ( umax - umin ) * ( vmax - v    ) * ( wmax - w    ) * yr00
      + ( umax - umin ) * ( vmax - v    ) * ( wmax - wmin ) * yr0t
      - ( umax - umin ) * ( vmax - v    ) * ( w    - wmin ) * yr01
      + ( umax - umin ) * ( vmax - vmin ) * ( wmax - w    ) * yrs0
      + ( umax - umin ) * ( vmax - vmin ) * ( w    - wmin ) * yrs1
      - ( umax - umin ) * ( v    - vmin ) * ( wmax - w    ) * yr10
      + ( umax - umin ) * ( v    - vmin ) * ( wmax - wmin ) * yr1t
      - ( umax - umin ) * ( v    - vmin ) * ( w    - wmin ) * yr11
      + ( u    - umin ) * ( vmax - v    ) * ( wmax - w    ) * y100
      - ( u    - umin ) * ( vmax - v    ) * ( wmax - wmin ) * y10t
      + ( u    - umin ) * ( vmax - v    ) * ( w    - wmin ) * y101
      - ( u    - umin ) * ( vmax - vmin ) * ( wmax - w    ) * y1s0
      + ( u    - umin ) * ( vmax - vmin ) * ( wmax - wmin ) * y1st
      - ( u    - umin ) * ( vmax - vmin ) * ( w    - wmin ) * y1s1
      + ( u    - umin ) * ( v    - vmin ) * ( wmax - w    ) * y110
      - ( u    - umin ) * ( v    - vmin ) * ( wmax - wmin ) * y11t
      + ( u    - umin ) * ( v    - vmin ) * ( w    - wmin ) * y111 
    ) / ( ( umax - umin ) * ( vmax - vmin ) * ( wmax - wmin ) );

    z[i] = (
        ( umax - u    ) * ( vmax - v    ) * ( wmax - w    ) * z000
      - ( umax - u    ) * ( vmax - v    ) * ( wmax - wmin ) * z00t
      + ( umax - u    ) * ( vmax - v    ) * ( w    - wmin ) * z001
      - ( umax - u    ) * ( vmax - vmin ) * ( wmax - w    ) * z0s0
      + ( umax - u    ) * ( vmax - vmin ) * ( wmax - wmin ) * z0st
      - ( umax - u    ) * ( vmax - vmin ) * ( w    - wmin ) * z0s1
      + ( umax - u    ) * ( v    - vmin ) * ( wmax - w    ) * z010
      - ( umax - u    ) * ( v    - vmin ) * ( wmax - wmin ) * z01t
      + ( umax - u    ) * ( v    - vmin ) * ( w    - wmin ) * z011
      - ( umax - umin ) * ( vmax - v    ) * ( wmax - w    ) * zr00
      + ( umax - umin ) * ( vmax - v    ) * ( wmax - wmin ) * zr0t
      - ( umax - umin ) * ( vmax - v    ) * ( w    - wmin ) * zr01
      + ( umax - umin ) * ( vmax - vmin ) * ( wmax - w    ) * zrs0
      + ( umax - umin ) * ( vmax - vmin ) * ( w    - wmin ) * zrs1
      - ( umax - umin ) * ( v    - vmin ) * ( wmax - w    ) * zr10
      + ( umax - umin ) * ( v    - vmin ) * ( wmax - wmin ) * zr1t
      - ( umax - umin ) * ( v    - vmin ) * ( w    - wmin ) * zr11
      + ( u    - umin ) * ( vmax - v    ) * ( wmax - w    ) * z100
      - ( u    - umin ) * ( vmax - v    ) * ( wmax - wmin ) * z10t
      + ( u    - umin ) * ( vmax - v    ) * ( w    - wmin ) * z101
      - ( u    - umin ) * ( vmax - vmin ) * ( wmax - w    ) * z1s0
      + ( u    - umin ) * ( vmax - vmin ) * ( wmax - wmin ) * z1st
      - ( u    - umin ) * ( vmax - vmin ) * ( w    - wmin ) * z1s1
      + ( u    - umin ) * ( v    - vmin ) * ( wmax - w    ) * z110
      - ( u    - umin ) * ( v    - vmin ) * ( wmax - wmin ) * z11t
      + ( u    - umin ) * ( v    - vmin ) * ( w    - wmin ) * z111 
  ) / ( ( umax - umin ) * ( vmax - vmin ) * ( wmax - wmin ) );

  }

  file_out << x[0] << "  " << y[0] << "  " << z[0] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[1] << "  " << y[1] << "  " << z[1] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[2] << "  " << y[2] << "  " << z[2] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";

  file_out << x[1] << "  " << y[1] << "  " << z[1] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[0] << "  " << y[0] << "  " << z[0] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[3] << "  " << y[3] << "  " << z[3] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";

  file_out << x[2] << "  " << y[2] << "  " << z[2] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[3] << "  " << y[3] << "  " << z[3] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[0] << "  " << y[0] << "  " << z[0] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";

  file_out << x[3] << "  " << y[3] << "  " << z[3] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[2] << "  " << y[2] << "  " << z[2] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";
  file_out << x[1] << "  " << y[1] << "  " << z[1] << "  " 
           << 0.0  << "  " << 0.0  << "  " << 0.0  << "\n";

  return;
# undef NPOINT
# undef R0
# undef R1
# undef S0
# undef S1
# undef T0
# undef T1
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
