# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "polygon_monte_carlo.hpp"

//****************************************************************************80

bool between ( double xa, double ya, double xb, double yb, double xc, 
  double yc )

//****************************************************************************80
//
//  Purpose:
//
//    BETWEEN is TRUE if vertex C is between vertices A and B.
//
//  Discussion:
//
//    The points must be (numerically) collinear.
//
//    Given that condition, we take the greater of XA - XB and YA - YB
//    as a "scale" and check where C's value lies.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, the coordinates of 
//    the vertices.
//
//    Output, bool BETWEEN, is TRUE if C is between A and B.
//
{
  bool value;
  double xmax;
  double xmin;
  double ymax;
  double ymin;

  if ( ! collinear ( xa, ya, xb, yb, xc, yc ) )
  {
    value = false;
  }
  else if ( fabs ( ya - yb ) < fabs ( xa - xb ) )
  {
    xmax = r8_max ( xa, xb );
    xmin = r8_min ( xa, xb );
    value = ( xmin <= xc && xc <= xmax );
  }
  else
  {
    ymax = r8_max ( ya, yb );
    ymin = r8_min ( ya, yb );
    value = ( ymin <= yc && yc <= ymax );
  }

  return value;
}
//****************************************************************************80

bool collinear ( double xa, double ya, double xb, double yb, double xc, 
  double yc )

//****************************************************************************80
//
//  Purpose:
//
//    COLLINEAR returns a measure of collinearity for three points.
//
//  Discussion:
//
//    In order to deal with collinear points whose coordinates are not
//    numerically exact, we compare the area of the largest square
//    that can be created by the line segment between two of the points
//    to (twice) the area of the triangle formed by the points.
//
//    If the points are collinear, their triangle has zero area.
//    If the points are close to collinear, then the area of this triangle
//    will be small relative to the square of the longest segment.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, the coordinates of 
//    the vertices.
//
//    Output, bool COLLINEAR, is TRUE if the points are judged 
//    to be collinear.
//
{
  double area;
  const double r8_eps = 2.220446049250313E-016;
  double side_ab_sq;
  double side_bc_sq;
  double side_ca_sq;
  double side_max_sq;
  bool value;

  area = 0.5 * ( 
      ( xb - xa ) * ( yc - ya ) 
    - ( xc - xa ) * ( yb - ya ) );

  side_ab_sq = pow ( xa - xb, 2 ) + pow ( ya - yb, 2 );
  side_bc_sq = pow ( xb - xc, 2 ) + pow ( yb - yc, 2 );
  side_ca_sq = pow ( xc - xa, 2 ) + pow ( yc - ya, 2 );

  side_max_sq = r8_max ( side_ab_sq, r8_max ( side_bc_sq, side_ca_sq ) );

  if ( side_max_sq <= r8_eps )
  {
    value = true;
  }
  else if ( 2.0 * fabs ( area ) <= r8_eps * side_max_sq )
  {
    value = true;
  }
  else
  {
    value = false;
  }

  return value;
}
//****************************************************************************80

bool diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIAGONAL: VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int IM1, IP1, the indices of two vertices.
//
//    Input, int N, the number of vertices.
//
//    Input, int PREV[N], the previous neighbor of each vertex.
//
//    Input, int NEXT[N], the next neighbor of each vertex.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, bool DIAGONAL, the value of the test.
//
{
  bool value;
  bool value1;
  bool value2;
  bool value3;

  value1 = in_cone ( im1, ip1, n, prev, next, x, y );
  value2 = in_cone ( ip1, im1, n, prev, next, x, y );
  value3 = diagonalie ( im1, ip1, n, next, x, y );

  value = ( value1 && value2 && value3 );

  return value;
}
//****************************************************************************80

bool diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIAGONALIE is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int IM1, IP1, the indices of two vertices.
//
//    Input, int N, the number of vertices.
//
//    Input, int NEXT[N], the next neighbor of each vertex.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, bool DIAGONALIE, the value of the test.
//
{
  int first;
  int j;
  int jp1;
  bool value;
  bool value2;

  first = im1;
  j = first;
  jp1 = next[first];

  value = true;
//
//  For each edge VERTEX(J):VERTEX(JP1) of the polygon:
//
  while ( 1 )
  {
//
//  Skip any edge that includes vertex IM1 or IP1.
//
    if ( j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1 )
    {
    }
    else
    {
      value2 = intersect ( x[im1], y[im1], x[ip1], y[ip1], x[j], y[j], 
        x[jp1], y[jp1] );

      if ( value2 )
      {
        value = false;
        break;
      }
    }
    j = jp1;
    jp1 = next[j];

    if ( j == first )
    {
      break;
    }
  }

  return value;
}
//****************************************************************************80

bool in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int IM1, IP1, the indices of two vertices.
//
//    Input, int N, the number of vertices.
//
//    Input, int PREV[N], the previous neighbor of each vertex.
//
//    Input, int NEXT[N], the next neighbor of each vertex.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, bool IN_CONE, the value of the test.
//
{
  int i;
  int im2;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  bool value;

  im2 = prev[im1];
  i = next[im1];

  t1 = triangle_area ( x[im1], y[im1], x[i], y[i], x[im2], y[im2] );

  if ( 0.0 <= t1 )
  {
    t2 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2] );
    t3 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[i], y[i] );
    value = ( ( 0.0 < t2 ) && ( 0.0 < t3 ) );
  }
  else
  {
    t4 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[i], y[i] );
    t5 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2] );
    value = ! ( ( 0.0 <= t4 ) && ( 0.0 <= t5 ) );
  }
  return value;
}
//****************************************************************************80

bool intersect ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd )

//****************************************************************************80
//
//  Purpose:
//
//    INTERSECT is true if lines VA:VB and VC:VD intersect.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
//    coordinates of the four vertices.
//
//    Output, bool INTERSECT, the value of the test.
//
{
  bool value;

  if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xc, yd ) )
  {
    value = true;
  }
  else if ( between ( xa, ya, xb, yb, xc, yc ) )
  {
    value = true;
  }
  else if ( between ( xa, ya, xb, yb, xd, yd ) )
  {
    value = true;
  }
  else if ( between ( xc, yc, xd, yd, xa, ya ) )
  {
    value = true;
  }
  else if ( between ( xc, yc, xd, yd, xb, yb ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool intersect_prop ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd )

//****************************************************************************80
//
//  Purpose:
//
//    INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
//    coordinates of the four vertices.
//
//    Output, bool INTERSECT_PROP, the result of the test.
//
{
  double t1;
  double t2;
  double t3;
  double t4;
  bool value;
  bool value1;
  bool value2;
  bool value3;
  bool value4;

  if ( collinear ( xa, ya, xb, yb, xc, yc ) )
  {
    value = false;
  }
  else if ( collinear ( xa, ya, xb, yb, xd, yd ) )
  {
    value = false;
  }
  else if ( collinear ( xc, yc, xd, yd, xa, ya ) )
  {
    value = false;
  }
  else if ( collinear ( xc, yc, xd, yd, xb, yb ) )
  {
    value = false;
  }
  else
  {
    t1 = triangle_area ( xa, ya, xb, yb, xc, yc );
    t2 = triangle_area ( xa, ya, xb, yb, xd, yd );
    t3 = triangle_area ( xc, yc, xd, yd, xa, ya );
    t4 = triangle_area ( xc, yc, xd, yd, xb, yb );

    value1 = ( 0.0 < t1 );
    value2 = ( 0.0 < t2 );
    value3 = ( 0.0 < t3 );
    value4 = ( 0.0 < t4 );

    value = ( l4_xor ( value1, value2 ) ) && ( l4_xor ( value3, value4 ) );
  }
  return value;
}
//****************************************************************************80

bool l4_xor ( bool l1, bool l2 )

//****************************************************************************80
//
//  Purpose:
//
//    L4_XOR returns the exclusive OR of two L4's.
//
//  Discussion:
//
//    An L4 is a logical value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2014
//
//  Author:
//
//   John Burkardt
//
//  Parameters:
//
//    Input, bool L1, L2, two values whose exclusive OR is needed.
//
//    Output, bool L4_XOR, the exclusive OR of L1 and L2.
//
{
  bool value;
  bool value1;
  bool value2;

  value1 = (     l1   && ( ! l2 ) );
  value2 = ( ( ! l1 ) &&     l2   );

  value = ( value1 || value2 );

  return value;
}
//****************************************************************************80

double *monomial_value ( int m, int n, int e[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= i <= m ) x(i)^e(i)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, int E[M], the exponents.
//
//    Input, double X[M*N], the point coordinates.
//
//    Output, double MONOMIAL_VALUE[N], the value of the monomial.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    if ( 0 != e[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], e[i] );
      }
    }
  }

  return v;
}
//****************************************************************************80

double polygon_area ( int nv, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_AREA determines the area of a polygon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NV, the number of vertices of the polygon.
//
//    Input, double V[2*N], the vertex coordinates.
//
//    Output, double POLYGON_AREA, the area of the polygon.
//
{
  double area;
  int e[2];

  e[0] = 0;
  e[1] = 0;

  area = polygon_monomial_integral ( nv, v, e );

  return area;
}
//****************************************************************************80

double polygon_monomial_integral ( int nv, double v[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_INTEGRAL integrates a monomial over a polygon.
//
//  Discussion:
//
//    Nu(P,Q) = Integral ( x, y in polygon ) x^e(1) y^e(2) dx dy
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carsten Steger,
//    On the calculation of arbitrary moments of polygons,
//    Technical Report FGBV-96-05,
//    Forschungsgruppe Bildverstehen, Informatik IX,
//    Technische Universitaet Muenchen, October 1996.
//
//  Parameters:
//
//    Input, int NV, the number of vertices of the polygon.
//
//    Input, double V[2*N], the vertex coordinates.
//
//    Input, int E[2], the exponents of the monomial.
//
//    Output, double POLYGON_MONOMIAL_INTEGRAL, the unnormalized moment Nu(P,Q).
//
{
  int i;
  int k;
  int l;
  double nu_pq;
  int p;
  int q;
  double s_pq;
  double xi;
  double xj;
  double yi;
  double yj;

  p = e[0];
  q = e[1];

  nu_pq = 0.0;

  xj = v[0+(nv-1)*2];
  yj = v[1+(nv-1)*2];

  for ( i = 0; i < nv; i++ )
  {
    xi = v[0+i*2];
    yi = v[1+i*2];

    s_pq = 0.0;

    for ( k = 0; k <= p; k++ )
    {
      for ( l = 0; l <= q; l++ )
      {
        s_pq = s_pq 
          + r8_choose ( k + l, l ) * r8_choose ( p + q - k - l, q - l ) 
          * pow ( xi, k ) * pow ( xj, p - k ) 
          * pow ( yi, l ) * pow ( yj, q - l );
      }
    }

    nu_pq = nu_pq + ( xj * yi - xi * yj ) * s_pq;

    xj = xi;
    yj = yi;
  }

  nu_pq = nu_pq / ( double ) ( p + q + 2 ) 
    / ( double ) ( p + q + 1 ) 
    / r8_choose ( p + q, p );

  return nu_pq;
}
//****************************************************************************80

double *polygon_sample ( int nv, double v[], int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_SAMPLE uniformly samples a polygon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NV, the number of vertices.
//
//    Input, double V[2*NV], the vertices of the polygon, listed in
//    counterclockwise order.
//
//    Input, int N, the number of points to create.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, double POLYGON_SAMPLE[2*N], the points.
//
{
  double *area_cumulative;
  double area_polygon;
  double *area_relative;
  double *area_triangle;
  double area_percent;
  int i;
  int ip1;
  int j;
  int k;
  double *r;
  double *s;
  int *triangles;
  double *x;
  double *y;
//
//  Triangulate the polygon.
//
  x = new double[nv];
  y = new double[nv];
  for ( i = 0; i < nv; i++ )
  {
    x[i] = v[0+i*2];
    y[i] = v[1+i*2];
  }

  triangles = polygon_triangulate ( nv, x, y );
//
//  Determine the areas of each triangle.
//
  area_triangle = new double[nv-2];

  for ( i = 0; i < nv - 2; i++ )
  {
    area_triangle[i] = triangle_area ( 
      v[0+triangles[0+i*3]*2], v[1+triangles[0+i*3]*2], 
      v[0+triangles[1+i*3]*2], v[1+triangles[1+i*3]*2], 
      v[0+triangles[2+i*3]*2], v[1+triangles[2+i*3]*2] );
  }
//
//  Normalize the areas.
//
  area_polygon = r8vec_sum ( nv - 2, area_triangle );

  area_relative = new double[nv-2];

  for ( i = 0; i < nv - 2; i++ )
  {
    area_relative[i] = area_triangle[i] / area_polygon;
  }
//
//  Replace each area by the sum of itself and all previous ones.
//
  area_cumulative = new double[nv-2];

  area_cumulative[0] = area_relative[0];
  for ( i = 1; i < nv - 2; i++ )
  {
    area_cumulative[i] = area_relative[i] + area_cumulative[i-1];
  }

  s = new double[2*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Choose triangle I at random, based on areas.
//
    area_percent = r8_uniform_01 ( seed );

    for ( k = 0; k < nv - 2; k++ )
    {
      i = k;

      if ( area_percent <= area_cumulative[k] )
      {
        break;
      }
    }
//
//  Now choose a point at random in triangle I.
//
    r = r8vec_uniform_01_new ( 2, seed );

    if ( 1.0 < r[0] + r[1] )
    {
      r[0] = 1.0 - r[0];
      r[1] = 1.0 - r[1];
    }

    s[0+j*2] = ( 1.0 - r[0] - r[1] ) * v[0+triangles[0+i*3]*2]
                     + r[0]          * v[0+triangles[1+i*3]*2]
                            + r[1]   * v[0+triangles[2+i*3]*2];

    s[1+j*2] = ( 1.0 - r[0] - r[1] ) * v[1+triangles[0+i*3]*2]
                     + r[0]          * v[1+triangles[1+i*3]*2]
                            + r[1]   * v[1+triangles[2+i*3]*2];
    delete [] r;
  }

  delete [] area_cumulative;
  delete [] area_relative;
  delete [] area_triangle;
  delete [] triangles;
  delete [] x;
  delete [] y;

  return s;
}
//****************************************************************************80

int *polygon_triangulate ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_TRIANGULATE determines a triangulation of a polygon.
//
//  Discussion:
//
//    There are N-3 triangles in the triangulation.
//
//    For the first N-2 triangles, the first edge listed is always an
//    internal diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int N, the number of vertices.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, int TRIANGLES[3*(N-2)], the triangles of the 
//    triangulation.
//
{
  double area;
  bool *ear;
  int first;
  int i;
  int i0;
  int i1;
  int i2;
  int i3;
  int i4;
  int *next;
  int node;
  int node_m1;
  int *prev;
  int triangle_num;
  int *triangles;
//
//  We must have at least 3 vertices.
//
  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_TRIANGULATE - Fatal error!\n";
    cerr << "  N < 3.\n";
    exit ( 1 );
  }
//
//  Consecutive vertices cannot be equal.
//
  node_m1 = n - 1;
  for ( node = 0; node < n; node++ )
  {
    if ( x[node_m1] == x[node] && y[node_m1] == y[node] )
    {
      cerr << "\n";
      cerr << "POLYGON_TRIANGULATE - Fatal error!\n";
      cerr << "  Two consecutive nodes are identical.\n";
      exit ( 1 );
    }
    node_m1 = node;
  }
//
//  Area must be positive.
//
  area = 0.0;
  for ( node = 0; node < n - 2; node++ )
  {
    area = area + 0.5 * 
    ( 
        ( x[node+1] - x[node] ) * ( y[node+2] - y[node] ) 
      - ( x[node+2] - x[node] ) * ( y[node+1] - y[node] )
    );
  }

  if ( area <= 0.0 )
  {
    cerr << "\n";
    cerr << "POLYGON_TRIANGULATE - Fatal error!\n";
    cerr << "  Polygon has zero or negative area.\n";
    exit ( 1 );
  }

  triangles = new int[3*(n-2)];
//
//  PREV and NEXT point to the previous and next nodes.
//
  prev = new int[n];
  next = new int[n];

  i = 0;
  prev[i] = n - 1;
  next[i] = i + 1;

  for ( i = 1; i < n - 1; i++ )
  {
    prev[i] = i - 1;
    next[i] = i + 1;
  }

  i = n - 1;
  prev[i] = i - 1;
  next[i] = 0;
//
//  EAR indicates whether the node and its immediate neighbors form an ear
//  that can be sliced off immediately.
//
  ear = new bool[n];
  for ( i = 0; i < n; i++ )
  {
    ear[i] = diagonal ( prev[i], next[i], n, prev, next, x, y );
  }

  triangle_num = 0;

  i2 = 0;

  while ( triangle_num < n - 3 )
  {
//
//  If I2 is an ear, gather information necessary to carry out
//  the slicing operation and subsequent "healing".
//
    if ( ear[i2] )
    {
      i3 = next[i2];
      i4 = next[i3];
      i1 = prev[i2];
      i0 = prev[i1];
//
//  Make vertex I2 disappear.
//
      next[i1] = i3;
      prev[i3] = i1;
//
//  Update the earity of I1 and I3, because I2 disappeared.
//
      ear[i1] = diagonal ( i0, i3, n, prev, next, x, y );
      ear[i3] = diagonal ( i1, i4, n, prev, next, x, y );
//
//  Add the diagonal [I3, I1, I2] to the list.
//
      triangles[0+triangle_num*3] = i3;
      triangles[1+triangle_num*3] = i1;
      triangles[2+triangle_num*3] = i2;
      triangle_num = triangle_num + 1;
    }
//
//  Try the next vertex.
//
    i2 = next[i2];
  }
//
//  The last triangle is formed from the three remaining vertices.
//
  i3 = next[i2];
  i1 = prev[i2];

  triangles[0+triangle_num*3] = i3;
  triangles[1+triangle_num*3] = i1;
  triangles[2+triangle_num*3] = i2;
  triangle_num = triangle_num + 1;

  delete [] ear;
  delete [] next;
  delete [] prev;

  return triangles;
}
//****************************************************************************80

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  double value;

  if ( k < n - k )
  {
    mn = k;
    mx = n - k;
  }
  else
  {
    mn = n - k;
    mx = k;
  }

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
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
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
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
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
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
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  const int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
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
//    19 August 2004
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
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
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

double triangle_area ( double xa, double ya, double xb, double yb, double xc, 
  double yc )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA computes the signed area of a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
//    the vertices of the triangle, given in counterclockwise order.
//
//    Output, double TRIANGLE_AREA, the signed area of the triangle.
//
{
  double value;

  value = 0.5 * ( 
      ( xb - xa ) * ( yc - ya ) 
    - ( xc - xa ) * ( yb - ya ) );

  return value;
}
