# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "circle_segment.hpp"

//****************************************************************************80

double circle_segment_angle_from_chord ( double r, double c[2], double p1[2], 
  double p2[2] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_ANGLE_FROM_CHORD computes the angle of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double P1[2], P2[2], the ends of the chord.
//
//    Output, double CIRCLE_SEGMENT_ANGLE_FROM_CHORD, the value of THETA,
//    the angle of the circle segment.  0 <= THETA < 2 * PI.
//
{
  const double pi = 3.141592653589793;
  double theta;
  double v1[2];
  double v2[2];
//
//  Compute the radial vectors V1 and V2.
//
  v1[0] = p1[0] - c[0];
  v1[1] = p1[1] - c[1];
  v2[0] = p2[0] - c[0];
  v2[1] = p2[1] - c[1];
//
//  The arc cosine will only give us an answer between 0 and PI.
//
  theta = r8_atan ( v2[1], v2[0] ) - r8_atan ( v1[1], v1[0] );
//
//  Force 0 <= THETA < 2 * PI.
//
  while ( theta < 0.0 )
  {
    theta = theta + 2.0 * pi;
  }

  while ( 2.0 * pi <= theta )
  {
    theta = theta - 2.0 * pi;
  }

  return theta;
}
//****************************************************************************80

double circle_segment_angle_from_chord_angles ( double omega1, double omega2 )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_ANGLE_FROM_CHORD_ANGLES computes angle of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double OMEGA1, OMEGA2, the angles of the points P1 
//    and P2.  OMEGA1 <= OMEGA2.
//
//    Output, double CIRCLE_SEGMENT_ANGLE_FROM_CHORD_ANGLES, the angle THETA
//    of the circle segment.  Essentially, THETA = OMEGA2 - OMEGA1.
//
{
  const double pi = 3.141592653589793;
  double theta;

  while ( omega2 < omega1 )
  {
    omega2 = omega2 + 2.0 * pi;
  }

  theta = omega2 - omega1;

  return theta;
}
//****************************************************************************80

double circle_segment_angle_from_height ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT computes the angle of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double H, the "height" of the circle segment.
//    0 <= H <= 2 * R.
//
//    Output, double CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT, the angle THETA
//    of the circle segment.
//
{
  const double pi = 3.141592653589793;
  double theta;

  if ( h <= 0.0 )
  {
    theta = 0.0;
  }
  else if ( h <= r )
  {
    theta = 2.0 * r8_acos ( ( r - h ) / r );
    theta = 2.0 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
  }
  else if ( h <= 2.0 * r )
  {
    theta = 2.0 * r8_acos ( ( r - h ) / r );
    theta = 2.0 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
    theta = 2.0 * pi - theta;
  }
  else
  {
    theta = 2.0 * pi;
  }

  return theta;
}
//****************************************************************************80

double circle_segment_area_from_angle ( double r, double theta )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_AREA_FROM_ANGLE computes the area of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double THETA, the angle of the circle segment.
//
//    Output, double CIRCLE_SEGMENT_AREA_FROM_ANGLE, the area of the 
//    circle segment.
//
{
  double area;

  area = r * r * ( theta - sin ( theta ) ) / 2.0;

  return area;
}
//****************************************************************************80

double circle_segment_area_from_chord ( double r, double c[2], double p1[2], 
  double p2[2] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_AREA_FROM_CHORD computes the area of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double P1[2], P2[2], the ends of the chord.
//
//    Output, double CIRCLE_SEGMENT_AREA_FROM_CHORD, the area of the 
//    circle segment.
//
{
  double area;
  double theta;

  theta = circle_segment_angle_from_chord ( r, c, p1, p2 );

  area = r * r * ( theta - sin ( theta ) ) / 2.0;

  return area;
}
//****************************************************************************80

double circle_segment_area_from_height ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double H, the height of the circle segment.
//    0 <= H <= 2 * R.
//
//    Output, double CIRCLE_SEGMENT_AREA_FROM_HEIGHT, the area of the 
//    circle segment.
//
{
  double area;
  const double pi = 3.141592653589793;
  double theta;

  if ( h <= 0.0 )
  {
    area = 0.0;
  }
  else if ( h <= r )
  {
    theta = 2.0 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
    area = r * r * ( theta - sin ( theta ) ) / 2.0;
  }
  else if ( h <= 2.0 * r )
  {
    theta = 2.0 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
    theta = 2.0 * pi - theta;
    area = r * r * ( theta - sin ( theta ) ) / 2.0;
  }
  else
  {
    area = pi * r * r;
  }

  return area;
}
//****************************************************************************80

double circle_segment_area_from_sample ( double r, double c[2], double p1[2], 
  double p2[2], int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_AREA_FROM_SAMPLE computes the area of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double P1[2], P2[2], the ends of the chord.
//
//    Input, int N, the number of sample points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double CIRCLE_SEGMENT_AREA_FROM_SAMPLE, the area of the 
//    circle segment.
//
{
  double *angle;
  double area;
  int i;
  int m;
  double omega1;
  double omega2;
  double p[2];
  const double pi = 3.141592653589793;
  double *r2;
  double rmh;
  double *vdotp;
  double *x;
  double *y;
//
//  Determine the angles of the chord endpoints.
//
  omega1 = r8_atan ( p1[1] - c[1], p1[0] - c[0] );
  while ( omega1 < 0.0 )
  {
    omega1 = omega1 + 2.0 * pi;
  }

  omega2 = r8_atan ( p2[1] - c[1], p2[0] - c[0] );
  while ( omega2 < omega1 )
  {
    omega2 = omega2 + 2.0 * pi;
  }
//
//  Get N random points in the circle.
//  To simplify angle measurement, take OMEGA1 as your smallest angle.
//  That way, the check OMEGA1 <= ANGLE <= OMEGA2 will be legitimate.
//
  angle = r8vec_uniform_01_new ( n, seed );
  for ( i = 0; i < n; i++ )
  {
    angle[i] = omega1 + 2.0 * pi * angle[i];
  }
  r2 = r8vec_uniform_01_new ( n, seed );
  for ( i = 0; i < n; i++ )
  {
    r2[i] = sqrt ( r2[i] );
  }
  x = new double[n];
  y = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = c[0] + r2[i] * cos ( angle[i] );
    y[i] = c[0] + r2[i] * sin ( angle[i] );
  }
//
//  Determine the vector that touches the circle segment base.
//
  p[0] = 0.5 * ( p1[0] + p2[0] ) - c[0];
  p[1] = 0.5 * ( p1[1] + p2[1] ) - c[1];

  rmh = sqrt ( p[0] * p[0] + p[1] * p[1] );

  p[0] = p[0] / rmh;
  p[1] = p[1] / rmh;

  if ( pi < omega2 - omega1 )
  {
    p[0] = - p[0];
    p[1] = - p[1];
    rmh =  - rmh;
  }
//
//  Compute the projection of each point onto P.
//
  vdotp = new double[n];
  for ( i = 0; i < n; i++ )
  {
    vdotp[i] = ( x[i] - c[0] ) * p[0] + ( y[i] - c[1] ) * p[1];
  }
//
//  Points in the segment lie in the sector, and project at least RMH onto P.
//
  m = 0;
  for ( i = 0; i < n; i++ )
  {
    if ( omega1 < angle[i] && angle[i] < omega2 && rmh < vdotp[i] )
    {
      m = m + 1;
    }
  }
//
//  The area of the segment is its relative share of the circle area.
//
  area = pi * r * r * ( double ) ( m ) / ( double ) ( n );

  delete [] angle;
  delete [] r2;
  delete [] vdotp;
  delete [] x;
  delete [] y;

  return area;
}
//****************************************************************************80

double circle_segment_cdf ( double r, double h, double h2 )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_CDF computes a CDF related to a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//    Now, suppose we want to assign a cumulative density function or CDF
//    based on a variable H2 which measures the height of the circle segment
//    formed by an arbitrary point in the circle segment.  CDF(H2) will
//    measure the probability that a point drawn uniformly at random from
//    the circle segment defines a (smaller) circle segment of height H2.
//
//    If we can define this CDF, then we will be able to sample uniformly
//    from the circle segment, since our "Y" value can be determined from H2,
//    and our X value is chosen uniformly over the horizontal chord 
//    through (0,Y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double H, the "height" of the circle segment.
//    0 <= H <= 2 * R.
//
//    Input, double H2, the "height" of the new circle segment 
//    defined by a given point in the circle segment.  0 <= H2 <= H.
//
//    Output, double CDF, the cumulative density function for H2, 
//    the probability that a point chosen at random in the circle segment 
//    would define a smaller circle segment of height H2 or less.
//
{
  double a;
  double a2;
  double cdf;

  if ( h2 <= 0.0 )
  {
    cdf = 0.0;
  }
  else if ( h <= h2 )
  {
    cdf = 1.0;
  }
  else
  {
    a = circle_segment_area_from_height ( r, h  );
    a2 = circle_segment_area_from_height ( r, h2 );
    cdf = a2 / a;
  }

  return cdf;
}
//****************************************************************************80

double *circle_segment_centroid_from_chord ( double r, double c[2], 
  double p1[2], double p2[2] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_CENTROID_FROM_CHORD computes the centroid of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//    For this function, we assume that the center of the circle is at (0,0),
//    that the chord is horizontal, and that the circle segment is at the top.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double P1[2], P2[2], the coordinates of the endpoints 
//    of the chord.
//
//    Output, double CIRCLE_SEGMENT_CENTROID_FROM_CHORD[2], the coordinates 
//    of the centroid.
//
{
  double *d;
  double s;
  double theta;
  double thetah;
  double v1[2];
//
//  Get the angle subtended by P1:P2.
//
  theta = circle_segment_angle_from_chord ( r, c, p1, p2 );
//
//  Construct V1, the vector from C to P1.
//
  v1[0] = p1[0] - c[0];
  v1[1] = p1[1] - c[1];
//
//  Rotate V1 through THETA / 2.
//
  thetah = theta / 2.0;

  d = new double[2];
  d[0] = cos ( thetah ) * v1[0] - sin ( thetah ) * v1[1];
  d[1] = sin ( thetah ) * v1[0] + cos ( thetah ) * v1[1];
//
//  Scale this vector so it represents the distance to the centroid
//  relative to R.
//
  s = 4.0 * pow ( sin ( theta / 2.0 ), 3 ) 
    / 3.0 / ( theta - sin ( theta ) );

  d[0] = s * d[0];
  d[1] = s * d[1];
//
//  Add the center.
//
  d[0] = d[0] + c[0];
  d[1] = d[1] + c[1];

  return d;
}
//****************************************************************************80

double *circle_segment_centroid_from_height ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_CENTROID_FROM_HEIGHT computes centroid of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//    For this function, we assume that the center of the circle is at (0,0).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double H, the "height" of the circle segment.
//    0 <= H <= 2 * R.
//
//    Output, double CIRCLE_SEGMENT_CENTROID_FROM_HEIGHT[2], the coordinates 
//    of the centroid.
//
{
  double *d;
  double theta;

  theta = circle_segment_angle_from_height ( r, h );

  d = new double[2];

  d[0] = 0.0;
  d[1] = 4.0 * r * pow ( sin ( theta / 2.0 ), 3 ) / 3.0 
    / ( theta - sin ( theta ) );

  return d;
}
//****************************************************************************80

double *circle_segment_centroid_from_sample ( double r, double c[2], 
  double p1[2], double p2[2], int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE estimates a circle segment centroid.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double P1[2], P2[2], the ends of the chord.
//
//    Input, int N, the number of sample points.
//
//    Input/output, int *&EED, a seed for the random 
//    number generator.
//
//    Output, double CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE[2], the estimated 
//    centroid of the circle segment.
//
{
  double *d;
  double *x;
  double *y;

  x = new double[n];
  y = new double[n];

  circle_segment_sample_from_chord ( r, c, p1, p2, n, seed, x, y );

  d = new double[2];

  d[0] = r8vec_sum ( n, x ) / ( double ) ( n );
  d[1] = r8vec_sum ( n, y ) / ( double ) ( n );

  delete [] x;
  delete [] y;

  return d;
}
//****************************************************************************80

int circle_segment_contains_point ( double r, double c[2], double omega1, 
  double omega2, double xy[2] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_CONTAINS_POINT reports whether a point is in a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//    In this function, we allow the circle to have an arbitrary center C,
//    arbitrary radius R, and we describe the points P1 and P2 by specifying
//    their angles OMEGA1 and OMEGA2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double OMEGA1, OMEGA2, the angles of the two points on 
//    the circumference of the circle that define the circle segment.
//    OMEGA1 < OMEGA2 <= OMEGA1 + 2 * PI
//
//    Input, double XY[2], a point.
//
//    Output, int CIRCLE_SEGMENT_CONTAINS_POINT, is TRUE if the point is inside 
//    the circle segment.
//
{
  double h;
  double omegah;
  const double pi = 3.141592653589793;
  double theta;
  double v[2];
  double v_omega;
  double v_project;
  double v_r;
  int value;

  if ( r <= 0.0 )
  {
    cerr << "\n";
    cerr << "CIRCLE_SEGMENT_CONTAINS_POINT - Fatal error!\n";
    cerr << "  R <= 0.0.\n";
    exit ( 1 );
  }

  while ( omega2 < omega1 )
  {
    omega2 = omega2 + 2.0 * pi;
  }
//
//  Compute the vector V = XY - C:
//
  v[0] = xy[0] - c[0];
  v[1] = xy[1] - c[1];
//
//  a: Point must be inside the circle, so ||V|| <= R.
//
  v_r = sqrt ( v[0] * v[0] + v[1] * v[1] );

  if ( r < v_r )
  {
    value = 0;
    return value;
  }
//
//  b: Angle made by the vector V must be between OMEGA1 and OMEGA2.
//
  v_omega = atan2 ( v[1], v[0] );

  while ( omega1 <= v_omega + 2.0 * pi )
  {
    v_omega = v_omega - 2.0 * pi;
  }

  while ( v_omega + 2.0 * pi <= omega1 )
  {
    v_omega = v_omega + 2.0 * pi;
  }

  if ( omega2 < v_omega )
  {
    value = 0;
    return value;
  }
//
//  c: Projection of V onto unit centerline must be at least R-H.
//
  omegah = 0.5 * ( omega1 + omega2 );
  v_project = v[0] * cos ( omegah ) + v[1] * sin ( omegah );

  theta = omega2 - omega1;
  h = circle_segment_height_from_angle ( r, theta );

  if ( v_project < r - h )
  {
    value = 0;
    return value;
  }

  value = 1;
  
  return value;
}
//****************************************************************************80

double circle_segment_height_from_angle ( double r, double angle )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE: height of a circle segment from angle.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//    This function is given the radius R and angle of the segment, and
//    determines the corresponding height.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double ANGLE, the angle of the circle segment.
//    0 <= ANGLE <= 2.0 * PI.
//
//    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE, the height of the 
//    circle segment.
//
{
  double h;
  const double pi = 3.141592653589793;

  if ( angle < 0.0 )
  {
    cerr << "\n";
    cerr << "CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!\n";
    cerr << "  ANGLE < 0.0.\n";
    exit ( 1 );
  }

  if ( angle == 0.0 )
  {
    h = 0.0;
    return h;
  }

  if ( angle == 2.0 * pi )
  {
    h = 2.0 * r;
    return h;
  }

  if ( 2.0 * pi < angle )
  {
    cerr << "\n";
    cerr << "CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!\n";
    cerr << "  2.0 * pi < ANGLE.\n";
    exit ( 1 );
  }

  if ( r <= 0.0 )
  {
    cerr << "\n";
    cerr << "CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!\n";
    cerr << "  R <= 0.0.\n";
    exit ( 1 );
  }

  if ( angle <= pi )
  {
    h = r * ( 1.0 - cos (              angle   / 2.0 ) );
  }
  else
  {
    h = r * ( 1.0 + cos ( ( 2.0 * pi - angle ) / 2.0 ) );
  }

  return h;
}
//****************************************************************************80

double circle_segment_height_from_area ( double r, double area )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_HEIGHT_FROM_AREA: height of a circle segment from area.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//    This function is given the radius R and area of the segment, and
//    determines the corresponding height.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double AREA, the area of the circle segment.
//    0 <= AREA <= 2.0 * PI * R^2.
//
//    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_AREA, the height of the
//    circle segment.
//
{
  double a;
  double a1;
  double a2;
  double area_circle;
  double eps;
  double h;
  double h1;
  double h2;
  int it;
  const double pi = 3.141592653589793;

  if ( area < 0.0 )
  {
    cerr << "\n";
    cerr << "CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!\n";
    cerr << "  AREA < 0.0.\n";
    exit ( 1 );
  }

  area_circle = 2.0 * pi * r * r;

  if ( area == 0.0 )
  {
    h = 0.0;
    return h;
  }

  if ( area == area_circle )
  {
    h = 2.0 * r;
    return h;
  }

  if ( area_circle < area )
  {
    cerr << "\n";
    cerr << "CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!\n";
    cerr << "  2.0 * pi * r^2 < AREA.\n";
    exit ( 1 );
  }

  if ( r <= 0.0 )
  {
    cerr << "\n";
    cerr << "CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!\n";
    cerr << "  R <= 0.0.\n";
    exit ( 1 );
  }

  h1 = 0.0;
  a1 = circle_segment_area_from_height ( r, h1 );
  h2 = 2.0 * r;
  a2 = circle_segment_area_from_height ( r, h2 );

  it = 0;
  eps = r8_epsilon ( );

  while ( it < 30 )
  {
    h = 0.5 * ( h1 + h2 );
    a = circle_segment_area_from_height ( r, h );
    it = it + 1;

    if ( fabs ( a - area ) < sqrt ( eps ) * area )
    {
      break;
    }

    if ( a < area )
    {
      h1 = h;
      a1 = a;
    }
    else
    {
      h2 = h;
      a2 = a;
    }
  }

  return h;
}
//****************************************************************************80

double circle_segment_height_from_chord ( double r, double c[2], double p1[2], 
  double p2[2] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_HEIGHT_FROM_CHORD: height of a circle segment from chord.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the coordinates of the circle center.
//
//    Input, double P1[2], P2[2], the coordinates of the 
//    chord endpoints.
//
//    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_CHORD, the height of the circle segment.
//
{
  double h;
  double theta;

  theta = circle_segment_angle_from_chord ( r, c, p1, p2 );

  h = circle_segment_height_from_angle ( r, theta );

  return h;
}
//****************************************************************************80

double circle_segment_rotation_from_chord ( double r, double c[], double p1[], 
  double p2[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_ROTATION_FROM_CHORD computes the rotation of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double P1[2], P2[2], the ends of the chord.
//    Warning! If P1 = P2, we can't tell whether the segment is the whole
//    circle or none of it!
//
//    Output, double CIRCLE_SEGMENT_ROTATION_FROM_CHORD, the rotation of the 
//    circle segment.  0 <= ALPHA < 2 * PI.
//
{
  double alpha;
  double pi = 3.141592653589793;
  double rho1;
  double rho2;
  double theta;
  double v1[2];
  double v2[2];
//
//  Compute the radial vectors V1 and V2.
//
  v1[0] = p1[0] - c[0];
  v1[1] = p1[1] - c[1];
  v2[0] = p2[0] - c[0];
  v2[1] = p2[1] - c[1];
//
//  Use R8_ATAN to guarantee that 0 <= RHO1, RHO2 <= 2 * PI.
//
  rho1 = r8_atan ( v1[1], v1[0] );
  rho2 = r8_atan ( v2[1], v2[0] );
//
//  Force RHO2 to be bigger than RHO1.
//
  while ( rho2 <= rho1 )
  {
    rho2 = rho2 + 2.0 * pi;
  }
//
//  Compute THETA.
//
  theta = rho2 - rho1;
//
//  ALPHA is RHO1, plus half of the angular distance between P1 and P2.
//
  alpha = rho1 + 0.5 * theta;

  while ( 2.0 * pi <= alpha )
  {
    alpha = alpha - 2.0 * pi;
  }

  return alpha;
}
//****************************************************************************80

void circle_segment_sample_from_chord ( double r, double c[2], double p1[2], 
  double p2[2], int n, int &seed, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_SAMPLE_FROM_CHORD samples points from a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double C[2], the center of the circle.
//
//    Input, double P1[2], P2[2], the endpoints of the chord.
//
//    Input, int N, the number of sample points.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double X[N], Y[N], the sample points.
//
{
  double c2[2];
  double *eta;
  double h;
  int i;
  double t;
  double vc[2];
  double vr[2];
  double *xi;
//
//  Determine unit vectors VR and VC.
//  VR points to the center of the chord from the radius.
//  VC points along the chord, from P1 to P2.
//
  vr[0] = 0.5 * ( p1[0] + p2[0] ) - c[0];
  vr[1] = 0.5 * ( p1[1] + p2[1] ) - c[1];

  t = sqrt ( vr[0] * vr[0] + vr[1] * vr[1] );
  vr[0] = vr[0] / t;
  vr[1] = vr[1] / t;

  vc[0] = p2[0] - p1[0];
  vc[1] = p2[1] - p1[1];

  t = sqrt ( vc[0] * vc[0] + vc[1] * vc[1] );
  vc[0] = vc[0] / t;
  vc[1] = vc[1] / t;
//
//  Get the height of the circle segment.
//
  c2[0] = 0.0;
  c2[1] = 0.0;
  h = circle_segment_height_from_chord ( r, c2, p1, p2 );
//
//  Sample (xi,eta) in the reference coordinates, where the chord
//  is horizontal.
//
  xi = new double[n];
  eta = new double[n];
  circle_segment_sample_from_height ( r, h, n, seed, xi, eta );
//
//  XI is the left/right coordinate along VC.
//  ETA is the distance along VR.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = c[0] + eta[i] * vr[0] + xi[i] * vc[0];
    y[i] = c[1] + eta[i] * vr[1] + xi[i] * vc[1];
  }

  delete [] eta;
  delete [] xi;

  return;
}
//****************************************************************************80

void circle_segment_sample_from_height ( double r, double h, int n, int &seed, 
  double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples points from a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double H, the height of the circle segment.
//    0 <= H <= 2 * R.
//
//    Input, int N, the number of sample points.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double X[N], Y[N], the sample points.
//
{
  double area;
  double *area2;
  double *h2;
  int i;
  double *u;
  double *wh;

  area = circle_segment_area_from_height ( r, h );
//
//  Pick CDF's randomly.
//
  u = r8vec_uniform_01_new ( n, seed );
//
//  Choose points randomly by choosing ordered areas randomly.
//
  area2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    area2[i] = u[i] * area;
  }
//
//  Each area corresponds to a height H2.  Find it.
//
  h2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    h2[i] = circle_segment_height_from_area ( r, area2[i] );
  }
//
//  Determine the half-width WH of the segment for each H2.
//
  wh = new double[n];
  for ( i = 0; i < n; i++ )
  {
    wh[i] = sqrt ( h2[i] * ( 2.0 * r - h2[i] ) );
  }
//
//  Choose an X position randomly in [-WH,+WH].
//
  delete [] u;
  u = r8vec_uniform_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( 2.0 * u[i] - 1.0 ) * wh[i];
  }
//
//  Our circle center is at (0,0).  Our height of H2 is subtracted
//  from the height R at the peak of the circle.  Determine the Y
//  coordinate using this fact.
//
  for ( i = 0; i < n; i++ )
  {
    y[i] = r - h2[i];
  }

  delete [] area2;
  delete [] h2;
  delete [] u;
  delete [] wh;

  return;
}
//****************************************************************************80

double circle_segment_width_from_height ( double r, double h )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_WIDTH_FROM_HEIGHT computes the width of a circle segment.
//
//  Discussion:
//
//    Begin with a circle of radius R.  Choose two points P1 and P2 on the
//    circle, and draw the chord P1:P2.  This chord divides the circle
//    into two pieces, each of which is called a circle segment.
//    Consider one of the pieces.  The "angle" of this segment is the angle 
//    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
//    the chord P1:P2 which is closest to C.  The "height" of the segment
//    is the distance from Q to the perimeter of the circle.  The "width"
//    of the circle segment is the length of P1:P2.
//
//    This function is given the radius R and height H of the segment, and
//    determines the corresponding width W.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the circle.
//    0 < R.
//
//    Input, double H, the height of the circle segment.
//    0 <= H <= 2 * R.
//
//    Output, double CIRCLE_SEGMENT_WIDTH_FROM_HEIGHT, the width of the 
//    circle segment.
//
{
  double w;

  w = 2.0 * sqrt ( h * ( 2.0 * r - h ) );

  return w;
}
//****************************************************************************80

void filename_inc ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the name is empty, then the routine stops.
//
//    If the name contains no digits, the empty string is returned.
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
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_INC - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

void gauss ( int n, double alpha[], double beta[], double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GAUSS computes a Gauss quadrature rule.
//
//  Discussion:
//
//    Given a weight function W encoded by the first N recurrence coefficients 
//    ALPHA and BETA for the associated orthogonal polynomials, the call 
//      call gauss ( n, alpha, beta, x, w ) 
//    generates the nodes and weights of the N-point Gauss quadrature rule 
//    for the weight function W.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2013
//
//  Author:
//
//    Original MATLAB version by Walter Gautschi.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//
//  Parameters:
//
//    Input, int N, the order of the desired quadrature rule.
//
//    Input, double ALPHA[N], BETA[N], the alpha and beta recurrence 
//    coefficients for the othogonal polynomials associated with the
//    weight function.
//
//    Output, double X[N], W[N], the nodes and  weights of the desired 
//    quadrature rule.  The nodes are listed in increasing order.
//
{
  double *a;
  int i;
  int it_max;
  int it_num;
  int j;
  int rot_num;
  double *v;
//
//  Define the tridiagonal Jacobi matrix.
//
  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = alpha[i];
      }
      else if ( i == j - 1 )
      {
        a[i+j*n] = sqrt ( beta[j] );
      }
      else if ( i - 1 == j )
      {
        a[i+j*n] = sqrt ( beta[i] );
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }
//
//  Get the eigenvectors and eigenvalues.
//
  it_max = 100;

  v = new double[n*n];

  jacobi_eigenvalue ( n, a, it_max, v, x, it_num, rot_num );

  for ( j = 0; j < n; j++ )
  {
    w[j] = beta[0] * v[0+j*n] * v[0+j*n];
  }

  delete [] a;
  delete [] v;

  return;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
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
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], 
  double d[], int &it_num, int &rot_num )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
//
//  Discussion:
//
//    This function computes the eigenvalues and eigenvectors of a
//    real symmetric matrix, using Rutishauser's modfications of the classical
//    Jacobi rotation method with threshold pivoting. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix, which must be square, real,
//    and symmetric.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, double V[N*N], the matrix of eigenvectors.
//
//    Output, double D[N], the eigenvalues, in descending order.
//
//    Output, int &IT_NUM, the total number of iterations.
//
//    Output, int &ROT_NUM, the total number of rotations.
//
{
  double *bw;
  double c;
  double g;
  double gapq;
  double h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int p;
  int q;
  double s;
  double t;
  double tau;
  double term;
  double termp;
  double termq;
  double theta;
  double thresh;
  double w;
  double *zw;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        v[i+j*n] = 1.0;
      }
      else
      {
        v[i+j*n] = 0.0;
      }
    }
  }

  for ( i = 0; i < n; i++ )
  {
    d[i] = a[i+i*n];
  }

  bw = new double[n];
  zw = new double[n];

  for ( i = 0; i < n; i++ )
  {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  it_num = 0;
  rot_num = 0;

  while ( it_num < it_max )
  {
    it_num = it_num + 1;
//
//  The convergence threshold is based on the size of the elements in
//  the strict upper triangle of the matrix.
//
    thresh = 0.0;
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < j; i++ )
      {
        thresh = thresh + a[i+j*n] * a[i+j*n];
      }
    }

    thresh = sqrt ( thresh ) / ( double ) ( 4 * n );

    if ( thresh == 0.0 )
    {
      break;
    }

    for ( p = 0; p < n; p++ )
    {
      for ( q = p + 1; q < n; q++ )
      {
        gapq = 10.0 * fabs ( a[p+q*n] );
        termp = gapq + fabs ( d[p] );
        termq = gapq + fabs ( d[q] );
//
//  Annihilate tiny offdiagonal elements.
//
        if ( 4 < it_num &&
             termp == fabs ( d[p] ) &&
             termq == fabs ( d[q] ) )
        {
          a[p+q*n] = 0.0;
        }
//
//  Otherwise, apply a rotation.
//
        else if ( thresh <= fabs ( a[p+q*n] ) )
        {
          h = d[q] - d[p];
          term = fabs ( h ) + gapq;

          if ( term == fabs ( h ) )
          {
            t = a[p+q*n] / h;
          }
          else
          {
            theta = 0.5 * h / a[p+q*n];
            t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
            if ( theta < 0.0 )
            {
              t = - t;
            }
          }
          c = 1.0 / sqrt ( 1.0 + t * t );
          s = t * c;
          tau = s / ( 1.0 + c );
          h = t * a[p+q*n];
//
//  Accumulate corrections to diagonal elements.
//
          zw[p] = zw[p] - h;                 
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;

          a[p+q*n] = 0.0;
//
//  Rotate, using information from the upper triangle of A only.
//
          for ( j = 0; j < p; j++ )
          {
            g = a[j+p*n];
            h = a[j+q*n];
            a[j+p*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = p + 1; j < q; j++ )
          {
            g = a[p+j*n];
            h = a[j+q*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = q + 1; j < n; j++ )
          {
            g = a[p+j*n];
            h = a[q+j*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[q+j*n] = h + s * ( g - h * tau );
          }
//
//  Accumulate information in the eigenvector matrix.
//
          for ( j = 0; j < n; j++ )
          {
            g = v[j+p*n];
            h = v[j+q*n];
            v[j+p*n] = g - s * ( h + g * tau );
            v[j+q*n] = h + s * ( g - h * tau );
          }
          rot_num = rot_num + 1;
        }
      }
    }

    for ( i = 0; i < n; i++ )
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
//
//  Restore upper triangle of input matrix.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[i+j*n] = a[j+i*n];
    }
  }
//
//  Ascending sort the eigenvalues and eigenvectors.
//
  for ( k = 0; k < n - 1; k++ )
  {
    m = k;
    for ( l = k + 1; l < n; l++ )
    {
      if ( d[l] < d[m] )
      {
        m = l;
      }
    }

    if ( m != k )
    {
      t    = d[m];
      d[m] = d[k];
      d[k] = t;
      for ( i = 0; i < n; i++ )
      {
        w        = v[i+m*n];
        v[i+m*n] = v[i+k*n];
        v[i+k*n] = w;
      }
    }
  }

  delete [] bw;
  delete [] zw;

  return;
}
//****************************************************************************80

void r_jacobi ( int n, double a, double b, double alpha[], double beta[] )

//****************************************************************************80
//
//  Purpose:
//
//    R_JACOBI computes recurrence coefficients for monic Jacobi polynomials.
//
//  Discussion:
//
//    This function generates the first N recurrence coefficients for monic 
//    Jacobi polynomials with parameters A and B. 
//
//    These polynomials are orthogonal on [-1,1] relative to the weight
//
//      w(x) = (1.0-x)^A * (1.0+x)^B. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2013
//
//  Author:
//
//    Original MATLAB version by Dirk Laurie, Walter Gautschi.
//    C version by John Burkardt.
//
//  Reference:
//
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//
//  Parameters:
//
//    Input, int N, the number of coefficients desired.
//
//    Input, double A, B, the parameters for the Jacobi polynomial.
//    -1.0 < A, -1.0 < B.
//
//    Output, double ALPHA[N], BETA[N], the first N recurrence
//    coefficients.
//
{
  int i;
  double i_r8;
  double mu;
  double nab;
  double nu;

  if ( a <= -1.0 )
  {
    cerr << "\n";
    cerr << "R_JACOBI - Fatal error!\n";
    cerr << "  Illegal value of A.\n";
    exit ( 1 );
  }

  if ( b <= -1.0 )
  {
    cerr << "\n";
    cerr << "R_JACOBI - Fatal error!\n";
    cerr << "  Illegal value of B.\n";
    exit ( 1 );
  }

  nu = ( b - a ) / ( a + b + 2.0 );

  mu = pow ( 2.0, a + b + 1.0 )
    * r8_gamma ( a + 1.0 ) 
    * r8_gamma ( b + 1.0 ) 
    / r8_gamma ( a + b + 2.0 );

  alpha[0] = nu;
  beta[0] = mu;

  if ( n == 1 )
  {
    return;
  }

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    alpha[i] = ( b - a ) * ( b + a ) 
      / ( 2.0 * ( i_r8 - 1.0 ) + a + b ) 
      / ( 2.0 * i_r8 + a + b );
  }

  beta[1] = 4.0 * ( a + 1.0 ) * ( b + 1.0 ) 
    / ( a + b + 2.0 ) / ( a + b + 2.0 )
    / ( a + b + 3.0 );

  for ( i = 2; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    nab = 2.0 * ( i_r8 - 1.0 ) + a +  b;
    beta[i] = 4.0 * ( i_r8 - 1.0 + a ) * ( i_r8 - 1.0 + b ) 
      * ( i_r8 - 1.0 ) * ( i_r8 - 1.0 + a + b ) 
      / nab / nab
      / ( nab + 1.0 ) 
      / ( nab - 1.0 );
  }

  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
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
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_acos ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOS computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double R8_ACOS, an angle whose cosine is C.
//
{
  const double pi = 3.141592653589793;
  double value;

  if ( c <= -1.0 )
  {
    value = pi;
  }
  else if ( 1.0 <= c )
  {
    value = 0.0;
  }
  else
  {
    value = acos ( c );
  }
  return value;
}
//****************************************************************************80

double r8_asin ( double s )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ASIN computes the arc sine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ASIN routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the argument, the sine of an angle.
//
//    Output, double R8_ASIN, an angle whose sine is S.
//
{
  double angle;
  const double pi = 3.141592653589793;

  if ( s <= -1.0 )
  {
    angle = - pi / 2.0;
  }
  else if ( 1.0 <= s )
  {
    angle = pi / 2.0;
  }
  else
  {
    angle = asin ( s );
  }
  return angle;
}
//****************************************************************************80

double r8_atan ( double y, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ATAN computes the inverse tangent of the ratio Y / X.
//
//  Discussion:
//
//    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
//    the built in functions ATAN and ATAN2 already do.
//
//    However:
//
//    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
//      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
//      and [-PI,+PI] respectively;
//
//    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
//     function by contrast always returns an angle in the first or fourth
//     quadrants.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double Y, X, two quantities which represent the tangent of
//    an angle.  If Y is not zero, then the tangent is (Y/X).
//
//    Output, double R8_ATAN, an angle between 0 and 2 * PI, whose tangent is
//    (Y/X), and which lies in the appropriate quadrant so that the signs
//    of its cosine and sine match those of X and Y.
//
{
  double abs_x;
  double abs_y;
  const double pi = 3.141592653589793;
  double theta;
  double theta_0;
//
//  Special cases:
//
  if ( x == 0.0 )
  {
    if ( 0.0 < y )
    {
      theta = pi / 2.0;
    }
    else if ( y < 0.0 )
    {
      theta = 3.0 * pi / 2.0;
    }
    else if ( y == 0.0 )
    {
      theta = 0.0;
    }
  }
  else if ( y == 0.0 )
  {
    if ( 0.0 < x )
    {
      theta = 0.0;
    }
    else if ( x < 0.0 )
    {
      theta = pi;
    }
  }
//
//  We assume that ATAN2 is correct when both arguments are positive.
//
  else
  {
    abs_y = r8_abs ( y );
    abs_x = r8_abs ( x );

    theta_0 = atan2 ( abs_y, abs_x );

    if ( 0.0 < x && 0.0 < y )
    {
      theta = theta_0;
    }
    else if ( x < 0.0 && 0.0 < y )
    {
      theta = pi - theta_0;
    }
    else if ( x < 0.0 && y < 0.0 )
    {
      theta = pi + theta_0;
    }
    else if ( 0.0 < x && y < 0.0 )
    {
      theta = 2.0 * pi - theta_0;
    }
  }

  return theta;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for an R8.
//
//  Discussion:
//
//    The C MATH library includes a function GAMMA ( X ) which should be
//    invoked instead of this function.
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  const double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = false;
  fact = 1.0;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

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
  int i4_huge = 2147483647;
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

double *r8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
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
//    03 October 2005
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
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
//
//    In other words, the interval is divided into N-1 even subintervals,
//    and the endpoints of intervals are used as the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
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
  int i4_huge = 2147483647;
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

double *tridisolve ( int n, double a[], double b[], double c[], double d[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIDISOLVE solves a tridiagonal system of linear equations.
//
//  Discussion:
//
//    We can describe an NxN tridiagonal matrix by vectors A, B, and C, where
//    A and C are of length N-1.  In that case, a linear system can be
//    represented as
//                        b(1) * x(1) + c(1) * x(2)   = d(1),
//      a(j-1) * x(j-1) + b(j) * x(j) + c(j) * x(j+1) = d(j), j = 2:n-1,
//      a(n-1) * x(n-1) + b(n) * x(n)                 = d(n)
//
//    This function produces the solution vector X.
//
//    This function is derived from Cleve Moler's Matlab suite.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 July 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A(N-1), B(N), C(N-1), the matrix entries.
//
//    Input, double D(N), the right hand side.
//
//    Output, double TRIDISOLVE[N], the solution.
//
{
  double *bi;
  int j;
  double mu;
  double *x;

  x = new double[n];

  for ( j = 0; j < n; j++ )
  {
    x[j] = d[j];
  }

  bi = new double[n];
  for ( j = 0; j < n; j++ )
  {
    bi[j] = 1.0 / b[j];
  }

  for ( j = 0; j < n - 1; j++ )
  {
    mu = a[j] * bi[j];
    b[j+1] = b[j+1] - mu * c[j];
    x[j+1] = x[j+1] - mu * x[j];
  }

  x[n-1] = x[n-1] * bi[n-1];
  for ( j = n - 2; 0 <= j; j-- )
  {
    x[j] = ( x[j] - c[j] * x[j+1] ) * bi[j];
  }

  delete [] bi;

  return x;
}
