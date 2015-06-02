# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "circle_segment.hpp"

int main ( );
void test01 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test11 ( );
void test13 ( );
void test14 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CIRCLE_SEGMENT_PRB.
//
//  Discussion:
//
//    CIRCLE_SEGMENT_PRB tests the CIRCLE_SEGMENT library.
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
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CIRCLE_SEGMENT_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CIRCLE_SEGMENT library.\n";

  test01 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test11 ( );
  test13 ( );
  test14 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CIRCLE_SEGMENT_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
///
//  Purpose:
//
//    TEST01 tests CIRCLE_SEGMENT_AREA_FROM_HEIGHT.
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
{
  double area;
  double h;
  int i;
  double r;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.\n";

  cout << "\n";
  cout << "          R               H               Area\n";
  cout << "\n";
  r = 1.0;
  h = 1.0;
  for ( i = 0; i <= 10; i++ )
  {
    area = circle_segment_area_from_height ( r, h );
    cout << "  " << setw(14) << r
         << "  " << setw(14) << h
         << "  " << setw(14) << area << "\n";
    h = h / 2.0;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests the AREA and HEIGHT functions.
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
{
  double a;
  double a2;
  double h;
  double h2;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  int test;

  cout << "\n";
  cout << "CIRCLE_SEGMENT_TEST05\n";
  cout << "  For circle segment with a given radius R,\n";
  cout << "  CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area A, given the height.\n";
  cout << "  CIRCLE_SEGMENT_HEIGHT_FROM_AREA computes height H, given the area.\n";
  cout << "  Check that these functions are inverses of each other\n";
  cout << "  using random values of R, A, and H.\n";

  cout << "\n";
  cout << "        R             H      =>     A    =>       H2\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( seed );
    h = 2.0 * r * r8_uniform_01 ( seed );
    a = circle_segment_area_from_height ( r, h );
    h2 = circle_segment_height_from_area ( r, a );
    cout << "  " << setw(12) << r
         << "  " << setw(12) << h
         << "  " << setw(12) << a
         << "  " << setw(12) << h2 << "\n";
  }

  cout << "\n";
  cout << "        R             A      =>     H    =>       A2\n";
  cout << "\n";

  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( seed );
    a = pi * r * r * r8_uniform_01 ( seed );
    h = circle_segment_height_from_area ( r, a );
    a2 = circle_segment_area_from_height ( r, h );
    cout << "  " << setw(12) << r
         << "  " << setw(12) << a
         << "  " << setw(12) << h
         << "  " << setw(12) << a2 << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 samples using CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT.
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
{
  double *an;
  int an_num = 51;
  string boundary_filename = "sample00_boundary.txt";
  ofstream boundary_unit;
  double *boundary_x;
  double *boundary_y;
  string command_filename = "sample00_commands.txt";
  ofstream command_unit;
  string data_filename = "sample00_data.txt";
  ofstream data_unit;
  int data_num = 100;
  double *data_x;
  double *data_y;
  string graphics_filename = "sample00.png";
  double h;
  int i;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  int test;
  double theta;
  double thetah;

  seed = 123456789;

  cout << "\n";
  cout << "CIRCLE_SEGMENT_TEST06\n";
  cout << "  CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples a circle segment.\n";
  cout << "\n";
  cout << "  Plot " << data_num << " points from several segments.\n";
  cout << "\n";

  r = 1.0;
  theta = pi;

  for ( test = 1; test <= 4; test++ )
  {
    h = circle_segment_height_from_angle ( r, theta );

    thetah = theta / 2.0;
//
//  Create boundary.
//
    an = r8vec_linspace_new ( an_num, -thetah, +thetah );
    for ( i = 0; i < an_num; i++ )
    {
      an[i] = an[i] + 0.5 * pi;
    }

    boundary_x = new double[an_num+1];
    boundary_y = new double[an_num+1];

    for ( i = 0; i < an_num; i++ )
    {
      boundary_x[i] = r * cos ( an[i] );
      boundary_y[i] = r * sin ( an[i] );
    }
    boundary_x[an_num] = boundary_x[0];
    boundary_y[an_num] = boundary_y[0];

    filename_inc ( &boundary_filename );
    boundary_unit.open ( boundary_filename.c_str ( ) );
    for ( i = 0; i <= an_num; i++ )
    {
      boundary_unit << "  " << setw(14) << boundary_x[i]
                    << "  " << setw(14) << boundary_y[i] << "\n";
    }
    boundary_unit.close ( );
    cout << "\n";
    cout << "  Created boundary file \"" << boundary_filename << "\".\n";
//
//  Create data.
//
    data_x = new double[data_num+1];
    data_y = new double[data_num+1];

    circle_segment_sample_from_height ( r, h, data_num, seed, data_x, data_y );

    filename_inc ( &data_filename );
    data_unit.open ( data_filename.c_str ( ) );
    for ( i = 0; i < data_num; i++ )
    {
      data_unit << "  " << setw(14) << data_x[i]
                << "  " << setw(14) << data_y[i] << "\n";
    }
    data_unit.close ( );
    cout << "\n";
    cout << "  Created data file \"" << data_filename << "\".\n";
//
//  Create commands.
//
    filename_inc ( &command_filename );
    command_unit.open ( command_filename.c_str ( ) );
    command_unit << "# " << command_filename << "\n";
    command_unit << "#\n";
    command_unit << "# Usage:\n";
    command_unit << "#  gnuplot < " << command_filename << "\n";
    command_unit << "#\n";
    command_unit << "set term png\n";
    filename_inc ( &graphics_filename );
    command_unit << "set output '" << graphics_filename << "'\n";
    command_unit << "set xlabel '<--- X --->'\n";
    command_unit << "set ylabel '<--- Y --->'\n";
    command_unit << "set title 'Circle Segment Sample'\n";
    command_unit << "set grid\n";
    command_unit << "set key off\n";
    command_unit << "set size ratio -1\n";
    command_unit << "set style data lines\n";
    command_unit << "plot '" << data_filename 
                 << "' using 1:2 with points lt 3 pt 3,\\\n";
    command_unit << "    '" << boundary_filename 
                 << "' using 1:2 lw 3 linecolor rgb 'black'\n";
    command_unit << "quit\n";
    command_unit.close ( );

    cout << "  Created command file \"" << command_filename << "\".\n";

    theta = theta / 2.0;

    delete [] an;
    delete [] boundary_x;
    delete [] boundary_y;
    delete [] data_x;
    delete [] data_y;
  }
 
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests the ANGLE and HEIGHT functions.
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
{
  double h;
  double h2;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  double t;
  double t2;
  int test;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  For circle segment with a given radius R,\n";
  cout << "  CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT computes the angle THETA, given the height.\n";
  cout << "  CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE computes height H, given the angle.\n";
  cout << "  Check that these functions are inverses of each other\n";
  cout << "  using random values of R, T, and H.\n";
  cout << "\n";
  cout << "        R             H      =>     T    =>       H2\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( seed );
    h = 2.0 * r * r8_uniform_01 ( seed );
    t = circle_segment_angle_from_height ( r, h );
    h2 = circle_segment_height_from_angle ( r, t );
    cout << "  " << setw(12) << r
         << "  " << setw(12) << h
         << "  " << setw(12) << t
         << "  " << setw(12) << h2 << "\n";
  }

  cout << "\n";
  cout << "        R             T      =>     H    =>       T2\n";
  cout << "\n";
  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( seed );
    t = 2.0 * pi * r8_uniform_01 ( seed );
    h = circle_segment_height_from_angle ( r, t );
    t2 = circle_segment_angle_from_height ( r, h );
    cout << "  " << setw(12) << r
         << "  " << setw(12) << t
         << "  " << setw(12) << h
         << "  " << setw(12) << t2 << "\n";
  }

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests CIRCLE_SEGMENT_CONTAINS_POINT.
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
{
  double area;
  double area_est;
  double c[2];
  int i;
  int *inout;
  int j;
  int n = 1000;
  double omega1;
  double omega2;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  int test;
  double theta;
  double *xy;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  CIRCLE_SEGMENT_CONTAINS_POINT reports whether\n";
  cout << "  a circle segment contains a point.\n";
  cout << "\n";
  cout << "  Pick a circle segment at random.\n";
  cout << "  Compute " << n << " sample points in the surrounding box.\n";
  cout << "  Compare the area of the segment to the percentage of points\n";
  cout << "  contained in the circle segment.\n";
  cout << "\n";
  cout << "       N       Omega1          Omega2           Area         Estimate\n";
  cout << "\n";

  r = 1.0;
  c[0] = 0.0;
  c[1] = 0.0;
  seed = 123456789;
  inout = new int[n];

  for ( test = 1; test <= 5; test++ )
  {
    omega1 = 2.0 * pi * r8_uniform_01 ( seed );
    omega2 = 2.0 * pi * r8_uniform_01 ( seed );
  
    if ( omega2 < omega1 )
    {
      omega2 = omega2 + 2.0 * pi;
    }

    xy = r8mat_uniform_01_new ( 2, n, seed );
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        xy[i+j*2] = 2.0 * xy[i+j*2] - 1.0;
      }
    }

    for ( j = 0; j < n; j++ )
    {
      inout[j] = circle_segment_contains_point ( r, c, omega1, omega2, xy + j * 2 );
    }

    theta = circle_segment_angle_from_chord_angles ( omega1, omega2 );
    area = circle_segment_area_from_angle ( r, theta );
    area_est = 4.0 * ( double ) ( i4vec_sum ( n, inout ) ) / ( double ) ( n );

    cout << "  " << setw(6) << n
         << "  " << setw(14) << omega1
         << "  " << setw(14) << omega2
         << "  " << setw(14) << area 
         << "  " << setw(14) << area_est << "\n";

    delete [] xy;
  }

  delete [] inout;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SEGMENT_TEST09 looks at the area and centroid calculations.
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
{
  double a1;
  double a2;
  double a3;
  double c[2];
  double *d1;
  double *d2;
  double *d3;
  double h;
  int n;
  double omega1;
  double omega2;
  double p1[2];
  double p2[2];
  const double pi = 3.141592653589793;
  double r;
  int seed;
  double theta;

  cout << "\n";
  cout << "CIRCLE_SEGMENT_TEST09\n";
  cout << "  CIRCLE_SEGMENT_AREA_FROM_CHORD and\n";
  cout << "  CIRCLE_SEGMENT_CENTROID_FROM_CHORD evaluate the area\n";
  cout << "  and centroid of a circle segment, given R, C and P1:P2.\n";
  cout << "\n";
  cout << "  CIRCLE_SEGMENT_AREA_FROM_SAMPLE and\n";
  cout << "  CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE give us Monte Carlo estimates.\n";
  cout << "\n";
  cout << "  GQCIRCSEGM can estimate these values by quadrature.\n";
  cout << "\n";
  cout << "  Start easy, with R = 1, C = (0,0), and Theta centered.\n";

  seed = 123456789;
  r = 1.0;
  c[0] = 0.0;
  c[1] = 0.0;
  theta = pi / 4.0;
  h = circle_segment_height_from_angle ( r, theta );
  omega1 = - theta / 2.0;
  omega2 = + theta / 2.0;
  p1[0] = c[0] + r * cos ( omega1 );
  p1[1] = c[1] + r * sin ( omega1 );
  p2[0] = c[0] + r * cos ( omega2 );
  p2[1] = c[1] + r * sin ( omega2 );

  a1 = circle_segment_area_from_chord ( r, c, p1, p2 );
  d1 = circle_segment_centroid_from_chord ( r, c, p1, p2 );

  cout << "\n";
  cout << "         Area          CentroidX    CentroidY\n";
  cout << "\n";
  cout << "  " << setw(14) << a1
       << "  " << setw(14) << d1[0]
       << "  " << setw(14) << d1[1] << "\n";
//
//  This only works because the centroid of the height-based circle segment 
//  is easily transformed to the centroid of the chord based circle segment.
//
  a2 = circle_segment_area_from_height ( r, h );
  d2 = circle_segment_centroid_from_height ( r, h );
  cout << "  " << setw(14) << a2
       << "  " << setw(14) << d2[1]
       << "  " << setw(14) << -d2[0] << "\n";

  n = 10000;
  a3 = circle_segment_area_from_sample ( r, c, p1, p2, n, seed );
  d3 = circle_segment_centroid_from_sample ( r, c, p1, p2, n, seed );
  cout << "  " << setw(14) << a3
       << "  " << setw(14) << d3[0]
       << "  " << setw(14) << d3[1] << "\n";

  delete [] d1;
  delete [] d2;
  delete [] d3;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 demonstrates CIRCLE_SEGMENT_ROTATION_FROM_CHORD.
//
//  Discussion:
//
//    We make a table of all pairs of angles that are multiples of pi/12.
//
//    For each pair, we compute the rotation, that is, the angle of the
//    central radius of the circle segment.  We print out the result in
//    terms of multiples of pi/12.
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
{
  double alpha;
  double c[2];
  int i;
  int j;
  double p1[2];
  double p2[2];
  double pi = 3.141592653589793;
  double r;
  double rho1;
  double rho2;
  double t;

  cout << "\n";
  cout << "TEST11:\n";
  cout << "  CIRCLE_SEGMENT_ROTATION_FROM_CHORD is given the endpoints\n";
  cout << "  of a chord, and is asked to determine the angle of the\n";
  cout << "  central radius vector.\n";
  cout << "\n";
  cout << "  We make a table of all pairs of angles that are multiples\n";
  cout << "  of pi/12, determine the corresponding chord endpoints, and\n";
  cout << "  compute the rotation angle, also printed as a multiple of pi/12.\n";

  r = 2.0;
  c[0] = 3.0;
  c[1] = 5.0;
  cout << "\n";
  cout << "     0.0   1.0   2.0   3.0   4.0   5.0   6.0   7.0";
  cout << "   8.0   9.0  10.0  11.0  12.0\n";
  cout << "\n";
  for ( i = 0; i <= 12; i++ )
  {
    rho1 = ( double ) ( i ) * pi / 6.0;
    p1[0] = c[0] + r * cos ( rho1 );
    p1[1] = c[1] + r * sin ( rho1 );
    cout << setw(2) << i;
    for ( j = 0; j <= 12; j++ )
    {
      rho2 = ( double ) ( j ) * pi / 6.0;
      p2[0] = c[0] + r * cos ( rho2 );
      p2[1] = c[1] + r * sin ( rho2 );
      alpha = circle_segment_rotation_from_chord ( r, c, p1, p2 );
      t = 6.0 * alpha / pi;
      cout << "  " << setw(4) << t;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 demonstrates GAUSS for quadrature rules.
//
//  Discussion:
//
//    Some recursion coefficients ALPHA and BETA are listed in Kautsky
//    and Elhay.
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
//  Reference
//
//    Jaroslav Kautsky, Sylvan Elhay,
//    Calculation of the Weights of Interpolatory Quadratures,
//    Numerische Mathematik,
//    Volume 40, Number 3, October 1982, pages 407-422.
//
{
  double *alpha;
  double *beta;
  int i;
  int n;
  double pi = 3.141592653589793;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  GAUSS computes the points and weights for a\n";
  cout << "  Gauss quadrature rule, given the ALPHA and BETA\n";
  cout << "  recursion coefficients.\n";
//
//  Legendre rule.
//
  n = 10;

  alpha = new double[n];
  beta = new double[n];

  for ( i = 0; i < n; i++ )
  {
    alpha[i] = 0.0;
    if ( i == 0 )
    {
      beta[i] = 2.0;
    }
    else
    {
      beta[i] = 1.0 / ( 4.0 - 1.0 / ( double ) ( i * i ) );
    }
  }

  x = new double[n];
  w = new double[n];

  gauss ( n, alpha, beta, x, w );

  cout << "\n";
  cout << "  LEGENDRE RULE\n";
  cout << "  Point   Weight\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[i] << "\n";
  }

  delete [] alpha;
  delete [] beta;
  delete [] w;
  delete [] x;
//
//  Hermite rule.
//
  n = 10;

  alpha = new double[n];
  beta = new double[n];

  for ( i = 0; i < n; i++ )
  {
    alpha[i] = 0.0;
    if ( i == 0 )
    {
      beta[i] = sqrt ( pi );
    }
    else
    {
      beta[i] = ( double ) ( i ) / 2.0;
    }
  }

  x = new double[n];
  w = new double[n];

  gauss ( n, alpha, beta, x, w );

  cout << "\n";
  cout << "  HERMITE RULE\n";
  cout << "  Point   Weight\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[i] << "\n";
  }

  delete [] alpha;
  delete [] beta;
  delete [] w;
  delete [] x;
//
//  Laguerre rule.
//
  n = 10;

  alpha = new double[n];
  beta = new double[n];

  for ( i = 0; i < n; i++ )
  {
    alpha[i] = 2.0 * ( double ) ( i + 1 ) - 1.0;
    if ( i == 0 )
    {
      beta[i] = 1.0;
    }
    else
    {
      beta[i] = ( double ) ( i * i );
    }
  }

  x = new double[n];
  w = new double[n];

  gauss ( n, alpha, beta, x, w );

  cout << "\n";
  cout << "  LAGUERRE RULE\n";
  cout << "  Point   Weight\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << x[i]
         << "  " << setw(14) << w[i] << "\n";
  }

  delete [] alpha;
  delete [] beta;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 demonstrates R_JACOBI.
//
//  Discussion:
//
//    R_JACOBI returns recursion coefficients ALPHA and BETA for rules
//    using a Jacobi type weight w(x) = (1-x)^A * (1+x)^B.
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
//    John Burkardt
//
//  Reference
//
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//
{
  double a;
  double *alpha;
  double b;
  double *beta;
  int i;
  int n;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  R_JACOBI computes recursion coefficients ALPHA and BETA\n";
  cout << "  Gauss quadrature rule, given the ALPHA and BETA\n";
  cout << "  recursion coefficients.\n";
//
//  Legendre rule.
//
  n = 10;

  a = 0.0;
  b = 0.0;
  alpha = new double[n];
  beta = new double[n];

  r_jacobi ( n, a, b, alpha, beta );

  cout << "\n";
  cout << "  Legendre weight\n";
  cout << "  A = " << a << ",  B = " << b << "\n";
  cout << "  Alpha          Beta\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << alpha[i]
         << "  " << setw(14) << beta[i] << "\n";
  }
  delete [] alpha;
  delete [] beta;
//
//  Chebyshev Type 1 rule.
//
  n = 10;

  a = -0.5;
  b = -0.5;
  alpha = new double[n];
  beta = new double[n];

  r_jacobi ( n, a, b, alpha, beta );

  cout << "\n";
  cout << "  Chebyshev Type 1 weight\n";
  cout << "  A = " << a << ",  B = " << b << "\n";
  cout << "  Alpha          Beta\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << alpha[i]
         << "  " << setw(14) << beta[i] << "\n";
  }
  delete [] alpha;
  delete [] beta;
//
//  Chebyshev Type 2 rule.
//
  n = 10;

  a = +0.5;
  b = +0.5;
  alpha = new double[n];
  beta = new double[n];

  r_jacobi ( n, a, b, alpha, beta );

  cout << "\n";
  cout << "  Chebyshev Type 2 weight\n";
  cout << "  A = " << a << ",  B = " << b << "\n";
  cout << "  Alpha          Beta\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << alpha[i]
         << "  " << setw(14) << beta[i] << "\n";
  }
  delete [] alpha;
  delete [] beta;
//
//  General Jacobi rule.
//
  n = 10;

  a = +0.5;
  b = +1.5;
  alpha = new double[n];
  beta = new double[n];

  r_jacobi ( n, a, b, alpha, beta );

  cout << "\n";
  cout << "  General Jacobi weight\n";
  cout << "  A = " << a << ",  B = " << b << "\n";
  cout << "  Alpha          Beta\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << alpha[i]
         << "  " << setw(14) << beta[i] << "\n";
  }
  delete [] alpha;
  delete [] beta;

  return;
}
