# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "polygon_properties.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POLYGON_PROPERTIES_PRB.
//
//  Discussion:
//
//    POLYGON_PROPERTIES_PRB tests the POLYGON_PROPERTIES library.
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
{
  timestamp ( );
  cout << "\n";
  cout << "POLYGON_PROPERTIES_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the POLYGON_PROPERTIES library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );

  test11 ( );
  test12 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLYGON_PROPERTIES_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests POLYGON_ANGLES.
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
{
  double *angle;
  int i;
  int n = 6;
  double v[2*6] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 1.0,
    3.0, 0.0,
    3.0, 2.0,
    1.0, 2.0 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For a polygon:\n";
  cout << "  POLYGON_ANGLES computes the angles.\n";

  cout << "\n";
  cout << "  Number of polygonal vertices = " << n << "\n";

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  angle = polygon_angles ( n, v );

  cout << "\n";
  cout << "  Polygonal angles in degrees:\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << setw(8) << i << "  "
         << setw(14) << r8_degrees ( angle[i] ) << "\n";
  }

  delete [] angle;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests POLYGON_AREA.
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
{
  double area;
  double area_exact1 = 2.0;
  double area_exact2 = 6.0;
  int n1 = 4;
  int n2 = 8;
  double v1[2*4] = {
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };
  double v2[2*8] = {
        0.0, 0.0, 
        3.0, 0.0, 
        3.0, 3.0, 
        2.0, 3.0, 
        2.0, 1.0, 
        1.0, 1.0, 
        1.0, 2.0, 
        0.0, 2.0 };

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a polygon:\n";
  cout << "  POLYGON_AREA computes the area.\n";

  cout << "\n";
  cout << "  Number of polygonal vertices = " << n1 << "\n";
  r8mat_transpose_print ( 2, n1, v1, "  The polygon vertices:" );
  area = polygon_area ( n1, v1 );
  cout << "\n";
  cout << "  Exact area is        " << area_exact1 << "\n";
  cout << "  The computed area is " << area << "\n";

  cout << "\n";
  cout << "  Number of polygonal vertices = " << n2 << "\n";
  r8mat_transpose_print ( 2, n2, v2, "  The polygon vertices:" );
  area = polygon_area ( n2, v2 );
  cout << "\n";
  cout << "  Exact area is        " << area_exact2 << "\n";
  cout << "  The computed area is " << area << "\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests POLYGON_AREA_2.
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
{
  double area;
  double area_exact1 = 2.0;
  double area_exact2 = 6.0;
  int n1 = 4;
  int n2 = 8;
  double v1[2*4] = {
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };
  double v2[2*8] = {
        0.0, 0.0, 
        3.0, 0.0, 
        3.0, 3.0, 
        2.0, 3.0, 
        2.0, 1.0, 
        1.0, 1.0, 
        1.0, 2.0, 
        0.0, 2.0 };

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a polygon:\n";
  cout << "  POLYGON_AREA_2 computes the area.\n";

  cout << "\n";
  cout << "  Number of polygonal vertices = " << n1 << "\n";
  r8mat_transpose_print ( 2, n1, v1, "  The polygon vertices:" );
  area = polygon_area_2 ( n1, v1 );
  cout << "\n";
  cout << "  Exact area is        " << area_exact1 << "\n";
  cout << "  The computed area is " << area << "\n";

  cout << "\n";
  cout << "  Number of polygonal vertices = " << n2 << "\n";
  r8mat_transpose_print ( 2, n2, v2, "  The polygon vertices:" );
  area = polygon_area_2 ( n2, v2 );
  cout << "\n";
  cout << "  Exact area is        " << area_exact2 << "\n";
  cout << "  The computed area is " << area << "\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests POLYGON_CENTROID and POLYGON_CENTROID_2.
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
{
  double *centroid;
  int n = 4;
  double v[2*4] = {
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For a polygon:\n";
  cout << "  POLYGON_CENTROID computes the centroid.\n";
  cout << "  POLYGON_CENTROID_2 computes the centroid.\n";

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  centroid = polygon_centroid ( n, v );
  r8vec_print ( 2, centroid, "  POLYGON_CENTROID:" );
  delete [] centroid;

  centroid = polygon_centroid_2 ( n, v );
  r8vec_print ( 2, centroid, "  POLYGON_CENTROID_2:" );
  delete [] centroid;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests POLYGON_CONTAINS_POINT and POLYGON_CONTAINS_POINT_2.
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
{
  int i;
  int inside1;
  int inside2;
  int n = 5;
  double p[2];
  double p_test[2*4] = {
    1.0,  1.0, 
    3.0,  4.0, 
    0.0,  2.0, 
    0.5, -0.25 };
  int test;
  int test_num = 4;
  double v[2*5] = {
    0.0, 0.0, 
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 2.0 };
 
  cout << "\n";
  cout << "TEST05\n";
  cout << "  POLYGON_CONTAINS_POINT determines if\n";
  cout << "  a point is in a polygon.\n";
  cout << "  POLYGON_CONTAINS_POINT_2 determines if\n";
  cout << "  a point is in a polygon.\n";

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  cout << "\n";
  cout << "          P          In1  In2\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    p[0] = p_test[0+test*2];
    p[1] = p_test[1+test*2];
 
    inside1 = polygon_contains_point ( n, v, p );

    inside2 = polygon_contains_point_2 ( n, v, p );

    cout << setw(14) << p[0] << "  "
         << setw(14) << p[1] << "  "
         << setw(1) << inside1 << "  "
         << setw(1) << inside2 << "\n";
  } 

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests POLYGON_DIAMETER;
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
{
  double diameter;
  double diameter_exact = 2.0;
  int n = 4;
  double v[2*4] = {
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a polygon:\n";
  cout << "  POLYGON_DIAMETER computes the diameter;\n";

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  diameter = polygon_diameter ( n, v );

  cout << "\n";
  cout << "  Diameter ( computed ) " << diameter << "\n";
  cout << "  Diameter ( exact )    " << diameter_exact << "\n";
 
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests POLYGON_EXPAND;
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
{
  double h;
  int n = 4;
  double v[2*4] = {
    1.0, 1.0, 
    5.0, 1.0, 
    2.0, 4.0, 
    1.0, 3.0 };
  double *w;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  For a polygon:\n";
  cout << "  POLYGON_EXPAND expands it by an amount H.\n";

  h = 0.5;

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  cout << "\n";
  cout << "  The expansion amount H = " << h << "\n";

  w = polygon_expand ( n, v, h );

  r8mat_transpose_print ( 2, n, w, "  The expanded polygon:" );

  delete [] w;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests POLYGON_INRAD_DATA, POLYGON_OUTRAD_DATA, POLYGON_SIDE_DATA;
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
{
  double area;
  int n;
  double radin;
  double radout;
  double side;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For a REGULAR polygon:\n";
  cout << "  the inradius, outradius and side are related.\n";
  cout << "  POLYGON_INRAD_DATA uses the inradius;\n";
  cout << "  POLYGON_OUTRAD_DATA uses the inradius;\n";
  cout << "  POLYGON_SIDE_DATA uses the inradius;\n";

  for ( n = 3; n <= 5; n++ )
  {
    cout << "\n";
    cout << "  Number of polygonal sides = " << n << "\n";
    side = 1.0;
    cout << "\n";
    cout << "  Assuming SIDE = " << side << "\n";
    polygon_side_data ( n, side, area, radin, radout );
    cout << "    AREA =   " << area << "\n";
    cout << "    RADIN =  " << radin << "\n";
    cout << "    RADOUT = " << radout << "\n";
    cout << "\n";
    cout << "  Assuming RADIN = " << radin << "\n";
    polygon_inrad_data ( n, radin, area, radout, side );
    cout << "    AREA =   " << area << "\n";
    cout << "    RADOUT = " << radout << "\n";
    cout << "    SIDE =   " << side << "\n";
    cout << "\n";
    cout << "  Assuming RADOUT = " << radout << "\n";
    polygon_outrad_data ( n, radout, area, radin, side );
    cout << "    AREA =   " << area << "\n";
    cout << "    RADIN =  " << radin << "\n";
    cout << "    SIDE =   " << side << "\n";
  }
  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests POLYGON_INTEGRAL_*.
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
{
  int n = 4;
  double result;
  double v[2*4] = {
    0.0, 0.0, 
    1.0, 0.0, 
    1.0, 1.0, 
    0.0, 1.0 };

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For a polygon:\n";
  cout << "  POLYGON_INTEGRAL_1 integrates 1\n";
  cout << "  POLYGON_INTEGRAL_X integrates X\n";
  cout << "  POLYGON_INTEGRAL_Y integrates Y\n";
  cout << "  POLYGON_INTEGRAL_XX integrates X*X\n";
  cout << "  POLYGON_INTEGRAL_XY integrates X*Y\n";
  cout << "  POLYGON_INTEGRAL_YY integrates Y*Y\n";

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  cout << "\n";
  cout << "  F(X,Y)    Integral\n";
  cout << "\n";

  result = polygon_integral_1 ( n, v );
  cout << "    1    " << result << "\n";

  result = polygon_integral_x ( n, v );
  cout << "    X    " << result << "\n";

  result = polygon_integral_y ( n, v );
  cout << "    Y    " << result << "\n";

  result = polygon_integral_xx ( n, v );
  cout << "  X*X    " << result << "\n";

  result = polygon_integral_xy ( n, v );
  cout << "  X*Y    " << result << "\n";

  result = polygon_integral_yy ( n, v );
  cout << "  Y*Y    " << result << "\n";

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests POLYGON_IS_CONVEX.
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
{ 
  double angle;
  int i;
  int n;
  int n01 = 1;
  int n02 = 2;
  int n03 = 3;
  int n04 = 3;
  int n05 = 3;
  int n06 = 4;
  int n07 = 5;
  int n08 = 5;
  int n09 = 6;
  int n10 = 6;
  int n11 = 8;
  const double r8_pi = 3.141592653589793;
  int result;
  int test;
  int test_num = 11;
  string title;
  double *v;
  double v01[2*1] = {
    0.0, 0.0 };
  double v02[2*2] = {
    0.0, 0.0, 
    1.0, 2.0 };
  double v03[2*3] = {
    0.0, 0.0, 
    2.0, 0.0, 
    1.0, 0.0 };
  double v04[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.0, 2.0 };
  double v05[2*3] = {
    0.0, 0.0, 
    0.0, 2.0, 
    1.0, 0.0 };
  double v06[2*4] = {
    1.0, 0.0, 
    2.0, 0.0, 
    3.0, 1.0, 
    0.0, 1.0 };
  double v07[2*5] = {
    0.0, 0.0, 
    0.5, 0.5, 
    1.0, 0.0, 
    1.0, 1.0, 
    0.0, 1.0 };
  double *v08;
  double *v09;
  double v10[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    1.0, 1.0, 
    0.0, 0.0, 
    2.0, 0.0, 
    1.0, 1.0 };
  double v11[2*8] = { 
    1.0, 0.0, 
    3.0, 0.0, 
    3.0, 3.0, 
    0.0, 3.0, 
    0.0, 1.0, 
    2.0, 1.0, 
    2.0, 2.0, 
    1.0, 2.0 };

  cout << "\n";
  cout << "TEST10\n";
  cout << "  POLYGON_IS_CONVEX determines if a polygon\n";
  cout << "  is convex.\n";

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = n01;
      v = v01;
      title = "  A point:";
    }
    else if ( test == 2 )
    {
      n = n02;
      v = v02;
      title = "  A line:";
    }
    else if ( test == 3 )
    {
      n = n03;
      v = v03;
      title = "  A triangle:";
    }
    else if ( test == 4 )
    {
      n = n04;
      v = v04;
      title = "  A CCW triangle:";
    }
    else if ( test == 5 )
    {
      n = n05;
      v = v05;
      title = "  A CW triangle:";
    }
    else if ( test == 6 )
    {
      n = n06;
      v = v06;
      title = "  Polygon with large angle:";
    }
    else if ( test == 7 )
    {
      n = n07;
      v = v07;
      title = "  Polygon with huge angle:";
    }
    else if ( test == 8 )
    {
      n = n08;
      v08 = new double[2*n];
      for ( i = 0; i < n; i++ )
      {
        angle = ( double ) ( i ) * 4.0 * r8_pi / ( double ) ( n );
        v08[0+i*2] = cos ( angle );
        v08[1+i*2] = sin ( angle );
      }
      v = v08;
      title = "  A five-pointed star:";
    }
    else if ( test == 9 )
    {
      n = n09;
      v09 = new double[2*n];
      for ( i = 0; i < n; i++ )
      {
        angle = ( double ) ( i ) * 2.0 * r8_pi / ( double ) ( n );
        v09[0+i*2] = cos ( angle );
        v09[1+i*2] = sin ( angle );
      }
      v = v09;
      title = "  A hexagon:";
    }
    else if ( test == 10 )
    {
      n = n10;
      v = v10;
      title = "  A triangle twice:";
    }
    else if ( test == 11 )
    {
      n = n11;
      v = v11;
      title = "  Square knot:";
    }

    r8mat_transpose_print ( 2, n, v, title );
    result = polygon_is_convex ( n, v );

    if ( result == -1 )
    {
      cout << "  The polygon is not convex.\n";
    }
    else if ( result == 0 )
    {
      cout << "  The polygon is degenerate and convex.\n";
    }
    else if ( result == 1 )
    {
      cout << "  The polygon is convex and counterclockwise.\n";
    }
    else if ( result == 2 )
    {
      cout << "  The polygon is convex and clockwise.\n";
    }
  }

  delete [] v08;
  delete [] v09;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests POLYGON_LATTICE_AREA.
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
{
  double area;
  int b;
  int i;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  POLYGON_LATTICE_AREA returns the area\n";
  cout << "  of a polygon, measured in lattice points.\n";

  i = 5;
  b = 6;

  cout << "\n";
  cout << "  Number of interior lattice points = " << i << "\n";
  cout << "  Number of boundary lattice points = " << b << "\n";

  area = polygon_lattice_area ( i, b );

  cout << "  Area of polygon is " << area << "\n";

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests POLYGON_SAMPLE.
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
{
  int n = 20;
  int nv = 6;
  int seed;
  double v[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    2.0, 1.0, 
    1.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };
  double *x;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  POLYGON_SAMPLE samples a polygon.\n";

  seed = 123456789;

  x = polygon_sample ( nv, v, n, seed );

  r8mat_transpose_print ( 2, n, x, "  Sample points:" );

  delete [] x;

  return;
}
