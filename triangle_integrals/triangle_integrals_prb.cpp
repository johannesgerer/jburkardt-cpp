# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "triangle_integrals.hpp"

int main ( );
void i4_to_pascal_test ( );
void i4_to_pascal_degree_test ( );
void pascal_to_i4_test ( );
void poly_power_test ( );
void poly_power_linear_test ( );
void poly_print_test ( );
void poly_product_test ( );
void r8mat_print_test ( );
void r8mat_print_some_test ( );
void rs_to_xy_map_test ( );
void triangle_area_test ( );
void triangle_monomial_integral_test ( );
void triangle_poly_integral_test ( );
void triangle_xy_integral_test ( );
void triangle01_monomial_integral_test ( );
void triangle01_poly_integral_test ( );
void trinomial_test ( );
void xy_to_rs_map_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_INTEGRALS_PRB.
//
//  Discussion:
//
//    TRIANGLE_INTEGRALS_PRB tests the TRIANGLE_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_INTEGRALS_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the TRIANGLE_INTEGRALS library.\n";

  i4_to_pascal_test ( );
  i4_to_pascal_degree_test ( );
  pascal_to_i4_test ( );
  r8mat_print_test ( );
  r8mat_print_some_test ( );
  trinomial_test ( );

  rs_to_xy_map_test ( );
  xy_to_rs_map_test ( );

  poly_print_test ( );
  poly_power_linear_test ( );
  poly_power_test ( );
  poly_product_test ( );

  triangle01_monomial_integral_test ( );
  triangle01_poly_integral_test ( );
  triangle_area_test ( );
  triangle_xy_integral_test ( );
  triangle_monomial_integral_test ( );
  triangle_poly_integral_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_INTEGRALS_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void i4_to_pascal_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_PASCAL_TEST tests I4_TO_PASCAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "I4_TO_PASCAL_TEST\n";
  cout << "  I4_TO_PASCAL converts a linear index to\n";
  cout << "  Pascal triangle indices.\n";
  cout << "\n";
  cout << "     K  =>   I     J\n";
  cout << "\n";

  for ( k = 1; k <= 20; k++ )
  {
    i4_to_pascal ( k, i, j );
    cout << "  " << setw(4) << k
         << "    " << setw(4) << i
         << "  " << setw(4) << j << "\n";
  }

  return;
}
//****************************************************************************80

void i4_to_pascal_degree_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_PASCAL_DEGREE_TEST tests I4_TO_PASCAL_DEGREE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int k;

  cout << "\n";
  cout << "I4_TO_PASCAL_DEGREE_TEST\n";
  cout << "  I4_TO_PASCAL_DEGREE converts a linear index to\n";
  cout << "  the degree of the corresponding Pascal triangle indices.\n";
  cout << "\n";
  cout << "     K  =>   D\n";
  cout << "\n";

  for ( k = 1; k <= 20; k++ )
  {
    d = i4_to_pascal_degree ( k );
    cout << "  " << setw(4) << k
         << "    " << setw(4) << d << "\n";
  }

  return;
}
//****************************************************************************80

void pascal_to_i4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PASCAL_TO_I4_TEST tests PASCAL_TO_I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "PASCAL_TO_I4_TEST\n";
  cout << "  PASCAL_TO_I4 converts Pascal triangle indices to a\n";
  cout << "  linear index.\n";
  cout << "\n";
  cout << "     I     J =>    K\n";
  cout << "\n";

  for ( d = 0; d <= 4; d++ )
  {
    for ( i = d; 0 <= i; i-- )
    {
      j = d - i;
      k = pascal_to_i4 ( i, j );
      cout << "  " << setw(4) << i
           << "  " << setw(4) << j
           << "    "  << setw(4) << k << "\n";
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void poly_power_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_POWER_TEST tests POLY_POWER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n1 = 2;
  int n4 = 3;

  int d1 = 1;
  int d2;
  int d3 = 2;
  int d4 = 2;
  int d5;
  int d6 = 6;

  double p1[3] = { 1.0, 2.0, 3.0 };
  double *p2;
  double p3[6] = { 1.0, 4.0, 6.0, 4.0, 12.0, 9.0 };
  double p4[6] = { 1.0, -2.0, 3.0, -4.0, +5.0, -6.0 };
  double *p5;
  double p6[28] = {
      1.0, 
     -6.0,  9.0,  
      0.0, -21.0,    9.0, 
     40.0, -96.0,  108.0,  -81.0, 
      0.0,  84.0, -141.0,  171.0,  -54.0, 
    -96.0, 384.0, -798.0, 1017.0, -756.0, 324.0, 
    -64.0, 240.0, -588.0,  845.0, -882.0, 540.0, -216.0 };

  cout << "\n";
  cout << "POLY_POWER_TEST:\n";
  cout << "  POLY_POWER computes the N-th power of an X,Y polynomial.\n";
//
//  P1 = ( 1 + 2 x + 3 y )
//  P2 = P1^2 = 1 + 4x + 6y + 4x^2 + 12xy + 9y^2 
//  P3 = correct value
//
  cout << "\n";
  poly_print ( d1, p1, "  p1(x,y)" );

  d2 = n1 * d1;
  p2 = poly_power ( d1, p1, n1 );
  cout << "\n";
  poly_print ( d2, p2, "  p2(x,y) = p1(x,y)^2" );

  cout << "\n";
  poly_print ( d3, p3, "  p3(x,y)=correct answer" );
//
//  P4 = ( 1 - 2 x + 3 y - 4 x^2 + 5 xy - 6 y^2 )
//  P5 = P4^3 =
//    1
//    -6x +9y
//    +0x^2 - 21xy + 9y^2
//    +40x^3 - 96x^2y  + 108x^y2 - 81y^3
//    +0x^4 + 84x^3y - 141 x^2y^2 +171xy^3 - 54y^4
//    -96x^5 + 384x^4y -798x^3y^2 + 1017 x^2y^3 - 756 xy^4 + 324 y^5
//    -64x^6 + 240x^5y - 588x^4y^2 + 845 x^3y^3 - 882 x^2y^4 +540 xy^5 - 216y^6
//
  cout << "\n";
  poly_print ( d4, p4, "  p4(x,y)" );

  d5 = n4 * d4;
  p5 = poly_power ( d4, p4, n4 );
  cout << "\n";
  poly_print ( d5, p5, "  p5(x,y) = p1(x,y)^3" );

  cout << "\n";
  poly_print ( d6, p6, "  p6(x,y)=correct answer" );

  delete [] p2;
  delete [] p5;

  return;
}
//****************************************************************************80

void poly_power_linear_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_POWER_LINEAR_TEST tests POLY_POWER_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n1 = 2;
  int n4 = 3;

  int d1 = 1;
  int d2 = d1 * n1;
  int d3 = 2;
  int d4 = 1;
  int d5 = d4 * n4;
  int d6 = 3;

  int m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
  int m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
  int m3 = ( ( d3 + 1 ) * ( d3 + 2 ) ) / 2;
  int m4 = ( ( d4 + 1 ) * ( d4 + 2 ) ) / 2;
  int m5 = ( ( d5 + 1 ) * ( d5 + 2 ) ) / 2;
  int m6 = ( ( d6 + 1 ) * ( d6 + 2 ) ) / 2;

  double p1[3] = { 1.0, 2.0, 3.0 };
  double *p2;
  double p3[6] = { 1.0, 4.0, 6.0, 4.0, 12.0, 9.0 };
  double p4[] = { 2.0, -1.0, 3.0 };
  double *p5;
  double p6[10] = { 8.0, -12.0, 36.0, 6.0, -36.0, 54.0, -1.0, 9.0, -27.0, 27.0 };

  cout << "\n";
  cout << "POLY_POWER_LINEAR_TEST:\n";
  cout << "  POLY_POWER_LINEAR computes the N-th power of\n";
  cout << "  a linear polynomial in X and Y.\n";
//
//  P1 = ( 1 + 2 x + 3 y )
//  P2 = P1^2
//  P3 = correct value
//
  cout << "\n";
  poly_print ( d1, p1, "  p1(x,y)" );

  p2 = poly_power_linear ( d1, p1, n1 );
  cout << "\n";
  poly_print ( d2, p2, "  p2(x,y) = p1(x,y)^n" );

  cout << "\n";
  poly_print ( d3, p3, "  Correct answer" );

  delete [] p2;
//
//  P4 = ( 2 - x + 3 y )
//  P5 = P4^3
//  P6 = correct value
//
  cout << "\n";
  poly_print ( d4, p4, "  p4(x,y)" );

  p5 = poly_power_linear ( d4, p4, n4 );
  cout << "\n";
  poly_print ( d5, p5, "  p5(x,y) = p4(x,y)^3" );

  cout << "\n";
  poly_print ( d6, p6, "  Correct answer" );

  delete [] p5;

  return;
}
//****************************************************************************80

void poly_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_PRINT_TEST tests POLY_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d1 = 0;
  int d2 = 1;
  int d3 = 2;
  int d4 = 3;

  int m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
  int m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
  int m3 = ( ( d3 + 1 ) * ( d3 + 2 ) ) / 2;
  int m4 = ( ( d4 + 1 ) * ( d4 + 2 ) ) / 2;

  double p1[1] = { 12.34 };
  double p2[3] = { 1.0, 2.0, 3.0 };
  double p3[6] = { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
  double p4[10] = { 1.0, -2.1, +3.2, -4.3, +5.4, 
    -6.5, +7.6, -8.7, +9.8, -10.9 };

  cout << "\n";
  cout << "POLY_PRINT_TEST:\n";
  cout << "  POLY_PRINT can print a D-degree polynomial in X and Y.\n";
//
//  P1 = 12.34
//
  cout << "\n";
  cout << "  P1(x,y) = 12.34\n";
  poly_print ( d1, p1, "  p1(x,y)" );
//
//  P2 = 1.0 + 2.0 * x + 3.0 * Y
//
  cout << "\n";
  cout << "  P2(x,y) = 1 + 2 * x + 3 * Y\n";
  poly_print ( d2, p2, "  p2(x,y)" );
//
//  P3 = XY
//
  cout << "\n";
  cout << "  P3(x,y) = xy\n";
  poly_print ( d3, p3, "  p3(x,y) = xy" );
//
//  P4 = 1 - 2.1 * x + 3.2 * y - 4.3 * x^2 + 5.4 * xy - 6.5 * y^2
//    + 7.6 * x^3 - 8.7 * x^2y + 9.8 * xy^2 - 10.9 * y^3.
//
  cout << "\n";
  cout << "  P4(x,y) = 1.0 - 2.1 * x + 3.2 * y - 4.3 * x^2 \n";
  cout << "          + 5.4 * xy - 6.5 * y^2 + 7.6 * x^3 \n";
  cout << "          - 8.7 * x^2y + 9.8 * xy^2 - 10.9 * y^3.\n";
  poly_print ( d4, p4, "  p4(x,y)" );

  return;
}
//****************************************************************************80

void poly_product_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_PRODUCT_TEST tests POLY_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d1 = 1;
  int d2 = 1;
  int d3;
  int d4 = d1 + d2;

  int d5 = 2;
  int d6 = 2;
  int d7;
  int d8 = d5 + d6;

  double p1[3] = { 1.0, 2.0, 3.0 };
  double p2[3] = { 4.0, 5.0, 0.0 };
  double *p3;
  double p4[6] = { 4.0, 13.0, 12.0, 10.0, 15.0, 0.0 };
  double p5[6] = { 1.0, -2.0, 3.0, -4.0, +5.0, -6.0 };
  double p6[6] = { 7.0, 0.0, 0.0, 3.0, 0.0, 0.0 };
  double *p7;
  double p8[15] = {
    7.0, 
  -14.0,  21.0, 
  -25.0, +35.0, -42.0, 
   -6.0,   9.0,   0.0, 0.0, 
  -12.0, +15.0, -18.0, 0.0, 0.0 };

  cout << "\n";
  cout << "POLY_PRODUCT_TEST:\n";
  cout << "  POLY_PRODUCT computes the product of two X,Y polynomials.\n";
//
//  P1 = ( 1 + 2 x + 3 y )
//  P2 = ( 4 + 5 x )
//  P3 = P1 * P2
//  P4 = 4 + 13x + 12y + 10x^2 + 15xy + 0y^2 
//
  cout << "\n";
  poly_print ( d1, p1, "  p1(x,y)" );

  cout << "\n";
  poly_print ( d2, p2, "  p2(x,y)" );

  d3 = d1 + d2;
  p3 = poly_product ( d1, p1, d2, p2 );
  cout << "\n";
  poly_print ( d3, p3, "  p3(x,y) = p1(x,y) * p2(x,y)" );

  cout << "\n";
  poly_print ( d4, p4, "  p4(x,y) = correct answer" );
//
//  P5 = ( 1 - 2 x + 3 y - 4x^2 + 5xy - 6y^2)
//  P6 = ( 7 + 3x^2 )
//  P7 = P5 * P6
//  P8 =    7 
//       - 14x   + 21   y 
//       - 25x^2 + 35x  y - 42   y^2 
//       -  6x^3 +  9x^2y +  0x  y^2 + 0  y^3
//       - 12x^4 + 15x^3y - 18x^2y^2 + 0 xy^3 + 0y^4
//
  cout << "\n";
  poly_print ( d5, p5, "  p5(x,y)" );

  cout << "\n";
  poly_print ( d6, p6, "  p6(x,y)" );

  d7 = d5 + d6;
  p7 = poly_product ( d5, p5, d6, p6 );
  cout << "\n";
  poly_print ( d7, p7, "  p7(x,y) = p5(x,y) * p6(x,y)" );

  cout << "\n";
  poly_print ( d8, p8, "  p8(x,y) = Correct answer" );
//
//  Free memory.
//
  delete [] p3;
  delete [] p7;

  return;
}
//****************************************************************************80

void r8mat_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_TEST tests R8MAT_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_PRINT_TEST\n";
  cout << "  R8MAT_PRINT prints an R8MAT.\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print ( m, n, a, "  The R8MAT:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void r8mat_print_some_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "R8MAT_PRINT_SOME_TEST\n";
  cout << "  R8MAT_PRINT_SOME prints some of an R8MAT.\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print_some ( m, n, a, 2, 1, 4, 2, "  The R8MAT, rows 2:4, cols 1:2:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void rs_to_xy_map_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RS_TO_XY_MAP_TEST tests RS_TO_XY_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  int j;
  double t[2*3] = {
    2.0, 0.0, 
    3.0, 4.0, 
    0.0, 3.0 };
  double tr[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.0, 1.0 };
  double x;
  double y;

  cout << "\n";
  cout << "RS_TO_XY_MAP_TEST:\n";
  cout << "  RS_TO_XY_MAP determines the coefficients of\n";
  cout << "  the linear map from a the reference in RS coordinates\n";
  cout << "  to the physical triangle in XY coordinates:\n";
  cout << "    X = a + b * R + c * S\n";
  cout << "    Y = d + e * R + f * S\n";

  r8mat_print ( 2, 3, t, "  XY triangle vertices:" );

  rs_to_xy_map ( t, a, b, c, d, e, f );

  cout << "\n";
  cout << "  Mapping coefficients are:\n";
  cout << "\n";
  cout << "    X = " << a << " + " << b << " * R + " << c << " * S\n";
  cout << "    Y = " << d << " + " << e << " * R + " << f << " * S\n";

  cout << "\n";
  cout << "  Apply map to RS triangle vertices.\n";
  cout << "  Recover XY vertices (2,0), (3,4) and (0,3).\n";
  cout << "\n";
  for ( j = 0; j < 3; j++ )
  {
    x = a + b * tr[0+j*2] + c * tr[1+j*2];
    y = d + e * tr[0+j*2] + f * tr[1+j*2];
    cout << "  V(" << j << ") = (" << x << "," << y << ")\n";
  }

  return;
}
//****************************************************************************80

void triangle_area_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_TEST tests TRIANGLE_AREA_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double angled;
  double angler;
  double area;
  int i;
  double r;
  const double r8_pi = 3.141592653589793;
  double t[2*3] = {
    0.0, 0.0, 
    2.0, 0.0, 
    0.0, 1.0 };

  cout << "\n";
  cout << "TRIANGLE_AREA_TEST:\n";
  cout << "  TRIANGLE_AREA determines the (signed) area of a triangle.\n";

  cout << "\n";
  cout << "  Triangle vertices are:\n";
  cout << "    (X1,Y1) = (0,0)\n";
  cout << "    (X2,Y2) = 2*(cos(angle),sin(angle))\n";
  cout << "    (X3,Y3) = (0,1)\n";
  cout << "  where angle will sweep from 0 to 360 degrees.\n";

  r = 2.0;

  cout << "\n";
  cout << "   I      Angle         X2          Y2          Area\n";
  cout << "        (degrees)\n";
  cout << "\n";
  for ( i = 0; i <= 24; i++ )
  {
    angled = ( double ) ( i ) * 180.0 / 12.0;
    angler = ( double ) ( i ) * r8_pi / 12.0;
    t[0+1*2] = r * cos ( angler );
    t[1+1*2] = r * sin ( angler );
    area = triangle_area ( t );
    cout << "  " << setw(2) << i
         << "  " << setw(10) << angled
         << "  " << setw(10) << t[0+1*2]
         << "  " << setw(10) << t[1+1*2]
         << "  " << setw(10) << area << "\n";
  }

  return;
}
//****************************************************************************80

void triangle_monomial_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_MONOMIAL_INTEGRAL_TEST estimates integrals over a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  double q;
  double q2;
  double t1[2*3] = {
     0.0, 0.0, 
     1.0, 0.0, 
     0.0, 1.0 };
  double t2[2*3] = {
     0.0, 0.0, 
     1.0, 0.0, 
     1.0, 2.0 };
  double t3[2*3] = {
    -3.0, 0.0, 
     6.0, 0.0, 
     0.0, 3.0 };
  double t4[2*3] = {
     0.0, 0.0, 
     4.0, 0.0, 
     0.0, 1.0  };

  cout << "\n";
  cout << "TRIANGLE_MONOMIAL_INTEGRAL_TEST\n";
  cout << "  TRIANGLE_MONOMIAL_INTEGRAL returns the integral Q of\n";
  cout << "  a monomial X^I Y^J over the interior of a triangle.\n";
//
//  Test 1:
//
  i = 1;
  j = 0;

  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t1[0+0*2] << "," << t1[1+0*2] << ")\n";
  cout << "    (" << t1[0+1*2] << "," << t1[1+1*2] << ")\n";
  cout << "    (" << t1[0+2*2] << "," << t1[1+2*2] << ")\n";
  cout << "  Integrand = x^" << i << " * y^" << j << "\n";

  q = triangle_monomial_integral ( i, j, t1 );
  q2 = 1.0 / 6.0;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";
//
//  Test 2:
//
  i = 1;
  j = 1;

  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t2[0+0*2] << "," << t2[1+0*2] << ")\n";
  cout << "    (" << t2[0+1*2] << "," << t2[1+1*2] << ")\n";
  cout << "    (" << t2[0+2*2] << "," << t2[1+2*2] << ")\n";
  cout << "  Integrand = x^" << i << " * y^" << j << "\n";

  q = triangle_monomial_integral ( i, j, t2 );
  q2 = 0.5;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";
//
//  Test 3:
//
  i = 1;
  j = 0;

  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t3[0+0*2] << "," << t3[1+0*2] << ")\n";
  cout << "    (" << t3[0+1*2] << "," << t3[1+1*2] << ")\n";
  cout << "    (" << t3[0+2*2] << "," << t3[1+2*2] << ")\n";
  cout << "  Integrand = x^" << i << " * y^" << j << "\n";

  q = triangle_monomial_integral ( i, j, t3 );
  q2 = 13.5;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";
//
//  Test 4:
//
  i = 1;
  j = 1;

  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t4[0+0*2] << "," << t4[1+0*2] << ")\n";
  cout << "    (" << t4[0+1*2] << "," << t4[1+1*2] << ")\n";
  cout << "    (" << t4[0+2*2] << "," << t4[1+2*2] << ")\n";
  cout << "  Integrand = x^" << i << " * y^" << j << "\n";

  q = triangle_monomial_integral ( i, j, t4 );
  q2 = 2.0 / 3.0;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";

  return;
}
//****************************************************************************80

void triangle_poly_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_POLY_INTEGRAL_TEST estimates integrals over a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d1 = 1;
  int d2 = 2;
  int d3 = 2;
  int d4 = 2;

  int m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
  int m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
  int m3 = ( ( d3 + 1 ) * ( d3 + 2 ) ) / 2;
  int m4 = ( ( d4 + 1 ) * ( d4 + 2 ) ) / 2;

  double p1[3] = { 0.0, 1.0, 0.0 };
  double p2[6] = { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
  double p3[6] = { 2.0, -3.0, 0.0, 0.0, 1.0, 0.0 };
  double p4[6] = { 0.0, 0.0,-40.0, 6.0, 0.0, 0.0 };
  double q;
  double q2;
  double t1[2*3] = {
     0.0, 0.0, 
     1.0, 0.0, 
     0.0, 1.0  };
  double t2[2*3] = {
     0.0, 0.0, 
     1.0, 0.0, 
     1.0, 2.0  };
  double t3[2*3] = {
     0.0, 0.0, 
     1.0, 0.0, 
     1.0, 3.0  };
  double t4[2*3] = {
     0.0, 3.0, 
     1.0, 1.0, 
     5.0, 3.0  };

  cout << "\n";
  cout << "TRIANGLE_POLY_INTEGRAL_TEST\n";
  cout << "  TRIANGLE_POLY_INTEGRAL returns the integral Q of\n";
  cout << "  a polynomial over the interior of a triangle.\n";
//
//  Test 1:
//  Integrate x over reference triangle.
//
  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t1[0+0*2] << "," << t1[1+0*2] << ")\n";
  cout << "    (" << t1[0+1*2] << "," << t1[1+1*2] << ")\n";
  cout << "    (" << t1[0+2*2] << "," << t1[1+2*2] << ")\n";

  poly_print ( d1, p1, "  Integrand p1(x,y)" );

  q = triangle_poly_integral ( d1, p1, t1 );
  q2 = 1.0 / 6.0;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";
//
//  Test 2:
//  Integrate xy over a general triangle.
//
  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t2[0+0*2] << "," << t2[1+0*2] << ")\n";
  cout << "    (" << t2[0+1*2] << "," << t2[1+1*2] << ")\n";
  cout << "    (" << t2[0+2*2] << "," << t2[1+2*2] << ")\n";

  poly_print ( d2, p2, "  Integrand p2(x,y)" );

  q = triangle_poly_integral ( d2, p2, t2 );
  q2 = 0.5;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";
//
//  Test 3:
//  Integrate 2-3x+xy over a general triangle.
//
  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t3[0+0*2] << "," << t3[1+0*2] << ")\n";
  cout << "    (" << t3[0+1*2] << "," << t3[1+1*2] << ")\n";
  cout << "    (" << t3[0+2*2] << "," << t3[1+2*2] << ")\n";

  poly_print ( d3, p3, "  Integrand p3(x,y)" );

  q = triangle_poly_integral ( d3, p3, t3 );
  q2 = 9.0 / 8.0;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";
//
//  Test 4:
//  Integrate -40y + 6x^2 over a general triangle.
//
  cout << "\n";
  cout << "  Triangle vertices:\n";
  cout << "    (" << t4[0+0*2] << "," << t4[1+0*2] << ")\n";
  cout << "    (" << t4[0+1*2] << "," << t4[1+1*2] << ")\n";
  cout << "    (" << t4[0+2*2] << "," << t4[1+2*2] << ")\n";

  poly_print ( d4, p4, "  Integrand p4(x,y)" );

  q = triangle_poly_integral ( d4, p4, t4 );
  q2 = - 935.0 / 3.0;

  cout << "  Computed Q = " << q << "\n";
  cout << "  Exact Q =    " << q2 << "\n";

  return;
}
//****************************************************************************80

void triangle01_monomial_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE01_MONOMIAL_INTEGRAL_TEST estimates integrals over the unit triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int j;
  double q;
 
  cout << "\n";
  cout << "TRIANGLE01_MONOMIAL_INTEGRAL_TEST\n";
  cout << "  TRIANGLE01_MONOMIAL_INTEGRAL returns the integral Q of\n";
  cout << "  a monomial X^I Y^J over the interior of the unit triangle.\n";

  cout << "\n";
  cout << "   I   J         Q(I,J)\n";

  for ( d = 0; d <= 5; d++ )
  {
    cout << "\n";
    for ( i = 0; i <= d; i++ )
    {
      j = d - i;
      q = triangle01_monomial_integral ( i, j );
      cout << "  " << setw(2) << i
           << "  " << setw(2) << j
           << "  " << q << "\n";
    }
  }

  return;
}
//****************************************************************************80

void triangle01_poly_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE01_POLY_INTEGRAL_TEST: polynomial integrals over the unit triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d_max = 6;
  int d1 = 1;
  int d2 = 2;
  int d3 = 2;

  int i;
  int j;
  int k;
  int km1;
  int m_max = ( ( d_max + 1 ) * ( d_max + 2 ) ) / 2;
  int m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
  int m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
  int m3 = ( ( d3 + 1 ) * ( d3 + 2 ) ) / 2;
  double p1[3] = { 1.0, 2.0, 3.0 };
  double p2[6] = { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
  double p3[6] = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0 };
  double q;
  double q2;
  double qm[28];

  for ( k = 1; k <= m_max; k++ )
  {
    i4_to_pascal ( k, i, j );
    km1 = k - 1;
    qm[km1] = triangle01_monomial_integral ( i, j );
  }

  cout << "\n";
  cout << "TRIANGLE01_POLY_INTEGRAL_TEST\n";
  cout << "  TRIANGLE01_POLY_INTEGRAL returns the integral Q of\n";
  cout << "  a polynomial P(X,Y) over the interior of the unit triangle.\n";

  cout << "\n";
  poly_print ( d1, p1, "  p(x,y)" );
  q = triangle01_poly_integral ( d1, p1 );
  cout << "\n";
  cout << "  Q =         " << q << "\n";
  q2 = r8vec_dot_product ( m1, p1, qm );
  cout << "  Q (exact) = " << q2 << "\n";

  cout << "\n";
  poly_print ( d2, p2, "  p(x,y)" );
  q = triangle01_poly_integral ( d2, p2 );
  cout << "\n";
  cout << "  Q =         " << q << "\n";
  q2 = r8vec_dot_product ( m2, p2, qm );
  cout << "  Q (exact) = " << q2 << "\n";

  cout << "\n";
  poly_print ( d3, p3, "  p(x,y)" );
  q = triangle01_poly_integral ( d3, p3 );
  cout << "\n";
  cout << "  Q =         " << q << "\n";
  q2 = r8vec_dot_product ( m3, p3, qm );
  cout << "  Q (exact) = " << q2 << "\n";

  return;
}
//****************************************************************************80

void triangle_xy_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_XY_INTEGRAL_TEST tests TRIANGLE_XY_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double q;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  cout << "\n";
  cout << "TRIANGLE_XY_INTEGRAL_TEST\n";
  cout << "  TRIANGLE_XY_INTEGRAL determines Q, the integral of the\n";
  cout << "  monomial X*Y over a triangle (X1,Y1), (X2,Y2), (X3,Y3).\n";

  x1 = 0.0;
  y1 = 0.0;

  x2 = 1.0;
  y2 = 0.0;

  x3 = 1.0;
  y3 = 2.0;

  q = triangle_xy_integral ( x1, y1, x2, y2, x3, y3 );

  cout << "\n";
  cout << "  (X1,Y1) = (" << x1 << "," << y1 << ")\n";
  cout << "  (X2,Y2) = (" << x2 << "," << y2 << ")\n";
  cout << "  (X3,Y3) = (" << x3 << "," << y3 << ")\n";
  cout << "  Q = " << q << "\n";
  cout << "  (Expecting answer 1/2.\n";

  x1 = 0.0;
  y1 = 0.0;

  x2 = 4.0;
  y2 = 0.0;

  x3 = 0.0;
  y3 = 1.0;

  q = triangle_xy_integral ( x1, y1, x2, y2, x3, y3 );

  cout << "\n";
  cout << "  (X1,Y1) = (" << x1 << "," << y1 << ")\n";
  cout << "  (X2,Y2) = (" << x2 << "," << y2 << ")\n";
  cout << "  (X3,Y3) = (" << x3 << "," << y3 << ")\n";
  cout << "  Q = " << q << "\n";
  cout << "  (Expecting answer 2/3.\n";

  return;
}
//****************************************************************************80

void trinomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRINOMIAL_TEST tests TRINOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int t;

  cout << "\n";
  cout << "TRINOMIAL_TEST\n";
  cout << "  TRINOMIAL evaluates the trinomial coefficient:\n";
  cout << "\n";
  cout << "  T(I,J,K) = (I+J+K)! / I! / J! / K!\n";
  cout << "\n";
  cout << "     I     J     K    T(I,J,K)\n";
  cout << "\n";
 
  for ( k = 0; k <= 4; k++ )
  {
    for ( j = 0; j <= 4; j++ )
    {
      for ( i = 0; i <= 4; i++ )
      {
        t = trinomial ( i, j, k );
        cout << "  " << setw(4) << i
             << "  " << setw(4) << j
             << "  " << setw(4) << k
             << "  " << setw(8) << t << "\n";
      }
    }
  }
 
  return;
}
//****************************************************************************80

void xy_to_rs_map_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    XY_TO_RS_MAP_TEST tests XY_TO_RS_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  int j;
  double r;
  double s;
  double t[2*3] = {
    2.0, 0.0, 
    3.0, 4.0, 
    0.0, 3.0 };

  cout << "\n";
  cout << "XY_TO_RS_MAP_TEST:\n";
  cout << "  XY_TO_RS_MAP determines the coefficients of the linear\n";
  cout << "  map from a general triangle in XY coordinates\n";
  cout << "  to the reference triangle in RS coordinates:\n";
  cout << "    R = a + b * X + c * Y\n";
  cout << "    S = d + e * X + f * Y\n";

  r8mat_print ( 2, 3, t, "  XY triangle vertices:" );

  xy_to_rs_map ( t, a, b, c, d, e, f );

  cout << "\n";
  cout << "  Mapping coefficients are:\n";
  cout << "\n";
  cout << "    R = " << a << " + " << b << " * X + " << c << " * Y\n";
  cout << "    S = " << d << " + " << e << " * X + " << f << " * Y\n";

  cout << "\n";
  cout << "  Apply map to XY triangle vertices.\n";
  cout << "  Recover RS vertices (0,0), (1,0) and (0,1).\n";
  cout << "\n";
  for ( j = 0; j < 3; j++ )
  {
    r = a + b * t[0+j*2] + c * t[1+j*2];
    s = d + e * t[0+j*2] + f * t[1+j*2];
    cout << "  V[" << j << "] = (" << r << "," << s << ")\n";
  }

  return;
}
