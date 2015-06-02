# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "piecewise_linear_product_integral.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB.
//
//  Discussion:
//
//    PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB tests the 
//    PIECEWISE_LINEAR_PRODUCT_INTEGRAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the PIECEWISE_LINEAR_PRODUCT_INTEGRAL_INTEGRAL library.\n";
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB\n";
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
//    TEST01 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    For the first test, we use the same single "piece" for both F and G.
//    Hence, we are actually integrating X^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
{
# define F_NUM 2
# define G_NUM 2

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM] = { 0.0, 5.0 };
  double f_x[F_NUM] = { 0.0, 5.0 };
  int g_num = G_NUM;
  double g_v[G_NUM] = { 0.0, 5.0 };
  double g_x[G_NUM] = { 0.0, 5.0 };
  int i;
  double integral;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a very simple problem.\n";
  cout << "  F and G are both defined over a single common\n";
  cout << "  interval, so that F(X) = G(X) = X.\n";
  cout << "\n";
  cout << "           A           B      Integral        Exact\n";
  cout << "\n";

  a = 1.0;
  for ( i = 1; i <= 5; i++ )
  {
    b = ( double ) ( i );
    integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, 
      g_x, g_v );
    exact = ( b * b * b - a * a * a ) / 3.0;
    cout << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(14) << integral
         << "  " << setw(14) << exact << "\n";
  }

  return;
# undef F_NUM
# undef G_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    For this test, we use multiple "pieces" for both F and G,
//    but we define the values so that we are still actually integrating X^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
{
# define F_NUM 3
# define G_NUM 4

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM] = { 0.0, 2.0, 5.0 };
  double f_x[F_NUM] = { 0.0, 2.0, 5.0 };
  int g_num = G_NUM;
  double g_v[G_NUM] = { 0.0, 1.5, 3.0, 5.0 };
  double g_x[G_NUM] = { 0.0, 1.5, 3.0, 5.0 };
  int i;
  double integral;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.\n";
  cout << "  F and G are both defined over separate, multiple\n";
  cout << "  intervals, but still true that F(X) = G(X) = X.\n";
  cout << "\n";
  cout << "           A           B      Integral        Exact\n";
  cout << "\n";

  a = 1.0;
  for ( i = 1; i <= 5; i++ )
  {
    b = ( double ) ( i );
    integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, 
      g_x, g_v );
    exact = ( b * b * b - a * a * a ) / 3.0;
    cout << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(14) << integral
         << "  " << setw(14) << exact << "\n";
  }

  return;
# undef F_NUM
# undef G_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    For this test, F(X) and G(X) are piecewise linear interpolants to
//    SIN(X) and 2 * COS(X), so we know the exact value of the integral
//    of the product of the original functions, but this is only an estimate 
//    of the exact value of the integral of the product of the interpolants.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 April 2009
//
//  Author:
//
//    John Burkardt
//
{
# define F_NUM 11
# define G_NUM 31

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM];
  double f_x[F_NUM];
  int g_num = G_NUM;
  double g_v[G_NUM];
  double g_x[G_NUM];
  int i;
  double integral;
  double pi = 3.141592653589793;
  double quad;
  int quad_num;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.\n";
  cout << "  F and G are defined over separate, multiple\n";
  cout << "  intervals.\n";
  cout << "\n";
  cout << "  F(X) interpolates SIN(X),\n";
  cout << "  G(X) interpolates 2*COS(X).\n";
  cout << "\n";
  cout << "  We compare:\n";
  cout << "\n";
  cout << "  INTEGRAL, our value for the integral,\n";
  cout << "  QUAD, a quadrature estimate for the integral, and\n";
  cout << "  CLOSE, the value of the integral of 2*COS(X)*SIN(X)\n";
  cout << "\n";
  cout << "           A           B      Integral        Quad            Close\n";
  cout << "\n";

  for ( i = 0; i < f_num; i++ )
  {
    f_x[i] = ( ( f_num - i - 1 ) * 0.0
             + (         i     ) * pi )
             / ( f_num     - 1 );
    f_v[i] = sin ( f_x[i] );
  }

  for ( i = 0; i < g_num; i++ )
  {
    g_x[i] = ( ( g_num - i - 1 ) * 0.0
             + (         i     ) * pi )
             / ( g_num     - 1 );
    g_v[i] = 2.0 * cos ( g_x[i] );
  }

  a = 0.0;
  for ( i = 0; i <= 6; i++ )
  {
    b = ( double ) ( i ) * pi / 6.0;
    integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, 
      g_num, g_x, g_v );
    exact = - ( cos ( 2.0 * b ) - cos ( 2.0 * a ) ) / 2.0;
    quad_num = 2000;
    quad = piecewise_linear_product_quad ( a, b, f_num, f_x, f_v, g_num, 
      g_x, g_v, quad_num );
    cout << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(14) << integral
         << "  " << setw(14) << quad
         << "  " << setw(14) << exact << "\n";
  }

  return;
# undef F_NUM
# undef G_NUM
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    For this test, we compute the integrals of a hat function with itself,
//    and a hat function with its neighbor.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
{
# define F_NUM 3
# define G_NUM 3

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM] = { 0.0, 1.0, 0.0 };
  double f_x[F_NUM] = { 0.0, 1.0, 2.0 };
  int g_num = G_NUM;
  double g_v[G_NUM] = { 1.0, 0.0, 0.0 };
  double g_x[G_NUM] = { 0.0, 1.0, 2.0 };
  int i;
  double integral;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL.\n";
  cout << "  The nodes are at 0, 1, and 2.\n";
  cout << "  F(X) = ( 0, 1, 0 ).\n";
  cout << "  G(X) = ( 1, 0, 0 ).\n";
  cout << "\n";

  a = 0.0;
  b = 2.0;

  integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, f_num, 
    f_x, f_v );

  cout << "  Integral F(X) * F(X) dx = " << integral << "\n";

  integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, 
    g_x, g_v );

  cout << "  Integral F(X) * G(X) dx = " << integral << "\n";

  integral = piecewise_linear_product_integral ( a, b, g_num, g_x, g_v, g_num, 
    g_x, g_v );

  cout << "  Integral G(X) * G(X) dx = " << integral << "\n";

  return;
# undef F_NUM
# undef G_NUM
}
