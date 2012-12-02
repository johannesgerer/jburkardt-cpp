# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "c8lib.hpp"

int main ( );

void test0061 ( );
void test0062 ( );
void test0063 ( );
void test0064 ( );
void test0065 ( );
void test0066 ( );
void test0067 ( );
void test102 ( );
void test103 ( );
void test104 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for C8LIB_PRB.
//
//  Discussion:
//
//    C8LIB_PRB calls the C8LIB tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "C8LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the C8LIB library.\n";

  test0061 ( );
  test0062 ( );
  test0063 ( );
  test0064 ( );
  test0065 ( );
  test0066 ( );
  test0067 ( );

  test102 ( );
  test103 ( );
  test104 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "C8LIB_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test0061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0061 tests C8_ARGUMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double argument;
  complex <double> c1;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0061\n";
  cout << "  C8_ARGUMENT computes the argument of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            ARG=C8_ARGUMENT(C1)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( &seed );
    argument = c8_argument ( c1 );

    cout << "  (" << setw(8) << real ( c1 )
         << ",  " << setw(8) << imag ( c1 ) << ")"
         << "  " << setw(8) << argument << "\n";
  }

  return;
}
//****************************************************************************80

void test0062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0062 tests C8_CUBE_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0062\n";
  cout << "  C8_CUBE_ROOT computes the principal cube root of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            C2=C8_CUBE_ROOT(C1)         C3=C2*C2*C2\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_cube_root ( c1 );
    c3 = c2 * c2 * c2;

    cout << "  (" << setw(8) << real ( c1 )
         << ",  " << setw(8) << imag ( c1 ) << ")"
         << "  (" << setw(8) << real ( c2 )
         << ",  " << setw(8) << imag ( c2 ) << ")"
         << "  (" << setw(8) << real ( c3 )
         << ",  " << setw(8) << imag ( c3 ) << ")\n";
  }

  return;
}
//****************************************************************************80

void test0063 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0063 tests C8_MAGNITUDE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  double magnitude;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0063\n";
  cout << "  C8_MAGNITUDE computes the magnitude of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            MAG=C8_MAGNITUDE(C1)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( &seed );
    magnitude = c8_magnitude ( c1 );

    cout << "  (" << setw(8) << real ( c1 )
         << ",  " << setw(8) << imag ( c1 ) << ")"
         << "  " << setw(8) << magnitude << "\n";
  }

  return;
}
//****************************************************************************80

void test0064 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0064 tests C8_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 20;
  complex <double> x;

  cout << "\n";
  cout << "TEST0064\n";
  cout << "  C8_NORMAL_01 generates unit pseudonormal\n";
  cout << "    complex values.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = c8_normal_01 ( &seed );
    cout << "  " << setw(8) << real ( x ) 
         << "  " << setw(8) << imag ( x ) << "\n";
  }

  return;
}
//****************************************************************************80

void test0065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0065 tests C8_SQRT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> c1;
  complex <double> c2;
  complex <double> c3;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0065\n";
  cout << "  C8_SQRT computes the principal square root of a C8.\n";
  cout << "\n";
  cout << 
    "            C1=random            C2=C8_SQRT(C1)              C3=C2*C2\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_sqrt ( c1 );
    c3 = c2 * c2;

    cout << "  (" << setw(8) << real ( c1 )
         << ",  " << setw(8) << imag ( c1 ) << ")"
         << "  (" << setw(8) << real ( c2 )
         << ",  " << setw(8) << imag ( c2 ) << ")"
         << "  (" << setw(8) << real ( c3 )
         << ",  " << setw(8) << imag ( c3 ) << ")\n";
  }

  return;
}
//****************************************************************************80

void test0066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0066 tests C8MAT_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int m = 5;
  int n = 4;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST0066\n";
  cout << "  C8MAT_UNIFORM_01_NEW computes a random complex matrix.\n";

  a = c8mat_uniform_01_new ( m, n, &seed );

  c8mat_print ( m, n, a, "  The matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void test0067 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0067 tests C8VEC_INDICATOR_NEW;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *a;
  int n = 10;

  cout << "\n";
  cout << "TEST0067\n";
  cout << "  C8VEC_INDICATOR_NEW sets A = (1-1i,2-2i,...,N-Ni)\n";

  a = c8vec_indicator_new ( n );
 
  c8vec_print ( n, a, "  The indicator vector:" );

  delete [] a;

  return;
}
//****************************************************************************80

void test102 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST102 tests R8POLY2_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double a;
  double a_test[TEST_NUM] = { 2.0, 1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = { -2.0, -20.0, -2.0 };
  double c;
  double c_test[TEST_NUM] = { -24.0, 100.0, 10.0 };
  complex <double> r1;
  complex <double> r2;
  int test;

  cout << "\n";
  cout << "TEST102\n";
  cout << "  R8POLY2_ROOT finds quadratic equation roots.\n";
  cout << "\n";
  cout << "         A         B         C     R1         R2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];

    r8poly2_root ( a, b, c, &r1, &r2 );

    cout << "  "   << setw(8) << a
         << "  "   << setw(8) << b
         << "  "   << setw(8) << c
         << "  ("  << setw(8) << real ( r1 )
         << "  ,"  << setw(8) << imag ( r1 )
         << ")  (" << setw(8) << real ( r2 )
         << ",  "  << setw(8) << imag ( r2 ) << ")\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test103 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST103 tests R8POLY3_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double a;
  double a_test[TEST_NUM] = { 1.0, 9.0, 1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = { -6.0, -36.0, -5.0, -8.0 };
  double c;
  double c_test[TEST_NUM] = { 11.0, 54.0, 8.0, 25.0 };
  double d;
  double r8_test[TEST_NUM] = { -6.0, -27.0, -4.0, -26.0 };
  complex <double> r1;
  complex <double> r2;
  complex <double> r3;
  int test;
//
//  1: Three distinct real roots, 1, 2, 3.
//  2: One repeated real root, 1.5, 1.5, 1.5.
//  3: Two real roots, one repeated, 1, 2, 2.
//  4: One real root, a complex conjugate pair, 2, 3+2I, 3-2I.
//
  cout << "\n";
  cout << "TEST103\n";
  cout << "  R8POLY3_ROOT finds roots of cubic equations.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];

    cout << "\n";
    cout << "  Polynomial coefficients A, B, C, D:\n";
    cout << "\n";

    cout << "  A = " << setw(8) << a << "\n";
    cout << "  B = " << setw(8) << b << "\n";
    cout << "  C = " << setw(8) << c << "\n";
    cout << "  D = " << setw(8) << d << "\n";

    r8poly3_root ( a, b, c, d, &r1, &r2, &r3 );

    cout << "\n";
    cout << "  Roots:\n";
    cout << "\n";
    cout << "  (" << real ( r1 )
         << ",  " << imag ( r1 ) << ")\n";
    cout << "  (" << real ( r2 )
         << ",  " << imag ( r2 ) << ")\n";
    cout << "  (" << real ( r3 )
         << ",  " << imag ( r3 ) << ")\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test104 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST104 tests R8POLY4_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 7

  double a;
  double a_test[TEST_NUM] = {
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = {
    -10.0, -5.0, -22.0, -16.0, -20.0,
    2.0, 0.0 };
  double c;
  double c_test[TEST_NUM] = {
    35.0, 1.0, 141.0, 72.0, 150.0,
    1.0, 13.0 };
  double d;
  double r8_test[TEST_NUM] = {
    -50.0, 21.0, -220.0, -128.0, -500.0,
    8.0, 0.0 };
  double e;
  double e_test[TEST_NUM] = {
    24.0, -18.0, +100.0, 80.0, 625.0,
    -12.0, 36.0 };
  complex <double> r1;
  complex <double> r2;
  complex <double> r3;
  complex <double> r4;
  int test;
//
//  1: Four distinct real roots, 1, 2, 3, 4.
//  2: Three distinct real roots, 1, -2, 3, 3
//  3: Two distinct real roots, 1, 1, 10, 10.
//  4: Two distinct real roots, 2, 2, 2, 10
//  5: One real root, 5, 5, 5, 5
//  6: Two distinct real roots, one complex conjugate pair.
//  7: Two distinct complex conjugate pairs.
//
  cout << "\n";
  cout << "TEST104\n";
  cout << "  R8POLY4_ROOT finds roots of quartic equations.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];
    e = e_test[test];

    cout << "\n";
    cout << "  A = " << setw(8) << a << "\n";
    cout << "  B = " << setw(8) << b << "\n";
    cout << "  C = " << setw(8) << c << "\n";
    cout << "  D = " << setw(8) << d << "\n";
    cout << "  E = " << setw(8) << e << "\n";

    r8poly4_root ( a, b, c, d, e, &r1, &r2, &r3, &r4 );

    cout << "\n";
    cout << "  Roots:\n";
    cout << "\n";
    cout << "  (" << real ( r1 )
         << ",  " << imag ( r1 ) << ")\n";
    cout << "  (" << real ( r2 )
         << ",  " << imag ( r2 ) << ")\n";
    cout << "  (" << real ( r3 )
         << ",  " << imag ( r3 ) << ")\n";
    cout << "  (" << real ( r4 )
         << ",  " << imag ( r4 ) << ")\n";
  }

  return;
# undef TEST_NUM
}
