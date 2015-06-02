# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "r4lib.hpp"

int main ( );

void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test008 ( );
void test009 ( );

void test010 ( );
void test011 ( );
void test012 ( );
void test013 ( );
void test014 ( );
void test015 ( );
void test016 ( );
void test017 ( );
void test018 ( );
void test019 ( );

void test020 ( );
void test023 ( );
void test0235 ( );
void test026 ( );
void test027 ( );
void test028 ( );

void test12555 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for R4LIB_PRB.
//
//  Discussion:
//
//    R4LIB_PRB tests the R4LIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "R4LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the R4LIB library.\n";

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test009 ( );

  test010 ( );
  test011 ( );
  test012 ( );
  test013 ( );
  test014 ( );
  test015 ( );
  test016 ( );
  test017 ( );
  test018 ( );
  test019 ( );

  test020 ( );
  test023 ( );
  test0235 ( );
  test026 ( );
  test027 ( );
  test028 ( );

  test12555 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "R4LIB_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST001 tests R4_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2008
//
//  Author:
//
//    John Burkardt
//
{
  float r4;
  float r4_absolute;
  float r4_hi = 5.0;
  float r4_lo = -3.0;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  cout << "\n";
  cout << "TEST001\n";
  cout << "  R4_ABS returns the absolute value of an R4.\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    r4 = r4_uniform_ab ( r4_lo, r4_hi, seed );
    r4_absolute = r4_abs ( r4 );
    cout << "  " << setw(10) << r4
         << "  " << setw(10) << r4_absolute << "\n";
  }

  return;
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 tests R4_ATAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  float x;
  float y;

  cout << "\n";
  cout << "TEST002\n";
  cout << "  R4_ATAN computes the arc-tangent given Y and X;\n";
  cout << "  ATAN2 is the system version of this routine.\n";
  cout << "\n";
  cout << "           X           Y  ATAN2(Y,X)  ATAN4(Y,X)\n";
  cout << "\n";

  x = 1.0;
  y = 0.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  x = 1.0;
  y = 1.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  x = 0.0;
  y = 1.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  x = -1.0;
  y = 1.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  x = -1.0;
  y = 0.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  x = - 1.0;
  y = - 1.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  x =   0.0;
  y = - 1.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  x =   1.0;
  y = - 1.0;
  cout                               << "  "
       << setw(10) << x              << "  "
       << setw(10) << y              << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << r4_atan ( y, x ) << "\n";

  return;
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests R4_CAS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 12

  int test;
  float x;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  R4_CAS evaluates the casine of a number.\n";
  cout << "\n";
  cout << "        X           R4_CAS ( X )\n";
  cout << "\n";

  for ( test = 0; test <= TEST_NUM; test++ )
  {
    x = r4_pi ( ) * ( float ) ( test ) / ( float ) ( TEST_NUM );
    cout << "  " << setw(14) << x
         << "  " << setw(14) << r4_cas ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test004 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST004 tests R4_CEILING.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  float rval;
  float rval_rounded;

  cout << "\n";
  cout << "TEST004\n";
  cout << "  R4_CEILING rounds a value up.\n";
  cout << "\n";

  for ( i = -6; i <= 6; i++ )
  {
    rval = ( float ) ( i ) / 5.0;
    rval_rounded = r4_ceiling ( rval );
    cout << "  " << setw(14) << rval
         << "  " << setw(14) << rval_rounded << "\n";
  }

  return;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests R4_DIFF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 15

  int ndig = 3;
  int test;
  float x = 1.0;
  float y;
  float y_test[TEST_NUM] = {
    0.0625, 0.125, 0.25, 0.50,  0.874,
    0.876,  0.90,  0.95, 0.99,  1.0,
    1.01,   1.05,  1.10, 3.0,  10.0 };

  cout << "\n";
  cout << "TEST005\n";
  cout << "  R4_DIFF computes a difference X-Y to a given\n";
  cout << "  number of binary places.\n";
  cout << "\n";
  cout << "  For this test, we use " << ndig << " binary places.\n";
  cout << "\n";
  cout << "       X       Y       X-Y     R4_DIFF(X,Y)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    y = y_test[test];
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y
         << "  " << setw(10) << x-y
         << "  " << setw(10) << r4_diff ( x, y, ndig ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests R4_DIGIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define MAXDIG 20

  int idigit;
  float x;

  x = r4_pi ( );

  cout << "\n";
  cout << "TEST006\n";
  cout << "  R4_DIGIT extracts decimal digits.\n";
  cout << "\n";
  cout << "  Here, we get digits of " << x << "\n";
  cout << "\n";

  cout << "  ";
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    cout << setw(3) << idigit;
  }
  cout << "\n";

  cout << "  ";
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    cout << setw(3) << r4_digit ( x, idigit );
  }
  cout << "\n";

  return;
# undef MAXDIG
}
//****************************************************************************80

void test007 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests R4_EPSILON
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  float r;
  float r2;
  float s;
  float t;

  cout << "\n";
  cout << "TEST007\n";
  cout << "  R4_EPSILON produces the R4 roundoff unit.\n";
  cout << "\n";

  r = r4_epsilon ( );
  cout << "  R = R4_EPSILON()  = " << setw(10) << r << "\n";

  s = 1.0 + r;
  t = s - 1.0;
  cout << "  ( 1 + R ) - 1     = " << setw(10) << t << "\n";

  r2 = r / 2.0;
  s = 1.0 + r2;
  t = s - 1.0;
  cout << "  ( 1 + (R/2) ) - 1 = " << setw(10) << t << "\n";

  return;
}
//****************************************************************************80

void test008 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST008 tests R4_FRACTIONAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  float fractional;
  float r4;
  float r4_hi = 5.0;
  float r4_lo = -3.0;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST008\n";
  cout << "  R4_FRACTIONAL returns the fractional part of an R4.\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    r4 = r4_uniform_ab ( r4_lo, r4_hi, seed );
    fractional = r4_fractional ( r4 );
    cout << "  " << setw(10) << r4
         << "  " << setw(10) << fractional << "\n";
  }

  return;
}
//****************************************************************************80

void test009 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST009 tests R4_HUGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST009\n";
  cout << "  R4_HUGE returns a large R4 value;\n";
  cout << "\n";
  cout << "  R4_HUGE =   " << r4_huge ( ) << "\n";

  return;
}
//****************************************************************************80

void test010 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST010 tests R4_LOG_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 18

  int test;
  float x;
  float x_test[TEST_NUM] = {
    0.0,  1.0,  2.0,   3.0,  9.0,
   10.0, 11.0, 99.0, 101.0, -1.0,
   -2.0, -3.0, -9.0,   0.5,  0.33,
    0.25, 0.20, 0.01 };

  cout << "\n";
  cout << "TEST010\n";
  cout << "  R4_LOG_2: computes the logarithm base 2.\n";
  cout << "\n";
  cout << "  X       R4_LOG_2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setw(12) << x
         << "  " << setw(12) << r4_log_2 ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test011 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST011 tests R4_LOG_B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 10

  float b;
  float b_test[TEST_NUM] = {
    2.0, 3.0, 4.0, 5.0, 6.0,
    7.0, 8.0, 16.0, 32.0, 256.0 };
  int test;
  float x;

  x = 16.0;

  cout << "\n";
  cout << "TEST011\n";
  cout << "  R4_LOG_B computes the logarithm base B.\n";
  cout << "\n";
  cout << "  X, B, R4_LOG_B\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = b_test[test];

    cout << "  " << setw(12) << x
         << "  " << setw(12) << b
         << "  " << setw(12) << r4_log_b ( x, b ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test012 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST012 tests R4_MANT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int is;
  int l;
  float r;
  float x;

  x = -314.159;

  cout << "\n";
  cout << "TEST012\n";
  cout << "  R4_MANT decomposes a value.\n";
  cout << "\n";
  cout << "  Number to be decomposed: X = " << x << "\n";

  r4_mant ( x, &is, &r, &l );

  cout << "\n";
  cout << "  X = " << is << " * " << r << " * 2 ^ " << l << "\n";

  return;
}
//****************************************************************************80

void test013 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST013 tests R4_MOD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float x_hi = 10.0;
  float x_lo = -10.0;
  float y;
  float z1;
  float z2;

  cout << "\n";
  cout << "TEST013\n";
  cout << "  R4_MOD returns the remainder after division.\n";
  cout << "  R4_MOD ( X, Y ) has the same sign as X.\n";
  cout << "\n";
  cout << "      X         Y    FMOD(X,Y)    R4_MOD(X,Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform_ab ( x_lo, x_hi, seed );
    y = r4_uniform_ab ( x_lo, x_hi, seed );

    z1 =   fmod ( x, y );
    z2 = r4_mod ( x, y );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << y
         << "  " << setw(12) << z1
         << "  " << setw(12) << z2 << "\n";
  }

  return;
}
//****************************************************************************80

void test014 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST014 tests R4_MODP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float x_hi = 10.0;
  float x_lo = -10.0;
  float y;
  float z1;
  float z2;

  cout << "\n";
  cout << "TEST014\n";
  cout << "  R4_MODP returns the remainder after division.\n";
  cout << "  R4_MODP ( X, Y ) is positive.\n";
  cout << "\n";
  cout << "      X       Y     FMOD(X,Y)  R4_MODP(X,Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform_ab ( x_lo, x_hi, seed );
    y = r4_uniform_ab ( x_lo, x_hi, seed );

    z1 =   fmod  ( x, y );
    z2 = r4_modp ( x, y );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << y
         << "  " << setw(12) << z1
         << "  " << setw(12) << z2 << "\n";
  }

  return;
}
//****************************************************************************80

void test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests R4_NINT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  float b;
  float c;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;

  cout << "\n";
  cout << "TEST015\n";
  cout << "  R4_NINT produces the nearest integer.\n";
  cout << "\n";

  b = -10.0;
  c = +10.0;

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform_ab ( b, c, seed );
    cout << setw(10) << x << "  "
         << setw(6)  << r4_nint ( x ) << "\n";
  }

  return;
}
//****************************************************************************80

void test016 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST016 tests R4_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 20

  int seed = 123456789;
  int test;
  float x;

  cout << "\n";
  cout << "TEST016\n";
  cout << "  R4_NORMAL_01 generates normally distributed random values.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = r4_normal_01 ( seed );
    cout                  << "  "
         << setw(10) << x << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test017 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST017 tests R4_PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  float four;
  float one;
  float v1;
  float v2;

  four = ( float ) ( 4 );
  one = ( float ) ( 1 );

  cout << "\n";
  cout << "TEST017\n";
  cout << "  R4_PI returns the value of PI.\n";
  cout << "\n";
  v1 = r4_pi ( );
  cout << "  R4_PI =     " << v1 << "\n";
  v2 = four * atan ( one );
  cout << "  4*atan(1) = " << v2 << "\n";

  return;
}
//****************************************************************************80

void test018 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST018 tests R4_POWER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int p;
  float r;
  float value;

  cout << "\n";
  cout << "TEST018\n";
  cout << "  R4_POWER computes R^P\n";
  cout << "\n";
  cout << "      R          P       R**P\n";
  cout << "\n";

  for ( p = -5; p <= 5; p++ )
  {
    r = 2.0;
    value = r4_power ( r, p );
    cout << "  " << setw(12) << r
         << "  " << setw(6)  << p
         << "  " << setw(12) << value << "\n";
  }

  return;
}
//****************************************************************************80

void test019 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST019 tests R4_POWER_FAST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int mults;
  int p;
  float r;
  float rp;

  cout << "\n";
  cout << "TEST019\n";
  cout << "  R4_POWER_FAST computes R^P, economizing on\n";
  cout << "    multiplications.\n";
  cout << "\n";
  cout << "      R          P       R**P       Mults\n";
  cout << "\n";

  for ( i = -10; i <= 40; i++ )
  {
    r = 2.0;
    p = i;
    rp = r4_power_fast ( r, p, &mults );
    cout << "  " << setw(12) << r
         << "  " << setw(6)  << p
         << "  " << setw(12) << rp
         << "  " << setw(6)  << mults << "\n";
  }

  return;
}
//****************************************************************************80

void test020 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST020 tests R4_ROUND2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nplace;
  float x;
  float xround;

  x = r4_pi ( );

  cout << "\n";
  cout << "TEST020\n";
  cout << "  R4_ROUND2 rounds a number to a\n";
  cout << "    specified number of base 2 digits.\n";
  cout << "\n";
  cout << "  Test effect on PI:\n";
  cout << "  X = " << x << "\n";
  cout << "\n";
  cout << "  NPLACE  XROUND\n";
  cout << "\n";

  for ( i = 0; i <= 20; i++ )
  {
    nplace = i;
    xround = r4_round2 ( nplace, x );
    cout << "  " << setw(8) << i
         << "  " << setprecision(16) << setw(24) << xround << "\n";
  }

  return;
}
//****************************************************************************80

void test023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023 tests R4_SIGN and R4_SIGN3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  float r4;
  float r4_test[5] = { -1.25E+00, -0.25E+00, 0.0E+00, +0.5E+00, +9.0E+00 };
  float s1;
  float s2;
  int test;
  const int test_num = 5;

  cout << "\n";
  cout << "TEST023\n";
  cout << "  R4_SIGN returns the sign of an R4.\n";
  cout << "  R4_SIGN3 returns the three-way sign of an R4.\n";
  cout << "\n";
  cout << "      R4    R4_SIGN(R4)  R4_SIGN3(R4)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    r4 = r4_test[test];
    s1 = r4_sign ( r4 );
    s2 = r4_sign3 ( r4 );
    cout << setw(10) << r4 << "  "
         << setw(8) << s1 << "  "
         << setw(8) << s2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0235 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0235 tests R4_SWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  float x;
  float y;

  cout << "\n";
  cout << "TEST0235\n";
  cout << "  R4_SWAP swaps two reals.\n";

  x = 1.0;
  y = 3.14159;

  cout << "\n";
  cout << "  Before swapping: \n";
  cout << "\n";
  cout << "    X = " << x << "\n";
  cout << "    Y = " << y << "\n";

  r4_swap ( x, y );

  cout << "\n";
  cout << "  After swapping: \n";
  cout << "\n";
  cout << "    X = " << x << "\n";
  cout << "    Y = " << y << "\n";

  return;
}
//****************************************************************************80

void test026 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026 tests R4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  float a;
  float b;
  float c;
  int i;
  int seed;

  b = 10.0;
  c = 25.0;
  seed = 17;

  cout << "\n";
  cout << "TEST026\n";
  cout << "  R4_UNIFORM_AB produces a random real in a given range.\n";
  cout << "\n";
  cout << "  Using range " << b << " <= A <= " << c << ".\n";
  cout << "\n";

  cout << "\n";
  cout << "  I   A\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    a = r4_uniform_ab ( b, c, seed );
    cout << setw ( 6 )  << i << " "
         << setw ( 10 ) << a << "\n";
  }

  return;
}
//****************************************************************************80

void test027 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST027 tests R4_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  float x;

  cout << "\n";
  cout << "TEST027\n";
  cout << "  R4_UNIFORM_01 produces a sequence of random values.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  cout << "\n";
  cout << "  SEED   R4_UNIFORM_01(SEED)\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << setw(12) << seed << "  ";
    x = r4_uniform_01 ( seed );
    cout << setw(10) << x << "\n";
  }

  cout << "\n";
  cout << "  Verify that the sequence can be restarted.\n";
  cout << "  Set the seed back to its original value, and see that\n";
  cout << "  we generate the same sequence.\n";

  seed = 123456789;
  cout << "\n";
  cout << "  SEED   R4_UNIFORM_01(SEED)\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    cout << setw(12) << seed << "  ";
    x = r4_uniform_01 ( seed );
    cout << setw(10) << x << "\n";
  }

  return;
}
//****************************************************************************80

void test028 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST028 tests R4_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 1000

  int i;
  float max;
  float mean;
  float min;
  int n;
  int seed = 123456789;
  float x[N];
  float variance;

  cout << "\n";
  cout << "TEST028\n";
  cout << "  R4_UNIFORM_01 samples a uniform random distribution in [0,1].\n";
  cout << "  distributed random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = r4_uniform_01 ( seed );
  }

  cout << "\n";
  cout << "  First few values:\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  min = r4vec_min ( N, x );
  max = r4vec_max ( N, x );
  mean = r4vec_mean ( N, x );
  variance = r4vec_variance ( N, x );

  cout << "\n";
  cout << "  Number of samples was " << N << "\n";
  cout << "  Minimum value was " << min << "\n";
  cout << "  Maximum value was " << max << "\n";
  cout << "  Average value was " << mean << "\n";
  cout << "  Variance was      " << variance << "\n";

  return;
# undef N
}
//****************************************************************************80

void test12555 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12555 tests R4VEC_INDICATOR0_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  float *v;

  cout << "\n";
  cout << "TEST012555\n";
  cout << "  R4VEC_INDICATOR0_NEW returns an indicator vector.\n";

  n = 10;
  v = r4vec_indicator0_new ( n );
  r4vec_print ( n, v, "  Indicator0 vector:" );
  delete [] v;

  return;
}
