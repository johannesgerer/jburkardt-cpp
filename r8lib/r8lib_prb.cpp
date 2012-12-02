# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "r8lib.hpp"

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
void test021 ( );
void test022 ( );
void test023 ( );
void test0235 ( );
void test024 ( );
void test025 ( );
void test026 ( );
void test027 ( );
void test028 ( );
void test029 ( );
void test0295 ( );

void test031 ( );
void test032 ( );
void test033 ( );
void test034 ( );
void test035 ( );
void test036 ( );
void test0363 ( );
void test0365 ( );
void test037 ( );
void test038 ( );
void test0385 ( );
void test039 ( );
void test0393 ( );
void test0395 ( );
void test0397 ( );

void test040 ( );
void test041 ( );
void test0415 ( );
void test042 ( );
void test043 ( );
void test044 ( );
void test0442 ( );
void test0443 ( );
void test0445 ( );
void test045 ( );
void test046 ( );
void test047 ( );
void test048 ( );
void test049 ( );

void test050 ( );
void test051 ( );
void test052 ( );
void test053 ( );
void test054 ( );
void test055 ( );
void test0555 ( );
void test056 ( );
void test057 ( );
void test058 ( );
double test058_f ( int n, double x[] );
double *test058_hess ( int n, double x[] );
void test059 ( );

void test060 ( );
void test061 ( );
void test062 ( );
void test063 ( );
void test064 ( );
void test065 ( );
void test066 ( );
double *test067_f ( int m, int n, double x[] );
double *test067_jac ( int m, int n, double x[] );
void test067 ( );
void test068 ( );
void test069 ( );

void test070 ( );
void test071 ( );
void test072 ( );
void test073 ( );
void test0731 ( );
void test0732 ( );
void test0733 ( );
void test0734 ( );
void test0735 ( );
void test0736 ( );
void test07365 ( );
void test0737 ( );
void test074 ( );
void test075 ( );
void test076 ( );
void test0764 ( );
void test0766 ( );
void test077 ( );
void test0775 ( );
void test0776 ( );
void test078 ( );
void test079 ( );

void test080 ( );
void test081 ( );
void test082 ( );
void test083 ( );
void test084 ( );
void test085 ( );
void test086 ( );
void test087 ( );
void test088 ( );
void test089 ( );

void test090 ( );
void test091 ( );
void test092 ( );
void test093 ( );
void test094 ( );
void test095 ( );
void test098 ( );
void test099 ( );

void test100 ( );
void test100_f ( double x, double *y, double *yp, double *ypp );
void test101 ( );
void test105 ( );
void test106 ( );
void test107 ( );
void test108 ( );
void test109 ( );

void test110 ( );
void test111 ( );
void test112 ( );
void test113 ( );
void test114 ( );
void test1143 ( );
void test1145 ( );
void test1147 ( );
void test115 ( );
void test116 ( );
double test116_f ( double x );
void test1165 ( );
void test1166 ( );
void test117 ( );
void test118 ( );

void test120 ( );
void test121 ( );
void test122 ( );
void test123 ( );
void test124 ( );
void test125 ( );
void test1251 ( );
void test1252 ( );
void test1255 ( );
void test1256 ( );
void test1258 ( );
void test126 ( );
void test127 ( );
void test128 ( );
void test129 ( );

void test130 ( );
void test131 ( );
void test132 ( );
void test133 ( );
void test134 ( );
void test135 ( );
void test136 ( );
void test137 ( );
void test138 ( );
void test139 ( );

void test140 ( );
void test141 ( );
void test142 ( );
void test143 ( );
void test144 ( );
void test145 ( );
void test146 ( );
void test1465 ( );
void test147 ( );
void test1475 ( );
void test148 ( );
void test149 ( );

void test150 ( );
void test1504 ( );
void test1505 ( );
void test151 ( );
void test152 ( );
void test153 ( );
void test154 ( );
void test155 ( );
void test156 ( );
void test157 ( );
void test158 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for R8LIB_PRB.
//
//  Discussion:
//
//    R8LIB_PRB calls the R8LIB tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "R8LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the R8LIB library.\n";

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
  test021 ( );
  test022 ( );
  test023 ( );
  test0235 ( );
  test024 ( );
  test025 ( );
  test026 ( );
  test027 ( );
  test028 ( );
  test029 ( );
  test0295 ( );

  test031 ( );
  test032 ( );
  test033 ( );
  test034 ( );
  test035 ( );
  test036 ( );
  test0363 ( );
  test0365 ( );
  test037 ( );
  test038 ( );
  test0385 ( );
  test039 ( );
  test0393 ( );
  test0395 ( );
  test0397 ( );

  test040 ( );
  test041 ( );
  test0415 ( );
  test042 ( );
  test043 ( );
  test044 ( );
  test0442 ( );
  test0443 ( );
  test0445 ( );
  test045 ( );
  test046 ( );
  test047 ( );
  test048 ( );
  test049 ( );

  test050 ( );
  test051 ( );
  test052 ( );
  test053 ( );
  test054 ( );
  test055 ( );
  test0555 ( );
  test056 ( );
  test057 ( );
  test058 ( );
  test059 ( );

  test060 ( );
  test061 ( );
  test062 ( );
  test063 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test067 ( );
  test068 ( );
  test069 ( );

  test070 ( );
  test071 ( );
  test072 ( );
  test073 ( );
  test0731 ( );
  test0732 ( );
  test0733 ( );
  test0734 ( );
  test0735 ( );
  test0736 ( );
  test07365 ( );
  test0737 ( );
  test074 ( );
  test075 ( );
  test076 ( );
  test0764 ( );
  test0766 ( );
  test077 ( );
  test0775 ( );
  test0776 ( );
  test078 ( );
  test079 ( );

  test080 ( );
  test081 ( );
  test082 ( );
  test083 ( );
  test084 ( );
  test085 ( );
  test086 ( );
  test086 ( );
  test088 ( );
  test089 ( );

  test090 ( );
  test091 ( );
  test092 ( );
  test093 ( );
  test094 ( );
  test095 ( );
  test098 ( );
  test099 ( );

  test100 ( );
  test101 ( );
  test105 ( );
  test106 ( );
  test107 ( );
  test108 ( );
  test109 ( );

  test110 ( );
  test111 ( );
  test112 ( );
  test113 ( );
  test114 ( );
  test1143 ( );
  test1145 ( );
  test1147 ( );
  test115 ( );
  test116 ( );
  test1165 ( );
  test1166 ( );
  test117 ( );
  test118 ( );

  test120 ( );
  test121 ( );
  test122 ( );
  test123 ( );
  test124 ( );
  test125 ( );
  test1251 ( );
  test1252 ( );
  test1255 ( );
  test1256 ( );
  test1258 ( );
  test126 ( );
  test127 ( );
  test128 ( );
  test129 ( );

  test130 ( );
  test152 ( );
  test131 ( );
  test132 ( );
  test133 ( );
  test134 ( );
  test135 ( );
  test136 ( );
  test137 ( );
  test138 ( );
  test139 ( );

  test140 ( );
  test141 ( );
  test142 ( );
  test143 ( );
  test144 ( );
  test145 ( );
  test146 ( );
  test1465 ( );
  test147 ( );
  test1475 ( );
  test148 ( );
  test149 ( );

  test150 ( );
  test1504 ( );
  test1505 ( );
  test151 ( );
  test153 ( );
  test154 ( );
  test155 ( );
  test156 ( );
  test157 ( );
  test158 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "R8LIB_PRB\n";
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
//    TEST001 tests R8_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double r8;
  double r8_absolute;
  double r8_hi = 5.0;
  double r8_lo = -3.0;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  cout << "\n";
  cout << "TEST001\n";
  cout << "  R8_ABS returns the absolute value of an R8.\n";
  cout << "\n";
  cout << "      X         R8_ABS(X)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    r8 = r8_uniform_ab ( r8_lo, r8_hi, seed );
    r8_absolute = r8_abs ( r8 );
    cout << "  " << setw(10) << r8
         << "  " << setw(10) << r8_absolute << "\n";
  }

  return;
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 tests R8_ATAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 8

  int test;
  double x;
  double xtest[TEST_NUM] = {
     1.0,  1.0,  0.0, -1.0,
    -1.0, -1.0,  0.0,  1.0 };
  double y;
  double ytest[TEST_NUM] = {
     0.0,  1.0,  1.0,  1.0,
     0.0, -1.0, -1.0, -1.0 };

  cout << "\n";
  cout << "TEST002\n";
  cout << "  R8_ATAN computes the arc-tangent given Y and X;\n";
  cout << "  ATAN2 is the system version of this routine.\n";
  cout << "\n";
  cout << "       X             Y          ATAN2(Y,X)    R8_ATAN(Y,X)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = xtest[test];
    y = ytest[test];
    cout << "  " << setw(14) << x
         << "  " << setw(14) << y
         << "  " << setw(14) << atan2 ( y, x )
         << "  " << setw(14) << r8_atan ( y, x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests R8_CAS.
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
# define TEST_NUM 12

  int test;
  double x;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  R8_CAS evaluates the casine of a number.\n";
  cout << "\n";
  cout << "        X           R8_CAS ( X )\n";
  cout << "\n";

  for ( test = 0; test <= TEST_NUM; test++ )
  {
    x = r8_pi ( ) * ( double ) ( test ) / ( double ) ( TEST_NUM );
    cout << "  " << setw(14) << x
         << "  " << setw(14) << r8_cas ( x ) << "\n";
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
//    TEST004 tests R8_CEILING.
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
  int i;
  double rval;
  double rval_rounded;

  cout << "\n";
  cout << "TEST004\n";
  cout << "  R8_CEILING rounds a value up.\n";
  cout << "\n";

  for ( i = -6; i <= 6; i++ )
  {
    rval = ( double ) ( i ) / 5.0;
    rval_rounded = r8_ceiling ( rval );
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
//    TEST005 tests R8_DIFF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 15

  int ndig = 3;
  int test;
  double x = 1.0;
  double y;
  double y_test[TEST_NUM] = {
    0.0625, 0.125, 0.25, 0.50,  0.874,
    0.876,  0.90,  0.95, 0.99,  1.0,
    1.01,   1.05,  1.10, 3.0,  10.0 };

  cout << "\n";
  cout << "TEST005\n";
  cout << "  R8_DIFF computes a difference X-Y to a given\n";
  cout << "  number of binary places.\n";
  cout << "\n";
  cout << "  For this test, we use " << ndig << " binary places.\n";
  cout << "\n";
  cout << "       X       Y       X-Y     R8_DIFF(X,Y)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    y = y_test[test];
    cout << "  " << setw(10) << x
         << "  " << setw(10) << y
         << "  " << setw(10) << x-y
         << "  " << setw(10) << r8_diff ( x, y, ndig ) << "\n";
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
//    TEST006 tests R8_DIGIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define MAXDIG 20

  int idigit;
  double x;

  x = r8_pi ( );

  cout << "\n";
  cout << "TEST006\n";
  cout << "  R8_DIGIT extracts decimal digits.\n";
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
    cout << setw(3) << r8_digit ( x, idigit );
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
//    TEST007 tests R8_EPSILON and R8_EPSILON_COMPUTE.
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
{
  double r1;
  double r2;
  double s;
  double t;

  cout << "\n";
  cout << "TEST007\n";
  cout << "  R8_EPSILON produces the R8 roundoff unit.\n";
  cout << "  R8_EPSILON_COMPUTE computes the R8 roundoff unit.\n";
  cout << "\n";

  r1 = r8_epsilon ( );
  cout << "  R = R8_EPSILON()         = " << setprecision(16) << setw(24) << r1 << "\n";

  r2 = r8_epsilon_compute ( );
  cout << "  R = R8_EPSILON_COMPUTE() =  " << setprecision(16) << setw(24) << r2 << "\n";

  s = 1.0 + r2;
  t = s - 1.0;
  cout << "  ( 1 + R2 ) - 1           =     " << setprecision(16) << setw(24) << t << "\n";

  s = 1.0 + ( r2 / 2.0 );
  t = s - 1.0;
  cout << "  ( 1 + (R2/2) ) - 1       = " << setw(10) << t << "\n";

  return;
}
//****************************************************************************80

void test008 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST008 tests R8_FRACTIONAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  double fractional;
  double r8;
  double r8_hi = 5.0;
  double r8_lo = -3.0;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST008\n";
  cout << "  R8_FRACTIONAL returns the fractional part of an R8.\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    r8 = r8_uniform_ab ( r8_lo, r8_hi, seed );
    fractional = r8_fractional ( r8 );
    cout << "  " << setw(10) << r8
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
//    TEST009 tests R8_HUGE.
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
  cout << "  R8_HUGE returns a large R8 value;\n";
  cout << "\n";
  cout << "  R8_HUGE =   " << r8_huge ( ) << "\n";

  return;
}
//****************************************************************************80

void test010 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST010 tests R8_LOG_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 18

  int test;
  double x;
  double x_test[TEST_NUM] = {
    0.0,  1.0,  2.0,   3.0,  9.0,
   10.0, 11.0, 99.0, 101.0, -1.0,
   -2.0, -3.0, -9.0,   0.5,  0.33,
    0.25, 0.20, 0.01 };

  cout << "\n";
  cout << "TEST010\n";
  cout << "  R8_LOG_2: computes the logarithm base 2.\n";
  cout << "\n";
  cout << "  X       R8_LOG_2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setw(12) << x
         << "  " << setw(12) << r8_log_2 ( x ) << "\n";
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
//    TEST011 tests R8_LOG_B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 10

  double b;
  double b_test[TEST_NUM] = {
    2.0, 3.0, 4.0, 5.0, 6.0,
    7.0, 8.0, 16.0, 32.0, 256.0 };
  int test;
  double x;

  x = 16.0;

  cout << "\n";
  cout << "TEST011\n";
  cout << "  R8_LOG_B computes the logarithm base B.\n";
  cout << "\n";
  cout << "  X, B, R8_LOG_B\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = b_test[test];

    cout << "  " << setw(12) << x
         << "  " << setw(12) << b
         << "  " << setw(12) << r8_log_b ( x, b ) << "\n";
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
//    TEST012 tests R8_MANT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int is;
  int l;
  double r;
  double x;

  x = -314.159;

  cout << "\n";
  cout << "TEST012\n";
  cout << "  R8_MANT decomposes a value.\n";
  cout << "\n";
  cout << "  Number to be decomposed: X = " << x << "\n";

  r8_mant ( x, &is, &r, &l );

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
//    TEST013 tests R8_MOD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;
  double x_hi = 10.0;
  double x_lo = -10.0;
  double y;
  double z1;
  double z2;

  cout << "\n";
  cout << "TEST013\n";
  cout << "  R8_MOD returns the remainder after division.\n";
  cout << "  R8_MOD ( X, Y ) has the same sign as X.\n";
  cout << "\n";
  cout << "      X         Y    FMOD(X,Y)    R8_MOD(X,Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( x_lo, x_hi, seed );
    y = r8_uniform_ab ( x_lo, x_hi, seed );

    z1 =   fmod ( x, y );
    z2 = r8_mod ( x, y );

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
//    TEST014 tests R8_MODP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;
  double x_hi = 10.0;
  double x_lo = -10.0;
  double y;
  double z1;
  double z2;

  cout << "\n";
  cout << "TEST014\n";
  cout << "  R8_MODP returns the remainder after division.\n";
  cout << "  R8_MODP ( X, Y ) is positive.\n";
  cout << "\n";
  cout << "      X       Y     FMOD(X,Y)  R8_MODP(X,Y)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( x_lo, x_hi, seed );
    y = r8_uniform_ab ( x_lo, x_hi, seed );

    z1 =   fmod  ( x, y );
    z2 = r8_modp ( x, y );

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
//    TEST015 tests R8_NINT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double c;
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;

  cout << "\n";
  cout << "TEST015\n";
  cout << "  R8_NINT produces the nearest integer.\n";
  cout << "\n";

  b = -10.0;
  c = +10.0;

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( b, c, seed );
    cout << setw(10) << x << "  "
         << setw(6)  << r8_nint ( x ) << "\n";
  }

  return;
}
//****************************************************************************80

void test016 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST016 tests R8_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 20

  int seed = 123456789;
  int test;
  double x;

  cout << "\n";
  cout << "TEST016\n";
  cout << "  R8_NORMAL_01 generates normally distributed random values.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = r8_normal_01 ( seed );
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
//    TEST017 tests R8_PI.
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
  double four;
  double one;
  double v1;
  double v2;

  four = ( double ) ( 4 );
  one = ( double ) ( 1 );

  cout << "\n";
  cout << "TEST017\n";
  cout << "  R8_PI returns the value of PI.\n";
  cout << "\n";
  v1 = r8_pi ( );
  cout << "  R8_PI =     " << v1 << "\n";
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
//    TEST018 tests R8_POWER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int p;
  double r;
  double value;

  cout << "\n";
  cout << "TEST018\n";
  cout << "  R8_POWER computes R^P\n";
  cout << "\n";
  cout << "      R          P       R**P\n";
  cout << "\n";

  for ( p = -5; p <= 5; p++ )
  {
    r = 2.0;
    value = r8_power ( r, p );
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
//    TEST019 tests R8_POWER_FAST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int mults;
  int p;
  double r;
  double rp;

  cout << "\n";
  cout << "TEST019\n";
  cout << "  R8_POWER_FAST computes R^P, economizing on\n";
  cout << "    multiplications.\n";
  cout << "\n";
  cout << "      R          P       R**P       Mults\n";
  cout << "\n";

  for ( i = -10; i <= 40; i++ )
  {
    r = 2.0;
    p = i;
    rp = r8_power_fast ( r, p, &mults );
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
//    TEST020 tests R8_ROUND2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nplace;
  double x;
  double xround;

  x = r8_pi ( );

  cout << "\n";
  cout << "TEST020\n";
  cout << "  R8_ROUND2 rounds a number to a\n";
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
    xround = r8_round2 ( nplace, x );
    cout << "  " << setw(8) << i
         << "  " << setprecision(16) << setw(24) << xround << "\n";
  }

  return;
}
//****************************************************************************80

void test021 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST021 tests R8_ROUNDB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int base;
  int i;
  int nplace;
  double x;
  double xround;

  base = 3;
  x = r8_pi ( );

  cout << "\n";
  cout << "TEST021\n";
  cout << "  R8_ROUNDB rounds a number to a \n";
  cout << "  specified number of base IBASE digits.\n";
  cout << "\n";
  cout << "  Here, we will use IBASE = " << base << "\n";
  cout << "\n";
  cout << "  Test effect on PI:\n";
  cout << "  X = " << setprecision(16) << setw(24) << x << "\n";
  cout << "\n";
  cout << "  NPLACE  XROUND\n";
  cout << "\n";

  for ( i = 0; i <= 20; i++ )
  {
    nplace = i;
    xround = r8_roundb ( base, nplace, x );
    cout << "  " << setw(8) << i
         << "  " << setprecision(16) << setw(24) << xround << "\n";
  }

  cout << "\n";
  cout << "  Try with a negative base:\n";
  x = 121.0;
  base = -3;
  nplace = 3;
  cout << "\n";
  cout << "  Input quantity is X = " << x << "\n";
  cout << "  to be rounded in base " << base << "\n";

  for ( nplace = 1; nplace <= 5; nplace++ )
  {
    xround = r8_roundb ( base, nplace, x );

    cout << "\n";
    cout << "  Output value to " << nplace << " places is " << xround << "\n";
  }

  return;
}
//****************************************************************************80

void test022 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST022 tests R8_ROUNDX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nplace;
  int seed;
  double x;
  double xround;

  seed = 123456789;
  x = r8_pi ( );

  cout << "\n";
  cout << "TEST022\n";
  cout << "  R8_ROUNDX rounds a number to a \n";
  cout << "  specified number of decimal digits.\n";
  cout << "\n";
  cout << "  Test effect on PI:\n";
  cout << "  X = " << setprecision(16) << x << "\n";
  cout << "\n";
  cout << "  NPLACE  XROUND\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    nplace = i;
    xround = r8_roundx ( nplace, x );
    cout << "  " << setw(6) << i
         << "  " << setw(20) << xround << "\n";
  }

  cout << "\n";
  cout << "  Test effect on random values:\n";
  cout << "\n";
  cout << "  NPLACE  X     XROUND\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    x = r8_uniform_01 ( seed );

    cout << "\n";

    for ( nplace = 0; nplace <= 10; nplace = nplace + 2 )
    {
      xround = r8_roundx ( nplace, x );

      cout << "  " << setw(6)  << nplace
           << "  " << setw(16) << x
           << "  " << setw(20) << xround << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023 tests R8_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  int test;
  double x;
  double x_test[TEST_NUM] = { -1.25, -0.25, 0.0, +0.5, +9.0 };

  cout << "\n";
  cout << "TEST023\n";
  cout << "  R8_SIGN returns the sign of a number.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setprecision(6) << setw(8) << x
         << "  " << r8_sign ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0235 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0235 tests R8_SWAP.
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
  double x;
  double y;

  cout << "\n";
  cout << "TEST0235\n";
  cout << "  R8_SWAP swaps two reals.\n";

  x = 1.0;
  y = 3.14159;

  cout << "\n";
  cout << "  Before swapping: \n";
  cout << "\n";
  cout << "    X = " << x << "\n";
  cout << "    Y = " << y << "\n";

  r8_swap ( &x, &y );

  cout << "\n";
  cout << "  After swapping: \n";
  cout << "\n";
  cout << "    X = " << x << "\n";
  cout << "    Y = " << y << "\n";

  return;
}
//****************************************************************************80

void test024 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST024 tests R8_TO_R8_DISCRETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int ndx = 19;
  double r;
  double rd;
  double rhi = 10.0;
  double rhi2;
  double rlo = 1.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 15;

  cout << "\n";
  cout << "TEST024\n";
  cout << "  R8_TO_R8_DISCRETE maps numbers to a discrete set\n";
  cout << "  of equally spaced numbers in an interval.\n";
  cout << "\n";
  cout << "  Number of discrete values = " << ndx << "\n";
  cout << "  Real interval: " << rlo << "  " << rhi << "\n";
  cout << "\n";
  cout << "  R   RD\n";
  cout << "\n";

  seed = 123456789;

  rlo2 = rlo - 2.0;
  rhi2 = rhi + 2.0;

  for ( test = 0; test < test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, seed );
    rd = r8_to_r8_discrete ( r, rlo, rhi, ndx );
    cout << "  " << setw(14) << r
         << "  " << setw(14) << rd << "\n";
  }

  return;
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests R8_TO_I4.
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
  int ix;
  int ixmax;
  int ixmin;
  double x;
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  R8_TO_I4 finds an integer IX in [IXMIN,IXMAX]\n";
  cout << "  corresponding to X in [XMIN,XMAX].\n";

  xmin = 2.5;
  x = 3.5;
  xmax = 5.5;

  ixmin = 10;
  ixmax = 40;

  ix = r8_to_i4 ( x, xmin, xmax, ixmin, ixmax );

  cout << "\n";
  cout << "   XMIN " <<  xmin << "   X = " <<  x << "  XMAX = "
    <<  xmax << "\n";
  cout << "  IXMIN " << ixmin << "  IX = " << ix << " IXMAX = "
    << ixmax << "\n";

  return;
}
//****************************************************************************80

void test026 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026 tests R8_UNIFORM.
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
  double a;
  double b;
  double c;
  int i;
  int seed;

  b = 10.0;
  c = 25.0;
  seed = 17;

  cout << "\n";
  cout << "TEST026\n";
  cout << "  R8_UNIFORM produces a random real in a given range.\n";
  cout << "\n";
  cout << "  Using range " << b << " <= A <= " << c << ".\n";
  cout << "\n";

  cout << "\n";
  cout << "  I   A\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    a = r8_uniform_ab ( b, c, seed );
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
//    TEST027 tests R8_UNIFORM_01.
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
  double x;

  cout << "\n";
  cout << "TEST027\n";
  cout << "  R8_UNIFORM_01 produces a sequence of random values.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  cout << "\n";
  cout << "  SEED   R8_UNIFORM_01(SEED)\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << setw(12) << seed << "  ";
    x = r8_uniform_01 ( seed );
    cout << setw(10) << x << "\n";
  }

  cout << "\n";
  cout << "  Verify that the sequence can be restarted.\n";
  cout << "  Set the seed back to its original value, and see that\n";
  cout << "  we generate the same sequence.\n";

  seed = 123456789;
  cout << "\n";
  cout << "  SEED   R8_UNIFORM_01(SEED)\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    cout << setw(12) << seed << "  ";
    x = r8_uniform_01 ( seed );
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
//    TEST028 tests R8_UNIFORM_01.
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
  double max;
  double mean;
  double min;
  int n;
  int seed = 123456789;
  double x[N];
  double variance;

  cout << "\n";
  cout << "TEST028\n";
  cout << "  R8_UNIFORM_01 samples a uniform random distribution in [0,1].\n";
  cout << "  distributed random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = r8_uniform_01 ( seed );
  }

  cout << "\n";
  cout << "  First few values:\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  min = r8vec_min ( N, x );
  max = r8vec_max ( N, x );
  mean = r8vec_mean ( N, x );
  variance = r8vec_variance ( N, x );

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

void test029 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST029 tests R8_WALSH_1D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double w0;
  double wm1;
  double wm2;
  double wm3;
  double wp1;
  double wp2;
  double x;

  cout << "\n";
  cout << "TEST029\n";
  cout << "  R8_WALSH_1D evaluates 1D Walsh functions:\n";
  cout << "\n";
  cout << "  X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)\n";
  cout << "\n";

  for ( i = 0; i <= 32; i++ )
  {
    x = ( double ) ( i ) / 4.0;

    wp2 = r8_walsh_1d ( x,  2 );
    wp1 = r8_walsh_1d ( x,  1 );
    w0  = r8_walsh_1d ( x,  0 );
    wm1 = r8_walsh_1d ( x, -1 );
    wm2 = r8_walsh_1d ( x, -2 );
    wm3 = r8_walsh_1d ( x, -3 );

    cout << "  " << setw(10) << x
         << "  " << setw(2)  << wp2
         << "  " << setw(2)  << wp1
         << "  " << setw(2)  << w0
         << "  " << setw(2)  << wm1
         << "  " << setw(2)  << wm2
         << "  " << setw(2)  << wm3 << "\n";
  }

  return;
}
//****************************************************************************80

void test0295 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0295 tests R8_WRAP;
//
//  Discussion:
//
//    Apparently if you turn on high precision somewhere in COUT, you're
//    stuck with it forever...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a = - 2.0;
  double b = 12.0;
  double r;
  double r2;
  double rhi = 6.5;
  double rlo = 3.0;
  int seed;
  int test;
  int test_num = 20;

  cout << "\n";
  cout << "TEST0295\n";
  cout << "  R8_WRAP \"wraps\" an R8 to lie within an interval:\n";
  cout << "\n";
  cout << "  Wrapping interval is " << rlo << ", " << rhi << "\n";
  cout << "\n";
  cout << "      R      R8_WRAP ( R )\n";
  cout << "\n";
  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( a, b, seed );
    r2 = r8_wrap ( r, rlo, rhi );
    cout << "  " << setw(10) << r
         << "  " << setw(10) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void test031 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST031 tests R82POLY2_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 12

  double a;
  double a_test[TEST_NUM] = {
    9.0, 4.0, 9.0,  1.0, 0.0,
    1.0, 0.0, 0.0,  0.0, 0.0,
    0.0, 0.0 };
  double b;
  double b_test[TEST_NUM] = {
    -4.0, 1.0,  16.0,  1.0,  0.0,
     2.0, 1.0,   1.0,  1.0,  0.0,
     0.0, 0.0 };
  double c;
  double c_test[TEST_NUM] = {
     0.0, -4.0,   0.0,   0.0, 1.0,
     0.0,  0.0,   0.0,  0.0,  0.0,
     0.0,  0.0 };
  double d;
  double r8_test[TEST_NUM] = {
    -36.0,  3.0,  36.0,  -6.0, 3.0,
    -2.0,   0.0,   0.0,  0.0,  2.0,
     0.0, 0.0 };
  double e;
  double e_test[TEST_NUM] = {
    -24.0, -4.0, -32.0, -10.0, -1.0,
     16.0, -6.0, -6.0, -2.0, -1.0,
     0.0, 0.0 };
  double f;
  double f_test[TEST_NUM] = {
    -36.0,  1.0, -92.0, 115.0, -3.0,
     33.0, +8.0, 10.0,  +1.0,  1.0,
      0.0, 1.0 };
  int test;
  int type;

  cout << "\n";
  cout << "TEST031\n";
  cout << "  R82POLY2_TYPE determines the type of a second order\n";
  cout << "  equation in two variables.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];
    e = e_test[test];
    f = f_test[test];

    cout << "\n";

    r82poly2_print ( a, b, c, d, e, f );

    type = r82poly2_type ( a, b, c, d, e, f );

    cout << "  Type = " << type << "\n";

    r82poly2_type_print ( type );
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST032 tests R82VEC_ORDER_TYPE.
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
{
# define N 4
# define TEST_NUM 10

  int i;
  int j;
  int order;
  int seed = 123456789;
  int test;
  double *x;

  cout << "\n";
  cout << "TEST032\n";
  cout << "  R82VEC_ORDER_TYPE classifies an R8VEC as\n";
  cout << "  -1: no order\n";
  cout << "   0: all equal;\n";
  cout << "   1: ascending;\n";
  cout << "   2: strictly ascending;\n";
  cout << "   3: descending;\n";
  cout << "   4: strictly descending.\n";
  cout << "\n";

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    x = r8mat_uniform_01_new ( 2, N, seed );

    for ( j = 0; j < N; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        x[i+j*2] = ( double ) ( r8_nint ( 3.0 * x[i+j*2] ) );
      }
    }
    order = r82vec_order_type ( N, x );

    cout << "  Order type = " << order << "\n";

    r82vec_print ( N, x, " " );

    delete [] x;
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test033 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST033 tests R82VEC_PART_QUICK_A.
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
{
# define N 12

  double *a;
  double b = 0.0E+00;
  double c = 2.0E+00;
  int i;
  int l;
  int r;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST033\n";
  cout << "  R82VEC_PART_QUICK_A reorders an R82VEC\n";
  cout << "  as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_ab_new ( 2, N, b, c, seed );

  r82vec_print ( N, a, "  Before rearrangment:" );

  r82vec_part_quick_a ( N, a, &l, &r );

  cout << "\n";
  cout << "  Rearranged array\n";
  cout << "  Left index =  " << l << "\n";
  cout << "  Key index =   " << l+1 << "\n";
  cout << "  Right index = " << r << "\n";

  r82vec_print ( l,     a,         "  Left half:" );
  r82vec_print ( 1,     a+2*l,     "  Key:" );
  r82vec_print ( N-l-1, a+2*(l+1), "  Right half:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test034 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST034 tests R82VEC_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double *a;
  double b = 0.0;
  int base = 0;
  double c = 10.0;
  int i;
  int *indx;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST034\n";
  cout << "  R82VEC_SORT_HEAP_INDEX_A index sorts an R82VEC\n";
  cout << "  using heapsort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_ab_new ( 2, N, b, c, seed );
//
//  Give a few elements the same first component.
//
  a[0+2*2] = a[0+4*2];
  a[0+3*2] = a[0+11*2];
//
//  Give a few elements the same second component.
//
  a[1+5*2] = a[1+0*2];
  a[1+1*2] = a[1+8*2];
//
//  Make two entries equal.
//
  a[0+6*2] = a[0+10*2];
  a[1+6*2] = a[1+10*2];

  r82vec_print ( N, a, "  Before rearrangement:" );

  indx = r82vec_sort_heap_index_a ( N, base, a );

  cout << "\n";
  cout << "         I     Index   A(Index)\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << indx[i]
         << "  " << setw(12) << a[0+indx[i]*2]
         << "  " << setw(12) << a[1+indx[i]*2] << "\n";
  }

  r82vec_permute ( N, indx, base, a );

  r82vec_print ( N, a, "  After rearrangement by R82VEC_PERMUTE:" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void test035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST035 tests R82VEC_SORT_QUICK_A.
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
{
# define N 12

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST035\n";
  cout << "  R82VEC_SORT_QUICK_A sorts an R82VEC\n";
  cout << "  as part of a quick sort.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8mat_uniform_ab_new ( 2, N, b, c, seed );
//
//  For better testing, give a few elements the same first component.
//
  a[2*(3-1)+0] = a[2*(5-1)+0];
  a[2*(4-1)+0] = a[2*(12-1)+0];
//
//  Make two entries equal.
//
  a[2*(7-1)+0] = a[2*(11-1)+0];
  a[2*(7-1)+1] = a[2*(11-1)+1];

  r82vec_print ( N, a, "  Before sorting:" );

  r82vec_sort_quick_a ( N, a );

  r82vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test036 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST036 tests R8BLOCK_EXPAND_LINEAR.
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
# define L 4
# define M 3
# define N 2

  int l2;
  int lfat = 1;
  int m2;
  int mfat = 2;
  int n2;
  int nfat = 1;
  double x[L*M*N] = {
        1.0,  2.0,  3.0,   4.0,  1.0,
        4.0,  9.0, 16.0,   1.0,  8.0,
       27.0, 64.0,  2.0,   4.0,  6.0,
        8.0,  2.0,  8.0,  18.0, 32.0,
        2.0, 16.0, 54.0, 128.0 };
  double *xfat;

  l2 = ( L - 1 ) * ( lfat + 1 ) + 1;
  m2 = ( M - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( N - 1 ) * ( nfat + 1 ) + 1;

  cout << "\n";
  cout << "TEST036\n";
  cout << "  R8BLOCK_EXPAND_LINEAR linearly interpolates new data\n";
  cout << "  between old values in a 3D block.\n";

  r8block_print ( L, M, N, x, "  Original block:" );

  cout << "\n";
  cout << "  LFAT = " << lfat << "\n";
  cout << "  MFAT = " << mfat << "\n";
  cout << "  NFAT = " << nfat << "\n";

  xfat = r8block_expand_linear ( L, M, N, x, lfat, mfat, nfat );

  r8block_print ( l2, m2, n2, xfat, "  Fattened block:" );

  delete [] xfat;

  return;
# undef L
# undef M
# undef N
}
//****************************************************************************80

void test0363 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0363 tests R8BLOCK_NEW and R8BLOCK_DELETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double ***a;
  double ***b;
  int i;
  int j;
  int k;
  int l;
  int m;
  int n;

  cout << "\n";
  cout << "TEST0363:\n";
  cout << "  R8BLOCK_NEW dynamically creates a 3D array.\n";
  cout << "  R8BLOCK_DELETE deletes it.\n";
  cout << "  Array entries can be addressed using the\n";
  cout << "  notation \"a[i][j][k]\".\n";
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by R8BLOCK_DELETE and new memory
//  requested by R8BLOCK_NEW.
//
  l = 2;
  m = 3;
  n = 2;
//
//  Allocate memory.
//
  cout << "\n";
  cout << "  Allocating memory for array A of size " << l << " by " << m << " by " << n << ".\n";

  a = r8block_new ( l, m, n );

  cout << "\n";
  cout << "  Assigning values to A.\n";
//
//  Store values in A.
//
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        a[i][j][k] = ( double ) ( 100 * i + 10 * j + k );
      }
    }
  }
//
//  Print A.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix A:\n";
  cout << "\n";
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        cout << "  " << setw(8) << a[i][j][k];
      }
      cout << "\n";
    }
    cout << "\n";
  }
//
//  Free memory.
//
  r8block_delete ( a, l, m, n );

  return;
}
//****************************************************************************80

void test0365 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0365 tests R8BLOCK_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2012
//
//  Author:
//
//    John Burkardt
//
{
# define L 4
# define M 3
# define N 2

  double x[L*M*N] = {
        1.0,  2.0,  3.0,   4.0,  1.0,
        4.0,  9.0, 16.0,   1.0,  8.0,
       27.0, 64.0,  2.0,   4.0,  6.0,
        8.0,  2.0,  8.0,  18.0, 32.0,
        2.0, 16.0, 54.0, 128.0 };

  cout << "\n";
  cout << "TEST0365\n";
  cout << "  R8BLOCK_PRINT prints an R8BLOCK.\n";

  r8block_print ( L, M, N, x, "  The 3D array:" );

  return;
# undef L
# undef M
# undef N
}
//****************************************************************************80

void test037 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST037 tests R8COL_FIND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  int col;
  double dtab[M*N];
  double r8vec[M];
  int i;
  int j;
  int k;

  k = 1;

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      dtab[i+j*M] = ( double ) k;
      if ( j == 2 )
      {
        r8vec[i] = ( double ) k;
      }
      k = k + 1;
    }
  }

  col = r8col_find ( M, N, dtab, r8vec );

  cout << "\n";
  cout << "TEST037\n";
  cout << "  R8COL_FIND finds a column in a table matching\n";
  cout << "  a given set of data.\n";
  cout << "\n";
  cout << "  R8COL_FIND returns COL = " << col << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void test038 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST038 tests R8COL_INSERT and R8COL_SORT_HEAP_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N_MAX 10

  double a[M*N_MAX] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0 };
  int col;
  double r8vec1[M] = { 3.0, 7.0, 11.0 };
  double r8vec2[M] = { 3.0, 4.0, 18.0 };
  int n;

  cout << "\n";
  cout << "TEST038\n";
  cout << "  R8COL_SORT_HEAP_A ascending heap sorts a table of columns.\n";
  cout << "  R8COL_INSERT inserts new columns.\n";

  n = 4;

  r8mat_print ( M, n, a, "  The unsorted matrix:" );

  r8col_sort_heap_a ( M, n, a );

  r8mat_print ( M, n, a, "  The sorted matrix:" );

  r8vec_print ( M, r8vec1, "  New column:" );

  col = r8col_insert ( N_MAX, M, n, a, r8vec1 );

  if ( col < 0 )
  {
    cout << "\n";
    cout << "  The data was already in column " << abs ( col ) << "\n";
  }
  else
  {
    r8mat_print ( M, n, a, "  The updated matrix:" );
  }

  r8vec_print ( M, r8vec2, "  New column:" );

  col = r8col_insert ( N_MAX, M, n, a, r8vec2 );

  if ( col < 0 )
  {
    cout << "\n";
    cout << "  The data was already in column " << abs ( col ) << "\n";
  }
  else
  {
    r8mat_print ( M, n, a, "  The updated matrix:" );
  }

  return;
# undef M
# undef N_MAX
}
//****************************************************************************80

void test0385 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0385 tests R8COL_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 15

  double a[M*N] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    3.0,  4.0, 18.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int base = 0;
  int i;
  int *indx;
  int j;
  int j2;
  int m = M;
  int n = N;

  cout << "\n";
  cout << "TEST0385\n";
  cout << "  R8COL_SORT_HEAP_INDEX_A computes an index vector which\n";
  cout << "  ascending sorts an R8COL.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  indx = r8col_sort_heap_index_a ( m, n, base, a );

  cout << "\n";
  cout << "  The implicitly sorted R8COL (transposed)\n";
  cout << "\n";

  for ( j = 0; j < n; j++ )
  {
    j2 = indx[j];
    cout << "  " << setw(4) << j2 << ":";
    for ( i = 0; i < m; i++ )
    {
      cout << "  " << setw(10) << a[i+j2*m];
    }
    cout << "\n";
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test039 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST039 tests R8COL_SORT_QUICK_A.
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
{
# define M 3
# define N 10

  double *a;
  double b = 0.0;
  double c = 10.0;
  int seed;

  cout << "\n";
  cout << "TEST039\n";
  cout << "  R8COL_SORT_QUICK_A sorts a table of columns.\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  The unsorted matrix:" );

  r8col_sort_quick_a ( M, N, a );

  r8mat_print ( M, N, a, "  The sorted matrix:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0393 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0393 tests R8COL_SORTED_TOL_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  cout << "\n";
  cout << "TEST0393\n";
  cout << "  R8COL_SORTED_TOL_UNIQUE finds tolerably unique columns \n";
  cout << "  in a sorted R8COL.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  unique_num = r8col_sorted_tol_unique ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << unique_num << "\n";

  r8mat_transpose_print ( m, unique_num, a,
    "  The sorted tolerably unique R8COL (transposed):" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0395 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0395 tests R8COL_SORTED_UNIQUE_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  cout << "\n";
  cout << "TEST0395\n";
  cout << "  R8COL_SORTED_UNIQUE_COUNT counts tolerably unique columns \n";
  cout << "  in a sorted R8COL.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  unique_num = r8col_sorted_tol_unique_count ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << unique_num << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0397 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0397 tests R8COL_SORTED_TOL_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << "\n";
  cout << "TEST0397\n";
  cout << "  R8COL_SORTED_TOL_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the tolerably unique columns of a sorted R8COL,\n";
  cout << "  and a map from the original R8COL to the (implicit)\n";
  cout << "  R8COL of sorted tolerably unique elements.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  n_unique = r8col_sorted_tol_unique_count ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << n_unique << "\n";

  au = new double[m*n_unique];
  undx = new int[n_unique];
  xdnu = new int[n];

  r8col_sorted_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  cout << "\n";
  cout << "  XDNU points to the representative for each item.\n";
  cout << "  UNDX selects the representatives.\n";
  cout << "\n";
  cout << "     I  XDNU  UNDX\n";
  cout << "\n";
  for ( i = 0; i < n_unique; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << undx[i] << "\n";
  }
  for ( i = n_unique; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i] << "\n";
  }
  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au,
    "  The tolerably unique R8COL (transposed):" );

  delete [] au;
  delete [] undx;
  delete [] xdnu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test040 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST040 tests R8COL_MAX and R8COL_MIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *amax;
  double *amin;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST040\n";
  cout << "  R8COL_MAX computes maximums of an R8COL;\n";
  cout << "  R8COL_MIN computes minimums of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  amax = r8col_max ( M, N, a );

  amin = r8col_min ( M, N, a );

  cout << "\n";
  cout << "  Column, maximum, minimum:\n";
  cout << "\n";

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(3) << j+1
         << "  " << setw(10) << amax[j]
         << "  " << setw(10) << amin[j] << "\n";
  }

  delete [] amax;
  delete [] amin;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test041 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST041 tests R8COL_MEAN and R8COL_SUM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *colsum;
  int i;
  int j;
  int k;
  double *mean;

  cout << "\n";
  cout << "TEST041\n";
  cout << "  R8COL_MEAN computes means of an R8COL;\n";
  cout << "  R8COL_SUM computes sums of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  colsum = r8col_sum ( M, N, a );

  mean = r8col_mean ( M, N, a );

  cout << "\n";
  cout << "  Column  sum, mean:\n";
  cout << "\n";

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(3) << j+1
         << "  " << setw(10) << colsum[j]
         << "  " << setw(10) << mean[j] << "\n";
  }

  delete [] mean;
  delete [] colsum;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0415 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0415 tests R8VEC_PERMUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 5

  double a[M*N] = {
    11.0, 21.0, 31.0,
    12.0, 22.0, 32.0,
    13.0, 23.0, 33.0,
    14.0, 24.0, 34.0,
    15.0, 25.0, 35.0 };
  int base = 1;
  int perm[N] = { 2, 4, 5, 1, 3 };


  cout << "\n";
  cout << "TEST0415\n";
  cout << "  R8COL_PERMUTE permutes an R8COL in place.\n";

  r8mat_print ( M, N, a, "  A (unpermuted):" );

  i4vec_print ( N, perm, "  The (column) permutation vector:" );

  r8col_permute ( M, N, perm, base, a );

  r8mat_print ( M, N, a, "  A (permuted):" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test042 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST042 tests R8COL_SORTR_A.
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
{
# define M 10
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int key;
  int seed;

  cout << "\n";
  cout << "TEST042\n";
  cout << "  R8COL_SORTR_A is given an array, and reorders\n";
  cout << "  it so that a particular column is sorted.\n";

  key = 2;
  cout << "\n";
  cout << "  Here, the special column is " << key << "\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  Unsorted array:" );

  r8col_sortr_a ( M, N, a, key );

  r8mat_print ( M, N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test043 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST043 tests R8COL_SWAP;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int icol1;
  int icol2;
  int j;
  int k;

  cout << "\n";
  cout << "TEST043\n";
  cout << "  R8COL_SWAP swaps two columns of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) ( k );
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  icol1 = 1;
  icol2 = 3;

  cout << "\n";
  cout << "  Swap columns " << icol1 << " and " << icol2 << ":\n";

  r8col_swap ( M, N, a, icol1, icol2 );

  r8mat_print ( M, N, a, "  The updated matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test044 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST044 tests R8COL_TO_R8VEC.
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
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  double *x;

  cout << "\n";
  cout << "TEST044\n";
  cout << "  R8COL_TO_R8VEC converts an array of columns to a vector.\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = ( double ) ( 10 * i + j );
    }
  }

  r8mat_print ( M, N, a, "  The array of columns:" );

  x = r8col_to_r8vec ( M, N, a );

  r8vec_print ( M*N, x, "  The resulting vector of columns:" );

  delete [] x;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0442 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0442 tests R8COL_TOL_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << "\n";
  cout << "TEST0442\n";
  cout << "  R8COL_TOL_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the tolerably unique columns of an R8COL,\n";
  cout << "  and a map from the original R8COL to the (implicit)\n";
  cout << "  R8COL of sorted tolerably unique elements.\n";

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  tol = 0.25;

  cout << "\n";
  cout << "  Using tolerance = " << tol << "\n";

  n_unique = r8col_tol_unique_count ( m, n, a, tol );

  cout << "\n";
  cout << "  Number of tolerably unique columns is " << n_unique << "\n";

  au = new double[m*n_unique];
  undx = new int[n_unique];
  xdnu = new int[n];

  r8col_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  cout << "\n";
  cout << "  XDNU points to the representative for each item.\n";
  cout << "  UNDX selects the representatives.\n";
  cout << "\n";
  cout << "     I  XDNU  UNDX\n";
  cout << "\n";
  for ( i = 0; i < n_unique; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << undx[i] << "\n";
  }
  for ( i = n_unique; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i] << "\n";
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au,
    "  The tolerably unique R8COL (transposed):" );

  delete [] au;
  delete [] undx;
  delete [] xdnu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0443 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0443 tests R8COL_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << "\n";
  cout << "TEST0443\n";
  cout << "  R8COL_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique columns of an (unsorted) R8COL,\n";
  cout << "  and a map from the original R8COL to the (implicit)\n";
  cout << "  R8COL of sorted unique elements.\n";

  r8mat_transpose_print ( m, n, a, "  The R8COL (transposed):" );

  n_unique = r8col_unique_count ( m, n, a );

  cout << "\n";
  cout << "  Number of unique columns is " << n_unique << "\n";

  au = new double[m*n_unique];
  undx = new int[n_unique];
  xdnu = new int[n];

  r8col_undex ( m, n, a, n_unique, undx, xdnu );

  cout << "\n";
  cout << "  XDNU points to the representative for each item.\n";
  cout << "  UNDX selects the representatives.\n";
  cout << "\n";
  cout << "     I  XDNU  UNDX\n";
  cout << "\n";
  for ( i = 0; i < n_unique; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << undx[i] << "\n";
  }
  for ( i = n_unique; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i] << "\n";
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au, "  The Unique R8COL (transposed):" );

  delete [] au;
  delete [] undx;
  delete [] xdnu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0445 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0445 tests R8COL_UNIQUE_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  cout << "\n";
  cout << "TEST0445\n";
  cout << "  R8COL_UNIQUE_COUNT counts unique columns.\n";

  r8mat_transpose_print ( m, n, a, "  The R8COL (transposed):" );

  unique_num = r8col_unique_count ( m, n, a );

  cout << "\n";
  cout << "  Number of unique columns is " << unique_num << "\n";

  return;
# undef M
# undef N
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests R8COL_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *variance;

  cout << "\n";
  cout << "TEST045\n";
  cout << "  R8COL_VARIANCE computes variances of an R8COL;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) ( k );
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  variance = r8col_variance ( M, N, a );

  cout << "\n";
  cout << "  Column  variance:\n";
  cout << "\n";

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(6)  << j
         << "  " << setw(10) << variance[j] << "\n";
  }

  delete [] variance;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test046 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST046 tests R8R8VEC_INDEX_INSERT_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 30

  int i;
  int ierror;
  int indx[N_MAX];
  int ival;
  int n;
  int seed;
  double x[N_MAX];
  double x_max = 4.0;
  double x_min = 1.0;
  double xval;
  double y[N_MAX];
  double y_max = 3.0;
  double y_min = 1.0;
  double yval;

  n = 0;

  cout << "\n";
  cout << "TEST046\n";
  cout << "  R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n";
  cout << "  index sorted array.\n";
  cout << "\n";
  cout << "  Generate " << N_MAX << " random values:\n";
  cout << "\n";
  cout << "    XVAL    YVAL   Index\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( x_min, x_max, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    yval = r8_uniform_ab ( y_min, y_max, seed );
    yval = ( double ) ( r8_nint ( yval ) );

    r8r8vec_index_insert_unique ( N_MAX, &n, x, y, indx, xval, yval,
      &ival, &ierror );

    cout << "  " << setw(6)  << ival
         << "  " << setw(12) << xval
         << "  " << setw(12) << yval << "\n";
  }

  cout << "\n";
  cout << "  Vector of unique X Y values:\n";
  cout << "\n";
  cout << "  I  X(I)   Y(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6)  << i+1
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  cout << "\n";
  cout << "  X, Y sorted by index\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(INDX(I))  Y(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i+1
         << "  " << setw(6) << indx[i]
         << "  " << setw(12) << x[indx[i]-1]
         << "  " << setw(12) << y[indx[i]-1] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test047 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST047 tests R8R8R8VEC_INDEX_INSERT_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 30

  int i;
  int ierror;
  int indx[N_MAX];
  int ival;
  int n;
  int seed;
  double x[N_MAX];
  double xval;
  double y[N_MAX];
  double yval;
  double z[N_MAX];
  double zval;

  n = 0;

  cout << "\n";
  cout << "TEST047\n";
  cout << "  R8R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values into\n";
  cout << "  an index sorted array.\n";
  cout << "\n";
  cout << "  Number of random values to generate = " << N_MAX << "\n";
  cout << "\n";
  cout << "    XVAL    YVAL  ZVAL  Index\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( 1.0, 4.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    yval = r8_uniform_ab ( 1.0, 3.0, seed );
    yval = ( double ) ( r8_nint ( yval ) );
    zval = r8_uniform_ab ( 1.0, 4.0, seed );
    zval = ( double ) ( r8_nint ( zval ) );

    r8r8r8vec_index_insert_unique ( N_MAX, &n, x, y, z, indx,
      xval, yval, zval, &ival, &ierror );

    cout << "  " << setw(6) << xval
         << "  " << setw(6) << yval
         << "  " << setw(6) << zval
         << "  " << setw(6) << ival << "\n";
  }

  cout << "\n";
  cout << "  Vector of unique X Y Z values:\n";
  cout << "\n";
  cout << "  I  X(I)   Y(I)    Z(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << y[i]
         << "  " << setw(6) << z[i] << "\n";
  }

  cout << "\n";
  cout << "  X Y Z sorted by index:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[indx[i]-1]
         << "  " << setw(6) << y[indx[i]-1]
         << "  " << setw(6) << z[indx[i]-1] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test048 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST048 tests R8INT_TO_I4INT and I4INT_TO_R8INT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ihi = 11;
  int ilo = 1;
  int ir;
  double r;
  double r2;
  double rhi = 200.0;
  double rhi2;
  double rlo = 100.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST048\n";
  cout << "  For data in an interval,\n";
  cout << "  I4INT_TO_R8INT converts an integer to a real;\n";
  cout << "  R8INT_TO_I4INT converts a real to an integer.\n";
  cout << "\n";
  cout << "  Integer interval: [" << ilo << ", " << ihi << "]\n";
  cout << "  Real interval:    [" << rlo << ", " << rhi << "]\n";
  cout << "\n";
  cout << "  R   I(R)  R(I(R))\n";
  cout << "\n";

  seed = 123456789;

  rlo2 = rlo - 15.0;
  rhi2 = rhi + 15.0;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, seed );
    ir = r8int_to_i4int ( rlo, rhi, r, ilo, ihi );
    r2 = i4int_to_r8int ( ilo, ihi, ir, rlo, rhi );
    cout << "  " << setw(12) << r
         << "  " << setw(6)  << ir
         << "  " << setw(12) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void test049 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST049 tests R8MAT_CHOLESKY_FACTOR, R8MAT_CHORESKY_FACTOR and R8MAT_CHOLESKY_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N*N];
  double b[N];
  int flag;
  double *d;
  int i;
  int j;
  double *l;
  double *lt;
  double *r;
  double *rt;
  double *x;

  cout << "\n";
  cout << "TEST049\n";
  cout << "  For a positive definite symmetric matrix,\n";
  cout << "  R8MAT_CHOLESKY_FACTOR computes the lower\n";
  cout << "  triangular Cholesky factor;\n";
  cout << "  R8MAT_CHORESKY_FACTOR computes the upper\n";
  cout << "  triangular Cholesky factor;\n";
  cout << "  R8MAT_CHOLESKY_SOLVE solves a linear system\n";
  cout << "  using the Cholesky factorization.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i == j )
      {
        a[i+j*N] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a[i+j*N] = -1.0;
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  r8mat_print ( N, N, a, "  Matrix to be factored:" );
//
//  Compute L, the lower Cholesky factor.
//
  l = r8mat_cholesky_factor ( N, a, flag );

  if ( flag != 0 )
  {
    cout << "\n";
    cout << "  R8MAT_CHOLESKY_FACTOR failed.\n";
    return;
  }

  r8mat_print ( N, N, l, "  Cholesky factor L:" );

  lt = r8mat_transpose_new ( N, N, l );

  d = r8mat_mm_new ( N, N, N, l, lt );

  r8mat_print ( N, N, d, "  Product L * L':" );
//
//  Compute R, the upper Cholesky factor.
//
  r = r8mat_choresky_factor ( N, a, flag );

  if ( flag != 0 )
  {
    cout << "\n";
    cout << "  R8MAT_CHORESKY_FACTOR failed.\n";
    return;
  }

  r8mat_print ( N, N, r, "  Cholesky factor R:" );

  rt = r8mat_transpose_new ( N, N, r );

  d = r8mat_mm_new ( N, N, N, r, rt );

  r8mat_print ( N, N, d, "  Product R * R':" );
//
//  Solve a system.
//
  for ( i = 0; i < N-1; i++ )
  {
    b[i] = 0.0;
  }
  b[N-1] = ( double ) ( N + 1 );

  r8vec_print ( N, b, "  Right hand side:" );

  x = r8mat_cholesky_solve ( N, l, b );

  r8vec_print ( N, x, "  Computed solution:" );

  delete [] d;
  delete [] l;
  delete [] lt;
  delete [] r;
  delete [] rt;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test050 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST050 tests R8MAT_DET_2D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0 };

  cout << "\n";
  cout << "TEST050\n";
  cout << "  R8MAT_DET_2D: determinant of a 2 by 2 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_2d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_2D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test051 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST051 tests R8MAT_DET_3D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0 };

  cout << "\n";
  cout << "TEST051\n";
  cout << "  R8MAT_DET_3D: determinant of a 3 by 3 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_3d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_3D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test052 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST052 tests R8MAT_DET_4D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0, 2.0 };

  cout << "\n";
  cout << "TEST052\n";
  cout << "  R8MAT_DET_4D determinant of a 4 by 4 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_4d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_4D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test053 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST053 tests R8MAT_DET_5D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0, 2.0, 3.0 };

  cout << "\n";
  cout << "TEST053\n";
  cout << "  R8MAT_DET_5D determinant of a 5 by 5 matrix;\n";

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_5d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  cout << "\n";
  cout << "  R8MAT_DET_5D computes determinant:" << det << "\n";
//
//  Special formula for the determinant of a Vandermonde matrix:
//
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  cout << "  Exact determinant is " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test054 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST054 tests R8MAT_EXPAND_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 3

  int m2;
  int mfat = 2;
  int n2;
  int nfat = 1;
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double x[M*N] = {
    1.0, 2.0, 3.0, 4.0, 1.0,
    4.0, 9.0, 16.0, 1.0, 8.0,
    27.0, 64.0 };
  double *xfat;

  m2 = ( M - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( N - 1 ) * ( nfat + 1 ) + 1;

  cout << "\n";
  cout << "TEST054\n";
  cout << "  R8MAT_EXPAND_LINEAR linearly interpolates new data\n";
  cout << "  between old values in a matrix.\n";

  r8mat_print ( M, N, x, "  Original matrix:" );

  cout << "\n";
  cout << "  MFAT = " << mfat << "\n";
  cout << "  NFAT = " << nfat << "\n";

  xfat = r8mat_expand_linear ( M, N, x, mfat, nfat );

  r8mat_print ( m2, n2, xfat, "  Fattened matrix:" );

  delete [] xfat;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test055 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST055 tests R8MAT_EXPAND_LINEAR2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 2

  double a[M*N];
  double *a2;
  int i;
  int j;
  int m2 = 10;
  int n2 = 5;

  cout << "\n";
  cout << "TEST055\n";
  cout << "  R8MAT_EXPAND_LINEAR2 fills in a large array by\n";
  cout << "  interpolating data from a small array.\n";
  cout << "\n";
  cout << "  Original matrix has dimensions:\n";
  cout << "\n";
  cout << "  M = " << M << ", N = " << N << "\n";
  cout << "\n";
  cout << "  Expanded matrix has dimensions:\n";
  cout << "\n";
  cout << "  M2 = " << m2 << ", N2 = " << n2 << "\n";

  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      a[i-1+(j-1)*M] = 10.0 * ( double ) ( i ) + ( double ) ( j );
    }
  }

  r8mat_print ( M, N, a, "  The little matrix A:" );

  a2 = r8mat_expand_linear2 ( M, N, a, m2, n2 );

  r8mat_print ( m2, n2, a2, "  Expanded array A2:" );

  delete [] a2;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0555 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0555 tests R8MAT_FSS_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 November 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define NB 3

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int k;
  int n = N;
  int nb = NB;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST0555\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8MAT_FSS_NEW factors and solves multiple linear systems.\n";
  cout << "\n";
  cout << "  Matrix order N = " << n << "\n";
  cout << "  Number of systems NB = " << nb << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( n, n, seed );
//
//  Set the desired solutions.
//
  b = new double[n * nb];

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0;
  }
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  k = 1;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( i % 3 ) + 1;
  }
  k = 2;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
//
//  Factor and solve the system.
//
  delete [] x;

  x = r8mat_fss_new ( n, a, nb, b );
  
  r8mat_print ( n, nb, x, "  Solutions:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
# undef NB
}
//****************************************************************************80

void test056 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST056 tests R8MAT_GIVENS_POST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N*N];
  double *ag;
  int col;
  double *g;
  int i;
  int j;
  int row;

  cout << "\n";
  cout << "TEST056\n";
  cout << "  R8MAT_GIVENS_POST computes a Givens postmultiplier rotation matrix.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_power ( i + 1, j );
    }
  }

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 3;
  col = 2;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  g = r8mat_givens_post ( N, a, row, col );

  r8mat_print ( N, N, g, "  G" );

  ag = r8mat_mm_new ( N, N, N, a, g );

  r8mat_print ( N, N, ag, "  A*G" );

  delete [] ag;
  delete [] g;

  return;
# undef N
}
//****************************************************************************80

void test057 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST057 tests R8MAT_GIVENS_PRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N*N];
  int col;
  double *g;
  double *ga;
  int i;
  int j;
  int row;

  cout << "\n";
  cout << "TEST057\n";
  cout << "  R8MAT_GIVENS_PRE computes a Givens premultiplier rotation matrix.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_power ( i + 1, j );
    }
  }

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 3;
  col = 2;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  g = r8mat_givens_pre ( N, a, row, col );

  r8mat_print ( N, N, g, "  G" );

  ga = r8mat_mm_new ( N, N, N, g, a );

  r8mat_print ( N, N, ga, "  G*A" );

  delete [] g;
  delete [] ga;

  return;
# undef N
}
//****************************************************************************80

void test058 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST058 tests R8MAT_HESS.
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
# define N 3

  double *h;
  double x[N] = { 1.0, 2.0, 3.0 };

  cout << "\n";
  cout << "TEST058\n";
  cout << "  R8MAT_HESS estimates the Hessian matrix\n";
  cout << "  of a scalar function.\n";

  h = r8mat_hess ( test058_f, N, x );

  r8mat_print ( N, N, h, "  Estimated jacobian:" );

  delete [] h;

  h = test058_hess ( N, x );

  r8mat_print ( N, N, h, "  Exact jacobian:" );

  delete [] h;

  return;
# undef N
}
//****************************************************************************80

double test058_f ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST058_F is a sample nonlinear function for treatment by R8MAT_JAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double TEST058_F, the function value.
//
{
  double f;

  f = x[0] * x[0] + x[0] * x[1] + x[1] * cos ( 10.0 * x[2] );

  return f;
}
//****************************************************************************80

double *test058_hess ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST058_HESS is the exact Hessian of TEST058_F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double TEST058_H[N*N], the Hessian values.
//
{
  double *h;

  h = new double[n*n];

  h[0+0*3] = 2.0;
  h[0+1*3] = 1.0;
  h[0+2*3] = 0.0;

  h[1+0*3] = 1.0;
  h[1+1*3] = 0.0;
  h[1+2*3] = -10.0 * sin ( 10.0 * x[2] );

  h[2+0*3] = 0.0;
  h[2+1*3] = -10.0 * sin ( 10.0 * x[2] );
  h[2+2*3] = -100.0 * x[1] * cos ( 10.0 * x[2] );

  return h;
}
//****************************************************************************80

void test059 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST059 tests R8MAT_HOUSE_FORM and R8VEC_HOUSE_COLUMN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 5.0;
  double *h;
  double *ha;
  int i;
  int j;
  int k;
  int n = 4;
  int seed;
  double *v;

  cout << "\n";
  cout << "TEST059\n";
  cout << "  R8VEC_HOUSE_COLUMN returns the compact form of\n";
  cout << "  a Householder matrix that packs a column\n";
  cout << "  of a matrix.\n";
//
//  Get a random matrix.
//
  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, b, c, seed );

  r8mat_print ( n, n, a, "  Matrix A:" );

  for ( k = 1; k <= n-1; k++ )
  {
    cout << "\n";
    cout << "  Working on column K = " << k << "\n";

    v = r8vec_house_column ( n, a+(k-1)*n, k );

    h = r8mat_house_form ( n, v );

    r8mat_print ( n, n, h, "  Householder matrix H:" );

    ha = r8mat_mm_new ( n, n, n, h, a );

    r8mat_print ( n, n, ha, "  Product H*A:" );
//
//  If we set A := HA, then we can successively convert A to upper
//  triangular form.
//
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        a[i+j*n] = ha[i+j*n];
      }
    }

    delete [] h;
    delete [] ha;
    delete [] v;
  }

  delete [] a;

  return;
}
//****************************************************************************80

void test060 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST060 tests R8MAT_HOUSE_FORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *h;
  double v[N] = { 0.0, 0.0, 1.0, 2.0, 3.0 };

  cout << "\n";
  cout << "TEST060\n";
  cout << "  R8MAT_HOUSE_FORM forms a Householder\n";
  cout << "  matrix from its compact form.\n";

  r8vec_print ( N, v, "  Compact vector form V:" ) ;

  h = r8mat_house_form ( N, v );

  r8mat_print ( N, N, h, "  Householder matrix H:" );

  delete [] h;

  return;
# undef N
}
//****************************************************************************80

void test061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST061 tests R8MAT_HOUSE_POST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *ah;
  double b = 0.0;
  double c = 5.0;
  int col;
  double *h;
  int n = 5;
  int row;
  int seed;

  cout << "\n";
  cout << "TEST061\n";
  cout << "  R8MAT_HOUSE_POST computes a Householder postmultiplier;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, b, c, seed );

  r8mat_print ( n, n, a, "  Matrix A:" );

  row = 2;
  col = 3;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  h = r8mat_house_post ( n, a, row, col );

  r8mat_print ( n, n, h, "  Householder matrix H:" );

  ah = r8mat_mm_new ( n, n, n, a, h );

  r8mat_print ( n, n, ah, "  Product A*H:" );

  delete [] a;
  delete [] ah;
  delete [] h;

  return;
}
//****************************************************************************80

void test062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST062 tests R8MAT_HOUSE_PRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double b = 0.0;
  double c = 5.0;
  int col;
  double *h;
  double *ha;
  int row;
  int seed;

  cout << "\n";
  cout << "TEST062\n";
  cout << "  R8MAT_HOUSE_PRE computes a Householder premultiplier;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( N, N, b, c, seed );

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 2;
  col = 3;

  cout << "\n";
  cout << "  I = " << row << "  J = " << col << "\n";

  h = r8mat_house_pre ( N, a, row, col );

  r8mat_print ( N, N, h, "  Householder matrix H:" );

  ha = r8mat_mm_new ( N, N, N, h, a );

  r8mat_print ( N, N, ha, "  Product H*A:" );

  delete [] a;
  delete [] h;
  delete [] ha;

  return;
# undef N
}
//****************************************************************************80

void test063 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST063 tests R8MAT_MAX_INDEX and R8MAT_MIN_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST063\n";
  cout << "  R8MAT_MAX_INDEX locates the maximum entry of an R8MAT;\n";
  cout << "  R8MAT_MIN_INDEX locates the minimum entry of an R8MAT;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );

  r8mat_print ( M, N, a, "  Random array:" );

  r8mat_max_index ( M, N, a, &i, &j );

  cout << "\n";
  cout << "  Maximum I,J indices            " << i << "  " << j << "\n";
  r8mat_min_index ( M, N, a, &i, &j );
  cout << "  Minimum I,J indices            " << i << "  " << j << "\n";

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test064 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST064 tests R8MAT_INVERSE_2D.
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
# define N 2

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST064\n";
  cout << "  R8MAT_INVERSE_2D inverts a 2 by 2 matrix.\n";

  a[0+0*N] = 1.0;
  a[0+1*N] = 2.0;

  a[1+0*N] = 3.0;
  a[1+1*N] = 4.0;

  r8mat_print ( 2, 2, a, "  Matrix A:" );

  b = r8mat_inverse_2d ( a );

  r8mat_print ( 2, 2, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 2, 2, c, "  C = A * B:" );

  delete [] b;

  return;

# undef N
}
//****************************************************************************80

void test065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 tests R8MAT_INVERSE_3D.
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
# define N 3

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST065\n";
  cout << "  R8MAT_INVERSE_3D inverts a 3 by 3 matrix.\n";

  a[0+0*N] = 3.0;
  a[0+1*N] = 2.0;
  a[0+2*N] = 1.0;

  a[1+0*N] = 2.0;
  a[1+1*N] = 2.0;
  a[1+2*N] = 1.0;

  a[2+0*N] = 0.0;
  a[2+1*N] = 1.0;
  a[2+2*N] = 1.0;

  r8mat_print ( 3, 3, a, "  Matrix A:" );

  b = r8mat_inverse_3d ( a );

  r8mat_print ( 3, 3, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 3, 3, c, "  C = A * B:" );

  delete [] b;

  return;

# undef N
}
//****************************************************************************80

void test066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST066 tests R8MAT_INVERSE_4D.
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
# define N 4

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST066\n";
  cout << "  R8MAT_INVERSE_4D inverts a 4 x 4 matrix.\n";


  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i <= j )
      {
        a[i+j*N] = ( double ) ( N - j );
      }
      else if ( j == i - 1 )
      {
        a[i+j*N] = ( double ) ( N - j - 1 );
      }
      else {
        a[i+j*N] = 0.0;
      }
    }
  }

  r8mat_print ( 4, 4, a, "  Matrix A:" );

  b = r8mat_inverse_4d ( a );

  r8mat_print ( 4, 4, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 4, 4, c, "  C = A * B:" );

  delete [] b;

  return;

# undef N
}
//****************************************************************************80

void test067 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST067 tests R8MAT_JAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double eps = 0.00001;
  double *fprime;
  int m = 3;
  double x[N] = { 1.0, 2.0, 3.0, 4.0 };

  cout << "\n";
  cout << "TEST067\n";
  cout << "  R8MAT_JAC estimates the M by N jacobian matrix\n";
  cout << "  of a nonlinear function.\n";

  fprime = r8mat_jac ( m, N, eps, test067_f, x );

  r8mat_print ( m, N, fprime, "  Estimated jacobian:" );

  delete [] fprime;

  fprime = test067_jac ( m, N, x );

  r8mat_print (  m, N, fprime, "  Exact jacobian:" );

  delete [] fprime;

  return;
# undef N
}
//****************************************************************************80

double *test067_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST067_F is a sample nonlinear function for treatment by R8MAT_JAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of functions.
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double TEST067_F[M], the function values.
//
{
  double *f;

  f = new double[m];

  f[0] = sin ( x[0] * x[1] );
  f[1] = sqrt ( 1.0 + x[0] * x[0] ) + x[2];
  f[2] = x[0] + 2.0 * x[1] + 3.0 * x[2] + 4.0 * x[3];

  return f;
}
//****************************************************************************80

double *test067_jac ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST067_JAC is the exact jacobian of TEST067_F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of functions.
//
//    Input, int N, the number of parameters.
//
//    Input, double X[N], the parameter values.
//
//    Output, double FPRIME[M*N], the jacobian values.
//
{
  double *fprime;

  fprime = new double[m*n];

  fprime[0+0*3] = cos ( x[0] * x[1] ) * x[1];
  fprime[0+1*3] = cos ( x[0] * x[1] ) * x[0];
  fprime[0+2*3] = 0.0;
  fprime[0+3*3] = 0.0;

  fprime[1+0*3] = x[0] / sqrt ( 1.0 + x[0] * x[0] );
  fprime[1+1*3] = 0.0;
  fprime[1+2*3] = 1.0;
  fprime[1+3*3] = 0.0;

  fprime[2+0*3] = 1.0;
  fprime[2+1*3] = 2.0;
  fprime[2+2*3] = 3.0;
  fprime[2+3*3] = 4.0;

  return fprime;
}
//****************************************************************************80

void test068 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST068 tests R8MAT_L_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
    1.0, 2.0, 4.0,  7.0,
    0.0, 3.0, 5.0,  8.0,
    0.0, 0.0, 6.0,  9.0,
    0.0, 0.0, 0.0, 10.0 };
  double *b;
  double *c;

  cout << "\n";
  cout << "TEST068\n";
  cout << "  R8MAT_L_INVERSE inverts a lower triangular matrix.\n";

  r8mat_print ( N, N, a, "  Matrix A to be inverted:" );

  b = r8mat_l_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test069 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST069 tests R8MAT_L_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double a1[28] = {
    11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          22.0, 32.0, 42.0, 52.0, 62.0, 72.0,
                33.0, 43.0, 53.0, 63.0, 73.0,
                      44.0, 54.0, 64.0, 74.0,
                            55.0, 65.0, 75.0,
                                  66.0, 76.0,
                                        77.0 };
  double a2[18] = {
    11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          22.0, 32.0, 42.0, 52.0, 62.0, 72.0,
                33.0, 43.0, 53.0, 63.0, 73.0 };
  double a3[10] = {
    11.0, 21.0, 31.0, 41.0,
          22.0, 32.0, 42.0,
                33.0, 43.0,
                      44.0 };
  int m1 = 7;
  int m2 = 7;
  int m3 = 4;
  int n1 = 7;
  int n2 = 3;
  int n3 = 7;

  cout << "\n";
  cout << "TEST069\n";
  cout << "  R8MAT_L_PRINT prints a lower triangular matrix\n";
  cout << "  stored compactly.  Only the (possibly) nonzero\n";
  cout << "  elements are printed.\n";

  r8mat_l_print ( m1, n1, a1, "  A 7 by 7 matrix." );

  r8mat_l_print ( m2, n2, a2, "  A 7 by 3 matrix." );

  r8mat_l_print ( m3, n3, a3, "  A 4 by 7 matrix." );

  return;
}
//****************************************************************************80

void test070 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST070 tests R8MAT_L1_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
     1.0, 2.0, 0.0, 5.0, 0.0, 75.0,
     0.0, 1.0, 0.0, 0.0, 0.0,  0.0,
     0.0, 0.0, 1.0, 3.0, 0.0,  0.0,
     0.0, 0.0, 0.0, 1.0, 0.0,  6.0,
     0.0, 0.0, 0.0, 0.0, 1.0,  4.0,
     0.0, 0.0, 0.0, 0.0, 0.0,  1.0 };
  double *b;
  double *c;

  cout << "\n";
  cout << "TEST070\n";
  cout << "  R8MAT_L1_INVERSE inverts a unit lower triangular matrix.\n";

  r8mat_print ( N, N, a, "  Matrix A to be inverted:" );

  b = r8mat_l1_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test071 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST071 tests R8MAT_LU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5

  double *a;
  double l[M*M];
  double *lu;
  double p[M*M];
  double *plu;
  double u[M*N];
  double x[N] = { 1.0, 10.0, 4.0, 2.0, 3.0 };

  cout << "\n";
  cout << "TEST071\n";
  cout << "  R8MAT_LU computes the LU factors of a matrix.\n";

  a = r8mat_vand2 ( N, x );

  r8mat_print ( M, N, a, "  Matrix to be factored:" );

  r8mat_lu ( M, N, a, l, p, u );

  r8mat_print ( M, M, p, "  P factor:" );

  r8mat_print ( M, M, l, "  L factor:" );

  r8mat_print ( M, N, u, "  U factor:" );

  lu = r8mat_mm_new ( M, M, N, l, u );

  plu = r8mat_mm_new ( M, M, N, p, lu );

  r8mat_print ( M, N, plu, "  P*L*U:" );

  delete [] a;
  delete [] lu;
  delete [] plu;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test072 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST072 tests R8MAT_MAX and R8MAT_MIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  cout << "\n";
  cout << "TEST072\n";
  cout << "  For a real matrix,\n";
  cout << "  R8MAT_MAX computes the maximum value;\n";
  cout << "  R8MAT_MIN computes the minimum value;\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  temp1 = r8mat_min ( m, n, a );
  temp2 = r8mat_max ( m, n, a );

  cout << "\n";
  cout << "  Minimum value = " << temp1 << "\n";
  cout << "  Maximum value = " << temp2 << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void test073 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST073 tests R8MAT_MAXCOL_MINROW, and its variations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  cout << "\n";
  cout << "TEST073\n";
  cout << "  R8MAT_MAXCOL_MINROW computes the maximum over\n";
  cout << "  columns of the mininum over rows;\n";
  cout << "  R8MAT_MAXROW_MINCOL computes the maximum over\n";
  cout << "  rows of the mininum over columns;\n";
  cout << "  R8MAT_MINCOL_MAXROW computes the minimum over\n";
  cout << "  columns of the maxinum over rows;\n";
  cout << "  R8MAT_MINROW_MAXCOL computes the minimum over\n";
  cout << "  rows of the maxinum over columns;\n";
  cout << "\n";

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, seed );

  r8mat_print ( m, n, a, "  Random array:" );

  cout << "  MAXCOL_MINROW = " << r8mat_maxcol_minrow ( m, n, a ) << "\n";
  cout << "  MINROW_MAXCOL = " << r8mat_minrow_maxcol ( m, n, a ) << "\n";
  cout << "  MAXROW_MINCOL = " << r8mat_maxrow_mincol ( m, n, a ) << "\n";
  cout << "  MINCOL_MAXROW = " << r8mat_mincol_maxrow ( m, n, a ) << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void test0731 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0731 tests R8MAT_MM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
# define N3 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2*N3] = {
     1.0,  2.0,  3.0,
     4.0,  5.0,  6.0,
     0.0,  0.0,  1.0,
    -1.0,  2.0, -1.0 };
  double *c;

  cout << "\n";
  cout << "TEST0731\n";
  cout << "  R8MAT_MM multiplies two (rectangular) matrices\n";
  cout << "  and returns the result as the function value.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8mat_print ( N2, N3, b, "  Matrix B:" );

  c = r8mat_mm_new ( N1, N2, N3, a, b );

  r8mat_print ( N1, N3, c, "  Product C = A * B:" );

  delete [] c;

  return;
# undef N1
# undef N2
# undef N3
}
//****************************************************************************80

void test0732 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0732 tests R8MAT_MXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
# define N3 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2*N3] = {
     1.0,  2.0,  3.0,
     4.0,  5.0,  6.0,
     0.0,  0.0,  1.0,
    -1.0,  2.0, -1.0 };
  double c[N1*N3];

  cout << "\n";
  cout << "TEST0732\n";
  cout << "  R8MAT_MXM multiplies two (rectangular) matrices\n";
  cout << "  and returns the result as an argument.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8mat_print ( N2, N3, b, "  Matrix B:" );

  r8mat_mxm ( N1, N2, N3, a, b, c );

  r8mat_print ( N1, N3, c, "  Product C = A * B:" );

  return;
# undef N1
# undef N2
# undef N3
}
//****************************************************************************80

void test0733 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0733 tests R8MAT_MV_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2] = {
     1.0,  2.0,  3.0 };
  double *c;

  cout << "\n";
  cout << "TEST0733\n";
  cout << "  R8MAT_MV_NEW multiplies a (rectangular) matrix times a vector,\n";
  cout << "  and returns the result as the function value.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N2, b, "  Vector B:" );

  c = r8mat_mv_new ( N1, N2, a, b );

  r8vec_print ( N1, c, "  Product C = A * B:" );

  delete [] c;

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void test0734 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0734 tests R8MAT_MV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2] = {
     1.0,  2.0,  3.0 };
  double c[N1];

  cout << "\n";
  cout << "TEST0734\n";
  cout << "  R8MAT_MV multiplies a (rectangular) matrix times a vector,\n";
  cout << "  and returns the result as an argument.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N2, b, "  Vector B:" );

  r8mat_mv ( N1, N2, a, b, c );

  r8vec_print ( N1, c, "  Product C = A * B:" );

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void test0735 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0735 tests R8MAT_MTV_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N1] = {
     1.0,  2.0 };
  double *c;

  cout << "\n";
  cout << "TEST0735\n";
  cout << "  R8MAT_MTV_NEW multiplies a transposed matrix times a vector,\n";
  cout << "  and returns the result as the function value.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N1, b, "  Vector B:" );

  c = r8mat_mtv_new ( N1, N2, a, b );

  r8vec_print ( N2, c, "  Product C = A' * B:" );

  delete [] c;

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void test0736 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0736 tests R8MAT_MTV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N1 2
# define N2 3
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N1] = {
     1.0,  2.0 };
  double c[N2];

  cout << "\n";
  cout << "TEST0736\n";
  cout << "  R8MAT_MTV multiplies a transposed matrix times a vector,\n";
  cout << "  and returns the result as an argument.\n";

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N1, b, "  Vector B:" );

  r8mat_mtv ( N1, N2, a, b, c );

  r8vec_print ( N2, c, "  Product C = A' * B:" );

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void test07365 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07365 tests R8MAT_NEW and R8MAT_DELETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double **a;
  double **b;
  int i;
  int j;
  int k;
  int m;
  int n;

  cout << "\n";
  cout << "TEST07365:\n";
  cout << "  R8MAT_NEW dynamically creates a 2D array.\n";
  cout << "  R8MAT_DELETE deletes it.\n";
  cout << "  Array entries can be addressed using the\n";
  cout << "  notation \"a[i][j]\".\n";
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by R8MAT_DELETE and new memory
//  requested by R8MAT_NEW.
//
  m = 4;
  n = 5;
//
//  Allocate memory.
//
  cout << "\n";
  cout << "  Allocating memory for array A of size " << m << " by " << n << ".\n";

  a = r8mat_new ( m, n );

  cout << "\n";
  cout << "  Assigning values to A.\n";
//
//  Store values in A.
//
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = ( double ) ( 10 * i + j );
    }
  }
//
//  Print A.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix A:\n";
  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << a[i][j];
    }
    cout << "\n";
  }
//
//  Create a new matrix B to store A' * A.
//
  b = r8mat_new ( n, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        b[i][j] = b[i][j] + a[k][i] * a[k][j];
      }
    }
  }
//
//  Print the matrix.
//
  cout << "\n";
  cout << "  Dynamically allocated matrix B = A' * A:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << b[i][j];
    }
    cout << "\n";
  }
//
//  Free memory.
//
  r8mat_delete ( a, m, n );
  r8mat_delete ( b, n, n );

  return;
}
//****************************************************************************80

void test0737 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0737 tests R8MAT_NULLSPACE_SIZE and R8MAT_NULLSPACE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  double *ax;
  int m = M;
  int n = N;
  double *nullspace;
  int nullspace_size;

  cout << "\n";
  cout << "TEST0737\n";
  cout << "  R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.\n";
  cout << "  R8MAT_NULLSPACE computes the nullspace of a matrix.\n";

  r8mat_print ( m, n, a, "  Input A:" );

  nullspace_size = r8mat_nullspace_size ( m, n, a );

  cout << "\n";
  cout << "  Nullspace size is " << nullspace_size << "\n";

  nullspace = r8mat_nullspace ( m, n, a, nullspace_size );

  r8mat_print ( n, nullspace_size, nullspace, "  Nullspace vectors:" );

  ax = r8mat_mxm_new ( m, n, nullspace_size, a, nullspace );

  r8mat_print ( m, nullspace_size, ax, "  Product A * Nullspace vectors:" );

  delete [] ax;
  delete [] nullspace;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test074 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST074 tests R8MAT_ORTH_UNIFORM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *at;
  double *ata;
  int n = 5;
  int seed;

  cout << "\n";
  cout << "TEST074\n";
  cout << "  R8MAT_ORTH_UNIFORM_NEW computes a random orthogonal matrix.\n";

  seed = 123456789;

  a = r8mat_orth_uniform_new ( n, seed );

  r8mat_print ( n, n, a, "  Random orthogonal matrix A" );

  at = r8mat_transpose_new ( n, n, a );

  ata = r8mat_mm_new ( n, n, n, at, a );

  r8mat_print ( n, n, ata, "  AT*A" );

  delete [] a;
  delete [] at;
  delete [] ata;

  return;
}
//****************************************************************************80

void test075 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST075 tests R8MAT_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define N 100

  double a[M*N];
  int i;
  int im1;
  int ip1;

  r8mat_zero ( M, N, a );

  for ( i = 0; i < M; i++ )
  {
    a[i+i*M] = -2.0;

    if ( i+1 < N )
    {
      ip1 = i+1;
    }
    else
    {
      ip1 = 0;
    }

    a[i+ip1*M] = 1.0;

    if ( 0 <= i-1 )
    {
      im1 = i-1;
    }
    else
    {
      im1 = N-1;
    }
    a[i+im1*M] = 1.0;
  }

  cout << "\n";
  cout << "TEST075\n";
  cout << "  R8MAT_PLOT prints a symbolic picture of a matrix.\n";
  cout << "  Typically,\n";
  cout << "\n";
  cout << "    - for negative, \n";
  cout << "    0 for zero, and\n";
  cout << "    + for positive entries\n";
  cout << "\n";
  cout << "  or\n";
  cout << "\n";
  cout << "    X for nonzero and\n";
  cout << "    0 for zero.\n";
  cout << "\n";

  r8mat_plot ( M, N, a, "  A plot of the matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test076 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST076 tests R8MAT_POWER_METHOD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N*N];
  double *av;
  int i;
  int j;
  double r;
  double v[N];

  cout << "\n";
  cout << "TEST076\n";
  cout << "  R8MAT_POWER_METHOD applies the power method\n";
  cout << "  to a matrix.\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( j == i - 1 || j == i + 1 )
      {
        a[i+j*N] = -1.0;
      }
      else if ( j == i )
      {
        a[i+j*N] = 2.0;
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }
  r8vec_zero ( N, v );

  r8mat_power_method ( N, a, &r, v );

  cout << "\n";
  cout << "  Estimated eigenvalue = " << r << "\n";

  r8vec_print ( N, v, "  Estimated eigenvector V:" );

  av = r8mat_mv_new ( N, N, a, v );

  r8vec_print ( N, av, "  Value of A*V:" );

  delete [] av;

  return;
# undef N
}
//****************************************************************************80

void test0764 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST076 tests R8MAT_REF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;

  cout << "\n";
  cout << "TEST0764\n";
  cout << "  R8MAT_REF computes the row echelon form of a matrix.\n";

  r8mat_print ( m, n, a, "  Input A:" );

  r8mat_ref ( m, n, a );

  r8mat_print ( m, n, a, "  REF form:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test0766 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0766 tests R8MAT_RREF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;

  cout << "\n";
  cout << "TEST0766\n";
  cout << "  R8MAT_RREF computes the reduced row echelon form of a matrix.\n";

  r8mat_print ( m, n, a, "  Input A:" );

  r8mat_rref ( m, n, a );

  r8mat_print ( m, n, a, "  REF form:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test077 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST077 tests R8MAT_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define RHS_NUM 2

  double a[N*(N+RHS_NUM)] = {
     1.0,  4.0,  7.0,
     2.0,  5.0,  8.0,
     3.0,  6.0,  0.0,
    14.0, 32.0, 23.0,
     7.0, 16.0,  7.0 };
  int i;
  int info;
  int j;

  cout << "\n";
  cout << "TEST077\n";
  cout << "  R8MAT_SOLVE solves linear systems.\n";
//
//  Print out the matrix to be inverted.
//
  r8mat_print ( N, N+RHS_NUM, a, "  The linear system:" );
//
//  Solve the systems.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  The input matrix was singular.\n";
    cout << "  The solutions could not be computed.\n";
    return;
  }

  cout << "\n";
  cout << "  The computed solutions:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    for ( j = N; j < N+RHS_NUM; j++ )
    {
      cout << setw(10) << a[i+j*N] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef RHS_NUM
}
//****************************************************************************80

void test0775 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0775 tests R8MAT_SOLVE_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double det;
  int i;
  int n = 2;
  int seed;
  int test;
  int test_num = 5;
  double *x;
  double *x2;

  cout << "\n";
  cout << "TEST0775\n";
  cout << "  R8MAT_SOLVE_2D solves 2D linear systems.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8mat_uniform_01_new ( n, n, seed );
    x = r8vec_uniform_01_new ( n, seed );
    b = r8mat_mv_new ( n, n, a, x );

    x2 = r8mat_solve_2d ( a, b, &det );

    cout << "\n";
    cout << "  Solution / Computed:\n";
    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(14) << x[i]
           << "  " << setw(14) << x2[i] << "\n";
    }

    delete [] a;
    delete [] b;
    delete [] x;
    delete [] x2;
  }

  return;
}
//****************************************************************************80

void test0776 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0776 tests R8MAT_SOLVE_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double det;
  int i;
  int n = 3;
  int seed;
  int test;
  int test_num = 5;
  double *x;
  double *x2;

  cout << "\n";
  cout << "TEST0776\n";
  cout << "  R8MAT_SOLVE_3D solves 3D linear systems.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8mat_uniform_01_new ( n, n, seed );
    x = r8vec_uniform_01_new ( n, seed );
    b = r8mat_mv_new ( n, n, a, x );

    x2 = r8mat_solve_3d ( a, b, &det );

    cout << "\n";
    cout << "  Solution / Computed:\n";
    cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(14) << x[i]
           << "  " << setw(14) << x2[i] << "\n";
    }

    delete [] a;
    delete [] b;
    delete [] x;
    delete [] x2;
  }

  return;
}
//****************************************************************************80

void test078 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST078 tests R8MAT_SOLVE2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double *a;
  double a1[2*2] = {
    1.0, 3.0,
    2.0, 4.0 };
  double a2[3*3] = {
    2.0, 1.0, 1.0,
    1.0, 1.0, 0.0,
    1.0, 0.0, 1.0 };
  double a3[4*4] = {
    1.0, 2.0, 1.0, 3.0,
    0.0, 1.0, 2.0, 1.0,
    0.0, 0.0, 3.0, 2.0,
    1.0, 3.0, 0.0, 1.0 };
  double a4[3*3] = {
    2.0, 1.0, 3.0,
    4.0, 2.0, 6.0,
    1.0, 4.0, 5.0 };
  double *b;
  double b1[2] = { 5.0, 11.0 };
  double b2[3] = { 4.0, 2.0, 2.0 };
  double b3[4] = { 5.0, 11.0, 16.0, 15.0 };
  double b4[3] = { 13.0, 17.0, 20.0 };
  int ierror;
  int n;
  int n_test[TEST_NUM] = { 2, 3, 4, 3 };
  int test;
  double *x;

  cout << "\n";
  cout << "TEST078\n";
  cout << "  R8MAT_SOLVE2 is a linear solver.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];

    if ( test == 0 )
    {
      a = a1;
      b = b1;
    }
    else if ( test == 1 )
    {
      a = a2;
      b = b2;
    }
    else if ( test == 2 )
    {
      a = a3;
      b = b3;
    }
    else if ( test == 3 )
    {
      a = a4;
      b = b4;
    }

    r8vec_print ( n, b, "  Right hand side:" );

    x = r8mat_solve2 ( n, a, b, &ierror );

    cout << "\n";
    if ( ierror == 0 )
    {
      cout << "  The system is nonsingular.\n";
    }
    else if ( ierror == 1 )
    {
      cout << "  The system is singular, but consistent.\n";
    }
    else if ( ierror == 2 )
    {
      cout << "  The system is singular and inconsistent.\n";
    }

    r8vec_print ( n, x, "  Computed solution:" );

    delete [] x;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test079 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST079 tests R8MAT_SYMM_JACOBI;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int n = 5;
  double *q;
  int seed;
  double *x;

  cout << "\n";
  cout << "TEST079\n";
  cout << "  For a symmetric R8MAT:\n";
  cout << "  R8MAT_SYMM_JACOBI diagonalizes;\n";
//
//  Choose eigenvalues.
//
  x = r8vec_indicator_new ( n );
//
//  Choose eigenvectors.
//
  seed = 123456789;

  q = r8mat_orth_uniform_new ( n, seed );
//
//  Now get A = Q*X*Q.
//
  a = r8mat_symm_eigen ( n, x, q );

  r8mat_print ( n, n, a, "  Matrix to diagonalize:" );

  r8mat_symm_jacobi ( n, a );

  cout << "\n";
  cout << "  Computed Eigenvalues:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(12) << a[i+i*n] << "\n";
  }

  delete [] a;
  delete [] q;
  delete [] x;

  return;
}
//****************************************************************************80

void test080 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST080 tests R8MAT_TO_R8PLU and R8PLU_TO_R8MAT;
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
{
# define N 5

  double *a;
  double a2[N*N];
  double b = 0.0;
  double c = 1.0;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST080\n";
  cout << "  R8MAT_TO_R8PLU determines the compressed PLU factors\n";
  cout << "  of a real general matrix.\n";
  cout << "  R8PLU_TO_R8MAT determines the original matrix from\n";
  cout << "  the compressed PLU factors.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_ab_new ( N, N, b, c, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
  }
//
//  Display the gory details.
//
  i4vec_print ( N, pivot, "  The pivot vector P:" );

  r8mat_print ( N, N, lu, "  The compressed LU factors:" );
//
//  Recover the matrix from the PLU factors.
//
  r8plu_to_r8mat ( N, pivot, lu, a2 );

  r8mat_print ( N, N, a2, "  The recovered matrix A2:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test081 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST081 tests R8MAT_TRACE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N];
  int i;
  int j;
  double trace;

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i <= j )
      {
        a[i+j*N] = ( double ) ( N - j );
      }
      else if ( j == i - 1 )
      {
        a[i+j*N] = ( double ) ( N - j - 1 );
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  cout << "\n";
  cout << "TEST081\n";
  cout << "  R8MAT_TRACE computes the trace of a matrix\n";

  r8mat_print ( N, N, a, "  Matrix:" );

  trace = r8mat_trace ( N, a );

  cout << "\n";
  cout << "  Trace is " << trace << "\n";

  return;
# undef N
}
//****************************************************************************80

void test082 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST082 tests R8MAT_TRANSPOSE_PRINT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 7
# define N 12

  double a[M*N];
  int i;
  int j;

  cout << "\n";
  cout << "TEST082\n";
  cout << "  R8MAT_TRANSPOSE_PRINT prints an R8MAT,\n";
  cout << "  transposed.\n";
  cout << "\n";
  cout << "  Matrix row order M =    " << M << "\n";
  cout << "  Matrix column order N = " << N << "\n";
//
//  Set the matrix.
//
  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      a[i-1+(j-1)*M] = ( double ) ( i * 100 + j );
    }
  }

  r8mat_transpose_print ( M, N, a, "  The transposed matrix A:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test083 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST083 tests R8MAT_U_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
    1.0, 0.0, 0.0,  0.0,
    2.0, 3.0, 0.0,  0.0,
    4.0, 5.0, 6.0,  0.0,
    7.0, 8.0, 9.0, 10.0 };
  double *b;
  double *c;
  int i;

  cout << "\n";
  cout << "TEST083\n";
  cout << "  R8MAT_U_INVERSE inverts an upper triangular matrix.\n";

  r8mat_print ( N, N, a, "  Input matrix A" );

  b = r8mat_u_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test084 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST084 tests R8MAT_U1_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
//
//  Each row of this definition is a COLUMN of the matrix.
//
  double a[N*N] = {
    1.0, 0.0, 0.0, 0.0, 0.0,  0.0,
    2.0, 1.0, 0.0, 0.0, 0.0,  0.0,
    0.0, 0.0, 1.0, 0.0, 0.0,  0.0,
    5.0, 0.0, 3.0, 1.0, 0.0,  0.0,
    0.0, 0.0, 0.0, 0.0, 1.0,  0.0,
   75.0, 0.0, 0.0, 6.0, 4.0,  1.0 };
  double *b;
  double *c;

  cout << "\n";
  cout << "TEST084\n";
  cout << "  R8MAT_U1_INVERSE inverts a unit upper triangular matrix.\n";

  r8mat_print ( N, N, a, "  Input matrix A" );

  b = r8mat_u1_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests R8MAT_UNIFORM.
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
{
# define M 5
# define N 4

  double *a;
  double b = 2.0E+00;
  double c = 10.0E+00;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST085\n";
  cout << "  R8MAT_UNIFORM sets a matrix to random values.\n";
  cout << "\n";

  a = r8mat_uniform_ab_new ( M, N, b, c, seed );
//
//  Print out the matrix to be inverted.
//
  r8mat_print ( M, N, a, "  The random matrix:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test086 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST086 tests R8PLU_DET;
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
{
# define N 5

  double *a;
  double b = 0.0;
  double c = 1.0;
  double det;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST086\n";
  cout << "  R8PLU_DET determines the determinant of a matrix from its\n";
  cout << "  compressed PLU factors.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_ab_new ( N, N, b, c, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
  }
//
//  Compute the determinant.
//
  det = r8plu_det ( N, pivot, lu );

  cout << "\n";
  cout << "  The determinant = " << det << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test087 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST087 tests R8PLU_INVERSE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double b[N*N];
  double *c;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST087\n";
  cout << "  R8PLU_INVERSE determines the inverse of a matrix from its\n";
  cout << "  compressed PLU factors.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( N, N, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
  }
//
//  Compute the inverse.
//
  r8plu_inverse ( N, pivot, lu, b );

  r8mat_print ( N, N, b, "  The inverse B:" );
//
//  Compute the product C = A * B.
//
  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] a;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test088 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST088 tests R8PLU_MUL;
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
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;
  double x[N];

  cout << "\n";
  cout << "TEST088\n";
  cout << "  R8PLU_MUL computes the product A*x\n";
  cout << "  using the compressed PLU factors of A.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( N, N, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Set the right hand side B1.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) (i+1);
  }

  b = r8mat_mv_new ( N, N, a, x );

  r8vec_print ( N, b, "  The right hand side B (computed from A):" );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );
//
//  Solve the system.
//
  r8plu_mul ( N, pivot, lu, x, b );

  r8vec_print ( N, b, "  The right hand side B (computed from PLU):" );

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void test089 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST089 tests R8PLU_SOL;
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
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;
  double x[N];

  cout << "\n";
  cout << "TEST089\n";
  cout << "  R8PLU_SOL solves the linear system A*x=b\n";
  cout << "  using the compressed PLU factors of A.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8mat_uniform_01_new ( N, N, seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
//
//  Set the desired solution.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) (i+1);
  }
//
//  Set the right hand side.
//
  b = r8mat_mv_new ( N, N, a, x );

  r8vec_print ( N, b, "  The right hand side B (computed from A):" );
//
//  Destroy the desired solution (no cheating now!)
//
  r8vec_zero ( N, x );
//
//  Factor the matrix.
//
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  R8MAT_TO_R8PLU declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  r8plu_sol ( N, pivot, lu, b, x );

  r8vec_print ( N, x, "  The computed solution X:" );

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void test090 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST090 tests R8POLY_DERIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double *c;
  double *cp;
  int d;
  double *x;

  cout << "\n";
  cout << "TEST090\n";
  cout << "  R8POLY_DERIV computes the coefficients of\n";
  cout << "  the derivative of a polynomial.\n";

  x = r8vec_indicator_new ( N );

  c = roots_to_r8poly ( N, x );

  r8poly_print ( N, c, "  The initial polynomial" );

  for ( d = 0; d <= N; d++ )
  {
    cp = r8poly_deriv ( N, c, d );
    cout << "\n";
    cout << "  The derivative of order " << d << "\n";
    cout << "\n";
    r8poly_print ( N-d, cp, " " );
    delete [] cp;
  }

  delete [] c;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test091 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST091 tests R8POLY_LAGRANGE_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 3

  int i;
  int ipol;
  double *pcof;
  double *xpol;

  cout << "\n";
  cout << "TEST091\n";
  cout << "  R8POLY_LAGRANGE_COEF returns the coefficients for\n";
  cout << "  a Lagrange basis polynomial.\n";

  xpol = r8vec_indicator_new ( NPOL );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );

  for ( ipol = 1; ipol <= NPOL; ipol++ )
  {
    pcof = r8poly_lagrange_coef ( NPOL, ipol, xpol );

    cout << "\n";
    cout << "  Lagrange basis polynomial " << setw(4) << ipol << ":\n";
    cout << "\n";

    for ( i = 0; i < NPOL; i++ )
    {
      cout << setw(10) << pcof[i] << "  "
           << setw(4)  << i << "\n";
    }
    delete [] pcof;

  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void test092 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST092 tests R8POLY_LAGRANGE_COEF and R8POLY_DERIV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  int d;
  int ipol;
  double *pcof;
  double *pprime;
  double *xpol;

  cout << "\n";
  cout << "TEST092\n";
  cout << "  R8POLY_LAGRANGE_COEF returns the coefficients\n";
  cout << "  for a Lagrange basis polynomial.\n";
  cout << "  R8POLY_DERIV computes derivatives of a polynomial.\n";

  xpol = r8vec_indicator_new ( NPOL );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );

  for ( ipol = 1; ipol <= NPOL; ipol++ )
  {
    pcof = r8poly_lagrange_coef ( NPOL, ipol, xpol );

    r8poly_print ( NPOL-1, pcof, "  The Lagrange basis polynomial:" );

    for ( d = 1; d <= NPOL-1; d++ )
    {
      pprime = r8poly_deriv ( NPOL-1, pcof, d );
      cout << "\n";
      cout << "  The derivative of order " << d << "\n";
      r8poly_print ( NPOL-1-d, pprime, " " );
      delete [] pprime;
    }

    delete [] pcof;
  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void test093 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST093 tests R8POLY_LAGRANGE_0, R8POLY_LAGRANGE_1 and R8POLY_LAGRANGE_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  double dw2dx2;
  double dwdx;
  int ival;
  int nx;
  double wval;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  cout << "\n";
  cout << "TEST093\n";
  cout << "  R8POLY_LAGRANGE_0 evaluates the Lagrange\n";
  cout << "  factor W(X) at a point.\n";
  cout << "  R8POLY_LAGRANGE_1 evaluates the Lagrange\n";
  cout << "  factor W'(X) at a point.\n";
  cout << "  R8POLY_LAGRANGE_2 evaluates the Lagrange\n";
  cout << "  factor W''(X) at a point.\n";
  cout << "\n";
  cout << "  The number of data points is " << NPOL << "\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );

  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate W(X), W'(X), W''.
//
  cout << "\n";
  cout << "      X          W(X)          W'(X)        W''(X)\n";
  cout << "\n";

  nx = 4 * NPOL - 1;

  for ( ival = 1; ival <= nx; ival++ )
  {
    xval = r8vec_even_select ( nx, xlo, xhi, ival );

    wval = r8poly_lagrange_0 ( NPOL, xpol, xval );
    dwdx = r8poly_lagrange_1 ( NPOL, xpol, xval );
    dw2dx2 = r8poly_lagrange_2 ( NPOL, xpol, xval );

    cout << setw(12) << xval   << "  "
         << setw(12) << wval   << "  "
         << setw(12) << dwdx   << "  "
         << setw(12) << dw2dx2 << "\n";
  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void test094 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST094 tests R8POLY_LAGRANGE_FACTOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  double dwdx;
  int i;
  double wval;
  double xhi;
  double xlo;
  double xpol[NPOL];
  double xval;

  cout << "\n";
  cout << "TEST094\n";
  cout << "  R8POLY_LAGRANGE_FACTOR evaluates the Lagrange\n";
  cout << "  factor W(X) at a point.\n";
  cout << "\n";
  cout << "  For this test, we use " << NPOL << " functions.\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0;
  xhi = ( double ) ( NPOL - 1 );
  for ( i = 0; i < NPOL; i++ )
  {
    xpol[i] = ( ( double ) ( NPOL - i ) * xlo + ( double ) i * xhi )
      / ( double ) ( NPOL );
  }

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate W(X) and W'(X).
//
  cout << "\n";
  cout << "      X          W(X)          W'(X)\n";
  cout << "\n";

  for ( i = 0; i < 2 * NPOL - 2; i++ )
  {
    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, i );

    r8poly_lagrange_factor ( NPOL, xpol, xval, &wval, &dwdx );

    cout << setw ( 10 ) << xval << " "
         << setw ( 10 ) << wval << " "
         << setw ( 10 ) << dwdx << "\n";
  }

  return;
# undef NPOL
}
//****************************************************************************80

void test095 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST095 tests R8POLY_LAGRANGE_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NPOL 5

  int i;
  int ipol;
  int ival;
  double p1;
  double p2;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  cout << "\n";
  cout << "TEST095\n";
  cout << "  R8POLY_LAGRANGE_VAL evaluates a Lagrange\n";
  cout << "  interpolating polynomial at a point.\n";
  cout << "\n";
  cout << "  For this test, we use " << NPOL << " functions.\n";
//
//  Set the abscissas of the polynomials.
//
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );
  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
//
//  Evaluate the polynomials.
//
  cout << "\n";
  cout << "  Here are the values of the functions at\n";
  cout << "  several points:\n";
  cout << "\n";
  cout << "      X          L1          L2          L3      L4          L5\n";
  cout << "\n";

  for ( ival = 0; ival < 2 * NPOL - 1; ival++ )
  {

    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, ival );
    cout << setw(10) << xval << "  ";

    for ( ipol = 0; ipol < NPOL; ipol++ )
    {
      r8poly_lagrange_val ( NPOL, ipol, xpol, xval, &p1, &p2 );
      cout << setw(10) << p1 << "  ";
    }
    cout << "\n";
  }
  cout << "\n";
  cout << "  And the derivatives:\n";
  cout << "\n";
  cout << "      X          L'1         L'2         L'3     L'4         L'5\n";
  cout << "\n";

  for ( ival = 0; ival < 2 * NPOL - 1; ival++ )
  {
    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, ival );
    cout << setw(10) << xval << " ";

    for ( ipol = 0; ipol < NPOL; ipol++ )
    {
      r8poly_lagrange_val ( NPOL, ipol, xpol, xval, &p1, &p2 );
      cout << setw(10) << p2 << " ";
    }
    cout << "\n";
  }

  delete [] xpol;

  return;
# undef NPOL
}
//****************************************************************************80

void test098 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST098 tests R8POLY_VALUE_HORNER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  double p[N+1] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  double pval;
  double x;

  cout << "\n";
  cout << "TEST098\n";
  cout << "  R8POLY_VALUE_HORNER evaluates a polynomial at a\n";
  cout << "  point, using Horner's method.\n";

  r8poly_print ( N, p, "  The polynomial:" );

  cout << "\n";
  cout << "        X            P(X)\n";
  cout << "\n";

  for ( i = 0; i <= 15; i++ )
  {
    x = ( double ) ( i ) / 3.0;
    pval = r8poly_value_horner ( N, p, x );
    cout << "  " << setw(14) << x
         << "  " << setw(14) << pval << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test099 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST099 tests R8POLY2_EX and R8POLY2_EX2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  int ierror;
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  cout << "\n";
  cout << "TEST099\n";
  cout << "  R8POLY2_EX finds the extreme value\n";
  cout << "  of a parabola determined by three points.\n";
  cout << "  R8POLY2_EX2 finds the extreme value\n";
  cout << "  of a parabola determined by three points.\n";

  a =  2.0;
  b = -4.0;
  c = 10.0;

  x1 = 1.0;
  y1 = a * x1 * x1 + b * x1 + c;
  x2 = 2.0;
  y2 = a * x2 * x2 + b * x2 + c;
  x3 = 3.0;
  y3 = a * x3 * x3 + b * x3 + c;

  cout << "\n";
  cout << "  Parabolic coefficients A = "
    << a << ", B = " << b << ", c = " << c << "\n";
  cout << "\n";
  cout << "  X, Y data:\n";
  cout << "\n";
  cout << "  " << x1 << "  " << y1;
  cout << "  " << x2 << "  " << y2;
  cout << "  " << x3 << "  " << y3;

  a = 0.0;
  b = 0.0;
  c = 0.0;

  ierror = r8poly2_ex ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( ierror == 0 )
  {
    cout << "\n";
    cout << "  R8POLY2_EX returns XMIN = "
      << xmin << ", YMIN = " << ymin << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8POLY2_EX returns error code " << ierror << ".\n";
  }

  ierror = r8poly2_ex2 ( x1, y1, x2, y2, x3, y3, &xmin, &ymin, &a, &b, &c );

  if ( ierror == 0 )
  {
    cout << "\n";
    cout << "  R8POLY2_EX2 returns XMIN = "
      << xmin << ", YMIN = " << ymin << "\n";
    cout << "  and A = " << a << ", B = " << b << ", c = " << c << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8POLY2_EX2 returns error code " << ierror << ".\n";
  }

  return;
}
//****************************************************************************80

void test100 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST100 tests R8POLY2_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double x;
  double x1;
  double x2;
  double x3;
  double y;
  double y1;
  double y2;
  double y3;
  double yp;
  double ypp;

  cout << "\n";
  cout << "TEST100\n";
  cout << "  R8POLY2_VAL evaluates a parabola given\n";
  cout << "  3 data points.\n";
  cout << "\n";
  cout << "  Our parabola will be 2*x^2 + 3 * x + 1.\n";
  cout << "\n";
  cout << "  Case 1: 3 distinct data points:\n";
  cout << "\n";

  x1 = -1.0;
  x2 = 1.0;
  x3 = 3.0;

  test100_f ( x1, &y1, &yp, &ypp );
  test100_f ( x2, &y2, &yp, &ypp );
  test100_f ( x3, &y3, &yp, &ypp );

  cout << "  " << x1 << " " << y1 << "\n";
  cout << "  " << x2 << " " << y2 << "\n";
  cout << "  " << x3 << " " << y3 << "\n";

  cout << "\n";
  cout << "  Sampled data:\n";
  cout << "\n";
  cout << "  X, Y, Y', Y''\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    cout << "  " << x << "  " << y << "  " << yp << "  " << ypp << "\n";
  }

  cout << "\n";
  cout << "  Case 2: X1=X2, X3 distinct:\n";
  cout << "\n";

  x1 = -1.0;
  x2 = -1.0;
  x3 = 3.0;

  test100_f ( x1, &y1, &y2, &ypp );
  test100_f ( x3, &y3, &yp, &ypp );

  cout << "  " << x1 << "  " << y1 << "\n";
  cout << "  " << x2 << "  " << y2 << "\n";
  cout << "  " << x3 << "  " << y3 << "\n";

  cout << "\n";
  cout << "  Sampled data:\n";
  cout << "\n";
  cout << "   X, Y, Y', Y''\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    cout << "  " << x << "  " << y << "  " << yp << "  " << ypp << "\n";
  }

  cout << "\n";
  cout << "  Case 3: X1=X2=X3:\n";
  cout << "\n";

  x1 = -1.0;
  x2 = -1.0;
  x3 = -1.0;

  test100_f ( x1, &y1, &y2, &y3 );

  cout << "  " << x1 << "  " << y1 << "\n";
  cout << "  " << x2 << "  " << y2 << "\n";
  cout << "  " << x3 << "  " << y3 << "\n";

  cout << "\n";
  cout << "  Sampled data:\n";
  cout << "\n";
  cout << "  X, Y, Y', Y''\n";
  cout << "\n";
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    cout << "  " << x << "  " << y << "  " << yp << "  " << ypp << "\n";
  }

  return;
}
//****************************************************************************80

void test100_f ( double x, double *y, double *yp, double *ypp )

//****************************************************************************80
//
//  Purpose:
//
//    TEST100_F evaluates a parabola for us.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
{
  *y = 2.0 * x * x + 3.0 * x + 1.0;
  *yp = 4.0 * x + 3.0;
  *ypp = 4.0;

  return;
}
//****************************************************************************80

void test101 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST101 tests R8POLY2_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 5

  int i;
  int left;
  double xdata[NDATA];
  double xval;
  double ydata[NDATA];
  double yval;
  double zdata[NDATA];
  double zval;

  cout << "\n";
  cout << "TEST101\n";
  cout << "  R8POLY2_VAL2 evaluates parabolas through\n";
  cout << "  3 points in a table\n";
  cout << "\n";
  cout << "  Our data tables will actually be parabolas:\n";
  cout << "    A: 2*x^2 + 3 * x + 1.\n";
  cout << "    B: 4*x^2 - 2 * x + 5.\n";
  cout << "\n";

  for ( i = 0; i < NDATA; i++ )
  {
    xval = 2.0 * ( double ) ( i + 1 );
    xdata[i] = xval;
    ydata[i] = 2.0 * xval * xval + 3.0 * xval + 1.0;
    zdata[i] = 4.0 * xval * xval - 2.0 * xval + 5.0;
    cout << setw(6)  << i << " "
         << setw(10) << xdata[i] << "  "
         << setw(10) << ydata[i] << "  "
         << setw(10) << zdata[i] << "\n";
  }

  cout << "\n";
  cout << "  Interpolated data:\n";
  cout << "\n";
  cout << "  LEFT, X, Y1, Y2\n";
  cout << "\n";

  for ( i = 0; i <= 4; i++ )
  {
    xval = ( double ) ( 2 * i + 1 );
    left = i;
    if ( NDATA - 3 < left )
    {
      left = NDATA - 3;
    }
    if ( left < 0 )
    {
      left = 0;
    }
    r8poly2_val2 ( NDATA, xdata, ydata, left, xval, &yval );
    r8poly2_val2 ( NDATA, xdata, zdata, left, xval, &zval );

    cout << setw(6)  << left << "  "
         << setw(10) << xval << "  "
         << setw(10) << yval << "  "
         << setw(10) << zval << "\n";
  }

  return;
# undef NDATA
}
//****************************************************************************80

void test105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST105 tests R8ROW_MAX and R8ROW_MIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  double *amax;
  double *amin;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST105\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_MAX computes maximums;\n";
  cout << "  R8ROW_MIN computes minimums;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  amax = r8row_max ( M, N, a );

  amin = r8row_min ( M, N, a );

  cout << "\n";
  cout << "  Row maximum, minimum:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(10) << amax[i]
         << "  " << setw(10) << amin[i] << "\n";
  }

  delete [] amax;
  delete [] amin;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test106 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST106 tests R8ROW_MEAN and R8ROW_SUM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *mean;
  double *rowsum;

  cout << "\n";
  cout << "TEST106\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_MEAN computes means;\n";
  cout << "  R8ROW_SUM computes sums;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  rowsum = r8row_sum ( M, N, a );

  mean = r8row_mean ( M, N, a );

  cout << "\n";
  cout << "  Row sum, mean:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(10) << rowsum[i]
         << "  " << setw(10) << mean[i] << "\n";
  }

  delete [] rowsum;
  delete [] mean;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test107 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST107 tests R8ROW_SWAP;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int row1;
  int row2;
  int j;
  int k;

  cout << "\n";
  cout << "TEST107\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_SWAP swaps two rows;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  row1 = 1;
  row2 = 3;

  cout << "\n";
  cout << "  Swap rows " << row1 << " and " << row2 << "\n";

  r8row_swap ( M, N, a, row1, row2 );

  r8mat_print ( M, N, a, "  The modified matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test108 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST108 tests R8ROW_TO_R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *x;

  cout << "\n";
  cout << "TEST108\n";
  cout << "  R8ROW_TO_R8VEC converts an array of rows into a vector.\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) ( 10 * i + j );
    }
  }

  r8mat_print ( M, N, a, "  The array of rows:" );

  x = r8row_to_r8vec ( M, N, a );

  r8vec_print ( M*N, x, "  The resulting vector of rows:" );

  delete [] x;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test109 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST109 tests R8ROW_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *variance;

  cout << "\n";
  cout << "TEST109\n";
  cout << "  For an R8ROW (a matrix regarded as rows):\n";
  cout << "  R8ROW_VARIANCE computes variances;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  variance = r8row_variance ( M, N, a );

  cout << "\n";
  cout << "  Row variances:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(10) << variance[i] << "\n";
  }

  delete [] variance;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test110 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST110 tests R8SLMAT_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double *a;
  double a1[21] = {
    21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          32.0, 42.0, 52.0, 62.0, 72.0,
                43.0, 53.0, 63.0, 73.0,
                      54.0, 64.0, 74.0,
                            65.0, 75.0,
                                  76.0 };
  double a2[15] = {
    21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          32.0, 42.0, 52.0, 62.0, 72.0,
                43.0, 53.0, 63.0, 73.0 };
  double a3[6] = {
    21.0, 31.0, 41.0,
          32.0, 42.0,
                43.0 };
  int m;
  int m_test[TEST_NUM] = { 7, 7, 4 };
  int n;
  int n_test[TEST_NUM] = { 7, 3, 7 };
  int size;
  int size_test[TEST_NUM] = { 21, 15, 6 };
  int test;

  cout << "\n";
  cout << "TEST110\n";
  cout << "  R8SLMAT_PRINT prints a strictly lower triangular matrix\n";
  cout << "  stored compactly.  Only the (possibly) nonzero \n";
  cout << "  elements are printed.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    m = m_test[test];
    n = n_test[test];
    size = size_test[test];
    a = new double[size];

    if ( test == 0 )
    {
      r8vec_copy ( size, a1, a );
    }
    else if ( test == 1 )
    {
      r8vec_copy ( size, a2, a );
    }
    else if ( test == 2 )
    {
      r8vec_copy ( size, a3, a );
    }

    r8slmat_print ( m, n, a, "  R8SLMAT:" );

    delete [] a;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test111 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST111 tests R8VEC_AMAX and R8VEC_AMIN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "TEST111\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_AMAX:      maximum magnitude entry;\n";
  cout << "  R8VEC_AMIN:      minimum magnitude entry.\n";

  b = - ( double ) N;
  c =  ( double ) N;

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  aval = r8vec_amax ( N, a );
  cout << "  Maximum absolute:         " << aval << "\n";

  aval = r8vec_amin ( N, a );
  cout << "  Minimum absolute:         " << aval << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test112 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST112 tests R8VEC_BRACKET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define TEST_NUM 6

  int i;
  int left;
  int right;
  int test;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "TEST112\n";
  cout << "  R8VEC_BRACKET finds a pair of entries in a\n";
  cout << "  sorted real array which bracket a value.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!)" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    xval = xtest[test];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    r8vec_bracket ( N, x, xval, &left, &right );

    cout << "  X[" << left  << "-1] = " << x[left-1]  << "\n";
    cout << "  X[" << right << "-1] = " << x[right-1] << "\n";
  }

  return;

# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test113 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST113 tests R8VEC_BRACKET2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define TEST_NUM 6

  int i;
  int left;
  int right;
  int start;
  int test;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "TEST113\n";
  cout << "  R8VEC_BRACKET2 finds a pair of entries in a\n";
  cout << "  sorted R8VEC which bracket a value.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!)" );

  for ( test = 0; test < TEST_NUM; test++ )
  {

    xval = xtest[test];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    if ( 0 < left )
    {
      start = left;
    }
    else
    {
      start = ( N + 1 ) / 2;
    }

    cout << "  Start = " << start << "\n";

    r8vec_bracket2 ( N, x, xval, start, &left, &right );

    cout << "  Left =  " << left  << "\n";
    cout << "  Right = " << right << "\n";

    if ( 1 <= left )
    {
      cout << "  X[" << left  << "-1] = " << x[left-1]  << "\n";
    }

    if ( 1 <= right )
    {
      cout << "  X[" << right << "-1] = " << x[right-1] << "\n";
    }
  }
  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test114 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST114 tests R8VEC_BRACKET3.
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
# define N 10
# define TEST_NUM 6

  int i;
  int itest;
  int left;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "TEST114\n";
  cout << "  R8VEC_BRACKET3 finds a pair of entries in a\n";
  cout << "  sorted real array which bracket a value.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!):" );

  left = ( N - 1 ) / 2;

  for ( itest = 0; itest < TEST_NUM; itest++ )
  {
    xval = xtest[itest];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    cout << "  Starting guess for interval is = " << left << "\n";

    r8vec_bracket3 ( N, x, xval, &left );

    cout << "  Nearest interval:\n";
    cout << "   X[" << left   << "]= " << x[left  ] << "\n";
    cout << "   X[" << left+1 << "]= " << x[left+1] << "\n";
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test1143 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1143 tests R8VEC_BRACKET5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int left;
  int n = 10;
  int right;
  int test;
  int test_num = 6;
  double *x;
  double xtest[6] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "TEST1143\n";
  cout << "  R8VEC_BRACKET5 finds a pair of entries in a\n";
  cout << "  sorted R8VEC which bracket a value.\n";

  x = r8vec_indicator_new ( n );
  x[5] = x[4];

  r8vec_print ( n, x, "  Sorted array:" );

  cout << "\n";
  cout << "        LEFT                   RIGHT\n";
  cout << "      X(LEFT)       XVAL     X(RIGHT)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    xval = xtest[test];

    left = r8vec_bracket5 ( n, x, xval );

    if ( left == -1 )
    {
      cout << "  " << setw(10) <<left << "\n";
      cout << "              " << setw(10) << xval << "  (Not bracketed!)\n";
    }
    else
    {
      right = left + 1;
      cout << "  " << setw(10) << left
           << "  " << "          "
           << "  " << setw(10) << right << "\n";
      cout << "  " << setw(10) << x[left]
           << "  " << setw(10) << xval
           << "  " << setw(10) << x[right] << "\n";
    }
  }

  delete [] x;

  return;
}
//****************************************************************************80

void test1145 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1145 tests R8VEC_CHEBYSPACE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double *r;
  double r1;
  double r2;

  cout << "\n";
  cout << "TEST1145\n";
  cout << "  R8VEC_CHEBYSPACE_NEW computes N Chebyshev points in [R1,R2].\n";

  r1 = -1.0;
  r2 = +1.0;
  n = 5;

  r = r8vec_chebyspace_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n << ", R1 = " << r1 << ", R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Chebyshev points:" );

  delete [] r;

  r1 =   0.0;
  r2 = +10.0;
  n = 7;

  r = r8vec_chebyspace_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n << ", R1 = " << r1 << ", R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Chebyshev points:" );

  delete [] r;

  return;
}
//****************************************************************************80

void test1147 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1147 tests R8VEC_CONVOLUTION
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2012
//
//  Author:
//
//    John Burkardt
//
{
# define M 4
# define N 3

  int m = M;
  int n = N;

  double x[M] = { 1.0, 2.0, 3.0, 4.0 };
  double y[N] = { -1.0, 5.0, 3.0 };
  double *z;
  double z_correct[M+N-1] = { -1.0, 3.0, 10.0, 17.0, 29.0, 12.0 };

  cout << "\n";
  cout << "TEST1147\n";
  cout << "  R8VEC_CONVOLUTION computes the convolution\n";
  cout << "  of two vectors.\n";

  r8vec_print ( m, x, "  The factor X:" );
  r8vec_print ( n, y, "  The factor Y:" );

  z = r8vec_convolution ( m, x, n, y );

  r8vec_print ( m + n - 1, z, "  The convolution z = x star y:" );

  r8vec_print ( m + n - 1, z_correct, "  Correct answer:" );

  delete [] z;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST115 tests R8VEC_CONVOLUTION_CIRC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double x[N] = { 1.0, 2.0, 3.0, 4.0 };
  double y[N] = { 1.0, 2.0, 4.0, 8.0 };
  double *z;
  double z_correct[N] = { 37.0, 44.0, 43.0, 26.0 };

  cout << "\n";
  cout << "TEST115\n";
  cout << "  R8VEC_CONVOLUTION_CIRC computes the circular convolution\n";
  cout << "  of two vectors.\n";

  r8vec_print ( N, x, "  The factor X:" );
  r8vec_print ( N, y, "  The factor Y:" );

  z = r8vec_convolution_circ ( N, x, y );

  r8vec_print ( N, z, "  The circular convolution z = xCCy:" );

  r8vec_print ( N, z_correct, "  Correct answer:" );

  delete [] z;

  return;
# undef N
}
//****************************************************************************80

void test116 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST116 tests R8VEC_DIF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *cof;
  double fdif;
  double h = 0.01;
  int i;
  int n = 4;
  double x = 1.0;
  double xi;

  cout << "\n";
  cout << "TEST116\n";
  cout << "  R8VEC_DIF estimates derivatives.\n";
  cout << "\n";
  cout << "  Estimate the derivative of order N = " << n << "\n";
  cout << "  Using H = " << h << "\n";
  cout << "  at argument X = " << x << "\n";
//
//  Get the coefficients.
//
  cof = r8vec_dif ( n, h );

  r8vec_print ( n+1, cof, "  The difference coefficients:" );

  fdif = 0.0;
  for ( i = 0; i <= n; i++ )
  {
    xi = x + ( double ) ( 2 * i - n ) * h;
    fdif = fdif + cof[i] * test116_f ( xi );
  }

  cout << "\n";
  cout << "  Estimate is FDIF = " << fdif << "\n";

  delete [] cof;

  return;
}
//****************************************************************************80

double test116_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    TEST116_F evaluates the function used in TEST116.
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
  double value;

  value = exp ( x );

  return value;
}
//****************************************************************************80

void test1165 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1165 tests R8VEC_DIRECT_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int factor_num = 3;
  int point_num = 24;

  int factor_index;
  int factor_order;
  double *factor_value;
  int i;
  int j;
  double x[factor_num*point_num];

  cout << "\n";
  cout << "TEST1165\n";
  cout << "  R8VEC_DIRECT_PRODUCT forms the entries of a\n";
  cout << "  direct product of a given number of R8VEC factors.\n";

  for ( j = 0; j < point_num; j++ )
  {
    for ( i = 0; i < factor_num; i++ )
    {
      x[i+j*factor_num] = 0.0;
    }
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = new double[factor_order];
      factor_value[0] = 1.0;
      factor_value[1] = 2.0;
      factor_value[2] = 3.0;
      factor_value[3] = 4.0;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = new double[factor_order];
      factor_value[0] = 50.0;
      factor_value[1] = 60.0;
      factor_value[2] = 70.0;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = new double[factor_order];
      factor_value[0] = 800.0;
      factor_value[1] = 900.0;
    }

    r8vec_direct_product ( factor_index, factor_order, factor_value,
      factor_num, point_num, x );

    delete [] factor_value;
  }

  cout << "\n";
  cout << "     J         X(1)      X(2)      X(3)\n";
  cout << "\n";

  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j << "  ";
    for ( i = 0; i < factor_num; i++ )
    {
      cout << "  " << setw(8) << x[i+j*factor_num];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test1166 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1166 tests R8VEC_DIRECT_PRODUCT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int factor_num = 3;
  int point_num = 24;

  int factor_index;
  int factor_order;
  double *factor_value;
  int i;
  int j;
  double w[point_num];

  cout << "\n";
  cout << "TEST1166\n";
  cout << "  R8VEC_DIRECT_PRODUCT2 forms the entries of a\n";
  cout << "  direct product of a given number of R8VEC factors.\n";

  for ( j = 0; j  < point_num; j++ )
  {
    w[j] = 1.0;
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = new double[factor_order];
      factor_value[0] = 2.0;
      factor_value[1] = 3.0;
      factor_value[2] = 5.0;
      factor_value[3] = 7.0;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = new double[factor_order];
      factor_value[0] = 11.0;
      factor_value[1] = 13.0;
      factor_value[2] = 17.0;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = new double[factor_order];
      factor_value[0] = 19.0;
      factor_value[1] = 21.0;
    }

    r8vec_direct_product2 ( factor_index, factor_order, factor_value,
      factor_num, point_num, w );

    delete [] factor_value;
  }

  cout << "\n";
  cout << "     J         W(J)\n";
  cout << "\n";

  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j << "  "
         << "  " << setw(8) << w[j] << "\n";
  }

  return;
}
//****************************************************************************80

void test117 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST117 tests R8VEC_EVEN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *x;
  double xhi = 99.0;
  double xlo = 0.0;

  cout << "\n";
  cout << "TEST117\n";
  cout << "  R8VEC_EVEN computes N evenly spaced values\n";
  cout << "  between XLO and XHI.\n";
  cout << "\n";
  cout << "  XLO = " << xlo << "\n";
  cout << "  XHI = " << xhi << "\n";
  cout << "  while N = " << N << "\n";

  x = r8vec_even_new ( N, xlo, xhi );

  r8vec_print ( N, x, "  Resulting array:" );

  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test118 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST118 tests R8VEC_EVEN2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define NOLD 5
# define MAXVAL 20

  int i;
  int istar;
  int jstar;
  int nfill[NOLD-1] = { 4, 3, 5, 0 };
  int nval;
  double xold[NOLD] = { 0.0, 1.0, 5.0, 2.0, 0.0 };
  double xval[MAXVAL];

  cout << "\n";
  cout << "TEST118:\n";
  cout << "  R8VEC_EVEN2 interpolates a specified number of\n";
  cout << "  points pairs of values in a vector.\n";
  cout << "\n";
  cout << "  Input data:\n";
  cout << "\n";
  for ( i = 1; i <= NOLD; i++ )
  {
    cout << "  " << setw(12) << xold[i-1] << "\n";
    if ( i < NOLD )
    {
      cout << "  (" << nfill[i-1] << ")\n";
    }
  }

  r8vec_even2 ( MAXVAL, nfill, NOLD, xold, &nval, xval );

  cout << "\n";
  cout << "  Resulting vector:\n";
  cout << "\n";

  istar = 1;
  jstar = 1;
  for ( i = 1; i <= nval; i++ )
  {
    if ( i == istar )
    {
      cout << "  " << '*' << "  " << xval[i-1] << "\n";

      if ( jstar < NOLD )
      {
        istar = istar + nfill[jstar-1] + 1;
        jstar = jstar + 1;
      }
    }
    else
    {
      cout << "     " << setw(12) << xval[i-1] << "\n";
    }
  }

  return;
# undef MAXVAL
# undef NOLD
}
//****************************************************************************80

void test120 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST120 tests R8VEC_EXPAND_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int fat = 3;
  int nfat;
  double x[N] = { 16.0, 4.0, 0.0, 4.0, 16.0, 36.0 };
  double *xfat;

  cout << "\n";
  cout << "TEST120\n";
  cout << "  R8VEC_EXPAND_LINEAR linearly interpolates new data\n";
  cout << "  between old values.\n";
  cout << "\n";

  r8vec_print ( N, x, "  Original vector:" );

  cout << "\n";
  cout << "  Expansion factor is " << fat << "\n";

  xfat = r8vec_expand_linear ( N, x, fat );

  nfat = ( N - 1 ) * ( fat + 1 ) + 1;

  r8vec_print ( nfat, xfat, "  Fattened vector:" );

  delete [] xfat;

  return;
# undef N
}
//****************************************************************************80

void test121 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST121 tests R8VEC_FRAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double afrac;
  int k;
  int seed;

  cout << "\n";
  cout << "TEST121\n";
  cout << "  R8VEC_FRAC: K-th smallest R8VEC entry;\n";

  seed = 123456789;

  a = r8vec_uniform_01_new ( N, seed );

  r8vec_print ( N, a, "  Array to search:" );

  cout << "\n";
  cout << "  Fractile  Value\n";
  cout << "\n";

  for ( k = 1; k < N; k = k + N/2 )
  {
    afrac = r8vec_frac ( N, a, k );
    cout << "  " << setw(6) << k
         << "  " << setw(14) << afrac << "\n";
  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test122 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST122 tests R8VEC_HISTOGRAM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define HISTO_NUM 20
# define N 1000

  double *a;
  double a_hi;
  double a_lo;
  double bin_hi;
  double bin_lo;
  int *histo_gram;
  int i;
  int seed = 123456789;
  int test;
  int test_num = 2;

  cout << "\n";
  cout << "TEST122\n";
  cout << "  R8VEC_HISTOGRAM histograms a real vector.\n";

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      cout << "\n";
      cout << "  Uniform data:\n";

      a_lo =  0.0;
      a_hi = +1.0;
      a = r8vec_uniform_01_new ( N, seed );
    }
    else if ( test == 2 )
    {
      cout << "\n";
      cout << "  Normal data:\n";
      a_lo = -3.0;
      a_hi = +3.0;
      a = r8vec_normal_01_new ( N, seed );
    }

    histo_gram = r8vec_histogram ( N, a, a_lo, a_hi, HISTO_NUM );

    cout << "\n";
    cout << "  Histogram of data:\n";
    cout << "\n";

    for ( i = 0; i <= HISTO_NUM+1; i++ )
    {
      if ( i == 0 )
      {
        cout << "  " << "          "
             << "  " << setw(10) << a_lo
             << "  " << setw(6)  << histo_gram[i] << "\n";
      }
      else if ( i <= HISTO_NUM )
      {
        bin_lo = ( ( double ) ( HISTO_NUM - i + 1 ) * a_lo
                 + ( double ) (             i - 1 ) * a_hi )
                 / ( double ) ( HISTO_NUM         );

        bin_hi = ( ( double ) ( HISTO_NUM - i     ) * a_lo
                 + ( double ) (             i     ) * a_hi )
                 / ( double ) ( HISTO_NUM         );

        cout << "  " << setw(10) << bin_lo
             << "  " << setw(10) << bin_hi
             << "  " << setw(6)  << histo_gram[i] << "\n";
      }
      else if ( i == HISTO_NUM+1 )
      {
        cout << "  " << setw(10) << a_hi
             << "  " << "          "
             << "  " << setw(6)  << histo_gram[i] << "\n";
      }
    }
    delete [] a;
    delete [] histo_gram;
  }

  return;
# undef HISTO_NUM
# undef N
}
//****************************************************************************80

void test123 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST123 tests R8VEC_INDEX_INSERT, R8VEC_INDEX_DELETE_ALL, R8VEC_INDEX_DELETE_DUPES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "TEST123\n";
  cout << "  R8VEC_INDEX_INSERT inserts values into an\n";
  cout << "  index sorted array.\n";
  cout << "  R8VEC_INDEX_DELETE_ALL deletes all copies of a\n";
  cout << "  particular value.\n";
  cout << "  R8VEC_INDEX_DELETE_ONE deletes one copies of a\n";
  cout << "  particular value.\n";
  cout << "  R8VEC_INDEX_DELETE_DUPES deletes duplicates.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";
  cout << "\n";

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "  " << xval << "\n";
    r8vec_index_insert ( &n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call R8VEC_INDEX_DELETE_ONE to delete one value of 8:\n";

  xval = 8.0;
  r8vec_index_delete_one ( n, x, indx, xval, &n, x, indx );

  cout << "\n";
  cout << "  Call R8VEC_INDEX_DELETE_ALL to delete all values of 7:\n";

  xval = 7.0;
  r8vec_index_delete_all ( n, x, indx, xval, &n, x, indx );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call R8VEC_INDEX_DELETE_DUPES to delete duplicates:\n";

  r8vec_index_delete_dupes ( n, x, indx, &n, x, indx );

  cout << "\n";
  cout << "  Indexed list of unique entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test124 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST124 tests R8VEC_INDEX_INSERT_UNIQUE and R8VEC_INDEX_ORDER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  int i;
  int indx[N_MAX];
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "TEST124\n";
  cout << "  R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n";
  cout << "  index sorted array.\n";
  cout << "  R8VEC_INDEX_ORDER sorts an index sorted array.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "  " << setw(6) << xval << "\n";
    r8vec_index_insert_unique ( &n, x, indx, xval );
  }

  cout << "\n";
  cout << "  Indexed list of unique entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Now call R8VEC_INDEX_ORDER to carry out the sorting:\n";

  r8vec_index_order ( n, x, indx );

  r8vec_print ( n, x, "  X:" );

  return;
# undef N_MAX
}
//****************************************************************************80

void test125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST125 tests R8VEC_INDEX_INSERT_UNIQUE and R8VEC_INDEX_SEARCH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double b;
  double c;
  int equal;
  int i;
  int indx[N_MAX];
  int less;
  int more;
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  cout << "\n";
  cout << "TEST125\n";
  cout << "  R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n";
  cout << "  index sorted array.\n";
  cout << "  R8VEC_INDEX_SEARCH searches for an entry \n";
  cout << "  with a given value.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";
  cout << "\n";

  b = 0.0;
  c = 20.0;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( b, c, seed );
    xval = ( double ) ( r8_nint ( xval ) );
    cout << "    " << setw(6) << xval << "\n";
    r8vec_index_insert_unique ( &n, x, indx, xval );
  }

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(6) << x[i]
         << "  " << setw(6) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Results of search for given XVAL:\n";
  cout << "\n";
  cout << "  XVAL  Less Equal More\n";
  cout << "\n";

  for ( i = 0; i <= N_MAX; i++ )
  {
    xval = ( double ) ( i );
    r8vec_index_search ( n, x, indx, xval, &less, &equal, &more );
    cout << "  " << setw(6) << xval
         << "  " << setw(3) << less
         << "  " << setw(3) << equal
         << "  " << setw(3) << more << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test1251 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1251 tests R8VEC_INDEX_SORTED_RANGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_hi;
  int i_lo;
  int *indx;
  int n = 20;
  double r[20];
  double r_lo;
  double r_hi;
  int seed;
  double t;
  int test;

  cout << "\n";
  cout << "TEST1251\n";
  cout << "  R8VEC_INDEX_SORTED_RANGE seeks the range I_LO:I_HI\n";
  cout << "  of entries of sorted indexed R so that\n";
  cout << "  R_LO <= R(INDX(I)) <= R_HI for I_LO <= I <= I_HI.\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r8vec_uniform_01 ( n, seed, r );

    r8vec_print ( n, r, "  Array" );

    indx = r8vec_sort_heap_index_a_new ( n, r );

    cout << "\n";
    cout << "     I  INDX    R(INDX(I))\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(4) << i
           << "  " << setw(4) << indx[i]
           << "  " << setw(14) << r[indx[i]] << "\n";
    }

    r_lo = r8_uniform_01 ( seed );
    r_hi = r8_uniform_01 ( seed );

    if ( r_hi < r_lo )
    {
      t = r_lo;
      r_lo = r_hi;
      r_hi = t;
    }

    r8vec_index_sorted_range ( n, r, indx, r_lo, r_hi, &i_lo, &i_hi );

    cout << "\n";
    if ( i_hi < i_lo )
    {
      cout << "  " << "R_LO" << "      "
           << "  " << setw(14) << r_lo << "\n";
      cout << "  " << "R_HI" << "      "
           << "  " << setw(14) << r_hi << "\n";
      cout << "  Empty range in R.\n";
    }
    else
    {
      cout << "  " << "R_LO" << "      "
           << "  " << setw(14) << r_lo << "\n";
      for ( i = i_lo; i <= i_hi; i++ )
      {
        cout << "  " << setw(4) << i
             << "  " << setw(4) << indx[i]
             << "  " << setw(14) << r[indx[i]] << "\n";
      }
      cout << "  " << "R_HI" << "      "
           << "  " << setw(14) << r_hi << "\n";
    }
    delete [] indx;
  }

  return;
}
//****************************************************************************80

void test1252 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1252 tests R8VEC_INDEXED_HEAP_D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double a[20] = {
    101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0,
    111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0 };
  int i;
  int indx[10] = {
    0, 10, 16, 4, 6, 12, 14, 2, 18, 8 };
  int m = 20;
  int n = 10;

  cout << "\n";
  cout << "TEST1252\n";
  cout << "  R8VEC_INDEXED_HEAP_D creates a descending heap\n";
  cout << "  from an indexed vector.\n";
//
//  Print before.
//
  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Heap the data.
//
  r8vec_indexed_heap_d ( n, a, indx );
//
//  Print afterwards.  Only INDX should change.
//
  r8vec_print ( m, a, "  The data vector (should NOT change):" );
  i4vec_print ( n, indx, "  The index vector (may change):" );
  cout << "\n";
  cout << "  A(INDX) is now a descending heap:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }

  return;
}
//****************************************************************************80

void test1255 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1255 tests R8VEC_INDEXED_HEAP_D_EXTRACT and related routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  cout << "\n";
  cout << "TEST1255\n";
  cout << "  For an indexed R8VEC,\n";
  cout << "  R8VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.\n";
  cout << "  R8VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;\n";
  cout << "  R8VEC_INDEXED_HEAP_D_MAX reports the maximum value.\n";
  cout << "\n";
  cout << "  These 3 operations are enough to model a priority queue.\n";
//
//  Set the data array.  To keep things easy, we will use the indicator vector.
//
  a = r8vec_indicator_new ( m );
//
//  The index array will initially be a random subset of the numbers 1 to M,
//  in random order.
//
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Create a descending heap from the indexed array.
//
  r8vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  cout << "\n";
  cout << "  A(INDX) after heaping:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Insert five entries, and monitor the maximum.
//
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    cout << "\n";
    cout << "  Inserting value " << a[indx_insert] << "\n";

    r8vec_indexed_heap_d_insert ( &n, a, indx, indx_insert );

    indx_max = r8vec_indexed_heap_d_max ( n, a, indx );

    cout << "  Current maximum is " << a[indx_max] << "\n";
  }
  r8vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  cout << "\n";
  cout << "  A(INDX) after insertions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }
//
//  Extract the first 5 largest elements.
//
  cout << "\n";
  cout << "  Now extract the maximum several times.\n";
  cout << "\n";

  for ( i = 0; i < 5; i++ )
  {
    indx_extract = r8vec_indexed_heap_d_extract ( &n, a, indx );
    cout << "  Extracting maximum element A[" << indx_extract
         << "] = " << a[indx_extract] << "\n";
  }
  r8vec_print ( m, a, "  The data vector after extractions:" );
  i4vec_print ( n, indx, "  The index vector after extractions:" );
  cout << "\n";
  cout << "  A(INDX) after extractions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << a[indx[i]] << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void test1256 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1256 tests R8VEC_LEGENDRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double *r;
  double r1;
  double r2;

  cout << "\n";
  cout << "TEST1256\n";
  cout << "  R8VEC_LEGENDRE computes N Legendre points in [R1,R2].\n";

  r1 = -1.0;
  r2 = +1.0;
  n = 5;

  r = r8vec_legendre_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n
       << "  R1 = " << r1
       << "  R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Legendre points:" );

  delete [] r;

  r1 =   0.0;
  r2 = +10.0;
  n = 7;

  r = r8vec_legendre_new ( n, r1, r2 );

  cout << "\n";
  cout << "  N = " << n
       << "  R1 = " << r1
       << "  R2 = " << r2 << "\n";

  r8vec_print ( n, r, "  Legendre points:" );

  delete [] r;

  return;
}
//****************************************************************************80

void test1258 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1258 tests R8VEC_LINSPACE_NEW and R8VEC_MIDSPACE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int n = 5;
  double *x;

  cout << "\n";
  cout << "TEST1258\n";
  cout << "  For a R8VEC:\n";
  cout << "  R8VEC_LINSPACE_NEW: evenly spaced points between A and B;\n";
  cout << "  R8VEC_MIDSPACE_NEW: evenly spaced midpoints between A and B\n";

  a = 10.0;
  b = 20.0;

  x = r8vec_linspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_linspace ( 5, 10, 20 )" );
  delete [] x;

  x = r8vec_midspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_midspace ( 5, 10, 20 )" );
  delete [] x;

  return;
}
//****************************************************************************80

void test126 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST126 tests R8VEC_MAX and R8VEC_MIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  int i;
  double rmax;
  double rmin;
  int seed;

  cout << "\n";
  cout << "TEST126\n";
  cout << "  R8VEC_MAX produces the maximum entry in a real array.\n";
  cout << "  R8VEC_MIN produces the minimum entry.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = r8vec_uniform_01_new ( N, seed );

  r8vec_print ( N, a, "  The array:" );

  rmax = r8vec_max ( N, a );
  rmin = r8vec_min ( N, a );

  cout << "\n";
  cout << "  R8VEC_MAX reports the maximum value is " << rmax << ".\n";
  cout << "  R8VEC_MIN reports the minimum value is " << rmin << ".\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test127 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST127 tests R8VEC_MAX_INDEX and R8VEC_MIN_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int ival;
  int seed;

  cout << "\n";
  cout << "TEST127\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_MAX_INDEX: index of maximum entry;\n";
  cout << "  R8VEC_MIN_INDEX: index of minimum entry;\n";

  b = - ( double ) ( N );
  c =   ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  ival = r8vec_max_index ( N, a );
  cout << "  Maximum index:           " << ival << "\n";

  ival = r8vec_min_index ( N, a );
  cout << "  Minimum index:           " << ival << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test128 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST128 tests R8VEC_MEAN and R8VEC_MEDIAN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  double mean;
  double median;
  int seed;

  cout << "\n";
  cout << "TEST128\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_MEAN:      mean value;\n";
  cout << "  R8VEC_MEDIAN:    median value;\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  mean = r8vec_mean ( N, a );
  cout << "  Mean:    " << mean << "\n";
  median = r8vec_median ( N, a );
  cout << "  Median:  " << median << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test129 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST129 tests R8VEC_NORM_L1, R8VEC_NORM_L2, R8VEC_NORM_LI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "TEST129\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_NORM_L1:   L1 norm.\n";
  cout << "  R8VEC_NORM_L2:   L2 norm.\n";
  cout << "  R8VEC_NORM_LI:   L-infinity norm.\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  cout << "\n";
  cout << "  L1 norm:           " << r8vec_norm_l1 ( N, a ) << "\n";
  cout << "  L2 norm:           " << r8vec_norm_l2 ( N, a ) << "\n";
  cout << "  L-Infinity norm:   " << r8vec_norm_li ( N, a ) << "\n";;

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test130 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST130 tests R8VEC_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 1000

  int i;
  int n;
  int seed = 123456789;
  double *x;
  double x_max;
  double x_mean;
  double x_min;
  double x_var;

  cout << "\n";
  cout << "TEST130\n";
  cout << "  R8VEC_NORMAL_01 computes a vector of normally\n";
  cout << "  distributed random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
//
//  Test 1:
//  Simply call 5 times for 1 value, and print.
//
  cout << "\n";
  cout << "  Test #1: Call 5 times, 1 value each time.\n";
  cout << "\n";

  n = 1;
  for ( i = 0; i < 5; i++ )
  {
    x = r8vec_normal_01_new ( n, seed );
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[0] << "\n";
    delete [] x;
  }
//
//  Test 2:
//  Restore the random number seed, and repeat.
//
  cout << "\n";
  cout << "  Test #2: Restore the random number seed.\n";
  cout << "  Call 5 times, 1 value each time.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  n = -1;
  x = r8vec_normal_01_new ( n, seed );
  delete [] x;

  seed = 123456789;

  n = 1;
  for ( i = 0; i < 5; i++ )
  {
    x = r8vec_normal_01_new ( n, seed );
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[0] << "\n";
    delete [] x;
  }
//
//  Test 3:
//  Restore the random number seed, compute all 5 values at once.
//
  cout << "\n";
  cout << "  Test #3: Restore the random number seed.\n";
  cout << "  Call 1 time for 5 values.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  n = -1;
  x = r8vec_normal_01_new ( n, seed );
  delete [] x;

  seed = 123456789;

  n = 5;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;
//
//  Test 4:
//  Restore the random number seed, compute all 5 values at once.
//
  cout << "\n";
  cout << "  Test #4: Restore the random number seed.\n";
  cout << "  Call for 2, 1, and 2 values.\n";
  cout << "  The results should be identical.\n";
  cout << "\n";

  n = -1;
  x = r8vec_normal_01_new ( n, seed );
  delete [] x;

  seed = 123456789;

  n = 2;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;

  n = 1;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;

  n = 2;
  x = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  delete [] x;
//
//  Test 5:
//  Determine the minimum, maximum, mean and variance.
//
  n = N_MAX;
  x = r8vec_normal_01_new ( n, seed );
  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  x_mean = r8vec_mean ( n, x );
  x_var = r8vec_variance ( n, x );
  delete [] x;

  cout << "\n";
  cout << "  Test #5:\n";
  cout << "  Number of samples was " << n << "\n";
  cout << "  Minimum value was " << x_min << "\n";
  cout << "  Maximum value was " << x_max << "\n";
  cout << "  Average value was " << x_mean << "\n";
  cout << "  Variance was      " << x_var << "\n";
  cout << "  Expected average  = 0.0\n";
  cout << "  Expected variance = 1.0\n";

  return;
# undef N_MAX
}
//****************************************************************************80

void test152 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST152 tests R8VEC_NORMALIZE_L1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "TEST152\n";
  cout << "  For an R8VEC:\n";
  cout << "  R8VEC_NORMALIZE_L1:  make unit sum;\n";

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Input vector:" );

  r8vec_normalize_l1 ( N, a );

  r8vec_print ( N, a, "  After calling R8VEC_NORMALIZE_L1:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test131 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST131 tests R8VEC_ORDER_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define TEST_NUM 6

  int itest;
  int j;
  int order;
  double x[N];

  cout << "\n";
  cout << "TEST131\n";
  cout << "  R8VEC_ORDER_TYPE classifies a real vector as\n";
  cout << "  -1: no order\n";
  cout << "   0: all equal;\n";
  cout << "   1: ascending;\n";
  cout << "   2: strictly ascending;\n";
  cout << "   3: descending;\n";
  cout << "   4: strictly descending.\n";
  cout << "\n";

  for ( itest = 1; itest <= TEST_NUM; itest++ )
  {
    if ( itest == 1 )
    {
      x[0] = 1.0;
      x[1] = 3.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 2 )
    {
      x[0] = 2.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 2.0;
    }
    else if ( itest == 3 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 4 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 3.0;
      x[3] = 4.0;
    }
    else if ( itest == 5 )
    {
      x[0] = 4.0;
      x[1] = 4.0;
      x[2] = 3.0;
      x[3] = 1.0;
    }
    else if ( itest == 6 )
    {
      x[0] = 9.0;
      x[1] = 7.0;
      x[2] = 3.0;
      x[3] = 0.0;
    }

    order = r8vec_order_type ( N, x );

    cout << "\n";
    cout << "The following vector has order type " << order << ".\n";
    cout << "\n";

    r8vec_print ( N, x, "" );
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test132 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST132 tests R8VEC_PERMUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int base = 1;
  int i;
  int perm[N] = { 2, 4, 5, 1, 3 };
  double x[N] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

  cout << "\n";
  cout << "TEST132\n";
  cout << "  R8VEC_PERMUTE permutes an R8VEC in place.\n";
  cout << "\n";
  cout << "  I, Perm(I), X(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6)  << i+1
         << "  " << setw(6)  << perm[i]
         << "  " << setw(12) << x[i] << "\n";
  }

  r8vec_permute ( N, perm, base, x );

  r8vec_print ( N, x, "  Permuted array:" );

  return;
# undef N
}
//****************************************************************************80

void test133 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST133 tests R8VEC_POLARIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N] = { 1.0, 2.0,  3.0 };
  double a2[N];
  double a_normal[N];
  double a_parallel[N];
  double ap_norm;
  int i;
  double p[N] = { 3.0, 1.0, -2.0 };
  double p_norm;
  double pan;
  double pap;

  cout << "\n";
  cout << "TEST133\n";
  cout << "  R8VEC_POLARIZE decomposes a vector into\n";
  cout << "  components parallel and normal to a direction.\n";

  r8vec_print ( N, a, "  Original vector:" );

  r8vec_print ( N, p, "  Direction vector:" );

  r8vec_polarize ( N, a, p, a_normal, a_parallel );

  r8vec_print ( N, a_normal, "  Normal component:" );

  r8vec_print ( N, a_parallel, "  Parallel component:" );

  pan = r8vec_dot_product ( N, p, a_normal );
  p_norm = r8vec_norm ( N, p );
  ap_norm = r8vec_norm ( N, a_parallel );

  pap = r8vec_dot_product ( N, p, a_parallel ) / ( p_norm * ap_norm );

  cout << "\n";
  cout << "  Dot product of P and A_normal (should be 0) " << pan << "\n";
  cout << "  Cosine of angle between P and A_parallel (should be 1 or -1) "
       << pap << "\n";

  for ( i = 0; i < N; i++ )
  {
    a2[i] = a_normal[i] + a_parallel[i];
  }

  r8vec_print ( N, a2, "  Sum of components (should equal A):" );

  return;
# undef N
}
//****************************************************************************80

void test134 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST134 tests R8VEC_ROTATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  int m = 2;

  cout << "\n";
  cout << "TEST134\n";
  cout << "  R8VEC_ROTATE rotates an R8VEC in place.\n";
  cout << "\n";
  cout << "  Rotate entries " << m << " places to the right.\n";

  r8vec_print ( N, a, "  Original array:" );

  r8vec_rotate ( N, a, m );

  r8vec_print ( N, a, "  Rotated array:" );

  return;
# undef N
}
//****************************************************************************80

void test135 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST135 tests R8VEC_REVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;

  cout << "\n";
  cout << "TEST135\n";
  cout << "  R8VEC_REVERSE reverses an R8VEC.\n";

  a = r8vec_indicator_new ( N );

  r8vec_print ( N, a, "  Original array:" );

  r8vec_reverse ( N, a );

  r8vec_print ( N, a, "  Reversed array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test136 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST136 tests R8VEC_SEARCH_BINARY_A;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  int index;
  int nc;
  double search_val;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST136\n";
  cout << "  For ascending order:\n";
  cout << "  R8VEC_SEARCH_BINARY_A searches a sorted array;\n";
  cout << "\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8vec_uniform_01_new ( N, seed );

  search_val = a[0];

  r8vec_sort_heap_a ( N, a );

  r8vec_print ( N, a, "  Sorted vector A:" );
//
//  Now search the sorted array for a given value.
//
  cout << "\n";
  cout << "  Search the array for the value " << search_val << "\n";

  index = r8vec_search_binary_a ( N, a, search_val );

  cout << "\n";
  cout << "  SEARCH RESULT:\n";
  cout << "\n";

  if ( 0 < index )
  {
    cout << "    The value occurs in index " << index << "\n";
  }
  else
  {
    cout << "    The value does not occur in the array.\n";
  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test137 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST137 tests R8VEC_SORT_BUBBLE_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST137\n";
  cout << "  R8VEC_SORT_BUBBLE_A sorts a real array.\n";

  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = r8vec_uniform_01_new ( N, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_bubble_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test138 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST138 tests R8VEC_SORT_HEAP_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "TEST138\n";
  cout << "  R8VEC_SORT_HEAP_A ascending sorts an R8VEC.\n";
  cout << "\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Original array:" );

  r8vec_sort_heap_a ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Ascending sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test139 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST139 tests R8VEC_SORT_HEAP_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "TEST139\n";
  cout << "  R8VEC_SORT_HEAP_D descending sorts an R8VEC.\n";
  cout << "\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Original array:" );

  r8vec_sort_heap_d ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Descending sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test140 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST140 tests R8VEC_SORT_HEAP_INDEX_A_NEW and R8VEC_SORT_HEAP_INDEX_D_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  int base;
  double c;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "TEST140\n";
  cout << "  R8VEC_SORT_HEAP_INDEX_A_NEW creates an ascending\n";
  cout << "  sort index for an R8VEC.\n";
  cout << "  R8VEC_SORT_HEAP_INDEX_D_NEW creates a descending\n";
  cout << "  sort index for an R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  indx = r8vec_sort_heap_index_a_new ( N, a );

  cout << "\n";
  cout << "  After indexed ascending sort:\n";
  cout << "\n";
  cout << "         I INDX(I)       A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "   INDX(I)  A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[indx[i]] << "\n";
  }
  cout << "\n";
  cout << "  Call R8VEC_PERMUTE to carry out the permutation explicitly.\n";
  cout << "\n";

  base = 0;
  r8vec_permute ( N, indx, base, a );

  r8vec_print ( N, a, "  I, A(I)" );

  delete [] indx;

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  indx = r8vec_sort_heap_index_d_new ( N, a );

  cout << "\n";
  cout << "  After indexed descending sort:\n";
  cout << "\n";
  cout << "         I  INDX(I)  A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "   INDX(I)  ARRAY(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << indx[i]
         << "  " << setw(12) << a[indx[i]] << "\n";
  }

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void test141 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST141 tests R8VEC_SORT_HEAP_MASK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define MASK_NUM 10
# define N 20

  double *a;
  double b;
  int base;
  double c;
  int i;
  int *indx;
  int mask[MASK_NUM] = { 2, 4, 7, 8, 9, 12, 13, 16, 18, 19 };
  int seed;

  cout << "\n";
  cout << "TEST141\n";
  cout << "  R8VEC_SORT_HEAP_MASK_A creates an ascending\n";
  cout << "  sort index for a masked R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  i4vec_print ( MASK_NUM, mask, "  The mask array:" );

  r8vec_mask_print ( N, a, MASK_NUM, mask, "  The masked unsorted array:" );

  indx = r8vec_sort_heap_mask_a ( N, a, MASK_NUM, mask );

  cout << "\n";
  cout << "  After masked indexed ascending sort:\n";
  cout << "\n";
  cout << "  I, INDX(I), MASK(INDX(I)), A(MASK(INDX(I)))\n";
  cout << "\n";
  for ( i = 0; i < MASK_NUM; i++ )
  {
    cout << "  " << setw(6) << i+1
         << "  " << setw(6) << indx[i]
         << "  " << setw(6) << mask[indx[i]-1]
         << "  " << setw(14) << a[mask[indx[i]-1]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call I4VEC_PERMUTE to carry out the index permutation\n";
  cout << "  explicitly on the MASK vector.\n";
  cout << "\n";

  base = 1;
  i4vec_permute ( MASK_NUM, indx, base, mask );
//
//  Essentially, INDX becomes the identity vector now.
//
  delete [] indx;

  indx = i4vec_indicator_new ( MASK_NUM );

  i4vec_print ( MASK_NUM, mask, "  The reordered mask array:" );

  r8vec_mask_print ( N, a, MASK_NUM, mask,
    "  The reordered masked sorted array:" );

  delete [] a;
  delete [] indx;

  return;
# undef MASK_NUM
# undef N
}
//****************************************************************************80

void test142 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST142 tests R8VEC_SORT_INSERT_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "TEST142\n";
  cout << "  R8VEC_SORT_INSERT_A ascending sorts an R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Unsorted array:" );

  r8vec_sort_insert_a ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test143 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST143 tests R8VEC_SORT_INSERT_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  int base = 1;
  double c;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "TEST143\n";
  cout << "  R8VEC_SORT_INSERT_INDEX_A creates an ascending\n";
  cout << "  sort index for an R8VEC.\n";

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print_some ( N, a, 1, 10, "  Unsorted array:" );

  indx = r8vec_sort_insert_index_a ( N, a );

  cout << "\n";
  cout << "  After indexed ascending sort:\n";
  cout << "\n";
  cout << "  I, INDX(I), A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i+1
         << "  " << setw(6) << indx[i]
         << "  " << setw(12) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "  I, INDX(I), A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(6) << i+1
         << "  " << setw(6) << indx[i]
         << "  " << setw(12) << a[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call R8VEC_PERMUTE to carry out the permutation explicitly.\n";
  cout << "\n";

  r8vec_permute ( N, indx, base, a );

  r8vec_print_some ( N, a, 1, 10, "  Permuted data" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void test144 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST144 tests R8VEC_SORT_QUICK_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  cout << "\n";
  cout << "TEST144\n";
  cout << "  R8VEC_SORT_QUICK_A sorts an R8VEC\n";
  cout << "  using quick sort.\n";

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_quick_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test145 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST145 tests R8VEC_SORTED_MERGE_A;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double *c;
  int na = 10;
  int nb = 10;
  int nc;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST145\n";
  cout << "  For ascending order:\n";
  cout << "  R8VEC_SORTED_MERGE_A merges two sorted R8VEC's;\n";
  cout << "\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = r8vec_uniform_01_new ( na, seed );
  b = r8vec_uniform_01_new ( nb, seed );

  r8vec_sort_heap_a ( na, a );

  r8vec_sort_heap_a ( nb, b );

  r8vec_print ( na, a, "  Sorted vector A:" );

  r8vec_print ( nb, b, "  Sorted vector B:" );

  c = r8vec_sorted_merge_a ( na, a, nb, b, &nc );

  r8vec_print ( nc, c, "  Merged vector C:" );

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void test146 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST146 tests R8VEC_SORTED_NEAREST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double b;
  double c;
  int i;
  int j;
  int seed = 123456789;
  double *x;
  double xval;

  cout << "\n";
  cout << "TEST146\n";
  cout << "  R8VEC_SORTED_NEAREST finds the nearest entry\n";
  cout << "  in a sorted real array.\n";

  b = 0.0;
  c = 10.0;

  x = r8vec_uniform_ab_new ( N, b, c, seed );

  r8vec_sort_heap_a ( N, x );

  r8vec_print ( N, x, "  Sorted array:" );

  cout << "\n";
  cout << "     Test        Nearest\n";
  cout << "     Value    Index   Value\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    xval = r8_uniform_01 ( seed );

    j = r8vec_sorted_nearest ( N, x, xval );

    cout << "  "
         << setw(8) << xval   << "    "
         << setw(6) << j      << "  "
         << setw(8) << x[j-1] << "\n";
  }

  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test1465 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1465 tests R8VEC_SORTED_RANGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_hi;
  int i_lo;
  int n = 10;
  double r[10];
  double r_lo;
  double r_hi;
  int seed;
  double t;
  int test;

  cout << "\n";
  cout << "TEST1465\n";
  cout << "  R8VEC_SORTED_RANGE seeks the range of indices\n";
  cout << "  in a sorted vector R so that\n";
  cout << "  R_LO <= R(I_LO:I_HI) <= R_HI.\n";

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r8vec_uniform_01 ( n, seed, r );

    r8vec_sort_heap_a ( n, r );

    r8vec_print ( n, r, "  Sorted array R:" );

    r_lo = r8_uniform_01 ( seed );
    r_hi = r8_uniform_01 ( seed );

    if ( r_hi < r_lo )
    {
      t = r_lo;
      r_lo = r_hi;
      r_hi = t;
    }

    r8vec_sorted_range ( n, r, r_lo, r_hi, &i_lo, &i_hi );

    cout << "\n";
    if ( i_hi < i_lo )
    {
      cout << "  R_LO  " << setw(14) << r_lo << "\n";
      cout << "  R_HI  " << setw(14) << r_hi << "\n";
      cout << "  Empty range in R.\n";
    }
    else
    {
      cout << "  R_LO  " << setw(14) << r_lo << "\n";
      for ( i = i_lo; i <= i_hi; i++)
      {
        cout << "  " << setw(4) << i << "  " << setw(14) << r[i] << "\n";
      }
      cout << "  R_HI  " << setw(14) << r_hi << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test147 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST147 tests R8VEC_SORTED_SPLIT and R8VEC_SPLIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b;
  double c;
  int i;
  int i_gt;
  int i_lt;
  int isplit;
  int n = 25;
  int seed;
  double split;

  cout << "\n";
  cout << "TEST147\n";
  cout << "  R8VEC_SORTED_SPLIT splits a sorted vector into\n";
  cout << "  entries less than and greater than a\n";
  cout << "  splitting value.\n";
  cout << "  R8VEC_SPLIT splits an unsorted vector\n";
  cout << "  in the same way.\n";
  cout << "\n";

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, b, c, seed );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.5 * ( double ) ( r8_nint ( a[i] ) );
  }

  r8vec_sort_heap_a ( n, a );

  split = 0.5 * ( a[0] + a[n-1] );

  r8vec_print ( n, a, "  The sorted array:" );

  cout << "\n";
  cout << "  Splitting value is " << split << "\n";
  cout << "\n";

  r8vec_sorted_split ( n, a, split, &i_lt, &i_gt );

  cout << "  Lower index I_LT = " << i_lt << "\n";
  cout << "  Upper index I_GT = " << i_gt << "\n";

  cout << "\n";
  cout << "  Now repeat test with R8VEC_SPLIT.\n";
  cout << "\n";

  r8vec_permute_uniform ( n, a, seed );

  r8vec_print ( n, a, "  The shuffled array:" );

  isplit = r8vec_split ( n, a, split );

  r8vec_print ( n, a, "  The split array:" );

  cout << "\n";
  cout << "  Array entries <= SPLIT up to index " << isplit << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void test1475 ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8LIB_TEST1475 tests R8VEC_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2008
//
//  Author:
//
//    John Burkardt
//
{
# define X_NUM 9

  int i;
  double tol;
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  double x_val[X_NUM] = { 11.0, 11.0, 11.0, 22.0, 22.0, 33.0, 33.0, 55.0, 55.0 };
  int *xdnu;
  double *xu_val;

  cout << "\n";
  cout << "R8LIB_TEST1475\n";
  cout << "  R8VEC_SORTED_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique elements of a sorted R8VEC,\n";
  cout << "  and a map from the original vector to the (implicit)\n";
  cout << "  vector of sorted unique elements.\n";

  r8vec_print ( x_num, x_val, "  The vector X:" );

  tol = r8_epsilon ( );
  x_unique_num = r8vec_sorted_unique_count ( x_num, x_val, tol );

  undx = new int[x_unique_num];
  xu_val = new double[x_unique_num];

  xdnu = new int[x_num];

  cout << "\n";
  cout << "  Number of unique entries in X is " << x_unique_num << "\n";

  r8vec_sorted_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu );

  cout << "\n";
  cout << "  UNDX can be used to list the unique elements of X\n";
  cout << "  in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX   X(UNDX)\n";
  cout << "\n";

  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << x_val[undx[i]] << "\n";
  }

  for ( i = 0; i < x_unique_num; i++ )
  {
    xu_val[i] = x_val[undx[i]];
  }
  cout << "\n";
  cout << "  UNDX can be used to created XU, a copy of X\n";
  cout << "  containing only the unique elements, in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX     XU(I)\n";
  cout << "\n";
  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << xu_val[i] << "\n";
  }

  cout << "\n";
  cout << "  XDNU can be used to match each element of X with one of the\n";
  cout << "  unique elements\n";
  cout << "\n";
  cout << "     I  XDNU  X(I)   XU(XDNU(I))\n";
  cout << "\n";

  for ( i = 0; i < x_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << x_val[i]
         << "  " << setw(12) << xu_val[xdnu[i]] << "\n";
  }

  delete [] undx;
  delete [] xdnu;
  delete [] xu_val;

  return;
# undef X_NUM
}
//****************************************************************************80

void test148 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST148 tests R8VEC_SORTED_UNIQUE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double *a;
  double *a_nint;
  double *a_unique;
  double b;
  double c;
  int i;
  int seed;
  double tol = 0.25;
  int unique_num;

  cout << "\n";
  cout << "TEST148\n";
  cout << "  R8VEC_SORTED_UNIQUE finds unique entries in a sorted R8VEC;\n";

  b = 0.0;
  c = ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  a_nint = r8vec_nint ( N, a );

  r8vec_print_some ( N, a_nint, 1, 10, "  Unsorted array:" );

  r8vec_sort_heap_a ( N, a_nint );

  a_unique = r8vec_sorted_unique ( N, a_nint, tol, &unique_num );

  r8vec_print ( unique_num, a_unique, "  Unique entries" );

  delete [] a;
  delete [] a_nint;
  delete [] a_unique;

  return;
# undef N
}
//****************************************************************************80

void test149 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST149 tests R8VEC_SORTED_UNIQUE_COUNT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 30

  double *a;
  double b;
  double c;
  int i;
  int unique_num;
  int seed;
  double tol = 0.25;

  cout << "\n";
  cout << "TEST149\n";
  cout << "  R8VEC_SORTED_UNIQUE_COUNT counts unique entries in a sorted R8VEC;\n";

  b = 0.0;
  c = ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  for ( i = 0; i < N; i++ )
  {
    a[i] = ( double ) r8_nint ( a[i] );
  }

  unique_num = r8vec_sorted_unique_count ( N, a, tol );

  cout << "\n";
  cout << "  Using a tolerance of " << tol << "\n";
  cout << "  R8VEC_SORTED_UNIQUE_COUNT counts " << unique_num
       << " unique entries in A.\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test150 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST150 tests R8VEC_SORTED_UNIQUE_HIST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MAXUNIQ 30
# define N 30

  double *a;
  int acount[MAXUNIQ];
  double auniq[MAXUNIQ];
  double b;
  double c;
  int i;
  int unique_num;
  int seed;
  double tol = 0.25;

  cout << "\n";
  cout << "TEST150\n";
  cout << "  R8VEC_SORTED_UNIQUE_HIST stores the unique entries\n";
  cout << "  and their multiplicities.\n";

  b = 0.0;
  c = ( double ) N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = r8vec_uniform_ab_new ( N, b, c, seed );

  for ( i = 0; i < N; i++ )
  {
    a[i] = ( ( double ) ( ( int ) a[i] ) ) + 0.5;
  }

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_bubble_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  r8vec_sorted_unique_hist ( N, a, tol, MAXUNIQ, &unique_num, auniq, acount );

  cout << "\n";
  cout << "  R8VEC_SORTED_UNIQUE_HIST counts " << unique_num << " unique entries.\n";
  cout << "\n";
  cout << "  Value  Multiplicity\n";
  cout << "\n";
  for ( i = 0; i < unique_num; i++ )
  {
    cout << setw(6)  << i         << "  "
         << setw(12) << auniq[i]  << "  "
         << setw(6)  << acount[i] << "\n";
  }

  delete [] a;

  return;
# undef MAXUNIQ
# undef N
}
//****************************************************************************80

void test1504 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R8LIB_TEST1504 tests R8VEC_TRANSPOSE_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  int n = 12;
  int seed;
  double *x;

  seed = 123456789;

  cout << "\n";
  cout << "R8LIB_TEST1504\n";
  cout << "  R8VEC_TRANSPOSE_PRINT prints an R8VEC \"tranposed\",\n";
  cout << "  that is, placing multiple entries on a line.\n";

  x = r8vec_uniform_01_new ( n, seed );

  r8vec_transpose_print ( n, x, "  The vector X:" );

  delete [] x;

  return;
}
//****************************************************************************80

void test1505 ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8LIB_TEST1505 tests R8VEC_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define X_NUM 9

  int i;
  double tol;
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  double x_val[X_NUM] = { 33.0, 55.0, 11.0, 11.0, 55.0, 33.0, 22.0, 22.0, 11.0 };
  int *xdnu;
  double *xu_val;

  cout << "\n";
  cout << "R8LIB_TEST1505\n";
  cout << "  R8VEC_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique elements of an (unsorted) R8VEC,\n";
  cout << "  and a map from the original vector to the (implicit)\n";
  cout << "  vector of sorted unique elements.\n";

  r8vec_print ( x_num, x_val, "  The vector X:" );

  tol = r8_epsilon ( );
  x_unique_num = r8vec_unique_count ( x_num, x_val, tol );

  undx = new int[x_unique_num];
  xu_val = new double[x_unique_num];

  xdnu = new int[x_num];

  cout << "\n";
  cout << "  Number of unique entries in X is " << x_unique_num << "\n";

  r8vec_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu );

  cout << "\n";
  cout << "  UNDX can be used to list the unique elements of X\n";
  cout << "  in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX   X(UNDX)\n";
  cout << "\n";

  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << x_val[undx[i]] << "\n";
  }

  for ( i = 0; i < x_unique_num; i++ )
  {
    xu_val[i] = x_val[undx[i]];
  }
  cout << "\n";
  cout << "  UNDX can be used to created XU, a copy of X\n";
  cout << "  containing only the unique elements, in sorted order.\n";
  cout << "\n";
  cout << "     I  UNDX     XU(I)\n";
  cout << "\n";
  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(8) << xu_val[i] << "\n";
  }

  cout << "\n";
  cout << "  XDNU can be used to match each element of X with one of the\n";
  cout << "  unique elements\n";
  cout << "\n";
  cout << "     I  XDNU  X(I)   XU(XDNU(I))\n";
  cout << "\n";

  for ( i = 0; i < x_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << xdnu[i]
         << "  " << setw(4) << x_val[i]
         << "  " << setw(12) << xu_val[xdnu[i]] << "\n";
  }

  delete [] undx;
  delete [] xdnu;
  delete [] xu_val;

  return;
# undef X_NUM
}
//****************************************************************************80

void test151 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST151 tests R8VEC_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  double b = 10.0;
  double c = 20.0;
  int j;
  double *r;
  int seed;

  cout << "\n";
  cout << "TEST151\n";
  cout << "  R8VEC_UNIFORM returns a random real vector\n";
  cout << "  with entries in a given range [ B, C ]\n";
  cout << "\n";
  cout << "  For this problem:\n";
  cout << "  B = " << b << "\n";
  cout << "  C = " << c << "\n";
  cout << "\n";

  seed = 123456789;

  for ( j = 1; j <= 3; j++ )
  {
    cout << "\n";
    cout << "  Input SEED = " << seed << "\n";
    cout << "\n";

    r = r8vec_uniform_ab_new ( N, b, c, seed );

    r8vec_print_some ( N, r, 1, 10, "  Random vector:" );

    delete [] r;
  }

  return;
# undef N
}
//****************************************************************************80

void test153 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST153 tests R8VEC2_SORT_A and R8VEC2_SORT_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "TEST153\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SORT_A ascending sorts;\n";
  cout << "  R8VEC2_SORT_D descending sorts;\n";

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_a ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after ascending sort:" );

  r8vec2_sort_d ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after descending sort:" );

  delete [] a1;
  delete [] a2;

  return;
}
//****************************************************************************80

void test154 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST154 tests R8VEC2_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int base = 0;
  int i;
  int *indx;
  int seed = 123456789;
  double x[N];
  double y[N];

  cout << "\n";
  cout << "TEST154\n";
  cout << "  R8VEC2_SORT_HEAP_INDEX_A creates a sort index\n";
  cout << "  for an (X,Y) array.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) i4_uniform_ab ( 0, N, seed ) / ( double ) N;
    y[i] = ( double ) i4_uniform_ab ( 0, N, seed ) / ( double ) N;
  }

  cout << "\n";
  cout << "  The unsorted array:\n";
  cout << "\n";
  cout << "         I  X(I), Y(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  indx = r8vec2_sort_heap_index_a ( N, base, x, y );

  cout << "\n";
  cout << "  After sorting:\n";
  cout << "\n";
  cout << "         I  INDX(I), X(I), Y(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << indx[i]
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  cout << "\n";
  cout << "  Now use the index array to carry out the\n";
  cout << "  permutation implicitly.\n";
  cout << "\n";
  cout << "         I  INDX(I), X(INDX(I)), Y(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << indx[i]
         << "  " << setw(12) << x[indx[i]]
         << "  " << setw(12) << y[indx[i]] << "\n";
  }

  cout << "\n";
  cout << "  R8VEC_PERMUTE carries out the permutation.\n";

  r8vec_permute ( N, indx, base, x );
  r8vec_permute ( N, indx, base, y );

  cout << "\n";
  cout << "         I X(I), Y(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(12) << x[i]
         << "  " << setw(12) << y[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test155 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST155 tests R8VEC2_SORTED_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;
  int unique_num;

  cout << "\n";
  cout << "TEST155\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SORTED_UNIQUE counts unique entries.\n";

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_a ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after ascending sort:" );

  r8vec2_sorted_unique ( n, a1, a2, &unique_num );

  r8vec2_print ( unique_num, a1, a2, "  UNIQed array:" );

  delete [] a1;
  delete [] a2;

  return;
}
//****************************************************************************80

void test156 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST156 tests R8VEC2_SORTED_UNIQUE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a1;
  double *a2;
  double b;
  double c;
  int indx[N];
  int seed;
  int unique_num;

  cout << "\n";
  cout << "TEST156\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SORTED_UNIQUE_INDEX indexes unique entries.\n";

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( N, b, c, seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( N, b, c, seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( N, a1, a2, "  The pair of arrays:" );

  r8vec2_sorted_unique_index ( N, a1, a2, &unique_num, indx );

  cout << "\n";
  cout << "  The number of unique elements is " << unique_num << "\n";

  i4vec_print ( unique_num, indx, "  Index of Unique Elements:" );

  r8vec_index_order ( unique_num, a1, indx );
  r8vec_index_order ( unique_num, a2, indx );

  r8vec2_print ( unique_num, a1, a2, "  After Indexed Nonunique Deletion." );

  delete [] a1;
  delete [] a2;

  return;
# undef N
}
//****************************************************************************80

void test157 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST157 tests R8VEC2_SUM_MAX_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *a1;
  double *a2;
  double b;
  double c;
  int ival;
  int n = 10;
  int seed;

  cout << "\n";
  cout << "TEST157\n";
  cout << "  For a pair of R8VEC's:\n";
  cout << "  R8VEC2_SUM_MAX_INDEX: index of the sum vector\n";
  cout << "  with maximum value.\n";

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, seed );

  b = 0.0;
  c = 5.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, seed );

  r8vec2_print ( n, a1, a2, "  The pair of vectors:" );

  ival = r8vec2_sum_max_index ( n, a1, a2 );

  cout << "\n";
  cout << "  Index of maximum in A+B: " << ival << "\n";

  delete [] a1;
  delete [] a2;

  return;
}
//****************************************************************************80

void test158 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST158 tests R8VECS_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int m = 5;
  int na = 15;

  double a[15] = {
    11.0, 12.0, 13.0, 
    21.0, 22.0, 
    31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 
    41.0, 42.0, 
    51.0 };

  int nvec[6] = { 1, 4, 6, 13, 15, 16 };

  cout << "\n";
  cout << "TEST158\n";
  cout << "  R8VECS_PRINT prints a packed R8VEC.\n";

  r8vecs_print ( m, nvec, na, a, "  Packed R8VEC:" );

  return;
}
