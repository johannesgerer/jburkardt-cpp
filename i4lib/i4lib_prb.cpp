# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "i4lib.hpp"

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
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test243 ( );
void test245 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test29 ( );

void test30 ( );
void test31 ( );
void test32 ( );
void test33 ( );
void test335 ( );
void test34 ( );
void test35 ( );
void test36 ( );
void test37 ( );
void test38 ( );
void test39 ( );

void test40 ( );
void test41 ( );
void test42 ( );
void test43 ( );
void test44 ( );
void test45 ( );
void test46 ( );
void test47 ( );
void test48 ( );
void test49 ( );

void test50 ( );
void test51 ( );
void test52 ( );
void test53 ( );
void test54 ( );
void test55 ( );
void test56 ( );
void test57 ( );
void test58 ( );
void test59 ( );

void test60 ( );
void test602 ( );
void test605 ( );
void test61 ( );
void test62 ( );
void test63 ( );
void test64 ( );
void test65 ( );
void test66 ( );
void test67 ( );
void test68 ( );
void test69 ( );

void test70 ( );
void test71 ( );
void test72 ( );
void test73 ( );
void test74 ( );
void test75 ( );
void test76 ( );
void test77 ( );
void test78 ( );
void test79 ( );

void test80 ( );
void test81 ( );
void test82 ( );
void test83 ( );
void test84 ( );
void test85 ( );
void test86 ( );
void test87 ( );
void test88 ( );
void test89 ( );

void test90 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for I4LIB_PRB.
//
//  Discussion:
//
//    I4LIB_PRB calls the I4LIB tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "I4LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the I4LIB library.\n";

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
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test243 ( );
  test245 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test335 ( );
  test34 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );

  test40 ( );
  test41 ( );
  test42 ( );
  test43 ( );
  test44 ( );
  test45 ( );
  test46 ( );
  test47 ( );
  test48 ( );
  test49 ( );

  test50 ( );
  test51 ( );
  test52 ( );
  test53 ( );
  test54 ( );
  test55 ( );
  test56 ( );
  test57 ( );
  test58 ( );
  test59 ( );

  test60 ( );
  test602 ( );
  test605 ( );
  test61 ( );
  test62 ( );
  test63 ( );
  test64 ( );
  test65 ( );
  test66 ( );
  test67 ( );
  test68 ( );
  test69 ( );

  test70 ( );
  test71 ( );
  test72 ( );
  test73 ( );
  test74 ( );
  test75 ( );
  test76 ( );
  test77 ( );
  test78 ( );
  test79 ( );

  test80 ( );
  test81 ( );
  test82 ( );
  test83 ( );
  test84 ( );
  test85 ( );
  test86 ( );
  test87 ( );
  test88 ( );
  test89 ( );

  test90 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "I4LIB_PRB\n";
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
//    TEST01 tests I4_BIT_HI1.
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
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  I4_BIT_HI1 returns the location of the high 1 bit.\n";
  cout << "\n";
  cout << "       I  I4_BIT_HI1(I)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 0, 100, seed );
    j = i4_bit_hi1 ( i );
    cout << "  " << setw(6) << i
         << "  " << setw(6) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests I4_BIT_LO0.
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
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  I4_BIT_LO0 returns the location of the low 0 bit.\n";
  cout << "\n";
  cout << "       I  I4_BIT_LO0(I)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 0, 100, seed );
    j = i4_bit_lo0 ( i );
    cout << "  " << setw(6) << i
         << "  " << setw(6) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests I4_BIT_LO1.
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
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  I4_BIT_LO1 returns the location of the low 1 bit.\n";
  cout << "\n";
  cout << "       I  I4_BIT_LO1(I)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 0, 100, seed );
    j = i4_bit_lo1 ( i );
    cout << "  " << setw(6) << i
         << "  " << setw(6) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests I4_BIT_REVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_hi;
  int j;
  int k;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  I4_BIT_REVERSE bit reverses I with respect to 2^J\n";
  cout << "\n";
  cout << "         I         J  I4_BIT_REVERSE(I,J)\n";
  cout << "\n";

  for ( j = 0; j <= 4; j++ )
  {
    i_hi = i4_power ( 2, j ) - 1;
    for ( i = 0; i <= i_hi; i++ )
    {
      k = i4_bit_reverse ( i, j );
      cout << "  " << i
           << "  " << j
           << "  " << k << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests I4_CHARACTERISTIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "TEST05\n";
  cout << "  I4_CHARACTERISTIC computes the characteristic\n";
  cout << "  of an integer Q, which is  \n";
  cout << "    Q if Q is prime;\n";
  cout << "    P, if Q = P**N for some prime P;\n";
  cout << "    0, if Q is negative, 0, 1, or the product of \n";
  cout << "      more than 1 distinct prime.\n";
  cout << "\n";
  cout << "  I, I4_CHARACTERISTIC\n";
  cout << "\n";

  for ( i = 1; i <= 50; i++)
  {
    cout << "  " << setw(2) << i
         << "  " << setw(4) << i4_characteristic ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests I4_DIV_ROUNDED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int a_hi =  100;
  int a_lo = -100;
  int b;
  int b_hi =  10;
  int b_lo = -10;
  double c0;
  int c1;
  int c2;
  int c3;
  int c4;
  int seed;
  int test;
  int test_num = 20;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  I4_DIV_ROUNDED performs rounded integer division.\n";
  cout << "\n";
  cout << "  C0 = ( double ) ( a ) / ( double ) ( b )\n";
  cout << "  C1 = I4_DIV_ROUNDED ( A, B )\n";
  cout << "  C2 = r8_nint ( ( double ) ( a ) / ( double ) ( b ) )\n";
  cout << "  C3 = A / B\n";
  cout << "  C4 = ( int ) ( ( double ) ( a ) / ( double ) ( b ) )\n";
  cout << "\n";
  cout << "  C1 and C2 should be equal;\n";
  cout << "  C3 and C4 should be equal.\n";
  cout << "\n";
  cout << "     A     B           C0         C1    C2      C3    C4\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4_uniform_ab ( a_lo, a_hi, seed );
    b = i4_uniform_ab ( b_lo, b_hi, seed );
    if ( b == 0 )
    {
      b = 7;
    }
    c0 = ( double ) ( a ) / ( double ) ( b );
    c1 = i4_div_rounded ( a, b );
    c2 = r8_nint ( ( double ) ( a ) / ( double ) ( b ) );
    c3 = a / b;
    c4 = ( int ) ( ( double ) ( a ) / ( double ) ( b ) );
    cout << "  " << setw(4) << a
         << "  " << setw(4) << b
         << "  " << setw(14) << c0
         << "  " << setw(4) << c1
         << "  " << setw(4) << c2
         << "  " << setw(4) << c3
         << "  " << setw(4) << c4 << "\n";
  }
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests I4_DIVP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int a_hi =  100;
  int a_lo = -100;
  int b;
  int b_hi =  10;
  int b_lo = -10;
  int c;
  int d;
  int seed;
  int test;
  int test_num = 20;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  I4_DIVP(A,B) returns the smallest multiplier of J\n";
  cout << "  that is less than I\n";
  cout << "\n";
  cout << "     A     B     C     D\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4_uniform_ab ( a_lo, a_hi, seed );
    b = i4_uniform_ab ( b_lo, b_hi, seed );
    if ( b == 0 )
    {
      b = 7;
    }
    c = i4_divp ( a, b );
    d = c * b;
    cout << "  " << setw(4) << a
         << "  " << setw(4) << b
         << "  " << setw(4) << c
         << "  " << setw(4) << d << "\n";
  }

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests I4_GCD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 7

  int i;
  int i_test[TEST_NUM] = { 36, 49, 0, 12, 36, 1, 91 };
  int j;
  int j_test[TEST_NUM] = { 30, -7, 71, 12, 49, 42, 28 };
  int test;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  I4_GCD computes the greatest common factor,\n";
  cout << "\n";
  cout << "     I     J   I4_GCD\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i = i_test[test];
    j = j_test[test];
    cout << "  " << setw(6) << i
         << "  " << setw(6) << j
         << "  " << setw(6) << i4_gcd ( i, j ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests I4_HUGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST09\n";
  cout << "  I4_HUGE returns a huge integer.\n";
  cout << "\n";
  cout << "  I4_HUGE() = " << i4_huge ( ) << "\n";

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests I4_HUGE_NORMALIZER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i4;
  double r8;
  double value;

  i4 = i4_huge ( );
  r8 = i4_huge_normalizer ( );

  cout << "\n";
  cout << "TEST10\n";
  cout << "  I4_HUGE_NORMALIZER returns 1/(I4_HUGE+1).\n";
  cout << "\n";
  cout << "  I4_HUGE() = " << i4 << "\n";
  cout << "  I4_HUGE_NORMALIZER() = " << r8 << "\n";

  value = ( double ) ( i4 ) * r8;

  cout << "\n";
  cout << "  I4_HUGE * I4_HUGE_NORMALIZER = " << value << "\n";

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests I4_IS_PRIME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  I4_IS_PRIME reports whether an integer is prime.\n";
  cout << "\n";
  cout << "  I     I4_IS_PRIME(I)\n";
  cout << "\n";

  for ( i = -2; i <= 25; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(1) << i4_is_prime ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests I4_LCM.
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
# define TEST_NUM 7

  int i;
  int i_test[TEST_NUM] = { 36, 49,  0, 12, 36,  1, 91 };
  int j;
  int j_test[TEST_NUM] = { 30, -7, 71, 12, 49, 42, 28 };
  int test;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  I4_LCM computes the least common multiple.\n";
  cout << "\n";
  cout << "     I     J   I4_LCM\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i = i_test[test];
    j = j_test[test];
    cout << "  " << setw(6) << i
         << "  " << setw(6) << j
         << "  " << setw(6) << i4_lcm ( i, j ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests I4_LOG_10.
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
# define N 13

  int i;
  int x[N] = { 0, 1, 2, 3, 9, 10, 11, 99, 101, -1, -2, -3, -9 };

  cout << "\n";
  cout << "TEST13\n";
  cout << "  I4_LOG_10: whole part of log base 10,\n";
  cout << "\n";
  cout << "  X, I4_LOG_10\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout                                 << "  "
         << setw(6) << x[i]              << "  "
         << setw(6) << i4_log_10 ( x[i] ) << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests I4_LOG_2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 17

  int test;
  int x;
  int x_test[TEST_NUM] = {
      0,    1,    2,    3,    9,
     10,   11,   99,  101,   -1,
     -2,   -3,   -9, 1000, 1023,
   1024, 1025 };

  cout << "\n";
  cout << "TEST14\n";
  cout << "  I4_LOG_2: whole part of log base 2.\n";
  cout << "\n";
  cout << "       X     I_LOG_2\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setw(6) << x
         << "  " << setw(12) << i4_log_2 ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests I4_LOG_I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i4;
  int j4;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  I4_LOG_I4: whole part of log base B,\n";
  cout << "\n";
  cout << "        I4        J4 I4_LOG_I4\n";
  cout << "\n";

  for ( j4 = 2; j4 <= 5; j4++ )
  {
    for ( i4 = 0; i4 <= 10; i4++ )
    {
      cout << "  " << setw(8) << i4
           << "  " << setw(8) << j4
           << "  " << setw(8) << i4_log_i4 ( i4, j4 ) << "\n";
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests I4_LOG_R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 10

  double b;
  double b_test[TEST_NUM] = {
    2.0, 3.0,  4.0,  5.0,   6.0,
    7.0, 8.0, 16.0, 32.0, 256.0 };
  int test;
  int x;

  x = 16;

  cout << "\n";
  cout << "TEST16\n";
  cout << "  I4_LOG_R8: whole part of log base B,\n";
  cout << "\n";
  cout << "  X  B  I4_LOG_R8\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = b_test[test];

    cout << "  " << setw(6) << x
         << "  " << setw(14) << b
         << "  " << setw(12) << i4_log_r8 ( x, b ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests I4_MANT.
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
  int is;
  int j;
  int k;
  int l;
  double x;

  x = -314.159;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  I4_MANT decomposes an integer,\n";
  cout << "\n";
  cout << "  Number to be decomposed is X = " << x << "\n";

  i4_mant ( x, &is, &j, &k, &l );

  cout << "\n";
  cout << "  X = "    << is
       << " * ( "     << j
       << " / "       << k
       << " ) * 2 ^ " << l << "\n";

  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests I4_MODDIV;
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
# define TEST_NUM 4

  int ndivid[TEST_NUM] = { 50, -50, 50, -50 };
  int nmult;
  int nrem;
  int number[TEST_NUM] = { 107, 107, -107, -107 };
  int test;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  I4_MODDIV factors a number\n";
  cout << "  into a multiple and a remainder.\n";
  cout << "\n";
  cout << "    Number   Divisor  Multiple Remainder\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4_moddiv ( number[test], ndivid[test], &nmult, &nrem );

    cout << "  " << setw(10) << number[test]
         << "  " << setw(10) << ndivid[test]
         << "  " << setw(10) << nmult
         << "  " << setw(10) << nrem << "\n";
  }

  cout << "\n";
  cout << "  Repeat using C++ percent\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nrem = ( number[test] % ndivid[test] );
    nmult = number[test] / ndivid[test];

    cout << "  " << setw(10) << number[test]
         << "  " << setw(10) << ndivid[test]
         << "  " << setw(10) << nmult
         << "  " << setw(10) << nrem << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests I4_MODP.
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
# define TEST_NUM 4

  int ndivid[TEST_NUM] = { 50, -50, 50, -50 };
  int nmult;
  int nrem;
  int number[TEST_NUM] = { 107, 107, -107, -107 };
  int test;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  I4_MODP factors a number\n";
  cout << "  into a multiple and a remainder.\n";
  cout << "\n";
  cout << "    Number   Divisor  Multiple Remainder\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nrem = i4_modp ( number[test], ndivid[test] );
    nmult = number[test] / ndivid[test];

    cout << "  " << setw(10) << number[test]
         << "  " << setw(10) << ndivid[test]
         << "  " << setw(10) << nmult
         << "  " << setw(10) << nrem << "\n";
  }

  cout << "\n";
  cout << "  Repeat using C++ percent operator:\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nrem = ( number[test] % ndivid[test] );
    nmult = number[test] / ndivid[test];

    cout << "  " << setw(10) << number[test]
         << "  " << setw(10) << ndivid[test]
         << "  " << setw(10) << nmult
         << "  " << setw(10) << nrem << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests I4_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  int test;
  int x;
  int x_test[TEST_NUM] = { -10, -7, 0, 5, 9 };

  cout << "\n";
  cout << "TEST20\n";
  cout << "  I4_SIGN returns the sign of a number.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setw(6) << x
         << "  " << setw(6) << i4_sign ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests I4_SWAP.
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
  int i;
  int j;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  I4_SWAP swaps two integers.\n";

  i = 1;
  j = 202;

  cout << "\n";
  cout << "  Before swapping: \n";
  cout << "\n";
  cout << "    I = " << i << "\n";
  cout << "    J = " << j << "\n";

  i4_swap ( &i, &j );

  cout << "\n";
  cout << "  After swapping: \n";
  cout << "\n";
  cout << "    I = " << i << "\n";
  cout << "    J = " << j << "\n";

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests I4_WALSH_1D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int w0;
  int wm1;
  int wm2;
  int wm3;
  int wp1;
  int wp2;
  double x;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  I4_WALSH_1D evaluates 1D Walsh functions:\n";
  cout << "\n";
  cout << "X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)\n";
  cout << "\n";

  for ( i = 0; i <= 32; i++ )
  {
    x = ( double ) i / 4.0;

    wp2 = i4_walsh_1d ( x,  2 );
    wp1 = i4_walsh_1d ( x,  1 );
    w0  = i4_walsh_1d ( x,  0 );
    wm1 = i4_walsh_1d ( x, -1 );
    wm2 = i4_walsh_1d ( x, -2 );
    wm3 = i4_walsh_1d ( x, -3 );

    cout << "  " << setw(10) << x
         << "  " << setw(2) << wp2
         << "  " << setw(2) << wp1
         << "  " << setw(2) << w0
         << "  " << setw(2) << wm1
         << "  " << setw(2) << wm2
         << "  " << setw(2) << wm3 << "\n";
  }

  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests I4_WRAP.
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
  int i;
  int ihi = 8;
  int ilo = 4;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  I4_WRAP forces an integer to lie within given limits.\n";
  cout << "\n";
  cout << "  ILO = " << ilo << "\n";
  cout << "  IHI = " << ihi << "\n";
  cout << "\n";
  cout << "     I  I4_WRAP(I)\n";
  cout << "\n";

  for ( i = -10; i <= 20; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << i4_wrap ( i, ilo, ihi ) << "\n";
  }

  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests I4_XOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_lo = 0;
  int i_hi = 100;
  int j;
  int k;
  int l;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  I4_XOR returns the bitwise exclusive OR of\n";
  cout << "  two integers.\n";
  cout << "  The operator ^ should generally be used instead!\n";
  cout << "\n";
  cout << "       I       J  I4_XOR     I^J\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( i_lo, i_hi, seed );
    j = i4_uniform_ab ( i_lo, i_hi, seed );
    k = i4_xor ( i, j );
    l = i ^ j;

    cout << "  " << setw(6) << i
         << "  " << setw(6) << j
         << "  " << setw(6) << k
         << "  " << setw(6) << l << "\n";
  }

  return;
}
//****************************************************************************80

void test243 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST243 tests I4BLOCK_NEW and I4BLOCK_DELETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int ***a;
  int i;
  int j;
  int k;
  int l;
  int m;
  int n;

  cout << "\n";
  cout << "TEST243:\n";
  cout << "  I4BLOCK_NEW dynamically creates a 3D array.\n";
  cout << "  I4BLOCK_DELETE deletes it.\n";
  cout << "  Array entries can be addressed using the\n";
  cout << "  notation \"a[i][j][k]\".\n";
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by I4BLOCK_DELETE and new memory
//  requested by I4BLOCK_NEW.
//
  l = 2;
  m = 3;
  n = 2;
//
//  Allocate memory.
//
  cout << "\n";
  cout << "  Allocating memory for array A of size " << l << " by " << m << " by " << n << ".\n";

  a = i4block_new ( l, m, n );

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
        a[i][j][k] = 100 * i + 10 * j + k;
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
  i4block_delete ( a, l, m, n );

  return;
}
//****************************************************************************80

void test245 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST245 tests I4BLOCK_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int l = 4;
  int m = 3;
  int n = 2;
  int x[4*3*2] = {
        1,  2,  3,   4,  1, 
        4,  9, 16,   1,  8, 
       27, 64,  2,   4,  6, 
        8,  2,  8,  18, 32, 
        2, 16, 54, 128 };

  cout << "\n";
  cout << "TEST245\n";
  cout << "  I4BLOCK_PRINT prints an I4BLOCK.\n";

  i4block_print ( l, m, n, x, "  The 3D array:" );

  return;
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests I4COL_FIND_ITEM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 4
# define TEST_NUM 3

  int a[M*N];
  int col;
  int i;
  int item;
  int item_test[TEST_NUM] = { 34, 12, 90 };
  int j;
  int row;
  int test;

  cout << "\n";
  cout << "TEST25\n";
  cout << "  I4COL_FIND_ITEM finds the first occurrence of\n";
  cout << "  an item in an integer array of columns.\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }

  i4mat_print ( M, N, a, "  The matrix of columns:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    item = item_test[test];

    i4col_find_item ( M, N, a, item, &row, &col );

    cout << "  Item " << item
         << " occurs in row " << row
         << " and column " << col << "\n";
  }

  return;
# undef M
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests I4COL_FIND_PAIR_WRAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 4
# define TEST_NUM 5

  int a[M*N];
  int col;
  int i;
  int item1;
  int item1_test[TEST_NUM] = { 22, 32, 22, 54, 54 };
  int item2;
  int item2_test[TEST_NUM] = { 32, 22, 23, 14, 11 };
  int j;
  int row;
  int test;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  I4COL_FIND_PAIR_WRAP finds the first occurrence of\n";
  cout << "  a pair of item in an integer array of columns.\n";
  cout << "  Items in the array are ordered by column, and\n";
  cout << "  wraparound is allowed.\n";

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }

  i4mat_print ( M, N, a, "  The matrix of columns:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    item1 = item1_test[test];
    item2 = item2_test[test];

    i4col_find_pair_wrap ( M, N, a, item1, item2, &row, &col );

    cout << "  Item " << item1
         << " followed by item " << item2
         << " occurs in row " << row
         << " and column " << col << "\n";
  }

  return;
# undef M
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests I4COL_SORT_A and I4COL_SORT_D.
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
  int *a;
  int b = 1;
  int c = 10;
  int m = 5;
  int n = 4;
  int seed;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  I4COL_SORT_A ascending sorts an integer array\n";
  cout << "  as a table of columns.\n";
  cout << "  I4COL_SORT_D descending sorts an integer array\n";
  cout << "  as a table of columns.\n";

  seed = 123456789;

  a = i4mat_uniform_new ( m, n, b, c, seed );

  i4mat_print ( m, n, a, "  The original matrix:" );

  i4col_sort_a ( m, n, a );

  i4mat_print ( m, n, a, "  Ascending sorted:" );

  i4col_sort_d ( m, n, a );

  i4mat_print ( m, n, a, "  Descending sorted:" );

  delete [] a;

  return;
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests I4COL_SORT2_A;
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
  int *a;
  int b = 0;
  int c = 20;
  int m = 6;
  int n = 4;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST28\n";
  cout << "  For a rectangular integer matrix:\n";
  cout << "  I4COL_SORT2_D sorts the elements of the columns.\n";

  a = i4mat_uniform_new ( m, n, b, c, seed );

  i4mat_print ( m, n, a, "  The matrix:" );

  i4col_sort2_a ( m, n, a );

  i4mat_print ( m, n, a, "  The element-sorted column matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests I4COL_SORTED_SINGLETON_COUNT;
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
  int *a;
  int b;
  int c;
  int m = 3;
  int n = 10;
  int seed;
  int singleton_num;
  int test;
  int test_num = 2;

  cout << "\n";
  cout << "TEST29\n";
  cout << "  I4COL_SORTED_SINGLETON_COUNT counts singletons\n";
  cout << "  in a sorted ICOL;\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    b = 0;
    c = 3;

    a = i4mat_uniform_new ( m, n, b, c, seed );

    i4col_sort_a ( m, n, a );

    i4mat_print ( m, n, a, "  Ascending sorted ICOL:" );

    singleton_num = i4col_sorted_singleton_count ( m, n, a );

    cout << "\n";
    cout << "  Number of singletons = " << singleton_num << "\n";

    delete [] a;
  }

  return;
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests I4COL_SORTED_UNIQUE_COUNT;
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
  int *a;
  int b;
  int c;
  int m = 3;
  int n = 10;
  int seed;
  int unique_num;
  int test;
  int test_num = 2;

  cout << "\n";
  cout << "TEST30\n";
  cout << "  I4COL_SORTED_UNIQUE_COUNT counts the unique entries\n";
  cout << "  of a sorted ICOL;\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    b = 0;
    c = 3;

    a = i4mat_uniform_new ( m, n, b, c, seed );

    i4col_sort_a ( m, n, a );

    i4mat_print ( m, n, a, "  Ascending sorted ICOL:" );

    unique_num = i4col_sorted_unique_count ( m, n, a );

    cout << "\n";
    cout << "  Number of unique entries = " << unique_num << "\n";

    delete [] a;
  }

  return;
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests I4MAT_ELIM and I4MAT_RED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5

  int a[M*N];
  int col[N];
  int factor;
  int i;
  int j;
  int k;
  int row[M];
  int test;
  int test_num = 3;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  I4MAT_ELIM does exact Gauss elimination.\n";
  cout << "  I4MAT_RED divides common factors in a matrix;\n";

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      k = 0;
      for ( i = 0; i < M; i++ )
      {
        for ( j = 0; j < N; j++ )
        {
          k = k + 1;
          a[i+j*M] = k;
        }
      }
    }
    else if ( test == 2 )
    {
      factor = 8 * 7 * 6 * 5 * 4 * 3 * 2;

      for ( i = 0; i < M; i++ )
      {
        for ( j = 0; j < N; j++ )
        {
          a[i+j*M] = factor / ( i + j + 1 );
        }
      }
    }
    else if ( test == 3 )
    {
      for ( i = 0; i < M; i++ )
      {
        for ( j = 0; j < N; j++ )
        {
          a[i+j*M] = ( i + 1 ) * ( j + 1 );
        }
      }
    }

    i4mat_print ( M, N, a, "  The original matrix:" );

    i4mat_red ( M, N, a, row, col );

    cout << "\n";
    cout << "  The matrix, as returned by I4MAT_RED:\n";
    cout << "  (Factors are displayed in an extra row and column.)\n";
    cout << "\n";
    for ( i = 0; i < M; i++ )
    {
      for ( j = 0; j < N; j++ )
      {
        cout << "  " << setw(6) << a[i+j*M];
      }
      cout << "  " << setw(6) << row[i] << "\n";
    }
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(6) << col[j];
    }
    cout << "\n";

    i4mat_elim ( M, N, a );

    i4mat_print ( M, N, a, "  The matrix returned by I4MAT_ELIM:" );
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests I4MAT_MAX_INDEX and I4MAT_MIN_INDEX.
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
  int *a;
  int b = 0;
  int c = 10;
  int i;
  int j;
  int m = 5;
  int n = 7;
  int seed;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  I4MAT_MAX_INDEX locates the maximum;\n";
  cout << "  I4MAT_MIN_INDEX locates the minimum;\n";

  seed = 123456789;

  a = i4mat_uniform_new ( m, n, b, c, seed );

  i4mat_print ( m, n, a, "  Random array:" );

  cout << "\n";
  i4mat_max_index ( m, n, a, &i, &j );
  cout << "  Maximum I,J indices            " << i << "  " << j << "\n";
  i4mat_min_index ( m, n, a, &i, &j );
  cout << "  Minimum I,J indices            " << i << "  " << j << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests I4MAT_L1_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2005
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
  int a[N*N] = {
     1,  2,  0,  5,  0, 75,
     0,  1,  0,  0,  0,  0,
     0,  0,  1,  3,  0,  0,
     0,  0,  0,  1,  0,  6,
     0,  0,  0,  0,  1,  4,
     0,  0,  0,  0,  0,  1 };
  int *b;
  int *c;

  cout << "\n";
  cout << "TEST33\n";
  cout << "  I4MAT_L1_INVERSE inverts a unit lower triangular matrix.\n";

  i4mat_print ( N, N, a, "  The original matrix:" );

  b = i4mat_l1_inverse ( N, a );

  i4mat_print ( N, N, b, "  The inverse matrix:" );

  c = i4mat_mm ( N, N, N, a, b );

  i4mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test335 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST335 tests I4MAT_NEW and I4MAT_DELETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int **a;
  int **b;
  int i;
  int j;
  int k;
  int m;
  int n;

  cout << "\n";
  cout << "TEST335:\n";
  cout << "  I4MAT_NEW dynamically creates a 2D array.\n";
  cout << "  I4MAT_DELETE deletes it.\n";
  cout << "  Array entries can be addressed using the\n";
  cout << "  notation \"a[i][j]\".\n";
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by I4MAT_DELETE and new memory
//  requested by I4MAT_NEW.
//
  m = 4;
  n = 5;
//
//  Allocate memory.
//
  cout << "\n";
  cout << "  Allocating memory for array A of size " << m << " by " << n << ".\n";

  a = i4mat_new ( m, n );

  cout << "\n";
  cout << "  Assigning values to A.\n";
//
//  Store values in A.
//
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = 10 * i + j;
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
  b = i4mat_new ( n, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = 0;
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
  i4mat_delete ( a, m, n );
  i4mat_delete ( b, n, n );

  return;
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests I4MAT_PERM_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N*N];
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST34\n";
  cout << "  I4MAT_PERM_UNIFORM applies a random permutation\n";
  cout << "  to a square integer matrix.\n";

  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = 10 * ( i + 1 ) + j + 1;
    }
  }
  i4mat_print ( N, N, a, "  The original matrix:" );

  i4mat_perm_uniform ( N, a, seed );

  i4mat_print ( N, N, a, "  The permuted matrix:" );

  return;
# undef N
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 tests I4MAT_U1_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2005
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
  int a[N*N] = {
    1,  0,  0,  0,  0,  0,
    2,  1,  0,  0,  0,  0,
    0,  0,  1,  0,  0,  0,
    5,  0,  3,  1,  0,  0,
    0,  0,  0,  0,  1,  0,
   75,  0,  0,  6,  4,  1 };
  int *b;
  int *c;

  cout << "\n";
  cout << "TEST35\n";
  cout << "  I4MAT_U1_INVERSE inverts a unit upper triangular matrix.\n";

  i4mat_print ( N, N, a, "  The original matrix:" );

  b = i4mat_u1_inverse ( N, a );

  i4mat_print ( N, N, b, "  The inverse matrix:" );

  c = i4mat_mm ( N, N, N, a, b );

  i4mat_print ( N, N, c, "  Product C = A * B:" );

  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test36 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST36 tests I4ROW_MAX and I4ROW_MIN;
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
# define M 3
# define N 4

  int a[M*N];
  int *amax;
  int *amin;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST36\n";
  cout << "  I4ROW_MAX computes row maximums;\n";
  cout << "  I4ROW_MIN computes row minimums;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = k;
    }
  }

  i4mat_print ( M, N, a, "  The matrix:" );

  amax = i4row_max ( M, N, a );

  amin = i4row_min ( M, N, a );

  cout << "\n";
  cout << "  Maximum, minimum:\n";
  cout << "\n";

  for ( i = 0; i < M; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(6) << amax[i]
         << "  " << setw(6) << amin[i] << "\n";
  }

  delete [] amax;
  delete [] amin;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test37 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST37 tests I4ROW_MEAN and I4ROW_SUM;
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
# define M 3
# define N 4

  int a[M*N];
  int i;
  int j;
  int k;
  double *mean;
  int *rowsum;

  cout << "\n";
  cout << "TEST37\n";
  cout << "  I4ROW_MEAN computes row means;\n";
  cout << "  I4ROW_SUM computes row sums;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = k;
    }
  }

  i4mat_print ( M, N, a, "  The matrix:" );

  rowsum = i4row_sum ( M, N, a );

  mean = i4row_mean ( M, N, a );

  cout << "\n";
  cout << "  Sum, mean:\n";
  cout << "\n";
  for ( i = 0; i < M; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(6) << rowsum[i]
         << "  " << setw(10) << mean[i] << "\n";
  }

  delete [] rowsum;
  delete [] mean;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test38 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST38 tests I4ROW_SORT_A;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int b = 0;
  int c = 10;
  int m = 10;
  int n = 4;
  int seed;

  cout << "\n";
  cout << "TEST38\n";
  cout << "  For a rectangular integer matrix:\n";
  cout << "  I4ROW_SORT_A sorts the rows;\n";

  seed = 123456789;

  a = i4mat_uniform_new ( m, n, b, c, seed );

  i4mat_print ( m, n, a, "  The original matrix:" );

  i4row_sort_a ( m, n, a );

  i4mat_print ( m, n, a, "  The row-sorted matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests I4ROW_SORT_D and I4ROW_SORT2_D;
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
# define M 6
# define N 4

  int a[M*N];
  int i;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST39\n";
  cout << "  For a rectangular integer matrix:\n";
  cout << "  I4ROW_SORT_D sorts the rows;\n";
  cout << "  I4ROW_SORT2_D sorts the elements of the rows.\n";

  seed = 123456789;

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }

  i4mat_print ( M, N, a, "  The original matrix:" );

  i4mat_perm2_uniform ( M, N, a, seed );

  i4mat_print ( M, N, a, "  The matrix, permuted by I4MAT_PERM2_UNIFORM:" );

  i4row_sort_d ( M, N, a );

  i4mat_print ( M, N, a, "  The row-sorted matrix:" );

  i4row_sort2_d ( M, N, a );

  i4mat_print ( M, N, a, "  The element-sorted row-sorted matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test40 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST40 tests I4ROW_SWAP;
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
# define M 3
# define N 4

  int a[M*N];
  int i;
  int j;
  int k;
  int row1;
  int row2;

  cout << "\n";
  cout << "TEST40\n";
  cout << "  For an integer matrix of rows,\n";
  cout << "  I4ROW_SWAP swaps two rows;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = k;
    }
  }

  i4mat_print ( M, N, a, "  The matrix:" );

  row1 = 1;
  row2 = 3;

  cout << "\n";
  cout << "  Swap rows " << row1 << " and " << row2 << "\n";
  cout << "\n";

  i4row_swap ( M, N, a, row1, row2 );

  i4mat_print ( M, N, a, "  The new matrix:" );

  return;
# undef M
# undef N
}
//****************************************************************************80

void test41 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST41 tests I4ROW_VARIANCE.
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
# define M 3
# define N 4

  int a[M*N];
  int i;
  int j;
  int k;
  double *variance;

  cout << "\n";
  cout << "TEST41\n";
  cout << "  I4ROW_VARIANCE computes row variances;\n";

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = k;
    }
  }

  i4mat_print ( M, N, a, "  The matrix:" );

  variance = i4row_variance ( M, N, a );

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

void test42 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST42 tests I4VEC_AMAX;
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
# define N 10

  int *a;
  int aval;
  int b;
  int c;
  int ival;
  int seed;

  cout << "\n";
  cout << "TEST42\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_AMAX:   maximum absolute entry;\n";

  seed = 123456789;

  b = -N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  aval = i4vec_amax ( N, a );

  cout << "\n";
  cout << "  Maximum absolute value: " << aval << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test43 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST43 tests I4VEC_AMIN;
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
# define N 10

  int *a;
  int aval;
  int b;
  int c;
  int ival;
  int seed;

  cout << "\n";
  cout << "TEST43\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_AMIN:   minimum absolute entry;\n";

  seed = 123456789;

  b = -N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  aval = i4vec_amin ( N, a );

  cout << "\n";
  cout << "  Minimum absolute value: " << aval << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test44 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST44 tests I4VEC_AMINZ and I4VEC_AMINZ_INDEX;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int aval;
  int b;
  int c;
  int ival;
  int seed;

  cout << "\n";
  cout << "TEST44\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_AMINZ:  minimum nonzero absolute entry;\n";
  cout << "  I4VEC_AMINZ_INDEX: index of minimum nonzero absolute entry;\n";

  seed = 123456789;

  b = -N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  aval = i4vec_aminz ( N, a );
  ival = i4vec_aminz_index ( N, a );

  cout << "  Minimum abs nonzero:       " << aval << "\n";
  cout << "  Minimum abs nonzero index: " << ival << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test45 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST45 tests I4VEC_AMAX_INDEX and I4VEC_AMIN_INDEX;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int aval;
  int b;
  int c;
  int ival;
  int seed;

  cout << "\n";
  cout << "TEST45\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_AMAX_INDEX:  index of maximum absolute entry;\n";
  cout << "  I4VEC_AMIN_INDEX:  index minimum absolute entry;\n";

  seed = 123456789;

  b = -N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  ival = i4vec_amax_index ( N, a );

  cout << "  Maximum abs index:        " << ival << "\n";

  ival = i4vec_amin_index ( N, a );

  cout << "  Minimum abs index:	     " << ival << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test46 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST46 tests I4VEC_MAX_INDEX, I4VEC_MAX_INDEX_LAST and I4VEC_INDEX;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int aval;
  int b;
  int c;
  int ival;
  int j;
  int seed;

  cout << "\n";
  cout << "TEST46\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_MAX_INDEX:          a maximal index;\n";
  cout << "  I4VEC_MAX_INDEX_LAST:     last maximal index;\n";
  cout << "  I4VEC_INDEX:              first index of given value;\n";

  seed = 123456789;

  b = -N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );
  aval = a[N/2];

  i4vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  ival = i4vec_max_index ( N, a );
  cout << "  Maximum index:            " << ival<< "\n";
  ival = i4vec_max_index_last ( N, a );
  cout << "  Last maximum index:       " << ival << "\n";
  ival = i4vec_min_index ( N, a );
  cout << "  Minimum index:            " << ival << "\n";
  cout << "\n";
  j = i4vec_index ( N, a, aval );
  cout << "  Index of first occurrence of " << aval << " is " << j << "\n";
  aval = aval + 1;
  j = i4vec_index ( N, a, aval );
  cout << "  Index of first occurrence of " << aval << " is " << j << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test47 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST47 tests I4VEC_ASCEND_SUB
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 14

  int *a;
  int b = 1;
  int c = 10;
  int length;
  int seed = 123456789;
  int *sub;
  int test;
  int test_num = 6;

  cout << "\n";
  cout << "TEST47\n";
  cout << "  I4VEC_ASCEND_SUB computes a longest ascending\n";
  cout << "  subsequence of an integer vector.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4vec_uniform_new ( N, b, c, seed );
    i4vec_print ( N, a, "  The vector to be tested:" );
    sub = i4vec_ascend_sub ( N, a, &length );
    i4vec_print ( length, sub, "  A longest ascending subsequence:" );
    delete [] a;
    delete [] sub;
  }

  return;
# undef N
}
//****************************************************************************80

void test48 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST48 tests I4VEC_ASCENDS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define TEST_NUM 6

  int test;
  int *x;
//
//  Each ROW of this definition is a COLUMN of the matrix.
//
  int x_test[N*TEST_NUM] = {
    1, 3, 2, 4,
    2, 2, 2, 2,
    1, 2, 2, 4,
    1, 2, 3, 4,
    4, 4, 3, 1,
    9, 7, 3, 0 };

  cout << "\n";
  cout << "TEST48\n";
  cout << "  I4VEC_ASCENDS determines if an integer vector ascends.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test + test * N;

    i4vec_print ( N, x, "  Test vector:" );

    cout << "  I4VEC_ASCENDS =  " << i4vec_ascends ( N, x ) << "\n";
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test49 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST49 tests I4VEC_BRACKET and I4VEC_INSERT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20
# define TEST_NUM 6

  int a[N_MAX];
  int atest[TEST_NUM] = {
    -10, 2, 9, 10, 20, 24 };
  int aval;
  int i;
  int left;
  int n;
  int right;
  int test;

  cout << "\n";
  cout << "TEST49\n";
  cout << "  I4VEC_BRACKET finds a pair of entries in a\n";
  cout << "  sorted integer array which bracket a value.\n";
  cout << "  I4VEC_INSERT inserts a value into a vector.\n";
  cout << "\n";
  cout << "  We use these two routines to bracket a value,\n";
  cout << "  and then insert it.\n";

  n = 10;
  for ( i = 0; i < n; i++ )
  {
    a[i] = 2 * ( i + 1 );
  }
  a[5] = a[4];

  i4vec_print ( n, a, "  Sorted array:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    aval = atest[test];

    cout << "\n";
    cout << "  Search for AVAL = " << aval << "\n";

    i4vec_bracket ( n, a, aval, &left, &right );

    cout << "  Left =  " << left  << "\n";
    cout << "  Right = " << right << "\n";

    if ( 1 <= left )
    {
      cout << "  A[LEFT] =  " << a[left-1] << "\n";
    }

    if ( 1 <= right )
    {
      cout << "  A(RIGHT) = " << a[right-1] << "\n";
    }
//
//  Insert the value.
//
    if ( left == -1 )
    {
      left = 0;
    }

    if ( left == right )
    {
      cout << "\n";
      cout << "  No insertion necessary.\n";
    }
    else
    {
      i4vec_insert ( n, a, left+1, aval );
      n = n + 1;
      i4vec_print ( n, a, "  Sorted, augmented array:" );
    }
  }

  return;
# undef N_MAX
# undef TEST_NUM
}
//****************************************************************************80

void test50 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST50 tests I4VEC_CUM_NEW and I4VEC_CUM0_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 December 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int *a_cum;
  int *a_cum0;
  int b;
  int c;
  int seed;

  cout << "\n";
  cout << "TEST50\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_CUM_NEW:  cumulative sum;\n";
  cout << "  I4VEC_CUM0_NEW: cumulative sum, (zero based);\n";

  seed = 123456789;

  b = -N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  cout << "\n";

  a_cum = i4vec_cum_new ( N, a );

  i4vec_print ( N, a_cum, "  Cumulative sums:" );

  a_cum0 = i4vec_cum_new ( N, a );

  i4vec_print ( N+1, a_cum0, "  0-based cumulative sums:" );

  delete [] a;
  delete [] a_cum;
  delete [] a_cum0;

  return;
# undef N
}
//****************************************************************************80

void test51 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST51 tests I4VEC_DESCENDS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define TEST_NUM 6

  int test;
  int *x;
//
//  Each ROW of this definition is a COLUMN of the matrix.
//
  int x_test[N*TEST_NUM] = {
    1, 3, 2, 4,
    2, 2, 2, 2,
    1, 2, 2, 4,
    1, 2, 3, 4,
    4, 4, 3, 1,
    9, 7, 3, 0 };

  cout << "\n";
  cout << "TEST51\n";
  cout << "  I4VEC_DESCENDS determines if an integer vector descends.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test + test * N;

    i4vec_print ( N, x, "  Test vector:" );

    cout << "  I4VEC_DESCENDS = " << i4vec_descends ( N, x ) << "\n";;
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test52 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST52 tests I4VEC_DIRECT_PRODUCT.
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
  int *factor_value;
  int i;
  int j;
  int x[factor_num*point_num];

  cout << "\n";
  cout << "TEST52\n";
  cout << "  I4VEC_DIRECT_PRODUCT forms the entries of a\n";
  cout << "  direct product of a given number of I4VEC factors.\n";

  for ( j = 0; j  < point_num; j++ )
  {
    for ( i = 0; i < factor_num; i++ )
    {
      x[i+j*factor_num] = 0;
    }
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = new int[factor_order];
      factor_value[0] = 1;
      factor_value[1] = 2;
      factor_value[2] = 3;
      factor_value[3] = 4;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = new int[factor_order];
      factor_value[0] = 50;
      factor_value[1] = 60;
      factor_value[2] = 70;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = new int[factor_order];
      factor_value[0] = 800;
      factor_value[1] = 900;
    }

    i4vec_direct_product ( factor_index, factor_order, factor_value,
      factor_num, point_num, x );

    delete [] factor_value;
  }

  cout << "\n";
  cout << "     J     X(1)  X(2)  X(3)\n";
  cout << "\n";

  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j << "  ";
    for ( i = 0; i < factor_num; i++ )
    {
      cout << "  " << setw(4) << x[i+j*factor_num];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test53 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST tests I4VEC_DIRECT_PRODUCT2.
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
  int *factor_value;
  int i;
  int j;
  int w[point_num];

  cout << "\n";
  cout << "TEST53\n";
  cout << "  I4VEC_DIRECT_PRODUCT2 forms the entries of a\n";
  cout << "  direct product of a given number of I4VEC factors.\n";

  for ( j = 0; j  < point_num; j++ )
  {
    w[j] = 1;
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = new int[factor_order];
      factor_value[0] = 2;
      factor_value[1] = 3;
      factor_value[2] = 5;
      factor_value[3] = 7;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = new int[factor_order];
      factor_value[0] = 11;
      factor_value[1] = 13;
      factor_value[2] = 17;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = new int[factor_order];
      factor_value[0] = 19;
      factor_value[1] = 21;
    }

    i4vec_direct_product2 ( factor_index, factor_order, factor_value,
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

void test54 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST54 tests I4VEC_FRAC;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int afrac;
  int b = 1;
  int c = 2 * N;
  int k;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST54\n";
  cout << "  I4VEC_FRAC: K-th smallest integer vector entry.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  The array to search:" );

  cout << "\n";
  cout << "  Fractile    Value\n";
  cout << "\n";

  for ( k = 1; k <= N; k = k + (N/2) )
  {
    afrac = i4vec_frac ( N, a, k );
    cout << "  " << setw(6) << k
         << "  " << setw(6) << afrac << "\n";
  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test55 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST55 tests I4VEC_HEAP_A and I4VEC_HEAP_D;
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
# define N 30

  int *a;
  int b;
  int c;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST55\n";
  cout << "  I4VEC_HEAP_A turns an integer array into an ascending heap;\n";
  cout << "  I4VEC_HEAP_D turns an integer array into a descending heap;\n";

  b = 1;
  c = 40;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unheaped array:" );

  i4vec_heap_a ( N, a );

  i4vec_print ( N, a, "  Ascending heaped array:" );

  delete [] a;

  seed = 123456789;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unheaped array:" );

  i4vec_heap_d ( N, a );

  i4vec_print ( N, a, "  Descending heaped array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test56 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST56 tests I4VEC_HEAP_D_EXTRACT, I4VEC_HEAP_D_INSERT and I4VEC_HEAP_D_MAX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  int a[N_MAX];
  int b;
  int c;
  int i;
  int n;
  int seed;
  int value;

  cout << "\n";
  cout << "TEST56\n";
  cout << "  For a descending heap of integers,\n";
  cout << "  I4VEC_HEAP_D_INSERT inserts a value into the heap.\n";
  cout << "  I4VEC_HEAP_D_EXTRACT extracts the maximum value;\n";
  cout << "  I4VEC_HEAP_D_MAX reports the maximum value.\n";
  cout << "\n";
  cout << "  These 3 operations are enough to model a priority queue.\n";

  n = 0;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    b = 0;
    c = 10;

    value = i4_uniform_ab ( b, c, seed );

    i4vec_heap_d_insert ( &n, a, value );

    cout << "\n";
    cout << "  Inserting value          " << value << "\n";

    value = i4vec_heap_d_max ( n, a );

    cout <<"  Current maximum value is " << value << "\n";
  }

  i4vec_print ( n, a, "  Current heap as a vector:" );

  cout << "\n";
  cout << "  Now extract the maximum several times.\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    value = i4vec_heap_d_extract ( &n, a );
    cout << "  Extracting maximum element = " << value << "\n";
  }

  i4vec_print ( n, a, "  Current heap as a vector:" );

  return;
# undef N_MAX
}
//****************************************************************************80

void test57 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST57 tests I4VEC_HISTOGRAM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 1000

  int *a;
  int *histo_gram;
  int histo_num;
  int i;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST57\n";
  cout << "  I4VEC_HISTOGRAM histograms an integer vector.\n";

  a = i4vec_uniform_new ( N, 0, 25, seed );

  histo_num = 20;

  histo_gram = i4vec_histogram ( N, a, histo_num );

  cout << "\n";
  cout << "  Histogram of data from 0 to " << histo_num << "\n";
  cout << "\n";

  for ( i = 0; i <= histo_num; i++ )
  {
    if ( 0 < histo_gram[i] )
    {
      cout << "  " << setw(6) << i
           << "  " << setw(6) << histo_gram[i] << "\n";
    }
  }

  delete [] a;
  delete [] histo_gram;

  return;
# undef N
}
//****************************************************************************80

void test58 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST58 tests I4VEC_INDEX_INSERT_UNIQUE and I4VEC_INDEX_SEARCH.
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

  int b;
  int c;
  int equal;
  int i;
  int indx[N_MAX];
  int less;
  int more;
  int n;
  int seed;
  int x[N_MAX];
  int xval;

  n = 0;

  cout << "\n";
  cout << "TEST58\n";
  cout << "  I4VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n";
  cout << "  index sorted array.\n";
  cout << "  I4VEC_INDEX_SEARCH searches for an entry with\n";
  cout << "  a given value.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";

  b = 0;
  c = N_MAX;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = i4_uniform_ab ( b, c, seed );
    i4vec_index_insert_unique ( &n, x, indx, xval );
  }
  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(3) << x[i]
         << "  " << setw(3) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Results of search for given XVAL:\n";
  cout << "\n";
  cout << "  XVAL  Less Equal More\n";
  cout << "\n";

  for ( xval = 0; xval <= N_MAX; xval++ )
  {
    i4vec_index_search ( n, x, indx, xval, &less, &equal, &more );
    cout << "  " << setw(3) << xval
         << "  " << setw(3) << less
         << "  " << setw(3) << equal
         << "  " << setw(3) << more << "\n";
  }
  return;
# undef N_MAX
}
//****************************************************************************80

void test59 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST59 tests I4VEC_INDEX_INSERT, I4VEC_INDEX_DELETE_DUPES, I4VEC_INDEX_DELETE_ALL.
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
# define N_MAX 25

  int b;
  int c;
  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  int x[N_MAX];
  int xval;

  n = 0;

  cout << "\n";
  cout << "TEST59\n";
  cout << "  I4VEC_INDEX_INSERT inserts values into an\n";
  cout << "  index sorted array of integers.\n";
  cout << "  I4VEC_INDEX_DELETE_ALL deletes all copies of a\n";
  cout << "  particular value.\n";
  cout << "  I4VEC_INDEX_DELETE_ONE deletes one copies of a\n";
  cout << "  particular value.\n";
  cout << "  I4VEC_INDEX_DELETE_DUPES deletes duplicates.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";
  cout << "\n";

  xval = 8;
  i4vec_index_insert ( &n, x, indx, xval );

  xval = 7;
  i4vec_index_insert ( &n, x, indx, xval );

  b = 0;
  c = 20;
  seed = 123456789;

  for ( i = 0; i < 20; i++ )
  {
    xval = i4_uniform_ab ( b, c, seed );
    cout << "  " << xval << "\n";
    i4vec_index_insert ( &n, x, indx, xval );
  }

  xval = 7;
  i4vec_index_insert ( &n, x, indx, xval );

  xval = 8;
  i4vec_index_insert ( &n, x, indx, xval );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(3) << x[i]
         << "  " << setw(3) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call I4VEC_INDEX_DELETE_ONE to delete a value of 8:\n";

  xval = 8;
  i4vec_index_delete_one ( n, x, indx, xval, &n, x, indx );

  cout << "\n";
  cout << "  Call I4VEC_INDEX_DELETE_ALL to delete values of 7:\n";

  xval = 7;
  i4vec_index_delete_all ( n, x, indx, xval, &n, x, indx );

  cout << "\n";
  cout << "  Indexed list of entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(3) << x[i]
         << "  " << setw(3) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Call I4VEC_INDEX_DELETE_DUPES to delete duplicates:\n";

  i4vec_index_delete_dupes ( n, x, indx, &n, x, indx );

  cout << "\n";
  cout << "  Indexed list of unique entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i+1
         << "  " << setw(3) << indx[i]
         << "  " << setw(3) << x[i] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test60 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST60 tests I4VEC_INDEX_INSERT_UNIQUE and I4VEC_INDEX_ORDER.
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
# define N_MAX 20

  int b;
  int c;
  int i;
  int indx[N_MAX];
  int n;
  int seed;
  int x[N_MAX];
  int xval;

  n = 0;

  cout << "\n";
  cout << "TEST60\n";
  cout << "  I4VEC_INDEX_INSERT_UNIQUE inserts unique values into\n";
  cout << "  an index sorted array.\n";
  cout << "  I4VEC_INDEX_ORDER sorts an index sorted array.\n";
  cout << "\n";
  cout << "  Generate some random values:\n";
  cout << "\n";

  b = 0;
  c = 20;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = i4_uniform_ab ( b, c, seed );
    cout << "  " << setw(3) << xval << "\n";
    i4vec_index_insert_unique ( &n, x, indx, xval );
  }

  cout << "\n";
  cout << "  Indexed list of unique entries:\n";
  cout << "\n";
  cout << "  I  INDX(I)  X(I)  X(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i
         << "  " << setw(3) << indx[i]
         << "  " << setw(3) << x[i]
         << "  " << setw(3) << x[indx[i]-1] << "\n";
  }

  cout << "\n";
  cout << "  Now call I4VEC_INDEX_ORDER to carry out the sorting:\n";

  i4vec_index_order ( n, x, indx );

  i4vec_print ( n, x, "  X:" );

  return;
# undef N_MAX
}
//****************************************************************************80

void test602 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST602 tests I4VEC_INDEXED_HEAP_D;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int a[20] = {
    101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
    111, 112, 113, 114, 115, 116, 117, 118, 119, 120 };
  int i;
  int indx[10] = {
    0, 10, 16, 4, 6, 12, 14, 2, 18, 8 };
  int m = 20;
  int n = 10;

  cout << "\n";
  cout << "TEST602\n";
  cout << "  I4VEC_INDEXED_HEAP_D creates a descending heap\n";
  cout << "  from an indexed vector.\n";
//
//  Print before.
//
  i4vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << a[indx[i]] << "\n";
  }
//
//  Heap the data.
//
  i4vec_indexed_heap_d ( n, a, indx );
//
//  Print afterwards.  Only INDX should change.
//
  i4vec_print ( m, a, "  The data vector (should NOT change):" );
  i4vec_print ( n, indx, "  The index vector (may change):" );
  cout << "\n";
  cout << "  A(INDX) is now a descending heap:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << a[indx[i]] << "\n";
  }

  return;
}
//****************************************************************************80

void test605 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST605 tests I4VEC_INDEXED_HEAP_D_EXTRACT and related routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  cout << "\n";
  cout << "TEST605\n";
  cout << "  For an indexed I4VEC,\n";
  cout << "  I4VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.\n";
  cout << "  I4VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;\n";
  cout << "  I4VEC_INDEXED_HEAP_D_MAX reports the maximum value.\n";
  cout << "\n";
  cout << "  These 3 operations are enough to model a priority queue.\n";
//
//  Set the data array.  To keep things easy, we will use the indicator vector.
//
  a = i4vec_indicator_new ( m );
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

  i4vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  cout << "\n";
  cout << "  A(INDX):\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << a[indx[i]] << "\n";
  }
//
//  Create a descending heap from the indexed array.
//
  i4vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  cout << "\n";
  cout << "  A(INDX) after heaping:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << a[indx[i]] << "\n";
  }
//
//  Insert five entries, and monitor the maximum.
//
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    cout << "\n";
    cout << "  Inserting value " << a[indx_insert] << "\n";

    i4vec_indexed_heap_d_insert ( &n, a, indx, indx_insert );

    indx_max = i4vec_indexed_heap_d_max ( n, a, indx );

    cout << "  Current maximum is " << a[indx_max] << "\n";
  }
  i4vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  cout << "\n";
  cout << "  A(INDX) after insertions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << a[indx[i]] << "\n";
  }
//
//  Extract the first 5 largest elements.
//
  cout << "\n";
  cout << "  Now extract the maximum several times.\n";
  cout << "\n";

  for ( i = 0; i < 5; i++ )
  {
    indx_extract = i4vec_indexed_heap_d_extract ( &n, a, indx );
    cout << "  Extracting maximum element A[" << indx_extract
         << "] = " << a[indx_extract] << "\n";
  }
  i4vec_print ( m, a, "  The data vector after extractions:" );
  i4vec_print ( n, indx, "  The index vector after extractions:" );
  cout << "\n";
  cout << "  A(INDX) after extractions:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << a[indx[i]] << "\n";
  }

  delete [] a;

  return;
}
//****************************************************************************80

void test61 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST61 tests I4VEC_INDICATOR_NEW;
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
# define N 10

  int *a;

  cout << "\n";
  cout << "TEST61\n";
  cout << "  I4VEC_INDICATOR_NEW sets A[0:N-1] = 1...N;\n";

  a = i4vec_indicator_new ( N );

  i4vec_print ( N, a, "  The indicator vector:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test62 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST62 tests I4VEC_MAX and I4VEC_MIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int a_max;
  int a_min;
  int b;
  int c;
  int seed;

  cout << "\n";
  cout << "TEST62\n";
  cout << "  I4VEC_MAX produces the maximum entry in an integer array.\n";
  cout << "  I4VEC_MIN produces the minimum entry.\n";

  b = 1;
  c = 30;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  The array:" );

  a_max = i4vec_max ( N, a );
  a_min = i4vec_min ( N, a );

  cout << "\n";
  cout << "  Maximum " << a_max << ".\n";
  cout << "  Minimum " << a_min << ".\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test63 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST63 tests I4VEC_MEAN and I4VEC_MEDIAN;
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

  int *a;
  int b;
  int c;
  int j;
  double mean;
  int median;
  int seed;

  cout << "\n";
  cout << "TEST63\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_MEAN:          mean value;\n";
  cout << "  I4VEC_MEDIAN:        median value;\n";

  seed = 123456789;

  b = -N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  mean = i4vec_mean ( N, a );
  median = i4vec_median ( N, a );

  cout << "\n";
  cout << "  Mean:    " << mean   << "\n";
  cout << "  Median:  " << median << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test64 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST64 tests I4VEC_MERGE_A;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N1 10
# define N2 10

  int *a1;
  int *a2;
  int *a3;
  int b;
  int c;
  int seed;

  cout << "\n";
  cout << "TEST64\n";
  cout << "  For ascending order:\n";
  cout << "  I4VEC_MERGE_A merges two sorted integer arrays;\n";

  seed = 123456789;

  b = 0;
  c = N1;

  a1 = i4vec_uniform_new ( N1, b, c, seed );

  i4vec_sort_heap_a ( N1, a1 );

  b = 0;
  c = N2;

  a2 = i4vec_uniform_new ( N2, b, c, seed );

  i4vec_sort_heap_a ( N2, a2 );

  i4vec_print ( N1, a1, "  Input vector A1:" );

  i4vec_print ( N2, a2, "  Input vector A2:" );

  a3 = i4vec_merge_a ( N1, a1, N2, a2 );

  i4vec_print ( N1+N2, a3, "  Merged vector A3:" );

  delete [] a1;
  delete [] a2;
  delete [] a3;

  return;
# undef N1
# undef N2
}
//****************************************************************************80

void test65 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST65 tests I4VEC_NONZERO_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int b;
  int c;
  int nonzero;
  int seed;

  cout << "\n";
  cout << "TEST65\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_NONZERO_COUNT: number of nonzeroes;\n";

  seed = 123456789;

  b = -2;
  c = 3;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  nonzero = i4vec_nonzero_count ( N, a );

  cout << "\n";
  cout << "  Number of nonzeroes :     " << nonzero << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test66 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST66 tests I4VEC_NONZERO_FIRST.
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

  int *a;
  int a_save[N];
  int i;
  int ihi = 2;
  int ilo = -1;
  int indx[N];
  int nz;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "TEST66\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_NONZERO_FIRST left shifts the nonzero entries\n";
  cout << "  of an I4VEC so they appear first.\n";
  cout << "\n";
  cout << "  ----------Before--------------    ----------After---------------\n";
  cout << "\n";
  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4vec_uniform_new ( N, ilo, ihi, seed );
    i4vec_copy ( N, a, a_save );
    i4vec_nonzero_first ( N, a, &nz, indx );
    cout << "  ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(3) << a_save[i];
    }
    cout << "    ";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(3) << a[i];
    }
    cout << "\n";
    delete [] a;
  }

  cout << "\n";
  cout << "  The value NZ counts the nonzeros, and\n";
  cout << "  the vector INDX indicates the original positions:\n";
  cout << "\n";

  a = i4vec_uniform_new ( N, ilo, ihi, seed );
  i4vec_copy ( N, a, a_save );
  i4vec_nonzero_first ( N, a, &nz, indx );

  cout << "\n";
  cout << "  Original vector:\n";
  cout << "\n";
  cout << "  ";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(3) << a_save[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "  Number of nonzeros NZ = " << nz << "\n";
  cout << "\n";
  cout << "  Shifted vector:\n";
  cout << "\n";
  cout << "  ";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(3) << a[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "  Index vector:\n";
  cout << "\n";
  cout << "  ";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(3) << indx[i];
  }
  cout << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test67 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST67 tests I4VEC_ORDER_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define TEST_NUM 6

  int j;
  int order;
  int test;
  int *x;
//
//  Each ROW of the definition is a COLUMN of the matrix.
//
  int x_test[N*TEST_NUM] = {
    1, 3, 2, 4,
    2, 2, 2, 2,
    1, 2, 2, 4,
    1, 2, 3, 4,
    4, 4, 3, 1,
    9, 7, 3, 0 };

  cout << "\n";
  cout << "TEST67\n";
  cout << "  I4VEC_ORDER_TYPE classifies an integer vector as\n";
  cout << "  -1: no order\n";
  cout << "   0: all equal;\n";
  cout << "   1: ascending;\n";
  cout << "   2: strictly ascending;\n";
  cout << "   3: descending;\n";
  cout << "   4: strictly descending.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test + N * test;

    order = i4vec_order_type ( N, x );

    cout << "\n";
    cout << "  The following vector has order type " << order << "\n";
    cout << "\n";
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(6) << j+1
           << "  " << setw(6) << x[j] << "\n";
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test68 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST68 tests I4VEC_PAIRWISE_PRIME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define TEST_NUM 6

  int i;
  int test;
  int *x;
//
//  Each ROW of the definition is a COLUMN of the matrix.
//
  int x_test[N*TEST_NUM] = {
     1,  3,  2,  4,
     2,  2,  2,  2,
     5,  7, 12, 29,
     1, 13,  1, 11,
     1,  4,  9, 16,
     6, 35, 13, 77 };

  cout << "\n";
  cout << "TEST68\n";
  cout << "  I4VEC_PAIRWISE_PRIME determines if a vector of\n";
  cout << "  integers is pairwise prime.\n";
  cout << "\n";
  cout << "              Pairwise\n";
  cout << "  Row Vector     Prime?\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test + N * test;

    for ( i = 0; i < N; i++ )
    {
      cout << "  " << setw(3) << x[i];
    }
    cout << "  " << setw(1) << i4vec_pairwise_prime ( N, x ) << "\n";
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test69 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST69 tests I4VEC_PART.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int a[N];
  int nval;

  cout << "\n";
  cout << "TEST69\n";
  cout << "  I4VEC_PART partitions an integer.\n";

  nval = 17;
  cout << "\n";
  cout << "  NVAL = " << nval << "\n";

  i4vec_part ( N, nval, a );

  i4vec_print ( N, a, "  Partitioned:" );

  nval = -49;
  cout << "\n";
  cout << "  NVAL = " << nval << "\n";

  i4vec_part ( N, nval, a );

  i4vec_print ( N, a, "  Partitioned:" );

  return;
# undef N
}
//****************************************************************************80

void test70 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST70 tests I4VEC_PART_QUICK_A.
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
# define N 16

  int *a;
  int b;
  int c;
  int i;
  int l;
  int r;
  int seed;

  cout << "\n";
  cout << "TEST70\n";
  cout << "  I4VEC_PART_QUICK_A reorders an int vector\n";
  cout << "  as part of a quick sort.\n";

  b = 0;
  c = N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Before rearrangement:" );

  i4vec_part_quick_a ( N, a, &l, &r );

  cout << "\n";
  cout << "  Rearranged array\n";
  cout << "  Left = " << l << "\n";
  cout << "  Right = " << r << "\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {

    if ( i == l && 0 < l )
    {
      cout << "\n";
    }

    if ( i == r-1 )
    {
      cout << "\n";
    }

    cout                 << "  "
         << setw(6) << i << "  "
         << setw(6) << a[i] << "\n";

  }

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test71 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST71 tests I4VEC_PERMUTE.
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
# define N 12

  int *a;
  int b = 0;
  int base = 1;
  int c = N;
  int *p;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST71\n";
  cout << "  I4VEC_PERMUTE reorders an integer vector\n";
  cout << "  according to a given permutation.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  A, before rearrangement:" );

  p = perm_uniform_new ( N, base, seed );

  i4vec_print ( N, p, "  Permutation vector P:" );

  i4vec_permute ( N, p, base, a );

  i4vec_print ( N, a, "  A, after rearrangement:" );

  delete [] a;
  delete [] p;

  return;
# undef N
}
//****************************************************************************80

void test72 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST72 tests I4VEC_REVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int b = 0;
  int c = 3 * N;
  int seed;

  cout << "\n";
  cout << "TEST72\n";
  cout << "  I4VEC_REVERSE reverses a list of integers.\n";

  seed = 123456789;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Original vector:" );

  i4vec_reverse ( N, a );

  i4vec_print ( N, a, "  Reversed:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test73 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST73 tests I4VEC_RUN_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int j;
  int n = 20;
  int run_count;
  int seed;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST73\n";
  cout << "  I4VEC_RUN_COUNT counts runs in an I4VEC\n";
  cout << "\n";
  cout << " Run Count        Sequence\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4vec_uniform_new ( n, 0, 1, seed );

    run_count = i4vec_run_count ( n, a );

    cout << "  " << setw(8) << run_count;
    cout << "        ";
    for ( j = 0; j < n; j++ )
    {
      cout << setw(2) << a[j];
    }
    cout << "\n";
    delete [] a;
  }

  return;
}
//****************************************************************************80

void test74 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST74 tests I4VEC_SEARCH_BINARY_A;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int b;
  int c;
  int index;
  int search_val;
  int seed;

  cout << "\n";
  cout << "TEST74\n";
  cout << "  For ascending order:\n";
  cout << "  I4VEC_SEARCH_BINARY_A searchs an array for a value;\n";

  seed = 123456789;

  b = 0;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  search_val = a[0];

  i4vec_sort_heap_a ( N, a );

  i4vec_print ( N, a, "  Input vector A:" );

  cout << "\n";
  cout << "  Search the array A for the value " << search_val << "\n";

  index = i4vec_search_binary_a ( N, a, search_val );

  cout << "\n";
  cout << "  SEARCH RESULT:\n";
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

void test75 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST75 tests I4VEC_SORT_BUBBLE_A and I4VEC_SORTED_UNIQUE_HIST.
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
# define MAXUNIQ 30
# define N 30

  int *a;
  int acount[MAXUNIQ];
  int auniq[MAXUNIQ];
  int b;
  int c;
  int i;
  int unique_num;
  int seed;

  cout << "\n";
  cout << "TEST75\n";
  cout << "  I4VEC_SORT_BUBBLE_A sorts an integer array;\n";
  cout << "  I4VEC_SORTED_UNIQUE_HIST stores the unique entries\n";
  cout << "  and their multiplicities.\n";

  b = 0;
  c = N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  i4vec_sort_bubble_a ( N, a );

  i4vec_print ( N, a, "  Sorted array:" );

  i4vec_sorted_unique_hist ( N, a, MAXUNIQ, &unique_num, auniq, acount );

  cout << "\n";
  cout << "  I4VEC_SORTED_UNIQUE_HIST counts " << unique_num
       << " unique entries.\n";
  cout << "\n";

  i4vec2_print ( unique_num, auniq, acount, "  Value and Multiplicity" );

  delete [] a;

  return;
# undef MAXUNIQ
# undef N
}
//****************************************************************************80

void test76 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST76 tests I4VEC_SORT_HEAP_A;
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
# define N 30

  int *a;
  int b;
  int c;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST76\n";
  cout << "  I4VEC_SORT_HEAP_A sorts an integer array;\n";

  b = 0;
  c = N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  i4vec_sort_heap_a ( N, a );

  i4vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test77 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST77 tests I4VEC_SORT_HEAP_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int b = 0;
  int c = 3 * N;
  int seed;

  cout << "\n";
  cout << "TEST77\n";
  cout << "  For a vector of integers,\n";
  cout << "  I4VEC_SORT_HEAP_D descending sorts.\n";

  seed = 123456789;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted:" );

  i4vec_sort_heap_d ( N, a );

  i4vec_print ( N, a, "  Descending sorted:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test78 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST78 tests I4VEC_SORT_HEAP_INDEX_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int b;
  int base = 0;
  int c;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "TEST78\n";
  cout << "  I4VEC_SORT_HEAP_INDEX_A creates an ascending\n";
  cout << "  sort index for an I4VEC.\n";
  cout << "  I4VEC_SORT_HEAP_INDEX_D creates a descending\n";
  cout << "  sort index for an I4VEC.\n";

  b = 0;
  c = 3 * N;
  seed = 123456789;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  indx = i4vec_sort_heap_index_a ( N, base, a );

  cout << "\n";
  cout << "  After indexed ascending sort:\n";
  cout << "\n";
  cout << "         I   INDX(I)    A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout                          << "  "
         << setw(8) <<        i   << "  "
         << setw(8) <<   indx[i]  << "  "
         << setw(8) <<      a[i]  << "\n";
  }

  cout << "\n";
  cout << "  The index array to carries out the permutation implicitly.\n";
  cout << "\n";
  cout << "       I   INDX(I)  A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout                          << "  "
         << setw(8) <<        i   << "  "
         << setw(8) <<   indx[i]  << "  "
         << setw(8) << a[indx[i]] << "\n";
  }

  cout << "\n";
  cout << "  I4VEC_PERMUTE carries out the permutation explicitly.\n";
  cout << "\n";

  i4vec_permute ( N, indx, base, a );

  i4vec_print ( N, a, "  I, A(I)" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void test79 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST79 tests I4VEC_SORT_HEAP_INDEX_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int base;
  int b;
  int c;
  int i;
  int *indx;
  int seed;

  cout << "\n";
  cout << "TEST79\n";
  cout << "  I4VEC_SORT_HEAP_INDEX_D creates a descending\n";
  cout << "  sort index for an I4VEC.\n";

  seed = 123456789;

  b = 0;
  c = 3 * N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  indx = i4vec_sort_heap_index_d ( N, a );

  cout << "\n";
  cout << "  After indexed descending sort:\n";
  cout << "\n";
  cout << "         I   INDX(I)      A(I)\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << indx[i]
         << "  " << setw(8) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  The index array carries out the permutation implicitly.\n";
  cout << "\n";
  cout << "         I   INDX(I)  A(INDX(I))\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << indx[i]
         << "  " << setw(8) << a[indx[i]] << "\n";
  }

  cout << "\n";
  cout << "  I4VEC_PERMUTE carries out the permutation explicitly.\n";
  cout << "\n";

  base = 0;
  i4vec_permute ( N, indx, base, a );

  i4vec_print ( N, a, "  I, A(I)" );

  delete [] a;
  delete [] indx;

  return;
# undef N
}
//****************************************************************************80

void test80 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST80 tests I4VEC_SORT_INSERT_A;
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
# define N 30

  int *a;
  int b;
  int c;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST80\n";
  cout << "  I4VEC_SORT_INSERT_A sorts an integer array;\n";

  b = 0;
  c = N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  i4vec_sort_insert_a ( N, a );

  i4vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test81 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST81 tests I4VEC_SORT_QUICK_A;
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
# define N 30

  int *a;
  int b;
  int c;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST81\n";
  cout << "  I4VEC_SORT_QUICK_A sorts an integer array;\n";

  b = 0;
  c = N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  i4vec_sort_quick_a ( N, a );

  i4vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test82 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST82 tests I4VEC_SORT_SHELL_A;
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
# define N 30

  int *a;
  int b;
  int c;
  int i;
  int seed;

  cout << "\n";
  cout << "TEST82\n";
  cout << "  I4VEC_SORT_SHELL_A sorts an integer array;\n";

  b = 0;
  c = N;
  seed = 123456789;

  cout << "\n";
  cout << "  Using random seed " << seed << ".\n";

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  i4vec_sort_shell_a ( N, a );

  i4vec_print ( N, a, "  Sorted array:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test83 ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4LIB_TEST83 tests I4VEC_SORTED_UNDEX.
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
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  int x_val[X_NUM] = { 11, 11, 11, 22, 22, 33, 33, 55, 55 };
  int *xdnu;
  int *xu_val;

  cout << "\n";
  cout << "I4LIB_TEST83\n";
  cout << "  I4VEC_SORTED_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique elements of a sorted I4VEC,\n";
  cout << "  and a map from the original vector to the (implicit)\n";
  cout << "  vector of sorted unique elements.\n";

  i4vec_print ( x_num, x_val, "  The vector X:" );

  x_unique_num = i4vec_sorted_unique_count ( x_num, x_val );

  undx = new int[x_unique_num];
  xu_val = new int[x_unique_num];

  xdnu = new int[x_num];

  cout << "\n";
  cout << "  Number of unique entries in X is " << x_unique_num << "\n";

  i4vec_sorted_undex ( x_num, x_val, x_unique_num, undx, xdnu );

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
  cout << "     I  UNDX XU(I)\n";
  cout << "\n";
  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(4) << xu_val[i] << "\n";
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

void test84 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST84 tests I4VEC_SORTED_UNIQUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int *a;
  int b = 0;
  int c = N;
  int seed;
  int unique_num;

  cout << "\n";
  cout << "TEST84\n";
  cout << "  I4VEC_SORTED_UNIQUE finds unique entries in a sorted array.\n";

  seed = 123456789;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_sort_heap_a ( N, a );

  i4vec_print ( N, a, "  Input vector:" );

  unique_num = i4vec_sorted_unique ( N, a );

  i4vec_print ( unique_num, a, "  Unique entries:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test85 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST85 tests I4VEC_TRANSPOSE_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  int *a;

  cout << "\n";
  cout << "TEST85\n";
  cout << "  I4VEC_TRANSPOSE_PRINT prints an integer vector\n";
  cout << "  with 5 entries to a row, and an optional title.\n";

  a = i4vec_indicator_new ( N );

  i4vec_print ( N, a, "  Output from I4VEC_PRINT:" );

  cout << "\n";
  cout << "  Now call I4VEC_TRANSPOSE_PRINT with a short title:\n";
  cout << "\n";

  i4vec_transpose_print ( N, a, "  My array:  " );

  cout << "\n";
  cout << "  Now call I4VEC_TRANSPOSE_PRINT with no title:\n";
  cout << "\n";

  i4vec_transpose_print ( N, a, " " );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test86 ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4LIB_TEST86 tests I4VEC_UNDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define X_NUM 9

  int i;
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  int x_val[X_NUM] = { 33, 55, 11, 11, 55, 33, 22, 22, 11 };
  int *xdnu;
  int *xu_val;

  cout << "\n";
  cout << "I4LIB_TEST86\n";
  cout << "  I4VEC_UNDEX produces index vectors which create a sorted\n";
  cout << "  list of the unique elements of an (unsorted) I4VEC,\n";
  cout << "  and a map from the original vector to the (implicit)\n";
  cout << "  vector of sorted unique elements.\n";

  i4vec_print ( x_num, x_val, "  The vector X:" );

  x_unique_num = i4vec_unique_count ( x_num, x_val );

  undx = new int[x_unique_num];
  xu_val = new int[x_unique_num];

  xdnu = new int[x_num];

  cout << "\n";
  cout << "  Number of unique entries in X is " << x_unique_num << "\n";

  i4vec_undex ( x_num, x_val, x_unique_num, undx, xdnu );

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
  cout << "     I  UNDX XU(I)\n";
  cout << "\n";
  for ( i = 0; i < x_unique_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(4) << undx[i]
         << "  " << setw(4) << xu_val[i] << "\n";
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

void test87 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST87 tests I4VEC_UNIQUE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *a;
  int b;
  int c;
  int i;
  int n = 20;
  int seed;
  int *unique_index;

  seed = 123456789;

  cout << "\n";
  cout << "TEST87\n";
  cout << "  I4VEC_UNIQUE_INDEX, for each entry in an I4VEC\n";
  cout << "  indexes the unique elements.\n";

  b = 1;
  c = 5;

  a = i4vec_uniform_new ( n, b, c, seed );

  unique_index = i4vec_unique_index ( n, a );

  cout << "\n";
  cout << "         I      A(I)    UNIQUE\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << a[i]
         << "  " << setw(8) << unique_index[i] << "\n";
  }

  delete [] a;
  delete [] unique_index;

  return;
}
//****************************************************************************80

void test88 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST88 tests I4VEC_VALUE_INDEX.
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
  int *a;
  int b;
  int c;
  int max_index = 3;
  int n = 25;
  int n_index;
  int seed = 123456789;
  int value = 3;
  int *value_index;

  cout << "\n";
  cout << "TEST88\n";
  cout << "  I4VEC_VALUE_INDEX indexes entries equal to\n";
  cout << "  a given value.\n";
  cout << "\n";
  cout << "  The desired value is " << value << "\n";
  cout << "  Maximum number of indices to find is " << max_index << "\n";

  b = 1;
  c = 5;

  a = i4vec_uniform_new ( n, b, c, seed );

  i4vec_print ( n, a, "  Input vector A:" );

  value_index = i4vec_value_index ( n, a, value, max_index, &n_index );

  i4vec_print ( n_index, value_index,
    "  Indices of entries equal to given value: " );

  delete [] a;
  delete [] value_index;

  return;
}
//****************************************************************************80

void test89 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST89 tests I4VEC_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int *a;
  int b;
  int c;
  int seed;
  double variance;

  cout << "\n";
  cout << "TEST89\n";
  cout << "  For an integer vector:\n";
  cout << "  I4VEC_VARIANCE:      variance.\n";

  seed = 123456789;

  b = - N;
  c = N;

  a = i4vec_uniform_new ( N, b, c, seed );

  i4vec_print ( N, a, "  Input vector:" );

  variance = i4vec_variance ( N, a );

  cout << "\n";
  cout << "  Variance: " << variance << "\n";

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test90 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST90 tests I4VEC2_SORT_A, I4VEC2_SORT_D, and I4VEC2_SORTED_UNIQUE.
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
  int b;
  int c;
  int *ivec;
  int *jvec;
  int n = 10;
  int unique_num;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST90\n";
  cout << "  For a pair of integer vectors:\n";
  cout << "  I4VEC2_SORT_A ascending sorts;\n";
  cout << "  I4VEC2_SORT_D descending sorts;\n";
  cout << "  I4VEC2_SORTED_UNIQUE counts unique entries.\n";

  b = 1;
  c = 3;

  ivec = i4vec_uniform_new ( n, b, c, seed );

  jvec = i4vec_uniform_new ( n, b, c, seed );

  ivec[2] = ivec[0];
  jvec[2] = jvec[0];

  ivec[4] = ivec[1];
  jvec[4] = jvec[1];

  ivec[8] = ivec[0];
  jvec[8] = jvec[0];

  i4vec2_print ( n, ivec, jvec, "  The array:" );

  i4vec2_sort_a ( n, ivec, jvec );

  i4vec2_print ( n, ivec, jvec, "  After ascending sort:" );

  i4vec2_sort_d ( n, ivec, jvec );

  i4vec2_print ( n, ivec, jvec, "  After descending sort:" );

  i4vec2_sorted_unique ( n, ivec, jvec, &unique_num );

  i4vec2_print ( unique_num, ivec, jvec, "  After UNIQ:" );

  delete [] ivec;
  delete [] jvec;

  return;
}
