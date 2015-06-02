# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "blas0.hpp"

int main ( );
void test01 ( );
void test015 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS0_PRB.
//
//  Discussion:
//
//    BLAS0_PRB tests the BLAS0 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BLAS0_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BLAS0 library.\n";

  test01 ( );
  test015 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BLAS0_PRB\n";
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
//    TEST01 tests R4_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2014
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
  cout << "TEST01\n";
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

void test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests R4_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2014
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5

  int test;
  float x;
  float x_test[TEST_NUM] = { -1.25, -0.25, 0.0, +0.5, +9.0 };

  cout << "\n";
  cout << "TEST015\n";
  cout << "  R4_SIGN returns the sign of a number.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setw(8) << x
         << "  " << setw(8) << r4_sign ( x ) << "\n";
  }

  return;
# undef TEST_N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8_ABS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2014
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
  cout << "TEST02\n";
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

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R8_SIGN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2014
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
  cout << "TEST03\n";
  cout << "  R8_SIGN returns the sign of a number.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    cout << "  " << setw(8) << x
         << "  " << setw(8) << r8_sign ( x ) << "\n";
  }

  return;
# undef TEST_NUM
}
