# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>

using namespace std;

# include "bisection_integer.hpp"

int main ( );
void test01 ( );
int f01 ( int x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BISECTION_INTEGER_PRB.
//
//  Discussion:
//
//    BISECTION_INTEGER_PRB tests the BISECTION_INTEGER library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "BISECTION_INTEGER_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BISECTION_INTEGER library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BISECTION_INTEGER_PRB\n";
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
//    TEST01 tests BISECTION_INTEGER;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int c;
  int fc;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  BISECTION_INTEGER attempts to locate an integer root C\n";
  cout << "  of an equation F(C) = 0.\n";
  cout << "  The user supplies a change of sign interval [A,B].\n";
  cout << "  The function considered here has two real roots\n";
  cout << "  as well as an integer root, so the algorithm can\n";
  cout << "  fail depending on how the change of sign interval is chosen.\n";

  a = 4;
  b = 100;

  cout << "\n";
  cout << "  The initial change of sign interval is:\n";
  cout << "  F(" << a << ") = " << f01 ( a ) << "\n";
  cout << "  F(" << b << ") = " << f01 ( b ) << "\n";

  bisection_integer ( f01, a, b, c, fc );

  if ( fc == 0 )
  {
    cout << "\n";
    cout << "  An exact root was found at C = " << c << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  An exact root was NOT found.\n";
    cout << "  The change of sign interval is now:\n";
    cout << "  F(" << a << ") = " << f01 ( a ) << "\n";
    cout << "  F(" << b << ") = " << f01 ( b ) << "\n";
  }

  a = -10;
  b = 15;

  cout << "\n";
  cout << "  The initial change of sign interval is:\n";
  cout << "  F(" << a << ") = " << f01 ( a ) << "\n";
  cout << "  F(" << b << ") = " << f01 ( b ) << "\n";

  bisection_integer ( f01, a, b, c, fc );

  if ( fc == 0 )
  {
    cout << "\n";
    cout << "  An exact root was found at C = " << c << "\n";
  }
  else
  {
    cout << "\n";
    cout << "  An exact root was NOT found.\n";
    cout << "  The change of sign interval is now:\n";
    cout << "  F(" << a << ") = " << f01 ( a ) << "\n";
    cout << "  F(" << b << ") = " << f01 ( b ) << "\n";
  }

  return;
}
//****************************************************************************80

int f01 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    F01 is a test function.
//
//  Discussion:
//
//    The polynomial has roots 1/2, 7/2, and 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument.
//
//    Output, int F01, the function value.
//
{
  int value;

  value = ( 2 * n - 7 ) * ( 2 * n - 1 ) * ( n - 10 );

  return value;
}
