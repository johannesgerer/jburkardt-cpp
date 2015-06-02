# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "cycle_floyd.hpp"

int main ( );

void test01 ( );
int f1 ( int i );
void test02 ( );
int f2 ( int i );
void test03 ( );
int f3 ( int i );
void test04 ( );
int f4 ( int i );
void test05 ( );
int f5 ( int i );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CYCLE_FLOYD_PRB.
//
//  Discussion:
//
//    CYCLE_FLOYD_PRB tests the CYCLE_FLOYD library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "CYCLE_FLOYD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CYCLE_FLOYD library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CYCLE_FLOYD_PRB\n";
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
//    TEST01 tests CYCLE_FLOYD for a tiny example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int lam;
  int mu;
  int x0;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test CYCLE_FLOYD on F1().\n";
  cout << "  f1(0) = 6.\n";
  cout << "  f1(1) = 6.\n";
  cout << "  f1(2) = 0.\n";
  cout << "  f1(3) = 1.\n";
  cout << "  f1(4) = 4.\n";
  cout << "  f1(5) = 3.\n";
  cout << "  f1(6) = 3.\n";
  cout << "  f1(7) = 4.\n";
  cout << "  f1(8) = 0.\n";

  x0 = 2;
  cout << "\n";
  cout << "  Starting argument X0 = " << x0 << "\n";

  cycle_floyd ( f1, x0, lam, mu );

  cout << "\n";
  cout << "  Reported cycle length is " << lam << "\n";
  cout << "  Expected value is 3\n";
  cout << "\n";
  cout << "  Reported distance to first cycle element is " << mu << "\n";
  cout << "  Expected value is 2\n";

  return;
}
//****************************************************************************80

int f1 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    F1 is the iteration function for example 1.
//
//  Discussion:
//
//    This function has two cycles:
//
//    6, 3, 1, of length 3
//    4, of length 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the argument of the function.
//
//    Output, int F1, the value of the function.
//
{
  static int f_table[9] = { 6, 6, 0, 1, 4, 3, 3, 4, 0 };
  int value;

  value = f_table[i];

  return value;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CYCLE_FLOYD for F2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int lam;
  int mu;
  int x0;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test CYCLE_FLOYD for F2().\n";
  cout << "  f2(i) = mod ( 22 * i + 1, 72 ).\n";

  x0 = 0;
  cout << "\n";
  cout << "  Starting argument X0 = " << x0 << "\n";

  cycle_floyd ( f2, x0, lam, mu );

  cout << "\n";
  cout << "  Reported cycle length is " << lam << "\n";
  cout << "  Expected value is 9\n";
  cout << "\n";
  cout << "  Reported distance to first cycle element is " << mu << "\n";
  cout << "  Expected value is 3\n";

  return;
}
//****************************************************************************80

int f2 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    F2 is the iteration function for example 2.
//
//  Discussion:
//
//    This function has a cycle
//
//    3, 67, 35, 51, 43, 11, 27, 19, 59, of length 9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the argument of the function.
//
//    Output, int F2, the value of the function.
//
{
  int value;

  value = ( 22 * i + 1 ) % 72;

  return value;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests CYCLE_FLOYD for F3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int lam;
  int mu;
  int x0;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Test CYCLE_FLOYD for F3().\n";
  cout << "  f3(i) = mod ( 123 * i + 456, 100000 ).\n";

  x0 = 789;
  cout << "\n";
  cout << "  Starting argument X0 = " << x0 << "\n";

  cycle_floyd ( f3, x0, lam, mu );

  cout << "\n";
  cout << "  Reported cycle length is " << lam << "\n";
  cout << "  Expected value is 50000\n";
  cout << "\n";
  cout << "  Reported distance to first cycle element is " << mu << "\n";
  cout << "  Expected value is 0\n";

  return;
}
//****************************************************************************80

int f3 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    F3 is the iteration function for example 3.
//
//  Discussion:
//
//    This function has a cycle of length 50000
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the argument of the function.
//
//    Output, int F3, the value of the function.
//
{
  int value;

  value = ( 123 * i + 456 ) % 1000000;

  return value;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests CYCLE_FLOYD for F4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int lam;
  int mu;
  int x0;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Test CYCLE_FLOYD for F4().\n";
  cout << "  f4(i) = mod ( 31421 * i + 6927, 65536 ).\n";

  x0 = 1;
  cout << "\n";
  cout << "  Starting argument X0 = " << x0 << "\n";

  cycle_floyd ( f4, x0, lam, mu );

  cout << "\n";
  cout << "  Reported cycle length is " << lam << "\n";
  cout << "  Expected value is 65536\n";
  cout << "\n";
  cout << "  Reported distance to first cycle element is " << mu << "\n";
  cout << "  Expected value is 0\n";

  return;
}
//****************************************************************************80

int f4 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    F4 is the iteration function for example 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the argument of the function.
//
//    Output, int F4, the value of the function.
//
{
  int value;

  value = ( 31421 * i + 6927 ) % 65536;

  return value;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests CYCLE_FLOYD for F5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int lam;
  int mu;
  int x0;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Test CYCLE_FLOYD for F5().\n";
  cout << "  f5(i) = mod ( 16383 * i + 1, 65536 ).\n";

  x0 = 1;
  cout << "\n";
  cout << "  Starting argument X0 = " << x0 << "\n";

  cycle_floyd ( f5, x0, lam, mu );

  cout << "\n";
  cout << "  Reported cycle length is " << lam << "\n";
  cout << "  Expected value is 8\n";
  cout << "\n";
  cout << "  Reported distance to first cycle element is " << mu << "\n";
  cout << "  Expected value is 0\n";

  i = 0;
  x0 = 1;
  cout << "  " << i << "  " << x0 << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    x0 = f5 ( x0 );
    cout << "  " << i << "  " << x0 << "\n";
  }

  return;
}
//****************************************************************************80

int f5 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    F5 is the iteration function for example 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the argument of the function.
//
//    Output, int F5, the value of the function.
//
{
  int value;

  value = ( 16383 * i + 1 ) % 65536;

  return value;
}
