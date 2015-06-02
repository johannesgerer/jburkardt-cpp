# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "combination_lock.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for COMBINATION_LOCK_PRB.
//
//  Discussion:
//
//    COMBINATION_LOCK_PRB tests the COMBINATION_LOCK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "COMBINATION_LOCK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the COMBINATION_LOCK libary.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "COMBINATION_LOCK_PRB\n";
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
//    TEST01 tests BICYCLE_LOCK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int m = 3;
  int n = 10;
  int seed;
  int step;
  
  cout << "\n";
  cout << "TEST01\n";
  cout << "  A bicycle combination lock consists of 3 dials,\n";
  cout << "  each having 10 symbols, 0 through 9.\n";
  cout << "  We seek to determine the combination C.\n";
//
//  Report on the problem data.
//
  cout << "\n";
  cout << "  The number of dials is M = " << m << "\n";
  cout << "  The number of symbols is N = " << n << "\n";
  cout << "  The number of possible combinations is M^N = " << i4_power ( n, m ) << "\n";

  seed = get_seed ( );
  c = i4_uniform ( 0, 999, seed );

  cout << "  The \"secret\" combination is " << c << "\n";

  step = bicycle_lock ( c );

  if ( step == -1 )
  {
    cout << "\n";
    cout << "  The combination was not found!\n";
  }
  else
  {
    cout << "\n";
    cout << "  The combination was found on step " << step << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests COMBINATION_LOCK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of dials.
//
//    Input, int N, the number of symbols on each dial.
//    We assume the symbols are the integers 0 to N-1.
//
//    Input, int C[M], the combination.
//
{
  int c[4] = { 1, 2, 3, 4 };
  int m = 4;
  int n = 5;
  int step;
  
  cout << "\n";
  cout << "TEST02\n";
  cout << "  A combination lock consists of M dials,\n";
  cout << "  each having N symbols.\n";
  cout << "  We seek to determine the combination C.\n";
//
//  Report on the problem data.
//
  cout << "\n";
  cout << "  The number of dials is M = " << m << "\n";
  cout << "  The number of symbols is N = " << n << "\n";
  cout << "  The number of possible combinations is M^N = " << i4_power ( n, m ) << "\n";

  i4vec_print ( m, c, "  The \"secret\" combination:" );

  step = combination_lock ( m, n, c );

  if ( step == -1 )
  {
    cout << "\n";
    cout << "  The combination was not found!\n";
  }
  else
  {
    cout << "\n";
    cout << "  The combination was found on step " << step << "\n";
  }

  return;
}
