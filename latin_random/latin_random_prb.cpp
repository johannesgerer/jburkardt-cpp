# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "latin_random.hpp"

int main ( );
void test01 ( int &seed );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LATIN_RANDOM_PRB.
//
//  Discussion:
//
//    LATIN_RANDOM_PRB tests the LATIN_RANDOM library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int seed;
  int test;

  timestamp ( );
  cout << "\n";
  cout << "LATIN_RANDOM_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LATIN_RANDOM library.\n";

  seed = 123456789;

  for ( test = 0; test < 3; test++ )
  {
    test01 ( seed );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "LATIN_RANDOM_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests LATIN_RANDOM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int m = 2;
  int i;
  int j;
  int k;
  int kk;
  int n = 10;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  LATIN_RANDOM chooses a Latin Square cell arrangement,\n";
  cout << "  and then chooses a random point from each cell.\n";
  cout << "\n";
  cout << "  Spatial dimension = " << m << "\n";
  cout << "  Number of points =  " << n << "\n";
  cout << "  Initial seed for UNIFORM = " << seed << "\n";

  x = latin_random_new ( m, n, seed );

  r8mat_transpose_print ( m, n, x, "  Latin Random Square:" );

  delete [] x;

  return;
}
