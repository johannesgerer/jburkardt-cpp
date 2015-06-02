# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pyramid_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PYRAMID_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    PYRAMID_MONTE_CARLO_PRB tests the PYRAMID_MONTE_CARLO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PYRAMID_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the PYRAMID_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PYRAMID_MONTE_CARLO_PRB\n";
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
//    TEST01 compares exact and estimated monomial integrals.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 3
# define TEST_NUM 10

  int e[M];
  int e_test[M*TEST_NUM] = {
    0, 0, 0, 
    0, 0, 1, 
    2, 0, 0, 
    0, 2, 0, 
    0, 0, 2, 
    2, 0, 1, 
    0, 2, 1, 
    0, 0, 3, 
    2, 2, 0, 
    2, 0, 2 };
  double error;
  double exact;
  int i;
  int j;
  int m = M;
  int n;
  double result;
  int seed;
  int test_num = TEST_NUM;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use PYRAMID01_SAMPLE to estimate integrals\n";
  cout << "  over the interior of the unit pyramid in 3D.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N";
  cout << "        1";
  cout << "               Z";
  cout << "             X^2";
  cout << "             Y^2";
  cout << "             Z^2";
  cout << "            X^2Z";
  cout << "            Y^2Z";
  cout << "             Z^3";
  cout << "          X^2Y^2";
  cout << "          X^2Z^2\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    cout << "  " << setw(8) << n;

    x = pyramid01_sample ( n, seed );

    for ( j = 0; j < test_num; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = pyramid01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      cout << "  " << setw(14) << result;
    }
    cout << "\n";

    delete [] value;
    delete [] x;

    n = 2 * n;
  }

  cout << "\n";
  cout << "     Exact";

  for ( j = 0; j < 10; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = pyramid01_integral ( e );
    cout << "  " << setw(14) << result;
  }
  cout << "\n";

  return;
}
