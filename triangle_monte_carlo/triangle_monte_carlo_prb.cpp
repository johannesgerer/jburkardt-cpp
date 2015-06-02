# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "triangle_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    TRIANGLE_MONTE_CARLO_PRB tests the TRIANGLE_MONTE_CARLO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_MONTE_CARLO_PRB\n";
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
//    TEST01 uses TRIANGLE_SAMPLE_01 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e[2];
  int e_test[2*7] = {
    0, 0, 
    1, 0, 
    0, 1, 
    2, 0, 
    1, 1, 
    0, 2, 
    3, 0 };
  double error;
  double exact;
  int i;
  int j;
  int m = 2;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use TRIANGLE01_SAMPLE for a Monte Carlo estimate of an\n";
  cout << "  integral over the interior of the unit triangle in 2D.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N        1               X               Y ";
  cout << "             X^2               XY             Y^2             X^3\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = triangle01_sample ( n, seed );
    cout << "  " << setw(8) << n;

    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = triangle01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      cout << "  " << setw(14) << result;

      delete [] value;
    }

    cout << "\n";

    delete [] x;

    n = 2 * n;
  }

  cout << "\n";
  cout << "     Exact";
  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = triangle01_monomial_integral ( e );
    cout << "  " << setw(14) << result;
  }

  cout << "\n";

  return;
}
