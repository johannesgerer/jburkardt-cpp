# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "ball_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BALL_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    BALL_MONTE_CARLO_PRB tests the BALL_MONTE_CARLO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BALL_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BALL_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BALL_MONTE_CARLO_PRB\n";
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
//    TEST01 uses BALL01_SAMPLE with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  int e_test[3*7] = {
    0, 0, 0, 
    2, 0, 0, 
    0, 2, 0, 
    0, 0, 2, 
    4, 0, 0, 
    2, 2, 0, 
    0, 0, 4 };
  double error;
  double exact;
  int i;
  int j;
  int n;
  double result[7];
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Estimate integrals over the interior of the unit ball\n";
  cout << "  using the Monte Carlo method.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N        1              X^2             Y^2 ";
  cout << "             Z^2             X^4           X^2Y^2           Z^4\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = ball01_sample ( n, seed );

    cout << "  " << setw(8) << n;
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        e[i] = e_test[i+j*3];
      }
      value = monomial_value ( 3, n, e, x );

      result[j] = ball01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      cout << "  " << setw(14) << result[j];

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
    for ( i = 0; i < 3; i++ )
    {
      e[i] = e_test[i+j*3];
    }
    result[j] = ball01_monomial_integral ( e );
    cout << "  " << setw(14) << result[j];
  }
  cout << "\n";

  return;
}
