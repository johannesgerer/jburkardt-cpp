# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "circle_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CIRCLE_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    CIRCLE_MONTE_CARLO_PRB tests the CIRCLE_MONTE_CARLO library.
//    
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CIRCLE_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CIRCLE_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CIRCLE_MONTE_CARLO_PRB\n";
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
//    TEST01 uses CIRCLE01_SAMPLE with an increasing number of points.
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
  int e[2];
  int e_test[2*7] = {
    0, 0, 
    2, 0, 
    0, 2, 
    4, 0, 
    2, 2, 
    0, 4,
    6, 0 };
  double error;
  double exact;
  int i;
  int j;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use CIRCLE01_SAMPLE to estimate integrals\n";
  cout << "  along the circumference of the unit circle in 2D.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N        1              X^2             Y^2";
  cout << "             X^4           X^2Y^2          Y^4          X^6\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = circle01_sample ( n, seed );
    cout << "  " << setw(8) << n;
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        e[i] = e_test[i+j*2];
      }

      value = monomial_value ( 2, n, e, x );

      result = circle01_length ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
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
    for ( i = 0; i < 2; i++ )
    {
      e[i] = e_test[i+j*2];
    }
    exact = circle01_monomial_integral ( e );
    cout << "  " << setw(14) << exact;
  }
  cout << "\n";

  return;
}
