# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "wedge_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WEDGE_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    WEDGE_MONTE_CARLO_PRB tests the WEDGE_MONTE_CARLO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "WEDGE_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the WEDGE_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WEDGE_MONTE_CARLO_PRB\n";
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
//    TEST01 uses WEDGE01_SAMPLE with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  int e_test[3*8] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    0, 0, 1, 
    2, 0, 0, 
    1, 1, 0, 
    0, 0, 2, 
    3, 0, 0 };
  int i;
  int j;
  int m = 3;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use WEDGE01_SAMPLE for a Monte Carlo estimate of an\n";
  cout << "  integral over the interior of the unit wedge in 3D.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N        1               X               Y ";
  cout << "              Z                X^2            XY              Z^2    ";
  cout << "        X^3\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = wedge01_sample ( n, seed );

    cout << "  " << setw(8) << n;

    for ( j = 0; j < 8; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
 
      value = monomial_value ( m, n, e, x );

      result = wedge01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      cout << "  " << setw(14 ) << result;
      delete [] value;
    }

    cout << "\n";

    delete [] x;

    n = 2 * n;
  }

  cout << "     Exact";

  for ( j = 0; j < 8; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = wedge01_integral ( e );
    cout << "  " << setw(14 ) << result;
  }
  cout << "\n";

  return;
}
