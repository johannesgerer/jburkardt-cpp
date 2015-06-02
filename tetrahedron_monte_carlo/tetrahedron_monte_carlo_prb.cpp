# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "tetrahedron_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TETRAHEDRON_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    TETRAHEDRON_MONTE_CARLO_PRB tests the TETRAHEDRON_MONTE_CARLO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TETRAHEDRON_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TETRAHEDRON_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TETRAHEDRON_MONTE_CARLO_PRB\n";
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
//    TEST01 uses TETRAHEDRON_SAMPLE_01 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  int e_test[3*10] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    0, 0, 1, 
    2, 0, 0, 
    1, 1, 0, 
    1, 0, 1, 
    0, 2, 0, 
    0, 1, 1, 
    0, 0, 2 };
  double error;
  double exact;
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
  cout << "  Use TETRAHEDRON01_SAMPLE for a Monte Carlo estimate of an\n";
  cout << "  integral over the interior of the unit tetrahedron in 3D.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N        1               X               Y ";
  cout << "              Z               X^2              XY             XZ";
  cout << "              Y^2             YZ               Z^2\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = tetrahedron01_sample ( n, seed );
    cout << "  " << setw(8) << n;

    for ( j = 0; j < 10; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = tetrahedron01_volume ( ) * r8vec_sum ( n, value ) 
        / ( double ) ( n );
      cout << "  " << setw(14) << result;

      delete [] value;
    }

    cout << "\n";

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
    result = tetrahedron01_monomial_integral ( e );
    cout << "  " << setw(14) << result;
  }

  cout << "\n";

  return;
}

