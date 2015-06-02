# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "cube_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CUBE_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    CUBE_MONTE_CARLO_PRB tests the CUBE_MONTE_CARLO library.
//    
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CUBE_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CUBE_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CUBE_MONTE_CARLO_PRB\n";
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
//    TEST01 estimates integrals over the unit cube in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 2014
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
  cout << "  Use CUBE01_SAMPLE to estimate integrals\n";
  cout << "  over the interior of the unit cube in 3D.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N";
  cout << "        1";
  cout << "               X";
  cout << "               Y ";
  cout << "              Z";
  cout << "               X^2";
  cout << "              XY";
  cout << "             XZ";
  cout << "              Y^2";
  cout << "             YZ";
  cout << "               Z^2\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = cube01_sample ( n, seed );
    cout << "  " << setw(8) << n;
    for ( j = 0; j < 10; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }

      value = monomial_value ( m, n, e, x );

      result = cube01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
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
    exact = cube01_monomial_integral ( e );
    cout << "  " << setw(14) << exact;
  }
  cout << "\n";

  return;
}
