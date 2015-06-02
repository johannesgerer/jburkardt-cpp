# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "cube_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CUBE_INTEGRALS_PRB.
//
//  Discussion:
//
//    CUBE_INTEGRALS_PRB tests the CUBE_INTEGRALS library.
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
  cout << "CUBE_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CUBE_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CUBE_INTEGRALS_PRB\n";
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
  int *e;
  double error;
  double exact;
  int i;
  int j;
  int m = 3;
  int n = 4192;
  double result;
  int seed;
  int test;
  int test_num = 20;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compare exact and estimated integrals\n";
  cout << "  over the interior of the unit cube in 3D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = cube01_sample ( n, seed );
  cout << "\n";
  cout << "  Number of sample points is " << n << "\n";
//
//  Randomly choose exponents.
//
  cout << "\n";
  cout << "  Ex  Ey  Ez     MC-Estimate           Exact      Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 7, seed );

    value = monomial_value ( m, n, e, x );

    result = cube01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = cube01_monomial_integral ( e );
    error = fabs ( result - exact );

    cout << "  " << setw(2) << e[0]
         << "  " << setw(2) << e[1]
         << "  " << setw(2) << e[2]
         << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(14) << error << "\n";

    delete [] value;
  }

  delete [] x;

  return;
}
