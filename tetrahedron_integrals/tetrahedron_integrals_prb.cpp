# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "tetrahedron_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TETRAHEDRON_INTEGRALS_PRB.
//
//  Discussion:
//
//    TETRAHEDRON_INTEGRALS_PRB tests the TETRAHEDRON_INTEGRALS library.
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
  cout << "TETRAHEDRON_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TETRAHEDRON_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TETRAHEDRON_INTEGRALS_PRB\n";
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
//    TEST01 uses TETRAHEDRON_SAMPLE_01 to compare exact and estimated integrals.
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
  double error;
  double exact;
  int i;
  int j;
  int k;
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
  cout << "  Estimate monomial integrals using Monte Carlo\n";
  cout << "  over the interior of the unit tetrahedron in 3D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = tetrahedron01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
// Run through the exponents.
//
  cout << "\n";
  cout << "  Ex  Ey  Ez     MC-Estimate      Exact           Error\n";
  cout << "\n";

  for ( i = 0; i <= 3; i++ )
  {
    e[0] = i;
    for ( j = 0; j <= 3; j++ )
    {
      e[1] = j;
      for ( k = 0; k <= 3; k++ )
      {
        e[2] = k;

        value = monomial_value ( m, n, e, x );

        result = tetrahedron01_volume ( ) * r8vec_sum ( n, value ) 
          / ( double ) ( n );
        exact = tetrahedron01_monomial_integral ( e );
        error = fabs ( result - exact );

        cout << "  " << setw(2) << e[0]
             << "  " << setw(2) << e[1]
             << "  " << setw(2) << e[2]
             << "  " << setw(14) << result
             << "  " << setw(14) << exact
             << "  " << setw(10) << error << "\n";

        delete [] value;
      }
    }
  }

  delete [] x;

  return;
}

