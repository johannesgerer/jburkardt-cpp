# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "triangle01_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE01_INTEGRALS_PRB.
//
//  Discussion:
//
//    TRIANGLE01_INTEGRALS_PRB tests the TRIANGLE01_INTEGRALS library.
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
  cout << "TRIANGLE01_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE01_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE01_INTEGRALS_PRB\n";
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
//    TEST01 uses TRIANGLE01_SAMPLE to compare exact and estimated monomial integrals.
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
  double error;
  double exact;
  int i;
  int j;
  int m = 2;
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
  cout << "  over the interior of the unit triangle in 2D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = triangle01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose X, Y exponents.
//
  cout << "\n";
  cout << "  We will restrict this test to randomly chosen even exponents.\n";
  cout << "\n";
  cout << "  Ex  Ey     MC-Estimate      Exact           Error\n";
  cout << "\n";

  for ( i = 0; i <= 4; i++ )
  {
    e[0] = i;
    for ( j = 0; j <= 4; j++ )
    {
      e[1] = j;

      value = monomial_value ( m, n, e, x );

      result = triangle01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      exact = triangle01_monomial_integral ( e );
      error = fabs ( result - exact );

      cout << "  " << setw(2) << e[0]
           << "  " << setw(2) << e[1]
           << "  " << setw(14) << result
           << "  " << setw(14) << exact
           << "  " << setw(10) << error << "\n";

      delete [] value;
    }
  }

  delete [] x;

  return;
}
