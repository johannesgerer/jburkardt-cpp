# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "disk_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DISK_INTEGRALS_PRB.
//
//  Discussion:
//
//    DISK_INTEGRALS_PRB tests the DISK_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "DISK_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the DISK_INTEGRALS library.\n";

  test01 ( );
/*
  Terminate.
*/
  cout << "\n";
  cout << "DISK_INTEGRALS_PRB\n";
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
//    TEST01 uses DISK01_SAMPLE to estimate a number of integrals.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2014
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
  cout << "  over the interior of the unit disk in 2D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = disk01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose X,Y exponents between 0 and 8.
//
  cout << "\n";
  cout << "  If any exponent is odd, the integral is zero.\n";
  cout << "  We will restrict this test to randomly chosen even exponents.\n";
  cout << "\n";
  cout << "  Ex  Ey     MC-Estimate           Exact      Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 4, seed );
    for ( i = 0; i < m; i++ )
    {
      e[i] = e[i] * 2;
    }
    value = monomial_value ( m, n, e, x );

    result = disk01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = disk01_monomial_integral ( e );
    error = fabs ( result - exact );

    cout << "  " << setw(2) << e[0]
         << "  " << setw(2) << e[1]
         << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(10) << error << "\n";

    delete [] e;
    delete [] value;
  }

  delete [] x;

  return;
}
