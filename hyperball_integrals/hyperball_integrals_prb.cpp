# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hyperball_integrals.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HYPERBALL_INTEGRALS_PRB.
//
//  Discussion:
//
//    HYPERBALL_INTEGRALS_PRB tests the HYPERBALL_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HYPERBALL_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the HYPERBALL_INTEGRALS library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HYPERBALL_INTEGRALS_PRB\n";
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
//    TEST01 uses HYPERBALL01_SAMPLE to compare exact and estimated integrals in 3D.
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
  int *e;
  double error;
  double exact;
  int i;
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
  cout << "  Use the Monte Carlo method to estimate integrals over\n";
  cout << "  the interior of the unit hyperball in M dimensions.\n";
  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = hyperball01_sample ( m, n, seed );
  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose exponents between 0 and 8.
//
  cout << "\n";
  cout << "  If any exponent is odd, the integral is zero.\n";
  cout << "  We will restrict this test to randomly chosen even exponents.\n";
  cout << "\n";
  cout << "  Ex  Ey  Ez     MC-Estimate           Exact      Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 4, seed );

    for ( i = 0; i < m; i++ )
    {
      e[i] = e[i] * 2;
    }
    value = monomial_value ( m, n, e, x );

    result = hyperball01_volume ( m ) * r8vec_sum ( n, value )
      / ( double ) ( n );
    exact = hyperball01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    for ( i = 0; i < m; i++ )
    {
      cout << "  " << setw(2) << e[i];
    }
    cout << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(10) << error << "\n";

    delete [] e;
    delete [] value;
  }

  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses HYPERBALL01_SAMPLE to compare exact and estimated integrals in 6D.
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
  int *e;
  double error;
  double exact;
  int i;
  int m = 6;
  int n = 4192;
  double result;
  int seed;
  int test;
  int test_num = 20;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use the Monte Carlo method to estimate integrals over\n";
  cout << "  the interior of the unit hyperball in M dimensions.\n";
  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = hyperball01_sample ( m, n, seed );
  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose exponents between 0 and 8.
//
  cout << "\n";
  cout << "  If any exponent is odd, the integral is zero.\n";
  cout << "  We will restrict this test to randomly chosen even exponents.\n";
  cout << "\n";
  cout << "  E1  E2  E3  E4  E5  E6     MC-Estimate           Exact      Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 4, seed );

    for ( i = 0; i < m; i++ )
    {
      e[i] = e[i] * 2;
    }
    value = monomial_value ( m, n, e, x );

    result = hyperball01_volume ( m ) * r8vec_sum ( n, value )
      / ( double ) ( n );
    exact = hyperball01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    for ( i = 0; i < m; i++ )
    {
      cout << "  " << setw(2) << e[i];
    }
    cout << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(10) << error << "\n";

    delete [] e;
    delete [] value;
  }

  delete [] x;

  return;
}
