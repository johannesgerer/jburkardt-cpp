# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "simplex_integrals.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SIMPLEX_INTEGRALS_PRB.
//
//  Discussion:
//
//    SIMPLEX_INTEGRALS_PRB tests the SIMPLEX_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SIMPLEX_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SIMPLEX_INTEGRALS library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SIMPLEX_INTEGRALS_PRB\n";
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
//    TEST01 compares exact and estimated integrals in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *e;
  double error;
  double exact;
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
  cout << "  over the interior of the unit simplex in M dimensions.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = simplex01_sample ( m, n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose exponents.
//
  cout << "\n";
  cout << "  We randomly choose the exponents.\n";
  cout << "\n";
  cout << "  Ex  Ey  Ez     MC-Estimate      Exact           Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 4, seed );

    value = monomial_value ( m, n, e, x );

    result = simplex01_volume ( m ) * r8vec_sum ( n, value ) 
      / ( double ) ( n );

    exact = simplex01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    cout << "  " << setw(2) << e[0]
         << "  " << setw(2) << e[1]
         << "  " << setw(2) << e[2]
         << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(14) << error << "\n";

    delete [] e;
    delete [] value;
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 compares exact and estimated integrals in 6D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *e;
  double error;
  double exact;
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
  cout << "  Estimate monomial integrals using Monte Carlo\n";
  cout << "  over the interior of the unit simplex in M dimensions.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = simplex01_sample ( m, n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose exponents.
//
  cout << "\n";
  cout << "  We randomly choose the exponents.\n";
  cout << "\n";
  cout << "  E1  E2  E3  E4  E5  E6     MC-Estimate      Exact           Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 4, seed );

    value = monomial_value ( m, n, e, x );

    result = simplex01_volume ( m ) * r8vec_sum ( n, value ) 
      / ( double ) ( n );

    exact = simplex01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    cout << "  " << setw(2) << e[0]
         << "  " << setw(2) << e[1]
         << "  " << setw(2) << e[2]
         << "  " << setw(2) << e[3]
         << "  " << setw(2) << e[4]
         << "  " << setw(2) << e[5] 
         << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(14) << error << "\n";

    delete [] e;
    delete [] value;
  }

  return;
}
