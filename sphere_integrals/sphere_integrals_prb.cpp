# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "sphere_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_INTEGRALS_PRB.
//
//  Discussion:
//
//    SPHERE_INTEGRALS_PRB tests the SPHERE_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPHERE_INTEGRALS_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_INTEGRALS_PRB\n";
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
//    TEST01 uses SPHERE01_SAMPLE to estimate monomial integrands.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2014
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
  int n;
  double result;
  int seed;
  int test;
  const int test_num = 20;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Estimate monomial integrands using Monte Carlo\n";
  cout << "  over the surface of the unit sphere in 3D.\n";
//
//  Get sample points.
//
  n = 8192;
  seed = 123456789;
  x = sphere01_sample ( n, seed );
  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose X,Y,Z exponents between (0,0,0) and (9,9,9).
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

    result = sphere01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = sphere01_monomial_integral ( e );
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

  delete [] x;

  return;
}
