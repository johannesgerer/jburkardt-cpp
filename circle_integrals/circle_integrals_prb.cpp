# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "circle_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CIRCLE_INTEGRALS_PRB.
//
//  Discussion:
//
//    CIRCLE_INTEGRALS_PRB tests the CIRCLE_INTEGRALS library.
//    
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CIRCLE_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CIRCLE_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CIRCLE_INTEGRALS_PRB\n";
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
//    TEST01 uses CIRCLE01_SAMPLE with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 January 2014
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
  cout << "  Use CIRCLE01_SAMPLE to compare exact and\n";
  cout << "  estimated integrals along the circumference\n";
  cout << "  of the unit circle in 2D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = circle01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose X, Y exponents.
//
  cout << "\n";
  cout << "  If any exponent is odd, the integral is zero.\n";
  cout << "  We restrict this test to randomly chosen even exponents.\n";
  cout << "\n";
  cout << "  Ex  Ey     MC-Estimate           Exact      Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 5, seed );

    for ( i = 0; i < m; i++ )
    {
      e[i] = e[i] * 2;
    }

    value = monomial_value ( m, n, e, x );

    result = circle01_length ( ) * r8vec_sum ( n, value ) 
      / ( double ) ( n );
    exact = circle01_monomial_integral ( e );
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
