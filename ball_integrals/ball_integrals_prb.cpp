# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "ball_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BALL_INTEGRALS_PRB.
//
//  Discussion:
//
//    BALL_INTEGRALS_PRB tests the BALL_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BALL_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BALL_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BALL_INTEGRALS_PRB\n";
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
//    TEST01 uses BALL01_SAMPLE to compare estimated and exact integrals.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2014
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
  cout << "  Estimate monomial integrals using Monte Carlo\n";
  cout << "  over the interior of the unit ball in 3D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = ball01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
//
//  Randomly choose X,Y exponents between 0 and 8.
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

    result = ball01_volume ( ) * r8vec_sum ( n, value ) 
      / ( double ) ( n );
    exact = ball01_monomial_integral ( e );
    error = fabs ( result - exact );

    cout << "  " << setw(2) << e[0]
         << "  " << setw(2) << e[1]
         << "  " << setw(2) << e[2]
         << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(10) << error << "\n";

    delete [] e;
    delete [] value;
  }
  delete [] x;

  return;
}
