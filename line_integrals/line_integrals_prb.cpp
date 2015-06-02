# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "line_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LINE_INTEGRALS_PRB.
//
//  Discussion:
//
//    LINE_INTEGRALS_PRB tests the LINE_INTEGRALS library.
//    
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LINE_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LINE_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LINE_INTEGRALS_PRB\n";
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
//    TEST01 compares exact and estimated monomial integrals.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double error;
  double exact;
  int m = 1;
  int n = 4192;
  double result;
  int seed;
  int test;
  int test_num = 11;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compare exact and estimated integrals \n";
  cout << "  over the length of the unit line in 1D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = line01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
  cout << "\n";
  cout << "   E     MC-Estimate      Exact           Error\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    e = test - 1;

    value = monomial_value_1d ( n, e, x );

    result = line01_length ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = line01_monomial_integral ( e );
    error = fabs ( result - exact );

    cout << "  " << setw(2) << e
         << "  " << setw(14) << result
         << "  " << setw(14) << exact
         << "  " << setw(10) << error << "\n";

    delete [] value;
  }

  delete [] x;

  return;
}
