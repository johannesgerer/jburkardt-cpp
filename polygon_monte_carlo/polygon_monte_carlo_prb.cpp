# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "polygon_monte_carlo.hpp"

int main ( );
void test01 ( int nv, double v[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POLYGON_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    POLYGON_MONTE_CARLO_PRB tests the POLYGON_MONTE_CARLO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  int nv1 = 4;
  double v1[2*4] = {
    -1.0, -1.0, 
     1.0, -1.0, 
     1.0,  1.0, 
    -1.0,  1.0 };

  timestamp ( );
  cout << "\n";
  cout << "POLYGON_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the POLYGON_MONTE_CARLO library.\n";

  test01 ( nv1, v1 );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLYGON_MONTE_CARLO_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int nv, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 estimates integrals over a polygon in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e[2];
  int e_test[2*7] = {
    0, 0, 
    2, 0, 
    0, 2, 
    4, 0, 
    2, 2, 
    0, 4, 
    6, 0 };
  double error;
  double exact;
  int j;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use POLYGON_SAMPLE to estimate integrals\n";
  cout << "  over the interior of a polygon in 2D.\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N";
  cout << "        1";
  cout << "              X^2 ";
  cout << "             Y^2";
  cout << "             X^4";
  cout << "           X^2Y^2";
  cout << "             Y^4";
  cout << "           X^6\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = polygon_sample ( nv, v, n, seed );

    cout << "  " << setw(8) << n;

    for ( j = 0; j < 7; j++ )
    {
      e[0] = e_test[0+j*2];
      e[1] = e_test[1+j*2];

      value = monomial_value ( 2, n, e, x );

      result = polygon_area ( nv, v ) * r8vec_sum ( n, value ) / ( double ) ( n );
      cout << "  " << setw(14) << result;
    }

    cout << "\n";

    delete [] value;
    delete [] x;

    n = 2 * n;

  }

  cout << "\n";
  cout << "     Exact";
  for ( j = 0; j < 7; j++ )
  {
    e[0] = e_test[0+j*2];
    e[1] = e_test[1+j*2];

    result = polygon_monomial_integral ( nv, v, e );
    cout << "  " << setw(14) << result;
  }

  cout << "\n";

  return;
}
