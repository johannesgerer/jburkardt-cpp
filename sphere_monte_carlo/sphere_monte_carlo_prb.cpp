# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "sphere_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_MONTE_CARLO_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPHERE_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_MONTE_CARLO_PRB\n";
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
//    TEST01 uses SPHERE_SAMPLE_01 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int e[3];
  int e_test[3*7] = {
    0, 0, 0, 
    2, 0, 0, 
    0, 2, 0, 
    0, 0, 2, 
    4, 0, 0, 
    2, 2, 0, 
    0, 0, 4 };
  double error;
  double exact;
  int i;
  int j;
  int n;
  double pi = 3.1415926535897932384626434;
  double result;
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N        1              X^2             Y^2";
  cout << "             Z^2             X^4            X^2Y^2          Z^4\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = sphere01_sample ( n, &seed );
    cout << "  " << setw(8) << n;
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        e[i] = e_test[i+j*3];
      }

      value = monomial_value ( 3, n, x, e );

      result = 4.0 * pi * r8vec_sum ( n, value ) / ( double ) ( n );
      cout << "  " << setprecision(10) << setw(14) << result;
    }
    cout << "\n";

    delete [] value;
    delete [] x;

    n = 2 * n;
  }

  cout << "\n";
  cout << "  " << "   Exact";
  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      e[i] = e_test[i+j*3];
    }
    exact = sphere01_monomial_integral ( e );
    cout << "  " << setprecision(10) << setw(14) << exact;
  }
  cout << "\n";

  return;
}
