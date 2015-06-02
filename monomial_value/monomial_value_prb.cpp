# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "monomial_value.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MONOMIAL_VALUE_PRB.
//
//  Discussion:
//
//    MONOMIAL_VALUE_PRB tests the MONOMIAL_VALUE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "MONOMIAL_VALUE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the MONOMIAL_VALUE library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "MONOMIAL_VALUE_PRB\n";
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
//    TEST01 tests MONOMIAL_VALUE on sets of data in various dimensions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *e;
  int e_max;
  int e_min;
  int i;
  int j;
  int m;
  int n;
  int seed;
  double *v;
  double *x;
  double x_max;
  double x_min;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Using monomial_value to evaluate monomials in\n";
  cout << "  dimensions 1 through 3.\n"; 

  e_min = -3;
  e_max = 6;
  n = 5;
  seed = 123456789;
  x_min = -2.0;
  x_max = +10.0;

  for ( m = 1; m <= 3; m++ )
  {
    cout << "\n";
    cout << "  Spatial dimension M = " << m << "\n";

    e = i4vec_uniform_ab_new ( m, e_min, e_max, seed );
    i4vec_transpose_print ( m, e, "  Exponents:" );
    x = r8mat_uniform_ab_new ( m, n, x_min, x_max, seed );
//
//  To make checking easier, make the X values integers.
//
    r8mat_nint ( m, n, x );
    v = monomial_value ( m, n, e, x );

    cout << "\n";
    cout << "   V(X)         ";
    for ( i = 0; i < m; i++ )
    {
      cout << "      X(" << i << ")";
    }
    cout << "\n";
    cout << "\n";
    for ( j = 0; j < n; j++ )
    {
      cout << setw(14) << v[j] << "  ";
      for ( i = 0; i < m; i++ )
      {
        cout << setw(10) << x[i+j*m];
      }
      cout << "\n";
    }

    delete [] e;
    delete [] v;
    delete [] x;
  }

  return;
}
