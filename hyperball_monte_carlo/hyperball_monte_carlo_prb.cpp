# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hyperball_monte_carlo.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HYPERBALL_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    HYPERBALL_MONTE_CARLO_PRB tests the HYPERBALL_MONTE_CARLO library.
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
  cout << "HYPERBALL_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the HYPERBALL_MONTE_CARLO library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HYPERBALL_MONTE_CARLO_PRB\n";
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
//    TEST01 uses HYPERBALL01_SAMPLE to estimate integrals in 3D.
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
  int e[3];
  int e_test[3*7] = {
    0, 0, 0, 
    2, 0, 0, 
    0, 2, 0, 
    0, 0, 2, 
    4, 0, 0, 
    2, 2, 0, 
    0, 0, 4 };
  int i;
  int j;
  int m = 3;
  int n;
  double result[7];
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use Monte Carlo to estimate integrals \n";
  cout << "  over the interior of the unit hyperball in M dimensions\n";
  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N        1              X^2             Y^2 ";
  cout << "             Z^2             X^4           X^2Y^2           Z^4\n";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = hyperball01_sample ( m, n, seed );

    cout << "  " << setw(8) << n;
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result[j] = hyperball01_volume ( m ) * r8vec_sum ( n, value ) 
        / ( double ) ( n );
      cout << "  " << setw(14) << result[j];

      delete [] value;
    }
    cout << "\n";

    delete [] x;

    n = 2 * n;
  }

  cout << "\n";
  cout << "     Exact";
  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result[j] = hyperball01_monomial_integral ( m, e );
    cout << "  " << setw(14) << result[j];
  }
  cout << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses HYPERBALL01_SAMPLE to estimate integrals in 6D.
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
  int e[6];
  int e_test[6*7] = {
    0, 0, 0, 0, 0, 0, 
    1, 0, 0, 0, 0, 0, 
    0, 2, 0, 0, 0, 0, 
    0, 2, 2, 0, 0, 0, 
    0, 0, 0, 4, 0, 0, 
    2, 0, 0, 0, 2, 2, 
    0, 0, 0, 0, 0, 6 };
  int i;
  int j;
  int m = 6;
  int n;
  double result[7];
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use Monte Carlo to estimate integrals \n";
  cout << "  over the interior of the unit hyperball in M dimensions\n";
  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";

  seed = 123456789;

  cout << "\n";
  cout << "         N";
  cout << "        1      ";
  cout << "        U      ";
  cout << "         V^2   ";
  cout << "         V^2W^2";
  cout << "         X^4   ";
  cout << "         Y^2Z^2";
  cout << "         Z^6";
  cout << "\n";

  n = 1;

  while ( n <= 65536 )
  {
    x = hyperball01_sample ( m, n, seed );

    cout << "  " << setw(8) << n;
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result[j] = hyperball01_volume ( m ) * r8vec_sum ( n, value ) 
        / ( double ) ( n );
      cout << "  " << setw(14) << result[j];

      delete [] value;
    }
    cout << "\n";

    delete [] x;

    n = 2 * n;
  }

  cout << "\n";
  cout << "     Exact";
  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result[j] = hyperball01_monomial_integral ( m, e );
    cout << "  " << setw(14) << result[j];
  }
  cout << "\n";

  return;
}
