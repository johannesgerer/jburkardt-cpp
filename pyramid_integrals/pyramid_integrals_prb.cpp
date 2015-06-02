# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pyramid_integrals.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PYRAMID_INTEGRALS_PRB.
//
//  Discussion:
//
//    PYRAMID_INTEGRALS_PRB tests the PYRAMID_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PYRAMID_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the PYRAMID_INTEGRALS library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PYRAMID_INTEGRALS_PRB\n";
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
//    14 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e_max = 6;
  int e1;
  int e2;
  int e3;
  int expon[3];
  double error;
  double exact;
  int m = 3;
  static int n = 500000;
  double q;
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compare exact and estimated integrals \n";
  cout << "  over the unit pyramid in 3D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = pyramid01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
  cout << "\n";
  cout << "   E1  E2  E3     MC-Estimate      Exact           Error\n";
  cout << "\n";
//
//  Check all monomials, with only even dependence on X or Y, 
//  up to total degree E_MAX.
//
  for ( e3 = 0; e3 <= e_max; e3++ )
  {
    expon[2] = e3;
    for ( e2 = 0; e2 <= e_max - e3; e2 = e2 + 2 )
    {
      expon[1] = e2;
      for ( e1 = 0; e1 <= e_max - e3 - e2; e1 = e1 + 2 )
      {
        expon[0] = e1;

        value = monomial_value ( m, n, expon, x );

        q = pyramid01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
        exact = pyramid01_integral ( expon );
        error = fabs ( q - exact );

        cout << "  " << setw(2) << expon[0]
             << "  " << setw(2) << expon[1]
             << "  " << setw(2) << expon[2]
             << "  " << setw(14) << q
             << "  " << setw(14) << exact
             << "  " << setw(12) << error << "\n";

        delete [] value;
      }
    }
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
//    TEST02 examines the sample points in the pyramid
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  int m = 3;
  int n = 20;
  int seed;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Print sample points in the unit pyramid in 3D.\n";
  seed = 123456789;
  x = pyramid01_sample ( n, seed );
  r8mat_transpose_print ( 3, n, x, "  Unit pyramid points" );

  delete [] x;

  return;
}
