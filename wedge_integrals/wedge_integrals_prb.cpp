# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "wedge_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WEDGE_INTEGRALS_PRB.
//
//  Discussion:
//
//    WEDGE_INTEGRALS_PRB tests the WEDGE_INTEGRALS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "WEDGE_INTEGRALS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the WEDGE_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WEDGE_INTEGRALS_PRB\n";
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
//    21 August 2014
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
  double error;
  double exact;
  int expon[3];
  int m = 3;
  int n = 500000;
  double q;
  int seed;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compare exact and estimated integrals \n";
  cout << "  over the unit wedge in 3D.\n";
//
//  Get sample points.
//
  seed = 123456789;
  x = wedge01_sample ( n, seed );

  cout << "\n";
  cout << "  Number of sample points used is " << n << "\n";
  cout << "\n";
  cout << "   E1  E2  E3     MC-Estimate      Exact           Error\n";
  cout << "\n";
//
//  Check all monomials up to total degree E_MAX.
//
  for ( e3 = 0; e3 <= e_max; e3++ )
  {
    expon[2] = e3;
    for ( e2 = 0; e2 <= e_max - e3; e2++ )
    {
      expon[1] = e2;
      for ( e1 = 0; e1 <= e_max - e3 - e2; e1++ )
      {
        expon[0] = e1;

        value = monomial_value ( m, n, expon, x );

        q = wedge01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
        exact = wedge01_integral ( expon );
        error = fabs ( q - exact );

        cout << setw(4) << expon[0] << "  "
             << setw(2) << expon[1] << "  "
             << setw(2) << expon[2] << "  "
             << setw(14) << q << "  "
             << setw(14) << exact << "  "
             << setw(14) << error << "\n";

        delete [] value;
      }
    }
  }

  delete [] x;

  return;
}
