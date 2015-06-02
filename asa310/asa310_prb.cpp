# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa310.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA310_PRB.
//
//  Discussion:
//
//    ASA310_PRB tests the ASA310 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA310_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA310 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA310_PRB:\n";
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
//    TEST01 demonstrates the use of NCBETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double errmax = 1.0E-10;
  double fx;
  double fx2;
  int ifault;
  double lambda;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  NCBETA computes the noncentral incomplete Beta function.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "      A        B     LAMBDA        X      "
       << "    FX                        FX2\n";
  cout << "                                          "
       << "    (Tabulated)               (NCBETA)            DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = ncbeta ( a, b, lambda, x, errmax, &ifault );

    cout << "  " << setprecision(2) << setw(7) << a
         << "  " << setprecision(2) << setw(7) << b
         << "  " << setprecision(3) << setw(7) << lambda
         << "  " << setprecision(4) << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
