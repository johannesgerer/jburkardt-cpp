# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "beta_nc.H"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BETA_NC_PRB.
//
//  Discussion:
//
//    BETA_NC_PRB tests the BETA_NC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BETA_NC_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the BETA_NC library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BETA_NC_PRB:\n";
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
//    TEST01 demonstrates the use of BETA_NONCENTRAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error_max;
  double fx;
  double fx2;
  double lambda;
  int n_data;
  double x;

  error_max = 1.0E-10;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  BETA_NONCENTRAL_CDF computes the noncentral incomplete \n";
  cout << "  Beta function.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "      A      B     LAMBDA        X      "
       << "    FX                        FX2\n";
  cout << "                                        "
       << "    (Tabulated)               (computed)          DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = beta_noncentral_cdf ( a, b, lambda, x, error_max );

    cout << "  " << setprecision(2) << setw(5) << a
         << "  " << setprecision(2) << setw(5) << b
         << "  " << setprecision(3) << setw(7) << lambda
         << "  " << setprecision(4) << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
