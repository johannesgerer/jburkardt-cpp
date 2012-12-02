# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "toms179.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS179_PRB.
//
//  Discussion:
//
//    TOMS179_PRB calls the TOMS179 routines.
//
//  Modified:
//
//    30 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TOMS179_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the TOMS179 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TOMS179_PRB:\n";
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
//    TEST01 demonstrates the use of ALOGAM.
//
//  Modified:
//
//    30 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  ALOGAM computes the logarithm of the \n";
  cout << "  Gamma function.  We compare the result\n";
  cout << "  to tabulated values.\n";
  cout << "\n";
  cout << "          X                     "
       << "FX                        FX2\n";
  cout << "                                "
       << "(Tabulated)               (ALOGAM)                DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = alogam ( x, &ifault );

    cout << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates the use of MDBETA.
//
//  Modified:
//
//    30 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int ier;
  int n_data;
  double p;
  double q;
  double x;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  MDBETA estimates the value of th modified Beta function.\n";
  cout << "  Compare with tabulated values.\n";
  cout << "\n";
  cout << "         X         P         Q         "
       << "Beta                       Beta                  DIFF\n";
  cout << "                                       "
       << "(Tabulated)                (MDBETA)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_cdf_values ( &n_data, &p, &q, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = mdbeta ( x, p, q, &ier );

    cout << "  " << setprecision(4)  << setw(8)  << x
         << "  " << setprecision(4)  << setw(8)  << p
         << "  " << setprecision(4)  << setw(8)  << q
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4)  << setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
