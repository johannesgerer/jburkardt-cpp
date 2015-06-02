# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "toms291.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS291_PRB.
//
//  Discussion:
//
//    TOMS291_PRB tests the TOMS291 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TOMS291_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the TOMS291 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TOMS291_PRB:\n";
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2008
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
