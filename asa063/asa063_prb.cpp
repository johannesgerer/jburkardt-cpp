# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa063.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA063_PRB.
//
//  Discussion:
//
//    ASA063_PRB calls the ASA063 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA063_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA063 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA063_PRB:\n";
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
//    TEST01 demonstrates the use of BETAIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double beta_log;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  BETAIN computes the incomplete Beta function.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "           A           B           X      "
       << "    FX                        FX2\n";
  cout << "                                          "
       << "    (Tabulated)               (BETAIN)            DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = alogam ( a, &ifault )
             + alogam ( b, &ifault )
             - alogam ( a + b, &ifault );

    fx2 = betain ( x, a, b, beta_log, &ifault );

    cout << "  " << setprecision(4) << setw(10) << a
         << "  " << setprecision(4) << setw(10) << b
         << "  " << setprecision(4) << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
