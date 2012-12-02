# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa239.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA239_PRB.
//
//  Discussion:
//
//    ASA239_PRB calls the ASA239 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA239_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA239 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA239_PRB:\n";
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
//    TEST01 demonstrates the use of GAMMAD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  GAMMAD computes the incomplete Gamma function.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "             A             X      "
       << "    FX                        FX2\n";
  cout << "                                  "
       << "    (Tabulated)               (GAMMAD)            DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gammad ( x, a, &ifault );

    cout << "  " << setprecision(4) << setw(12) << a
         << "  " << setprecision(4) << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
