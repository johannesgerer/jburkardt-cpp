# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa121.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA121_PRB.
//
//  Discussion:
//
//    ASA121_PRB tests the ASA121 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA121_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA121 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA121_PRB:\n";
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
//    TEST01 demonstrates the use of TRIGAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2008
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
  cout << "  TRIGAMMA computes the trigamma function. \n";
  cout << "  We compare the result to tabulated values.\n";
  cout << "\n";
  cout << "          X                     "
       << "FX                        FX2\n";
  cout << "                                "
       << "(Tabulated)               (TRIGAMMA)                DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    trigamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = trigamma ( x, &ifault );

    cout << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
