# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa103.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA103_PRB.
//
//  Discussion:
//
//    ASA103_PRB calls the ASA103 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA103_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA103 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA103_PRB:\n";
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
//    TEST01 demonstrates the use of DIGAMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
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
  cout << "  DIGAMA computes the Digama or Psi function. \n";
  cout << "  Compare the result to tabulated values.\n";
  cout << "\n";
  cout << "          X       "
       << "FX                        FX2\n";
  cout << "                  "
       << "(Tabulated)               (DIGAMA)                DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = digama ( x, &ifault );

    cout << "  " << setprecision(4) << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
