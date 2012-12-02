# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa111.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    ASA111_PRB tests ASA111.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "ASA111_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA111 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA111_PRB:\n";
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
//    TEST01 compares PPND against tabulated values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int ifault;
  int n_data;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PPND computes percentage points of the normal distribution.\n";
  cout << "  Compare against tabulated values.\n";
  cout << "\n";
  cout << "         CDF        X                           X  "
       << "                  DIFF\n";
  cout << "                 (tabulated)                   (PPND)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = ppnd ( fx, &ifault );

    cout << "  " << setprecision(4)  << setw(10) << fx
         << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << x2
         << "  " << setprecision(4)  << setw(10) << r8_abs ( x - x2 ) << "\n";
  }

  return;
}
