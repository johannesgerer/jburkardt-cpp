# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa109.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA109_PRB.
//
//  Discussion:
//
//    ASA109_PRB calls the ASA109 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA109_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA109 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA109_PRB:\n";
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
//    TEST01 demonstrates the use of XINBTA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
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
  int ifault;
  int n_data;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  XINBTA inverts the incomplete Beta function.\n";
  cout << "  Given CDF, it computes an X.\n";
  cout << "\n";
  cout << "           A           B           CDF    "
       << "    X                         X\n";
  cout << "                                          "
       << "    (Tabulated)               (XINBTA)            DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = alngam ( a, &ifault )
             + alngam ( b, &ifault )
             - alngam ( a + b, &ifault );

    x2 = xinbta ( a, b, beta_log, fx, &ifault );

    cout << "  " << setprecision(4) << setw(10) << a
         << "  " << setprecision(4) << setw(10) << b
         << "  " << setprecision(4) << setw(10) << fx
         << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << x2
         << "  " << setprecision(4) << setw(10) << r8_abs ( x - x2 ) << "\n";
  }

  return;
}
