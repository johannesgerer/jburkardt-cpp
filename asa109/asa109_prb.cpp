# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa109.hpp"

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
//    MAIN is the main program for ASA109_PRB.
//
//  Discussion:
//
//    ASA109_PRB tests the ASA109 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2014
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
  test02 ( );
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
//    28 April 2013
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
    beta_inc_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = lgamma ( a )
             + lgamma ( b )
             - lgamma ( a + b );

    x2 = xinbta ( a, b, beta_log, fx, ifault );

    cout << "  " << setprecision(4) << setw(10) << a
         << "  " << setprecision(4) << setw(10) << b
         << "  " << setprecision(4) << setw(10) << fx
         << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << x2
         << "  " << setprecision(4) << setw(10) << fabs ( x - x2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BETA_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2014
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
  cout << "TEST02:\n";
  cout << "  BETA_INC_VALUES stores values of\n";
  cout << "  the incomplete Beta function.\n";
  cout << "\n";
  cout << "      A            B            X            BETA_INC(A,B)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = lgamma ( a )
             + lgamma ( b )
             - lgamma ( a + b );

    ifault = 0;
    fx2 = betain ( x, a, b, beta_log, ifault );

    cout << "  " << setprecision(4) << setw(10) << a
         << "  " << setprecision(4) << setw(10) << b
         << "  " << setprecision(4) << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }
  return;
}
