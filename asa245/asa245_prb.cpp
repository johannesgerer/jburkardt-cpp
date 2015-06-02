# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa245.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA245_PRB.
//
//  Discussion:
//
//    ASA245_PRB tests the ASA245 library.
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
  cout << "ASA245_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA245 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA245_PRB:\n";
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
//    TEST01 demonstrates the use of ALNGAM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2008
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
  cout << "  ALNGAM computes the logarithm of the \n";
  cout << "  Gamma function.  We compare the result\n";
  cout << "  to tabulated values.\n";
  cout << "\n";
  cout << "          X                     "
       << "FX                        FX2\n";
  cout << "                                "
       << "(Tabulated)               (ALNGAM)                DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = alngam ( x, ifault );

    cout << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates the use of LNGAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2008
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
  double x;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  LNGAMMA computes the logarithm of the \n";
  cout << "  Gamma function.  We compare the result\n";
  cout << "  to tabulated values.\n";
  cout << "\n";
  cout << "          X                     "
      << "FX                        FX2\n";
  cout << "                                "
       << "(Tabulated)               (LNGAMMA)                DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lngamma ( x, ier );

    cout << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 demonstrates the use of LGAMMA.
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
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  LGAMMA computes the logarithm of the \n";
  cout << "  Gamma function.\n";
  cout << "  LGAMMA is available with the G++ compiler.\n";
  cout << "  Compare the result to tabulated values.\n";
  cout << "\n";
  cout << "          X                     "
      << "FX                        FX2\n";
  cout << "                                "
       << "(Tabulated)               (LNGAMMA)                DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lgamma ( x );

    cout << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
