# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa005.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA005_PRB.
//
//  Discussion:
//
//    ASA005_PRB tests the ASA005 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA005_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA005 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA005_PRB:\n";
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
//    TEST01 demonstrates the use of PRNCST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double fx;
  double fx2;
  int ifault;
  double lambda;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PRNCST computes the noncentral Student's\n";
  cout << "  Cumulative Density Function.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "      X   LAMBDA  DF     "
       << " CDF                       CDF                     DIFF\n";
  cout << "                         "
       << " Tabulated                 PRNCST\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    student_noncentral_cdf_values ( &n_data, &df, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = prncst ( x, df, lambda, &ifault );

    cout << "  " << setprecision(2) << setw(6) << x
         << "  " << setprecision(2) << setw(6) << lambda
         << "  "                    << setw(2)  << df
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
