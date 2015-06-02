# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa243.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA243_PRB.
//
//  Discussion:
//
//    ASA243_PRB tests the ASA243 library.
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
  timestamp ( );
  cout << "\n";
  cout << "ASA243_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA243 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA243_PRB:\n";
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
//    TEST01 demonstrates the use of TNC.
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
  double delta;
  int df;
  double df_real;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  TNC computes the noncentral Student T\n";
  cout << "  Cumulative Density Function.\n";
  cout << "  Compare with tabulated values.\n";
  cout << "\n";
  cout << "        X         LAMBDA        DF     "
       << " CDF             CDF           DIFF\n";
  cout << "                                       "
       << " Tabulated       PRNCST\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    student_noncentral_cdf_values ( &n_data, &df, &delta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    df_real = ( double ) ( df );

    fx2 = tnc ( x, df_real, delta, &ifault );

    cout << "  " << setw(10) << setprecision(4) << x
         << "  " << setw(10) << setprecision(4) << delta
         << "  " << setw(8)                     << df
         << "  " << setw(24) << setprecision(16) << fx
         << "  " << setw(24) << setprecision(16) << fx2
         << "  " << setw(10) << setprecision(4) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
