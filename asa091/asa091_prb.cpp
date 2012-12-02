# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa091.hpp"

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
//    MAIN is the main program for ASA091_PRB.
//
//  Discussion:
//
//    ASA091_PRB calls the ASA091 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA091_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA091 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA091_PRB:\n";
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
//    TEST01 makes a single simple calculation with PPCHI2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double g;
  int ifault;
  double p;
  double v;
  double value;
  double value_correct = 0.4;

  p = 0.017523;
  v = 4.0;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Perform a simple sample calculation using\n";
  cout << "  PPCHI2 to invert the Chi-Squared CDF.\n";

  g = alngam ( v / 2.0, &ifault );

  cout << "\n";
  cout << "  P =                  "
       << setw(24) << setprecision(16) << p << "\n";
  cout << "  V =                  "
       << setw(24) << setprecision(16) << v << "\n";
  cout << "  G Log(Gamma(V/2)) =  "
       << setw(24) << setprecision(16) << g << "\n";

  value = ppchi2 ( p, v, g, &ifault );

  cout << "  VALUE =              "
       << setw(24) << setprecision(16) << value << "\n";
  cout << "  VALUE (correct) =    "
       << setw(24) << setprecision(16) << value_correct << "\n";

  cout << "\n";
  cout << "  Error flag IFAULT = " << ifault << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 compares PPCHI2 against tabulated values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double fx;
  double g;
  int ifault;
  int n_data;
  double v;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  PPCHI2 computes percentage points of the Chi-Square CDF.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "         N        CDF       X                        "
       << " X2                    DIFF\n";
  cout << "                           (tabulated)               "
       << "(PPCHI2)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    v = ( double ) ( a );

    g = alngam ( v / 2.0, &ifault );

    x2 = ppchi2 ( fx, v, g, &ifault );

    cout << "  " << setprecision(4) << setw(10) << a
         << "  " << setprecision(4) << setw(10) << fx
         << "  " << setprecision(16) << setw(24) << x
         << "  " << setprecision(16) << setw(24) << x2
         << "  " << setprecision(4) << setw(10) << r8_abs ( x - x2 ) << "\n";
  }

  return;
}
