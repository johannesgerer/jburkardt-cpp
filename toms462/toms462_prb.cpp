# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "toms462.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    TOMS462_PRB tests BIVNOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TOMS462_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TOMS462 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TOMS462_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************807
//
//  Purpose:
//
//    TEST01 tests BIVNOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  double expect;
  double r;
  double x;
  double y;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compare BIVNOR with some simple data\n";
  cout << "  with 3 digit accuracy.\n";
  cout << "\n";
  cout << "       X         Y          R          P               P\n";
  cout << "                                      (Tabulated)     (BIVNOR)\n";
  cout << "\n";

  x =  0.8;
  y = -1.5;
  r =  -0.9;
  expect = 0.148;
  cdf = bivnor ( x, y, r );
  cout << "  " << setw(9) << x
       << "  " << setw(8) << y
       << "  " << setw(8) << r
       << "  " << setw(14) << expect
       << "  " << setw(14) << cdf << "\n";

  x =  0.6;
  y = -1.4;
  r =  -0.7;
  expect = 0.208;
  cdf = bivnor ( x, y, r );
  cout << "  " << setw(9) << x
       << "  " << setw(8) << y
       << "  " << setw(8) << r
       << "  " << setw(14) << expect
       << "  " << setw(14) << cdf << "\n";

  x =  0.2;
  y = -1.0;
  r =  -0.5;
  expect = 0.304;
  cdf = bivnor ( x, y, r );
  cout << "  " << setw(9) << x
       << "  " << setw(8) << y
       << "  " << setw(8) << r
       << "  " << setw(14) << expect
       << "  " << setw(14) << cdf << "\n";

  x = -1.2;
  y =  0.1;
  r =   0.0;
  expect = 0.407;
  cdf = bivnor ( x, y, r );
  cout << "  " << setw(9) << x
       << "  " << setw(8) << y
       << "  " << setw(8) << r
       << "  " << setw(14) << expect
       << "  " << setw(14) << cdf << "\n";

  x = -1.2;
  y = -0.1;
  r =   0.3;
  expect = 0.501;
  cdf = bivnor ( x, y, r );
  cout << "  " << setw(9) << x
       << "  " << setw(8) << y
       << "  " << setw(8) << r
       << "  " << setw(14) << expect
       << "  " << setw(14) << cdf << "\n";

  x = -0.4;
  y = -0.9;
  r =   0.6;
  expect = 0.601;
  cdf = bivnor ( x, y, r );
  cout << "  " << setw(9) << x
       << "  " << setw(8) << y
       << "  " << setw(8) << r
       << "  " << setw(14) << expect
       << "  " << setw(14) << cdf << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BIVNOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fxy1;
  double fxy2;
  int n_data;
  double r;
  double x;
  double y;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Compare BIVNOR with some tabulated data.\n";
  cout << "\n";
  cout << "      X          Y          ";
  cout << "R           P                         P";
  cout << "                      DIFF\n";
  cout << "                                ";
  cout << "       (Tabulated)               (BIVNOR)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( n_data, x, y, r, fxy1 );

    if ( n_data == 0 )
    {
      break;
    }
//
//  BIVNOR computes the "tail" of the probability, and we want the
//  initial part//
//
    fxy2 = bivnor ( - x, - y, r );

    cout << "  " << setw(8) << x
         << "  " << setw(8) << y
         << "  " << setw(8) << r
         << "  " << setw(24) << fxy1
         << "  " << setw(24) << fxy2
         << "  " << setw(10) << r8_abs ( fxy1 - fxy2 ) << "\n";
  }
  return;
}
