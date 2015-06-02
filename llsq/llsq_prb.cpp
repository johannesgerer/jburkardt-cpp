# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "llsq.hpp"

using namespace std;

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LLSQ_PRB.
//
//  Discussion:
//
//    LLSQ_PRB tests the LLSQ library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LLSQ_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LLSQ library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LLSQ_PRB\n";
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
//    TEST01 calls LLSQ to match 15 data values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  int i;
  int n = 15;
  double x[15] = { 
    1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 
    1.73, 1.75, 1.78, 1.80, 1.83 };
  double y[15] = {
    52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47,
    66.28, 68.10, 69.92, 72.19, 74.46 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  LLSQ can compute the formula for a line of the form\n";
  cout << "    y = A * x + B\n";
  cout << "  which minimizes the RMS error to a set of N data values.\n";

  llsq ( n, x, y, a, b );

  cout << "\n";
  cout << "  Estimated relationship is y = " << a << " * x + " << b << "\n";
  cout << "  Expected value is         y = 61.272 * x - 39.062\n";
  cout << "\n";
  cout << "     I      X       Y      B+A*X    |error|\n";
  cout << "\n";
  error = 0.0;
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(7) << x[i]
         << "  " << setw(7) << y[i]
         << "  " << setw(7) << b + a * x[i]
         << "  " << setw(7) << b + a * x[i] - y[i] << "\n";
    error = error + pow ( b + a * x[i] - y[i], 2 );
  }
  error = sqrt ( error / ( double ) n );
  cout << "\n";
  cout << "  RMS error =                      " << error << "\n";

  return;
}
