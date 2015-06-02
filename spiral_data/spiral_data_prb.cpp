# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "spiral_data.hpp"

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
//    MAIN is the main program for SPIRAL_DATA_PRB.
//
//  Discussion:
//
//    SPIRAL_DATA_PRB tests the SPIRAL_DATA library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPIRAL_DATA_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPIRAL_DATA library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPIRAL_DATA_PRB\n";
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
//    TEST01 generates a field and estimates its range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double c;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Sample a spiral velocity field and estimate\n";
  cout << "  the range of the solution values.\n";

  n = 1000;

  xy_lo = +0.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  c = 1.0;

  u = new double[n];
  v = new double[n];
  uv_spiral ( n, x, y, c, u, v );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:  " << "  " << setw(14) << r8vec_min ( n, u )
                   << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:  " << "  " << setw(14) << r8vec_min ( n, v )
                   << "  " << setw(14) << r8vec_max ( n, v ) << "\n";

  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
/*
  Purpose:

    TEST02 generates a field and samples its residuals.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt
*/
{
  int n;
  double c;
  double *pr;
  int seed;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Sample a spiral velocity field and estimate the\n";
  cout << "  range of residuals in the continuity equation.\n";

  n = 1000;

  xy_lo = +0.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  c = 1.0;

  pr = new double[n];
  resid_spiral ( n, x, y, c, pr );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  Pr:  " << "  " << setw(14) << r8vec_amin ( n, pr )
                    << "  " << setw(14) << r8vec_amax ( n, pr ) << "\n";

  delete [] pr;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 generates a field on a regular grid and plots it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  string header;
  int n;
  double s;
  int seed;
  double *u;
  double *v;
  double *x;
  double x_hi;
  double x_lo;
  int x_num = 21;
  double *y;
  double y_hi;
  double y_lo;
  int y_num = 21;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Generate a spiral velocity field on a regular grid.\n";
  cout << "  Store in GNUPLOT data and command files.\n";

  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = new double[x_num*y_num];
  y = new double[x_num*y_num];

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;
  c = 1.0;

  u = new double[x_num*y_num];
  v = new double[x_num*y_num];

  uv_spiral ( n, x, y, c, u, v );

  header = "spiral";
  s = 0.05;
  spiral_gnuplot ( header, n, x, y, u, v, s );

  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}