# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fd1d_heat_steady.hpp"

int main ( );
double k1 ( double x );
double f1 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int n = 11;
  double *u;
  string u_file;
  double ua;
  double ub;
  double *x;
  string x_file;

  timestamp ( );
  cout << "\n";
  cout << "PROBLEM1:\n";
  cout << "  C++ version\n";
  cout << "  A test problem for FD1D_HEAT_STEADY.\n";

  a = 0.0;
  b = 1.0;

  x = r8vec_even ( n, a, b );

  ua = 0.0;
  ub = 1.0;

  u = fd1d_heat_steady ( n, a, b, ua, ub, k1, f1, x );

  x_file = "problem1_nodes.txt";
  r8mat_write ( x_file, 1, n, x );

  cout << "\n";
  cout << "  X data written to \"" << x_file << "\".\n";

  u_file = "problem1_values.txt";
  r8mat_write ( u_file, 1, n, u );

  cout << "  U data written to \"" << u_file << "\".\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "PROBLEM1:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  delete [] u;
  delete [] x;

  return 0;
}
//****************************************************************************80

double k1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    K1 evaluates the heat transfer coefficient K(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the position.
//
//    Output, double K1, the value of K(X).
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double f1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F1 evaluates the right hand side of the steady state heat equation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the position.
//
//    Output, double F1, the value of F(X).
//
{
  double value;

  value = 0.0;

  return value;
}


