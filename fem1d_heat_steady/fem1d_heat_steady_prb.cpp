# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem1d_heat_steady.hpp"

int main ( );
void fem1d_heat_steady_test01 ( );
double k1 ( double x );
double f1 ( double x );
double exact1 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_HEAT_STEADY_PRB tests the routines in FEM1D_HEAT_STEADY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FEM1D_HEAT_STEADY library.\n";

  fem1d_heat_steady_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM1D_HEAT_STEADY_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void fem1d_heat_steady_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_HEAT_STEADY_TEST01 carries out test case #1.
//
//  Discussion:
//
//    Use K1, F1, EXACT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int n = 11;
  double *u;
  double ua;
  double ub;
  double uexact;
  double *x;

  cout << "\n";
  cout << "FEM1D_HEAT_STEADY_TEST01\n";
  cout << "  K1(X)  = 1.0\n";
  cout << "  F1(X)  = X * ( X + 3 ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
//
//  Geometry definitions.
//
  a = 0.0;
  b = 1.0;
  ua = 0.0;
  ub = 0.0;
  x = r8vec_even_new ( n, a, b );

  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
  cout << "  Left endpoint A = " << a << "\n";
  cout << "  Right endpoint B = " << b << "\n";
  cout << "  Prescribed U(A) = " << ua << "\n";
  cout << "  Prescribed U(B) = " << ub << "\n";

  u = fem1d_heat_steady ( n, a, b, ua, ub, k1, f1, x );

  cout << "\n";
  cout << "     I         X          U                Uexact      Error\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    cout << "  " << setw(4) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(14) << u[i]
         << "  " << setw(14) << uexact
         << "  " << setw(14) << r8_abs ( u[i] - uexact ) << "\n";
  }

  delete [] u;
  delete [] x;

  return;
}
//****************************************************************************80

double k1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    K1 evaluates K function #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
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
//    F1 evaluates right hand side function #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double F1, the value of F(X).
//
{
  double value;

  value = x * ( x + 3.0 ) * exp ( x );

  return value;
}
//****************************************************************************80

double exact1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT1 evaluates exact solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double EXACT1, the value of U(X).
//
{
  double value;

  value = x * ( 1.0 - x ) * exp ( x );

  return value;
}
