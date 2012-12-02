# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem1d_bvp_linear.hpp"

int main ( );
void fem1d_bvp_linear_test01 ( );
void fem1d_bvp_linear_test02 ( );
void fem1d_bvp_linear_test03 ( );
void fem1d_bvp_linear_test04 ( );
void fem1d_bvp_linear_test05 ( );
void fem1d_bvp_linear_test06 ( );
double a1 ( double x );
double a2 ( double x );
double a3 ( double x );
double c1 ( double x );
double c2 ( double x );
double c3 ( double x );
double exact1 ( double x );
double exact_ux1 ( double x );
double exact2 ( double x );
double exact3 ( double x );
double exact4 ( double x );
double exact_ux4 ( double x );
double f1 ( double x );
double f2 ( double x );
double f3 ( double x );
double f4 ( double x );
double f5 ( double x );
double f6 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR_PRB tests the routines in FEM1D_BVP_LINEAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
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
  cout << "  Test the FEM1D_BVP_LINEAR library.\n";

  fem1d_bvp_linear_test01 ( );
  fem1d_bvp_linear_test02 ( );
  fem1d_bvp_linear_test03 ( );
  fem1d_bvp_linear_test04 ( );
  fem1d_bvp_linear_test05 ( );
  fem1d_bvp_linear_test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void fem1d_bvp_linear_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR_TEST01 carries out test case #1.
//
//  Discussion:
//
//    Use A1, C1, F1, EXACT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 11;
  double l2_error;
  double seminorm_error;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_TEST01\n";
  cout << "  A1(X)  = 1.0\n";
  cout << "  C1(X)  = 0.0\n";
  cout << "  F1(X)  = X * ( X + 3 ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
//
//  Geometry definitions.
//
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a1, c1, f1, x );

  cout << "\n";
  cout << "     I    X         U         Uexact    Error\n";
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

  l2_error = compute_l2_error ( n, x, u, exact1 );
  seminorm_error = compute_seminorm_error ( n, x, u, exact_ux1 );
  cout << "\n";
  cout << "  L2 norm of error  = " << l2_error << "\n";
  cout << "  Seminorm of error = " << seminorm_error << "\n";

  delete [] u;
  delete [] x;

  return;
}
//****************************************************************************80

void fem1d_bvp_linear_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR_TEST02 carries out test case #2.
//
//  Discussion:
//
//    Use A1, C2, F2, EXACT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 11;
  double l2_error;
  double seminorm_error;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_TEST02\n";
  cout << "  A1(X)  = 1.0\n";
  cout << "  C2(X)  = 2.0\n";
  cout << "  F2(X)  = X * ( 5 - X ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
//
//  Geometry definitions.
//
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a1, c2, f2, x );

  cout << "\n";
  cout << "     I    X         U         Uexact    Error\n";
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

  l2_error = compute_l2_error ( n, x, u, exact1 );
  seminorm_error = compute_seminorm_error ( n, x, u, exact_ux1 );
  cout << "\n";
  cout << "  L2 norm of error  = " << l2_error << "\n";
  cout << "  Seminorm of error = " << seminorm_error << "\n";

  delete [] u;
  delete [] x;

  return;
}
//****************************************************************************80

void fem1d_bvp_linear_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR_TEST03 carries out test case #3.
//
//  Discussion:
//
//    Use A1, C3, F3, EXACT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 11;
  double l2_error;
  double seminorm_error;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_TEST03\n";
  cout << "  A1(X)  = 1.0\n";
  cout << "  C3(X)  = 2.0 * X\n";
  cout << "  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
//
//  Geometry definitions.
//
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a1, c3, f3, x );

  cout << "\n";
  cout << "     I    X         U         Uexact    Error\n";
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

  l2_error = compute_l2_error ( n, x, u, exact1 );
  seminorm_error = compute_seminorm_error ( n, x, u, exact_ux1 );
  cout << "\n";
  cout << "  L2 norm of error  = " << l2_error << "\n";
  cout << "  Seminorm of error = " << seminorm_error << "\n";

  delete [] u;
  delete [] x;

  return;
}
//****************************************************************************80

void fem1d_bvp_linear_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR_TEST04 carries out test case #4.
//
//  Discussion:
//
//    Use A2, C1, F4, EXACT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 11;
  double l2_error;
  double seminorm_error;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_TEST04\n";
  cout << "  A2(X)  = 1.0 + X * X\n";
  cout << "  C1(X)  = 0.0\n";
  cout << "  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
//
//  Geometry definitions.
//
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a2, c1, f4, x );

  cout << "\n";
  cout << "     I    X         U         Uexact    Error\n";
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

  l2_error = compute_l2_error ( n, x, u, exact1 );
  seminorm_error = compute_seminorm_error ( n, x, u, exact_ux1 );
  cout << "\n";
  cout << "  L2 norm of error  = " << l2_error << "\n";
  cout << "  Seminorm of error = " << seminorm_error << "\n";

  delete [] u;
  delete [] x;

  return;
}
//****************************************************************************80

void fem1d_bvp_linear_test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR_TEST05 carries out test case #5.
//
//  Discussion:
//
//    Use A3, C1, F5, EXACT1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 11;
  double l2_error;
  double seminorm_error;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_TEST05\n";
  cout << "  A3(X)  = 1.0 + X * X for X <= 1/3\n";
  cout << "         = 7/9 + X     for      1/3 < X\n";
  cout << "  C1(X)  = 0.0\n";
  cout << "  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n";
  cout << "                       for X <= 1/3\n";
  cout << "         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) .* exp ( X )\n";
  cout << "                       for      1/3 <= X\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
//
//  Geometry definitions.
//
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_linear ( n, a3, c1, f5, x );

  cout << "\n";
  cout << "     I    X         U         Uexact    Error\n";
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

  l2_error = compute_l2_error ( n, x, u, exact1 );
  seminorm_error = compute_seminorm_error ( n, x, u, exact_ux1 );
  cout << "\n";
  cout << "  L2 norm of error  = " << l2_error << "\n";
  cout << "  Seminorm of error = " << seminorm_error << "\n";

  delete [] u;
  delete [] x;

  return;
}
/******************************************************************************/

void fem1d_bvp_linear_test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_BVP_LINEAR_TEST06 does an error analysis.
//
//  Discussion:
//
//    Use A1, C1, F6, EXACT4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double l2_error;
  double seminorm_error;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  cout << "\n";
  cout << "FEM1D_BVP_LINEAR_TEST06\n";
  cout << "  A1(X)  = 1.0n";
  cout << "  C1(X)  = 0.0\n";
  cout << "  F6(X)  = pi*pi*sin(pi*X)\n";
  cout << "  U4(X)  = sin(pi*x)\n";
  cout << "\n";
  cout << "  Compute L2 norm and seminorm of error for various N.\n";
  cout << "\n";
  cout << "     N        L2 error      Seminorm error\n";
  cout << "\n";

  n = 11;
  for ( i = 0; i <= 4; i++ )
  {
//
//  Geometry definitions.
//
    x_first = 0.0;
    x_last = 1.0;
    x = r8vec_even_new ( n, x_first, x_last );

    u = fem1d_bvp_linear ( n, a1, c1, f6, x );

    l2_error = compute_l2_error ( n, x, u, exact4 );
    seminorm_error = compute_seminorm_error ( n, x, u, exact_ux4 );

    cout << "  " << setw(4) << n
         << "  " << setw(14) << l2_error
         << "  " << setw(14) << seminorm_error << "\n";

    n = 2 * ( n - 1 ) + 1;

    delete [] u;
    delete [] x;
  }

  return;
}
//****************************************************************************80

double a1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    A1 evaluates A function #1.
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
//    Output, double A1, the value of A(X).
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double a2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    A2 evaluates A function #2.
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
//    Output, double A2, the value of A(X).
//
{
  double value;

  value = 1.0 + x * x;

  return value;
}
//****************************************************************************80

double a3 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    A3 evaluates A function #3.
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
//    Output, double A3, the value of A(X).
//
{
  double value;

  if ( x <= 1.0 / 3.0 )
  {
    value = 1.0 + x * x;
  }
  else
  {
    value = x + 7.0 / 9.0;
  }

  return value;
}
//****************************************************************************80

double c1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    C1 evaluates C function #1.
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
//    Output, double C1, the value of C(X).
//
{
  double value;

  value = 0.0;

  return value;
}
//****************************************************************************80

double c2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    C2 evaluates C function #2.
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
//    Output, double C2, the value of C(X).
//
{
  double value;

  value = 2.0;

  return value;
}
//****************************************************************************80

double c3 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    C3 evaluates C function #3.
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
//    Output, double C3, the value of C(X).
//
{
  double value;

  value = 2.0 * x;

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
//****************************************************************************80

double exact_ux1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT_UX1 evaluates the derivative of exact solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double EXACT_UX1, the value of dUdX(X).
//
{
  double value;

  value = ( 1.0 - x - x * x ) * exp ( x );

  return value;
}
//****************************************************************************80

double exact2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT2 returns exact solution #2.
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
//    Output, double EXACT2, the value of U(X).
//
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value =  x * ( 1.0 - x ) * exp ( x );
  }
  else
  {
    value = x * ( 1.0 - x )  * exp ( 2.0 / 3.0 );
  }
  return value;
}
//****************************************************************************80

double exact3 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT3 returns exact solution #3.
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
//    Output, double EXACT3, the value of U(X).
//
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = x * ( 1.0 - x ) * exp ( x );
  }
  else
  {
    value = x * ( 1.0 - x );
  }
  return value;
}
//****************************************************************************80

double exact4 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT4 returns exact solution #4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double EXACT4, the value of U(X).
//
{
  static double pi = 3.141592653589793;
  double value;

  value = sin ( pi * x );

  return value;
}
//****************************************************************************80

double exact_ux4 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT_UX4 returns the derivative of exact solution #4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double EXACT_UX4, the value of U(X).
//
{
  static double pi = 3.141592653589793;
  double value;

  value = pi * cos ( pi * x );

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

double f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F2 evaluates right hand side function #2.
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
//    Output, double F2, the value of F(X).
//
{
  double value;

  value = x * ( 5.0 - x ) * exp ( x );

  return value;
}
//****************************************************************************80

double f3 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F3 evaluates right hand side function #3.
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
//    Output, double F3, the value of F(X).
//
{
  double value;

  value = - x * ( 2.0 * x * x - 3.0 * x - 3.0 ) * exp ( x );

  return value;
}
//****************************************************************************80

double f4 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F4 evaluates right hand side function #4.
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
//    Output, double F4, the value of F(X).
//
{
  double value;

  value = ( x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x ) * exp ( x );

  return value;
}
//****************************************************************************80

double f5 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F5 evaluates right hand side function #5.
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
//    Output, double F5, the value of F(X).
//
{
  double value;

  if ( x <= 1.0 / 3.0 )
  {
    value = ( x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x ) * exp ( x );
  }
  else
  {
    value = ( - 1.0 + ( 10.0 / 3.0 ) * x 
      + ( 43.0 / 9.0 ) * x * x + x * x * x ) * exp ( x );
  }

  return value;
}
//****************************************************************************80

double f6 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F6 evaluates right hand side function #6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double F6, the value of F(X).
//
{
  static double pi = 3.141592653589793;
  double value;

  value = pi * pi * sin ( pi * x );

  return value;
}
