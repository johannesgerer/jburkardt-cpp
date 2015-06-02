# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "fd1d_bvp.hpp"

int main ( );
void fd1d_bvp_test01 ( );
void fd1d_bvp_test02 ( );
void fd1d_bvp_test03 ( );
void fd1d_bvp_test04 ( );
void fd1d_bvp_test05 ( );
double a1 ( double x );
double a1prime ( double x );
double a2 ( double x );
double a2prime ( double x );
double a3 ( double x );
double a3prime ( double x );
double c1 ( double x );
double c2 ( double x );
double c3 ( double x );
double f1 ( double x );
double f2 ( double x );
double f3 ( double x );
double f4 ( double x );
double f5 ( double x );
double *exact1 ( int n, double x[] );
double *exact2 ( int n, double x[] );
double *exact3 ( int n, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FD1D_BVP_PRB.
//
//  Discussion:
//
//    FD1D_BVP_PRB tests the FD1D_BVP library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FD1D_BVP_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FD1D_BVP library.\n";

  fd1d_bvp_test01 ( );
  fd1d_bvp_test02 ( );
  fd1d_bvp_test03 ( );
  fd1d_bvp_test04 ( );
  fd1d_bvp_test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FD1D_BVP_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void fd1d_bvp_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_BVP_TEST01 carries out test case #1.
//
//  Discussion:
//
//    Use A1, C1, F1, EXACT1.
//
//    Repeat using a nonuniform mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int i;
  int n = 21;
  double *u;
  double *u2;
  double *uexact;
  double *x;
  double x1 = 0.0;
  double x2 = 1.0;

  cout << "\n";
  cout << "FD1D_BVP_TEST01\n";
  cout << "  A1(X)  = 1.0\n";
  cout << "  A1'(X) = 0.0\n";
  cout << "  C1(X)  = 0.0\n";
  cout << "  F1(X)  = X * ( X + 3 ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
  cout << "  X1 = " << x1 << "\n";
  cout << "  X2 = " << x2 << "\n";

  x = r8vec_even ( n, x1, x2 );

  u = fd1d_bvp ( n, a1, a1prime, c1, f1, x );

  uexact = exact1 ( n, x );

  cout << "\n";
  cout << "     I         X        U             Uexact         Error\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(12) << u[i]
         << "  " << setw(12) << uexact[i]
         << "  " << setw(14) << r8_abs ( u[i] - uexact[i] ) << "\n";
  }

  delete [] u;
  delete [] uexact;
  delete [] x;

  cout << "\n";
  cout << "  Repeat, using a nonuniform mesh.\n";

  x = r8vec_even ( n, x1, x2 );

  for ( i = 0; i < n; i++ )
  {
    x[i] = sqrt ( x[i] );
  }

  u = fd1d_bvp ( n, a1, a1prime, c1, f1, x );

  uexact = exact1 ( n, x );

  cout << "\n";
  cout << "     I         X        U             Uexact         Error\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(12) << u[i]
         << "  " << setw(12) << uexact[i]
         << "  " << setw(14) << r8_abs ( u[i] - uexact[i] ) << "\n";
  }
//
//  Write the data to files.
//
  filename = "fd1d_bvp_test01_nodes.txt";
  r8mat_write ( filename, n, 1, x );

  u2 = new double[n*2];
  for ( i = 0; i < n; i++ )
  {
    u2[i+0*n] = u[i];
  }
  for ( i = 0; i < n; i++ )
  {
    u2[i+1*n] = uexact[i];
  }

  filename = "fd1d_bvp_test01_values.txt";
  r8mat_write ( filename, n, 2, u2 );

  delete [] u;
  delete [] u2;
  delete [] uexact;
  delete [] x;

  return;
}
//****************************************************************************80

void fd1d_bvp_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_BVP_TEST02 carries out test case #2.
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
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int i;
  int n = 11;
  double *u;
  double *u2;
  double *uexact;
  double *x;
  double x1 = 0.0;
  double x2 = 1.0;

  cout << "\n";
  cout << "FD1D_BVP_TEST02\n";
  cout << "  A1(X)  = 1.0\n";
  cout << "  A1''(X) = 0.0\n";
  cout << "  C2(X)  = 2.0\n";
  cout << "  F2(X)  = X * ( 5 - X ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
  cout << "  X1 = " << x1 << "\n";
  cout << "  X2 = " << x2 << "\n";

  x = r8vec_even ( n, x1, x2 );

  u = fd1d_bvp ( n, a1, a1prime, c2, f2, x );

  uexact = exact1 ( n, x );

  cout << "\n";
  cout << "     I         X        U             Uexact         Error\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(12) << u[i]
         << "  " << setw(12) << uexact[i]
         << "  " << setw(14) << r8_abs ( u[i] - uexact[i] ) << "\n";
  }
//
//  Write the data to files.
//
  filename = "fd1d_bvp_test02_nodes.txt";
  r8mat_write ( filename, n, 1, x );

  u2 = new double[n*2];
  for ( i = 0; i < n; i++ )
  {
    u2[i+0*n] = u[i];
  }
  for ( i = 0; i < n; i++ )
  {
    u2[i+1*n] = uexact[i];
  }

  filename = "fd1d_bvp_test02_values.txt";
  r8mat_write ( filename, n, 2, u2 );

  delete [] u;
  delete [] u2;
  delete [] uexact;
  delete [] x;

  return;
}
//****************************************************************************80

void fd1d_bvp_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_BVP_TEST03 carries out test case #3.
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
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int i;
  int n = 11;
  double *u;
  double *u2;
  double *uexact;
  double *x;
  double x1 = 0.0;
  double x2 = 1.0;

  cout << "\n";
  cout << "FD1D_BVP_TEST03\n";
  cout << "  A1(X)  = 1.0\n";
  cout << "  A1''(X) = 0.0\n";
  cout << "  C3(X)  = 2.0 * X\n";
  cout << "  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
  cout << "  X1 = " << x1 << "\n";
  cout << "  X2 = " << x2 << "\n";

  x = r8vec_even ( n, x1, x2 );

  u = fd1d_bvp ( n, a1, a1prime, c3, f3, x );

  uexact = exact1 ( n, x );

  cout << "\n";
  cout << "     I         X        U             Uexact         Error\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(12) << u[i]
         << "  " << setw(12) << uexact[i]
         << "  " << setw(14) << r8_abs ( u[i] - uexact[i] ) << "\n";
  }
//
//  Write the data to files.
//
  filename = "fd1d_bvp_test03_nodes.txt";
  r8mat_write ( filename, n, 1, x );

  u2 = new double[n*2];
  for ( i = 0; i < n; i++ )
  {
    u2[i+0*n] = u[i];
  }
  for ( i = 0; i < n; i++ )
  {
    u2[i+1*n] = uexact[i];
  }

  filename = "fd1d_bvp_test03_values.txt";
  r8mat_write ( filename, n, 2, u2 );

  delete [] u;
  delete [] u2;
  delete [] uexact;
  delete [] x;

  return;
}
//****************************************************************************80

void fd1d_bvp_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_BVP_TEST04 carries out test case #4.
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
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int i;
  int n = 11;
  double *u;
  double *u2;
  double *uexact;
  double *x;
  double x1 = 0.0;
  double x2 = 1.0;

  cout << "\n";
  cout << "FD1D_BVP_TEST04\n";
  cout << "  A2(X)  = 1.0 + X * X\n";
  cout << "  A2''(X) = 2.0 * X\n";
  cout << "  C1(X)  = 0.0\n";
  cout << "  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
  cout << "  X1 = " << x1 << "\n";
  cout << "  X2 = " << x2 << "\n";

  x = r8vec_even ( n, x1, x2 );

  u = fd1d_bvp ( n, a2, a2prime, c1, f4, x );

  uexact = exact1 ( n, x );

  cout << "\n";
  cout << "     I         X        U             Uexact         Error\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(12) << u[i]
         << "  " << setw(12) << uexact[i]
         << "  " << setw(14) << r8_abs ( u[i] - uexact[i] ) << "\n";
  }
//
//  Write the data to files.
//
  filename = "fd1d_bvp_test04_nodes.txt";
  r8mat_write ( filename, n, 1, x );

  u2 = new double[n*2];
  for ( i = 0; i < n; i++ )
  {
    u2[i+0*n] = u[i];
  }
  for ( i = 0; i < n; i++ )
  {
    u2[i+1*n] = uexact[i];
  }

  filename = "fd1d_bvp_test04_values.txt";
  r8mat_write ( filename, n, 2, u2 );

  delete [] u;
  delete [] u2;
  delete [] uexact;
  delete [] x;

  return;
}
//****************************************************************************80

void fd1d_bvp_test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_BVP_TEST05 carries out test case #5.
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
//    15 February 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int i;
  int n = 11;
  double *u;
  double *u2;
  double *uexact;
  double *x;
  double x1 = 0.0;
  double x2 = 1.0;

  cout << "\n";
  cout << "FD1D_BVP_TEST05\n";
  cout << "  A3(X)  = 1.0 + X * X for X <= 1/3\n";
  cout << "         = 7/9 + X     for      1/3 < X\n";
  cout << "  A3''(X) = 2.0 * X     for X <= 1/3\n";
  cout << "           1           for      1/3 < X\n";
  cout << "  C1(X)  = 0.0\n";
  cout << "  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n";
  cout << "                       for X <= 1/3\n";
  cout << "         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) * exp ( X )\n";
  cout << "                       for      1/3 <= X\n";
  cout << "  U1(X)  = X * ( 1 - X ) * exp ( X )\n";
  cout << "\n";
  cout << "  Number of nodes = " << n << "\n";
  cout << "  X1 = " << x1 << "\n";
  cout << "  X2 = " << x2 << "\n";

  x = r8vec_even ( n, x1, x2 );

  u = fd1d_bvp ( n, a3, a3prime, c1, f5, x );

  uexact = exact1 ( n, x );

  cout << "\n";
  cout << "     I         X        U             Uexact         Error\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(12) << u[i]
         << "  " << setw(12) << uexact[i]
         << "  " << setw(14) << r8_abs ( u[i] - uexact[i] ) << "\n";
  }
//
//  Write the data to files.
//
  filename = "fd1d_bvp_test05_nodes.txt";
  r8mat_write ( filename, n, 1, x );

  u2 = new double[n*2];
  for ( i = 0; i < n; i++ )
  {
    u2[i+0*n] = u[i];
  }
  for ( i = 0; i < n; i++ )
  {
    u2[i+1*n] = uexact[i];
  }

  filename = "fd1d_bvp_test05_values.txt";
  r8mat_write ( filename, n, 2, u2 );

  delete [] u;
  delete [] u2;
  delete [] uexact;
  delete [] x;

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
//    17 May 2009
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

double a1prime ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    A1PRIME evaluates A' function #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double A1PRIME, the value of A'(X).
//
{
  double value;

  value = 0.0;

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
//    17 May 2009
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

double a2prime ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    A2PRIME evaluates A' function #2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double A2PRIME, the value of A'(X).
//
{
  double value;

  value = 2.0 * x;

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
//    17 May 2009
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

double a3prime ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    A3PRIME evaluates A' function #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double A3PRIME, the value of A'(X).
//
{
  double value;

  if ( x <= 1.0 / 3.0 )
  {
    value = 2.0 * x;
  }
  else
  {
    value = 1.0;
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
//    17 May 2009
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
//    17 May 2009
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
//    17 May 2009
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
//    17 May 2009
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
//    17 May 2009
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
//    17 May 2009
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
//    17 May 2009
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
//    17 May 2009
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

double *exact1 ( int n, double x[] )

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
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double EXACT1[N], the values of U(X(1:N)).
//
{
  int i;
  double *uexact;

  uexact = new double[n];

  for ( i = 0; i < n; i++ )
  {
    uexact[i] = x[i] * ( 1.0 - x[i] ) * exp ( x[i] );
  }

  return uexact;
}
//****************************************************************************80

double *exact2 ( int n, double x[] )

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
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double EXACT2[N], the values of U(X(1:N)).
//
{
  int i;
  double *uexact;

  uexact = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= 2.0 / 3.0 )
    {
      uexact[i] = x[i] * ( 1.0 - x[i] ) * exp ( x[i] );
    }
    else
    {
      uexact[i] = x[i] * ( 1.0 - x[i] ) * exp ( 2.0 / 3.0 );
    }
  }

  return uexact;
}
//****************************************************************************80

double *exact3 ( int n, double x[] )

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
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double EXACT3[N], the values of U(X(1:N)).
//
{
  int i;
  double *uexact;

  uexact = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= 2.0 / 3.0 )
    {
      uexact[i] = x[i] * ( 1.0 - x[i] ) * exp ( x[i] );
    }
    else
    {
      uexact[i] = x[i] * ( 1.0 - x[i] );
    }
  }

  return uexact;
}
