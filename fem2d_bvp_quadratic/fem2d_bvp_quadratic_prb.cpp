
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem2d_bvp_quadratic.hpp"

int main ( );

void test01 ( );
double a1 ( double x, double y );
double c1 ( double x, double y );
double exact1 ( double x, double y );
double exact_ux1 ( double x, double y );
double exact_uy1 ( double x, double y );
double f1 ( double x, double y );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_BVP_QUADRATIC_PRB.
//
//  Discussion:
//
//    FEM2D_BVP_QUADRATIC_PRB tests the FEM2D_BVP_QUADRATIC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FEM2D_BVP_QUADRATIC_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FEM2D_BVP_QUADRATIC library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM2D_BVP_QUADRATIC_PRB\n";
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
//    TEST01 carries out test case #1.
//
//  Discussion:
//
//    Use A1, C1, F1, EXACT1, EXACT_UX1, EXACT_UX2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double e1;
  double e2;
  double h1s;
  int i;
  int j;
  int mn;
  int nx = 3;
  int ny = 3;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;
  double *y;
  double y_first;
  double y_last;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Solve - del ( A del U ) + C U = F \n";
  cout << "  on the unit square with zero boundary conditions.\n";
  cout << "  A1(X,Y) = 1.0\n";
  cout << "  C1(X,Y) = 0.0\n";
  cout << "  F1(X,Y) = 2*X*(1-X)+2*Y*(1-Y)\n";
  cout << "  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )\n";
  cout << "\n";
  cout << "  Number of X grid values NX = " << nx << "\n";
  cout << "  Number of Y grid values NY = " << ny << "\n";
  mn = nx * ny;
//
//  Geometry definitions.
//
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( nx, x_first, x_last );

  y_first = 0.0;
  y_last = 1.0;
  y = r8vec_even_new ( ny, y_first, y_last );

  u = fem2d_bvp_quadratic ( nx, ny, a1, c1, f1, x, y );

  cout << "\n";
  cout << "     I     J    X         Y         U         Uexact    Error\n";
  cout << "\n";

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      uexact = exact1 ( x[i], y[j] );
      cout << setw(6) << i << "  "
           << setw(4) << j << "  "
           << setw(8) << x[i] << "  "
           << setw(8) << y[j] << "  "
           << setw(14) << u[i+j*nx] << "  "
           << setw(14) << uexact << "  "
           << setw(14) << fabs ( u[i+j*nx] - uexact ) << "\n";
    }
  }

  e1 = fem2d_l1_error ( nx, ny, x, y, u, exact1 );
  e2 = fem2d_l2_error_quadratic ( nx, ny, x, y, u, exact1 );
  h1s = fem2d_h1s_error_quadratic ( nx, ny, x, y, u, exact_ux1, exact_uy1 );
  cout << "\n";
  cout << "  l1 norm of error  = " << e1 << "\n";
  cout << "  L2 norm of error  = " << e2 << "\n";
  cout << "  Seminorm of error = " << h1s << "\n";

  delete [] u;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

double a1 ( double x, double y )

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
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double A1, the value of A(X,Y).
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double c1 ( double x, double y )

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
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double C1, the value of C(X,Y).
//
{
  double value;

  value = 0.0;

  return value;
}
//****************************************************************************80

double exact1 ( double x, double y )

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
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double EXACT1, the value of U(X,Y).
//
{
  double value;

  value = x * ( 1.0 - x ) * y * ( 1.0 - y );

  return value;
}
//****************************************************************************80

double exact_ux1 ( double x, double y )

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
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double EXACT_UX1, the value of dUdX(X,Y).
//
{
  double value;

  value = ( 1.0 - 2.0 * x ) * ( y - y * y );

  return value;
}
//****************************************************************************80

double exact_uy1 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT_UY1 evaluates the derivative of exact solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double EXACT_UY1, the value of dUdY(X,Y).
//
{
  double value;

  value = ( x - x * x ) * ( 1.0 - 2.0 * y );

  return value;
}
//****************************************************************************80

double f1 ( double x, double y )

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
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double F1, the value of F(X,Y).
//
{
  double value;

  value = 2.0 * x * ( 1.0 - x ) 
        + 2.0 * y * ( 1.0 - y );

  return value;
}

