// File recommented by recomment.C
// on Jul  7 2014 at 20:37:06.
//
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem2d_bvp_serene.hpp"

int main ( );

void test01 ( );
double a1 ( double x, double y );
double c1 ( double x, double y );
double exact1 ( double x, double y );
double exact_ux1 ( double x, double y );
double exact_uy1 ( double x, double y );
double f1 ( double x, double y );
void test02 ( );
void test03 ( );
double a3 ( double x, double y );
double c3 ( double x, double y );
double exact3 ( double x, double y );
double exact_ux3 ( double x, double y );
double exact_uy3 ( double x, double y );
double f3 ( double x, double y );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_BVP_SERENE_PRB.
//
//  Discussion:
//
//    FEM2D_BVP_SERENE_PRB tests the FEM2D_BVP_SERENE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FEM2D_BVP_SERENE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FEM2D_BVP_SERENE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM2D_BVP_SERENE_PRB\n";
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
//    Use A1, C1, F1, EXACT1, EXACT_UX1, EXACT_UX1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
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
  int inc;
  int j;
  int k;
  int mn;
  int nx = 5;
  int ny = 5;
  int show11;
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
  x = r8vec_linspace_new ( nx, x_first, x_last );

  y_first = 0.0;
  y_last = 1.0;
  y = r8vec_linspace_new ( ny, y_first, y_last );

  show11 = 0;
  u = fem2d_bvp_serene ( nx, ny, a1, c1, f1, x, y, show11 );

  if ( nx * ny <= 25 )
  {
    cout << "\n";
    cout << "     I     J    X         Y               U               Uexact     Error\n";
    cout << "\n";

    k = 0;

    for ( j = 0; j < ny; j++ )
    {
      if ( ( j % 2 ) == 0 )
      {
        inc = 1;
      }
      else
      {
        inc = 2;
      }
      for ( i = 0; i < nx; i = i + inc )
      {
        uexact = exact1 ( x[i], y[j] );
        cout << setw(4) << i << "  "
             << setw(4) << j << "  "
             << setw(8) << x[i] << "  "
             << setw(8) << y[j] << "  "
             << setw(14) << u[k] << "  "
             << setw(14) << uexact << "  "
             << setw(14) << fabs ( u[k] - uexact ) << "\n";
        k = k + 1;
      }
    }
  }

  e1 = fem2d_l1_error_serene ( nx, ny, x, y, u, exact1 );
  e2 = fem2d_l2_error_serene ( nx, ny, x, y, u, exact1 );
  h1s = fem2d_h1s_error_serene ( nx, ny, x, y, u, exact_ux1, exact_uy1 );

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
//    07 July 2014
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
//    07 July 2014
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
//    07 July 2014
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
//    07 July 2014
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
//    07 July 2014
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
//    07 July 2014
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
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 checks the basis functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int seed;
  double *v;
  double *vx;
  double *vy;
  double xe;
  double xq;
  double xw;
  double xx[8] = { 2.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0 };
  double yn;
  double yq;
  double ys;
  double yy[8] = { 5.0, 5.0, 5.0, 4.0, 3.0, 3.0, 3.0, 4.0 };

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Basis function checks.\n";
//
//  Check that V is identity matrix at nodes.
//
  cout << "\n";
  cout << "  The matrix Aij = V(j)(X(i),Y(i)) should be the identity.\n";
  cout << "\n";

  xw = 0.0;
  ys = 3.0;
  xe = 2.0;
  yn = 5.0;

  for ( i = 0; i < 8; i++ )
  {
    xq = xx[i];
    yq = yy[i];
    v = basis_serene ( xq, yq, xw, ys, xe, yn, xx, yy );
    for ( j = 0; j < 8; j++ )
    {
      cout << setw(10) << v[j];
    }
    cout << "\n";
    delete [] v;
  }
//
//  Check that VX and VY sum to zero anywhere.
//
  cout << "\n";
  cout << "  The vectors dVdX(1:8)(X,Y) and dVdY(1:8)(X,Y)\n";
  cout << "  should both sum to zero for any (X,Y).\n";

  seed = 123456789;
  xq = 2.0 * r8_uniform_01 ( seed );
  yq = 3.0 + 2.0 * r8_uniform_01 ( seed );
  xw = 0.0;
  ys = 3.0;
  xe = 2.0;
  yn = 5.0;

  vx = basis_dx_serene ( xq, yq, xw, ys, xe, yn, xx, yy );
  vy = basis_dy_serene ( xq, yq, xw, ys, xe, yn, xx, yy );

  cout << "\n";
  cout << "  Random evaluation point is (" << xq << "," << yq << ")\n";
  cout << "\n";
  cout << "              dVdX        dVdY\n";
  cout << "\n";
  for ( i = 0; i < 8; i++ )
  {
    cout << setw(6) << i << "  "
         << setw(10) << vx[i] << "  "
         << setw(10) << vy[i] << "\n";
  }
  cout << "\n";
  cout << "  Sum:  "
       << setw(10) << r8vec_sum ( 8, vx ) << "  "
       << setw(10) << r8vec_sum ( 8, vy ) << "\n";

  delete [] vx;
  delete [] vy;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 carries out test case #3.
//
//  Discussion:
//
//    Use A3, C3, F3, EXACT3, EXACT_UX3, EXACT_UX3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *amat;
  double e1;
  double e2;
  double h1s;
  int i;
  int inc;
  int j;
  int k;
  int mn;
  int nx = 5;
  int ny = 5;
  double scale;
  int show11;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;
  double *y;
  double y_first;
  double y_last;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Solve - del ( A del U ) + C U = F \n";
  cout << "  on the unit square with zero boundary conditions.\n";
  cout << "  A1(X,Y) = 0.0\n";
  cout << "  C1(X,Y) = 1.0\n";
  cout << "  F1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )\n";
  cout << "  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )\n";
  cout << "\n";
  cout << "  This example is contrived so that the system matrix\n";
  cout << "  is the WATHEN matrix.\n";
  cout << "\n";
  cout << "  Number of X grid values NX = " << nx << "\n";
  cout << "  Number of Y grid values NY = " << ny << "\n";
  mn = nx * ny;
//
//  Geometry definitions.
//
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_linspace_new ( nx, x_first, x_last );

  y_first = 0.0;
  y_last = 1.0;
  y = r8vec_linspace_new ( ny, y_first, y_last );

  show11 = 1;
  u = fem2d_bvp_serene ( nx, ny, a3, c3, f3, x, y, show11 );

  if ( nx * ny <= 25 )
  {
    cout << "\n";
    cout << "     I     J    X         Y               U               Uexact     Error\n";
    cout << "\n";

    k = 0;

    for ( j = 0; j < ny; j++ )
    {
      if ( ( j % 2 ) == 0 )
      {
        inc = 1;
      }
      else
      {
        inc = 2;
      }
      for ( i = 0; i < nx; i = i + inc )
      {
        uexact = exact3 ( x[i], y[j] );
        cout << setw(4) << i << "  "
             << setw(4) << j << "  "
             << setw(8) << x[i] << "  "
             << setw(8) << y[j] << "  "
             << setw(14) << u[k] << "  "
             << setw(14) << uexact << "  "
             << setw(14) << fabs ( u[k] - uexact ) << "\n";
        k = k + 1;
      }
    }
  }

  e1 = fem2d_l1_error_serene ( nx, ny, x, y, u, exact3 );
  e2 = fem2d_l2_error_serene ( nx, ny, x, y, u, exact3 );
  h1s = fem2d_h1s_error_serene ( nx, ny, x, y, u, exact_ux3, exact_uy3 );

  cout << "\n";
  cout << "  l1 norm of error  = " << e1 << "\n";
  cout << "  L2 norm of error  = " << e2 << "\n";
  cout << "  Seminorm of error = " << h1s << "\n";

  delete [] u;
  delete [] x;
  delete [] y;
//
//  Pull out the Wathen matrix from MATLAB.
//  It will have been multiplied by a random scale factor.
//  While my numbering scheme is
//    3  2  1
//    4     8
//    5  6  7
//  the numbering scheme used here is 
//    1  2  3
//    4     5
//    6  7  8
//    
  amat = wathen ( 1, 1, 8 );
 
  scale = 0.5 * amat[0+2*8];
  for ( j = 0; j < 8; j++ )
  {
    for ( i = 0; i < 8; i++ )
    {
      amat[i+j*8] = amat[i+j*8] / scale;
    }
  }
 
  r8mat_print ( 8, 8, amat, "  Wathen matrix:" );

  delete [] amat;

  return;
}
//****************************************************************************80

double a3 ( double x, double y )

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
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double A3, the value of A(X,Y).
//
{
  double value;

  value = 0.0;

  return value;
}
//****************************************************************************80

double c3 ( double x, double y )

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
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double C3, the value of C(X,Y).
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double exact3 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT3 evaluates exact solution #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double EXACT3, the value of U(X,Y).
//
{
  double value;

  value = x * ( 1.0 - x ) * y * ( 1.0 - y );

  return value;
}
//****************************************************************************80

double exact_ux3 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT_UX3 evaluates the derivative of exact solution #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double EXACT_UX3, the value of dUdX(X,Y).
//
{
  double value;

  value = ( 1.0 - 2.0 * x ) * ( y - y * y );

  return value;
}
//****************************************************************************80

double exact_uy3 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT_UY3 evaluates the derivative of exact solution #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double EXACT_UY3, the value of dUdY(X,Y).
//
{
  double value;

  value = ( x - x * x ) * ( 1.0 - 2.0 * y );

  return value;
}
//****************************************************************************80

double f3 ( double x, double y )

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
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double F3, the value of F(X,Y).
//
{
  double value;

  value = x * ( 1.0 - x ) * y * ( 1.0 - y );

  return value;
}
