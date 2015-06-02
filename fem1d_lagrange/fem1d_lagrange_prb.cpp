# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <iostream>

using namespace std;

# include "fem1d_lagrange.hpp"

int main ( );
void legendre_set_test ( );
void lagrange_value_test ( );
void lagrange_derivative_test ( );
void fem1d_lagrange_stiffness_test ( int x_num, int q_num );
double f ( double x );
double *exact ( int x_num, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_LAGRANGE_PRB.
//
//  Discussion:
//
//    FEM1D_LAGRANGE_PRB tests FEM1D_LAGRANGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  int q_num;
  int x_num;

  timestamp ( );
  cout << "\n";
  cout << "FEM1D_LAGRANGE_PRB\n";
  cout << "  C version.\n";
  cout << "  Test the FEM1D_LAGRANGE library.\n";

  legendre_set_test ( );
  lagrange_value_test ( );
  lagrange_derivative_test ( );

  x_num = 11;
  q_num = 5;
  fem1d_lagrange_stiffness_test ( x_num, q_num );

  x_num = 11;
  q_num = 10;
  fem1d_lagrange_stiffness_test ( x_num, q_num );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM1D_LAGRANGE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void legendre_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_TEST tests LEGENDRE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  double e1;
  double e2;
  double e3;
  int i;
  int n;
  double *w;
  double *x;

  cout << "\n";
  cout << "LEGENDRE_SET_TEST\n";
  cout << "  LEGENDRE_SET returns points and weights of\n";
  cout << "  Gauss-Legendre quadrature rules.\n";
  cout << "\n";
  cout << "   N               1             X^4           Runge\n";
  cout << "\n";

  for ( n = 1; n <= 10; n++ )
  {
    x = new double[n];
    w = new double[n];

    legendre_set ( n, x, w );
    e1 = 0.0;
    e2 = 0.0;
    e3 = 0.0;
    for ( i = 0; i < n; i++ )
    {
      e1 = e1 + w[i];
      e2 = e2 + w[i] * pow ( x[i], 4 );
      e3 = e3 + w[i] / ( 1.0 + 25.0 * x[i] * x[i] );
    }
    cout << "  " << setw(2) << n
         << "  " << setw(14) << e1
         << "  " << setw(14) << e2
         << "  " << setw(14) << e3 << "\n";
    free ( w );
    free ( x );
  }

  return;
}
//****************************************************************************80

void lagrange_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_VALUE_TEST tests LAGRANGE_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  double *li;
  int nd;
  int ni;
  double *xd;
  double xhi;
  double *xi;
  double xlo;

  cout << "\n";
  cout << "LAGRANGE_VALUE_TEST\n";
  cout << "  LAGRANGE_VALUE evaluates the Lagrange basis polynomials.\n";

  nd = 5;
  xlo = 0.0;
  xhi = ( double ) ( nd - 1 ) ;
  xd = r8vec_linspace_new ( nd, xlo, xhi );

  r8vec_print ( nd, xd, "  Lagrange basis points:" );
//
//  Evaluate the polynomials.
//
  cout << "\n";
  cout << "   I      X          L1(X)       L2(X)       L3(X)";
  cout << "       L4(X)       L5(X)\n";
  cout << "\n";
 
  ni = 2 * nd - 1;
  xi = r8vec_linspace_new ( ni, xlo, xhi );

  li = lagrange_value ( nd, xd, ni, xi );

  for ( i = 0; i < ni; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << xi[i];
    for ( j = 0; j < nd; j++ )
    {
      cout << "  " << setw(10) << li[i+j*ni];
    }
    cout << "\n";
  }

  delete [] li;
  delete [] xd;
  delete [] xi;

  return;
}
//****************************************************************************80

void lagrange_derivative_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_DERIVATIVE_TEST tests LAGRANGE_DERIVATIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  double *lpi;
  int nd;
  int ni;
  double *xd;
  double xhi;
  double *xi;
  double xlo;

  cout << "\n";
  cout << "LAGRANGE_DERIVATIVE_TEST\n";
  cout << "  LAGRANGE_DERIVATIVE evaluates the Lagrange basis derivative.\n";

  nd = 5;
  xlo = 0.0;
  xhi = ( double ) ( nd - 1 );
  xd = r8vec_linspace_new ( nd, xlo, xhi );

  r8vec_print ( nd, xd, "  Lagrange basis points:" );
//
//  Evaluate the polynomials.
//
  cout << "\n";
  cout << "   I      X         L1'(X)      L2'(X)      L3'(X)";
  cout << "      L4'(X)      L5'(X)\n";
  cout << "\n";
 
  ni = 2 * nd - 1;
  xi = r8vec_linspace_new ( ni, xlo, xhi );
  lpi = lagrange_derivative ( nd, xd, ni, xi );

  for ( i = 0; i < ni; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << xi[i];
    for ( j = 0; j < nd; j++ )
    {
      cout << "  " << setw(10) << lpi[i+j*ni];
    }
    cout << "\n";
  }

  delete [] lpi;
  delete [] xd;
  delete [] xi;

  return;
}
//****************************************************************************80

void fem1d_lagrange_stiffness_test ( int x_num, int q_num )

//****************************************************************************80
//
//  Purpose:
//
//    FEM1D_LAGRANGE_STIFFNESS_TEST tests FEM1D_LAGRANGE_STIFFNESS.
//
//  Discussion:
//
//    The results are very sensitive to the quadrature rule.
//
//    In particular, if X_NUM points are used, the mass matrix will
//    involve integrals of polynomials of degree 2*(X_NUM-1), so the
//    quadrature rule should use at least Q_NUM = X_NUM - 1 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, int Q_NUM, the number of quadrature points.
//
{
  double *a;
  double *b;
  int i;
  int info;
  int j;
  double *k;
  double *m;
  double *u;
  double *u_e;
  double *x;
  double x_hi;
  double x_lo;

  cout << "\n";
  cout << "FEM1D_LAGRANGE_STIFFNESS_TEST\n";
  cout << "  FEM1D_LAGRANGE_STIFFNESS computes the stiffness matrix,\n";
  cout << "  the mass matrix, and right hand side vector for a\n";
  cout << "  finite element problem using Lagrange interpolation\n";
  cout << "  basis polynomials.\n";
  cout << "\n";
  cout << "  Solving:\n";
  cout << "    -u''+u=x on 0 < x < 1\n";
  cout << "    u(0) = u(1) = 0\n";
  cout << "  Exact solution:\n";
  cout << "    u(x) = x - sinh(x)/sinh(1)\n";
  cout << "\n";
  cout << "  Number of mesh points = " << x_num << "\n";
  cout << "  Number of quadrature points = " << q_num << "\n";

  x_lo = 0.0;
  x_hi = 1.0;
  x = r8vec_linspace_new ( x_num, x_lo, x_hi );

  a = new double[x_num*x_num];
  m = new double[x_num*x_num];
  b = new double[x_num];

  fem1d_lagrange_stiffness ( x_num, x, q_num, f, a, m, b );

  k = new double[x_num*x_num];

  for ( j = 0; j < x_num; j++ )
  {
    for ( i = 0; i < x_num; i++ )
    {
      k[i+j*x_num] = a[i+j*x_num] + m[i+j*x_num];
    }
  }
  for ( j = 0; j < x_num; j++ )
  {
    k[0+j*x_num] = 0.0;
  }
  k[0+0*x_num] = 1.0;
  b[0] = 0.0;

  for ( j = 0; j < x_num; j++ )
  {
    k[x_num-1+j*x_num] = 0.0;
  }
  k[x_num-1+(x_num-1)*x_num] = 1.0;
  b[x_num-1] = 0.0;

  u = r8mat_fs_new ( x_num, k, b );

  u_e = exact ( x_num, x );

  cout << "\n";
  cout << "   I      X             U              U(exact)         Error\n";
  cout << "\n";

  for ( i = 0; i < x_num; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(14) << u[i]
         << "  " << setw(14) << u_e[i]
         << "  " << setw(14) << fabs ( u[i] - u_e[i] ) << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] k;
  delete [] m;
  delete [] u;
  delete [] u_e;
  delete [] x;

  return;
}
//****************************************************************************80

double f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F evaluates the right hand side function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double F, the value of the right hand side at X.
//
{
  double value;

  value = x;

  return value;
}
//****************************************************************************80

double *exact ( int x_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT returns the exact solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X[X_NUM], the nodes.
//
//    Output, double UE[X_NUM], the exact solution at the nodes.
//
{
  double *ue;
  int x_i;

  ue = new double[x_num];

  for ( x_i = 0; x_i < x_num; x_i++ )
  {
    ue[x_i] = x[x_i] - sinh ( x[x_i] ) / sinh ( 1.0 );
  }
  return ue;
}

