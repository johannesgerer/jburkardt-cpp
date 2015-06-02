# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fd1d_heat_explicit.hpp"

int main ( );
void fd1d_heat_explicit_test01 ( );
void bc_test01 ( int x_num, double x[], double t, double h[] );
double *ic_test01 ( int x_num, double x[], double t );
double *rhs_test01 ( int x_num, double x[], double t );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FD1D_HEAT_EXPLICIT_PRB.
//
//  Discussion:
//
//    FD1D_HEAT_EXPLICIT_TEST tests the FD1D_HEAT_EXPLICIT library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FD1D_HEAT_EXPLICIT_TEST:\n";
  cout << "  C++ version.\n";
  cout << "  Test the FD1D_HEAT_EXPLICIT library.\n";

  fd1d_heat_explicit_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FD1D_HEAT_EXPLICIT_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void fd1d_heat_explicit_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    FD1D_HEAT_EXPLICIT_TEST01 does a simple test problem
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double cfl;
  double dt;
  double *h;
  double *h_new;
  double *hmat;
  int i;
  int j;
  double k;
  double *t;
  double t_max;
  double t_min;
  int t_num;
  double *x;
  double x_max;
  double x_min;
  int x_num;

  cout << "\n";
  cout << "FD1D_HEAT_EXPLICIT_TEST01:\n";
  cout << "  Compute an approximate solution to the time-dependent\n";
  cout << "  one dimensional heat equation:\n";
  cout << "\n";
  cout << "    dH/dt - K * d2H/dx2 = f(x,t)\n";
  cout << "\n";
  cout << "  Run a simple test case.\n";
//
//  Heat coefficient.
//
  k = 0.002;
//
//  X_NUM is the number of equally spaced nodes to use between 0 and 1.
//
  x_num = 21;
  x_min = 0.0;
  x_max = 1.0;
  x = r8vec_linspace_new ( x_num, x_min, x_max );
//
//  T_NUM is the number of equally spaced time points between 0 and 10.0.
//
  t_num = 201;
  t_min = 0.0;
  t_max = 80.0;
  dt = ( t_max - t_min ) / ( double ) ( t_num - 1 );
  t = r8vec_linspace_new ( t_num, t_min, t_max );
//
//  Get the CFL coefficient.
//
  cfl = fd1d_heat_explicit_cfl ( k, t_num, t_min, t_max, x_num, x_min, x_max );
//
//  Running the code produces an array H of temperatures H(t,x),
//  and vectors x and t.
//
  h = ic_test01 ( x_num, x, t[0] );
  bc_test01 ( x_num, x, t[0], h );

  hmat = new double[x_num*t_num];

  j = 0;
  for ( i = 0; i < x_num; i++ )
  {
    hmat[i+j*x_num] = h[i];
  }

  for ( j = 1; j < t_num; j++ )
  {
    h_new = fd1d_heat_explicit ( x_num, x, t[j-1], dt, cfl, rhs_test01, bc_test01, h );

    for ( i = 0; i < x_num; i++ )
    {
      hmat[i+j*x_num] = h_new[i];
      h[i] = h_new[i];
    }
    delete [] h_new;
  }
//
//  Write the data to files.
//
  r8mat_write ( "h_test01.txt", x_num, t_num, hmat );
  r8vec_write ( "t_test01.txt", t_num, t );
  r8vec_write ( "x_test01.txt", x_num, x );

  delete [] h;
  delete [] hmat;
  delete [] t;
  delete [] x;

  return;
}
//****************************************************************************80

void bc_test01 ( int x_num, double x[], double t, double h[] )

//****************************************************************************80
//
//  Purpose:
//
//    BC_TEST01 evaluates the boundary conditions for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X[X_NUM], the node coordinates.
//
//    Input, double T, the current time.
//
//    Input, double H[X_NUM], the current heat values.
//
//    Output, double H[X_NUM], the current heat values, after boundary
//    conditions have been imposed.
//
{
  h[0]  = 90.0;
  h[x_num-1] = 70.0;

  return;
}
//****************************************************************************80

double *ic_test01 ( int x_num, double x[], double t )

//****************************************************************************80
//
//  Purpose:
//
//    IC_TEST01 evaluates the initial condition for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X[X_NUM], the node coordinates.
//
//    Input, double T, the initial time.
//
//    Output, double H[X_NUM], the heat values at the initial time.
//
{
  int j;
  double *h;

  h = new double[x_num];

  for ( j = 0; j < x_num; j++ )
  {
    h[j] = 50.0;
  }
  return h;
}
//****************************************************************************80

double *rhs_test01 ( int x_num, double x[], double t )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_TEST01 evaluates the right hand side for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X[X_NUM], the node coordinates.
//
//    Input, double T, the current time.
//
//    Output, double RHS_TEST01[X_NUM], the source term.
//
{
  int i;
  double *value;

  value = new double[x_num];

  for ( i = 0; i < x_num; i++ )
  {
    value[i] = 0.0;
  }
  return value;
}
