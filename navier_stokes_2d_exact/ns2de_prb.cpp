# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "ns2de.hpp"

int main ( );
void uvp_taylor_test ( );
void uvp_taylor_test2 ( );
void rhs_taylor_test ( );
void resid_taylor_test ( );
void gnuplot_taylor_test ( );
void parameter_taylor_test ( );
void uvp_spiral_test ( );
void uvp_spiral_test2 ( );
void rhs_spiral_test ( );
void resid_spiral_test ( );
void gnuplot_spiral_test ( );
void parameter_spiral_test ( );
void uvp_lucas_test ( );
void uvp_lucas_test2 ( );
void rhs_lucas_test ( );
void resid_lucas_test ( );
void gnuplot_lucas_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    NS2DE_TEST tests the NS2DE library.
//
//  Location:
//
//    http://people.sc.fsu.edu/~jburkardt/cpp_src/navier_stokes_2d_exact/ns2de_prb.cpp
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "NS2DE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the NS2DE library.\n";
//
//  Taylor Vortex Flow.
//
  uvp_taylor_test ( );
  uvp_taylor_test2 ( );
  rhs_taylor_test ( );
  resid_taylor_test ( );
  gnuplot_taylor_test ( );
  parameter_taylor_test ( );
//
//  Spiral Flow.
//
  uvp_spiral_test ( );
  uvp_spiral_test2 ( );
  rhs_spiral_test ( );
  resid_spiral_test ( );
  gnuplot_spiral_test ( );
  parameter_spiral_test ( );
//
//  Lucas Bystricky Flow.
//
  uvp_lucas_test ( );
  uvp_lucas_test2 ( );
  rhs_lucas_test ( );
  resid_lucas_test ( );
  gnuplot_lucas_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "NS2DE_TEST\n";
  cout << "  Normal end of execution.\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void uvp_taylor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_TAYLOR_TEST samples the solution at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "UVP_TAYLOR_TEST\n";
  cout << "  Taylor Vortex Flow:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  at the initial time T = 0, using a region that is\n";
  cout << "  the square centered at (1.5,1.5) with 'radius' 1.0,\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  p = new double[n];
  u = new double[n];
  v = new double[n];

  r8_lo = 0.5;
  r8_hi = +2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:"
       << "  " << setw(14) << r8vec_min ( n, u )
       << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:"
       << "  " << setw(14) << r8vec_min ( n, v )
       << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:"
       << "  " << setw(14) << r8vec_min ( n, p )
       << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void uvp_taylor_test2 ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_TAYLOR_TEST2 samples the solution on the boundary at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double nu;
  double *p;
  double r8_hi;
  double r8_lo;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double *y;

  nu = 1.0;
  rho = 1.0;
  t = 0.0;

  cout << "\n";
  cout << "UVP_TAYLOR_TEST2\n";
  cout << "  Taylor Vortex Flow:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  on the boundary\n";
  cout << "  at the initial time T = 0, using a region that is\n";
  cout << "  the square centered at (1.5,1.5) with 'radius' 1.0,\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 400;

  p = new double[n];
  u = new double[n];
  v = new double[n];
  x = new double[n];
  y = new double[n];

  r8_lo = 0.5;
  r8_hi = +2.5;

  r8vec_linspace ( 100, r8_lo, r8_hi, x );
  for ( i = 0; i < 100; i++ )
  {
    y[i] = r8_lo;
  }

  for ( i = 100; i < 200; i++ )
  {
    x[i] = r8_hi;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+100 );

  r8vec_linspace ( 100, r8_hi, r8_lo, x+200 );
  for ( i = 200; i < 300; i++ )
  {
    y[i] = r8_hi;
  }

  for ( i = 300; i < 400; i++ )
  {
    x[i] = r8_lo;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+300 );

  uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:"
       << "  " << setw(14) << r8vec_min ( n, u )
       << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:"
       << "  " << setw(14) << r8vec_min ( n, v )
       << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:"
       << "  " << setw(14) << r8vec_min ( n, p )
       << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void rhs_taylor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_TAYLOR_TEST samples the right hand side at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *f;
  double *g;
  double *h;
  int n;
  double nu;
  double rho;
  int seed;
  double t;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "RHS_TAYLOR_TEST\n";
  cout << "  Taylor Vortex Flow:\n";
  cout << "  Sample the Navier-Stokes right hand sides\n";
  cout << "  at the initial time T = 0, using a region that is\n";
  cout << "  the square centered at (1.5,1.5) with 'radius' 1.0,\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  f = new double[n];
  g = new double[n];
  h = new double[n];

  r8_lo = 0.5;
  r8_hi = +2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  rhs_taylor ( nu, rho, n, x, y, t, f, g, h );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  F:"
       << "  " << setw(14) << r8vec_min ( n, f )
       << "  " << setw(14) << r8vec_max ( n, f ) << "\n";
  cout << "  G:"
       << "  " << setw(14) << r8vec_min ( n, g )
       << "  " << setw(14) << r8vec_max ( n, g ) << "\n";
  cout << "  H:"
       << "  " << setw(14) << r8vec_min ( n, h )
       << "  " << setw(14) << r8vec_max ( n, h ) << "\n";

  delete [] f;
  delete [] g;
  delete [] h;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void resid_taylor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_TAYLOR_TEST samples the residual at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double nu;
  double *pr;
  double rho;
  int seed;
  double t;
  double *ur;
  double *vr;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "RESID_TAYLOR_TEST\n";
  cout << "  Taylor Vortex Flow:\n";
  cout << "  Sample the Navier-Stokes residuals\n";
  cout << "  at the initial time T = 0, using a region that is\n";
  cout << "  the square centered at (1.5,1.5) with 'radius' 1.0,\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  pr = new double[n];
  ur = new double[n];
  vr = new double[n];

  r8_lo = 0.5;
  r8_hi = +2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  resid_taylor ( nu, rho, n, x, y, t, ur, vr, pr );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  Ur:"
       << "  " << setw(14) << r8vec_amin ( n, ur )
       << "  " << setw(14) << r8vec_amax ( n, ur ) << "\n";
  cout << "  Vr:"
       << "  " << setw(14) << r8vec_amin ( n, vr )
       << "  " << setw(14) << r8vec_amax ( n, vr ) << "\n";
  cout << "  Pr:"
       << "  " << setw(14) << r8vec_amin ( n, pr )
       << "  " << setw(14) << r8vec_amax ( n, pr ) << "\n";

  delete [] pr;
  delete [] ur;
  delete [] vr;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void gnuplot_taylor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GNUPLOT_TAYLOR_TEST generates a field on a regular grid and plots it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  string header;
  int n;
  double nu;
  double *p;
  double rho;
  double s;
  int seed;
  double t;
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
  cout << "GNUPLOT_TAYLOR_TEST:\n";
  cout << "  Taylor Vortex Flow:\n";
  cout << "  Generate a Taylor vortex velocity field on a regular grid.\n";
  cout << "  Store in GNUPLOT data and command files.\n";

  x_lo = 0.5;
  x_hi = 2.5;

  y_lo = 0.5;
  y_hi = 2.5;

  x = new double[x_num*y_num];
  y = new double[x_num*y_num];

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  nu = 1.0;
  rho = 1.0;
  n = x_num * y_num;
  t = 0.0;

  u = new double[x_num*y_num];
  v = new double[x_num*y_num];
  p = new double[x_num*y_num];

  uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

  header = "taylor";
  s = 0.10;
  ns2de_gnuplot ( header, n, x, y, u, v, s );

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void parameter_taylor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PARAMETER_TAYLOR_TEST monitors solution norms for various values of NU, RHO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int n;
  double nu;
  double *p;
  double p_norm;
  double rho;
  int seed;
  double t;
  double *u;
  double u_norm;
  double *v;
  double v_norm;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  cout << "\n";
  cout << "PARAMETER_TAYLOR_TEST\n";
  cout << "  Taylor Vortex Flow:\n";
  cout << "  Monitor solution norms over time for various\n";
  cout << "  values of NU, RHO.\n";

  n = 1000;

  u = new double[n];
  v = new double[n];
  p = new double[n];

  r8_lo = 0.5;
  r8_hi = +2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
//
//  Vary RHO.
//
  cout << "\n";
  cout << "  RHO affects the pressure scaling.\n";
  cout << "\n";
  cout << "     RHO         NU           T     ||U||       ||V||       ||P||\n";
  cout << "\n";

  nu = 1.0;
  rho = 1.0;

  for ( j = 1; j <= 3; j++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      cout << "  " << setw(10) << rho
           << "  " << setw(10) << nu
           << "  " << setw(8) << t
           << "  " << setw(10) << u_norm
           << "  " << setw(10) << v_norm
           << "  " << setw(10) << p_norm << "\n";
    }
    cout << "\n";
    rho = rho / 100.0;
  }
/*
  Vary NU.
*/
  cout << "\n";
  cout << "  NU affects the time scaling.\n";
  cout << "\n";
  cout << "     RHO         NU           T     ||U||       ||V||       ||P||\n";
  cout << "\n";

  nu = 1.0;
  rho = 1.0;

  for ( i = 1; i <= 4; i++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      cout << "  " << setw(10) << rho
           << "  " << setw(10) << nu
           << "  " << setw(8) << t
           << "  " << setw(10) << u_norm
           << "  " << setw(10) << v_norm
           << "  " << setw(10) << p_norm << "\n";
    }

    cout << "\n";

    nu = nu / 10.0;
  }

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void uvp_spiral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_SPIRAL_TEST samples the Spiral Flow at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "UVP_SPIRAL_TEST\n";
  cout << "  Spiral Flow:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  at the initial time T = 0, over the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  p = new double[n];
  u = new double[n];
  v = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:"
       << "  " << setw(14) << r8vec_min ( n, u )
       << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:"
       << "  " << setw(14) << r8vec_min ( n, v )
       << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:"
       << "  " << setw(14) << r8vec_min ( n, p )
       << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void uvp_spiral_test2 ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_SPIRAL_TEST2 samples the Spiral Flow on the boundary.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double nu;
  double *p;
  double r8_hi;
  double r8_lo;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double *y;

  nu = 1.0;
  rho = 1.0;
  t = 0.0;

  cout << "\n";
  cout << "UVP_SPIRAL_TEST2\n";
  cout << "  Spiral Flow:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  on the boundary\n";
  cout << "  at the initial time T = 0, over the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 400;

  p = new double[n];
  u = new double[n];
  v = new double[n];
  x = new double[n];
  y = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;

  r8vec_linspace ( 100, r8_lo, r8_hi, x );
  for ( i = 0; i < 100; i++ )
  {
    y[i] = r8_lo;
  }

  for ( i = 100; i < 200; i++ )
  {
    x[i] = r8_hi;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+100 );

  r8vec_linspace ( 100, r8_hi, r8_lo, x+200 );
  for ( i = 200; i < 300; i++ )
  {
    y[i] = r8_hi;
  }

  for ( i = 300; i < 400; i++ )
  {
    x[i] = r8_lo;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+300 );

  uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:"
       << "  " << setw(14) << r8vec_min ( n, u )
       << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:"
       << "  " << setw(14) << r8vec_min ( n, v )
       << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:"
       << "  " << setw(14) << r8vec_min ( n, p )
       << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void rhs_spiral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_SPIRAL_TEST samples the right hand side at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *f;
  double *g;
  double *h;
  int n;
  double nu;
  double rho;
  int seed;
  double t;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "RHS_SPIRAL_TEST\n";
  cout << "  Spiral Flow:\n";
  cout << "  Sample the Navier-Stokes right hand sides\n";
  cout << "  at the initial time T = 0, using the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  f = new double[n];
  g = new double[n];
  h = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  rhs_spiral ( nu, rho, n, x, y, t, f, g, h );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  F:"
       << "  " << setw(14) << r8vec_min ( n, f )
       << "  " << setw(14) << r8vec_max ( n, f ) << "\n";
  cout << "  G:"
       << "  " << setw(14) << r8vec_min ( n, g )
       << "  " << setw(14) << r8vec_max ( n, g ) << "\n";
  cout << "  H:"
       << "  " << setw(14) << r8vec_min ( n, h )
       << "  " << setw(14) << r8vec_max ( n, h ) << "\n";

  delete [] f;
  delete [] g;
  delete [] h;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void resid_spiral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_SPIRAL_TEST samples the residual at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double nu;
  double *pr;
  double rho;
  int seed;
  double t;
  double *ur;
  double *vr;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "RESID_SPIRAL_TEST\n";
  cout << "  Spiral Flow:\n";
  cout << "  Sample the Navier-Stokes residuals\n";
  cout << "  at the initial time T = 0, on the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  pr = new double[n];
  ur = new double[n];
  vr = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  resid_spiral ( nu, rho, n, x, y, t, ur, vr, pr );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  Ur:"
       << "  " << setw(14) << r8vec_amin ( n, ur )
       << "  " << setw(14) << r8vec_amax ( n, ur ) << "\n";
  cout << "  Vr:"
       << "  " << setw(14) << r8vec_amin ( n, vr )
       << "  " << setw(14) << r8vec_amax ( n, vr ) << "\n";
  cout << "  Pr:"
       << "  " << setw(14) << r8vec_amin ( n, pr )
       << "  " << setw(14) << r8vec_amax ( n, pr ) << "\n";

  delete [] pr;
  delete [] ur;
  delete [] vr;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void gnuplot_spiral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GNUPLOT_SPIRAL_TEST generates a field on a regular grid and plots it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  string header;
  int n;
  double nu;
  double *p;
  double rho;
  double s;
  int seed;
  double t;
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
  cout << "GNUPLOT_SPIRAL_TEST:\n";
  cout << "  Spiral Flow:\n";
  cout << "  Generate a velocity field on a regular grid.\n";
  cout << "  Store in GNUPLOT data and command files.\n";

  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = new double[x_num*y_num];
  y = new double[x_num*y_num];

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  nu = 1.0;
  rho = 1.0;
  n = x_num * y_num;
  t = 0.0;

  u = new double[x_num*y_num];
  v = new double[x_num*y_num];
  p = new double[x_num*y_num];

  uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

  header = "spiral";
  s = 5.00;
  ns2de_gnuplot ( header, n, x, y, u, v, s );

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void parameter_spiral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PARAMETER_SPIRAL_TEST monitors solution norms for various values of NU, RHO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int n;
  double nu;
  double *p;
  double p_norm;
  double rho;
  int seed;
  double t;
  double *u;
  double u_norm;
  double *v;
  double v_norm;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  cout << "\n";
  cout << "PARAMETER_SPIRAL_TEST\n";
  cout << "  Spiral Flow:\n";
  cout << "  Monitor solution norms over time for various\n";
  cout << "  values of NU, RHO.\n";

  n = 1000;

  u = new double[n];
  v = new double[n];
  p = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
//
//  Vary RHO.
//
  cout << "\n";
  cout << "  RHO affects the pressure scaling.\n";
  cout << "\n";
  cout << "     RHO         NU           T     ||U||       ||V||       ||P||\n";
  cout << "\n";

  nu = 1.0;
  rho = 1.0;

  for ( j = 1; j <= 3; j++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      cout << "  " << setw(10) << rho
           << "  " << setw(10) << nu
           << "  " << setw(8) << t
           << "  " << setw(10) << u_norm
           << "  " << setw(10) << v_norm
           << "  " << setw(10) << p_norm << "\n";
    }
    cout << "\n";
    rho = rho / 100.0;
  }
/*
  Vary NU.
*/
  cout << "\n";
  cout << "  NU affects the time scaling.\n";
  cout << "\n";
  cout << "     RHO         NU           T     ||U||       ||V||       ||P||\n";
  cout << "\n";

  nu = 1.0;
  rho = 1.0;

  for ( i = 1; i <= 4; i++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      cout << "  " << setw(10) << rho
           << "  " << setw(10) << nu
           << "  " << setw(8) << t
           << "  " << setw(10) << u_norm
           << "  " << setw(10) << v_norm
           << "  " << setw(10) << p_norm << "\n";
    }

    cout << "\n";

    nu = nu / 10.0;
  }

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void uvp_lucas_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_LUCAS_TEST samples the Lucas Bystricky Flow at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "UVP_LUCAS_TEST\n";
  cout << "  Lucas Bystricky Flow:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  at the initial time T = 0, over the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  p = new double[n];
  u = new double[n];
  v = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  uvp_lucas ( nu, rho, n, x, y, t, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:"
       << "  " << setw(14) << r8vec_min ( n, u )
       << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:"
       << "  " << setw(14) << r8vec_min ( n, v )
       << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:"
       << "  " << setw(14) << r8vec_min ( n, p )
       << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void uvp_lucas_test2 ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_LUCAS_TEST2 samples the Lucas Bystricky Flow on the boundary.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double nu;
  double *p;
  double r8_hi;
  double r8_lo;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double *y;

  nu = 1.0;
  rho = 1.0;
  t = 0.0;

  cout << "\n";
  cout << "UVP_LUCAS_TEST2\n";
  cout << "  Lucas Bystricky Flow:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  on the boundary\n";
  cout << "  at the initial time T = 0, over the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 400;

  p = new double[n];
  u = new double[n];
  v = new double[n];
  x = new double[n];
  y = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;

  r8vec_linspace ( 100, r8_lo, r8_hi, x );
  for ( i = 0; i < 100; i++ )
  {
    y[i] = r8_lo;
  }

  for ( i = 100; i < 200; i++ )
  {
    x[i] = r8_hi;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+100 );

  r8vec_linspace ( 100, r8_hi, r8_lo, x+200 );
  for ( i = 200; i < 300; i++ )
  {
    y[i] = r8_hi;
  }

  for ( i = 300; i < 400; i++ )
  {
    x[i] = r8_lo;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+300 );

  uvp_lucas ( nu, rho, n, x, y, t, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:"
       << "  " << setw(14) << r8vec_min ( n, u )
       << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:"
       << "  " << setw(14) << r8vec_min ( n, v )
       << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:"
       << "  " << setw(14) << r8vec_min ( n, p )
       << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void rhs_lucas_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_LUCAS_TEST samples the right hand side at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *f;
  double *g;
  double *h;
  int n;
  double nu;
  double rho;
  int seed;
  double t;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "RHS_LUCAS_TEST\n";
  cout << "  Lucas Bystricky Flow:\n";
  cout << "  Sample the Navier-Stokes right hand sides\n";
  cout << "  at the initial time T = 0, using the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  f = new double[n];
  g = new double[n];
  h = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  rhs_lucas ( nu, rho, n, x, y, t, f, g, h );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  F:"
       << "  " << setw(14) << r8vec_min ( n, f )
       << "  " << setw(14) << r8vec_max ( n, f ) << "\n";
  cout << "  G:"
       << "  " << setw(14) << r8vec_min ( n, g )
       << "  " << setw(14) << r8vec_max ( n, g ) << "\n";
  cout << "  H:"
       << "  " << setw(14) << r8vec_min ( n, h )
       << "  " << setw(14) << r8vec_max ( n, h ) << "\n";

  delete [] f;
  delete [] g;
  delete [] h;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void resid_lucas_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_LUCAS_TEST samples the residual at the initial time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double nu;
  double *pr;
  double rho;
  int seed;
  double t;
  double *ur;
  double *vr;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  cout << "\n";
  cout << "RESID_LUCAS_TEST\n";
  cout << "  Lucas Bystricky Flow:\n";
  cout << "  Sample the Navier-Stokes residuals\n";
  cout << "  at the initial time T = 0, on the unit square.\n";
  cout << "  Kinematic viscosity NU = " << nu << "\n";
  cout << "  Fluid density RHO = " << rho << "\n";

  n = 1000;

  pr = new double[n];
  ur = new double[n];
  vr = new double[n];

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, seed );
  t = 0.0;

  resid_lucas ( nu, rho, n, x, y, t, ur, vr, pr );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  Ur:"
       << "  " << setw(14) << r8vec_amin ( n, ur )
       << "  " << setw(14) << r8vec_amax ( n, ur ) << "\n";
  cout << "  Vr:"
       << "  " << setw(14) << r8vec_amin ( n, vr )
       << "  " << setw(14) << r8vec_amax ( n, vr ) << "\n";
  cout << "  Pr:"
       << "  " << setw(14) << r8vec_amin ( n, pr )
       << "  " << setw(14) << r8vec_amax ( n, pr ) << "\n";

  delete [] pr;
  delete [] ur;
  delete [] vr;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void gnuplot_lucas_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GNUPLOT_LUCAS_TEST generates a field on a regular grid and plots it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string header;
  int n;
  double nu;
  double *p;
  double rho;
  double s;
  int seed;
  double t;
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
  cout << "GNUPLOT_LUCAS_TEST:\n";
  cout << "  Lucas Bystricky Flow:\n";
  cout << "  Generate a velocity field on a regular grid.\n";
  cout << "  Store in GNUPLOT data and command files.\n";

  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = new double[x_num*y_num];
  y = new double[x_num*y_num];

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  nu = 1.0;
  rho = 1.0;
  n = x_num * y_num;
  t = 0.0;

  u = new double[x_num*y_num];
  v = new double[x_num*y_num];
  p = new double[x_num*y_num];

  uvp_lucas ( nu, rho, n, x, y, t, u, v, p );

  header = "lucas";
  s = 0.25;
  ns2de_gnuplot ( header, n, x, y, u, v, s );

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
