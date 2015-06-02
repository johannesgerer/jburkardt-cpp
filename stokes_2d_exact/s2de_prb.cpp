# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "s2de.hpp"

int main ( );
void uvp_stokes1_test ( );
void resid_stokes1_test ( );
void gnuplot_stokes1_test ( );
void uvp_stokes2_test ( );
void resid_stokes2_test ( );
void gnuplot_stokes2_test ( );
void uvp_stokes3_test ( );
void resid_stokes3_test ( );
void gnuplot_stokes3_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    S2DE_PRB tests the S2DE library.
//
//  Location:
//
//    http://people.sc.fsu.edu/~jburkardt/cpp_src/stokes_2d_exact/s2de_prb.cpp
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "S2DE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the S2DE library.\n";

  uvp_stokes1_test ( );
  resid_stokes1_test ( );
  gnuplot_stokes1_test ( );

  uvp_stokes2_test ( );
  resid_stokes2_test ( );
  gnuplot_stokes2_test ( );

  uvp_stokes3_test ( );
  resid_stokes3_test ( );
  gnuplot_stokes3_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "S2DE_PRB\n";
  cout << "  Normal end of execution.\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void uvp_stokes1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_STOKES1_TEST samples the solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n = 1000;
  double *p;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "UVP_STOKES1_TEST\n";
  cout << "  Exact Stokes solution #1:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  using a region that is the unit square.\n";

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );

  u = new double[n];
  v = new double[n];
  p = new double[n];

  uvp_stokes1 ( n, x, y, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:" << "  " << setw(14) << r8vec_min ( n, u )
                 << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:" << "  " << setw(14) << r8vec_min ( n, v )
                 << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:" << "  " << setw(14) << r8vec_min ( n, p )
                 << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void resid_stokes1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_STOKES1_TEST samples the residual for solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n = 1000;
  double *pr;
  int seed;
  double *ur;
  double *vr;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "RESID_STOKES1_TEST\n";
  cout << "  Exact Stokes solution #1:\n";
  cout << "  Sample the Stokes residuals\n";
  cout << "  using a region that is the unit square.\n";

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );

  ur = new double[n];
  vr = new double[n];
  pr = new double[n];

  resid_stokes1 ( n, x, y, ur, vr, pr );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  Ur:" << "  " << setw(14) << r8vec_amin ( n, ur )
                  << "  " << setw(14) << r8vec_amax ( n, ur ) << "\n";
  cout << "  Vr:" << "  " << setw(14) << r8vec_amin ( n, vr )
                  << "  " << setw(14) << r8vec_amax ( n, vr ) << "\n";
  cout << "  Pr:" << "  " << setw(14) << r8vec_amin ( n, pr )
                  << "  " << setw(14) << r8vec_amax ( n, pr ) << "\n";
  delete [] pr;
  delete [] ur;
  delete [] vr;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void gnuplot_stokes1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GNUPLOT_STOKES1_TEST plots solution #1 on a regular grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  string header;
  int n;
  double *p;
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
  cout << "GNUPLOT_STOKES1_TEST:\n";
  cout << "  Exact Stokes solution #1:\n";
  cout << "  Generate a Stokes velocity field on a regular grid.\n";
  cout << "  Store in GNUPLOT data and command files.\n";
   
  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = new double[x_num*y_num];
  y = new double[x_num*y_num];

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;

  u = new double[n];
  v = new double[n];
  p = new double[n];

  uvp_stokes1 ( n, x, y, u, v, p );

  header = "stokes1";
  s = 4.0;
  stokes_gnuplot ( header, n, x, y, u, v, s );

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void uvp_stokes2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_STOKES2_TEST samples the solution #2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n = 1000;
  double *p;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "UVP_STOKES2_TEST\n";
  cout << "  Exact Stokes solution #2:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  using a region that is the unit square.\n";

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );

  u = new double[n];
  v = new double[n];
  p = new double[n];

  uvp_stokes2 ( n, x, y, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:" << "  " << setw(14) << r8vec_min ( n, u )
                 << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:" << "  " << setw(14) << r8vec_min ( n, v )
                 << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:" << "  " << setw(14) << r8vec_min ( n, p )
                 << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void resid_stokes2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_STOKES2_TEST samples the residual for solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n = 1000;
  double *pr;
  int seed;
  double *ur;
  double *vr;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "RESID_STOKES2_TEST\n";
  cout << "  Exact Stokes solution #2:\n";
  cout << "  Sample the Stokes residuals\n";
  cout << "  using a region that is the unit square.\n";

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );

  ur = new double[n];
  vr = new double[n];
  pr = new double[n];

  resid_stokes2 ( n, x, y, ur, vr, pr );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  Ur:" << "  " << setw(14) << r8vec_amin ( n, ur )
                  << "  " << setw(14) << r8vec_amax ( n, ur ) << "\n";
  cout << "  Vr:" << "  " << setw(14) << r8vec_amin ( n, vr )
                  << "  " << setw(14) << r8vec_amax ( n, vr ) << "\n";
  cout << "  Pr:" << "  " << setw(14) << r8vec_amin ( n, pr )
                  << "  " << setw(14) << r8vec_amax ( n, pr ) << "\n";

  delete [] pr;
  delete [] ur;
  delete [] vr;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void gnuplot_stokes2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GNUPLOT_STOKES2_TEST plots solution #2 on a regular grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  string header;
  int n;
  double *p;
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
  cout << "GNUPLOT_STOKES2_TEST:\n";
  cout << "  Exact Stokes solution #2:\n";
  cout << "  Generate a Stokes velocity field on a regular grid.\n";
  cout << "  Store in GNUPLOT data and command files.\n";
   
  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = new double[x_num*y_num];
  y = new double[x_num*y_num];

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;

  u = new double[n];
  v = new double[n];
  p = new double[n];

  uvp_stokes2 ( n, x, y, u, v, p );

  header = "stokes2";
  s = 0.05;
  stokes_gnuplot ( header, n, x, y, u, v, s );

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void uvp_stokes3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_STOKES3_TEST samples the solution #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n = 1000;
  double *p;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "UVP_STOKES3_TEST\n";
  cout << "  Exact Stokes solution #3:\n";
  cout << "  Estimate the range of velocity and pressure\n";
  cout << "  using a region that is [-1,+1]x[-1,+1].\n";

  xy_lo = -1.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );

  u = new double[n];
  v = new double[n];
  p = new double[n];

  uvp_stokes3 ( n, x, y, u, v, p );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  U:" << "  " << setw(14) << r8vec_min ( n, u )
                 << "  " << setw(14) << r8vec_max ( n, u ) << "\n";
  cout << "  V:" << "  " << setw(14) << r8vec_min ( n, v )
                 << "  " << setw(14) << r8vec_max ( n, v ) << "\n";
  cout << "  P:" << "  " << setw(14) << r8vec_min ( n, p )
                 << "  " << setw(14) << r8vec_max ( n, p ) << "\n";

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void resid_stokes3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_STOKES3_TEST samples the residual for solution #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n = 1000;
  double *pr;
  int seed;
  double *ur;
  double *vr;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  cout << "\n";
  cout << "RESID_STOKES3_TEST\n";
  cout << "  Exact Stokes solution #3:\n";
  cout << "  Sample the Stokes residuals\n";
  cout << "  using a region that is [-1,+1]x[-1,+1].\n";

  xy_lo = -1.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, seed );

  ur = new double[n];
  vr = new double[n];
  pr = new double[n];

  resid_stokes3 ( n, x, y, ur, vr, pr );

  cout << "\n";
  cout << "           Minimum       Maximum\n";
  cout << "\n";
  cout << "  Ur:" << "  " << setw(14) << r8vec_amin ( n, ur )
                  << "  " << setw(14) << r8vec_amax ( n, ur ) << "\n";
  cout << "  Vr:" << "  " << setw(14) << r8vec_amin ( n, vr )
                  << "  " << setw(14) << r8vec_amax ( n, vr ) << "\n";
  cout << "  Pr:" << "  " << setw(14) << r8vec_amin ( n, pr )
                  << "  " << setw(14) << r8vec_amax ( n, pr ) << "\n";
  delete [] pr;
  delete [] ur;
  delete [] vr;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void gnuplot_stokes3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GNUPLOT_STOKES3_TEST plots solution #3 on a regular grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  string header;
  int n;
  double *p;
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
  cout << "GNUPLOT_STOKES3_TEST:\n";
  cout << "  Exact Stokes solution #3:\n";
  cout << "  Generate a Stokes velocity field on [-1,+1]x[-1,+1].\n";
  cout << "  Store in GNUPLOT data and command files.\n";
   
  x_lo = -1.0;
  x_hi = +1.0;

  y_lo = -1.0;
  y_hi = +1.0;

  x = new double[x_num*y_num];
  y = new double[x_num*y_num];

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;

  u = new double[n];
  v = new double[n];
  p = new double[n];

  uvp_stokes3 ( n, x, y, u, v, p );

  header = "stokes3";
  s = 0.05;
  stokes_gnuplot ( header, n, x, y, u, v, s );

  delete [] p;
  delete [] u;
  delete [] v;
  delete [] x;
  delete [] y;

  return;
}
